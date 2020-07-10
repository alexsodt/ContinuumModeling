#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include <sys/time.h>
#include "interp.h"
#include "3d_interp.h"
#include "fitRho.h"
#include "simulation.h"
#include "parallel.h"
#include "mutil.h"
#include "M_matrix.h"

// fitRho globals
double fitCoupling = 1.0;
double fitThickness = 15.0;

// fitRho statics
static double eps_f_min = 1.0;
int fitRho_activated = 0;
void Simulation::setupDensity( char *fileName, int shiftRho )
{

	FILE *rhoFile = fopen(fileName, "r");
	char *buffer = (char *)malloc( sizeof(char) * 100000 );
	if( !rhoFile )
	{
		printf("Cannot open density file '%s'.\n", fileName );
		exit(1);
	}

	double La,Lb,Lc;
	int nx, ny, nz;
	
	fscanf(rhoFile, "%lf %lf %lf\n", &La, &Lb, &Lc );
	fscanf(rhoFile, "%d %d %d\n", &nx, &ny, &nz );

	double *rho = (double *)malloc( sizeof(double) * nx * ny * nz );


	for( int ix = 0; ix < nx; ix++ )
	for( int iy = 0; iy < ny; iy++ )
	{
		getLine( rhoFile, buffer );

//		for( int iz = 0; iz < nz; iz++ )
			readNDoubles( buffer, rho+ix*ny*nz+iy*nz, nz );
	}

	int cdist[3] = {0,0,0};
	int lims[3] = {nx,ny,nz};

	if( shiftRho )
		rhoShifter( rho, nx, ny, nz );

/*	{

		double distx[nx];
		memset( distx, 0, sizeof(double) * nx );
		double disty[ny];
		memset( disty, 0, sizeof(double) * ny );
		double distz[nz];
		memset( distz, 0, sizeof(double) * nz );

		double *dists[3] = { distx, disty, distz };

		for( int ix = 0; ix < nx; ix++ )
		for( int iy = 0; iy < ny; iy++ )
		for( int iz = 0; iz < nz; iz++ )
		{
			int bins[3] = { ix, iy, iz };

			distx[ix] += rho[ix*ny*nz+iy*nz+iz];	
			disty[iy] += rho[ix*ny*nz+iy*nz+iz];	
			distz[iz] += rho[ix*ny*nz+iy*nz+iz];	
		}

		for( int c = 0; c < 3; c++ )
		{
			double best_chi2 = 1e10;
			double best_cen = 0;
			double *dist = dists[c];

			for( int cb = 0; cb < lims[c]; cb++ )
			{
				double sigma2 = 0;

				for( int b = 0; b < lims[c]; b++ )
				{
					double dc = (b - cb);
					if( dc > lims[c]/2 ) dc -= lims[c]/2;
					if( dc < -lims[c]/2 ) dc += lims[c]/2;
					sigma2 += dist[b] * dc*dc;				
				}

				if( sigma2 < best_chi2 )
				{
					best_cen = cb;
					best_chi2 = sigma2;
				}
			}

			cdist[c] = best_cen; 	
		}

		printf("Cdist: %d %d %d\n", cdist[0], cdist[1], cdist[2] );

#if 1 
		double shift[3] = { cdist[0] * La / nx,
				    cdist[1] * Lb / ny, 
				    cdist[2] * Lc / nz }; 
		for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
		{
			surface *theSurface = sRec->theSurface;
			int nv = theSurface->nv;
			for( int v = 0; v < nv; v++ )
			{
				sRec->r[3*v+0] += shift[0];
				sRec->r[3*v+1] += shift[1];
				sRec->r[3*v+2] += shift[2];
			
				theSurface->theVertices[v].r[0] += shift[0];	
				theSurface->theVertices[v].r[1] += shift[1];	
				theSurface->theVertices[v].r[2] += shift[2];	
			}
		} 

#else
		// shift
		double *rhoCopy = (double *)malloc( sizeof(double) * nx * ny * nz );

		for( int ix = 0; ix < nx; ix++ )	
		for( int iy = 0; iy < ny; iy++ )	
		for( int iz = 0; iz < nz; iz++ )
		{
			int sx = ix - cdist[0];
			int sy = iy - cdist[1];
			int sz = iz - cdist[2];
		
			if( sx < 0 ) sx += nx;
			if( sx >= nx ) sx -= nx;
			if( sy < 0 ) sy += ny;
			if( sy >= ny ) sy -= ny;
			if( sz < 0 ) sz += nz;
			if( sz >= nz ) sz -= nz;

			rhoCopy[sx*nx*ny+sy*nz+sz] = rho[ix*ny*nz+iy*nz+iz];
		}

		memcpy( rho, rhoCopy, sizeof(double) * nx * ny *nz );	
#endif
	}
*/

	fclose(rhoFile);
	

	// #DBG BAD
	
	double LXS = PBC_vec[0][0];
	double LYS = PBC_vec[1][1];
	double LZS = PBC_vec[2][2];

	double scale_x = La / LXS;
	double scale_y = Lb / LYS;
	double scale_z = Lc / LYS;

	PBC_vec[0][0] = La;
	PBC_vec[1][1] = Lb;
	PBC_vec[2][2] = Lc;

	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		for( int i = 0; i < sRec->theSurface->nv; i++ )
		{
			sRec->theSurface->theVertices[i].r[0] *= scale_x;
			sRec->theSurface->theVertices[i].r[1] *= scale_y;
			sRec->theSurface->theVertices[i].r[2] *= scale_z;
		}
	}
	double *rho_smooth = (double *)malloc( sizeof(double) * nx * ny * nz );
	
	double xbw = La/nx;
	double ybw = Lb/ny;
	double zbw = Lc/nz;

	double sigma =  sqrt(xbw*xbw+ybw*ybw+zbw*zbw)*3;

#define RUN_SMOOTHER

	#define N_SMOOTH	10
	printf("Smoothing.\n");
	for( int siter = 0; siter < N_SMOOTH; siter++ )
	{
		for( int ix = 0; ix < nx; ix++ )
		for( int iy = 0; iy < ny; iy++ )
		for( int iz = 0; iz < nz; iz++ )
		{
			double r = 0;

			double sum = 0;

			for( int dx = -1; dx <= 1; dx++ )
			for( int dy = -1; dy <= 1; dy++ )
			for( int dz = -1; dz <= 1; dz++ )
			{
				double dr[3] = { dx * xbw, dy * ybw, dz * zbw };
				double r2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
					
				int rx = ix+dx;
				int ry = iy+dy;
				int rz = iz+dz;

				if( rx < 0 ) rx += nx;
				if( rx >= nx ) rx -= nx;

				if( ry < 0 )   ry += ny;
				if( ry >= ny ) ry -= ny;
				
				if( rz < 0 )   rz += nz;
				if( rz >= nz ) rz -= nz;

				int alt_bin = rx*ny*nz+ry*nz+rz;

				sum += exp( -r2 / sigma );	
			}

			for( int dx = -1; dx <= 1; dx++ )
			for( int dy = -1; dy <= 1; dy++ )
			for( int dz = -1; dz <= 1; dz++ )
			{
				double dr[3] = { dx * xbw, dy * ybw, dz * zbw };
				double r2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
					
				int rx = ix+dx;
				int ry = iy+dy;
				int rz = iz+dz;

				if( rx < 0 ) rx += nx;
				if( rx >= nx ) rx -= nx;

				if( ry < 0 )   ry += ny;
				if( ry >= ny ) ry -= ny;
				
				if( rz < 0 )   rz += nz;
				if( rz >= nz ) rz -= nz;

				int alt_bin = rx*ny*nz+ry*nz+rz;

				r += rho[alt_bin] * exp( -r2 / sigma ) / sum;	
			}

			rho_smooth[ix*ny*nz+iy*nz+iz] = r;
		}
		memcpy( rho, rho_smooth, sizeof(double) * nx * ny *nz );
	}
	// #DBG BAD
	printf("Done smoothing.\n");
	setup_rho( rho, nx, ny, nz );
	free(buffer);	
	

	fitRho_activated = 1;
}

/* This returns the overlap energy and computes the overlap gradient of the surface with a 3d interpolated density of neutral surface atoms. */

double surface::rhoEnergy( double *r, double PBC_vec[3][3], double thickness_inner, double thickness_outer )
{
#ifdef PARALLEL
	if( par_info.my_id != BASE_TASK )
		return 0;
#endif
	if( ! fitRho_activated ) return 0;
	if( fabs(fitCoupling) < 1e-30 ) return 0;
	return rhoWorker( r, NULL, PBC_vec, 0, thickness_inner, thickness_outer, NULL, NULL );
}

double surface::rhoGrad( double *r, double *gr, double PBC_vec[3][3], double thickness_inner, double thickness_outer, double *tDerInner, double *tDerOuter )
{
#ifdef PARALLEL
	if( par_info.my_id != BASE_TASK )
		return 0;
#endif

	if( ! fitRho_activated ) return 0;
	if( fabs(fitCoupling) < 1e-30 ) return 0;

	return rhoWorker( r, gr, PBC_vec, 1, thickness_inner, thickness_outer, tDerInner, tDerOuter );
}

double surface::rhoWorker( double * r, double *gr, double PBC_vec[3][3], int do_grad, double thickness_inner, double thickness_outer, double *tDerInner, double *tDerOuter)
{
	if( !theFormulas )
		generatePlan();

	int static print_trigger = 0;

	double alpha_x = r[3*nv];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];

	double e = 0;

	int pgrid = 1;

	int use_max = 20;

	double *rGrad = (double *)malloc( sizeof(double) * 3 * use_max );
	double *nGrad = (double *)malloc( sizeof(double) * 3 * 3 * use_max );
	double *hGrad = (double *)malloc( sizeof(double) * 3 * use_max );
	double *kGrad = (double *)malloc( sizeof(double) * 3 * use_max );
	int *coor_list = (int *)malloc( sizeof(int) * use_max );
	int nCoor;

	double der_thickness[2] = { 0, 0 };

	int ngrid_pts = pgrid * (pgrid+1)/2;
	
	FILE *outFile = NULL;		
//	if( print_trigger % 10 == 0 && do_grad )
//		outFile = fopen("tempout.xyz","w");

	double ATOT = 0;
	for( int f = 0; f < nt; f++ )
	{
		if( f < nf_faces )
		{
			for( int q = 0; q < nf_g_q_p; q++ )
				ATOT += theFormulas[f*nf_g_q_p+q].g0*0.5;
		}
		else
		{
			for( int q = 0; q < nf_irr_pts; q++ )
				ATOT += theIrregularFormulas[(f-nf_faces)*nf_irr_pts+q].g0 * theIrregularFormulas[(f-nf_faces)*nf_irr_pts+q].weight;
		}
	}

	int nuv = pgrid*pgrid; 
	
	//
	double uv_array[2*nuv];
	memset(uv_array, 0, sizeof(double) * 2 * nuv );
	double l = 1.0 / (double)pgrid;

	// advances per row:
	double d_row =  l * sin(60.0*M_PI/180.0);
	// advances per triangle;
	double d_row_half =  l * sin(30.0*M_PI/180.0);
	
	double uoff = d_row_half;
	int uv_pts = 0;
	for( int row = 0; row < pgrid; row++ )
	{
		double voff = d_row_half;

		for( int j = 0; j < pgrid-row; j++ )
		{
			uv_array[2*uv_pts+0] = row * l + (1.0/3.0) * l;
			uv_array[2*uv_pts+1] = j * l + (1.0/3.0) * l;
//			printf("%lf %lf\n", uv_array[2*uv_pts+0], uv_array[2*uv_pts+1]);
			uv_pts++;

			if( row > 0 )
			{
				uv_array[2*uv_pts+0] = row * l - (1.0/3.0) * l;
				uv_array[2*uv_pts+1] =   j * l + (2.0/3.0) * l;
//				printf("%lf %lf\n", uv_array[2*uv_pts+0], uv_array[2*uv_pts+1]);
				uv_pts++;
			}
		} 
	}


	int max_f=0;
	double max_u=1e10,max_v=1e10;
	double max_strain = 0;
	double max_c=1e10, max_k=1e10;

	for( int f = 0; f < nt; f++ )
	{
		double g0 = 0;
	
		if( f < nf_faces )
		{
			for( int q = 0; q < nf_g_q_p; q++ )
				g0 += theFormulas[f*nf_g_q_p+q].g0*0.5;

			double *rf = r + 3 * theFormulas[f*nf_g_q_p].cp[0];
		
			double dr[3] = { rf[0],rf[1],rf[2]-40 };

			double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

			if( r < 40.0 )
			{
//				printf("Break here.\n");
			}
		}
		else
		{
			for( int q = 0; q < nf_irr_pts; q++ )
				g0 += theIrregularFormulas[(f-nf_faces)*nf_irr_pts+q].g0 * theIrregularFormulas[(f-nf_faces)*nf_irr_pts+q].weight;
		}

		g0 /= ATOT;


		for( int tuv = 0; tuv < uv_pts; tuv++ )
		{
			double u = uv_array[2*tuv+0];
			double v = uv_array[2*tuv+1];

			double ctot=0;
			double k = 0;
			if( do_grad )
			{
				ctot = gradFetch( f, u, v, r,
				rGrad,
				nGrad,
				hGrad,
				kGrad,
				&nCoor,
				coor_list, &k );
			}
			else
				ctot = c( f, u, v, r, &k );
//#define FDIFF_CHECK		
#ifdef FDIFF_CHECK
{
			double eps = 1e-4;

			for( int cx = 0; cx < nCoor; cx++ )
			{
				for( int cart = 0; cart < 3; cart++ )
				{
					double dc[2];
					double dk[2];
					double dn[2][3];
					double dr[2][3];

					for( int pm = 0; pm < 2; pm++ )
					{
						r[3*coor_list[cx]+cart] += eps * (pm == 0 ? -1 : 1 );
						double k=0;
						dc[pm] = c(f,u,v,r,&k);
						dk[pm] = k;
						evaluateRNRM( f, u ,v, dr[pm], dn[pm], r ); 

						r[3*coor_list[cx]+cart] -= eps * (pm == 0 ? -1 : 1 );
					}
	
					double fd_k = (dk[1]-dk[0])/(2*eps);	
					double fd_c = (dc[1]-dc[0])/(2*eps);
					double fd_r[3] = { 
						(dr[1][0]-dr[0][0])/(2*eps),
						(dr[1][1]-dr[0][1])/(2*eps),
						(dr[1][2]-dr[0][2])/(2*eps) };
					double fd_n[3] = { 
						(dn[1][0]-dn[0][0])/(2*eps),
						(dn[1][1]-dn[0][1])/(2*eps),
						(dn[1][2]-dn[0][2])/(2*eps) };

					printf("f %d cart %d cx %d fdc: %.14le dc: %.14le fdk: %.14le dk: %.14le ac: %.14le fdr: %.14le ar: %.14le , fdn: %.14le %.14le %.14le an: %.14le %.14le %.14le\n",
						f,cart,cx,
						fd_c, 
						hGrad[3*cx+cart],
						fd_k,
						kGrad[3*cx+cart],
						fd_r[cart],
						rGrad[3*cx+cart],

						fd_n[0], fd_n[1], fd_n[2],
						nGrad[9*cx+cart+3*0], nGrad[9*cx+cart+3*1], nGrad[9*cx+cart+3*2] );
					
				}
			}
}			
#endif
			double rpt[3],nrm[3];

			evaluateRNRM( f, u, v, rpt, nrm, r );

			double use_thickness[2] = { thickness_inner, thickness_outer };
			
			double signs[2] = {-1,1};
			for( int leaflet = 0; leaflet < 2; leaflet++ )
			{
				double use_sign = signs[leaflet];

				double thickness = use_thickness[leaflet];
				double strain_f = exp( use_sign * thickness * ctot + thickness*thickness*k );

				double R[3] = { 
					rpt[0] + use_sign * nrm[0] * thickness* strain_f,
					rpt[1] + use_sign * nrm[1] * thickness* strain_f,
					rpt[2] + use_sign * nrm[2] * thickness* strain_f};

				double strain = strain_f - 1;

				if( strain_f > 200 )
				{
					e += 1e10;
					continue;
				}

				if( fabs(strain) > fabs(max_strain) )
				{
					max_f = f;
					max_u = u;
					max_v = v;
					max_c = ctot;
					max_k = k;
					max_strain = strain;
				}
				double d_R_d_thickness[3] = { 
						use_sign * nrm[0] * strain_f + use_sign * nrm[0] * thickness * strain_f * ( use_sign * ctot + 2 * thickness * k),
						use_sign * nrm[1] * strain_f + use_sign * nrm[1] * thickness * strain_f * ( use_sign * ctot + 2 * thickness * k),
						use_sign * nrm[2] * strain_f + use_sign * nrm[2] * thickness * strain_f * ( use_sign * ctot + 2 * thickness * k) };
//						use_sign * nrm[0] * (1+use_sign * thickness*ctot) + nrm[0] * ctot * thickness + use_sign * nrm[0] * 3 * thickness * thickness * k,
//						use_sign * nrm[1] * (1+use_sign * thickness*ctot) + nrm[1] * ctot * thickness + use_sign * nrm[1] * 3 * thickness * thickness * k,
//						use_sign * nrm[2] * (1+use_sign * thickness*ctot) + nrm[2] * ctot * thickness + use_sign * nrm[2] * 3 * thickness * thickness * k };

				double fractional[3] = { R[0] / PBC_vec[0][0], R[1] / PBC_vec[1][1], R[2] / PBC_vec[2][2] };

				e += -(g0/ngrid_pts) * fitCoupling * log( eps_f_min + eval_rho( fractional[0], fractional[1], fractional[2] ) );

				if( print_trigger % 100 == 0 && outFile )
				{
					fprintf(outFile, "%c %lf %lf %lf\n", (leaflet == 0 ? 'O' : 'C' ), R[0], R[1],R[2] );

				}
				if( do_grad )
				{
					double arg = eps_f_min + eval_rho( fractional[0], fractional[1], fractional[2] );

					double rho_g[3] = {0,0,0};
					eval_drho( fractional[0], fractional[1], fractional[2], rho_g );
	
					for( int cx = 0; cx < nCoor; cx++ )
					{
						double d_Rx_d_vx = rGrad[3*cx+0];
						double d_Rx_d_vy = 0;
						double d_Rx_d_vz = 0;
	
						double d_Ry_d_vx = 0;
						double d_Ry_d_vy = rGrad[3*cx+1];
						double d_Ry_d_vz = 0;
						
						double d_Rz_d_vx = 0;
						double d_Rz_d_vy = 0;
						double d_Rz_d_vz = rGrad[3*cx+2];
						d_Rx_d_vx += use_sign * nGrad[9*cx+0*3+0] * thickness * strain_f;
						d_Ry_d_vx += use_sign * nGrad[9*cx+1*3+0] * thickness * strain_f; 
						d_Rz_d_vx += use_sign * nGrad[9*cx+2*3+0] * thickness * strain_f; 
						
						d_Rx_d_vy += use_sign * nGrad[9*cx+0*3+1] * thickness * strain_f; 
						d_Ry_d_vy += use_sign * nGrad[9*cx+1*3+1] * thickness * strain_f; 
						d_Rz_d_vy += use_sign * nGrad[9*cx+2*3+1] * thickness * strain_f; 
						
						d_Rx_d_vz += use_sign * nGrad[9*cx+0*3+2] * thickness * strain_f; 
						d_Ry_d_vz += use_sign * nGrad[9*cx+1*3+2] * thickness * strain_f; 
						d_Rz_d_vz += use_sign * nGrad[9*cx+2*3+2] * thickness * strain_f; 
						
						d_Rx_d_vx += use_sign * nrm[0] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+0]; 
						d_Ry_d_vx += use_sign * nrm[1] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+0]; 
						d_Rz_d_vx += use_sign * nrm[2] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+0]; 
						                                                                             
						d_Rx_d_vy += use_sign * nrm[0] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+1]; 
						d_Ry_d_vy += use_sign * nrm[1] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+1]; 
						d_Rz_d_vy += use_sign * nrm[2] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+1]; 
						                                                                             
						d_Rx_d_vz += use_sign * nrm[0] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+2]; 
						d_Ry_d_vz += use_sign * nrm[1] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+2]; 
						d_Rz_d_vz += use_sign * nrm[2] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+2]; 
						
						d_Rx_d_vx += use_sign * nrm[0] * thickness * strain_f * thickness * thickness * kGrad[3*cx+0]; 
						d_Ry_d_vx += use_sign * nrm[1] * thickness * strain_f * thickness * thickness * kGrad[3*cx+0]; 
						d_Rz_d_vx += use_sign * nrm[2] * thickness * strain_f * thickness * thickness * kGrad[3*cx+0]; 
						                                                                                             
						d_Rx_d_vy += use_sign * nrm[0] * thickness * strain_f * thickness * thickness * kGrad[3*cx+1]; 
						d_Ry_d_vy += use_sign * nrm[1] * thickness * strain_f * thickness * thickness * kGrad[3*cx+1]; 
						d_Rz_d_vy += use_sign * nrm[2] * thickness * strain_f * thickness * thickness * kGrad[3*cx+1]; 
						                                                                                             
						d_Rx_d_vz += use_sign * nrm[0] * thickness * strain_f * thickness * thickness * kGrad[3*cx+2]; 
						d_Ry_d_vz += use_sign * nrm[1] * thickness * strain_f * thickness * thickness * kGrad[3*cx+2]; 
						d_Rz_d_vz += use_sign * nrm[2] * thickness * strain_f * thickness * thickness * kGrad[3*cx+2]; 

						gr[3*coor_list[cx]+0] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[0] / PBC_vec[0][0] * alpha_x * d_Rx_d_vx;	
						gr[3*coor_list[cx]+1] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[0] / PBC_vec[0][0] * alpha_x * d_Rx_d_vy;	
						gr[3*coor_list[cx]+2] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[0] / PBC_vec[0][0] * alpha_x * d_Rx_d_vz;	
						
						gr[3*coor_list[cx]+0] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[1] / PBC_vec[1][1] * alpha_y * d_Ry_d_vx;	
						gr[3*coor_list[cx]+1] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[1] / PBC_vec[1][1] * alpha_y * d_Ry_d_vy;	
						gr[3*coor_list[cx]+2] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[1] / PBC_vec[1][1] * alpha_y * d_Ry_d_vz;	
						
						gr[3*coor_list[cx]+0] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[2] / PBC_vec[2][2] * alpha_z * d_Rz_d_vx;	
						gr[3*coor_list[cx]+1] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[2] / PBC_vec[2][2] * alpha_z * d_Rz_d_vy;	
						gr[3*coor_list[cx]+2] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[2] / PBC_vec[2][2] * alpha_z * d_Rz_d_vz;	

					}
						
					der_thickness[leaflet] -=  fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[0] / PBC_vec[0][0] * alpha_x  * d_R_d_thickness[0];
					der_thickness[leaflet] -=  fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[1] / PBC_vec[1][1] * alpha_y  * d_R_d_thickness[1];
					der_thickness[leaflet] -=  fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[2] / PBC_vec[2][2] * alpha_z  * d_R_d_thickness[2];
				}
			}	
		}	
	} 

	free(rGrad);
	free(nGrad);
	free(hGrad);
	// thickness penalty to prohibit locking onto a single leaflet.

#if 1	
	double massive_k = 1e5;
	double min_thresh = 8.0;
	double max_thresh = 25.0;
	if( thickness_inner < min_thresh )
	{
		double dh = (thickness_inner-min_thresh);
		e += massive_k * (dh*dh);	
		if( do_grad )
			*tDerInner += 2 * massive_k * dh;
	}
	
	if( thickness_outer < min_thresh )
	{
		double dh = (thickness_outer-min_thresh);
		e += massive_k * (dh*dh);	

		if( do_grad )
			*tDerOuter += 2 * massive_k * dh;
	}
	
	if( thickness_inner > max_thresh )
	{
		double dh = (max_thresh-thickness_inner);
		e += massive_k * (dh*dh);	
		if( do_grad )
			*tDerInner -= 2 * massive_k * dh;
	}
	
	if( thickness_outer > max_thresh )
	{
		double dh = (max_thresh-thickness_outer);
		e += massive_k * (dh*dh);	

		if( do_grad )
			*tDerOuter -= 2 * massive_k * dh;
	}
#endif
	if( do_grad )//&& max_strain > 0 )
	{
		printf("max strain %lf %lf %lf f %d u %lf v %lf\n", max_strain, max_c, max_k, max_f, max_u, max_v );
		*tDerInner += der_thickness[0];
		*tDerOuter += der_thickness[1];
	
		print_trigger++;

	}

	if( outFile )
	{
//		fclose(outFile);
	}

	return e;	
}

/*
 * finds the best translational shift of the density by doing 1-dimensional projection fits.
 *
 * */

void Simulation::rhoShifter( double *rho, int nx, int ny, int nz )
{
	double distx[nx];
	memset( distx, 0, sizeof(double) * nx );
	double disty[ny];
	memset( disty, 0, sizeof(double) * ny );
	double distz[nz];
	memset( distz, 0, sizeof(double) * nz );

	double *dists[3] = { distx, disty, distz };

	for( int ix = 0; ix < nx; ix++ )
	for( int iy = 0; iy < ny; iy++ )
	for( int iz = 0; iz < nz; iz++ )
	{
		int bins[3] = { ix, iy, iz };

		distx[ix] += rho[ix*ny*nz+iy*nz+iz];	
		disty[iy] += rho[ix*ny*nz+iy*nz+iz];	
		distz[iz] += rho[ix*ny*nz+iy*nz+iz];	
	}

	double La = PBC_vec[0][0];
	double Lb = PBC_vec[1][1];
	double Lc = PBC_vec[2][2];

	double meshx[nx];
	double meshy[ny];
	double meshz[nz];

	memset( meshx, 0, sizeof(double)*nx);
	memset( meshy, 0, sizeof(double)*ny);
	memset( meshz, 0, sizeof(double)*nz);
	
	double *meshes[3] = { meshx, meshy, meshz };
	int lims[3] = {nx,ny,nz};

	// loop over the surfaces and their triangles, make a rough estimate of the density.
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{	
		surface *theSurface = sRec->theSurface;

		int grain = 10;
		for( int t = 0; t < theSurface->nt; t++ )
		{
			for( int iu = 0; iu < grain; iu ++ )
			for( int iv = 0; iv < grain-iu; iv++ )
			{
				double rpt[3],npt[3];
	
				double u = (iu+0.5)/grain;
				double v = (iv+0.5)/grain;

				theSurface->evaluateRNRM( t, u, v, rpt, npt, sRec->r );
	
				double lipds = theSurface->g( t, u, v, sRec->r );

				while( rpt[0] < 0 ) rpt[0] += PBC_vec[0][0];				
				while( rpt[0] >= PBC_vec[0][0] ) rpt[0] -= PBC_vec[0][0];				
				while( rpt[1] < 0 ) rpt[1] += PBC_vec[1][1];				
				while( rpt[1] >= PBC_vec[1][1] ) rpt[1] -= PBC_vec[1][1];				
				while( rpt[2] < 0 ) rpt[2] += PBC_vec[2][2];				
				while( rpt[2] >= PBC_vec[2][2] ) rpt[2] -= PBC_vec[2][2];				
	
				int  ix = nx * rpt[0] / PBC_vec[0][0];
				int  iy = ny * rpt[1] / PBC_vec[1][1];
				int  iz = nz * rpt[2] / PBC_vec[2][2];
	
				meshx[ix] += lipds;
				meshy[iy] += lipds;
				meshz[iz] += lipds;
			}
		}
	}
	
	int cdist[3] = {0,0,0};

	// for each dimension, find the shift that leads to the best match of the density profile.
	for( int c = 0; c < 3; c++ )
	{
//		printf("DIMENSION %d\n", c );
		double best_chi2 = 1e100;
		double best_cen = 0;
		double *dist = dists[c];
		double *mesh = meshes[c];
	
		// cb is the center of the box.
		for( int cb = 0; cb < lims[c]; cb++ )
		{
			double sigma2 = 0;

			// find the scaling of the mesh density that gives the best match.
			
			double xy = 0;
			double xx = 0;

			for( int b = 0; b < lims[c]; b++ )
			{
				int tb = b + cb;
				while( tb < 0 ) tb += lims[c];
				while( tb >= lims[c] ) tb -= lims[c];
				xx += mesh[tb]*mesh[tb];
				xy += mesh[tb]*dist[b];
			}

			double m = xy/xx;

			for( int b = 0; b < lims[c]; b++ )
			{
				int tb = b + cb;
				while( tb < 0 ) tb += lims[c];
				while( tb >= lims[c] ) tb -= lims[c];
				
				double del = dist[b] - m * mesh[tb];

				if( cb == 0 ) printf("%d %lf %lf\n", b, dist[b], mesh[tb] );

				sigma2 += del*del;	
			}

			if( sigma2 < best_chi2 )
			{
				best_cen = cb;
				best_chi2 = sigma2;
			}
		}

		cdist[c] = best_cen; 	
	}

	printf("Rho shifter Cdist: %d %d %d\n", -cdist[0], -cdist[1], -cdist[2] );

	double shift[3] = { -cdist[0] * La / nx,
			    -cdist[1] * Lb / ny, 
			    -cdist[2] * Lc / nz }; 

	double PBCs[3] = { La, Lb, Lc };
	for( int c = 0; c < 3; c++ )
	{
		while( shift[c] >  PBCs[c]/2 ) shift[c] -= PBCs[c];
		while( shift[c] < -PBCs[c]/2 ) shift[c] += PBCs[c];
	}
	
	printf("Rho shifter Cdist: %d %d %d, shift: %le %le %le\n", cdist[0], cdist[1], cdist[2], shift[0], shift[1], shift[2] );
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		surface *theSurface = sRec->theSurface;
		int nv = theSurface->nv;
		for( int v = 0; v < nv; v++ )
		{
			sRec->r[3*v+0] += shift[0];
			sRec->r[3*v+1] += shift[1];
			sRec->r[3*v+2] += shift[2];
		
			theSurface->theVertices[v].r[0] += shift[0];	
			theSurface->theVertices[v].r[1] += shift[1];	
			theSurface->theVertices[v].r[2] += shift[2];	
		}
	} 
}

// FIXED CUTS: some point on the surface must go through the set of points.

static int n_cut_points = 5;

int surface::movePointToCut( int *f_in, double *u_in, double *v_in, int cartesian_component, double value, double *r)
{
	int f = *f_in;
	double u = *u_in;	
	double v = *v_in;

	int found_pt = 0;

	do {
		double drdu[3], drdv[3];

		double rpt[3],nrm[3];
		// find the point
		evaluateRNRM( f, u, v, rpt, nrm,  r );			
		// find the tangent plane
		ru( f, u, v, r, drdu );
		rv( f, u, v, r, drdv );
		// find the uv-vector that moves us in the proper direction.

		double dsign = 1;

		if( rpt[cartesian_component] < value ) // needs to move up
			dsign = 1;
		else
			dsign = -1;	

		// find drdu/drdv combination that maximizes dcart.
		double kv_on_ku = 1;

		double rux = drdu[0], ruy = drdu[1], ruz = drdu[2];
		double rvx = drdv[0], rvy = drdv[1], rvz = drdv[2];

		switch( cartesian_component )
		{
			case 0:	
				kv_on_ku = -(( ruy*ruy*rvx  + ruz*ruz*rvx - rux*ruy*rvy - rux*ruz*rvz)/(  ruy*rvx*rvy -  rux*rvy*rvy + ruz*rvx*rvz - rux*rvz*rvz));
				break;	
			case 1:
				kv_on_ku = -((-(rux*ruy*rvx) + rux*rux*rvy + ruz*ruz*rvy - ruy*ruz*rvz)/(-(ruy*rvx*rvx) + rux*rvx*rvy + ruz*rvy*rvz - ruy*rvz*rvz));	
				break;	
			case 2:
				kv_on_ku = -((-(rux*ruz*rvx) - ruy*ruz*rvy + rux*rux*rvz + ruy*ruy*rvz)/(-(ruz*rvx*rvx) - ruz*rvy*rvy + rux*rvx*rvz + ruy*rvy*rvz));
				break;	
		}

		double duv[2] = { 1, kv_on_ku };

		double move[3] = { duv[0] * drdu[0] + duv[1] * drdv[0],
				   duv[0] * drdu[1] + duv[1] * drdv[1],
				   duv[0] * drdu[2] + duv[1] * drdv[2] };

		double ln = sqrt(duv[0]*duv[0]+duv[1]*duv[1]);

		duv[0] /= ln;
		duv[1] /= ln;

		double move_expec = normalize(move);
		if( fabs(move[cartesian_component])  < 0.05 )
		{
			found_pt = -1;		
		} 
		else
		{
//			printf("Going for it.\n");
		}
		if( move[cartesian_component] * (rpt[cartesian_component]-value) > 0 )
		{
			duv[0] *= -1;
			duv[1] *= -1;
		}

		int f_o = f;
		double uv_o[2] = { u,v};
		double duv_o[2] = { duv[0], duv[1] };

		double alpha_low  = 0.0;
		double alpha_high = 1.0;

		// check to see if we end up past the point
		
		int nf;
		do {
			nf = f;
			f = nextFace( f, &u, &v, duv+0, duv+1, r );
		} while( nf != f );	

		double pt2[3];
		evaluateRNRM( f, u, v, pt2, nrm, r );

		double move_length = sqrt( 
				(pt2[0]-rpt[0])*(pt2[0]-rpt[0])+
				(pt2[1]-rpt[1])*(pt2[1]-rpt[1])+
				(pt2[2]-rpt[2])*(pt2[2]-rpt[2]) );
		double move_expec_o = move_expec;
		if( (pt2[cartesian_component] - value) * ( rpt[cartesian_component] - value) < 0 )
		{
			// moved past, bisect to find the point.
			
			double alpha_low  = 0.0;
			double alpha_high = 1.0;
			double alpha_test = 0.5;

			int niters = 0;
			do { 
				duv[0] = duv_o[0] * alpha_test;						
				duv[1] = duv_o[1] * alpha_test;						
				f = f_o;
				u = uv_o[0];		
				v = uv_o[1];
	
				int nf;
				do {
					nf = f;
					f = nextFace( f, &u, &v, duv+0, duv+1, r );
				} while( nf != f );	
			
				evaluateRNRM( f, u, v, pt2, nrm, r );	
		

//				printf("alphas: %.14le %.14le %.14le cur: %.14le\n", alpha_low, alpha_test, alpha_high, pt2[cartesian_component] );
	
				if( (pt2[cartesian_component] - value) * ( rpt[cartesian_component] - value) < 0 ) // moved past
				{
					alpha_high = alpha_test;
					alpha_test = (alpha_low+alpha_test)/2;
				}
				else 
				{
					alpha_low = alpha_test;
					alpha_test = (alpha_high+alpha_test)/2;
				}
	
				if( fabs(pt2[cartesian_component]-value) < 1e-3)
					found_pt = 1;
				niters++;
			} while( found_pt != 1 && niters < 100 );
		}
		else 
		{	
			double scalef = 0.5;
			while( move_length/move_expec < 0.5 )
			{
				duv[0] = duv_o[0] * scalef;
				duv[1] = duv_o[1] * scalef;
				
				double move[3] = { duv[0] * drdu[0] + duv[1] * drdv[0],
						   duv[0] * drdu[1] + duv[1] * drdv[1],
						   duv[0] * drdu[2] + duv[1] * drdv[2] };
				move_expec = normalize(move);

				f = f_o;
				u = uv_o[0];
				v = uv_o[1];

				int nf;
				do {
					nf = f;
					f = nextFace( f, &u, &v, duv+0, duv+1, r );
				} while( nf != f );	
		
				evaluateRNRM( f, u, v, pt2, nrm, r );

				move_length = sqrt( 
					(pt2[0]-rpt[0])*(pt2[0]-rpt[0])+
					(pt2[1]-rpt[1])*(pt2[1]-rpt[1])+
					(pt2[2]-rpt[2])*(pt2[2]-rpt[2]) );

				scalef *= 0.5;
			}
		}					
	} while( found_pt == 0 );
	
	if( found_pt == 1 )
	{
		*f_in = f;
		*u_in = u;
		*v_in = v;
	}


	return found_pt;
}

void surface::get_cut_points( int cartesian_component, double value, int *f_pts, double *uv_pts, double *rall, int n_cut_points, double *r)
{
	for( int p = 0; p < n_cut_points; p++ )
	{
		int done = 0;

		while( !done )
		{
			int f;
			double u,v;
			randomPointOnSurface( &f, &u, &v );
	
			int found_pt = movePointToCut(  &f, &u, &v, cartesian_component, value, r );

			if( found_pt == 1 )
			{
				f_pts[p] = f;
				uv_pts[p*2+0] = u;
				uv_pts[p*2+1] = v;
				done = 1;

				double rpt[3], nrm[3];
				evaluateRNRM( f, u, v, rpt, nrm, r );


				rall[p*3+0] = rpt[0];
				rall[p*3+1] = rpt[1];
				rall[p*3+2] = rpt[2];
			}
		}
	}		

	int n_mc = 30 * n_cut_points;

	for( int mc = 0; mc < n_mc; mc++ )
	{
		double cure = 0;

		for( int p = 0; p < n_cut_points; p++ )
		for( int p2 = p+1; p2 < n_cut_points; p2++ )
		{
			double dr[3] = { rall[3*p+0] - rall[3*p2+0], rall[3*p+1] - rall[3*p2+1], rall[3*p+2] - rall[3*p2+2] };

			double l = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

			cure += 1.0/l/l;	
		}

		int to_move = rand() % n_cut_points;	
	
		double rsave[3] = { rall[3*to_move+0], rall[3*to_move+1], rall[3*to_move+2] };
		int f_save = f_pts[to_move];
		double u_save = uv_pts[2*to_move];
		double v_save = uv_pts[2*to_move+1];

		double sigma = 5.0;

		double frc_duv[2] = { 0,0};
		double dt = 1;
		double fstep;

		int f = f_save;
		double u = u_save;
		double v = v_save;

		localMove( &f, &u, &v, sigma, r, frc_duv, dt, &fstep, 0 ); 

		int found_pt = movePointToCut(  &f, &u, &v, cartesian_component, value, r );
		
		if( found_pt == 1 )
		{
			double rpt[3], nrm[3];
			evaluateRNRM( f, u, v, rpt, nrm, r );

			rall[3*to_move+0] = rpt[0]; 
			rall[3*to_move+1] = rpt[1]; 
			rall[3*to_move+2] = rpt[2]; 
			double newe = 0;

			for( int p = 0; p < n_cut_points; p++ )
			for( int p2 = p+1; p2 < n_cut_points; p2++ )
			{
				double dr[3] = { rall[3*p+0] - rall[3*p2+0], rall[3*p+1] - rall[3*p2+1], rall[3*p+2] - rall[3*p2+2] };
	
				double l = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
	
				newe += 1.0/l/l;	
			}

			if( newe < cure )
			{
				f_pts[to_move] = f;
				uv_pts[2*to_move] = u;
				uv_pts[2*to_move+1] = v;

				cure = newe;
			}
			else
			{
				rall[3*to_move+0] = rsave[0];
				rall[3*to_move+1] = rsave[1];
				rall[3*to_move+2] = rsave[2];
			}
		}

	}
	
}

void surface::addFixedPoint( double *r_fixed )
{
	fixed_cut_point *new_point = (fixed_cut_point *)malloc( sizeof(fixed_cut_point) );
	new_point->rpt[0] = r_fixed[0];
	new_point->rpt[1] = r_fixed[1];
	new_point->rpt[2] = r_fixed[2];
	new_point->k = 10; // 10 kcal/mol/a^2

	new_point->next = cutPoints;
	cutPoints = new_point; 
}

void surface::setupCut( int cartesian_component, double value, double *r)
{
	// here's how we'll get them. find n_cut_points random points on the surface at the cut value
	// use a pseudopotential to optimize their position at the cut value to try to make them evenly spaced.		
	int *f_pts = (int *)malloc( sizeof(int) * n_cut_points );
	double *uv_pts = (double *)malloc( sizeof(double) * 2 * n_cut_points );	
	double *rall = (double *)malloc( sizeof(double) * 3 * n_cut_points );

	get_cut_points( cartesian_component, value, f_pts, uv_pts, rall, n_cut_points, r );

#if 0
	for( int p = 0; p < n_cut_points; p++ )
	{
		int done = 0;

		while( !done )
		{
			int f;
			double u,v;
			randomPointOnSurface( &f, &u, &v );
	
			int found_pt = movePointToCut(  &f, &u, &v, cartesian_component, value, r );

			if( found_pt == 1 )
			{
				f_pts[p] = f;
				uv_pts[p*2+0] = u;
				uv_pts[p*2+1] = v;
				done = 1;

				double rpt[3], nrm[3];
				evaluateRNRM( f, u, v, rpt, nrm, r );


				rall[p*3+0] = rpt[0];
				rall[p*3+1] = rpt[1];
				rall[p*3+2] = rpt[2];
			}
		}
	}		

	int n_mc = 1000;

	for( int mc = 0; mc < n_mc; mc++ )
	{
		double cure = 0;

		for( int p = 0; p < n_cut_points; p++ )
		for( int p2 = p+1; p2 < n_cut_points; p2++ )
		{
			double dr[3] = { rall[3*p+0] - rall[3*p2+0], rall[3*p+1] - rall[3*p2+1], rall[3*p+2] - rall[3*p2+2] };

			double l = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

			cure += 1.0/l/l;	
		}

		int to_move = rand() % n_cut_points;	
	
		double rsave[3] = { rall[3*to_move+0], rall[3*to_move+1], rall[3*to_move+2] };
		int f_save = f_pts[to_move];
		double u_save = uv_pts[2*to_move];
		double v_save = uv_pts[2*to_move+1];

		double sigma = 5.0;

		double frc_duv[2] = { 0,0};
		double dt = 1;
		double fstep;

		int f = f_save;
		double u = u_save;
		double v = v_save;

		localMove( &f, &u, &v, sigma, r, frc_duv, dt, &fstep, 0 ); 

		int found_pt = movePointToCut(  &f, &u, &v, cartesian_component, value, r );
		
		if( found_pt == 1 )
		{
			double rpt[3], nrm[3];
			evaluateRNRM( f, u, v, rpt, nrm, r );

			rall[3*to_move+0] = rpt[0]; 
			rall[3*to_move+1] = rpt[1]; 
			rall[3*to_move+2] = rpt[2]; 
			double newe = 0;

			for( int p = 0; p < n_cut_points; p++ )
			for( int p2 = p+1; p2 < n_cut_points; p2++ )
			{
				double dr[3] = { rall[3*p+0] - rall[3*p2+0], rall[3*p+1] - rall[3*p2+1], rall[3*p+2] - rall[3*p2+2] };
	
				double l = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
	
				newe += 1.0/l/l;	
			}

			if( newe < cure )
			{
				f_pts[to_move] = f;
				uv_pts[2*to_move] = u;
				uv_pts[2*to_move+1] = v;

				cure = newe;
			}
			else
			{
				rall[3*to_move+0] = rsave[0];
				rall[3*to_move+1] = rsave[1];
				rall[3*to_move+2] = rsave[2];
			}
		}

		printf("cure: %le\n", cure );
	}

//	for( int t = 0; t < n_cut_points; t++ )
//		printf("r: %le %le %le\n", rall[3*t+0], rall[3*t+1], rall[3*t+2] );

#endif

	for( int t = 0; t < n_cut_points; t++ )
	{
		fixed_cut_point *new_point = (fixed_cut_point *)malloc( sizeof(fixed_cut_point) );
		new_point->rpt[0] = rall[3*t+0];
		new_point->rpt[1] = rall[3*t+1];
		new_point->rpt[2] = rall[3*t+2];
		new_point->k = 10; // 10 kcal/mol/a^2

		new_point->next = cutPoints;
		cutPoints = new_point; 
	}

	free(f_pts);
	free(uv_pts);
	free(rall);	
}

double surface::cutEnergy( double *r )
{
	double **M;
	int mlow,mhigh;
	getM(&M,&mlow,&mhigh);

	double cute = 0;

	for( fixed_cut_point *thePt = cutPoints; thePt; thePt = thePt->next )
	{
		double r_compr[3] = { thePt->rpt[0], thePt->rpt[1], thePt->rpt[2] };

		int f;
		double u, v;
		double distance;
		nearPointOnBoxedSurface( r_compr, &f, &u, &v, M, mlow, mhigh, &distance );

		double rpt[3],npt[3];
		evaluateRNRM( f, u, v, rpt, npt, r );

		double dr[3] = { rpt[0] - r_compr[0], rpt[1] - r_compr[1], rpt[2] - r_compr[2] };

		cute += (thePt->k/2) * (dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);	
	}

	return cute;	
}

double surface::cutGrad( double *r, double *g )
{
	double **M;
	int mlow,mhigh;
	getM(&M,&mlow,&mhigh);

	double cute = 0;

	for( fixed_cut_point *thePt = cutPoints; thePt; thePt = thePt->next )
	{
		double r_compr[3] = { thePt->rpt[0], thePt->rpt[1], thePt->rpt[2] };

		int f;
		double u, v;
		double distance;
		nearPointOnBoxedSurface( r_compr, &f, &u, &v, M, mlow, mhigh, &distance );

		double rpt[3],npt[3];

		evaluateRNRM( f, u, v, rpt, npt, r );

		double dr[3] = { rpt[0] - r_compr[0], rpt[1] - r_compr[1], rpt[2] - r_compr[2] };

		// dr d mesh points.

		int ncoords;	
		double coeffs[12+MAX_VALENCE];
		int coord_list[12+MAX_VALENCE];

		get_pt_coeffs( f, u, v,	coeffs, coord_list, &ncoords ); 	

		for( int t = 0; t < ncoords; t++ )
		{
			int c = coord_list[t];

			g[3*c+0] += (thePt->k) * dr[0] * coeffs[t];
			g[3*c+1] += (thePt->k) * dr[1] * coeffs[t];
			g[3*c+2] += (thePt->k) * dr[2] * coeffs[t];
		}

		cute += (thePt->k/2) * (dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);	
	}

	return cute;
}






