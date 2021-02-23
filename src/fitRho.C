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
#include "meshCollisionLib.h"
#include "io_mol_read.h"
#include "input.h"
#include "mutil.h"

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
//		printf("max strain %lf %lf %lf f %d u %lf v %lf\n", max_strain, max_c, max_k, max_f, max_u, max_v );
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
	int noiters=0;

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

		noiters++;
		if( noiters > 100 )
			found_pt = -1; 

	} while( found_pt == 0 );
	

	if( found_pt == 1 )
	{
		*f_in = f;
		*u_in = u;
		*v_in = v;
	}


	return found_pt;
}


void surface::mc_cut_worker( double value, int cartesian_component, double r_value, int *f_pts, double *uv_pts, double *rall, double *r, int n_cut_points, double *cen, int use_cen )
{
	int n_temp_stages = 1;
	int n_mc = 30 * n_cut_points;
	int mode = 0;
	double k_r = 1;

	if( r_value > 0 )
		mode = 1;

	if( mode == 1 )
		n_mc *= 100;

	double bad_tries = 0;

	double nacc = 0;
	double nrej = 0;

	double T = 0.1;

	for( int temp_stage = 0; temp_stage < n_temp_stages; temp_stage++ )
	{
		double ring_target = 2*M_PI * r_value / n_cut_points;
		ring_target /= (1 + 0.1 * (n_temp_stages-1-temp_stage)/(n_temp_stages) );
 
		for( int mc = 0; mc < n_mc; mc++ )
		{
			double cure = 0;
	
			if( mode == 0 )
			{
				for( int p = 0; p < n_cut_points; p++ )
				{
					double rp1[3] = { rall[3*p+0], rall[3*p+1], rall[3*p+2] };
					if( use_cen )
					{
						double dr[3];
						dr[0] = rp1[0] - cen[0];
						dr[1] = rp1[1] - cen[1];
						dr[2] = rp1[2] - cen[2];
		
						wrapPBC( dr );
		
						rp1[0] = cen[0] + dr[0];
						rp1[1] = cen[1] + dr[1];
						rp1[2] = cen[2] + dr[2];
					}
					for( int p2 = p+1; p2 < n_cut_points; p2++ )
					{
						double rp2[3] = { rall[3*p2+0], rall[3*p2+1], rall[3*p2+2] };
						if( use_cen )
						{
							double dr[3];
							dr[0] = rp2[0] - cen[0];
							dr[1] = rp2[1] - cen[1];
							dr[2] = rp2[2] - cen[2];
			
							wrapPBC( dr );
			
							rp2[0] = cen[0] + dr[0];
							rp2[1] = cen[1] + dr[1];
							rp2[2] = cen[2] + dr[2];
						}
						double dr[3] = { rp2[0]-rp1[0], rp2[1]-rp1[1], rp2[2]-rp1[2] };
			
						double l = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
			
						cure += 1.0/l/l;	
					}
				}
			}
			else if( mode == 1 )
			{
				for( int p = 0; p < n_cut_points; p++ )
				{
					double rp[3] = { rall[3*p+0], rall[3*p+1], rall[3*p+2] };
					if( use_cen )
					{
						double dr[3];
						dr[0] = rp[0] - cen[0];
						dr[1] = rp[1] - cen[1];
						dr[2] = rp[2] - cen[2];
	
						wrapPBC( dr );
	
						rp[0] = cen[0] + dr[0];
						rp[1] = cen[1] + dr[1];
						rp[2] = cen[2] + dr[2];
					}
					
					double dr[3] = { rp[0]-cen[0], rp[1]-cen[1], rp[2]-cen[2] };
	
					if( cartesian_component >= 0 && cartesian_component <= 2 )
						dr[cartesian_component] = dr[2];
	
					dr[2] = 0;
	
					double l = normalize(dr);
	
					cure += k_r * (l - r_value) * (l - r_value)/r_value/r_value; 
				}
				
				for( int p = 0; p < n_cut_points; p++ )
				{
					double rp1[3] = { rall[3*p+0], rall[3*p+1], rall[3*p+2] };
					if( use_cen )
					{
						double dr[3];
						dr[0] = rp1[0] - cen[0];
						dr[1] = rp1[1] - cen[1];
						dr[2] = rp1[2] - cen[2];
		
						wrapPBC( dr );
		
						rp1[0] = cen[0] + dr[0];
						rp1[1] = cen[1] + dr[1];
						rp1[2] = cen[2] + dr[2];
					}
					int p2 = p+1;
	
					if( p2 >= n_cut_points ) p2 -= n_cut_points;
	
					double rp2[3] = { rall[3*p2+0], rall[3*p2+1], rall[3*p2+2] };
	
					if( use_cen )
					{
						double dr[3];
						dr[0] = rp2[0] - cen[0];
						dr[1] = rp2[1] - cen[1];
						dr[2] = rp2[2] - cen[2];
			
						wrapPBC( dr );
			
						rp2[0] = cen[0] + dr[0];
						rp2[1] = cen[1] + dr[1];
						rp2[2] = cen[2] + dr[2];
					}
					double dr[3] = { rp2[0]-rp1[0], rp2[1]-rp1[1], rp2[2]-rp1[2] };
			
					double l = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
			
					double strain = (l-ring_target)/(ring_target);
					cure += (k_r/10) * strain*strain;				
				}
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
	
			
			int found_pt;
	
			if( mode == 0 )
				found_pt = movePointToCut(  &f, &u, &v, cartesian_component, value, r );
			
			if( found_pt == -1 )
			{
				bad_tries += 1;
			}
	
			if( mode == 1 || found_pt == 1 )
			{
	
				double rpt[3], nrm[3];
				evaluateRNRM( f, u, v, rpt, nrm, r );
	
				rall[3*to_move+0] = rpt[0]; 
				rall[3*to_move+1] = rpt[1]; 
				rall[3*to_move+2] = rpt[2]; 
				double newe = 0;
	
				if( mode == 0 )
				{
					for( int p = 0; p < n_cut_points; p++ )
					{
						double rp1[3] = { rall[3*p+0], rall[3*p+1], rall[3*p+2] };
						if( use_cen )
						{
							double dr[3];
							dr[0] = rp1[0] - cen[0];
							dr[1] = rp1[1] - cen[1];
							dr[2] = rp1[2] - cen[2];
		
							wrapPBC( dr );
		
							rp1[0] = cen[0] + dr[0];
							rp1[1] = cen[1] + dr[1];
							rp1[2] = cen[2] + dr[2];
						}
						for( int p2 = p+1; p2 < n_cut_points; p2++ )
						{
							double rp2[3] = { rall[3*p2+0], rall[3*p2+1], rall[3*p2+2] };
							if( use_cen )
							{
								double dr[3];
								dr[0] = rp2[0] - cen[0];
								dr[1] = rp2[1] - cen[1];
								dr[2] = rp2[2] - cen[2];
			
								wrapPBC( dr );
			
								rp2[0] = cen[0] + dr[0];
								rp2[1] = cen[1] + dr[1];
								rp2[2] = cen[2] + dr[2];
							}
							double dr[3] = { rp2[0]-rp1[0], rp2[1]-rp1[1], rp2[2]-rp1[2] };
			
							double l = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
			
							newe += 1.0/l/l;	
						}
					}
				}
				else if( mode == 1 )
				{
					for( int p = 0; p < n_cut_points; p++ )
					{
						double rp[3] = { rall[3*p+0], rall[3*p+1], rall[3*p+2] };
						if( use_cen )
						{
							double dr[3];
							dr[0] = rp[0] - cen[0];
							dr[1] = rp[1] - cen[1];
							dr[2] = rp[2] - cen[2];
	
							wrapPBC( dr );
	
							rp[0] = cen[0] + dr[0];
							rp[1] = cen[1] + dr[1];
							rp[2] = cen[2] + dr[2];
						}
						
						double dr[3] = { rp[0]-cen[0], rp[1]-cen[1], rp[2]-cen[2] };
	
						if( cartesian_component >= 0 && cartesian_component <= 2 )
							dr[cartesian_component] = dr[2];
	
						dr[2] = 0;
	
						double l = normalize(dr);
	
						newe += k_r * (l - r_value) * (l - r_value)/r_value/r_value; 
					}
					for( int p = 0; p < n_cut_points; p++ )
					{
						double rp1[3] = { rall[3*p+0], rall[3*p+1], rall[3*p+2] };
						if( use_cen )
						{
							double dr[3];
							dr[0] = rp1[0] - cen[0];
							dr[1] = rp1[1] - cen[1];
							dr[2] = rp1[2] - cen[2];
			
							wrapPBC( dr );
			
							rp1[0] = cen[0] + dr[0];
							rp1[1] = cen[1] + dr[1];
							rp1[2] = cen[2] + dr[2];
						}
						int p2 = p+1;
		
						if( p2 >= n_cut_points ) p2 -= n_cut_points;
		
						double rp2[3] = { rall[3*p2+0], rall[3*p2+1], rall[3*p2+2] };
		
						if( use_cen )
						{
							double dr[3];
							dr[0] = rp2[0] - cen[0];
							dr[1] = rp2[1] - cen[1];
							dr[2] = rp2[2] - cen[2];
				
							wrapPBC( dr );
				
							rp2[0] = cen[0] + dr[0];
							rp2[1] = cen[1] + dr[1];
							rp2[2] = cen[2] + dr[2];
						}
						double dr[3] = { rp2[0]-rp1[0], rp2[1]-rp1[1], rp2[2]-rp1[2] };
				
						double l = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
				
						double strain = (l-ring_target)/(ring_target);
						newe += (k_r/10) * strain*strain;				
					}
				}
			
				double pr;

//				if( mode == 1 ) printf("T: %lf dE: %lf\n", T, newe-cure );
				if( temp_stage == n_temp_stages-1 )
				{
					pr =0;
					if( newe < cure ) pr = 1;
				}
				else
					pr = exp( - (newe-cure)/T );

				double p = rand()/(double)RAND_MAX;

				if( p < pr  )
				{
					nacc++;
					f_pts[to_move] = f;
					uv_pts[2*to_move] = u;
					uv_pts[2*to_move+1] = v;
	
					cure = newe;
				}
				else
				{
					nrej++;
					rall[3*to_move+0] = rsave[0];
					rall[3*to_move+1] = rsave[1];
					rall[3*to_move+2] = rsave[2];
				}
			}
	
			if( mc > 20 && bad_tries / (double)mc > 0.5 )
				break;
		}

		T *= 0.8;
	}
}

void surface::spread_evenly_at_cut( double value, int cartesian_component, int *f_pts, double *uv_pts, double *rall, double *r, int np, double *cen, int use_cen )
{
	mc_cut_worker( value, cartesian_component, -1, f_pts, uv_pts, rall, r, np, cen, use_cen );
}

void surface::fit_to_circle( int cartesian_component, double r_off, int *f_pts, double *uv_pts, double *rall, double *r, int np, double *cen, int use_cen )
{
	mc_cut_worker( 0, cartesian_component, r_off, f_pts, uv_pts, rall, r, np, cen, use_cen);
}



int surface::get_cut_points( int cartesian_component, double value, int *f_pts, double *uv_pts, double *rall, int n_cut_points, double *r, int *nconvex_out, double cen[3], int use_cen, int reduce_convex )
{
	
	for( int p = 0; p < n_cut_points; p++ )
	{
		int done = 0;
		int ncyc = 0;

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
		
			if( ncyc > 1000 ) 
				return 0;

			ncyc++;
		}
	}		

	spread_evenly_at_cut( value, cartesian_component, f_pts, uv_pts, rall, r, n_cut_points, cen, use_cen );

#if 0
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
#endif
	
	if( reduce_convex )	
	{
		void giftwrap( double *pre_pts_in, int *ptsOrdered, int npts, int *nconvex, int expand );
	
		int ptsOrdered[n_cut_points];
		double r_xy[2*n_cut_points];
		int nconvex = 0;
		
		for( int i = 0; i < n_cut_points; i++ )
		{
			double t[3] = { rall[3*i+0], rall[3*i+1], rall[3*i+2] };
	
			if( use_cen )
			{
				double dr[3] = { t[0] - cen[0], t[1] - cen[1], t[2] - cen[2] };
	
				wrapPBC( dr );
	
				t[0] = cen[0] + dr[0];	
				t[1] = cen[1] + dr[1];	
				t[2] = cen[2] + dr[2];	
			}
	
			// this just replaces x/y with z if we are looking normal to that component.
			if( cartesian_component >=0 && cartesian_component <= 2 )
				t[cartesian_component] = t[2];
		
			r_xy[2*i+0] = t[0];
			r_xy[2*i+1] = t[1];
		}
	
		giftwrap( r_xy, ptsOrdered, n_cut_points, &nconvex, 0 );
	
		int f_out[nconvex];
		double uv_out[nconvex*2];
		double rallout[nconvex*3];
	
		for( int t = 0; t < nconvex; t++ )
		{
			f_out[t] = f_pts[ptsOrdered[t]];
			uv_out[2*t+0] = uv_pts[ptsOrdered[t]*2+0];	
			uv_out[2*t+1] = uv_pts[ptsOrdered[t]*2+1];	
	
			rallout[3*t+0] = rall[3*ptsOrdered[t]+0];
			rallout[3*t+1] = rall[3*ptsOrdered[t]+1];
			rallout[3*t+2] = rall[3*ptsOrdered[t]+2];
		}	
	
		memcpy( f_pts, f_out, sizeof(int) * nconvex );	
		memcpy( uv_pts, uv_out, sizeof(double) * nconvex*2 );	
		memcpy( rall, rallout, sizeof(double) * nconvex*3 );	
	
	
		*nconvex_out = nconvex;
	}
	else *nconvex_out = n_cut_points; 
	return 1;
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

int surface::setupCut( int cartesian_component, double value, double *r)
{
	// here's how we'll get them. find n_cut_points random points on the surface at the cut value
	// use a pseudopotential to optimize their position at the cut value to try to make them evenly spaced.		
	int *f_pts = (int *)malloc( sizeof(int) * n_cut_points );
	double *uv_pts = (double *)malloc( sizeof(double) * 2 * n_cut_points );	
	double *rall = (double *)malloc( sizeof(double) * 3 * n_cut_points );

	int nconvex;
	double cen[3]={0,0,0};
	int ok = get_cut_points( cartesian_component, value, f_pts, uv_pts, rall, n_cut_points, r, &nconvex, cen, 0, 0 );

	if( !ok ) { n_cut_points = 0; return 0; }

	n_cut_points = nconvex;

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
		thePt->frc[0] = 0;
		thePt->frc[1] = 0;
		thePt->frc[2] = 0;
	}
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

			thePt->frc[0] -= (thePt->k) * dr[0] * coeffs[t];	
			thePt->frc[1] -= (thePt->k) * dr[1] * coeffs[t];	
			thePt->frc[2] -= (thePt->k) * dr[2] * coeffs[t];	
		}

		cute += (thePt->k/2) * (dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);	
	}

	return cute;
}

double ring_perimeter( surface *theSurface, double trial_z, double *rsurf, double *cen, int use_cen)
{
	int n_cut_points = 20;
	int *f_pts = (int *)malloc( sizeof(int) * n_cut_points );
	double *uv_pts = (double *)malloc( sizeof(double) * 2 * n_cut_points );
	double *rall = (double *)malloc( sizeof(double) * 3 * n_cut_points );
	
	int nconvex;

	int ok = theSurface->get_cut_points( 2, trial_z, f_pts, uv_pts, rall, n_cut_points, rsurf, &nconvex, cen, use_cen, 1 );
	if( !ok ) return -1;

	double perimeter = 0;

	for( int i = 0; i < nconvex; i++ )
	{
		int ip1 = i+1;
		while( ip1 >= nconvex )
			ip1 -= nconvex;

		double dr[3] = { rall[3*i+0] - rall[3*ip1+0],
				 rall[3*i+1] - rall[3*ip1+1],
				 rall[3*i+2] - rall[3*ip1+2] };
		theSurface->wrapPBC(dr);

		//printf("pt %lf %lf %lf\n", rall[3*i+0], rall[3*i+1], rall[3*i+2] );

		perimeter += normalize(dr);
	}

	free(f_pts);
	free(uv_pts);
	free(rall);

	return perimeter;
}



void surface::GetFusionPoreRegionStats( double *rsurf, parameterBlock *block )
{
	put( rsurf );

	int default_dz = 0;
	double outer_R = block->pore_outer_cut;
	double pore_dz_from_center = block->pore_dz_from_center;

	if( pore_dz_from_center < 0 )
		default_dz = 1;

	double zm[3] = {-25, 0, 25};	
	double pm[3];
	double cen[3] = {0,0,0};

	for( int t= 0; t < 3; t++ )
	{
		pm[t] = ring_perimeter(this, zm[t], rsurf, cen, 0 ); 
		printf("Perimeter %d %lf\n", t, pm[t] );
	}

	if( pm[0] == -1 && pm[1] == -1 && pm[2] == -1 )
	{
		printf("Can't find the pore center.\n");
		exit(1);
	}
	
	if( pm[0] == -1 || pm[1] == -1 || pm[2] == -1 )
	{
		printf("Can't find the pore center, the code needs to be improved.\n");
		exit(1);
	}


	int nop = 20;

	for( int t = 0; t < nop; t++ )
	{
		double trial_low = (zm[0]+zm[1])/2;

		double p = ring_perimeter(this, trial_low, rsurf, cen, 0 ); 
	
		if( p < pm[1] )
		{
//			printf("Low: low value: %lf at %lf\n", p, trial_low );
			zm[2] = zm[1];
			pm[2] = pm[1];
			zm[1] = trial_low;
			pm[1] = p;	
		}
		else
		{
//			printf("Low: higher value: %lf at %lf\n", p, trial_low );
			zm[0] = trial_low;
			pm[0] = p;
		}

		double trial_high = (zm[1]+zm[2])/2;

		p = ring_perimeter(this, trial_high, rsurf, cen, 0 ); 
	
		if( p < pm[1] )
		{
//			printf("High: low value: %lf at %lf\n", p, trial_high );
			zm[0] = zm[1];
			pm[0] = pm[1];
			zm[1] = trial_high;
			pm[1] = p;	
		}
		else
		{
//			printf("High: high value: %lf at %lf\n", p, trial_high );
			zm[2] = trial_high;
			pm[2] = p;	
		}

//		printf(" %lf %lf %lf <-> %lf %lf %lf\n", zm[0], zm[1], zm[2], pm[0], pm[1], pm[2] );

	} 

	int n_cut_points = 20;
	int *f_pts = (int *)malloc( sizeof(int) * n_cut_points );
	double *uv_pts = (double *)malloc( sizeof(double) * 2 * n_cut_points );
	double *rall = (double *)malloc( sizeof(double) * 3 * n_cut_points );
	int nconvex = 0;	
	int ok = get_cut_points( 2, zm[1], f_pts, uv_pts, rall, n_cut_points, rsurf, &nconvex, cen, 0, 1 );

	double c_sign_factor = 1;


	// this defines a set of points in the interior.
	// we now want to define circular cuts with other criteria.
	
	double pore_center[3] = {0,0,0};
	double w = 0;

	for( int i = 0; i < nconvex; i++ )
	{
		int ip1 = i+1;
		while( ip1 >= nconvex ) ip1 -= nconvex;
		double dr[3] = { rall[3*ip1+0] - rall[3*i+0],
				 rall[3*ip1+1] - rall[3*i+1],
				 rall[3*ip1+2] - rall[3*i+2] }; 
		double l = normalize(dr);
		
		pore_center[0] += rall[3*i+0] * l/2; 
		pore_center[1] += rall[3*i+1] * l/2; 
		pore_center[2] += rall[3*i+2] * l/2; 
		
		pore_center[0] += rall[3*ip1+0] * l/2; 
		pore_center[1] += rall[3*ip1+1] * l/2; 
		pore_center[2] += rall[3*ip1+2] * l/2; 

		w += l;
	}

	pore_center[0] /= w;
	pore_center[1] /= w;
	pore_center[2] /= w;

	printf("Pore center: %lf %lf %lf\n", pore_center[0], pore_center[1], zm[1] );
	
	double av_c_pore = 0;
	double av_c_pore_perp = 0;
	double av_dp = 0;
	for( int p = 0; p < nconvex; p++ )
	{
		int f = f_pts[p];
		double u = uv_pts[2*p+0];
		double v = uv_pts[2*p+1];

		double rpt[3],npt[3];	
		evaluateRNRM( f, u, v, rpt, npt, rsurf );

		double K;
		double c_vec_1[2], c_vec_2[2], c_val_1, c_val_2;
		c( f, u, v, rsurf, &K, c_vec_1, c_vec_2, & c_val_1, & c_val_2 );	 
		double drdu[3],drdv[3];
		ru( f, u, v, rsurf, drdu );
		rv( f, u, v, rsurf, drdv );

		double dr[3] = { pore_center[0] - rpt[0], pore_center[1] - rpt[1], pore_center[2] - rpt[2] };
		wrapPBC(dr);

		normalize(dr);
		double dp = dr[0]*npt[0] + dr[1]*npt[1]+dr[2]*npt[2];

		av_dp += dp;

		double r_c_vec_1[3] = { drdu[0] * c_vec_1[0] + drdv[0] * c_vec_1[1],
					drdu[1] * c_vec_1[0] + drdv[1] * c_vec_1[1],
					drdu[2] * c_vec_1[0] + drdv[2] * c_vec_1[1] };
		double r_c_vec_2[3] = { drdu[0] * c_vec_2[0] + drdv[0] * c_vec_2[1],
					drdu[1] * c_vec_2[0] + drdv[1] * c_vec_2[1],
					drdu[2] * c_vec_2[0] + drdv[2] * c_vec_2[1] };

		normalize(r_c_vec_1);
		normalize(r_c_vec_2);

		if( fabs(r_c_vec_1[2]) > fabs(r_c_vec_2[2]) )
		{
			av_c_pore_perp += c_val_1;
			av_c_pore += c_val_2;
		}	
		else
		{
			av_c_pore_perp += c_val_2;
			av_c_pore += c_val_1;
		}	
	}

	av_dp /= nconvex;
	av_c_pore_perp /= nconvex;
	av_c_pore /= nconvex;

	if( av_dp < 0 )
	{
		printf("Switching c sign (%le)\n", av_dp );
		c_sign_factor *= -1;	
	}
	av_c_pore_perp *= c_sign_factor;
	av_c_pore *= c_sign_factor;

	printf("Pore curvature %le perp %le\n", av_c_pore, av_c_pore_perp );
	// translate these point directly up in z, fix that z, then minimize the perimeter?

	// angstroms, spacing
	double dz = 1;
	// numer of discrete cuts we look at
	int scan_num = 50;	

	// this scans the fusion pore from min to max.
	double min_layer = -scan_num*dz;
	double max_layer = scan_num*dz;


	if( outer_R < 0 )
		outer_R = PBC_vec[0][0]/2 * 0.9;

	int best_lower_cut = -1;
	int best_upper_cut = -1;
	double best_lower_chi2 = 1e10;
	double best_upper_chi2 = 1e10;

	// loop over the cuts.
	for( int iz = -scan_num*dz; iz <= scan_num*dz; iz++ )
	{
		// this is the z value of the cut.
		double z_cut = zm[1] + dz * iz;

		// the perimeter at this cut, implies the radius.
		double p = ring_perimeter(this, z_cut, rsurf, pore_center, 1 ); 

		if( p < 0 ) // p < 0 means the cut was beyond our pore.
		{
			if( iz < 0 )
			{
				if( z_cut + dz > min_layer )
					min_layer = z_cut + dz;
			}
			else 
			{
				if( z_cut - dz < max_layer )
					max_layer = z_cut - dz;
			}
			// this is outside what we'd care about.
		}	
		else
		{
			// if the radius matches where we establish our boundary, start there.
			double R = p / (2*M_PI);
			double chi2 = (R-outer_R)*(R-outer_R);

			if( iz < 0 )
			{
				if( chi2 < best_lower_chi2 )
				{
					best_lower_chi2 = chi2;
					best_lower_cut = z_cut;	
				}
			}
			else if( chi2 < best_upper_chi2 )
			{
				best_upper_chi2 = chi2;
				best_upper_cut = z_cut;	
			}
		}

		printf("z: %lf perimeter: %lf\n", z_cut, p );
		//get_cut_points( 2, z_cut, f_pts, uv_pts, rall, n_cut_points, rsurf, &nconvex );						
	}

	if( default_dz )
	{
		pore_dz_from_center = (max_layer-min_layer)/2 - 20.0;
		printf("Using default_dz calc: %le\n", pore_dz_from_center );
	} 

	int *face_labels = (int *)malloc( sizeof(int) * nt );

	for( int f = 0; f < nt; f++ )
	{
		double rp[3],np[3];

		evaluateRNRM( f, 1.0/3.0, 1.0/3.0, rp, np, rsurf );

		double dr_cen[3] = { rp[0] - pore_center[0], rp[1] - pore_center[1], rp[2] - pore_center[2] };

		wrapPBC(dr_cen);

		double rxy = sqrt(dr_cen[0]*dr_cen[0]+dr_cen[1]*dr_cen[1]);

		if( fabs(dr_cen[2]) < pore_dz_from_center )
			face_labels[f] = 2; // inner pore
		else if( rxy < outer_R )
			face_labels[f] = 1;
		else
			face_labels[f] = 0; 
	}

	typedef struct lipid_data
	{
		int num;
		char *resname;
		struct lipid_data *next;
	} lipid_data;

	typedef struct
	{
		double r[3];

		lipid_data *theDataInside;
		lipid_data *theDataOutside;
		double sumJ,sumK;
		double sumJ2,sumK2;
		double area;
		double av_psm_orientation;
		double nav_psm;
	} face_data;	

	face_data *theData = (face_data *)malloc( sizeof(face_data) * nt );

	for( int t = 0; t < nt; t++ )
	{
		double nrm[3];
		evaluateRNRM( t, 1.0/3.0 ,1.0/3.0, theData[t].r, nrm, rsurf ); 

		theData[t].theDataInside = NULL;
		theData[t].theDataOutside = NULL;
		theData[t].av_psm_orientation = 0;
		theData[t].nav_psm = 0;
		theData[t].sumJ=0;
		theData[t].sumJ2=0;
		theData[t].sumK=0;
		theData[t].sumK2=0;
	}

	int plim = 5;

	int ntri_set = -1;

	for( int f = 0; f < nt; f++ )
	{
		int ntri = 0;
		for( int iu = 0; iu < plim; iu++ )
		{
			for( int iv = 0; iv <= iu; iv++ )
			{
				double fu = (iu) / plim; 
				double fup1 = (iu+1) / plim; 
				double fv = (iv) / plim; 
				double fvp1 = (iv+1) / plim; 

				double sqrs[2][3][2] = { { { iu, iv}, {iu+1,iv}, { iu, iv+1} },
						         { { iu+1, iv+1}, {iu+1,iv}, { iu, iv+1} } };
				for( int sqr = 0; sqr < 2; sqr++ )
				{
					int ok = 1;
					double cen[2] = { 0,0};
					for( int p =0; p < 3; p++ )
					{
						double fu = sqrs[sqr][p][0] / (double)plim;
						double fv = sqrs[sqr][p][1] / (double)plim;


						cen[0] += fu/3;
						cen[1] += fv/3;
					}

					if( cen[0] > 1 || cen[1] > 1 || cen[0]+cen[1] > 1 )
						continue;

					ntri += 1;

					double K;
					double c_vec_1[3], c_vec_2[3], c_val_1, c_val_2;
					c( f, cen[0], cen[1], rsurf, &K, c_vec_1, c_vec_2, & c_val_1, & c_val_2 );
	
					double J = c_sign_factor * c_val_1 + c_sign_factor * c_val_2;
					double l_A = g( f, cen[0], cen[1], rsurf )/2;

					theData[f].sumJ += J * l_A;
					theData[f].sumJ2 += J*J * l_A;
					theData[f].sumK += K * l_A;
					theData[f].sumK2 += K*K * l_A;
					theData[f].area += l_A;	
	 
				}
			}
		}		

		if( ntri_set < 0 ) ntri_set = ntri;
	}




	double av_c = 0;
	double av_k = 0;
	double sum_g = 0;	

	for( int f = 0; f < nt; f++ )
	{
		double k;
		double gv =  g(f, 1.0/3.0, 1.0/3.0, rsurf ); 
		av_c += gv * c(f, 1.0/3.0, 1.0/3.0, rsurf, &k ); 
		av_k += k * gv;
		sum_g += gv;
	}

	av_c /= sum_g;
	av_k /= sum_g;

	double **M;
	int mlow,mhigh;

	getM( &M, &mlow, &mhigh );

        char *gathered = (char *)malloc( sizeof(char) * (1+strlen(block->jobName) + strlen("_gathered.xyz")) );
	sprintf(gathered, "%s_gathered.xyz", block->jobName );
	FILE *gatheredXYZ = fopen(gathered,"w");
	if( gatheredXYZ )		
	{
		writeLimitingSurface(gatheredXYZ);
		fclose(gatheredXYZ);
	}
	if( !block->structureName )
	{
		printf("Gathering requires a structure (e.g., PSF) file.\n");
		exit(1);
	}

	io_readStructureFile( block->structureName );	
	io_initialize_read( block->dcdName );

	const char *append_string = "_fusion_pore_gather.txt";
	char *gatherName = (char *)malloc( sizeof(char) * ( 1 + strlen(block->jobName) + strlen(append_string) ) );

	sprintf(gatherName, "%s%s", block->jobName, append_string );

	FILE *gatherFile = fopen(gatherName, "w");

	if( !gatherFile )
	{
		printf("Couldn't open gather file '%s' for writing.\n", gatherName );
		exit(1);
	}

	printf("Starting gather.\n");

	char assignmentName[256];

	sprintf(assignmentName, "%s_assignment.txt", block->jobName );
	FILE*assignmentFile = fopen(assignmentName, "w");

	int nframes = io_nframes();

	double MIN_K=-1e-3;
	double MAX_K=1e-3;

	int N_K_BINS = 20;
	int N_ORDER_BINS = 20;
	
	struct order_record
	{	
		char *resName;
		double *order_p_dist;
		struct order_record *next;
	};

	struct order_record *allRecords = NULL;
	
	for( int f = 0; f < nframes; f++ )
	{
		io_readFrame( block->dcdName );
		printf("Gathering from frame %d.\n", f );
		fflush(stdout);
		int ns = io_nNSAtoms();
		io_align();
 
		struct atom_rec *at = io_getCurrentFrame();

		for( int ax = 0; ax < ns; ax ++ )
		{
			int a = io_getNSAtom(ax);
			int astart=-1;
			int astop=-1;
			io_getNSAtomStartStop(ax,&astart,&astop); 
			double r[3] = { at[a].x, at[a].y, at[a].z };
			double cout,kout,dp_out,dz_out;
			int leaflet_out;
			double nrm[3];
			int f;
			double u,v;
			double distance;	
			nearPointOnBoxedSurface( r, &f, &u, &v, M, mlow, mhigh, &distance );					

			double rp[3];
			evaluateRNRM( f, u, v, rp, nrm, rsurf );

			// get nrm from polarity.

			int new_leaflet_out = -1;		

			double psm_chain_1[3]={0,0,0};
			double psm_chain_2[3]={0,0,0};
			double npsm1=0,npsm2=0;
			
			double pc_chain_1[3]={0,0,0};
			double pc_chain_2[3]={0,0,0};
			double npc1=0,npc2=0;

			double chain_1[3] = {0,0,0};
			double chain_2[3] = {0,0,0};
			int nc1=0;
			int nc2=0;

			if( astart >= 0 && astop >= 0 )
			{
				struct order_record *theRec = NULL;

				for( struct order_record *tRec = allRecords; tRec; tRec = tRec->next )
				{
					if( !strcasecmp( tRec->resName, at[astart].resname) ) 
						theRec = tRec;
				}

				if( !theRec )
				{
					theRec = (struct order_record *)malloc( sizeof(order_record) );
					theRec->order_p_dist = (double *)malloc( sizeof(double) * N_K_BINS * N_ORDER_BINS );
					memset( theRec->order_p_dist, 0, sizeof(double) * N_K_BINS * N_ORDER_BINS ); 
					theRec->next = allRecords;
					allRecords = theRec;
					theRec->resName = (char *)malloc( sizeof(char)*(1+strlen(at[astart].resname)) );
					strcpy( theRec->resName, at[astart].resname );
				}


				for( int tx = astart; tx < astop; tx++ )
				{
					int wc = 0;

					if( !strcasecmp( at[tx].resname, "PSM" ) || !strncasecmp( at[tx].resname, "CER1", 4 ) )
					{
						if( at[tx].atname[0] == 'C' && at[tx].atname[strlen(at[tx].atname)-1] == 'F' )
							wc = 1;
						if( at[tx].atname[0] == 'C' && at[tx].atname[strlen(at[tx].atname)-1] == 'S' )
							wc = 2;
					}
					else if( !strcasecmp( at[tx].resname, "CHL1" ) )
					{
						if( !strcasecmp( at[tx].atname, "C4" ) || !strcasecmp(at[tx].atname,"C6") ||
						    !strcasecmp( at[tx].atname, "C7" ) || !strcasecmp(at[tx].atname,"C15") || 
						    !strcasecmp( at[tx].atname, "C16" ) )
							wc = 1;
						if( !strcasecmp( at[tx].atname, "C2" ) || !strcasecmp(at[tx].atname,"C1") ||
						    !strcasecmp( at[tx].atname, "C11" ) || !strcasecmp(at[tx].atname,"C12") )
							wc = 2; 
					}
					else
					{
						if( at[tx].atname[0] == 'C' && at[tx].atname[1] == '2' )
							wc = 1;
						if( at[tx].atname[0] == 'C' && at[tx].atname[1] == '3' )
							wc = 2;
					}
			
					if( wc == 0 ) continue;

					if( wc == 1 )
					{
						chain_1[0] += at[tx].x;
						chain_1[1] += at[tx].y;
						chain_1[2] += at[tx].z;
						nc1++;
					}
					else	
					{
						chain_2[0] += at[tx].x;
						chain_2[1] += at[tx].y;
						chain_2[2] += at[tx].z;
						nc2++;
					}	
				}

				chain_1[0] /= nc1;
				chain_1[1] /= nc1;
				chain_1[2] /= nc1;
				
				chain_2[0] /= nc2;
				chain_2[1] /= nc2;
				chain_2[2] /= nc2;


				double polarity = 0;
	
				double AVP = 0;
				for( int tx = astart; tx <= astop; tx++ )
				{
					double PM = -1;
					if( at[tx].atname[0] == 'C' || at[tx].atname[0] == 'H' )
						PM = 1;
					AVP += PM;
				}
				
				AVP /= (astop-astart+1);
				for( int tx = astart; tx <= astop; tx++ )
				{
					double PM = -1;
	
					if( at[tx].atname[0] == 'C' || at[tx].atname[0] == 'H' )
						PM = 1;
					polarity += (PM-AVP) * (at[tx].x * nrm[0] + at[tx].y * nrm[1] + at[tx].z * nrm[2]);
				}

				if( polarity < 0 ) 
					new_leaflet_out = 1;
				else
					new_leaflet_out = 0;

			//	if( !strcasecmp( at[a].resname, "PSM" ) || !strcasecmp( at[a].resname, "PLPC" ) )
				if( nc1 > 0 && nc2 > 0 )
				{
					double lipid_dr[3];

					int do_psm = 0;
					if( !strcasecmp( at[a].resname, "PSM" ) ) 
						do_psm = 1;
					lipid_dr[0] = chain_2[0] - chain_1[0];
					lipid_dr[1] = chain_2[1] - chain_1[1];
					lipid_dr[2] = chain_2[2] - chain_1[2];

					double dp = lipid_dr[0] * nrm[0] + lipid_dr[1] * nrm[1] + lipid_dr[2] * nrm[2];
					lipid_dr[0] -= dp * nrm[0];
					lipid_dr[1] -= dp * nrm[1];
					lipid_dr[2] -= dp * nrm[2];
					normalize(lipid_dr);

					double K;
					double c_vec_1[2], c_vec_2[2], c_val_1, c_val_2;
					c( f, u, v, rsurf, &K, c_vec_1, c_vec_2, & c_val_1, & c_val_2 );	 
					double drdu[3],drdv[3];
					ru( f, u, v, rsurf, drdu );
					rv( f, u, v, rsurf, drdv );

					double r_c_vec_1[3] = { drdu[0] * c_vec_1[0] + drdv[0] * c_vec_1[1],
								drdu[1] * c_vec_1[0] + drdv[1] * c_vec_1[1],
								drdu[2] * c_vec_1[0] + drdv[2] * c_vec_1[1] };
					double r_c_vec_2[3] = { drdu[0] * c_vec_2[0] + drdv[0] * c_vec_2[1],
								drdu[1] * c_vec_2[0] + drdv[1] * c_vec_2[1],
								drdu[2] * c_vec_2[0] + drdv[2] * c_vec_2[1] };
			
					normalize(r_c_vec_1);
					normalize(r_c_vec_2);

					c_val_1 *= c_sign_factor;	
					c_val_2 *= c_sign_factor;	

					if( new_leaflet_out == 1 )
					{
						c_val_1 *= -1;
						c_val_2 *= -1;
					}

					double p_c_dp=0;
					if( c_val_2 > c_val_1 )
						p_c_dp = r_c_vec_2[0] * lipid_dr[0] + r_c_vec_2[1] * lipid_dr[1] + r_c_vec_2[2] * lipid_dr[2];
					else
						p_c_dp = r_c_vec_1[0] * lipid_dr[0] + r_c_vec_1[1] * lipid_dr[1] + r_c_vec_1[2] * lipid_dr[2];


					int kbin = N_K_BINS*(c_val_1*c_val_2-MIN_K)/(MAX_K-MIN_K);
					if( kbin < 0 ) kbin = 0;
					if( kbin >= N_K_BINS ) kbin = N_K_BINS-1;

					int ob = p_c_dp*p_c_dp * N_ORDER_BINS;
					int ab = N_ORDER_BINS * acos(p_c_dp)/M_PI;
					if( ab < 0 ) ab = 0; if( ab >= N_ORDER_BINS) ab = N_ORDER_BINS-1;

					theRec->order_p_dist[kbin*N_ORDER_BINS+ab] += 1;

					if( do_psm )
					{
						theData[f].av_psm_orientation += p_c_dp*p_c_dp;				
						theData[f].nav_psm += 1;				
					}
				}

			}

			lipid_data *got= NULL;

			lipid_data **theList = &theData[f].theDataInside;
		
			if( new_leaflet_out == 1 )
				theList = &theData[f].theDataOutside;
	

			for( lipid_data *aLipid = (*theList); aLipid; aLipid = aLipid->next )
			{
				if( !strcasecmp( aLipid->resname, at[astart].resname) )
					got = aLipid;
			}

			if( !got )
			{
				got = (lipid_data *)malloc( sizeof(lipid_data) );
				got->resname = (char *)malloc( sizeof(char) * (1+strlen(at[astart].resname) ) );
				strcpy( got->resname, at[astart].resname );
				got->num=0;
				got->next = *theList;
				*theList = got;
			}
	
			got->num += 1;			
		}	

	}

	fprintf(gatherFile, "min_layer %lf\n", min_layer );
	fprintf(gatherFile, "max_layer %lf\n", max_layer );
	fprintf(gatherFile, "center %lf %lf %lf\n", pore_center[0], pore_center[1], pore_center[2] );
	fprintf(gatherFile, "nframes %d\n", nframes );	
	fprintf(gatherFile, "pore_curvature %le\n", av_c_pore );
	fprintf(gatherFile, "perp_pore_curvature %le\n", av_c_pore_perp );


	for( int f = 0; f < nt; f++ )
	{
		fprintf(gatherFile, "%d J: %le J2: %le K: %le K2: %le area: %le r: %le %le %le", f, theData[f].sumJ/ntri_set, theData[f].sumJ2/ntri_set, theData[f].sumK/ntri_set, theData[f].sumK2/ntri_set, theData[f].area/ntri_set, theData[f].r[0], theData[f].r[1], theData[f].r[2] );
		fprintf(gatherFile, " inside");
		for( lipid_data *aLipid = theData[f].theDataInside; aLipid; aLipid = aLipid->next )
			fprintf(gatherFile, " %s:%d", aLipid->resname, aLipid->num);
		fprintf(gatherFile, " outside");
		for( lipid_data *aLipid = theData[f].theDataOutside; aLipid; aLipid = aLipid->next )
			fprintf(gatherFile, " %s:%d", aLipid->resname, aLipid->num);
		fprintf(gatherFile, " psm_orientation");
			fprintf(gatherFile, " %lf:%lf", theData[f].av_psm_orientation, theData[f].nav_psm );
		fprintf(gatherFile,"\n");
	}
		
	printf("Finished with gather.\n" );
	fflush(stdout);

	char fileNamePD[256];
	sprintf(fileNamePD, "%s_pdist.txt", block->jobName );
	FILE *orderFile = fopen(fileNamePD, "w");

	for( struct order_record *theRec = allRecords; theRec; theRec = theRec->next )
	{
		fprintf(orderFile, "residue %s\n", theRec->resName );

		double *order_p_dist = theRec->order_p_dist;
		for( int k = 0; k < N_K_BINS; k++ )
		{
			for( int o = 0; o < N_ORDER_BINS; o++ )
			{
				fprintf(orderFile, "%lf %lf %lf\n", MIN_K + (k+0.5)*(MAX_K-MIN_K)/N_K_BINS, M_PI*(o+0.5)/N_ORDER_BINS, order_p_dist[k*N_ORDER_BINS+o] );
			}
			fprintf(orderFile,"\n");
		}
	}
	fclose(orderFile);
}



