#include "simulation.h"
#include "interp.h"
#include "pcomplex.h"
#include <string.h>
#include <math.h>
#include "l-bfgs.h"
#include "p_p.h"
#include "parallel.h"
#include "fitRho.h"
#include "globals.h"
#include "clathrin.h"


extern int enable_elastic_interior;
static int min_freeze_clathrin = 0;
static int min_z_only = 0;
static int min_ncomplex = 0;
static int min_nparams = 0;
static int min_nsurfaceparams = 0;
static int do_freeze_membrane = 0;
static pcomplex **min_complexes;
static int fix_thickness = 0;
extern double VA,VC;
double VV = 0;
extern double water_KV;
Simulation *min_simulation = NULL;
static surface_record one_srec;
static double cur_rho_thickness[2] = { 15.0, 15.0 };

#define FIT_RHO_ONE_THICKNESS

double surface_f( double *p )
{
	int offset = 3; 
	
	surface_record *header = NULL;

	int do_sim = 1;

	if( min_simulation )
		header = min_simulation->allSurfaces;
	else
	{
		do_sim = 0;
		header = &one_srec;
	}
	for( surface_record *sRec = header; sRec; sRec = sRec->next )
	{
		int nv = sRec->theSurface->nv;
		p[offset+3*nv+0] = p[0];
		p[offset+3*nv+1] = p[1];
		p[offset+3*nv+2] = p[2];
		
		sRec->r = p+offset;		

//		memcpy( sRec->r, p+offset, sizeof(double) * (3*nv+3) );
		offset += sRec->theSurface->nv*3+3;
	}

#ifdef PARALLEL
	ParallelSyncComplexes( min_complexes, min_ncomplex );
#endif

	VV = 0;
	VA = 0;
	VC = 0;
	double v = 0;
	for( surface_record *sRec = header; sRec; sRec = sRec->next )
	{
		v += sRec->theSurface->energy( p+sRec->temp_min_offset, NULL );
	}
	ParallelSum(&VA,1);
	ParallelSum(&VC,1);
	int nparams = min_nsurfaceparams;

	for( int c = 0; c < min_ncomplex; c++ )
	{
		min_complexes[c]->applyParamsToComplex( p + nparams );
		nparams += min_complexes[c]->nparams();
	}
		
	if( do_sim && fitRho_activated && fabs(fitCoupling) > 1e-30 )
	{
#ifdef FIT_RHO_ONE_THICKNESS
		double thick_inner = p[nparams]; nparams++;
		double thick_outer = thick_inner;
#else
		double thick_inner = p[nparams]; nparams++;
		double thick_outer = p[nparams]; nparams++;
#endif
		for( surface_record *sRec = header; sRec; sRec = sRec->next )
			v += sRec->theSurface->rhoEnergy( p+sRec->temp_min_offset, min_simulation->PBC_vec, thick_inner, thick_outer );
	}	

	if( global_block->clathrinStructure && ! min_freeze_clathrin )
	{
		rotate_pts( min_simulation, p+nparams, p + nparams+3 );	
		nparams += 6;
	}

	for( surface_record *sRec = header; sRec; sRec = sRec->next )
		v += sRec->theSurface->cutEnergy( p+sRec->temp_min_offset );

	if( do_sim )
	{
		for( int cx = 0; cx < par_info.nc; cx++ )
		{
			int c = par_info.complexes[cx];
	
			v += min_complexes[c]->V( min_simulation );
			v += min_complexes[c]->AttachV( min_simulation );
		}
		v += Boxed_PP_V( min_simulation );
	}
	
	ParallelSum(&v,1);

	// AFTER PARALLEL SUM:

	if( do_sim && water_KV > 0 )
	{
		for( surface_record *sRec = header; sRec; sRec = sRec->next )
		{
			double Vi,Vo;
			sRec->theSurface->new_volume( &Vi, &Vo, sRec->r, min_simulation->alpha, NULL );
			double strain_i = (Vi-sRec->V0_i)/sRec->V0_i;
			double strain_o = (Vo-sRec->V0_o)/sRec->V0_o;
			double vi = 0.5 * water_KV * sRec->V0_i * strain_i * strain_i;
			double vo = 0.5 * water_KV * sRec->V0_o * strain_o * strain_o;
//			printf("en vol %le %le en: %le %le\n", Vi,Vo,vi,vo );

			if( global_block->restrain_volume_outside )
			{
				VV += vo;
				v += vo;
			}
			if( global_block->restrain_volume_inside )
			{
				v += vi;
				VV += vi;
			}	
//			printf("vol: %le vol_v: %le\n",  Vi, 0.5 * water_KV * sRec->V0 * pow( (Vi-sRec->V0)/sRec->V0, 2 ) );
		}
	}

#ifdef PARALLEL
	MPI_Bcast( &v, 1, MPI_DOUBLE, BASE_TASK, MPI_COMM_WORLD );
#endif
//	printf("v: %le\n", v);
	return v;
}



double surface_fdf( double *p, double *g)
{
#ifdef PARALLEL
	ParallelSyncComplexes( min_complexes, min_ncomplex );
#endif

	memset( g, 0, sizeof(double) * min_nparams );
	double v = 0;	
	int offset = 3; 
	
	int do_sim = 1;

	surface_record *header = NULL;
	if( min_simulation )
		header = min_simulation->allSurfaces;
	else
	{
		do_sim = 0;
		header = &one_srec;
	}
	
	for( surface_record *sRec = header; sRec; sRec = sRec->next )
	{
		int nv = sRec->theSurface->nv;
		p[offset+3*nv+0] = p[0];
		p[offset+3*nv+1] = p[1];
		p[offset+3*nv+2] = p[2];
		sRec->g = g+offset;
		sRec->r = p+offset;		

		sRec->theSurface->grad( p+offset, g+offset );
		offset += sRec->theSurface->nv*3+3;
	}
	int nparams = min_nsurfaceparams;	
	for( int c = 0; c < min_ncomplex; c++ )
	{
		int np = min_complexes[c]->nparams();
		min_complexes[c]->applyParamsToComplex( p + nparams );
		nparams += np;
	}
	
	//if( fitRho_activated  && fabs(fitCoupling) > 1e-30 )
	{
		double thick_inner, thick_outer;
		double *rho_g_i, *rho_g_o;
		int do_rho = 0;

		if( fitRho_activated && fabs(fitCoupling) > 1e-30 )
		{
			do_rho = 1;
#ifdef FIT_RHO_ONE_THICKNESS
			thick_inner = p[nparams];
			rho_g_i = g+nparams;
			*rho_g_i = 0;
	
			thick_outer = p[nparams]; 
			rho_g_o = g+nparams;
			*rho_g_o = 0;
	
			nparams++;
#else
			thick_inner = p[nparams];
			rho_g_i = g+nparams;
			*rho_g_i = 0;
			nparams++;
	
			thick_outer = p[nparams]; 
			rho_g_o = g+nparams;
			*rho_g_o = 0;
			nparams++;
#endif
		}

		int offset = 3; 

		if( global_block->clathrinStructure && ! min_freeze_clathrin )
		{
			// when the gradient is computed it sets the frc points for subsequent iteration. we need to update the base pdb as well at this point.
			set_saved_transform( min_simulation, p+nparams, p+nparams+3 );
			// this I'm sure hurts the BFGS. resets the angles
			p[nparams+0] = 0;
			p[nparams+1] = 0;
			p[nparams+2] = 0;
			p[nparams+3] = 0;
			p[nparams+4] = 0;
			p[nparams+5] = 0;
			min_simulation->clathrinGrad( g+nparams );
			nparams += 6;
		}

		for( surface_record *sRec = header; sRec; sRec = sRec->next )
		{
			int nv = sRec->theSurface->nv;
			p[offset+3*nv+0] = p[0];
			p[offset+3*nv+1] = p[1];
			p[offset+3*nv+2] = p[2];
			sRec->g = g+offset;
			sRec->r = p+offset;		


			if(do_rho )
			{
				v += sRec->theSurface->rhoGrad( p+offset, g+offset, min_simulation->PBC_vec, thick_inner, thick_outer, rho_g_i, rho_g_o );

				if( fix_thickness )
				{
					*rho_g_i = 0;
					*rho_g_o = 0;
				}
			}
			v += sRec->theSurface->cutGrad( p+offset, g+offset );
			offset += sRec->theSurface->nv*3+3;
		}
	}	
	
	for( int c = 0; c < min_ncomplex; c++ )
		memset( min_complexes[c]->save_grad, 0, sizeof(double) * 3 * min_complexes[c]->nsites );

#ifdef PARALLEL
	ParallelSyncComplexes( min_complexes, min_ncomplex );
#endif

	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];


		int nsites = min_complexes[c]->nsites;

		double rg[3*nsites];
		double ng[3*nsites];

		memset( rg, 0, sizeof(double) * 3 * nsites );
		memset( ng, 0, sizeof(double) * 3 * nsites );
		// derivative wrt positions, normals
		double le = min_complexes[c]->grad( min_simulation, rg, ng );
			
		v += le;
	
		if( min_complexes[c]->bound )
		{
			// attachment coordinates need gradients.
			
			for( int a = 0; a < min_complexes[c]->nattach; a++ )
			{
				double point_grad[2]={0,0};
				double point_rgrad[3] = { 0,0,0};
				surface_record *sRec = min_simulation->fetch( min_complexes[c]->sid[a] );
				
				// function updates the gradient of the energy for the surface coordinates.

				sRec->theSurface->pointGradient( min_complexes[c]->grad_fs[a], min_complexes[c]->grad_puv[2*a+0], 
										     min_complexes[c]->grad_puv[2*a+1],
										     sRec->r, sRec->g, point_grad, rg+3*a, ng+3*a );  			
				min_complexes[c]->save_grad[2*a+0] += point_grad[0];		
				min_complexes[c]->save_grad[2*a+1] += point_grad[1];		
			}

			for( int a = min_complexes[c]->nattach; a < min_complexes[c]->nsites; a++ )
			{
				min_complexes[c]->save_grad[3*a+0] += rg[3*a+0];	
				min_complexes[c]->save_grad[3*a+1] += rg[3*a+1];	
				min_complexes[c]->save_grad[3*a+2] += rg[3*a+2];	
			}			
		}
	
		v += min_complexes[c]->AttachG( min_simulation, min_complexes[c]->save_grad );
	}

	nparams = min_nsurfaceparams;

	for( int c = 0; c < min_ncomplex; c++ )
	{
		int np = min_complexes[c]->nparams();
		if( min_complexes[c]->bound )
		{
			// attachment coordinates need gradients.
		
			for( int a = 0; a < min_complexes[c]->nattach; a++ )
			{
				int np = min_complexes[c]->nparams();

				double M[4] = { min_complexes[c]->coord_transform[4*a+0],
						min_complexes[c]->coord_transform[4*a+1],
						min_complexes[c]->coord_transform[4*a+2],
						min_complexes[c]->coord_transform[4*a+3] };
				double det = M[0]*M[3]-M[1]*M[2];
				double MINV[4] = { M[3]/det, -M[1]/det, -M[2]/det, M[0]/det };

#if 1 
				g[nparams+2*a+0] += min_complexes[c]->save_grad[2*a+0] * M[0];
				g[nparams+2*a+0] += min_complexes[c]->save_grad[2*a+1] * M[2];
				g[nparams+2*a+1] += min_complexes[c]->save_grad[2*a+0] * M[1];
				g[nparams+2*a+1] += min_complexes[c]->save_grad[2*a+1] * M[3];
#else
				g[nparams+2*a+0] += min_complexes[c]->save_grad[2*a+0];
				g[nparams+2*a+1] += min_complexes[c]->save_grad[2*a+1];

#endif

			}
	 	
			int poff = 2 * min_complexes[c]->nattach;
			for( int a = min_complexes[c]->nattach; a < min_complexes[c]->nsites; a++ )
			{
				int ax = a - min_complexes[c]->nattach;
				g[nparams+poff+3*ax+0] += min_complexes[c]->save_grad[3*a+0];
				g[nparams+poff+3*ax+1] += min_complexes[c]->save_grad[3*a+1];
				g[nparams+poff+3*ax+2] += min_complexes[c]->save_grad[3*a+2];
			}			
		}

		nparams += np;
	}	
	
	for( int c = 0; c < min_ncomplex; c++ )
		memset( min_complexes[c]->save_grad, 0, sizeof(double) * 3 * min_complexes[c]->nsites );
		
	v += Boxed_PP_G( min_simulation);
	
	nparams = min_nsurfaceparams;
	for( int c = 0; c < min_ncomplex; c++ )
	{
		int np = min_complexes[c]->nparams();

		for( int a = 0; a < min_complexes[c]->nattach;  a++ )
		{
				double M[4] = { min_complexes[c]->coord_transform[4*a+0],
						min_complexes[c]->coord_transform[4*a+1],
						min_complexes[c]->coord_transform[4*a+2],
						min_complexes[c]->coord_transform[4*a+3] };
				double det = M[0]*M[3]-M[1]*M[2];
				double MINV[4] = { M[3]/det, -M[1]/det, -M[2]/det, M[0]/det };

				g[nparams+2*a+0] += min_complexes[c]->save_grad[2*a+0] * M[0];
				g[nparams+2*a+0] += min_complexes[c]->save_grad[2*a+1] * M[2];
				g[nparams+2*a+1] += min_complexes[c]->save_grad[2*a+0] * M[1];
				g[nparams+2*a+1] += min_complexes[c]->save_grad[2*a+1] * M[3];
		}
	 	int poff = 2 * min_complexes[c]->nattach;
		for( int a = min_complexes[c]->nattach; a < min_complexes[c]->nsites; a++ )
		{
			int ax = a - min_complexes[c]->nattach;
			g[nparams+poff+3*ax+0] += min_complexes[c]->save_grad[3*a+0];
			g[nparams+poff+3*ax+1] += min_complexes[c]->save_grad[3*a+1];
			g[nparams+poff+3*ax+2] += min_complexes[c]->save_grad[3*a+2];
		}
		nparams += np;
	}

	// for now zero box dimension gradient.
	
	g[0] = 0;
	g[1] = 0;
	g[2] = 0;

	for( surface_record *sRec = header; sRec; sRec = sRec->next )
	{
		g[sRec->temp_min_offset+3*sRec->theSurface->nv+0] = 0;
		g[sRec->temp_min_offset+3*sRec->theSurface->nv+1] = 0;
		g[sRec->temp_min_offset+3*sRec->theSurface->nv+2] = 0;

		if( min_z_only )
		{
			for( int x = 0; x < sRec->theSurface->nv; x++ )
			{
				g[sRec->temp_min_offset+3*x+0] = 0;
				g[sRec->temp_min_offset+3*x+1] = 0;
			}
		}
	}
	

	if( do_freeze_membrane )
	{
		for( int x = 0; x < min_nsurfaceparams; x++ )
			g[x] = 0;
	}

	for( surface_record *sRec = header; sRec; sRec = sRec->next )
		v += sRec->theSurface->energy( p+sRec->temp_min_offset, NULL );
	ParallelSum( &v, 1 );
	
	ParallelSum( g, min_nparams );

	if( water_KV > 0 )
	{
		for( surface_record *sRec = header; sRec; sRec = sRec->next )
		{
			double * dvol = (double *)malloc( sizeof(double) * 3 * sRec->theSurface->nv );
			memset( dvol, 0, sizeof(double) * 3 * sRec->theSurface->nv );
			double Vi,Vo;
			sRec->theSurface->new_volume( &Vi, &Vo, sRec->r, min_simulation->alpha, dvol );
				
			double strain_i = (Vi-sRec->V0_i)/sRec->V0_i;
			double strain_o = (Vo-sRec->V0_o)/sRec->V0_o;

			double vi = 0.5 * water_KV * sRec->V0_i * strain_i * strain_i;
			double vo = 0.5 * water_KV * sRec->V0_o * strain_o * strain_o;
			if( global_block->restrain_volume_inside )
			for( int p = 0; p < sRec->theSurface->nv; p++ )
			{
				g[sRec->temp_min_offset+3*p+0] += water_KV * strain_i * dvol[3*p+0] ;		
				g[sRec->temp_min_offset+3*p+1] += water_KV * strain_i * dvol[3*p+1] ;		
				g[sRec->temp_min_offset+3*p+2] += water_KV * strain_i * dvol[3*p+2] ;		
			}

			
			if( global_block->restrain_volume_outside )
			for( int p = 0; p < sRec->theSurface->nv; p++ )
			{
				g[sRec->temp_min_offset+3*p+0] -= water_KV * strain_o * dvol[3*p+0] ;		
				g[sRec->temp_min_offset+3*p+1] -= water_KV * strain_o * dvol[3*p+1] ;		
				g[sRec->temp_min_offset+3*p+2] -= water_KV * strain_o * dvol[3*p+2] ;		
			}
//			printf("grad vol %le %le en: %le %le\n", Vi,Vo,vi,vo );

			free(dvol);

			if( global_block->restrain_volume_inside )
			{
				VV += vi;
				v += vi;
			}
			if( global_block->restrain_volume_outside )
			{
				VV += vo;
				v += vo;
			}
		}
	}

#if 0 	
	MPI_Barrier(MPI_COMM_WORLD);
	if( par_info.my_id == BASE_TASK )
	{	
		FILE *gradf = fopen("grad.txt","w");
		for( int c = 0; c < min_nparams; c++ )
			fprintf(gradf, "%d %.14le\n", c, g[c] );
		fclose(gradf);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	exit(1);
#endif

#ifdef PARALLEL
	MPI_Bcast( &v, 1, MPI_DOUBLE, BASE_TASK, MPI_COMM_WORLD );
#endif
	return v;

}

void full_fd_test( double *p )
{
	double deps = 1e-8;

	double *g = (double *)malloc( sizeof(double) * min_nparams );
	memset( g, 0, sizeof(double) * min_nparams );

	double e0 = surface_fdf( p, g );
	printf("Finite difference test.\n");

	for( int xp  = 0; xp < min_nparams; xp++ )
	{
		double de_pm[2];
	
		for( int im = 0; im < 2; im++ )
		{
			p[xp] += deps * (im == 0 ? 1 : -1);	
		
			double v = surface_f( p );
			
			de_pm[im] = v;

			p[xp] -= deps * (im == 0 ? 1 : -1);	
		}

		double fr_error = fabs(((de_pm[0]-de_pm[1])/(2*deps) - g[xp])/((de_pm[0]-de_pm[1])/(2*deps)));	
	
		if( fabs( (de_pm[0]-de_pm[1])/(2*deps) - g[xp]) > 1e-5 && fr_error > 1e-5 )
			printf("parm %d fd_der %.14le g %.14le del %.14le BAD\n", xp, (de_pm[0]-de_pm[1])/(2*deps), g[xp],  (de_pm[0]-de_pm[1])/(2*deps) - g[xp] );
		else if( fabs( (de_pm[0]-de_pm[1])/(2*deps)) > 1e-10 )	
			printf("parm %d fd_der %.14le g %.14le del %.14le OK\n", xp, (de_pm[0]-de_pm[1])/(2*deps), g[xp],  (de_pm[0]-de_pm[1])/(2*deps) - g[xp] );
	}
	

	free(g);

//	exit(1);
}

void fd_test( double *p )
{
	double deps = 1e-6;

	double *g = (double *)malloc( sizeof(double) * min_nparams );
	memset( g, 0, sizeof(double) * min_nparams );

	double e0 = surface_fdf( p, g );
	printf("Finite difference test.\n");
	for( double eps = 0; eps < 20 * deps; eps += deps)
	{
		double expec = 0;
		for( int x = 0; x < min_nparams;x++ )
		{
			p[x] -= eps * g[x];
			expec += -g[x] * g[x] * eps;
		}

		double v = surface_f( p );
		printf("%le %.14le %.14le del %.14le\n", eps, v-e0, expec, (v-(e0+expec)) );
		
		for( int x = 0; x < min_nparams;x++ )
			p[x] += eps * g[x];
	}

	free(g);
}

void Simulation::minimize( int freeze_membrane, int freeze_clathrin  )
{
	do_freeze_membrane = freeze_membrane;
	min_freeze_clathrin = freeze_clathrin;
	int prev_enable = enable_elastic_interior; 

	enable_elastic_interior = 1;


	min_simulation = this;
	min_ncomplex = ncomplex;	
	min_complexes = allComplexes;
	min_z_only = 0;	
	if( global_block->z_only ) min_z_only = 1;

	if( global_block->fitThickness >= 0 ) {
		fix_thickness = 1;
		cur_rho_thickness[0] = global_block->fitThickness;
		cur_rho_thickness[1] = global_block->fitThickness;
	}
#ifdef PARALLEL
	ParallelSyncComplexes( min_complexes, min_ncomplex );
#endif
	
	int num_params = 3; // alphas

	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
		num_params += sRec->theSurface->nv*3+3;

	min_nsurfaceparams = num_params;

	for( int c = 0; c < ncomplex; c++ )
		num_params += allComplexes[c]->nparams();

	if( fitRho_activated && fabs(fitCoupling) > 1e-30  )
	{
#ifdef FIT_RHO_ONE_THICKNESS
		num_params += 1; // same thickness
#else
		num_params += 2; // inner and outer thickness.
#endif
	}

	if( global_block->clathrinStructure && !min_freeze_clathrin)
		num_params += 6;


	min_nparams = num_params;
	
	double *p = (double *)malloc( sizeof(double) * num_params );
	double *g = (double *)malloc( sizeof(double) * num_params );

	p[0] = alpha[0];
	p[1] = alpha[1];
	p[2] = alpha[2];
	int offset = 3;
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		int nv = sRec->theSurface->nv;

		memcpy( p+offset, sRec->r, sizeof(double) * 3 * nv );
		sRec->temp_min_offset = offset;
		sRec->temp_r = sRec->r;
		sRec->temp_g = sRec->g;
		sRec->r = p + offset;
		sRec->g = g + offset;

		p[offset+3*nv+0] = alpha[0];
		p[offset+3*nv+1] = alpha[1];
		p[offset+3*nv+2] = alpha[2];

		offset += 3*nv+3;
	}

	int tp = min_nsurfaceparams;

	for( int c = 0; c < ncomplex; c++ )
	{
		allComplexes[c]->getParamsFromComplex( p + tp );
		tp += allComplexes[c]->nparams();
	}	

	int thickness_ptr = 0;

	if( fitRho_activated && fabs(fitCoupling) > 1e-30  )
	{
		printf("Current thickness: %lf %lf\n", cur_rho_thickness[0], cur_rho_thickness[1] );
		thickness_ptr = tp;

#ifdef FIT_RHO_ONE_THICKNESS
		p[tp] =cur_rho_thickness[0]; tp++; 
#else
		p[tp] =cur_rho_thickness[0]; tp++; 
		p[tp] =cur_rho_thickness[1]; tp++; 
#endif
	}
	
	if( global_block->clathrinStructure && !min_freeze_clathrin )
	{
		p[tp] = 0; tp++;
		p[tp] = 0; tp++;
		p[tp] = 0; tp++;
		p[tp] = 0; tp++;
		p[tp] = 0; tp++;
		p[tp] = 0; tp++;
	}
	// derivative might be zero (absolutely)

	surface_fdf(p,g);
	double mag_init = 0;
	for( int p = 0; p < num_params; p++ )
		mag_init += g[p]*g[p];
	

	int nsteps = 50;
	
	int use_m = nsteps;
	if( use_m > num_params )
		use_m = num_params;
	double e_init = surface_f(p);
	if( global_block->fdiff_check )
		full_fd_test(p);
	printf("Entering minimize with e_init: %le\n", e_init );
	l_bfgs_setup( use_m, num_params, p, 1.0, surface_f, surface_fdf); 

	printf("Minimize IN: V: %.14le VA: %lf VC: %lf VV: %lf\n", e_init, VA, VC, VV );
	if( mag_init > 1e-20 )
	{
		for( int x = 0; x < nsteps; x++ )
		{
			if( ! l_bfgs_iteration( p ) ) { break; }
	
//		if( x %10 == 0 )
//			printf("Sub iteration %d, V: %le\n", x, surface_f(p) );
		}
	}	
	else
	{
		printf("Initial gradient zero.\n");
	}	
	l_bfgs_clear();
	double v =surface_fdf(p,g);
	double e = surface_f(p);

	double rms = 0;
	for( int x = 0; x < num_params; x++ )
		rms += g[x]*g[x];
	rms /= num_params;
	rms = sqrt(rms);

//	full_fd_test(p);

	printf("Minimize OUT: VG: %.14le VV: %.14le VA: %lf VC: %lf VV: %lf grad rms %le\n", v, e, VA, VC, VV, rms );

	enable_elastic_interior = prev_enable;

	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		int nv = sRec->theSurface->nv;

		memcpy( sRec->temp_r, sRec->r, sizeof(double) * 3 * nv );
		sRec->temp_r[3*nv+0] = p[0];
		sRec->temp_r[3*nv+1] = p[1];
		sRec->temp_r[3*nv+2] = p[2];
		sRec->r = sRec->temp_r;
		sRec->g = sRec->temp_g;
	}
	
	for( int c = 0; c < ncomplex; c++ )
		allComplexes[c]->refresh(this);

	if( fitRho_activated && fabs(fitCoupling) > 1e-30  )
	{
#ifdef FIT_RHO_ONE_THICKNESS
		cur_rho_thickness[0] = p[thickness_ptr];		
		cur_rho_thickness[1] = p[thickness_ptr];		
#else
		cur_rho_thickness[0] = p[thickness_ptr];		
		cur_rho_thickness[1] = p[thickness_ptr+1];		
#endif
	}
	
	free(p);
	free(g);
}

/*
 * Minimizes the curvature energy of the surface only, no frills.
 * */

void surface::minimize( double *r )
{
	min_simulation = NULL;
	min_ncomplex = 0;	
	min_complexes = NULL;

	int num_params = 3 + (nv * 3+3); // alphas

	min_nsurfaceparams = num_params;
	min_nparams = num_params;
	
	min_z_only = 0;	
	if( global_block->z_only ) min_z_only = 1;
	double *p = (double *)malloc( sizeof(double) * num_params );
	double *g = (double *)malloc( sizeof(double) * num_params );
	
	p[0] = r[3*nv+0];
	p[1] = r[3*nv+1];
	p[2] = r[3*nv+2];

	int offset = 3;
	
	one_srec.theSurface = this;
	one_srec.temp_r = (double *)malloc( sizeof(double ) * 3 * nv );
	one_srec.r = p + offset;
	one_srec.g = g + offset;
	one_srec.temp_g = (double *)malloc( sizeof(double) * 3 * nv );
	one_srec.next = NULL;
	
	memcpy( p + offset, r, sizeof(double ) * (3 * nv+3) );

	offset += 3*nv+3;

	int tp = min_nsurfaceparams;

	// derivative might be zero (absolutely)

	surface_fdf(p,g);
	double mag_init = 0;
	for( int p = 0; p < num_params; p++ )
		mag_init += g[p]*g[p];
	
	int nsteps = 100;
	
	int use_m = nsteps;
	if( use_m > num_params )
		use_m = num_params;
	double e_init = surface_f(p);
	//full_fd_test(p);
	printf("Entering minimize with e_init: %le\n", e_init );
	l_bfgs_setup( use_m, num_params, p, 1.0, surface_f, surface_fdf); 

	printf("Minimize IN: V: %.14le VA: %lf VC: %lf VV: %lf\n", e_init, VA, VC, VV );
	if( mag_init > 1e-20 )
	{
		for( int x = 0; x < nsteps; x++ )
		{
			if( ! l_bfgs_iteration( p ) ) { break; }
	
//		if( x %10 == 0 )
//			printf("Sub iteration %d, V: %le\n", x, surface_f(p) );
		}
	}	
	else
	{
		printf("Initial gradient zero.\n");
	}	
	l_bfgs_clear();
	double v =surface_fdf(p,g);
	double e = surface_f(p);

	double rms = 0;
	for( int x = 0; x < num_params; x++ )
		rms += g[x]*g[x];
	rms /= num_params;
	rms = sqrt(rms);

	printf("Minimize OUT: VG: %.14le VV: %.14le VA: %lf VC: %lf VV: %lf grad rms %le\n", v, e, VA, VC, VV, rms );
	fflush(stdout);
	memcpy( r, one_srec.r, sizeof(double)* 3* nv );

	r[3*nv+0] = p[0];
	r[3*nv+1] = p[1];
	r[3*nv+2] = p[2];


	free( one_srec.temp_g );	
	free( one_srec.temp_r );	
	
	free(p);
	free(g);
}
	


