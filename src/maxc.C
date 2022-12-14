#include "interp.h"
#include "uv_map.h"
#include "l-bfgs.h"
#include <math.h>
#include "maxc.h"

// a simple utility to find saddle curvature on a surface.

// these are the initial points from which we diverge.
static int src_f, out_f;
static double src_uv[2], out_uv[2];
static double *rsurf = NULL;
static surface *minSurface = NULL;
static double irreg_penalty = 1e5;
static int c_mode = 0;

double maxc_f( double *parms )
{
	int curf=src_f;

	double uv1[2] = { src_uv[0], src_uv[1] };
	double duv1[2] = { parms[0], parms[1] };
	int f_1,nf=src_f;

	do {
		f_1 = nf;
		nf = minSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf  ); 
	} while( nf != f_1 );

	double k;
	double cvec1[2],cvec2[2],c1,c2;
//	printf("f Eval at %d %.14le %.14le\n", f_1, uv1[0], uv1[1] );
	minSurface->c( f_1, uv1[0], uv1[1], rsurf,  &k, cvec1, cvec2, &c1, &c2 );
	double rpt[3],nrm[3];
	minSurface->evaluateRNRM( f_1, uv1[0], uv1[1], rpt, nrm, rsurf );	
//	printf("H %lf %lf %lf c1: %le c2: %le\n", rpt[0], rpt[1], rpt[2], c1, c2 );

	out_f = f_1;
	out_uv[0] = uv1[0];
	out_uv[1] = uv1[1];

	double pot = 0;

	switch( c_mode )
	{
		case MAX_C_GAUSS:
			pot += c1*c2;
			break;
		case MAX_C_POS:
			// max positive.
			if( c1 > c2 )
				pot += -c1;
			else
				pot += -c2;
			break;
		case MAX_C_NEG:
			// max neg
			if( c1 < c2 )
				pot += c1;
			else
				pot += c2;
			break;
		case MAX_C_TOT_POS:
			pot +=  c1+c2;
			break;
		case MAX_C_TOT_NEG:
			pot += -(c1+c2);
			break;
		case MAX_C_MIN:
			pot += (c1+c2)*(c1+c2);
			break;
		default:
			printf("Min C mode %d not implemented.\n", c_mode );
			break;
	}

	if( f_1 >= minSurface->nf_faces )
	{
		pot += irreg_penalty * pow(1-uv1[0]-uv1[1],2.0);	
	} 

	

//	printf("f: %.14le --- %d nf_reg: %d %.14le %.14le\n", c1*c2, f_1, minSurface->nf_faces, uv1[0], uv1[1] );
	return pot;	
}

double maxc_fdf( double *parms, double *parms_g )
{
	int curf=src_f;

	double uv1[2] = { src_uv[0], src_uv[1] };
	double duv1[2] = { parms[0], parms[1] };
	double M[4] = { 1,0,0,1};
	int f_1,nf=src_f;
	double null_mom[2]={0,0};
	do {
		f_1 = nf;
		nf = minSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf, null_mom, M ); 
	} while( nf != f_1 );

	double k;

	double cvec1[2],cvec2[2],c1,c2;
	minSurface->c( f_1, uv1[0], uv1[1], rsurf,  &k, cvec1, cvec2, &c1, &c2 );

	double d_c1_duv[2];
	double d_c2_duv[2];

	// these are in terms of the local coordinate system. ``M'' transforms between the original and the local.

//	printf("fdf Eval at %d %.14le %.14le\n", f_1, uv1[0], uv1[1] );
	minSurface->d_c_duv( rsurf, f_1, uv1[0], uv1[1], d_c1_duv, d_c2_duv );

	double rpt[3], nrm[3];
	minSurface->evaluateRNRM( f_1, uv1[0], uv1[1], rpt, nrm, rsurf );
//	printf("N %lf %lf %lf c1: %le c2: %le d_c1_duv: %le %le\n", rpt[0], rpt[1], rpt[2], c1, c2, d_c1_duv[0], d_c1_duv[1] );
	double pot = 0;

	parms_g[0] = 0;
	parms_g[1] = 0;

	double de_dc1 = 0;
	double de_dc2 = 0;

	switch( c_mode )
	{
		case MAX_C_GAUSS:
			de_dc1 = c2;
			de_dc2 = c1;
			pot += c1*c2;
			break;
		case MAX_C_POS:
			// max positive.
			if( c1 > c2 )
				de_dc1 = -1;
			else
				de_dc2 = -1;
			if( c1 > c2 )
				pot += -c1;
			else
				pot += -c2;
			break;
		case MAX_C_NEG:
			// max neg
			if( c1 < c2 )
				de_dc1 = 1;
			else
				de_dc2 = 1;
			if( c1 < c2 )
				pot += c1;
			else
				pot += c2;
			break;
		case MAX_C_TOT_POS:
			de_dc1 = 1;
			de_dc2 = 1;
			pot +=  c1+c2;
			break;
		case MAX_C_TOT_NEG:
			de_dc1 = -1;
			de_dc2 = -1;
			pot += -(c1+c2);
			break;
		case MAX_C_MIN:
			de_dc1 = 2*(c1+c2);
			de_dc2 = 2*(c1+c2);
			pot += (c1+c2)*(c1+c2);
			break;
		default:
			printf("Min C mode %d not implemented.\n", c_mode );
			break;
	}


	parms_g[0] += M[0*2+0] * d_c1_duv[0] * de_dc1;	
	parms_g[0] += M[1*2+0] * d_c1_duv[1] * de_dc1;	
	
	parms_g[0] += M[0*2+0] * d_c2_duv[0] * de_dc2;	
	parms_g[0] += M[1*2+0] * d_c2_duv[1] * de_dc2;	
	
	parms_g[1] += M[0*2+1] * d_c1_duv[0] * de_dc1;	
	parms_g[1] += M[1*2+1] * d_c1_duv[1] * de_dc1;	
	
	parms_g[1] += M[0*2+1] * d_c2_duv[0] * de_dc2;	
	parms_g[1] += M[1*2+1] * d_c2_duv[1] * de_dc2;	
	

#if 0
	printf("g: %le %le\n", parms_g[0], parms_g[1] );
	for( double deps = -0.02; deps <= 0.02; deps += 0.001 )
	{
		double uv1[2] = { src_uv[0], src_uv[1] };
		double M[4] = { 1,0,0,1};
		int f_1,nf=src_f;
		double null_mom[2]={0,0};
		double duv[2] = { deps * parms_g[0], deps * parms_g[1] };

		do {
			f_1 = nf;
			nf = minSurface->nextFace( f_1, uv1+0, uv1+1, duv+0, duv+1, rsurf, null_mom, M ); 
		} while( nf != f_1 );

		double k, cvec1[2], cvec2[2], c1,c2;
		minSurface->c( f_1, uv1[0], uv1[1], rsurf,  &k, cvec1, cvec2, &c1, &c2 );

		printf("deps: %le f: %d u %lf v %lf c1 %.14le c2 %.14le\n", deps, f_1, uv1[0], uv1[1],
										c1, c2 );
	}
#endif

//	printf("fdf: %.14le --- %d nf_reg: %d %.14le %.14le %.14le %.14le\n", c1*c2, f_1, minSurface->nf_faces, uv1[0], uv1[1], parms_g[0], parms_g[1] );

	if( f_1 >= minSurface->nf_faces )
	{
		pot += irreg_penalty * pow(1-uv1[0]-uv1[1],2.0);	

		parms_g[0] += irreg_penalty * ( 1-uv1[0]-uv1[1]);// * M[0*2+0];
		parms_g[0] += irreg_penalty * ( 1-uv1[0]-uv1[1]);// * M[1*2+0];
		parms_g[1] += irreg_penalty * ( 1-uv1[0]-uv1[1]);// * M[0*2+1];
		parms_g[1] += irreg_penalty * ( 1-uv1[0]-uv1[1]);// * M[1*2+1];
	} 

	return pot;	
} 

/*
	// mode 0: max gauss
	// mode 1; max individual pos
	// mode 2: max individual neg
	// mode 3: min C
	// mode 4: max C
	// mode 5: min C^2 (flat)

	curvature isn't differentiable along a face: makes minimization challenging.

*/

void max_c( surface *theSurface, int *f, double *u, double *v, int niter, double *rsurf_in, int mode )
{
	c_mode = mode;
	disable_random_uv_step();
	minSurface = theSurface;
	rsurf = rsurf_in;

	double *duv = (double *)malloc( sizeof(double) * 2 );
	
	double outk = 0;
	
	for( int it = 0; it < niter; it++ )
	{
		// where we evaluate the gradient.

		src_f = *f;
		src_uv[0] = *u;
		src_uv[1] = *v;
		duv[0] = 0;
		duv[1] = 0;
		double g[2];
		double k0 = maxc_fdf( duv, g );
		double best_f = k0;
		double deps = 1;	

		l_bfgs_setup( 4, 2, duv, 0.1, maxc_f, maxc_fdf); 



	
		
		int nsteps = 200;

		for( int x = 0; x < nsteps; x++ )
		{
			if( ! l_bfgs_iteration( duv ) ) { break; }
		}

		outk = maxc_f( duv );
				

		// do a line minimization just to escape very small wells that the BFGS routines aren't made for.
		
		
		k0 = maxc_fdf( duv, g );
		best_f = k0;
		
		src_f = *f = out_f;
		src_uv[0] = *u = out_uv[0];				
		src_uv[1] = *v = out_uv[1];			
		
		int line_found = 0;
		int do_hop = 0;
#if 0 

//		printf("FEPS %d %d\n", out_f, minSurface->nf_faces );	

		for( int ieps = 0; ieps < 100; ieps++ )
		{
		 	double trial[2] = {  duv[0] + deps * ieps * g[0],
					     duv[1] + deps * ieps * g[1] };
			double keval = maxc_f(trial);
//			printf("eps: %le eval: %le expec: %le\n", ieps * deps,
//				keval-k0, - deps * ieps * pow(g[0],2.0) - deps * ieps * pow(g[1],2.0)  );  

			if( ieps == 1 )
			{
				double rat = (keval-k0)/( - deps * ieps * pow(g[0],2.0) - deps * ieps * pow(g[1],2.0));

				if( rat < -0.5 && rat >= -2 )
				{
					do_hop = 1;
//					printf("Reverse gradient. check here!\n");
//					double rev_grad = maxc_fdf( duv, g );
				}
//				break;
			}

			if( keval < best_f )
			{
				line_found = 1;
				best_f = keval;
				*f = out_f;
				*u = out_uv[0];				
				*v = out_uv[1];				
			}
		}
#endif
		
		// if we stopped on an edge, do a monte carlo optimization to try to get out of it.

		
	//	if( !line_found && (fabs(*u) < 1e-2 || fabs(*v) < 1e-2 || fabs(1-*u-*v) < 1e-2) )	
		//if( do_hop && !line_found )
		{
//			printf("Hopping!\n");
			int n_mc_hops = 100;

			double duv[2] = { 4 * (0.5 - rand() / (double)(RAND_MAX)),
					  4 * (0.5 - rand() / (double)(RAND_MAX)) };
			double trial = maxc_f(duv);

			if( trial < best_f )
			{
				*f = src_f = out_f;
				*u = src_uv[0] = out_uv[0];
				*v = src_uv[1] = out_uv[1];
				best_f = trial;
			} 
		}

		l_bfgs_clear();
	}

	free(duv);
	
	
	double k, cvec1[2],cvec2[2],c1,c2;
	minSurface->c( *f, *u, *v, rsurf,  &k, cvec1, cvec2, &c1, &c2 );

	printf("Final c_1/c_2 : %le %le f: %d uv: %le %le \n", c1, c2, *f, *u, *v );
	enable_random_uv_step();
}
