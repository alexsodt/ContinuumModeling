#define __gatherc__
/*
 * Routines for gathering nearest-point curvature from a surface.
 * Initial purpose is to gather lipid's sampled curvature given a fit to the surface.
 * */
#include "simulation.h"
#include "interp.h"
#include "gather.h"
#include "M_matrix.h"
#include "meshCollisionLib.h"
#include "io_mol_read.h"
#include "input.h"
#include "mutil.h"
#include "pcomplex.h"
#include <math.h>
#include "maxc.h"

extern double SUMK;

double Simulation::nearCurvature( double *rpt, double *c_out, double *k_out, double *dp_out, double *dz_out, int *leaflet_out, double *out_nrm )
{
	double **M;
	int mlow,mhigh;

	getM( &M, &mlow, &mhigh );

	double min_dist = 1e10;

	double min_c = 0;
	double min_k = 0;
	double min_dp = 0;
	double min_dz = 0;
	int min_leaflet = 0; // 0 is lower.
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		int f;
		double u, v;
		double distance;

		sRec->theSurface->nearPointOnBoxedSurface( rpt, &f, &u, &v, M, mlow, mhigh, &distance );					

		if( distance < min_dist )
		{
			min_c = sRec->theSurface->c( f, u, v, sRec->r, &min_k ); 
			min_dist = distance;

			// evaluate upper/lower

			double fpt[3],nrm[3];
			sRec->theSurface->evaluateRNRM( f, u, v, fpt, nrm, sRec->r );		

			if( sRec->gather_flip )
			{
				min_c *= -1;
				nrm[0] *= -1;
				nrm[1] *= -1;
				nrm[2] *= -1;
			}

			double dr[3] = { rpt[0] - fpt[0], rpt[1] - fpt[1], rpt[2] - fpt[2] };
			wrapPBC( dr, alpha );	

			double l = normalize(dr);
			min_dp = dr[0]*nrm[0]+dr[1]*nrm[1]+dr[2]*nrm[2];
			min_dz = l;
	
			if( min_dp < 0 )
			{
				min_leaflet = 0;
				min_c *= -1;
			}
			else
				min_leaflet = 1;	 

			out_nrm[0] = nrm[0];
			out_nrm[1] = nrm[1];
			out_nrm[2] = nrm[2];
		}
	}	

	if( c_out ) *c_out = min_c;
	if( k_out ) *k_out = min_k;
	if( dp_out ) *dp_out = min_dp;
	if( dz_out ) *dz_out = min_dz;
	if( leaflet_out ) *leaflet_out = min_leaflet;

	return min_c;
}

static double fusion_pore_center[3] = {0,0,0};
static double fusion_pore_r 	    = 0;
static double fusion_pore_h	    = 0;
static int is_fusion_pore = 0;

int Simulation::getFusionPoreData( double *center, double *r, double *h )
{
	double vec1[3],vec2[3];	
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		double av_c = 0;
		double av_k = 0;
		double sum_g = 0;	
	
		double uv_max_neg[2];
		int f_max_neg;	
		double max_neg = 1e10;
	
		double uv_max_pos[2];
		int f_max_pos;	
		double max_pos = -1e10;
		
		double uv_max_k[2];
		int f_max_k;	
		double max_k = 1e10;

		// compute the energy to get sum K.

		sRec->theSurface->energy( sRec->r,NULL);

		double theK = SUMK;

		for( int f = 0; f < sRec->theSurface->nt; f++ )
		{
			double k;
			double g = sRec->theSurface->g(f, 1.0/3.0, 1.0/3.0, sRec->r ); 
			double c1,c2;
			double c = sRec->theSurface->c(f, 1.0/3.0, 1.0/3.0, sRec->r, &k, vec1, vec2, &c1, &c2 );
			av_c += g * c; 
			av_k += k * g;
			sum_g += g;

			if( c1 > max_pos )
			{
				max_pos = c1;
				uv_max_pos[0] = 1.0/3.0;
				uv_max_pos[1] = 1.0/3.0;
				f_max_pos = f;
			}
			else if( c2 > max_pos )
			{
				max_pos = c2;
				uv_max_pos[0] = 1.0/3.0;
				uv_max_pos[1] = 1.0/3.0;
				f_max_pos = f;
			}
			
			if( c1 < max_neg )
			{
				max_neg = c1;
				uv_max_neg[0] = 1.0/3.0;
				uv_max_neg[1] = 1.0/3.0;
				f_max_neg = f;
			}
			if( c2 < max_neg )
			{
				max_neg = c2;
				uv_max_neg[0] = 1.0/3.0;
				uv_max_neg[1] = 1.0/3.0;
				f_max_neg = f;
			}

			if( k < max_k )
			{
				max_k = c;
				uv_max_k[0] = 1.0/3.0;
				uv_max_k[1] = 1.0/3.0;
				f_max_k = f;
			}
		}

		double *rsurf = sRec->r;
		double chi_check = theK / (-4*M_PI);
	
		printf("K_tot: %lf\n", theK );	
		if( fabs(chi_check-1) < 0.1 )
	 	{
			printf("Fusion pore detected. Setting sign...");
			double rpt[3], rnrm[3];
			sRec->theSurface->evaluateRNRM( 0, 1.0/3.0, 1.0/3.0, rpt, rnrm, sRec->r );

			sRec->gather_flip = 0;

			if( fabs(rnrm[2]) > 0.7 )
			{  
				if( rnrm[2] < -0.7 && rpt[2] < 0 )
					sRec->gather_flip  = 1;
				else if( rnrm[2] > 0.7 && rpt[2] > 0 )
					sRec->gather_flip  = 1;
			}
			else 
			{
				double dp = rpt[0] * rnrm[0] + rpt[1] * rnrm[1] + rpt[2] * rnrm[2];

				if( dp < 0 )
					sRec->gather_flip = 1;
			}
			
			double nrm_junk[3];
			double max_saddle[3];
			double max_pos[3];
			double max_neg[3];
			double k;
			double cpos[2];
			double cneg[2];
			double ck[2];

			max_c( sRec->theSurface, &f_max_k, uv_max_k+0, uv_max_k+1, 1000, rsurf, MAX_C_GAUSS );     
			sRec->theSurface->c(f_max_k, uv_max_k[0], uv_max_k[1], rsurf, &k, vec1,vec2, ck+0, ck+1 );
		
			double du[3], dv[3];
			sRec->theSurface->ru( f_max_k, uv_max_k[0], uv_max_k[1], rsurf, du );	
			sRec->theSurface->rv( f_max_k, uv_max_k[0], uv_max_k[1], rsurf, dv );	

			double cvec1[3] = { vec1[0] * du[0] + vec1[1] * dv[0],
					    vec1[0] * du[1] + vec1[1] * dv[1], 
					    vec1[0] * du[2] + vec1[1] * dv[2] };

			normalize(cvec1);
			printf("ck[0]: %lf cvec1: %lf %lf %lf\n", ck[0], cvec1[0], cvec1[1], cvec1[2] );

			if( ck[0] < 0 && fabs(cvec1[2]) > 0.5 )
				sRec->gather_flip = 0;
			else if( ck[0] > 0 && fabs(cvec1[2]) > 0.5 )
				sRec->gather_flip = 1;
			else if( ck[0] < 0 && fabs(cvec1[2]) <= 0.5 )
				sRec->gather_flip = 1;
			else if( ck[0] > 0 && fabs(cvec1[2]) <= 0.5 )
				sRec->gather_flip = 0;
			 
			if( sRec->gather_flip )
				printf("Flipping.\n");
			else
				printf("Not flipping.\n");

			// printf("Getting fusion pore info.

			if( sRec->gather_flip )
			{
				double uvt[2] = { uv_max_pos[0], uv_max_pos[1] };
				int ft = f_max_pos;

				uv_max_pos[0] = uv_max_neg[0];
				uv_max_pos[1] = uv_max_neg[1];
				f_max_pos = f_max_neg;

				uv_max_neg[0] = uvt[0];
				uv_max_neg[1] = uvt[1];
				f_max_neg = ft;	
			}

                               
			max_c( sRec->theSurface, &f_max_k, uv_max_k+0, uv_max_k+1, 1000, rsurf, MAX_C_GAUSS );     
			max_c( sRec->theSurface, &f_max_pos, uv_max_pos+0, uv_max_pos+1, 1000, rsurf, MAX_C_POS );     
			max_c( sRec->theSurface, &f_max_neg, uv_max_neg+0, uv_max_neg+1, 1000, rsurf, MAX_C_NEG );     
   
			sRec->theSurface->c(f_max_k, uv_max_k[0], uv_max_k[1], rsurf, &k, vec1,vec2, ck+0, ck+1 );
			sRec->theSurface->c(f_max_pos, uv_max_pos[0], uv_max_pos[1], rsurf, &k, vec1,vec2,cpos+0, cpos+1 );
			sRec->theSurface->c(f_max_neg, uv_max_neg[0], uv_max_neg[1], rsurf, &k, vec1,vec2,cneg+0, cneg+1 );
			

			

			sRec->theSurface->evaluateRNRM( f_max_k, uv_max_k[0], uv_max_k[1], max_saddle, nrm_junk, rsurf ); 
			sRec->theSurface->evaluateRNRM( f_max_pos, uv_max_pos[0], uv_max_pos[1], max_pos, nrm_junk, rsurf ); 
			sRec->theSurface->evaluateRNRM( f_max_neg, uv_max_neg[0], uv_max_neg[1], max_neg, nrm_junk, rsurf ); 

			printf("Max saddle at %lf %lf %lf c: %lf %lf\n", max_saddle[0], max_saddle[1], max_saddle[2] , ck[0], ck[1] ); 
			printf("Max pos at %lf %lf %lf c: %lf %lf\n", max_pos[0], max_pos[1], max_pos[2], cpos[0], cpos[1] ); 
			printf("Max neg at %lf %lf %lf c: %lf %lf\n", max_neg[0], max_neg[1], max_neg[2], cneg[0], cneg[1]  ); 
		
			double cen[3]={0,0,0};

			
			int nrim = 50;
			int nrim_out = 50;
			int f_ring[nrim];
			double uv_ring[nrim*2];
			double rim_pts[nrim*3];
			sRec->theSurface->get_cut_points( 2, max_pos[2], f_ring, uv_ring, rim_pts, nrim, rsurf, &nrim_out, cen, 0, 2  );
			cen[0]=0;
			cen[1]=0;
			cen[2]=0;
			double w = 0;
			for( int i = 0; i < nrim_out; i++ )
			{
				int ip1 = i+1;
				if( ip1 >= nrim_out ) ip1--;
				double dr[3] = { rim_pts[3*ip1+0] - rim_pts[3*i+0], rim_pts[3*ip1+1] - rim_pts[3*i+1], rim_pts[3*ip1+2] - rim_pts[3*i+2] };
				wrapPBC(dr,alpha);
				double l = normalize(dr); 
				//printf("C %lf %lf %lf\n", rim_pts[3*i+0], rim_pts[3*i+1], rim_pts[3*i+2] );
				cen[0] += rim_pts[3*i+0] * l;
				cen[1] += rim_pts[3*i+1] * l;
				cen[2] += rim_pts[3*i+2] * l;
				w += l;
			}

			cen[0]/=w;
			cen[1]/=w;
			cen[2]/=w;

			printf("the center: %lf %lf %lf\n", cen[0], cen[1], cen[2] );
			
			double av_r = 0;
			double av_r2 = 0;
			double com[3]={0,0,0};
			for( int i = 0; i < nrim_out; i++ )
			{
				com[0] += rim_pts[3*i+0] / nrim_out;
				com[1] += rim_pts[3*i+1] / nrim_out;
				com[2] += rim_pts[3*i+2] / nrim_out;
			} 

			printf("com: %lf %lf %lf\n", com[0], com[1], com[2] );	 
			
			for( int i = 0; i < nrim_out; i++ )
			{
				double dr[3] = { rim_pts[3*i+0] - com[0],
						 rim_pts[3*i+1] - com[1],
						 rim_pts[3*i+2] - com[2] };
				double r =normalize(dr);
				av_r += r;
				av_r2 += r*r;
			}
			av_r /= nrim_out;
			av_r2 /= nrim_out;
			double var = av_r2 -av_r*av_r;
	
			printf("r: %lf stdd %lf\n", av_r, sqrt(var) );
	
			fusion_pore_center[0] = com[0];
			fusion_pore_center[1] = com[1];
			fusion_pore_center[2] = com[2];
			fusion_pore_r = av_r;
			fusion_pore_h = 0;

			is_fusion_pore = 1;
		}
		else
		{	

			av_c /= sum_g;
			av_k /= sum_g;
			printf("Surface %d has <c> %lf <k> %lf ", sRec->id, av_c, av_k );
	
	
			if( av_c < 0 )
				printf("(flipping).\n");
			else
				printf("(not flipping).\n");
	
			if( av_c < 0 )
				sRec->gather_flip = 1;
			else
				sRec->gather_flip = 0;
		}
	}	

	if( is_fusion_pore )
	{
		center[0] = fusion_pore_center[0];
		center[1] = fusion_pore_center[1];
		center[2] = fusion_pore_center[2];
		*r = fusion_pore_r;
		*h = fusion_pore_h;

		return 1;
	}
	
	return 0;
}

void Simulation::gather( parameterBlock *block )
{
			double vec1[3], vec2[3];
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		double av_c = 0;
		double av_k = 0;
		double sum_g = 0;	
	
		double uv_max_neg[2];
		int f_max_neg;	
		double max_neg = 1e10;
	
		double uv_max_pos[2];
		int f_max_pos;	
		double max_pos = -1e10;
		
		double uv_max_k[2];
		int f_max_k;	
		double max_k = 1e10;

		// compute the energy to get sum K.

		sRec->theSurface->energy( sRec->r,NULL);

		double theK = SUMK;

		for( int f = 0; f < sRec->theSurface->nt; f++ )
		{
			double k;
			double g = sRec->theSurface->g(f, 1.0/3.0, 1.0/3.0, sRec->r ); 
			double c1,c2;
			double c = sRec->theSurface->c(f, 1.0/3.0, 1.0/3.0, sRec->r, &k, vec1, vec2, &c1, &c2 );
			av_c += g * c; 
			av_k += k * g;
			sum_g += g;

			if( c1 > max_pos )
			{
				max_pos = c1;
				uv_max_pos[0] = 1.0/3.0;
				uv_max_pos[1] = 1.0/3.0;
				f_max_pos = f;
			}
			else if( c2 > max_pos )
			{
				max_pos = c2;
				uv_max_pos[0] = 1.0/3.0;
				uv_max_pos[1] = 1.0/3.0;
				f_max_pos = f;
			}
			
			if( c1 < max_neg )
			{
				max_neg = c1;
				uv_max_neg[0] = 1.0/3.0;
				uv_max_neg[1] = 1.0/3.0;
				f_max_neg = f;
			}
			if( c2 < max_neg )
			{
				max_neg = c2;
				uv_max_neg[0] = 1.0/3.0;
				uv_max_neg[1] = 1.0/3.0;
				f_max_neg = f;
			}

			if( k < max_k )
			{
				max_k = c;
				uv_max_k[0] = 1.0/3.0;
				uv_max_k[1] = 1.0/3.0;
				f_max_k = f;
			}
		}

		double *rsurf = sRec->r;
		double chi_check = theK / (-4*M_PI);
	
		printf("K_tot: %lf\n", theK );	
		if( fabs(chi_check-1) < 0.1 )
	 	{
			printf("Fusion pore detected. Setting sign...");
			double rpt[3], rnrm[3];
			sRec->theSurface->evaluateRNRM( 0, 1.0/3.0, 1.0/3.0, rpt, rnrm, sRec->r );

			sRec->gather_flip = 0;

			if( fabs(rnrm[2]) > 0.7 )
			{  
				if( rnrm[2] < -0.7 && rpt[2] < 0 )
					sRec->gather_flip  = 1;
				else if( rnrm[2] > 0.7 && rpt[2] > 0 )
					sRec->gather_flip  = 1;
			}
			else 
			{
				double dp = rpt[0] * rnrm[0] + rpt[1] * rnrm[1] + rpt[2] * rnrm[2];

				if( dp < 0 )
					sRec->gather_flip = 1;
			}
			
			double nrm_junk[3];
			double max_saddle[3];
			double max_pos[3];
			double max_neg[3];
			double k;
			double cpos[2];
			double cneg[2];
			double ck[2];

			max_c( sRec->theSurface, &f_max_k, uv_max_k+0, uv_max_k+1, 1000, rsurf, MAX_C_GAUSS );     
			sRec->theSurface->c(f_max_k, uv_max_k[0], uv_max_k[1], rsurf, &k, vec1,vec2, ck+0, ck+1 );
		
			double du[3], dv[3];
			sRec->theSurface->ru( f_max_k, uv_max_k[0], uv_max_k[1], rsurf, du );	
			sRec->theSurface->rv( f_max_k, uv_max_k[0], uv_max_k[1], rsurf, dv );	

			double cvec1[3] = { vec1[0] * du[0] + vec1[1] * dv[0],
					    vec1[0] * du[1] + vec1[1] * dv[1], 
					    vec1[0] * du[2] + vec1[1] * dv[2] };

			normalize(cvec1);
			printf("ck[0]: %lf cvec1: %lf %lf %lf\n", ck[0], cvec1[0], cvec1[1], cvec1[2] );

			if( ck[0] < 0 && fabs(cvec1[2]) > 0.5 )
				sRec->gather_flip = 0;
			else if( ck[0] > 0 && fabs(cvec1[2]) > 0.5 )
				sRec->gather_flip = 1;
			else if( ck[0] < 0 && fabs(cvec1[2]) <= 0.5 )
				sRec->gather_flip = 1;
			else if( ck[0] > 0 && fabs(cvec1[2]) <= 0.5 )
				sRec->gather_flip = 0;
			 
			if( sRec->gather_flip )
				printf("Flipping.\n");
			else
				printf("Not flipping.\n");

			// printf("Getting fusion pore info.

			if( sRec->gather_flip )
			{
				double uvt[2] = { uv_max_pos[0], uv_max_pos[1] };
				int ft = f_max_pos;

				uv_max_pos[0] = uv_max_neg[0];
				uv_max_pos[1] = uv_max_neg[1];
				f_max_pos = f_max_neg;

				uv_max_neg[0] = uvt[0];
				uv_max_neg[1] = uvt[1];
				f_max_neg = ft;	
			}

                               
			max_c( sRec->theSurface, &f_max_k, uv_max_k+0, uv_max_k+1, 1000, rsurf, MAX_C_GAUSS );     
			max_c( sRec->theSurface, &f_max_pos, uv_max_pos+0, uv_max_pos+1, 1000, rsurf, MAX_C_POS );     
			max_c( sRec->theSurface, &f_max_neg, uv_max_neg+0, uv_max_neg+1, 1000, rsurf, MAX_C_NEG );     
   
			sRec->theSurface->c(f_max_k, uv_max_k[0], uv_max_k[1], rsurf, &k, vec1,vec2, ck+0, ck+1 );
			sRec->theSurface->c(f_max_pos, uv_max_pos[0], uv_max_pos[1], rsurf, &k, vec1,vec2,cpos+0, cpos+1 );
			sRec->theSurface->c(f_max_neg, uv_max_neg[0], uv_max_neg[1], rsurf, &k, vec1,vec2,cneg+0, cneg+1 );
			

			

			sRec->theSurface->evaluateRNRM( f_max_k, uv_max_k[0], uv_max_k[1], max_saddle, nrm_junk, rsurf ); 
			sRec->theSurface->evaluateRNRM( f_max_pos, uv_max_pos[0], uv_max_pos[1], max_pos, nrm_junk, rsurf ); 
			sRec->theSurface->evaluateRNRM( f_max_neg, uv_max_neg[0], uv_max_neg[1], max_neg, nrm_junk, rsurf ); 

			printf("Max saddle at %lf %lf %lf c: %lf %lf\n", max_saddle[0], max_saddle[1], max_saddle[2] , ck[0], ck[1] ); 
			printf("Max pos at %lf %lf %lf c: %lf %lf\n", max_pos[0], max_pos[1], max_pos[2], cpos[0], cpos[1] ); 
			printf("Max neg at %lf %lf %lf c: %lf %lf\n", max_neg[0], max_neg[1], max_neg[2], cneg[0], cneg[1]  ); 
		
			double cen[3]={0,0,0};

			
			int f_ring[20];
			double uv_ring[40];
			double rim_pts[60];
			int nrim = 20;
			int nrim_out = 20;
			sRec->theSurface->get_cut_points( 2, max_pos[2], f_ring, uv_ring, rim_pts, nrim, rsurf, &nrim_out, cen, 0, 0  );

			printf("center: %lf %lf %lf\n", cen[0], cen[1], cen[2] );
			
			double av_r = 0;
			double av_r2 = 0;
			double com[3]={0,0,0};
			for( int i = 0; i < nrim_out; i++ )
			{
				com[0] += rim_pts[3*i+0] / nrim_out;
				com[1] += rim_pts[3*i+1] / nrim_out;
				com[2] += rim_pts[3*i+2] / nrim_out;
			} 

			printf("com: %lf %lf %lf\n", com[0], com[1], com[2] );	 
			
			for( int i = 0; i < nrim_out; i++ )
			{
				double dr[3] = { rim_pts[3*i+0] - com[0],
						 rim_pts[3*i+1] - com[1],
						 rim_pts[3*i+2] - com[2] };
				double r =normalize(dr);
				av_r += r;
				av_r2 += r*r;
			}
			av_r /= nrim_out;
			av_r2 /= nrim_out;
			double var = av_r2 -av_r*av_r;
	
			printf("r: %lf stdd %lf\n", av_r, sqrt(var) );
	
			fusion_pore_center[0] = com[0];
			fusion_pore_center[1] = com[1];
			fusion_pore_center[2] = com[2];
			fusion_pore_r = av_r;
			fusion_pore_h = 0;

			is_fusion_pore = 1;
		}
		else
		{	

			av_c /= sum_g;
			av_k /= sum_g;
			printf("Surface %d has <c> %lf <k> %lf ", sRec->id, av_c, av_k );
	
	
			if( av_c < 0 )
				printf("(flipping).\n");
			else
				printf("(not flipping).\n");
	
			if( av_c < 0 )
				sRec->gather_flip = 1;
			else
				sRec->gather_flip = 0;
		}
	}	


        char *gathered = (char *)malloc( sizeof(char) * (1+strlen(block->jobName) + strlen("_gathered.xyz")) );
	sprintf(gathered, "%s_gathered.xyz", block->jobName );
	FILE *gatheredXYZ = fopen(gathered,"w");
	if( !block->structureName )
	{
		printf("Gathering requires a structure (e.g., PSF) file.\n");
		exit(1);
	}

	io_readStructureFile( block->structureName );	
	io_initialize_read( block->dcdName );

	char *gatherName = (char *)malloc( sizeof(char) * ( 1 + strlen(block->jobName) + strlen("_gather.txt") ) );

	sprintf(gatherName, "%s_gather.txt", block->jobName );

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
	for( int f = 0; f < nframes; f++ )
	{
		io_readFrame( block->dcdName );
		printf("Gathering from frame %d.\n", f );
		fflush(stdout);
		int ns = io_nNSAtoms();
		io_align();
		double Lx,Ly,Lz;
		io_get_PBC(&Lx,&Ly,&Lz);	
		printf("Gather PBC: %le %le %le\n", Lx, Ly, Lz );
 
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
			double c = nearCurvature( r,&cout,&kout,&dp_out,&dz_out,&leaflet_out,nrm );
	

			// get nrm from polarity.

			int new_leaflet_out = -1;		
			if( astart >= 0 && astop >= 0 )
			{
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

			}
			else
			{
			}

//			printf("leaflet_out: %d new_leaflet_out: %d\n", leaflet_out, new_leaflet_out );

			if( new_leaflet_out >= 0 ) 
			{
				if( new_leaflet_out != leaflet_out )
				{
					leaflet_out = new_leaflet_out;
					cout *= -1;
				}
			}

			if( f == 0 && assignmentFile )
				fprintf(assignmentFile, "%s %s %d %s\n", at[a].resname, at[a].segid, at[a].res, (leaflet_out == 0 ? "lower" : "upper" ) );

			fprintf(gatherFile, "%s %s %d %le %le %lf %lf %s r: %le %le %le\n", at[a].segid, at[a].resname, at[a].res, cout, kout, dp_out, dz_out, (leaflet_out == 0 ? "lower" : "upper" ), at[a].x, at[a].y, at[a].z );	
			fflush(gatherFile);
		}	
		
		if( ncomplex > 0 )
		{
			int cur_complex = 0;

			int a = 0;
			char pseg[256];

			pseg[0] = '\0';

			int nat = io_nAtoms();
			int seg_start = 0;
			
			for( int a = 0; a <= nat; a++ )
			{
				if( a == nat || strcasecmp(at[a].segid, pseg ) )
				{
					int seg_stop = a-1;

					// a new segment, is it a protein?
					if( strlen(pseg)>1 && pseg[0] == 'P' && pseg[1] >= '0' && pseg[1] <= '9' )
					{
						if( cur_complex  < ncomplex )
						{
							allComplexes[cur_complex]->get( this, at, seg_start, seg_stop, nat ); 
							allComplexes[cur_complex]->refresh(this);
							cur_complex++;
						}
					}
					
					if( !strcasecmp( pseg, "PROA" ) )
					{
						if( cur_complex  < ncomplex )
						{
							allComplexes[cur_complex]->get( this, at, seg_start, seg_stop, nat ); 
							allComplexes[cur_complex]->refresh(this);
							cur_complex++;
						}
					}

					seg_start = a;

					if( a < nat )
						strcpy( pseg, at[a].segid );
				}
					 
			}
		}
		if( f== 0 )
			fclose(assignmentFile);
		if( gatheredXYZ )		
			writeLimitingSurface(gatheredXYZ);
	}	

	if( gatheredXYZ )
	{
		fclose(gatheredXYZ);
		gatheredXYZ=NULL;
	}	
	
	printf("Finished with gather.\n" );
	fflush(stdout);
} 

