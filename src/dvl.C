#include "interp.h"
#include "pcomplex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mutil.h"
#include "gsl_random_globals.h"
#include "alignSet.h"
#include "pdbFetch.h"
#include "face_mask.h"
#include "aa_build_util.h"
#include "uv_map.h"
#include "M_matrix.h"
#include "lapack_we_use.h"
#include "gather.h"
#include "io_mol_read.h"

void DVL::get( 
		Simulation *theSimulation, // this is the surface/PBC
		struct atom_rec *at,
		int dvl_start, int dvl_stop, int nat_tot )
{
	int nsegments = 2;
	
	surface *theSurface;
	double *rsurf;
	theSimulation->fetch(sid[0],&theSurface,&rsurf);
	surface_record *sRec = theSimulation->fetch(sid[0]);

	// angle data for creating a histogram.	

	getLocalLipidComposition( theSimulation, at, dvl_start, dvl_stop, nat_tot );
	

	// CLASSIFICATION

	// * C2A near fusion pore
	
	// * C2B near fusion pore

	// * Both near fusion pore
	
	// Both in bulk region
}

void DVL::getLocalLipidComposition( 
		Simulation *theSimulation, // this is the surface/PBC
		struct atom_rec *at,
		int dvl_start, int dvl_stop, int nat_tot )
{
	// for near pointing
	double **M;
	int mlow,mhigh;

	getM( &M, &mlow, &mhigh );
	int nsegments = 2;
	
	surface *theSurface;
	double *rsurf;
	theSimulation->fetch(sid[0],&theSurface,&rsurf);
	surface_record *sRec = theSimulation->fetch(sid[0]);

	int flip_sign = sRec->gather_flip;
	double Ls[3] = {-1,-1,-1};

	/*
 *		For each ``face'' of the protein that we're interested in, get the distance from the membrane and establish a local coordinate system based on our orientation. Then, output lipid centers-of-mass in the lipid coordinate system.
 *
 * 	the distance will be defined by a point-of-contact, and the orientation by the vector from the point-of-contact to the center of mass of the domain.
 *
 * 	*/

	int ndomains = 1;

	// STEP 1: get domain COMs	
	
	double alphas[3] = {1,1,1};
	double domain_com[1][3] = { {0,0,0} };
	double ndom[1] = {0};

        int domain_start[1] = { 261 };
        int domain_stop[1] = { 355 };

	
	for( int domain = 0; domain < ndomains; domain++ )
	{
		for( int a = dvl_start; a < dvl_stop ; a++ )
		{
			if( at[a].res >= domain_start[domain] && at[a].res <= domain_stop[domain] )
			{
				double cur_com[3] = {0,0,0};
		
				if( ndom[domain] > 0 )
				{
					cur_com[0] = domain_com[domain][0] / ndom[domain];
					cur_com[1] = domain_com[domain][1] / ndom[domain];
					cur_com[2] = domain_com[domain][2] / ndom[domain];
				}
		
				
				double dr[3] = { at[a].x - cur_com[0],
						 at[a].y - cur_com[1],
						 at[a].z - cur_com[2] };
					
				theSimulation->wrapPBC( dr, alphas );
				
				domain_com[domain][0] += cur_com[0] + dr[0];
				domain_com[domain][1] += cur_com[1] + dr[1];
				domain_com[domain][2] += cur_com[2] + dr[2];
				ndom[domain] += 1;	
			}
		}
	}

	for( int d = 0; d < ndomains; d++ )
	{	
		domain_com[d][0] /= ndom[d];
		domain_com[d][1] /= ndom[d];
		domain_com[d][2] /= ndom[d];
	}

	printf("DOMAIN COM %lf %lf %lf\n", domain_com[0][0], domain_com[0][1], domain_com[0][2] );

	// STEP 2: get all lipid COMs.

	// neutral surface atoms.

	int ns = io_nNSAtoms();
	struct atom_rec *at_lipid = io_getCurrentFrame();

	double *lcom = (double *)malloc( sizeof(double) * 3 * ns );
	for( int ax = 0; ax < ns; ax ++ )
	{
		int a = io_getNSAtom(ax);
			
		int astart=-1;
		int astop=-1;
		io_getNSAtomStartStop(ax,&astart,&astop); 

		lcom[3*ax+0]=0;
		lcom[3*ax+1]=0;
		lcom[3*ax+2]=0;
			
		for( int x = 0; x < astop-astart; x++ )
		{
			lcom[3*ax+0] += at[astart+x].x;
			lcom[3*ax+1] += at[astart+x].y;
			lcom[3*ax+2] += at[astart+x].z;
		}

		lcom[3*ax+0]/=(astop-astart);
		lcom[3*ax+1]/=(astop-astart);
		lcom[3*ax+2]/=(astop-astart);
	}
	
	// STEP 3: get points-of-contact.

	// number of points-of-contact:
	int nPofCon[1] = { 5 };
	// # of residues in PofCon:
	int nRes[1][5] = {
		{ 1, 1, 1, 1, 1 }
	};

	int pocRes[1][5][1] =
	{
		{
			{290},
			{310},
			{274},
			{301},
			{335}
		}
	};
	
	int nfaces = 0;
	for( int domain = 0; domain < ndomains; domain++ )
		nfaces += nPofCon[domain];

	double pcom[1][5][3];	
	int f_ind_for_d[1][5];
	int domain_for_face[nfaces];
	int face_index[nfaces];
	double orientation[nfaces*3];
	double poc[nfaces*3];
	double near_point[nfaces*3];
	double tangent_plane[nfaces*9];		
	double dist[nfaces];

	int f_ind = 0;
	for( int domain = 0; domain < ndomains; domain++ )
	{
		for( int poc = 0; poc < nPofCon[domain]; poc++ )
		{
			domain_for_face[f_ind] = domain;
			face_index[f_ind] = poc;
			

			pcom[domain][poc][0] = 0;
			pcom[domain][poc][1] = 0;
			pcom[domain][poc][2] = 0;
			
			f_ind_for_d[domain][poc] = f_ind;
	
			f_ind++;
		}
	}

	for( int a = dvl_start; a < dvl_stop ; a++ )
	{
		int f_ind = 0;
		for( int domain = 0; domain < ndomains; domain++ )
		{
			for( int poc = 0; poc < nPofCon[domain]; poc++ )
			{
				for( int r = 0; r < nRes[domain][poc]; r++ )
				{
					if( at[a].res == pocRes[domain][poc][r] && !strcasecmp(at[a].atname, "CA") )
					{
						pcom[domain][poc][0] += at[a].x;
						pcom[domain][poc][1] += at[a].y;
						pcom[domain][poc][2] += at[a].z;
					} 
				}
				f_ind++;
			}
		}
	}

	f_ind=0;	
	for( int domain = 0; domain < ndomains; domain++ )
	{
		for( int tpoc = 0; tpoc < nPofCon[domain]; tpoc++ )
		{
			pcom[domain][tpoc][0] /= nRes[domain][tpoc];
			pcom[domain][tpoc][1] /= nRes[domain][tpoc];
			pcom[domain][tpoc][2] /= nRes[domain][tpoc];

			poc[f_ind*3+0] = pcom[domain][tpoc][0];	
			poc[f_ind*3+1] = pcom[domain][tpoc][1];	
			poc[f_ind*3+2] = pcom[domain][tpoc][2];	

			f_ind++;
		}
	}

	// Step 4: get each face's orientation and tangent plane coordinate system.
	
	int f_ref_CBL[1][5][2] =
	{
		{	{0,1}, {0,0}, {0, 3}, {0,2}, {0,0} } // domain 0, either 1 or 0.
	};
	
	f_ind=0;	
	for( int domain = 0; domain < ndomains; domain++ )
	{
		for( int tpoc = 0; tpoc < nPofCon[domain]; tpoc++ )
		{
			// pof - domain com
			orientation[f_ind*3+0] = poc[f_ind*3+0] - domain_com[domain][0]; 
			orientation[f_ind*3+1] = poc[f_ind*3+1] - domain_com[domain][1]; 
			orientation[f_ind*3+2] = poc[f_ind*3+2] - domain_com[domain][2]; 

			// normalized.
			normalize(orientation+f_ind*3);

			// near point on surface:
	
			// find near points for the attachment sites.
			double the_dist;
			int f;
			double col_u, col_v;
			double pt[3] = { poc[3*f_ind+0], poc[3*f_ind+1], poc[3*f_ind+2] };
			theSurface->nearPointOnBoxedSurface( pt, &f, &col_u, &col_v, M, mlow, mhigh, &the_dist );					
			double rpt[3], npt[3];
			theSurface->evaluateRNRM( f, col_u, col_v, rpt, npt, rsurf );

			double dr[3] = { 
				rpt[0] - poc[3*f_ind+0],
				rpt[1] - poc[3*f_ind+1],
				rpt[2] - poc[3*f_ind+2] };

			theSimulation->wrapPBC(dr,alphas);

			// near point , wrapped.
			near_point[3*f_ind+0] = poc[3*f_ind+0] + dr[0];
			near_point[3*f_ind+1] = poc[3*f_ind+1] + dr[1];
			near_point[3*f_ind+2] = poc[3*f_ind+2] + dr[2];

#define CBL_ORIENTATION

#ifdef CBL_ORIENTATION 
			int np_f = f_ind_for_d[f_ref_CBL[domain][tpoc][0]][f_ref_CBL[domain][tpoc][1]];
			double xvec[3] = { poc[3*f_ind+0] - poc[3*np_f+0],
					   poc[3*f_ind+1] - poc[3*np_f+1],
					   poc[3*f_ind+2] - poc[3*np_f+2] };
#else

			double xvec[3] = { orientation[f_ind*3+0], orientation[f_ind*3+1], orientation[3*f_ind+2] };
#endif
			double dp = dot( xvec, npt );

			xvec[0] -= dp * npt[0];
			xvec[1] -= dp * npt[1];
			xvec[2] -= dp * npt[2];

			normalize(xvec);

			double yvec[3];
			cross( xvec, npt, yvec );

			normalize(yvec );

			double check[3];

			cross( xvec, yvec, check );

			dp = dot( check, dr );
	
			// dr is near point - poc (points toward surface away from protein )
			// if x cross y dotted into dr is positive, the orientation is wrong (y points "down" if x to the right is positive.)

			if( dp > 0 )
			{
				yvec[0] *= -1;
				yvec[1] *= -1;
				yvec[2] *= -1;
			
			} 

			dp = dot( dr, npt );

			if( dp > 0 )
			{
				// make normal point to the protein.
				//
				npt[0] *= -1;
				npt[1] *= -1;
				npt[1] *= -1;
			}

			normalize(xvec);
			normalize(yvec);

			tangent_plane[9*f_ind+0] = xvec[0];
			tangent_plane[9*f_ind+1] = xvec[1];
			tangent_plane[9*f_ind+2] = xvec[2];

			tangent_plane[9*f_ind+3] = yvec[0];
			tangent_plane[9*f_ind+4] = yvec[1];
			tangent_plane[9*f_ind+5] = yvec[2];

			tangent_plane[9*f_ind+6] = npt[0];
			tangent_plane[9*f_ind+7] = npt[1];
			tangent_plane[9*f_ind+8] = npt[2];

			// orientation is POC-Domain_COM, 
			// dr is          pt-poc  
			double dp_dist_check = dot( orientation + 3*f_ind, dr );
	
			if( dp_dist_check < 0 )
				dist[f_ind] = -the_dist; 
			else
				dist[f_ind] = the_dist; 

			f_ind++;
		}
	}
	
	double max_d = 30.0; // 2*max_d by 2*max_d grid
	double max_dx = 60.0; // 2*max_d by 2*max_d grid

	// Step 5: print the positions of nearby lipids in the tangent-plane coordinate system.

	for( int f1 = 0; f1 < nfaces; f1++ )
	{
		int d1 = domain_for_face[f1];
		int p1 = face_index[f1];

		for( int f2 = 0; f2 < nfaces; f2++ )
		{
			int d2 = domain_for_face[f2];
			int p2 = face_index[f2];
			
			double dr[3] = { near_point[f2*3+0] - near_point[3*f1+0],
					 near_point[f2*3+1] - near_point[3*f1+1],
					 near_point[f2*3+2] - near_point[3*f1+2] };
			theSimulation->wrapPBC(dr,alphas);	
			
			double x_tangent_plane = dot( dr, tangent_plane+9*f1 );
			double y_tangent_plane = dot( dr, tangent_plane+9*f1 + 3 );
			double z_tangent_plane = dot( dr, tangent_plane+9*f1 + 6 ); 

			if( z_tangent_plane < 0 )
			{
				// opposite leaflet..
				continue;
			}
					
			if( z_tangent_plane > 30 )
			{
				// probably on the other side of the protein?
				continue;
			}

			if( fabs(x_tangent_plane) > max_dx ) continue;	
			if( fabs(y_tangent_plane) > max_d ) continue;	


			double dist_to_CBL1, dist_to_CBL2;

			int f_CBL1 = f_ind_for_d[0][0]; 
			int f_CBL2 = f_ind_for_d[1][0]; 
			double dr_check[3] = {  near_point[3*f1+0] - near_point[3*f_CBL1+0],
					        near_point[3*f1+1] - near_point[3*f_CBL1+1],
						near_point[3*f1+2] - near_point[3*f_CBL1+2] };
			theSimulation->wrapPBC( dr_check, alphas );
			
			double r_check1 = normalize(dr_check);
			dr_check[0] = near_point[3*f1+0] - near_point[3*f_CBL2+0];
			dr_check[1] = near_point[3*f1+1] - near_point[3*f_CBL2+1];
			dr_check[2] = near_point[3*f1+2] - near_point[3*f_CBL2+2];
			theSimulation->wrapPBC( dr_check, alphas );
			double r_check2 = normalize(dr_check);

			printf("DVL-DVL %d domain %d POC %d dist %lf domain2 %d POC %d dist2 %lf x %lf y %lf z %lf NP: %lf %lf %lf rCBL1: %lf rCBL2: %lf\n", my_id, d1, p1, dist[f1], d2, p2, dist[f2], x_tangent_plane, y_tangent_plane, z_tangent_plane, near_point[3*f1+0], near_point[3*f1+1], near_point[3*f1+2], r_check1, r_check2  ); 
		}
	}


	for( int l = 0; l < ns; l++ )
	{
		int a = io_getNSAtom(l);

		for( int f = 0; f < nfaces; f++ )
		{
			int domain = domain_for_face[f];
			int p = face_index[f];

			// lipid minus near point.
			double dr_l[3] = { lcom[3*l+0] - near_point[3*f+0],
					   lcom[3*l+1] - near_point[3*f+1],
					   lcom[3*l+2] - near_point[3*f+2] };
			theSimulation->wrapPBC(dr_l,alphas);	

			double x_tangent_plane = dot( dr_l, tangent_plane+9*f );
			double y_tangent_plane = dot( dr_l, tangent_plane+9*f + 3 );
			double z_tangent_plane = dot( dr_l, tangent_plane+9*f + 6 ); 

			if( z_tangent_plane < 0 )
			{
				// opposite leaflet..
				continue;
			}
					
			if( z_tangent_plane > 30 )
			{
				// probably on the other side of the protein?
				continue;
			}

			if( fabs(x_tangent_plane) > max_dx ) continue;	
			if( fabs(y_tangent_plane) > max_d ) continue;	


			double dist_to_CBL1, dist_to_CBL2;

			int f_CBL1 = f_ind_for_d[0][0]; 
			int f_CBL2 = f_ind_for_d[1][0]; 
			double dr_check[3] = {  near_point[3*f+0] - near_point[3*f_CBL1+0],
					        near_point[3*f+1] - near_point[3*f_CBL1+1],
						near_point[3*f+2] - near_point[3*f_CBL1+2] };
			theSimulation->wrapPBC( dr_check, alphas );
			
			double r_check1 = normalize(dr_check);
			dr_check[0] = near_point[3*f+0] - near_point[3*f_CBL2+0];
			dr_check[1] = near_point[3*f+1] - near_point[3*f_CBL2+1];
			dr_check[2] = near_point[3*f+2] - near_point[3*f_CBL2+2];
			theSimulation->wrapPBC( dr_check, alphas );
			double r_check2 = normalize(dr_check);

			printf("DVL %d domain %d POC %d dist %lf lipid %s segid %s res %d x %lf y %lf z %lf NP: %lf %lf %lf rCBL1: %lf rCBL2: %lf\n", my_id, domain, p, dist[f], at[a].resname, at[a].segid, at[a].res, x_tangent_plane, y_tangent_plane, z_tangent_plane, near_point[3*f+0], near_point[3*f+1], near_point[3*f+2], r_check1, r_check2  ); 
				 
		}
	}

	// Last: clean up.
	
	free(lcom);
}
