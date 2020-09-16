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

static double monomer_MW = 10*26492.14; //amu 

void ifitm3::init( double *r )
{
	// aqueous initialization.

	// syt7 model has six sites
	
	base_init();

	nsites = 1;
	nattach = 0;

	alloc();

	double vdw_r = 40.0;
	sigma[0] = vdw_r;

	mass[0] = monomer_MW*100;

	bound = 0;
	
	// initial geometry.

}

// initialize the BAR domain on the membrane.

void ifitm3::init( Simulation *theSimulation, surface *theSurface, double *rsurf, int f, double u, double v )
{
	// assume for now this is one of the points on the membrane neck.

	base_init();

	nsites = 1;
	nattach = 1;

	alloc();

	double vdw_r = 40.0;
	sigma[0] = vdw_r;

	mass[0] = monomer_MW*100;

	bound = 1;
	
	double rpt_attach1[3], nrm_attach1[3];
	theSurface->evaluateRNRM( f, u, v, rpt_attach1, nrm_attach1, rsurf );
	sid[0] = theSurface->surface_id;
	
	//
	// This is the neck-attachment point of the first protein.
	// goes into site-position "0"

	puv[0] = u;
	puv[1] = v;
	fs[0] = f;


	rall[0] = rpt_attach1[0];	
	rall[1] = rpt_attach1[1];	
	rall[2] = rpt_attach1[2];	
	


	memcpy( grad_fs, fs, sizeof(int) * nattach );
	memcpy( grad_puv, puv, sizeof(double) * nattach*2 );

	memset( PBC_ext, 0, sizeof(double) * 3*nsites );
	
	setrall(theSimulation);
}

// custom orient procedure.



void ifitm3::bind( int f, double u, double v)
{
	bound = 1;
}

void ifitm3::unbind( void )
{
	bound = 0;
}


void ifitm3::loadParams( parameterBlock *block )
{
}


int ifitm3::getNBonds( void )
{
	return 0;
}

void ifitm3::putBonds( int *bond_list )
{
}

double ifitm3::V( Simulation *theSimulation  )
{
	// could do something here with Gaussian curvature to test its hypothesis

	double *alphas = theSimulation->alpha;

	double r[3*nsites];
	double n[3*nsites];

	if( bound )
	{
		// evaluate the real-space coordinates and normals based on the membrane surface coordinates.
		for( int s = 0; s < nattach; s++ )
		{
			surface_record *sRec = theSimulation->fetch(sid[s]);
			surface *theSurface = sRec->theSurface;
			double *rsurf = sRec->r;
			int f_1 = fs[s], nf = fs[s];
			double uv1[2] = { 0.33, 0.33 };
			double duv1[2] = { puv[2*s+0]-uv1[0], puv[2*s+1]-uv1[1] };
				
			
			if( puv[2*s+0] <= 0 || puv[2*s+1] <= 0 || puv[2*s+0]+puv[2*s+1] >= 1 )
			{
				double ro[3];
				theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], ro, n+3*s, rsurf );  

				do {
					f_1 = nf;
					nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf ); 
				} while( nf != f_1 );
	
				uv1[0] += duv1[0];		
				uv1[1] += duv1[1];		
			
				grad_fs[s] = f_1;
				grad_puv[2*s+0] = uv1[0];
				grad_puv[2*s+1] = uv1[1];

				theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], r+3*s, n+3*s, rsurf );  

				double dr[3] = { r[3*s+0] - ro[0], r[3*s+1] - ro[1], r[3*s+2] - ro[2] };
				double del[3];
				MinImage3D( dr, theSurface->PBC_vec, del, rsurf+theSurface->nv*3 );

				r[3*s+0] = ro[0] + dr[0];
				r[3*s+1] = ro[1] + dr[1];
				r[3*s+2] = ro[2] + dr[2];
			}
			else
			{
				theSurface->evaluateRNRM( fs[s], puv[2*s+0], puv[2*s+1], r+3*s, n+3*s, rsurf );
			
	
				grad_fs[s] = fs[s];
				grad_puv[2*s+0] = puv[2*s+0];
				grad_puv[2*s+1] = puv[2*s+1];
			
			}			
				
			r[3*s+0] += theSurface->PBC_vec[0][0] * PBC_ext[3*s+0] *alphas[0] + theSurface->PBC_vec[1][0] * PBC_ext[3*s+1] *alphas[0] + theSurface->PBC_vec[2][0] * PBC_ext[3*s+2] *alphas[0];
			r[3*s+1] += theSurface->PBC_vec[0][1] * PBC_ext[3*s+0] *alphas[1] + theSurface->PBC_vec[1][1] * PBC_ext[3*s+1] *alphas[1] + theSurface->PBC_vec[2][1] * PBC_ext[3*s+2] *alphas[1];
			r[3*s+2] += theSurface->PBC_vec[0][2] * PBC_ext[3*s+0] *alphas[2] + theSurface->PBC_vec[1][2] * PBC_ext[3*s+1] *alphas[2] + theSurface->PBC_vec[2][2] * PBC_ext[3*s+2] *alphas[2];

		}
		
		memcpy(r+nattach*3 , rall+nattach*3, sizeof(double) * 3*(nsites-nattach) );
	}
	else
	{
		memcpy(r , rall, sizeof(double) * 3*nsites );
		memset(n, 0, sizeof(double) * 3*nsites ); 
	}

	double pot = 0;

	return pot;
}

// gets derivative of internal energy relative to position (surfacer_g) and the normal (surfacen_g).

double ifitm3::grad( Simulation *theSimulation,  double *surfacer_g, double *surfacen_g )
{
	double *alphas = theSimulation->alpha;
	double r[3*nsites];
	double n[3*nsites];

	if( bound )
	{
		for( int s = 0; s < nattach; s++ )
		{
			surface_record *sRec = theSimulation->fetch(sid[s]);
			surface *theSurface = sRec->theSurface;
			double *rsurf = sRec->r;
			int f_1 = fs[s], nf = fs[s];
			double uv1[2] = { 0.33, 0.33 };
			double duv1[2] = { puv[2*s+0]-uv1[0], puv[2*s+1]-uv1[1] };
			
			double null_mom[2] = {0,0};

			coord_transform[4*s+0] = 1;
			coord_transform[4*s+1] = 0;
			coord_transform[4*s+2] = 0;
			coord_transform[4*s+3] = 1;

			if( puv[2*s+0] <= 0 || puv[2*s+1] <= 0 || puv[2*s+0]+puv[2*s+1] >= 1 )
			{
				double ro[3];
				theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], ro, n+3*s, rsurf );  

				do {
					f_1 = nf;
					nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf, null_mom, coord_transform+4*s  ); 
				} while( nf != f_1 );
	
				uv1[0] += duv1[0];		
				uv1[1] += duv1[1];		
			
				grad_fs[s] = f_1;
				grad_puv[2*s+0] = uv1[0];
				grad_puv[2*s+1] = uv1[1];

				theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], r+3*s, n+3*s, rsurf );  

				double dr[3] = { r[3*s+0] - ro[0], r[3*s+1] - ro[1], r[3*s+2] - ro[2] };
				double del[3];
				MinImage3D( dr, theSurface->PBC_vec, del, rsurf+theSurface->nv*3 );

				r[3*s+0] = ro[0] + dr[0];
				r[3*s+1] = ro[1] + dr[1];
				r[3*s+2] = ro[2] + dr[2];
			}
			else
			{
				grad_fs[s] = fs[s];
				grad_puv[2*s+0] = puv[2*s+0];
				grad_puv[2*s+1] = puv[2*s+1];
				theSurface->evaluateRNRM( fs[s], puv[2*s+0], puv[2*s+1], r+3*s, n+3*s, rsurf );
	
			}			
				
			r[3*s+0] += theSurface->PBC_vec[0][0] * PBC_ext[3*s+0] * alphas[0] + theSurface->PBC_vec[1][0] * PBC_ext[3*s+1] * alphas[0] + theSurface->PBC_vec[2][0] * PBC_ext[3*s+2] * alphas[0];
			r[3*s+1] += theSurface->PBC_vec[0][1] * PBC_ext[3*s+0] * alphas[1] + theSurface->PBC_vec[1][1] * PBC_ext[3*s+1] * alphas[1] + theSurface->PBC_vec[2][1] * PBC_ext[3*s+2] * alphas[1];
			r[3*s+2] += theSurface->PBC_vec[0][2] * PBC_ext[3*s+0] * alphas[2] + theSurface->PBC_vec[1][2] * PBC_ext[3*s+1] * alphas[2] + theSurface->PBC_vec[2][2] * PBC_ext[3*s+2] * alphas[2];

		}
		
		memcpy(r+nattach*3 , rall+nattach*3, sizeof(double) * 3*(nsites-nattach) );
	}

	double pot = 0;

	return pot;
}

void ifitm3::move_inside( void )
{
	pcomplex::move_inside();
}

void ifitm3::move_outside( void )
{
	pcomplex::move_outside();
}

void ifitm3::writeStructure( Simulation *theSimulation, 
		surface_mask *upperSurfaceMask, 
		surface_mask *lowerSurfaceMask, 
		struct atom_rec **at_out, int *nat_out, char **sequence, struct ion_add **ions, int *nions, aa_build_data *buildData )
{
	struct atom_rec *IFI = NULL;
	int nIFI=0;
	
	surface *theSurface;
	double *rsurf;
	theSimulation->fetch(sid[0],&theSurface,&rsurf);
	
	double io_sign = (is_inside ? -1 : 1 );

	int flipped = is_inside;
	int f_attach = fs[0];
	double u_attach = puv[2*0+0];
	double v_attach = puv[2*0+1];

	// get the direction of max neg curavture and align the helix along this?
	
	int align_on_max_neg = 1;
	
	double cvec1[2] = {0,0}, cvec2[2]={0,0};
	double c1=0,c2=0;
	double k;
	double ctot = theSurface->c(f_attach,u_attach,v_attach,rsurf,&k,cvec1,cvec2,&c1,&c2);

	if( io_sign < 0 )
	{
		ctot *= -1;
		c1 *= -1;
		c2 *= -1;
	}

	double align_vec[3];
	
	double drdu[3]={0,0,0}, drdv[3]={0,0,0};
	theSurface->ru( f_attach, u_attach, v_attach,   rsurf, drdu ); 
	theSurface->rv( f_attach, u_attach, v_attach,   rsurf, drdv ); 

	if( (c1 < c2 &&  align_on_max_neg) || (c1 > c2 && !align_on_max_neg) )
	{
		align_vec[0] = cvec1[0] * drdu[0] + cvec1[1] * drdv[0];
		align_vec[1] = cvec1[0] * drdu[1] + cvec1[1] * drdv[1];
		align_vec[2] = cvec1[0] * drdu[2] + cvec1[1] * drdv[2];
	}	
	else
	{
		align_vec[0] = cvec2[0] * drdu[0] + cvec2[1] * drdv[0];
		align_vec[1] = cvec2[0] * drdu[1] + cvec2[1] * drdv[1];
		align_vec[2] = cvec2[0] * drdu[2] + cvec2[1] * drdv[2];
	}	

	normalize(align_vec);

	// r_C2A_CBL determines the attachment point, nrm_* does a vector orientation.
	double r[3], nrm[3];
	theSurface->evaluateRNRM( f_attach, u_attach, v_attach, r, nrm, rsurf );

	nrm[0] *= io_sign;
	nrm[1] *= io_sign;
	nrm[2] *= io_sign;
	
	// structures loaded into the pool can be re-used and don't need to be freed here. cleanup should be through the pool mechanism.
	
	// chooses which residues to place at this site.
	
	// these are centered at the bilayer midplane.
	int pool_code = pdbFetch( &IFI, &nIFI, "ifitm3_helix", "ifitm3_helix", addToPool );

	const char *segid_search = "PROA";

	int grab = 0;
	for( int x = 0; x < nIFI; x++ )
	{
		if( !strcasecmp(IFI[x].segid, segid_search ) )
			grab++;		
	}
	
	// four extra coordinates for orientation purposes.
	
	double *pcopy = (double *)malloc( sizeof(double) * 3 * (grab+4) );

	grab = 0;
	
	int set_IFI[3] = {0,1};
	int set_IFI_pdb[3] = {-1,-1,-1};
	int res_IFI_pdb[3] = { 2, 17 };
	const char *at_IFI_pdb[3] = { "CA", "CA" };
	double src_xy[2], src_orientation[2];

	for( int x = 0; x < nIFI; x++ )
	{
		if( !strcasecmp(IFI[x].segid, segid_search ) )
		{
			pcopy[3*grab+0] = IFI[x].x * (flipped ? -1 : 1 );  
			pcopy[3*grab+1] = IFI[x].y; 
			pcopy[3*grab+2] = IFI[x].z * (flipped ? -1 : 1 ); 

			for( int t = 0; t < 2; t++ )
			{
				if( !strcasecmp( IFI[x].atname, at_IFI_pdb[t] ) && IFI[x].res == res_IFI_pdb[t] )
					set_IFI_pdb[t] = grab;
			}
			grab++;		
		}
	}
	
	
	if( set_IFI_pdb[0] == -1 || set_IFI_pdb[1] == -1 )
	{
		printf("ERROR could not grab correct atoms from library PDB for Syt7 IFI.\n");
		exit(1);
	}

	double cen[2] = { 
			(pcopy[3*set_IFI_pdb[0]+0]+pcopy[3*set_IFI_pdb[1]+0])/2,
			(pcopy[3*set_IFI_pdb[0]+1]+pcopy[3*set_IFI_pdb[1]+1])/2 };

	// sets the mid-plane position directly below the protein.
	pcopy[3*grab+0] = cen[0]; 	
	pcopy[3*grab+1] = cen[1]; 	
	pcopy[3*grab+2] = 0; 	
	
	// one angstrom above, the normal	
	pcopy[3*(grab+1)+0] = cen[0]; 
	pcopy[3*(grab+1)+1] = cen[1];
	pcopy[3*(grab+1)+2] = (flipped ? -1 : 1) * 1; 


	double d_mod[3] = { pcopy[3*set_IFI_pdb[1]+0] - pcopy[3*set_IFI_pdb[0]+0],
			    pcopy[3*set_IFI_pdb[1]+1] - pcopy[3*set_IFI_pdb[0]+1],
			    0 };

	normalize(d_mod);

	// the alignment direction (with align_vec)
	pcopy[3*(grab+2)+0] = cen[0] + d_mod[0]; 	
	pcopy[3*(grab+2)+1] = cen[1] + d_mod[1]; 	
	pcopy[3*(grab+2)+2] = 0; // midplane	
	
	// the local x-direction.
	pcopy[3*(grab+3)+0] = cen[0] + 1.0;	
	pcopy[3*(grab+3)+1] = cen[1]; 	
	pcopy[3*(grab+3)+2] = 0; // midplane	
	
	// the center that gets mapped to the attachment point.
	src_xy[0] = cen[0]; 
	src_xy[1] = cen[1]; 


	double align[9] = {
		r[0], r[1], r[2], // attachment point, bilayer midplane.
		r[0]+nrm[0], r[1]+nrm[1], r[2]+nrm[2], // attachment point plus normal.
		r[0] + align_vec[0], r[1] + align_vec[1], r[2] + align_vec[2], // align vector.
	};

	int main_align_set[3] = { grab, grab+1, grab+2 };
	int surf_align_set[3] = { 0, 1, 2 };

	alignStructuresOnAtomSet( align, surf_align_set, pcopy, main_align_set, 3, grab+4 ); 

	double x_dir[3] = { 
			pcopy[(grab+3)*3+0] - pcopy[(grab)*3+0],
			pcopy[(grab+3)*3+1] - pcopy[(grab)*3+1],
			pcopy[(grab+3)*3+2] - pcopy[(grab)*3+2] };

	double rr = x_dir[0]*x_dir[0]+x_dir[1]*x_dir[1]+x_dir[2]*x_dir[2];
	double ruru = drdu[0]*drdu[0]+drdu[1]*drdu[1]+drdu[2]*drdu[2];
	double rvrv = drdv[0]*drdv[0]+drdv[1]*drdv[1]+drdv[2]*drdv[2];
	double rurv = drdu[0]*drdv[0]+drdu[1]*drdv[1]+drdu[2]*drdv[2];
	double rru = drdu[0]*x_dir[0]+drdu[1]*x_dir[1]+drdu[2]*x_dir[2];
	double rrv = x_dir[0]*drdv[0]+x_dir[1]*drdv[1]+x_dir[2]*drdv[2];
	
	// this is the du,dv, vector that should be matched to the x-dir of the simulation
	src_orientation[0] = -(-rrv * rurv + 2 * rru * rvrv)/(rurv*rurv-4*ruru*rvrv);
	src_orientation[1] = -(-rru * rurv + 2 * rrv * ruru)/(rurv*rurv-4*ruru*rvrv);	
	
	// x_dir = ku drdu + kv drdv
	(*at_out) = (struct atom_rec *)realloc( *at_out, sizeof(struct atom_rec) * (*nat_out + grab ) );
	
	int a_start = *nat_out;

	int naddSpace = *nions+10;

	*ions = (ion_add *)realloc( *ions,  sizeof(ion_add) * naddSpace );

	int tx = 0;
	for( int a = 0; a < nIFI; a++ )
	{
		if( !strcasecmp(IFI[a].segid, segid_search ) )
		{
			if( !strncasecmp( IFI[a].atname, "CA",2 ) && fabs(IFI[a].charge-2) < 1e-4 )
			{
				if( *nions == naddSpace )
				{
					naddSpace *= 2;
					*ions = (ion_add *)realloc( *ions, sizeof(ion_add) * naddSpace );
				}

				(*ions)[*nions].type = ION_CAL;
				(*ions)[*nions].x = pcopy[3*tx+0];
				(*ions)[*nions].y = pcopy[3*tx+1];
				(*ions)[*nions].z = pcopy[3*tx+2];
		
				(*nions) += 1;
			}
			else
			{
				(*at_out)[*nat_out] = IFI[a];
				// the oriented structure.
				(*at_out)[*nat_out].x = pcopy[3*tx+0];
				(*at_out)[*nat_out].y = pcopy[3*tx+1];
				(*at_out)[*nat_out].z = pcopy[3*tx+2];
				(*nat_out) += 1;
			}

			tx++;
		}
	}

	free(pcopy);	

	// write sequence

	struct atom_rec *at = *at_out;
	int pres = at[a_start].res;
	
	int seq_space = 10;
	int nseq = 0;

	char *seq = (char *)malloc( sizeof(char) * seq_space );

	seq[nseq] = threeToOne( at[a_start].resname );
	nseq++;

	for( int a = a_start; a < *nat_out; a++ )
	{
		if( at[a].res != pres )
		{
			int doff = at[a].res - pres;

			if( nseq-1 + doff >= seq_space )
			{
				seq_space *= 2;
				seq_space += doff;	
			

				seq = (char *)realloc( seq, sizeof(char) * (seq_space+1) );
			}

			for( int t = nseq; t < nseq-1+doff; t++ )
				seq[t] = 'X';		

			nseq = nseq-1+doff;

			seq[nseq] = threeToOne( at[a].resname ); 
			nseq++;
		}

		pres = at[a].res;
	}	

	seq[nseq] = '\0';
	*sequence = seq;
	

	if( io_sign > 0 )
		upperSurfaceMask->modifyMaskWithPoolAtPoint( theSurface, rsurf, theSimulation->PBC_vec, pool_code,
				f_attach, u_attach, v_attach, 20.0, 
					src_xy,
					src_orientation, flipped );				 
	else
		lowerSurfaceMask->modifyMaskWithPoolAtPoint( theSurface, rsurf, theSimulation->PBC_vec, pool_code,
				f_attach, u_attach, v_attach, 20.0, 
					src_xy,
					src_orientation, flipped );				 
	
/*
	printf("SEQUENCE:");
	for( int t = 0; t < nseq; t++ )
		printf("%c", seq[t] );
	printf("\n");
*/
}
