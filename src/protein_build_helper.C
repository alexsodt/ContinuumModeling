#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mutil.h"
#include "pcomplex.h"
#include "simulation.h"
#include "io_mol_read.h"
#include "face_mask.h"
#include "alignSet.h"
#include "pdb.h"
#include "aa_build_util.h"
#include "pdbFetch.h"

int pcomplex::addCurvatureOrientedPeripheralProteinHelper( Simulation *theSimulation, surface_mask *upperMask, surface_mask *lowerMask,
		struct atom_rec **at_out,
		int *nat_out,
		struct ion_add **ions,
		int *nions,
		int pool_code, // the pool code 
		const char *segid, // the segid of the protein to extract
		int orient_res[3], // the residues in the protein used to define the orientation.
		int attach_site, // the sites in the complex used to define the attachment.
		int align_on_neg,
		aa_build_data *buildData, // data structure to put build information (PSF-derived info and atom placements for detecting clashes). 
		double noise // for removing clashes
		)
{
	int attach_sites[2] = { attach_site, -1 };
	int aqueous_sites[2] = { -1, -1 };
	return addGeneralProteinHelper( theSimulation, upperMask, lowerMask,
		at_out, nat_out,
		ions, nions,
		pool_code,
		segid,
		orient_res,
		attach_sites,
		aqueous_sites,
		align_on_neg,
		ALIGN_TYPE_CURVATURE,
		buildData, noise);
}


int pcomplex::addPeripheralProteinHelper( Simulation *theSimulation, surface_mask *upperMask, surface_mask *lowerMask,
		struct atom_rec **at_out,
		int *nat_out,
		struct ion_add **ions,
		int *nions,
		int pool_code, // the pool code 
		const char *segid, // the segid of the protein to extract
		int orient_res[3], // the residues in the protein used to define the orientation.
		int attach_sites[2], // the sites in the complex used to define the attachment.
		int aqueous_sites[2], // the sites in the complex used to define the attachment.
		aa_build_data *buildData, // data structure to put build information (PSF-derived info and atom placements for detecting clashes). 
		double noise
		)
{
	return addGeneralProteinHelper( theSimulation, upperMask, lowerMask,
		at_out, nat_out,
		ions, nions,
		pool_code,
		segid,
		orient_res,
		attach_sites,
		aqueous_sites,
		0, // unused, curvature
		ALIGN_TYPE_SURF_AQ,
		buildData, noise );

}
;
int pcomplex::addAlignedProteinHelper( Simulation *theSimulation, surface_mask *upperMask, surface_mask *lowerMask,
		struct atom_rec **at_out,
		int *nat_out,
		struct ion_add **ions,
		int *nions,
		int pool_code, // the pool code 
		const char *segid, // the segid of the protein to extract
		int nalign,	// number of alignments.
		int *pdb_residues, // which residues
		int *complex_sites,	// which complex sites to align to those residues
		aa_build_data *buildData, // data structure to put build information (PSF-derived info and atom placements for detecting clashes). 
		double noise
		)
{
	// this routine should be more fundamental than the ``general'' one below ... combine them?
	//
	surface *theSurface;
	double *rsurf;
	theSimulation->fetch(sid[0],&theSurface,&rsurf);

	pool_structure *thePool = getPool(pool_code);

	int nat = thePool->nat;
	struct atom_rec *at = thePool->at;

	int grab = 0;


	for( int x = 0; x < nat; x++ )
	{
		if( !strcasecmp(at[x].segid, segid ) )
			grab++;		
	}
	
	double *pcopy = (double *)malloc( sizeof(double) * 3 * (grab+2) );

	grab = 0;
	
	int ngrab = nalign;
	int set_pdb[nalign];

	for( int t = 0; t < nalign; t++ )
		set_pdb[t] = -1;

	for( int x = 0; x < nat; x++ )
	{
		if( !strcasecmp(at[x].segid, segid ) )
		{
			pcopy[3*grab+0] = at[x].x;  
			pcopy[3*grab+1] = at[x].y;
			pcopy[3*grab+2] = at[x].z; 

			for( int t = 0; t < nalign; t++ )
			{
				if( (!strcasecmp( at[x].atname, "CA" ) || !strcasecmp( at[x].atname, "BB" )) && at[x].res == pdb_residues[t] )
					set_pdb[t] = grab;
			}
			grab++;		
		}
	}
	
	pcopy[grab*3+0] = 1;
	pcopy[grab*3+1] = 0;
	pcopy[grab*3+2] = 0;
	grab++;
	pcopy[grab*3+0] = 0;
	pcopy[grab*3+1] = 0;
	pcopy[grab*3+2] = 0;
	grab++;
	
	for( int t = 0; t < nalign; t++ )
	{
		if( set_pdb[t] == -1 )
		{
			printf("ERROR could not grab correct atoms for system '%s'.\n", thePool->fileNameStruct);
			exit(1);
		}
	}
	
	double align[3*nalign];
	int surf_align_set[nalign];
	
	for( int t = 0; t < nalign; t++ )
	{
		align[3*t+0] = rall[3*complex_sites[t]+0];
		align[3*t+1] = rall[3*complex_sites[t]+1];
		align[3*t+2] = rall[3*complex_sites[t]+2];
		surf_align_set[t] = t;
	} 


	alignStructuresOnAtomSet( align, surf_align_set, pcopy, set_pdb, nalign, grab ); 

	// CHECK FOR UNRECOVERABLE CLASHES HERE

	int clash =0;
	
	int *output_map = (int *)malloc( sizeof(int) *nat );
	double *flat_coords = (double *)malloc( sizeof(double) * 3 * nat );
	int tx = 0;
	int flat_map = 0;

	// first, get the map and coords.

	for( int a = 0; a < nat; a++ )
	{
		output_map[a] = -1;

		if( !strcasecmp(at[a].segid, segid ) )
		{
			if( strncasecmp( at[a].atname, "CAL",3 ) || fabs(at[a].charge-2) > 1e-4 )
			{
				// loop over its bonds.

				flat_coords[3*flat_map+0] = pcopy[3*tx+0];
				flat_coords[3*flat_map+1] = pcopy[3*tx+1];
				flat_coords[3*flat_map+2] = pcopy[3*tx+2];
				output_map[a] = flat_map;
				flat_map++;
				// the map needs to go from the pool'd index (a) to the output index (*nat_out).
			}

			tx++;
		}
	}
	
	clash = buildData->checkMappedCycles( flat_coords, output_map, nat, getPool(pool_code)->cycles, getPool(pool_code)->cycle_lengths, getPool(pool_code)->ncycles ); 
	if( ! clash )
		clash = buildData->checkMappedBonds( flat_coords, output_map, nat, getPool(pool_code)->bonds, getPool(pool_code)->bond_offsets, getPool(pool_code)->nbonds ); 

/*
	if( clash )
	{
		printf("CLASH! redoing.\n");
		free(output_map);
		free(flat_coords);
		free(pcopy);
		return 1;	
	}
*/
	// END CHECK.
	
	// the x dir check, unused now.
	grab -= 2;
	
	double x_dir[3] = { 
			pcopy[(grab+1)*3+0] - pcopy[(grab)*3+0],
			pcopy[(grab+1)*3+1] - pcopy[(grab)*3+1],
			pcopy[(grab+1)*3+2] - pcopy[(grab)*3+2] };

	(*at_out) = (struct atom_rec *)realloc( *at_out, sizeof(struct atom_rec) * (*nat_out + grab ) );
	
	int a_start = *nat_out;

	int naddSpace = *nions+10;

	*ions = (ion_add *)realloc( *ions,  sizeof(ion_add) * naddSpace );

	int place_offset = buildData->curPlace();

	flat_map = 0;
	tx = 0;

	for( int a = 0; a < nat; a++ )
	{
		output_map[a] = -1;

		if( !strcasecmp(at[a].segid, segid ) )
		{
			if( !strncasecmp( at[a].atname, "CAL",3 ) && fabs(at[a].charge-2) < 1e-4 )
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
				(*at_out)[*nat_out] = at[a];
				// the oriented structure.
				(*at_out)[*nat_out].x = pcopy[3*tx+0];
				(*at_out)[*nat_out].y = pcopy[3*tx+1];
				(*at_out)[*nat_out].z = pcopy[3*tx+2];

				// the map needs to go from the pool'd index (a) to the output index (*nat_out).

				flat_coords[3*flat_map+0] = pcopy[3*tx+0];
				flat_coords[3*flat_map+1] = pcopy[3*tx+1];
				flat_coords[3*flat_map+2] = pcopy[3*tx+2];
				output_map[a] = flat_map; // flat_map is the index relative to place_offset below, this is correct.
			
				buildData->addAtom( flat_coords+3*flat_map );
				
				flat_map++;
				(*nat_out) += 1;
			}

			tx++;
		}
	}


	buildData->addMappedCycles( place_offset, flat_coords, output_map, nat, getPool(pool_code)->cycles, getPool(pool_code)->cycle_lengths, getPool(pool_code)->ncycles ); 
	buildData->addMappedBonds( place_offset, output_map, nat, getPool(pool_code)->bonds, getPool(pool_code)->bond_offsets, getPool(pool_code)->nbonds ); 

	free(output_map);
	free(flat_coords);
	free(pcopy);	

	return 0;
}


int pcomplex::addGeneralProteinHelper( Simulation *theSimulation, surface_mask *upperMask, surface_mask *lowerMask,
		struct atom_rec **at_out,
		int *nat_out,
		struct ion_add **ions,
		int *nions,
		int pool_code, // the pool code 
		const char *segid, // the segid of the protein to extract
		int orient_res[3], // the residues in the protein used to define the orientation.
		int attach_sites[2], // the sites in the complex used to define the attachment.
		int aqueous_sites[2], // the sites in the complex used to define the attachment.
		int align_on_neg,
		int align_type,
		aa_build_data *buildData, // data structure to put build information (PSF-derived info and atom placements for detecting clashes). 
		double noise
		)
{
	surface *theSurface;
	double *rsurf;
	theSimulation->fetch(sid[0],&theSurface,&rsurf);

	double io_sign = (is_inside ? -1 : 1 ); // is_inside means put it opposite the surface normal instead of along it.
	
	int flipped = is_inside;
	
	int f_attach = fs[attach_sites[0]];
	double u_attach = puv[attach_sites[0]*2+0]; 
	double v_attach = puv[attach_sites[0]*2+1];

	int noise_done = 0;

	while( !noise_done )
	{
		double u_trial = u_attach;
		double v_trial = v_attach;

		double duv[2] = { noise * 2*(rand()/(double)RAND_MAX - 0.5),
				  noise * 2*(rand()/(double)RAND_MAX - 0.5) };
	
		int nf;
		do {
			nf = f_attach;
			f_attach = theSurface->nextFace( f_attach, &u_trial, &v_trial, duv+0, duv+1, rsurf );
		} while( nf != f_attach );	

		u_attach = u_trial;
		v_attach = v_trial;
	
		noise_done = 1;
	} 

	int f_orient;
	double u_orient,v_orient;
	
	if( align_type == ALIGN_TYPE_SURF_SURF )
	{
		f_orient = fs[attach_sites[1]];
		u_orient = puv[attach_sites[1]*2+0];
		v_orient = puv[attach_sites[1]*2+1];

		noise_done = 0;
		while( !noise_done )
		{
			double u_trial = u_orient;
			double v_trial = v_orient;
	
			double duv[2] = { noise * 2*(rand()/(double)RAND_MAX - 0.5),
					  noise * 2*(rand()/(double)RAND_MAX - 0.5) };
		
			int nf;
			do {
				nf = f_orient;
				f_orient = theSurface->nextFace( f_orient, &u_trial, &v_trial, duv+0, duv+1, rsurf );
			} while( nf != f_orient );	
	
			u_orient = u_trial;
			v_orient = v_trial;
		} 
	}
	
	double drdu[3]={0,0,0}, drdv[3]={0,0,0};
	theSurface->ru( f_attach, u_attach, v_attach,   rsurf, drdu ); 
	theSurface->rv( f_attach, u_attach, v_attach,   rsurf, drdv ); 
	
	double r_attach[3], nrm_attach[3];
	theSurface->evaluateRNRM( f_attach, u_attach, v_attach, r_attach, nrm_attach, rsurf );
	nrm_attach[0] *= io_sign; 
	nrm_attach[1] *= io_sign; 
	nrm_attach[2] *= io_sign; 
	
	double r_orient[3], nrm_orient[3];

	if( align_type == ALIGN_TYPE_SURF_SURF )
	{
		theSurface->evaluateRNRM( f_orient, u_orient, v_orient, r_orient, nrm_orient, rsurf );
		nrm_orient[0] *= io_sign; 
		nrm_orient[1] *= io_sign; 
		nrm_orient[2] *= io_sign; 
	}
	else if( align_type == ALIGN_TYPE_SURF_AQ )
	{
		r_orient[0] = rall[3*aqueous_sites[0]+0] + 20*noise * 2*rand()/(double)RAND_MAX;
		r_orient[1] = rall[3*aqueous_sites[0]+1] + 20*noise * 2*rand()/(double)RAND_MAX;
		r_orient[2] = rall[3*aqueous_sites[0]+2] + 20*noise * 2*rand()/(double)RAND_MAX;
	}	

	pool_structure *thePool = getPool(pool_code);

	int nat = thePool->nat;
	struct atom_rec *at = thePool->at;

	int grab = 0;

	double zcom = 0;

	for( int x = 0; x < nat; x++ )
	{
		if( !strcasecmp(at[x].segid, segid ) )
		{
			grab++;		
			zcom += at[x].z;
		}
	}

	zcom /= grab;
	
	if( zcom < 0 ) // if the peripheral protein happened to be simulated on the lower (-z) leaflet, rotate around y 180 degrees. we do the same thing to the lipids later (this info is added to the mask below).
		flipped = 1;
	else
		flipped = 0;
	
	double *pcopy = (double *)malloc( sizeof(double) * 3 * (grab+4) );

	grab = 0;
	
	const char *at_grab[] = { "CA", "CA", "CA" };
	
	int ngrab = 3;
	int set_pdb[]={-1,-1,-1};
	double src_xy[2], src_orientation[2];

	for( int x = 0; x < nat; x++ )
	{
		if( !strcasecmp(at[x].segid, segid ) )
		{
			pcopy[3*grab+0] = at[x].x * (flipped ? -1 : 1 );  
			pcopy[3*grab+1] = at[x].y; 
			pcopy[3*grab+2] = at[x].z * (flipped ? -1 : 1 ); 

			for( int t = 0; t < ngrab; t++ )
			{
				if( !strcasecmp( at[x].atname, at_grab[t] ) && at[x].res == orient_res[t] )
					set_pdb[t] = grab;
			}
			grab++;		
		}
	}
	
	if( set_pdb[0] == -1 || set_pdb[1] == -1 || (align_type == ALIGN_TYPE_SURF_AQ && set_pdb[2] == -1) )
	{
		printf("ERROR could not grab correct atoms for system '%s'.\n", thePool->fileNameStruct);
		exit(1);
	}
	
	double cen[2] = { 
			pcopy[3*set_pdb[0]+0],
			pcopy[3*set_pdb[0]+1] };

	// sets the mid-plane position directly below the protein.
	pcopy[3*grab+0] = cen[0]; 	
	pcopy[3*grab+1] = cen[1]; 	
	pcopy[3*grab+2] = 0; 	
	
	// one angstrom above, the normal	
	pcopy[3*(grab+1)+0] = cen[0]; 
	pcopy[3*(grab+1)+1] = cen[1];
	pcopy[3*(grab+1)+2] = (flipped ? -1 : 1) * 1; 


	double d_mod[3];

	if( align_type == ALIGN_TYPE_SURF_SURF || align_type == ALIGN_TYPE_CURVATURE )
	{
		d_mod[0] = pcopy[3*set_pdb[1]+0] - pcopy[3*set_pdb[0]+0];	
		d_mod[1] = pcopy[3*set_pdb[1]+1] - pcopy[3*set_pdb[0]+1];	
		d_mod[2] = 0;	
	}
	else if( align_type == ALIGN_TYPE_SURF_AQ )
	{
		d_mod[0] = pcopy[3*set_pdb[2]+0] - pcopy[3*set_pdb[1]+0];
		d_mod[1] = pcopy[3*set_pdb[2]+1] - pcopy[3*set_pdb[1]+1];
		d_mod[2] = pcopy[3*set_pdb[2]+2] - pcopy[3*set_pdb[1]+2];
	}
	else
	{
		printf("Align type error.\n");
		exit(1);
	}

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

	double align_vec[3];

	if( align_type == ALIGN_TYPE_SURF_SURF )
	{
		align_vec[0] = r_orient[0] - r_attach[0];
		align_vec[1] = r_orient[1] - r_attach[1];
		align_vec[2] = r_orient[2] - r_attach[2];
	}	
	else if( align_type == ALIGN_TYPE_SURF_AQ )
	{
		align_vec[0] = rall[3*aqueous_sites[1]+0] - rall[3*aqueous_sites[0]+0];
		align_vec[1] = rall[3*aqueous_sites[1]+1] - rall[3*aqueous_sites[0]+1];
		align_vec[2] = rall[3*aqueous_sites[1]+2] - rall[3*aqueous_sites[0]+2];
	}	
	else if( align_type == ALIGN_TYPE_CURVATURE )
	{
		double cvec1[3], cvec2[3];
		double c1, c2, k;

		double ctot = theSurface->c(f_attach,u_attach,v_attach,rsurf,&k,cvec1,cvec2,&c1,&c2);

		if( io_sign < 0 )
		{
			ctot *= -1;
			c1 *= -1;
			c2 *= -1;
		}

		if( (c1 < c2 &&  align_on_neg) || (c1 > c2 && !align_on_neg) )
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
	}
	else
	{
		printf("Unknown alignment type: %d.\n", align_type );
		exit(1);
	}

	normalize(align_vec);

	double align[9] = {
		r_attach[0], r_attach[1], r_attach[2], // attachment point, bilayer midplane.
		r_attach[0]+nrm_attach[0], r_attach[1]+nrm_attach[1], r_attach[2]+nrm_attach[2], // attachment point plus normal.
		r_attach[0] + align_vec[0], r_attach[1] + align_vec[1], r_attach[2] + align_vec[2], // align vector.
	};

	int main_align_set[3] = { grab, grab+1, grab+2 };
	int surf_align_set[3] = { 0, 1, 2 };

	alignStructuresOnAtomSet( align, surf_align_set, pcopy, main_align_set, 3, grab+4 ); 

	// CHECK FOR UNRECOVERABLE CLASHES HERE

	int clash =0;
	
	int *output_map = (int *)malloc( sizeof(int) *nat );
	double *flat_coords = (double *)malloc( sizeof(double) * 3 * nat );
	int tx = 0;
	int flat_map = 0;

	// first, get the map and coords.

	for( int a = 0; a < nat; a++ )
	{
		output_map[a] = -1;

		if( !strcasecmp(at[a].segid, segid ) )
		{
			if( strncasecmp( at[a].atname, "CA",2 ) || fabs(at[a].charge-2) > 1e-4 )
			{
				// loop over its bonds.

				flat_coords[3*flat_map+0] = pcopy[3*tx+0];
				flat_coords[3*flat_map+1] = pcopy[3*tx+1];
				flat_coords[3*flat_map+2] = pcopy[3*tx+2];
				output_map[a] = flat_map;
				flat_map++;
				// the map needs to go from the pool'd index (a) to the output index (*nat_out).
			}

			tx++;
		}
	}
	
	clash = buildData->checkMappedCycles( flat_coords, output_map, nat, getPool(pool_code)->cycles, getPool(pool_code)->cycle_lengths, getPool(pool_code)->ncycles ); 
	if( ! clash )
		clash = buildData->checkMappedBonds( flat_coords, output_map, nat, getPool(pool_code)->bonds, getPool(pool_code)->bond_offsets, getPool(pool_code)->nbonds ); 

	if( clash )
	{
		printf("CLASH! redoing.\n");
		free(output_map);
		free(flat_coords);
		free(pcopy);
		return 1;	
	}

	// END CHECK.

	
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
	
	(*at_out) = (struct atom_rec *)realloc( *at_out, sizeof(struct atom_rec) * (*nat_out + grab ) );
	
	int a_start = *nat_out;

	int naddSpace = *nions+10;

	*ions = (ion_add *)realloc( *ions,  sizeof(ion_add) * naddSpace );


	int place_offset = buildData->curPlace();

	flat_map = 0;
	tx = 0;

	for( int a = 0; a < nat; a++ )
	{
		output_map[a] = -1;

		if( !strcasecmp(at[a].segid, segid ) )
		{
			if( !strncasecmp( at[a].atname, "CA",2 ) && fabs(at[a].charge-2) < 1e-4 )
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
				(*at_out)[*nat_out] = at[a];
				// the oriented structure.
				(*at_out)[*nat_out].x = pcopy[3*tx+0];
				(*at_out)[*nat_out].y = pcopy[3*tx+1];
				(*at_out)[*nat_out].z = pcopy[3*tx+2];

				// the map needs to go from the pool'd index (a) to the output index (*nat_out).

				flat_coords[3*flat_map+0] = pcopy[3*tx+0];
				flat_coords[3*flat_map+1] = pcopy[3*tx+1];
				flat_coords[3*flat_map+2] = pcopy[3*tx+2];
				output_map[a] = flat_map; // flat_map is the index relative to place_offset below, this is correct.
			
				buildData->addAtom( flat_coords+3*flat_map );
				
				flat_map++;
				(*nat_out) += 1;
			}

			tx++;
		}
	}


	buildData->addMappedCycles( place_offset, flat_coords, output_map, nat, getPool(pool_code)->cycles, getPool(pool_code)->cycle_lengths, getPool(pool_code)->ncycles ); 
	buildData->addMappedBonds( place_offset, output_map, nat, getPool(pool_code)->bonds, getPool(pool_code)->bond_offsets, getPool(pool_code)->nbonds ); 

	free(output_map);
	free(flat_coords);
	free(pcopy);	

	if( io_sign > 0 )
		upperMask->modifyMaskWithPoolAtPoint( theSurface, rsurf, theSurface->PBC_vec, pool_code,
		         f_attach, u_attach, v_attach, 20.0, 
		         	src_xy,
		         	src_orientation, flipped );				 
	else
		lowerMask->modifyMaskWithPoolAtPoint( theSurface, rsurf, theSurface->PBC_vec, pool_code,
				f_attach, u_attach, v_attach, 20.0, 
					src_xy,
					src_orientation, flipped );				 
	return 0;
}

