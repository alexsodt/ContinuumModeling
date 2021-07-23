
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
#include "io_mol_read.h"
#include "globals.h"

void CNT::writeStructure( Simulation *theSimulation, surface_mask *upperSurfaceMask, surface_mask *lowerSurfaceMask, struct atom_rec **at_out, int *nat_out, char ***seq, int *nseq, int **seq_at_array,  char ***patches, ion_add **ions, int *nions, struct aa_build_data *buildData, int *add_type )
{
	// the CNT is a transmembrane pore oriented along the bilayer normal.

	*add_type = BUILD_GENERAL; // includes its own RTF.

	int nsegments = 2;
	
	surface *theSurface;
	double *rsurf;
	theSimulation->fetch(sid[0],&theSurface,&rsurf);

	struct atom_rec *at_cnt;
	int nat;

	int pool_code;

	if( doMartini() )
		pool_code  = pdbFetch( &at_cnt, &nat, "CNT", "CNT_martini", addToPool );
	else
		pool_code  = pdbFetch( &at_cnt, &nat, "CNT", "CNT_aa", addToPool );



	char *rtf;	

	if( doMartini() )
		rtfFetch( &rtf, (const char *)"CNT", (const char *)"CNT_martini", 1 ); 
	else
		rtfFetch( &rtf, (const char *)"CNT", (const char *)"CNT_aa", 1 ); 


	int len = strlen("PCOMPLEX_RTF") + 128 /* pcomplex # */ + strlen("str");

	char *fileName = (char *)malloc( sizeof(char) * (1+len) );
	sprintf(fileName, "PCOMPLEX_RTF_%d.str", my_id ); 
	FILE *rtf_local = fopen(fileName,"w");
	free(fileName);

	if( rtf )
	{	
		char *t = rtf;

		while( *t )
		{
			fputc( *t, rtf_local );
			t += 1;
		}

		free(rtf);
	}
	else
	{
		fprintf(rtf_local, "! empty on purpose, rtf already loaded.\n");
		fprintf(rtf_local, "\n");
	}
	fclose(rtf_local);	

	double pt[3];
	double nrm[3];

	theSurface->evaluateRNRM( fs[0], puv[0], puv[1], pt, nrm, rsurf );

	// put all the atoms in the CNT pdb into the atom set. it's aligned on Z.
	
	double *pcopy = (double *)malloc( sizeof(double) * 3 * nat );

	double cog[3] = {0,0,0};
	for( int a = 0; a < nat; a++ )
	{
		pcopy[3*a+0] = at_cnt[a].x;
		pcopy[3*a+1] = at_cnt[a].y;
		pcopy[3*a+2] = at_cnt[a].z;

		cog[0] += pcopy[3*a+0];
		cog[1] += pcopy[3*a+1];
		cog[2] += pcopy[3*a+2];
	}

	cog[0] /= nat;
	cog[1] /= nat;
	cog[2] /= nat;

	double avr=0;

	for( int a = 0; a < nat; a++ )
	{
		avr += sqrt( (pcopy[3*a+0] - cog[0])*(pcopy[3*a+0] - cog[0]) + (pcopy[3*a+1] - cog[1])*(pcopy[3*a+1] - cog[1]) );

		pcopy[3*a+0] += (pt[0] - cog[0]);	
		pcopy[3*a+1] += (pt[1] - cog[1]);	
		pcopy[3*a+2] += (pt[2] - cog[2]);	
	}

	avr /= nat;

	double z_axis[3] = { 0,0,1};
	double rotor_axis[3];

	cross( nrm, z_axis, rotor_axis );

	double dp = dot( nrm, z_axis );
	double theta = acos(dp);	

	if( theta > 1e-4 )
		rotateArbitrary( pcopy, rotor_axis, pt, nat, theta ); 
	
	(*at_out) = (struct atom_rec *)realloc( *at_out, sizeof(struct atom_rec) * (*nat_out + nat ) );

	int a_start = *nat_out;
	int place_offset = buildData->curPlace();

	double *flat_coords = (double *)malloc( sizeof(double) * 3 * nat );

	int flat_map = 0;
	int tx = 0;

	int *output_map = (int *)malloc( sizeof(int) * nat );

	for( int a = 0; a < nat; a++ )
	{
		output_map[a] = -1;

		(*at_out)[*nat_out] = at_cnt[a];
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
	
//		printf("C %lf %lf %lf\n", flat_coords[3*flat_map+0], flat_coords[3*flat_map+1], flat_coords[3*flat_map+2] );
	
		flat_map++;
		(*nat_out) += 1;

		tx++;
	}

	int *bonds = NULL;
	int nb = 0;
	autoBonds( at_cnt, nat, &bonds, &nb );
	
	int **basis, *basis_length, nbasis;
	fetchCycleBasis( &basis, &basis_length, &nbasis, bonds, nb, nat ); 
		
	int * bond_offsets = (int *)malloc( sizeof(int) * nat );
	int * nbonds = (int *)malloc( sizeof(int) * nat);
	memset( nbonds, 0, sizeof(int) * nat );
	int off = 0;		

	for( int b = 0; b < nb; b++ )
	{
		nbonds[bonds[2*b+0]] += 1;	
		nbonds[bonds[2*b+1]] += 1;	
	}

	for( int a = 0; a < curNAtoms(); a++ )
	{
		bond_offsets[a]=off;
		off += nbonds[a];
	}

	buildData->addMappedCycles( place_offset, flat_coords, output_map, nat, 
			basis, basis_length, nbasis );
	buildData->addMappedBonds( place_offset, output_map, nat, 
			bonds, bond_offsets, nbonds );	

	// add virtual sites inside the channel so that it kills anything that might accidentally be put inside.
	
	int nvirt = 200;
	double spread = 30.0;
	double noise_mag = avr/sqrt(2.0);

	for( int v = 0; v <= nvirt; v++ )
	{
		double displ[3] = { 
				noise_mag * 2*((double)rand()/(double)RAND_MAX-0.5),
				noise_mag * 2*((double)rand()/(double)RAND_MAX-0.5),
				noise_mag * 2*((double)rand()/(double)RAND_MAX-0.5) };
		
		double dp = dot(displ, nrm);
		
		displ[0] -= nrm[0] * dp;
		displ[1] -= nrm[1] * dp;
		displ[2] -= nrm[2] * dp;

		double p[3] = { pt[0] - nrm[0] * spread/2 + nrm[0] * spread * (v)/nvirt + displ[0], 
				pt[1] - nrm[1] * spread/2 + nrm[1] * spread * (v)/nvirt + displ[1],
				pt[2] - nrm[2] * spread/2 + nrm[2] * spread * (v)/nvirt + displ[2] };


		buildData->addAtom( p );
		
//		printf("O %lf %lf %lf\n", p[0], p[1], p[2]  );
	}	
	

	free(output_map);
	free(flat_coords);
	free( bond_offsets);
	free( nbonds);
	free( bonds);	
	free( basis_length );

	for( int i = 0; i < nbasis; i++ )
	{
		if( basis[i] )
			free(basis[i]);
	}

	*nseq = 0;

}
