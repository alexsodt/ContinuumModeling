#define __io_mol_readc__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
#include "library.h"
#include "io_mol_read.h"
#include "mutil.h"

// trajectory reading provided by either LOOS or some legacy functions I wrote -AJS
// header file is the same, regardless.

//void io_readStructureFile( char *fileName );
//void io_readFrame( char *fileName );

#include "pdb.h"

#ifdef HAVE_LIBLOOS
#include "loos.hpp"
#else
#include "dcd.h"
#endif

typedef struct
{
	FILE *file_handle;
	


	struct atom_rec *cur_frame;
	struct atom_rec *last_frame;
	int natoms;
	double cur_PBC[3][3];	
	int nframes;
	int populated;
	int read_init;
	int n_ns_atoms;
	int *ns_atoms;
	int *ns_atoms_start;
	int *ns_atoms_stop;
	void align_on_ns();
} open_frame_file;

static open_frame_file *theFrame= NULL;

// only handles psf/pdb for structure and dcd/pdb for coords
// not particularly robust or well-tested.
// LOOS up to date and tested by broader community.

int io_nframes( void )
{
	if( !theFrame )
	{
		printf("Internal error: Read hasn't started for trajectory file yet.\n");
		exit(1);
	}

	return theFrame->nframes;
}

void io_readStructureFile( char *fileName )
{
	int natoms = 0;
#ifdef HAVE_LIBLOOS

#else 
	if( strlen(fileName) < 4 || (strcasecmp(fileName+strlen(fileName)-3,"pdb") && strcasecmp(fileName+strlen(fileName)-3,"psf")) )
	{	
		printf("File must have either pdb or psf extension.\n");
		printf("For a wide variety of I/O, build with LOOS.\n");
		exit(1);
	}

	FILE *theFile = fopen(fileName,"r");

	if( !theFile ) 
	{
		printf("Couldn't open structure file '%s'.\n", fileName );
		exit(1);
	}

	if( !strcasecmp(fileName+strlen(fileName)-3,"pdb") )
		loadPSFfromPDB(theFile);
	else
		loadPSF(theFile);
	natoms = curNAtoms();
	fclose(theFile);
#endif
	if( !theFrame )
	{
		theFrame = (open_frame_file *)malloc( sizeof(open_frame_file) );
		theFrame->file_handle = NULL;
		theFrame->natoms = natoms;
		theFrame->cur_frame = (struct atom_rec *)malloc( sizeof(struct atom_rec) * theFrame->natoms );
		theFrame->last_frame = NULL;
		
		theFrame->n_ns_atoms = 0;	
		theFrame->ns_atoms = NULL;
		theFrame->ns_atoms_start = NULL;
		theFrame->ns_atoms_stop = NULL;

		theFrame->cur_PBC[0][0] = 0;
		theFrame->cur_PBC[0][1] = 0;
		theFrame->cur_PBC[0][2] = 0;
		
		theFrame->cur_PBC[1][0] = 0;
		theFrame->cur_PBC[1][1] = 0;
		theFrame->cur_PBC[1][2] = 0;
		
		theFrame->cur_PBC[2][0] = 0;
		theFrame->cur_PBC[2][1] = 0;
		theFrame->cur_PBC[2][2] = 0;

		theFrame->read_init = 0;
		theFrame->populated = 0;
		theFrame->nframes = 0;
	}	
}

void io_initialize_read( char *fileName )
{
#ifndef HAVE_LIBLOOS
	if( !theFrame->file_handle )
	{
		if( strlen(fileName) < 4 || (strcasecmp(fileName+strlen(fileName)-3,"pdb") && strcasecmp(fileName+strlen(fileName)-3,"dcd")) )
		{	
			printf("File must have either pdb or psf extension.\n");
			printf("For a wide variety of I/O, build with LOOS.\n");
			exit(1);
		}

		theFrame->file_handle = fopen(fileName,"r");

		if( !theFrame->file_handle )
		{
			printf("Couldn't open file '%s'.\n", fileName );
			exit(1);
		}

		if( !strcasecmp(fileName+strlen(fileName)-3,"pdb") )
			theFrame->nframes = 1;
		else
		{
			readDCDHeader( theFrame->file_handle );	
			setAligned();
			theFrame->nframes = curNFrames();
		}
	}
	theFrame->read_init = 1;
#endif	
}

void io_readFrame( char *fileName )
{
	if( !theFrame )
	{
		printf("Internal error: Cannot read a trajectory frame without structural information.\n");
		exit(1);
	}	
	
	if( !theFrame->read_init )
		io_initialize_read(fileName);

	int init = theFrame->populated;

	if( theFrame->populated )
	{	// if the last frame is populated, swap.
		if( !theFrame->last_frame )
			theFrame->last_frame = (struct atom_rec *)malloc( sizeof(struct atom_rec) * theFrame->natoms );
		else
		{
			for( int a = 0; a < theFrame->natoms; a++ )
				theFrame->last_frame[a].zap();
		}

		struct atom_rec *temp = theFrame->cur_frame;	
		theFrame->cur_frame = theFrame->last_frame;
		theFrame->last_frame = temp;	
	}

#ifdef HAVE_LIBLOOS

#else
	if( !strcasecmp(fileName+strlen(fileName)-3,"dcd") )
		loadFrame(theFrame->file_handle, theFrame->cur_frame);  
	else
		loadPDB(theFrame->file_handle, theFrame->cur_frame, theFrame->natoms);
	double La,Lb,Lc,alpha,beta,gamma;
	PBCD(&La,&Lb,&Lc,&alpha,&beta,&gamma);
	theFrame->cur_PBC[0][0] = La;
	theFrame->cur_PBC[1][1] = Lb;
	theFrame->cur_PBC[2][2] = Lc;
#endif

	if( !init )
	{
		int nlibrary = nLipidsInLibrary();
		int nspace = 10;
		theFrame->ns_atoms = (int *)malloc( sizeof(int) * nspace );
		theFrame->ns_atoms_start = (int *)malloc( sizeof(int) * nspace );
		theFrame->ns_atoms_stop = (int *)malloc( sizeof(int) * nspace );
		int pres = -1;
		int added = 0;
		char pseg[256];
		pseg[0] = '\0';

		struct atom_rec *at = theFrame->cur_frame;

		int cur_start = 0;
		int cur_stop  = -1;
		int cur_ns    = -1;
		
		for( int a = 0; a < theFrame->natoms; a++ )
		{
			if( (at[a].segid && strcasecmp( pseg, at[a].segid )) || at[a].res != pres || !added )
			{
				if( at[a].res != pres || (at[a].segid && strcasecmp( pseg, at[a].segid )) )
				{
					if( added && cur_ns >= 0 )
					{
						theFrame->ns_atoms_start[cur_ns] = cur_start;
						theFrame->ns_atoms_stop[cur_ns] = a-1;
					}
					cur_start = a;
					added = 0;
				}
				for( int xl = 0; xl < nlibrary; xl++ )
				{
					if( !lipidLibrary[xl].ns_atom ) continue;

					if( !strcasecmp( lipidLibrary[xl].name, at[a].resname) &&
					    !strcasecmp( lipidLibrary[xl].ns_atom, at[a].atname) )
					{
						added = 1;
						if( theFrame->n_ns_atoms == nspace )
						{
							nspace *= 2;	
							theFrame->ns_atoms = (int *)realloc( theFrame->ns_atoms, sizeof(int) * nspace );
							theFrame->ns_atoms_start = (int *)realloc( theFrame->ns_atoms_start, sizeof(int) * nspace );
							theFrame->ns_atoms_stop  = (int *)realloc( theFrame->ns_atoms_stop, sizeof(int) * nspace );
						}
						theFrame->ns_atoms[theFrame->n_ns_atoms] = a;
						theFrame->ns_atoms_start[theFrame->n_ns_atoms] = -1;
						theFrame->ns_atoms_stop[theFrame->n_ns_atoms] = -1;
						cur_ns = theFrame->n_ns_atoms;
						theFrame->n_ns_atoms+=1;
					} 
				}
			}
		
			if( strlen(at[a].segid) < 255 )
				strcpy( pseg, at[a].segid );
			else
				pseg[0] = '\0';

			pres = at[a].res;		
		}
	}

	theFrame->populated = 1;
}

int io_nAtoms( void )
{
	return theFrame->natoms;
}

int io_nNSAtoms( void )
{
	if( !theFrame )
	{
		printf("Internal error: Read hasn't started for trajectory file yet.\n");
		exit(1);
	}

	return theFrame->n_ns_atoms;	
}

void io_getNSAtomStartStop( int index, int *start, int *stop )
{
	if( index < 0 || index >= theFrame->n_ns_atoms )
	{
		printf("ERROR: requested atom index out of bounds.\n");
		exit(1);
	}

	*start = theFrame->ns_atoms_start[index];
	*stop  = theFrame->ns_atoms_stop[index];
}

int io_getNSAtom( int index )
{
	if( index < 0 || index >= theFrame->n_ns_atoms )
	{
		printf("ERROR: requested atom index out of bounds.\n");
		exit(1);
	}

	return theFrame->ns_atoms[index];
}

struct atom_rec *io_getCurrentFrame( void )
{
	return theFrame->cur_frame;
}

// align on the neutral surface atoms.

void io_align( void )
{
	theFrame->align_on_ns();
}

void open_frame_file::align_on_ns( void )
{
	struct atom_rec *at = cur_frame;
	if( n_ns_atoms == 0 ) return;
	if( last_frame == NULL ) return;

	double com_prev[3] = { 0,0,0};
	for( int t = 0; t < n_ns_atoms; t++ )
	{
		com_prev[0] += at[ns_atoms[t]].x;
		com_prev[1] += at[ns_atoms[t]].y;
		com_prev[2] += at[ns_atoms[t]].z;
	}	

	com_prev[0] /= n_ns_atoms;
	com_prev[1] /= n_ns_atoms;
	com_prev[2] /= n_ns_atoms;

	double cur_com[3] = { 0,0,0};
	
	for( int t = 0; t < n_ns_atoms; t++ )
	{
		double dr[3] = { 
			at[ns_atoms[t]].x - last_frame[ns_atoms[t]].x,
			at[ns_atoms[t]].y - last_frame[ns_atoms[t]].y,
			at[ns_atoms[t]].z - last_frame[ns_atoms[t]].z };
		double put[3];
		MinImage3D( dr, cur_PBC, put );

		cur_com[0] += last_frame[ns_atoms[t]].x + dr[0];
		cur_com[1] += last_frame[ns_atoms[t]].y + dr[1];
		cur_com[2] += last_frame[ns_atoms[t]].z + dr[2];
	}	

	cur_com[0] /= n_ns_atoms;
	cur_com[1] /= n_ns_atoms;
	cur_com[2] /= n_ns_atoms;

	for( int t = 0; t < n_ns_atoms; t++ )
	{
		cur_frame[t].x -= (cur_com[0] - com_prev[0]);
		cur_frame[t].y -= (cur_com[1] - com_prev[1]);
		cur_frame[t].z -= (cur_com[2] - com_prev[2]);
	}
}

int addStructureToPool( const char *fileNameStruct, const char *fileNamePSF )
{
	for (struct pool_structure *t = thePool; t; t = t->next )
	{
		if( !strcmp( fileNameStruct, t->fileNameStruct) && (fileNamePSF && !strcmp(fileNamePSF, t->fileNamePSF )) )
			return t->id; 
	}

	struct pool_structure *pool = (struct pool_structure*)malloc( sizeof( struct pool_structure ) );

	FILE *structFile = fopen(fileNameStruct,"r");	
	
	if( !structFile )
	{
		printf("Couldn't open file \"%s\"\n", fileNameStruct );
		exit(1);
	}

	int loaded_psf = 0;
	
	int **basis=NULL;
	int *basis_length=NULL;
	int nbasis=0;

	if( fileNamePSF )
	{
		FILE *psfFile = fopen(fileNamePSF,"r");
	
		if( !psfFile )
		{
			printf("Couldn't open file \"%s\"\n", psfFile );
			exit(1);
		}
		
		loadPSF( psfFile );
		loaded_psf = 1;

		fclose(psfFile);
	}
	else
	{
		loadPSFfromPDB( structFile ); 		
		if( curNAtoms() == 0 )
		{
			printf("ERROR. Failed to read any atoms from file '%s'.\n", fileNameStruct );
			exit(1); 
		}
	}
	int nat = curNAtoms();


	struct atom_rec *at = (struct atom_rec *)malloc( sizeof(struct atom_rec ) * curNAtoms() );

	rewind(structFile);
	int nread = loadPDB(structFile, at ); 

	if( nread != curNAtoms() )
	{
		printf("ERROR reading atoms from pdb '%s'. Expected %d but found %d.\n",
			fileNameStruct, curNAtoms(), nread );
		exit(1);			
	}

	fclose(structFile);			
	
	double Lx,Ly,Lz;
	double a,b,g;

	PBCD( &Lx, &Ly, &Lz, &a, &b, &g );

	pool->Lx = Lx;
	pool->Ly = Ly;
	pool->Lz = Lz;

	int *bond_offsets = NULL;
	int *nbonds = NULL;
	int *bond_list = NULL; 	


	if( loaded_psf ) // get rings of the molecules for penetration search.
	{
		fetchCycleBasis( &basis, &basis_length, &nbasis );

		int nb = getNBonds();
		int *bonds = (int *)malloc( sizeof(int) * nb * 2 );
		getBonds(bonds);

		bond_offsets = (int *)malloc( sizeof(int) * nat );
		nbonds = (int *)malloc( sizeof(int) * nat);
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

		bond_list = (int *)malloc( sizeof(int) * off );
		memset( nbonds, 0, sizeof(int) * nat );

		for( int b = 0; b < nb; b++ )
		{
			int a1 = bonds[2*b+0];
			int a2 = bonds[2*b+1];
			bond_list[bond_offsets[a1]+nbonds[a1]] = a2;
			nbonds[a1] += 1;	
			bond_list[bond_offsets[a2]+nbonds[a2]] = a1;
			nbonds[a2] += 1;	
		}						
	}

	double *lipid_xyz = pool->xyz = (double *)malloc(sizeof(double)*3*curNAtoms());
	int *leaflet = pool->leaflet = (int *)malloc( sizeof(int) * curNAtoms() );
	int *lipid_start = pool->lipid_start = (int *)malloc( sizeof(int) * curNAtoms() );
		// the index where a lipid's atoms stop
	int *lipid_stop = pool->lipid_stop  = (int *)malloc( sizeof(int) * curNAtoms() );

	pool->at = at;
	pool->nat = nat;

	pool->cycles = basis;
	pool->cycle_lengths = basis_length;
	pool->ncycles = nbasis;
	
	pool->bonds = bond_list;
	pool->nbonds = nbonds;
	pool->bond_offsets = bond_offsets;
	
	
	char p_segid[256];
	int pres = -1;
	int seg_continuity_mode = 0;
	p_segid[0] = '\0';
	// fill the index

	int protein_mode = !strncasecmp( at[0].segid, "PRO", 3);	
	int nlipids=-1;
	for( int a = 0; a < curNAtoms(); a++ )
	{
		if( !strcasecmp( at[a].resname, "TIP3") ) continue;
		if( !strcasecmp( at[a].resname, "W") ) continue;
		if( !strcasecmp( at[a].resname, "SOD") ) continue;
		if( !strcasecmp( at[a].resname, "POT") ) continue;
		if( !strcasecmp( at[a].resname, "CLA") ) continue;
		if( !strcasecmp( at[a].resname, "CAL") ) continue;

		if( !strncasecmp( at[a].segid, "GLP", 3) ||
	  	    !strncasecmp( at[a].segid, "PRO", 3) )
			seg_continuity_mode = 1;
		else
			seg_continuity_mode = 0;

		if( protein_mode )
		{
			// protein_mode untested. remove this comment when you've tested it.
		}
		else if( (!seg_continuity_mode && (at[a].res != pres)) || strcasecmp( at[a].segid, p_segid) )
		{  
			nlipids++;
			// a new residue. is it a lipid, or protein?
			// for now, assume everything is going in, except water.
			
			lipid_start[nlipids] = a;
			lipid_stop[nlipids] = a;
				
		}
		else
		{
			if( nlipids >= 0 )
				lipid_stop[nlipids] = a;
		}
		strcpy( p_segid, at[a].segid ); 
		pres = at[a].res;
		
		if( !strncasecmp( at[a].segid, "PRO", 3) )
			protein_mode = 1;
		else
			protein_mode = 0;
	}

	// the last "lipid" doesn't get incremented.
	nlipids++;
		
	// get the bilayer center.
	
#define N_BINS_MOLDIST 100
	
        double best_chi2 = 1e10;
	double wrapto = 0;
 
        double moldist[N_BINS_MOLDIST];
        memset( moldist, 0, sizeof(double) * N_BINS_MOLDIST );


	for( int l = 0; l < nlipids; l++ )
	{
		for( int p = lipid_start[l]; p <= lipid_stop[l]; p++ )
		{
                      double tz = at[p].z;

                      while( tz < 0 ) tz += Lz;
                      while( tz >= Lz ) tz -= Lz;

                      int zb = N_BINS_MOLDIST * tz / Lz; // this is right
                      if( zb < 0 ) zb = 0;
                      if( zb >= N_BINS_MOLDIST ) zb = N_BINS_MOLDIST-1;
                      moldist[zb] += 1;
		}
	}

         for( int zb = 0; zb < N_BINS_MOLDIST; zb++ )
         {
                 double zv = Lz * (zb+0.5) / (double)N_BINS_MOLDIST;
 
                  int zlow  = zb- N_BINS_MOLDIST/2;
                  int zhigh = zlow + N_BINS_MOLDIST;
 
                  double Lzhi2 = 0;
                  for( int iz = zlow; iz < zhigh; iz++ )
                  {
                          double dz = Lz * (iz+0.5) / N_BINS_MOLDIST - zv;
 
                          int iiz = iz;
                          while( iiz < 0 ) iiz += N_BINS_MOLDIST;
                          while( iiz >= N_BINS_MOLDIST ) iiz -= N_BINS_MOLDIST;
 
                          Lzhi2 += moldist[iiz] * (dz) * (dz);
                  }
 
                  if( Lzhi2 < best_chi2 )
                  {
                          best_chi2 = Lzhi2;
                          wrapto = zv;
                  }
         }



	// wrap around z periodic dimension

	for( int l = 0; l < nlipids; l++ )
	{
		int midlipid = (lipid_start[l]+lipid_stop[l])/2;
		while( at[midlipid].z - wrapto < -Lz/2 ) at[midlipid].z += Lz;
		while( at[midlipid].z - wrapto > Lz/2 ) at[midlipid].z -= Lz;

		double lcom_z = 0;

		for( int p = lipid_start[l]; p <= lipid_stop[l]; p++ )
		{
			while( at[p].x - at[midlipid].x < -Lx/2 ) at[p].x += Lx;
			while( at[p].x - at[midlipid].x >  Lx/2 ) at[p].x -= Lx;
			
			while( at[p].y - at[midlipid].y < -Ly/2 ) at[p].y += Ly;
			while( at[p].y - at[midlipid].y >  Ly/2 ) at[p].y -= Ly;
			
			while( at[p].z - at[midlipid].z < -Lz/2 ) at[p].z += Lz;
			while( at[p].z - at[midlipid].z >  Lz/2 ) at[p].z -= Lz;
	
			lcom_z += (at[p].z - wrapto);
		}

		if( lcom_z > 0 )
			leaflet[l] = 1;
		else
			leaflet[l] = -1;
	} 

	// subtract off wrapto

	for( int a = 0; a < curNAtoms(); a++ )
		at[a].z -= wrapto;

	// wrap proteins around z.
	
	char pseg[256] = { '\0' };

	if( at[0].segid )
		strcpy( pseg, at[0].segid );

	int start_flag = 1;
	int protein_seg_start = -1;
	int protein_seg_stop  = -1;
	for( int a = 0; a <= curNAtoms(); a++ )
	{
		if( protein_seg_start >= 0 && ( a== curNAtoms() || strcasecmp( at[a].segid, pseg)) )
		{
			double com_z = 0;

			for( int a2 = protein_seg_start; a2 <= protein_seg_stop; a2++ )
				com_z += at[a2].z;
			com_z /= (protein_seg_stop-protein_seg_start+1);
			double shift_z = 0;
			while(com_z + shift_z < -Lz/2 ) shift_z += Lz;
			while(com_z + shift_z >  Lz/2 ) shift_z -= Lz;
			// this is what we're actually trying to accomplish:
			for( int a2 = protein_seg_start; a2 <= protein_seg_stop; a2++ )
				at[a2].z += shift_z;		

			protein_seg_start=-1;
			protein_seg_stop =-1;
		}

		if( a == curNAtoms() )
			break;

		int is_protein = !strncasecmp( at[a].segid, "PRO", 3);

		// terminating a structure?

		if( is_protein && (start_flag || strcasecmp(at[a].segid, pseg ) ) )
			protein_seg_start = a;
		else if( is_protein ) // continuation
			protein_seg_stop = a;	

		start_flag = 0;

		strcpy( pseg, at[a].segid );
	}


	// fill lipid positions

	for( int l = 0; l < nlipids; l++ )
	{
		lipid_xyz[3*l+0] = 0;	
		lipid_xyz[3*l+1] = 0;	
		lipid_xyz[3*l+2] = 0;
	
		for( int p = lipid_start[l]; p <= lipid_stop[l]; p++ )
		{
			lipid_xyz[3*l+0] += at[p].x;
			lipid_xyz[3*l+1] += at[p].y;
			lipid_xyz[3*l+2] += at[p].z;
		}

		lipid_xyz[3*l+0] /= (lipid_stop[l] - lipid_start[l] + 1 );
		lipid_xyz[3*l+1] /= (lipid_stop[l] - lipid_start[l] + 1 );
		lipid_xyz[3*l+2] /= (lipid_stop[l] - lipid_start[l] + 1 );
	}

	pool->nlipids = nlipids;
	
	int id = 0;

	for( pool_structure *tp = thePool; tp; tp = tp->next )
	{
		if( id <= tp->id ) id = tp->id+1;
	}

	pool->fileNameStruct = (char *)malloc( sizeof(char) *(1+strlen(fileNameStruct) ) );
	if( fileNamePSF )
		pool->fileNamePSF = (char *)malloc( sizeof(char) *(1+strlen(fileNamePSF) ) );
	else
		pool->fileNamePSF = NULL;
	strcpy(pool->fileNameStruct, fileNameStruct);
	if( fileNamePSF )
		strcpy(pool->fileNamePSF, fileNamePSF);

	pool->next = thePool;
	pool->id = id;
	thePool = pool;

	return pool->id;
}

struct pool_structure *getPool(int id)
{
	for( struct pool_structure *pool = thePool; pool; pool = pool->next )	
	{
		if( pool->id == id )
			return pool;
	}

	printf("Logical error. Couldn't find pool with id '%d'.\n", id );
	exit(1);

	return NULL;
}

void deleteFromPool( int id )
{

}


