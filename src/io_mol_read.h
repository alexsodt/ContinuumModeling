#ifndef __io_mol_readh__
#define __io_mol_readh__
#include "pdb.h"


int io_nframes( void);
void io_readStructureFile( char *fileName );
void io_readFrame( char *fileName );
struct atom_rec *io_getCurrentFrame();
int io_nNSAtoms( void );
int io_nAtoms( void );
int io_getNSAtom( int index );
void io_align( void );
void io_initialize_read( char *fileName );
void io_getNSAtomStartStop( int index, int *start, int *stop );
void io_get_PBC( double *Lx, double *Ly, double *Lz );

// the structure for how lipids are stored and extracted to build an all-atom sim.


struct pool_structure
{
	int id; // for loading this pool'd structure.

	char *fileNameStruct;
	char *fileNamePSF; // to see if it is already loaded.

	int nlipids;
	int *lipid_start;
	int *lipid_stop;
	int *leaflet;
	double *xyz;

	int **cycles;
	int *cycle_lengths;
	int ncycles;

	int *bonds;
	int *nbonds;
	int *bond_offsets;

	struct atom_rec *at;
	int nat;

	int has_pbc;

	double Lx,Ly,Lz;

	struct pool_structure *next;
};

#ifdef __io_mol_readc__
struct pool_structure *thePool=NULL;
#else
extern struct pool_structure *thePool;
#endif

int addStructureToPool( const char *fileNameStruct, const char *fileNamePSF );
struct pool_structure *getPool(int id);
void deleteFromPool(int id);
int poolIsMartini( int id );

#endif




