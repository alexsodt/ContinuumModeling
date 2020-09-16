#ifndef __pdbfetchh__
#define __pdbfetchh__

#include "pdb.h"
#include "dcd.h"

const int addToPool = 1;
const int doNotAddToPool = 0;

int pdbFetch( struct atom_rec **out_pdb, int *nout, const char *dir, const char *file, int addToPool=0 );
char threeToOne( const char *code );

#define ION_CAL		0	
#define ION_MG		1
#define ION_FE		2

struct ion_add
{
	int type;
	double x;
	double y;
	double z;
};

#endif
