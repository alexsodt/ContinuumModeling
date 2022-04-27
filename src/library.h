#ifndef __libraryh__
#define __libraryh__

struct fixed_lipid
{	
	const char *name;
	double c0;
	double APL;
	const char *ns_atom;
	const char *ns_atom_M;
};

struct lipid
{
	char *name;
	double c0;
	double APL;
	char *ns_atom;
	char *ns_atom_M;
};

#ifdef __lipidcompositionc__
lipid *lipidLibrary = NULL;
#else
extern lipid *lipidLibrary;
#endif
int nLipidsInLibrary(void);
#endif
