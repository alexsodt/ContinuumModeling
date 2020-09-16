#include "pdb.h"
#include "dcd.h"
#include "io_mol_read.h"

// need to set this with configure.
const char *libDir = "/Users/sodtaj/git_projects/HD/lib";

struct map31
{
        const char *threeCode;
        char oneCode;
        char cgCode;
};

static const map31 themap31[] =
{
        {"Ala",'A','V'},
        {"Arg",'R','L'},
        {"Asn",'N','L'},
        {"Asp",'D','L'},
        {"Cys",'C','B'},
        {"CysP",'1','B'},
        {"Glu",'E','L'},
        {"Gln",'Q','L'},
        {"Gly",'G','N'},
        {"His",'H','N'},
        {"Ile",'I','B'},
        {"Leu",'L','B'},
        {"Lys",'K','N'},
        {"Met",'M','B'},
        {"Phe",'F','B'},
        {"Pro",'P','N'},
        {"Ser",'S','N'},
        {"Thr",'T','N'},
        {"Trp",'W','B'},
        {"Tyr",'Y','B'},
        {"Val",'V','V'},
        {"UNK",'U','U'},
        {"HSE",'H','N'},
        {"HSD",'H','N'}
};	

char threeToOne( const char *code )
{
	int nl = sizeof(themap31)/sizeof(map31);

        for( int c = 0; c < nl; c++ )
        {   
                if( !strcasecmp( themap31[c].threeCode, code ) ) 
                        return themap31[c].oneCode;
        }    

        printf("code '%s' does not map.\n", code);
        exit(1);
}


int pdbFetch( struct atom_rec **out_pdb, int *nout, const char *dir, const char *file, int addToPool )
{
	int lenLIB = strlen(libDir);
	int lenDir = strlen(dir);
	int lenFile = strlen(file);
	int lenTot = lenLIB + lenDir + lenFile + 256;
	char fileName[lenTot];
	sprintf(fileName, "%s/%s/%s.pdb", libDir, dir, file );
	
	char fileName_PSF[lenTot];
	sprintf(fileName_PSF, "%s/%s/%s.psf", libDir, dir, file );
	
	int pool=-1;

	if( !addToPool ) 
	{
		FILE *theFile = fopen(fileName,"r");
	
		if( !theFile )
		{
			printf("Couldn't open file '%s' from pdbFetch.\n", fileName );
			exit(1);
		}
		int usePSF = 0;
		FILE *theFilePSF = fopen(fileName_PSF,"r");
	
		if( theFilePSF )
		{
			usePSF = 1;
			loadPSF( theFilePSF );
		}
		else
			loadPSFfromPDB( theFile );
	
		int nat = curNAtoms();
	
		(*out_pdb) = (struct atom_rec *)malloc( sizeof(struct atom_rec) * nat );
		*nout = nat;
	
		rewind(theFile);
	
		loadPDB( theFile, *out_pdb );	
	
		fclose(theFile);
	
		if( theFilePSF )
			fclose(theFilePSF);
	}
	else
	{	
		const char *usePSFFile = fileName_PSF;
		FILE *checkFile = fopen(fileName_PSF,"r");
		if( !checkFile )
			usePSFFile = NULL;
		else
			fclose(checkFile);	

		pool = addStructureToPool(fileName,usePSFFile); 
		pool_structure *thePool = getPool(pool);

		*nout = thePool->nat;
		*out_pdb = thePool->at;
	}
	return pool;
}

