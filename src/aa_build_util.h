#ifndef __aa_build_util_h__
#define __aa_build_util_h__

struct caa_box
{
	int np;
	int npSpace;
	int *plist;
};


// for clash workers.
//
#define ACTION_REPORT		0
#define ACTION_DELETE		1

struct aa_build_data
{
	int system_pool;

	int **global_cycles;
	int *global_cycle_len;
	int global_ncycles;
	int global_nspace;
	
	int *global_bonds;
	int *global_nbonds;
	int *global_bond_offsets;
	int global_bonds_tot;
	int global_bond_space;
	int global_atom_space;
	
	int nplaced;
	int nplaced_pcut;
	int nplacedSpace;
	double *placed_atoms;
	int *deleted;

	char *charmm_delete_buffer;
	int del_space;

	caa_box *theBoxes;
	caa_box *cycleBoxes;

	int add_mode;

	double PBC_vec[3][3];
	int nx, ny, nz;
	int nx_c, ny_c, nz_c;

	int halve[3];
	double frac_cen[3]; // the center point of the fractional grab.
	void init( void );
	void setupBoxing( double Lx, double Ly, double Lz, int nx, int ny, int nz);

	void addFromPool( int pool_code, double *pcopy, int *map );
	void addLipidFromPool( int pool_code, double *pcopy, int l );

	int curPlace(void);
	int addAtom( double *r );
	int checkMappedBonds( double *coords, int *map, int nmapped, int *bonds, int *bond_offsets, int *nbonds );
	int checkMappedCycles( double *coords, int *map, int nmapped, int **cycles, int *cycle_lens, int ncycles );
	void addMappedBonds( int offset, int *map, int nmapped, int *bonds, int *bond_offsets, int *nbonds );
	void addBondsInRun( int offset, int a_start, int a_stop, int *bonds, int *bond_offsets, int *nbonds );
	void addMappedCycles( int offset, double *coords, int *map, int nmap, int **cycles, int *cycle_lens, int ncycles );
	void addCyclesInRun( int offset, double *coords, int a_start, int a_stop, int **cycles, int *cycle_lens, int ncycles );
	void addSpecialCycle( int *cycle, int len, double *coords );
	void setupBoxing( double PBC_in[3][3], int nx_in, int ny_in, int nz_in);
	int cycleClash( double *coords, int a_start, int *cycle, int len );
	int bondClash( double *r1, double *r2 );
	int nclash_aa( double *coords, int lc, int is_mod, double cutoff=0.5, int *halve=NULL );
	
	int nclash_aa_worker( double *coords, int lc, int is_mod, double cutoff=0.5, int *halve=NULL, int action=ACTION_REPORT );
	int bondClash_worker( double *r1, double *r2, int action=ACTION_REPORT, int i=-1, int j=-1 );
	int cycleClash_worker( double *coords, int a_start, int *cycle, int len, int action=ACTION_REPORT );
	
	int getPcut( void ) { return nplaced_pcut; }
	void markPcut(void) { nplaced_pcut = nplaced; }
	void addCrossedBonds( int a_start, int a_stop ); // someone please notice I restrained myself from calling this fn hot crossed bonds
	
	void addSystemPool( int pool_code );
	int deleteResidue( int index );
	void deleteClashes( int start, int num );
	int inAddMode( void );
	int activateAddMode( void );
	int deactivateAddMode( void );
};


void boxit( double *r_in, int index, caa_box *theBoxes, double PBC_vec[3][3], int nx, int ny, int nz ); 
void autoBonds( struct atom_rec *at, int nat, int **bonds_out, int *nbonds_out);

#endif

