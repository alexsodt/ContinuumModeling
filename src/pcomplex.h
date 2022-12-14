#ifndef __pcomplexh__

#include "simulation.h"
#include "interp.h"
#include "input.h"
#include "pdbFetch.h"

struct atom_rec;
struct surface_mask;
struct aa_build_data;

#define __pcomplexh__

#define DEBUG_OFF	0
#define DEBUG_NO_V	1
#define DEBUG_NO_T	2

#define ALIGN_TYPE_SURF_SURF	0
#define ALIGN_TYPE_CURVATURE	1
#define ALIGN_TYPE_SURF_AQ	2

#define NO_BUILD		-1
#define BUILD_SEQUENCE		0	// build a pcomplex from its protein sequence
#define BUILD_GENERAL		1	// build a protein from its included RTF.

struct pcomplex
{
	int disabled; //check
	int nwatchers; //check

	int nmer_saved;
	int is_inside; // copyparent
	int nsites; //check
	int nattach; //check
	int bound; //check
	int debug; //check
	int do_bd; //check
	int my_id;

	int *rd_timestep_disabled; // set to true if we disable its diffusion this timestep.

	char *complex_name; //check
	
	// diffusion constants.
	double *DC; //check

	// all masses.
	double *mass; //check

	// all WCA/LJ radii
	double *sigma; //check
	double *att_eps; //check
	double *att_sigma; //check

	// all coordinates in 3D 
	double *rall; //check
	// attach coordinates in f/uv; can exceed face bounds
	int *sid; //check // id of surface we are on.
	int *stype; //check // if doing rxn/diff
	int *fs;	
	double *puv;

	// in the final face bounds.
	int *grad_fs;	
	double *grad_puv;
	
	double *p; // Hamiltonian conjugate momenta
	double *qdot; // time derivatives of generalized coordinates.
	double *save_grad;
	double *cache_grad; // if a particle move is rejected we may need to save the gradient
	double *cache_rall; // if a particle move is rejected we may need to save the gradient
	double *cache_puv; // if a particle move is rejected we may need to save the gradient
	int *cache_f; // if a particle move is rejected we may need to save the gradient
	
	double *PBC_ext; // the current PBC vector of each site.
	double *last_pos;

	double *p_area;
	double *p_c0;

	double *coord_transform;

// This is added by Kayla for keeping track of when to "delete" a complex
	int alive; //1: alive, 0:dead -- remove from memory during garbage cleanup  
// End of what Kayla added

	virtual pcomplex *clone( void ) { printf("Clone function not implemented for this pcomplex class.\n"); exit(1); }
	virtual void clone( Simulation *theSimulation, surface *theSurface, double *rsurf, pcomplex *takeFrom, int add_to);

	virtual void alloc( void );
	virtual void pfree( void );

	virtual int getNBonds( void ) { return 0; };
	virtual void putBonds( int *bond_list, double *bond_r = NULL, double *bond_k = NULL );

	virtual void init( Simulation *theSimulation, surface *theSurface, double *, int f, double u, double v, int nmer=1 ); 
	virtual void init( double *r ); 

	virtual void bind( int f, double u, double v );
	virtual void unbind( void );

	virtual void applyParamsToComplex( double *p );
	virtual void getParamsFromComplex( double *p );
	virtual int nparams( void );
	virtual void orient( surface *theSurface, double *rsurf ) { }

	virtual void refresh( Simulation *theSimulation  );
	virtual double V( Simulation *theSimulation );
	virtual double grad( Simulation *theSimulation, double *surface_g, double *particle_g ); 
	virtual void fd_grad_debug(surface *theSurface, double *rsurf ); 
	virtual void loadParams( parameterBlock *block );

	virtual void move_inside( void );
	virtual void move_outside( void );

	virtual int packLenF( void );
	virtual int packLenI( void );
	virtual void packageF( double * );
	virtual void unpackF( double * );
	virtual void packageI( int *);
	virtual void unpackI( int *); 
	
	virtual int isElastic(void) { return 0; }	

	virtual char getSiteCode( int ind );

	void activateBrownianDynamics(void);

	void printType( char **type );
	void cacheVelocities( void );
	void prepareForGradient( void );
	void loadCoords( surface *theSurface, double *rsurf, double *r, double *n );
	void setrall( Simulation *theSimulation  );
	void base_init( void );
	double T(Simulation *theSimulation, int subp=-1 );
	double update_dH_dq( Simulation *theSimulation, double time_step=-1, double timestep_total=0 );
	void propagate_p( Simulation *theSimulation, double dt );
	//void compute_qdot( surface *theSurface, double *rsurf, double *mesh_qdot0, double *mesh_qdot, double frac_mult=1.0 );
	void compute_qdot( Simulation *theSimulation,  double frac_mult=1.0 );
	void propagate_surface_q( Simulation *theSimulation, double dt );
	void copyParentParameters( pcomplex *parent );

	void debug_dPinv_dq( surface * theSurface, double *rsurf  );
	void getMeshQxdot( surface *theSurface, double *rsurf, double *Minv, double *mesh_p, double *mesh_qdot, double *mesh_qdot0, double *mesh_der_qdot );
	void applyLangevinFriction( Simulation *theSimulation, double dt, double gamma );
	void applyLangevinNoise( Simulation *theSimulation, double dt, double gamma, double temperature );
	int saveComplex( char *buffer, int *buflen, int bufLenMax );
	void saveComplex( FILE * theFile );
	int loadComplex( FILE * theFile, Simulation *theSimulation, int load_to );
	double AttachV( Simulation *theSimulation );
	double AttachG( Simulation *theSimulation, double *pg );
	void evaluate_momentum( surface *theSurface, double *rsurf, double *pout );
	double local_curvature( Simulation *theSimulation);
	void print_type( char **outp );
	virtual void writeStructure( Simulation *theSimulation, surface_mask *upperSurfaceMask, surface_mask *lowerSurfaceMask, struct atom_rec **at, int *nat, char ***seq, int *nseq, int **seq_at_array, char ***patches, ion_add **ions, int *nions, aa_build_data *buildData, int *build_type );
	virtual void get( 
		Simulation *theSimulation, // this is the surface/PBC
		struct atom_rec *at,
		int p_start, int p_stop, int nat_tot);


	void disable( void) { disabled = 1; }
	void destroy( void);
	void watch( void );
	void forget( void );

	void cache( void);
	void uncache(void);
	

	// for building all-atom.
	int addPeripheralProteinHelper(Simulation *theSimulation, surface_mask *upperMask, surface_mask *lowerMask,
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
		);
	
	int addCurvatureOrientedPeripheralProteinHelper(Simulation *theSimulation, surface_mask *upperMask, surface_mask *lowerMask,
		struct atom_rec **at_out,
		int *nat_out,
		struct ion_add **ions,
		int *nions,
		int pool_code, // the pool code 
		const char *segid, // the segid of the protein to extract
		int orient_res[3], // the residues in the protein used to define the orientation.
		int attach_site, // the site in the complex used to define the attachment.
		int align_on_neg, // align on negative curvature?
		aa_build_data *buildData, // data structure to put build information (PSF-derived info and atom placements for detecting clashes). 
		double noise
		);

	int addGeneralProteinHelper(Simulation *theSimulation, surface_mask *upperMask, surface_mask *lowerMask,
		struct atom_rec **at_out,
		int *nat_out,
		struct ion_add **ions,
		int *nions,
		int pool_code, // the pool code 
		const char *segid, // the segid of the protein to extract
		int orient_res[3], // the residues in the protein used to define the orientation.
		int attach_sites[2], // the sites in the complex used to define the attachment.
		int aqueous_sites[2], // the sites in the complex used to define the attachment.
		int align_on_neg, // align on neg curvature
		int align_type,
		aa_build_data *buildData, // data structure to put build information (PSF-derived info and atom placements for detecting clashes). 
		double noise
		);

	int addAlignedProteinHelper( Simulation *theSimulation, surface_mask *upperMask, surface_mask *lowerMask,
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
		double noise,
		int multi_segment=0,
		double nrm_displace=0
		);
};

struct actin : pcomplex
{
	int *interfaces; // array of avaiable interfaces by sub ID (for no branching - "0" and "nattach-1")
	double *r_int; // position of the available interfaces
	double k_on;
	double k_off;
};

struct simpleParticle : pcomplex
{
};

struct CNT : pcomplex
{
	void writeStructure( Simulation *theSimulation, surface_mask *upperSurfaceMask, surface_mask *lowerSurfaceMask, struct atom_rec **at, int *nat, char ***seq, int *nseq, int **seq_at_array,  char ***patches, ion_add **ions, int *nions, struct aa_build_data *buildData, int *add_type );
};

struct simpleLipid : pcomplex
{
	double c0_val;
	void loadParams( parameterBlock *block );
	virtual void init(  Simulation *theSimulation,surface *theSurface, double *, int f, double u, double v, int nmer=1 ); 
};

struct simpleBound : simpleLipid
{
	double bound_sigma;
	void loadParams( parameterBlock *block );
	virtual void init(  Simulation *theSimulation,surface *theSurface, double *, int f, double u, double v, int nmer=1 ); 
};

struct simpleDimer : simpleLipid
{
	void loadParams( parameterBlock *block );
	virtual void init(  Simulation *theSimulation,surface *theSurface, double *, int f, double u, double v, int nmer=1); 
};

struct NBAR : pcomplex
{
	double bond_length;
	double k_phi;
	double k_theta;
	double bond_k;
	double theta_0;
	double phi_0;

	virtual void getDirName( char put[256]);

	void init( Simulation *theSimulation, surface *, double *rsurf, int f, double u, double v, int nmer=1 ); 
	void init( double *r );
	void bind( int f, double u, double v);
	void unbind( void );
	void loadParams( parameterBlock *block );
	void orient( surface *theSurface, double *rsurf );
	
	int getNBonds( void );
	virtual void putBonds( int *bond_list, double *bond_r = NULL, double *bond_k = NULL );

	double V( Simulation *theSimulation );
	double grad( Simulation *theSimulation, double *surface_g, double *particle_g );
	
	void writeStructure( Simulation *theSimulation, surface_mask *upperSurfaceMask, surface_mask *lowerSurfaceMask, struct atom_rec **at, int *nat, char ***seq, int *nseq, int **seq_at_array,  char ***patches, ion_add **ions, int *nions, struct aa_build_data *buildData, int *build_type );

	void move_inside(void);
	void move_outside(void);	
};

struct DVL : pcomplex // disshevelled
{
	void get( 
		Simulation *theSimulation, // this is the surface/PBC
		struct atom_rec *at,
		int p_start, int p_stop, int nat_tot );
	void getLocalLipidComposition( 
		Simulation *theSimulation, 
		struct atom_rec *at,
		int dvl_start, int dvl_stop, int nat_tot );
};

struct syt7 : pcomplex
{

	void init( Simulation *theSimulation, surface *, double *rsurf, int f, double u, double v, int nmer=1 ); 
	void init( double *r );
	void bind( int f, double u, double v);
	void unbind( void );
	void loadParams( parameterBlock *block );
	
	int getNBonds( void );
	virtual void putBonds( int *bond_list, double *bond_r = NULL, double *bond_k = NULL );

	double V( Simulation *theSimulation );
	double grad( Simulation *theSimulation, double *surface_g, double *particle_g );
	char getSiteCode( int ind);	
	void get( 
		Simulation *theSimulation, // this is the surface/PBC
		struct atom_rec *at,
		int p_start, int p_stop, int nat_tot );
	void getLocalLipidComposition( 
		Simulation *theSimulation, 
		struct atom_rec *at,
		int syt7_start, int syt7_stop, int nat_tot );

	//void writeStructure( Simulation *theSimulation, struct atom_rec **at, int *nat );
	void writeStructure( Simulation *theSimulation, surface_mask *upperSurfaceMask, surface_mask *lowerSurfaceMask, struct atom_rec **at, int *nat, char ***seq, int *nseq, int **seq_at_array /*where the atoms of a segment start and stop*/,  char ***patches, ion_add **ions, int *nions, struct aa_build_data *buildData, int *build_type );
	void move_inside(void);
	void move_outside(void);	
};

struct C2Domain : pcomplex
{
	virtual int resStart( ) { return 1; } 
	virtual int resStop( ) { return 255; } 
	virtual int putCBLResidues( int put[5] ) { put[0] = 166; put[1] = 172; put[2] = 225; put[3] = 227; put[4] = 233; }

	void get( 
		Simulation *theSimulation, // this is the surface/PBC
		struct atom_rec *at,
		int p_start, int p_stop, int nat_tot );
	void getLocalLipidComposition( 
		Simulation *theSimulation, 
		struct atom_rec *at,
		int c2_start, int c2_stop, int nat_tot );
};

struct C2A : C2Domain
{
	virtual int resStart( ) { return 1; } 
	virtual int resStop( ) { return 255; } 
	virtual int putCBLResidues( int put[5] ) { put[0] = 166; put[1] = 172; put[2] = 225; put[3] = 227; put[4] = 233; }
	virtual int altContactSite( int put[5] ) { put[0] = 297; put[1] = 303; put[2] = 357; put[3] = 359; put[4] = 365; }
};

struct C2B : C2Domain
{
	virtual int resStart( ) { return 256; } 
	virtual int resStop( ) { return 500; } 
	virtual int putCBLResidues( int put[5] ) { put[0] = 297; put[1] = 303; put[2] = 357; put[3] = 359; put[4] = 365; }
	// for alignment:
	virtual int altContactSite( int put[5] ) { put[0] = 297; put[1] = 303; put[2] = 357; put[3] = 359; put[4] = 365; }
};


struct dynamin : pcomplex
{

	virtual pcomplex *clone( void ) { return new dynamin; }
	void init( Simulation *theSimulation, surface *, double *rsurf, int f, double u, double v, int nmer ); 
	void init( double *r, int nmer );
	void bind( int f, double u, double v);
	void unbind( void );
	void loadParams( parameterBlock *block );

	virtual void getDirName( char put[256]);
	virtual void getParams(void);
	virtual int getNDihe(void);
	virtual void putDihe( int *dihe_list, double *dihe_theta, double *dihe_k );
	virtual int getNAngles( void );
	virtual void putAngles( int *angle_list, double *angle_th, double *angle_k );
	
	virtual int getNBonds( void );
	virtual void putBonds( int *bond_list, double *bond_r = NULL, double *bond_k = NULL );

	double V( Simulation *theSimulation );
	double grad( Simulation *theSimulation, double *surface_g, double *particle_g );
	char getSiteCode( int ind);	

	void get( 
		Simulation *theSimulation, // this is the surface/PBC
		struct atom_rec *at,
		int p_start, int p_stop, int nat_tot );

	void placeSubunit( Simulation *theSimulation, surface *theSurface, double *rsurf, int m, double dynamin_axis[3] );
	void clone( Simulation *theSimulation, surface *theSurface, double *rsurf, pcomplex *takeFrom, int add_to);
	//void writeStructure( Simulation *theSimulation, struct atom_rec **at, int *nat );
	void writeStructure( Simulation *theSimulation, surface_mask *upperSurfaceMask, surface_mask *lowerSurfaceMask, struct atom_rec **at, int *nat, char ***seq, int *nseq, int **seq_at_array,  char ***patches, ion_add **ions, int *nions, struct aa_build_data *buildData, int *build_type );
	void move_inside(void);
	void move_outside(void);	
};

struct sdynamin : dynamin
{
	virtual pcomplex *clone( void ) { return new sdynamin; }
	virtual void getDirName( char put[256]);
	void getParams(void);
	int getNDihe(void);
	void putDihe( int *dihe_list, double *dihe_theta, double *dihe_k );
	int getNAngles( void );
	void putAngles( int *angle_list, double *angle_th, double *angle_k );
	
	int getNBonds( void );
	void putBonds( int *bond_list, double *bond_r = NULL, double *bond_k = NULL );

};

#define IFITM3_BASE_AH	0
#define	IFITM3_P2	1
#define IFITM3_P2_W60A  2
#define IFITM3_P2_F63Q  3
#define IFITM3_P2_F67Q  4

struct ifitm3 : pcomplex
{

	void init( Simulation *theSimulation, surface *, double *rsurf, int f, double u, double v, int nmer=1 ); 
	void init( double *r );
	void bind( int f, double u, double v);
	void unbind( void );
	void loadParams( parameterBlock *block );
	
	int getNBonds( void );
	virtual void putBonds( int *bond_list, double *bond_r = NULL, double *bond_k = NULL );

	double V( Simulation *theSimulation );
	double grad( Simulation *theSimulation, double *surface_g, double *particle_g );

	//void writeStructure( Simulation *theSimulation, struct atom_rec **at, int *nat );
	void writeStructure( Simulation *theSimulation, surface_mask *upperSurfaceMask, surface_mask *lowerSurfaceMask, struct atom_rec **at, int *nat, char ***seq, int *nseq, int **seq_at_array,  char ***patches, ion_add **ions, int *nions, struct aa_build_data *buildData, int *build_type );
	void move_inside(void);
	void move_outside(void);	
	virtual int getIFITM3Index() { return IFITM3_BASE_AH; };
};

struct P2 : ifitm3
{
	virtual int getIFITM3Index() { return IFITM3_P2; };
};

struct P2_w60a : P2
{
	virtual int getIFITM3Index() { return IFITM3_P2_W60A; };
};

struct P2_f63q : P2
{

	virtual int getIFITM3Index() { return IFITM3_P2_F63Q; };
};

struct P2_f67q : P2
{
	virtual int getIFITM3Index() { return IFITM3_P2_F67Q; };
};


struct dimer : pcomplex
{
	double bond_length;
	double bond_k;
	void init( Simulation *theSimulation, surface *, double *rsurf, int f, double u, double v, int nmer=1 ); 
	void init( double *r );
	void loadParams( parameterBlock *block );
	
	int getNBonds( void );
	virtual void putBonds( int *bond_list, double *bond_r = NULL, double *bond_k = NULL );
	
	 void bind( int f, double u, double v );
	 void unbind( void );

	virtual double V( Simulation *theSimulation );
	virtual double grad( Simulation*theSimulation, double *surface_g, double *particle_g );
};

struct MAB : dimer
{
	// dtheta_0 is the optimal angle between the normals, projected along the vector between the particles
	double dtheta_0;
	// k_theta is the force constant for dtheta_0.
	double k_theta;

	void move_inside( void );
	void move_outside( void );

	void loadParams( parameterBlock * block );
	virtual double V( Simulation *theSimulation );
	virtual double grad( Simulation *theSimulation, double *surface_g, double *particle_g );
};

struct elasticCrowder : pcomplex
{
	double hard_r;
	double d; // the distance from the surface to the sphere center.
	double bond_k;
	double att_eps_val;
	double att_sigma_val;
	double att_sigma_inner;
	double crowder_mass;

	void init( Simulation *theSimulation, surface *, double *rsurf, int f, double u, double v, int nmer=1 ); 
	void loadParams( parameterBlock * block );
	int isElastic(void) { return 1; }	
	void move_inside( void );
	void move_outside( void );
	int getNBonds( void ) { return 1; }
	virtual void putBonds( int *bond_list, double *bond_r = NULL, double *bond_k = NULL );
	virtual double V( Simulation *theSimulation );
	virtual double grad( Simulation *theSimulation, double *surface_g, double *particle_g );
};

void propagateSolutionParticles( Simulation *theSimulation, double dt );
pcomplex *loadComplex( const char *name, const char *mod );

#endif
