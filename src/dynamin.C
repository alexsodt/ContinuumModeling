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

//#define DEBUG_DYNAMIN

static int dihe_loop_only = 14;

static double bond_force = 1; //1
static double angle_force = 500; // 100
static double dihedral_force =  50;
static int aq_sites_per_nmer = 5;

#define SQRT_SAFETY (1e-7)

#define EPS_SMALL (1e-14)
#define WORKING

static double k_nrm = 50;
static double syt7_phi0 = 0; // M_PI/3; //M_PI;
static double bilayer_thickness = 50.0;
static double PH_width = 30.0;
static double PH_to_stalk = 30.0;
static double inner_to_outer = bilayer_thickness/2 + PH_width/2 + PH_to_stalk; 

static double dynamin_stride = 6.35;
static double dynamin_rotor  = 23.68; // the three-sites have to be this distance apart.
static double dynamin_R         = 85; // the approximate radius of the ideal dynamin tube
static double dynamin_outer_R = dynamin_R + inner_to_outer;
static double dynamin_angle = 100; // the angle that one dynamin dimer takes up.
static double dynamin_axis_angle = 10.0; // just a guess, replace with correct value.
static double axis_stride = 40.0; // approx.


static double monomer_MW = 1000 * 200000; //amu per dimer.
static double site_radius = 10.0;

// this is determined by making the sites line up through one rotation.

double elbow_1_to_membrane = 65.0;
double elbow_angle	 = 120.0;	// check
double elbow_to_elbow    = 70.0;
double plane_angle       = 21.094992;	// check
double el_el_G_angle     = 70.596031;	// check
double elbow_to_G	 = 84.113721;	// check
double el_G_el_angle  = (161.273414);	// check
double final_dihe	 = 11.862716;	// check

void sdynamin::getParams( void )
{
	PH_width = 30.0;
	PH_to_stalk = 30.0;
	inner_to_outer = bilayer_thickness/2 + PH_width/2 + PH_to_stalk;
	dynamin_stride = 14.63;
	dynamin_rotor  = 26.14;
	dynamin_angle = 100;
	dynamin_axis_angle = 10.0;

	dynamin_R         = 65; // the approximate radius of the ideal dynamin tube
	dynamin_outer_R = dynamin_R + inner_to_outer;
	
	elbow_1_to_membrane = 65.0;
	elbow_angle	 = 122.0;
	elbow_to_elbow    = 70.0;
	plane_angle       = 29.742715;
	el_el_G_angle     = 73.292420; // check
	elbow_to_G	 = 86.032050;
	el_G_el_angle  = 157.978543;
	final_dihe	 = 18.720648;
}

void dynamin::getParams( void )
{
	PH_width = 30.0;
	PH_to_stalk = 30.0;
	inner_to_outer = bilayer_thickness/2 + PH_width/2 + PH_to_stalk;
	dynamin_stride = 6.35;
	dynamin_rotor  = 23.68;
	dynamin_angle = 100;
	dynamin_axis_angle = 10.0;

	dynamin_R         = 85; // the approximate radius of the ideal dynamin tube
	dynamin_outer_R = dynamin_R + inner_to_outer;

	elbow_1_to_membrane = 65.0;
	elbow_angle	 = 120.0;
	elbow_to_elbow    = 70.0;
	plane_angle       = 21.094992;
	el_el_G_angle     = 70.596031;
	elbow_to_G	 = 84.113721;
	el_G_el_angle  = (161.273414);
	final_dihe	 = 11.862716;
} 

void dynamin::init( double *r, int nmer )
{
	// aqueous initialization.

	nmer_saved =nmer;

	printf("Currently dynamin can only be initialized on the membrane.\n");
	exit(1);

	base_init();

	nattach = 2 * nmer;
	nsites = nattach + aq_sites_per_nmer * nmer;

	alloc();

	double vdw_r = site_radius;

	for( int s = 0; s < nsites; s++ )
	{
		sigma[s] = vdw_r;
		mass[s] = monomer_MW;
	}

	bound = 0;
	
	// initial geometry.

}

// initialize the BAR domain on the membrane.


// this is just above the PHD:


void dynamin::clone( Simulation *theSimulation, surface *theSurface, double *rsurf, pcomplex *takeFrom, int add_to)
{
	getParams();

	int nmer_copy = takeFrom->nmer_saved;

	double io_sign = 1;
	if( is_inside )
		io_sign *= -1;

	if( add_to < 0 )
		nmer_copy += add_to;

#ifdef OOPS //delete me
	base_init();
	
	nmer_saved = takeFrom->nmer_saved + add_to;
	int nmer = nmer_saved;

	nattach = 2 * nmer;
	nsites = nattach + aq_sites_per_nmer * nmer;

	alloc();

	double vdw_r = site_radius;

	for( int s = 0; s < nsites; s++ )
	{
		sigma[s] = vdw_r;
		mass[s] = monomer_MW;
	}

	for( int a = 0; a < nattach; a++ )
		sid[a] = theSurface->surface_id;

	bound = 1;

	// we get the attachment point of one connection, then place the other site by attempting to wrap around the positive curvature direction.

	int init = 1;


	double p_origin[3]={0,0,0};
	
	double **M;
	int mlow,mhigh;

	getM( &M, &mlow, &mhigh );
#endif
	
	double dynamin_axis[3] = { 0,0,0}; // the current best estimate of the dynamin axis.
	

	for( int m = 0; m < nmer_copy; m++ )
	{
		memcpy( puv + m * 4, takeFrom->puv+ m * 4, sizeof(double) * 4 );
		memcpy( fs + 2 * m, takeFrom->fs + 2 * m, sizeof(int) * 2 );

		memcpy( rall           + nattach*3           + m * 3 * aq_sites_per_nmer, 
			takeFrom->rall + takeFrom->nattach*3 + m * 3 * aq_sites_per_nmer,
				sizeof(double) * 3 * aq_sites_per_nmer ); 

		double r_1[3], nrm_1[3];
		double r_2[3], nrm_2[3];

		theSurface->evaluateRNRM( fs[2*m+0], puv[4*m+0], puv[4*m+1], r_1, nrm_1, rsurf );	
		nrm_1[0] *= io_sign;
		nrm_1[1] *= io_sign;
		nrm_1[2] *= io_sign;
		theSurface->evaluateRNRM( fs[2*m+1], puv[4*m+2], puv[4*m+3], r_2, nrm_2, rsurf );	
		nrm_2[0] *= io_sign;
		nrm_2[1] *= io_sign;
		nrm_2[2] *= io_sign;
		
		double dr_u[3];
		theSurface->ru(fs[2*m+0], puv[4*m+0], puv[4*m+1], rsurf, dr_u );
		double drdu = normalize(dr_u);
		double dr_v[3];
		theSurface->rv(fs[2*m+0], puv[4*m+0], puv[4*m+1], rsurf, dr_v );
		double drdv = normalize(dr_v);
			
		double c_vec_1[3];
		double c_vec_2[3];
		double c_val1, c_val2;
		double k;
		theSurface->c( fs[2*m+0], puv[4*m+0], puv[4*m+1], rsurf, &k, c_vec_1, c_vec_2, &c_val1, &c_val2 ); 

		double positive_dir_uv[2];
		double negative_dir_uv[2];
		double pos_c;
		if( c_val1*io_sign > c_val2*io_sign )
		{
			pos_c = c_val1;
			memcpy( positive_dir_uv, c_vec_1, sizeof(double) * 2 );
			memcpy( negative_dir_uv, c_vec_2, sizeof(double) * 2 );
		}
		else
		{
			pos_c = c_val2;
			memcpy( positive_dir_uv, c_vec_2, sizeof(double) * 2 );
			memcpy( negative_dir_uv, c_vec_1, sizeof(double) * 2 );
		}
	
		double d_axis[2] = { negative_dir_uv[0], negative_dir_uv[1] };
		double d_axis_expec_1[3] = {0,0,0};
		
		d_axis_expec_1[0] += d_axis[0] * dr_u[0]*drdu;
		d_axis_expec_1[1] += d_axis[0] * dr_u[1]*drdu;
		d_axis_expec_1[2] += d_axis[0] * dr_u[2]*drdu;
		
		d_axis_expec_1[0] += d_axis[1] * dr_v[0]*drdv;
		d_axis_expec_1[1] += d_axis[1] * dr_v[1]*drdv;
		d_axis_expec_1[2] += d_axis[1] * dr_v[2]*drdv;

		double l_axis_expec = normalize(d_axis_expec_1);
		
		double dr_orient[3] = { r_2[0] - r_1[0], r_2[1] - r_1[1], r_2[2] - r_1[2] };
	
		double temp[3];
		cross( dr_orient, nrm_1, temp );

		double dp = dot( d_axis_expec_1, temp );

		if( dp < 0 )
		{
			d_axis_expec_1[0] *= -1;
			d_axis_expec_1[1] *= -1;
			d_axis_expec_1[2] *= -1;
		}

		double w = 0.5;
	
		if( m == 0 )
			memcpy( dynamin_axis, d_axis_expec_1, sizeof(double) * 3 );
		else
		{
			double dp1=1, dp2;
			dp1 = dot( d_axis_expec_1, dynamin_axis );
			
			if( dp1 < 0 )
			{
				d_axis_expec_1[0] *= -1;
				d_axis_expec_1[1] *= -1;
				d_axis_expec_1[2] *= -1;
			}
			dynamin_axis[0] = w * dynamin_axis[0] + (1-w) * d_axis_expec_1[0];
			dynamin_axis[1] = w * dynamin_axis[1] + (1-w) * d_axis_expec_1[1];
			dynamin_axis[2] = w * dynamin_axis[2] + (1-w) * d_axis_expec_1[2];
		}
	}	
	
	for( int m = nmer_copy; m < nmer_saved; m++ )
		placeSubunit( theSimulation, theSurface, rsurf, m, dynamin_axis );
	
	for( int a = nmer_copy*2; a < nsites; a++ )
	{
		PBC_ext[3*a+0] = 0;
		PBC_ext[3*a+1] = 0;
		PBC_ext[3*a+2] = 0;
	}
	
	setrall(theSimulation);
	
	int nbonds = getNBonds();
	int theBonds[2*nbonds];

	putBonds(theBonds);

	int lock[nsites];
	memset( lock, 0, sizeof(int) * nsites);
	lock[0] = 1;


//#define DEBUG_BLOCK
#ifdef DEBUG_BLOCK
	for( int a = 0; a < nsites; a++ )
	{
		if( a < nattach )
			printf("C %lf %lf %lf\n", rall[3*a+0], rall[3*a+1], rall[3*a+2] );
		else
			printf("O %lf %lf %lf\n", rall[3*a+0], rall[3*a+1], rall[3*a+2] );
	}
	exit(1);
#endif

	for( int s = 0; s < nattach; s++ )
	{
		double nrmx[3];
	//	theSurface->evaluateRNRM( fs[s], puv[2*s+0], puv[2*s+1], rall+3*s, nrmx, rsurf ); 
	}

	int done = 0;
	while( !done )
	{
		done = 1;
		for( int b = 0; b < nbonds; b++ )
		{
			int a1 = theBonds[2*b+0];
			int a2 = theBonds[2*b+1];
	
			if( lock[a2] && ! lock[a1] )
			{
				int t = a2;
				a2 = a1;
				a1 = t;
			}
	
			if( lock[a2] ) continue;
			if( !lock[a1] ) continue;

	
	
			double dr_r[3] = { rall[3*a2+0] - rall[3*a1+0], rall[3*a2+1] - rall[3*a1+1], rall[3*a2+2] - rall[3*a1+2] };
	
			double put[3];
	

			MinImage3D( dr_r, theSimulation->PBC_vec, put, theSimulation->alpha );
			
	
			if( a2 >= nattach )
			{
				rall[3*a2+0] = rall[3*a1+0] + dr_r[0];
				rall[3*a2+1] = rall[3*a1+1] + dr_r[1];
				rall[3*a2+2] = rall[3*a1+2] + dr_r[2];
	
				PBC_ext[3*a2+0] = 0;
				PBC_ext[3*a2+1] = 0;
				PBC_ext[3*a2+2] = 0;
			}	
			else
			{
				PBC_ext[3*a2+0] += put[0];
				PBC_ext[3*a2+1] += put[1];
				PBC_ext[3*a2+2] += put[2];
			}

			done = 0;	
			lock[a2] = 1;
		}
	}

	refresh(theSimulation);
}

void dynamin::placeSubunit( Simulation *theSimulation, surface *theSurface, double *rsurf, int m, double dynamin_axis[3]  )
{
	getParams();

	double r_1[3]={1e10,1e10,1e10};
	double r_2[3]={1e10,1e10,1e10};

	double nrm_1[3]={1e10,1e10,1e10};
	double nrm_2[3];

	// pick the next spot	

	double uv1[2], duv1[2];
	double dr_u[3];
	double dr_v[3];
//	double uv1[2]={1e10,1e10}, duv1[2]={1e10,1e10};
//	double dr_u[3]={1e10,1e10,1e10};
//	double dr_v[3]={1e10,1e10,1e10};
	int f_1=-1, nf=-1;
	double **M=NULL;
	int mlow=0,mhigh=0;

	getM( &M, &mlow, &mhigh );
		
	double io_sign = 1;
	if( is_inside )
		io_sign *= -1;

	if( m == 0 )
	{
	}
	else
	{
		double the_sign = 1;
		double the_sign2 = -1;
		double use_axis[3] = { dynamin_axis[0] * the_sign,
				       dynamin_axis[1] * the_sign,
				       dynamin_axis[2] * the_sign };

		for( int subp = 0; subp < 2; subp++ )
		{
			int p_f = fs[2*(m-1)+subp];
			double p_u = puv[4*(m-1)+2*subp+0],
				p_v = puv[4*(m-1)+2*subp+1];
		
			double r_start[3], n_start[3];
			theSurface->evaluateRNRM( p_f, p_u, p_v, r_start, n_start, rsurf );
			n_start[0] *= io_sign;
			n_start[1] *= io_sign;
			n_start[2] *= io_sign;
			double dr_u[3];
			theSurface->ru( p_f, p_u, p_v, rsurf, dr_u );
			double dr_v[3];
			theSurface->rv( p_f, p_u, p_v, rsurf, dr_v );
			
			// direction that aligns with the current dynamin axis.
			
			double duv1[2]; 
			best_align( duv1, dr_u, dr_v, use_axis );
		
			duv1[0] *= dynamin_stride;
			duv1[1] *= dynamin_stride;
		
			int f_1 = p_f, nf = p_f;
			double uv1[2] = { p_u, p_v };
			
			do {
				f_1 = nf;
				nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf ); 
			} while( nf != f_1 );
			
			// now for the angular move.	
		
			double rotor_stride = (2*M_PI*dynamin_R) * (dynamin_rotor / 360.0);
		
			double r_pt[3], npt[3];
		
			theSurface->evaluateRNRM(f_1, uv1[0], uv1[1], r_pt, npt, rsurf );
			npt[0] *= io_sign;
			npt[1] *= io_sign;
			npt[2] *= io_sign;
			double cdr[3] = { r_pt[0] - r_start[0], r_pt[1] - r_start[1], r_pt[2] - r_start[2] };
		
			double put_discard[3];
			MinImage3D( cdr, theSimulation->PBC_vec, put_discard, theSimulation->alpha );
		
		
			double perp_dir[3];
			cross( use_axis, npt, perp_dir );
	
			theSurface->ru( f_1, uv1[0], uv1[1], rsurf, dr_u );
			theSurface->rv( f_1, uv1[0], uv1[1], rsurf, dr_v );
		
			best_align( duv1, dr_u, dr_v, perp_dir );
		
			duv1[0] *= the_sign2*rotor_stride;
			duv1[1] *= the_sign2*rotor_stride;
		
			do {
				f_1 = nf;
				nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf ); 
			} while( nf != f_1 );
			  
			fs[2*m+subp] = f_1;
		
			puv[4*m+2*subp+0] = uv1[0];
			puv[4*m+2*subp+1] = uv1[1];
			
			double r_final[3], n_final[3];
			theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], r_final, n_final, rsurf );
			n_final[0] *= io_sign;
			n_final[1] *= io_sign;
			n_final[2] *= io_sign;
		}
	}

	if( m == 0 )
	{
		theSurface->evaluateRNRM( fs[2*m+0], puv[4*m+0], puv[4*m+1], r_1, nrm_1, rsurf );	
	
		nrm_1[0] *= io_sign;
		nrm_1[1] *= io_sign;
		nrm_1[2] *= io_sign;
	
		theSurface->ru(fs[2*m+0], puv[4*m+0], puv[4*m+1], rsurf, dr_u );
		double drdu = normalize(dr_u);
		theSurface->rv(fs[2*m+0], puv[4*m+0], puv[4*m+1], rsurf, dr_v );
		double drdv = normalize(dr_v);
			
		double c_vec_1[3];
		double c_vec_2[3];
		double c_val1, c_val2;
		double k;
		theSurface->c( fs[2*m+0], puv[4*m+0], puv[4*m+1], rsurf, &k, c_vec_1, c_vec_2, &c_val1, &c_val2 ); 
	
		double positive_dir_uv[2];
		double negative_dir_uv[2];
		double pos_c;
		if( c_val1*io_sign > c_val2*io_sign )
		{
			pos_c = c_val1;
			memcpy( positive_dir_uv, c_vec_1, sizeof(double) * 2 );
			memcpy( negative_dir_uv, c_vec_2, sizeof(double) * 2 );
		}
		else
		{
			pos_c = c_val2;
			memcpy( positive_dir_uv, c_vec_2, sizeof(double) * 2 );
			memcpy( negative_dir_uv, c_vec_1, sizeof(double) * 2 );
		}
		
		
		
		double d_axis[2] = { negative_dir_uv[0], negative_dir_uv[1] };
		double d_axis_expec_1[3] = {0,0,0};
		
		d_axis_expec_1[0] += d_axis[0] * dr_u[0]*drdu;
		d_axis_expec_1[1] += d_axis[0] * dr_u[1]*drdu;
		d_axis_expec_1[2] += d_axis[0] * dr_u[2]*drdu;
		
		d_axis_expec_1[0] += d_axis[1] * dr_v[0]*drdv;
		d_axis_expec_1[1] += d_axis[1] * dr_v[1]*drdv;
		d_axis_expec_1[2] += d_axis[1] * dr_v[2]*drdv;
		
		double l_axis_expec = normalize(d_axis_expec_1);
		
		
		if( m == 0 )
		{
			memcpy( dynamin_axis, d_axis_expec_1, sizeof(double) * 3 );
		}
		else
		{
			double w = 0.5;
			double check_dp = dot( dynamin_axis, d_axis_expec_1 );
		
			if( check_dp < 0 )
			{
				d_axis[0] *=-1;
				d_axis[1] *=-1;
			
				d_axis_expec_1[0] *= -1;
				d_axis_expec_1[1] *= -1;
				d_axis_expec_1[2] *= -1;
			}
	
	
			dynamin_axis[0] = dynamin_axis[0] * w + (1-w) * d_axis_expec_1[0];
			dynamin_axis[1] = dynamin_axis[1] * w + (1-w) * d_axis_expec_1[1];
			dynamin_axis[2] = dynamin_axis[2] * w + (1-w) * d_axis_expec_1[2];
		}
	

		double elbow_1_A[3] = { r_1[0] +  nrm_1[0] * elbow_1_to_membrane,
				        r_1[1] +  nrm_1[1] * elbow_1_to_membrane,
				        r_1[2] +  nrm_1[2] * elbow_1_to_membrane };
	
		double elbow_2_A[3] = { elbow_1_A[0] +  nrm_1[0] * elbow_to_elbow,			
					elbow_1_A[1] +  nrm_1[1] * elbow_to_elbow,
					elbow_1_A[2] +  nrm_1[2] * elbow_to_elbow };
	
		// bent straight across the tube by about 60 degrees:
		rotateArbitrary( elbow_2_A, dynamin_axis, elbow_1_A, 1, (180.0-elbow_angle)*(M_PI/180.0) );	
	
		// we now bend "back" along the dynamin axis by about 21 degrees. 
	
		double d_elbow[3] = { elbow_2_A[0] - elbow_1_A[0], elbow_2_A[1] - elbow_1_A[1], elbow_2_A[2] - elbow_1_A[2] };
		normalize(d_elbow);
	
		double G_center[3] = { elbow_2_A[0] + d_elbow[0] * elbow_to_G,
				             elbow_2_A[1] + d_elbow[1] * elbow_to_G,
					     elbow_2_A[2] + d_elbow[2] * elbow_to_G };
		//rotateArbitrary( G_center, dynamin_axis, elbow_2_A, 1, (180.0 - el_el_G_angle) * (M_PI/180.0) ); 
		rotateArbitrary( G_center, dynamin_axis, elbow_2_A, 1, (el_el_G_angle) * (M_PI/180.0) ); 
		// the dihedral rotation:
		rotateArbitrary( G_center, d_elbow, elbow_2_A, 1, -plane_angle * (M_PI/180.0) ); 
		
		// the off-axis plane:
		double dr_G_el[3] = { G_center[0] - elbow_2_A[0], G_center[1] - elbow_2_A[1], G_center[2] - elbow_2_A[2] };
		normalize( dr_G_el);
		double off_axis[3];
		cross( d_elbow, dr_G_el, off_axis );
		normalize( off_axis);
	
		double elbow_2_B[3] = { G_center[0] + dr_G_el[0] * elbow_to_G,
					G_center[1] + dr_G_el[1] * elbow_to_G,
					G_center[2] + dr_G_el[2] * elbow_to_G };
		rotateArbitrary( elbow_2_B, off_axis, G_center, 1, (180.0 - el_G_el_angle) * (M_PI/180.0) ); 		
	
		double b_G_el[3] = { elbow_2_B[0] - G_center[0],
				     elbow_2_B[1] - G_center[1],
				     elbow_2_B[2] - G_center[2] };
		normalize(b_G_el);
	
		double elbow_1_B[3] = { elbow_2_B[0] + b_G_el[0] * elbow_to_elbow,
					elbow_2_B[1] + b_G_el[1] * elbow_to_elbow,
					elbow_2_B[2] + b_G_el[2] * elbow_to_elbow };
	
		rotateArbitrary( elbow_1_B, off_axis, elbow_2_B, 1, (180.0-elbow_angle) * (M_PI/180.0) );	
		
		rotateArbitrary( elbow_1_B, b_G_el, elbow_2_B, 1, final_dihe * (M_PI/180.0) );	
	
		double elb_to_membrane_B[3] = { 
				elbow_1_B[0] - elbow_2_B[0],
				elbow_1_B[1] - elbow_2_B[1],
				elbow_1_B[2] - elbow_2_B[2] };
	
		normalize(elb_to_membrane_B);
					
	
		double expec_membrane_pt[3] =
		{
			elbow_1_B[0] + elb_to_membrane_B[0] * elbow_1_to_membrane, 
			elbow_1_B[1] + elb_to_membrane_B[1] * elbow_1_to_membrane, 
			elbow_1_B[2] + elb_to_membrane_B[2] * elbow_1_to_membrane
		};
	
		rotateArbitrary( expec_membrane_pt, off_axis, elbow_1_B, 1, (180.0-elbow_angle) * (M_PI/180.0) );	
	
		if( m == 0 )
		{
			double del[3] = { expec_membrane_pt[0] - r_1[0], expec_membrane_pt[1] - r_1[1], expec_membrane_pt[2] - r_1[2] };
			
			double du_move[2]; 
				
			theSurface->ru( fs[2*m+0], puv[4*m+0], puv[4*m+1], rsurf, dr_u );
			theSurface->rv( fs[2*m+0], puv[4*m+0], puv[4*m+1], rsurf, dr_v );
			best_align( du_move, dr_u, dr_v, del, 0 );
		 
			f_1 = fs[2*m+0];
			nf = fs[2*m+0];
			uv1[0] = puv[4*m+0];
			uv1[1] = puv[4*m+1];
			duv1[0] = du_move[0];
			duv1[1] = du_move[1];
			do {
				f_1 = nf;
				nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf ); 
		
			} while( nf != f_1 );
		
			fs[2*m+1] = f_1;
			puv[4*m+2] = uv1[0];
			puv[4*m+3] = uv1[1];	
			
			theSurface->evaluateRNRM( fs[2*m+1], puv[4*m+2], puv[4*m+3], r_2, nrm_2, rsurf );	
			nrm_2[0] *= io_sign;
			nrm_2[1] *= io_sign;
			nrm_2[2] *= io_sign;
		}
	
	
		memcpy( rall+3*nattach+m*3*(aq_sites_per_nmer)+0*3, elbow_1_A, sizeof(double)*3);
		memcpy( rall+3*nattach+m*3*(aq_sites_per_nmer)+1*3, elbow_2_A, sizeof(double)*3);
		memcpy( rall+3*nattach+m*3*(aq_sites_per_nmer)+2*3, G_center, sizeof(double)*3);
		memcpy( rall+3*nattach+m*3*(aq_sites_per_nmer)+3*3, elbow_2_B, sizeof(double)*3);
		memcpy( rall+3*nattach+m*3*(aq_sites_per_nmer)+4*3, elbow_1_B, sizeof(double)*3);
	
	#define ROTOR_PLACE
	#ifdef ROTOR_PLACE
		// rotate around the first point to be consistent with the second membrane-bound pt.
	
		double tpts[3*aq_sites_per_nmer];
	
		double *put = rall+3*nattach+m*3*(aq_sites_per_nmer);
		memcpy( tpts, put, sizeof(double)*3*aq_sites_per_nmer );
	
		double rotor_min = - M_PI/2;
		double rotor_max = M_PI/2;
	
		double *lastp = put + (aq_sites_per_nmer-1)*3;
	
		double best_pot=1e20;
		double best_rot=0;
	
		int n_rotor_check = 50;
	
		for( int ir = 0; ir <= n_rotor_check; ir++ )
		{
			double rot = rotor_min+ (rotor_max-rotor_min)*ir / (double)n_rotor_check;
			memcpy(put, tpts, sizeof(double) * 3 * aq_sites_per_nmer );
			rotateArbitrary( put, off_axis, tpts, 5, rot );
	
			double dr[3] = { lastp[0] - r_2[0], lastp[1]-r_2[1], lastp[2]-r_2[2] };
	
			double l = normalize(dr);
	
			double dp = dr[0] * nrm_2[0] + dr[1] * nrm_2[1] + dr[2] * nrm_2[2];
			
			double dp_target = 1;	
			if( is_inside )
				dp_target *= -1;
	
			double pot = k_nrm/2 * (dp-dp_target)*(dp-dp_target) + bond_force/2 * pow(l - elbow_1_to_membrane,2.0); 
	
			if( pot < best_pot )
			{
//				printf("bestPot: %lf rot: %lf\n", best_pot, rot );
				best_pot = pot;
				best_rot = rot;
			}
//			else
//				printf("less bestPot: %lf rot: %lf\n", best_pot, rot );
		}
			
		memcpy(put, tpts, sizeof(double) * 3 * aq_sites_per_nmer );
		rotateArbitrary( put, off_axis, tpts, aq_sites_per_nmer, best_rot );
	}
	else
	{

		double *prev = rall + nattach*3 + (m-1)*aq_sites_per_nmer*3;
		double *next = rall + nattach*3 + (m)*aq_sites_per_nmer*3;

		memcpy( next, prev, sizeof(double) * 3 * aq_sites_per_nmer );

		double p_r_1[3], p_r_2[3];
		double p_n_1[3], p_n_2[3];

		theSurface->evaluateRNRM( fs[2*(m-1)+0], puv[4*(m-1)+0], puv[4*(m-1)+1], p_r_1, p_n_1, rsurf );
		theSurface->evaluateRNRM( fs[2*(m-1)+1], puv[4*(m-1)+2], puv[4*(m-1)+3], p_r_2, p_n_2, rsurf );
		double av_nrm[3] = { p_n_1[0]+p_n_2[0], p_n_1[1]+p_n_2[1], p_n_1[2]+p_n_2[2] };
		normalize(av_nrm);

		double dynamin_origin[3] = {prev[3*(aq_sites_per_nmer/2)+0] - av_nrm[0] * dynamin_outer_R,
					    prev[3*(aq_sites_per_nmer/2)+1] - av_nrm[1] * dynamin_outer_R,
					    prev[3*(aq_sites_per_nmer/2)+2] - av_nrm[2] * dynamin_outer_R };
		dynamin_origin[0] += dynamin_axis[0] * dynamin_stride;
		dynamin_origin[1] += dynamin_axis[1] * dynamin_stride;
		dynamin_origin[2] += dynamin_axis[2] * dynamin_stride;
		rotateArbitrary( next, dynamin_axis, dynamin_origin, aq_sites_per_nmer, -(M_PI/180.0) * dynamin_rotor );
		
		for( int aq = 0; aq < aq_sites_per_nmer; aq++ )
		{
			rall[nattach*3+m*aq_sites_per_nmer+aq*3+0] += dynamin_axis[0] * dynamin_stride;
			rall[nattach*3+m*aq_sites_per_nmer+aq*3+1] += dynamin_axis[1] * dynamin_stride;
			rall[nattach*3+m*aq_sites_per_nmer+aq*3+2] += dynamin_axis[2] * dynamin_stride;
		} 
	}

#endif
	
}

void dynamin::init( Simulation *theSimulation, surface *theSurface, double *rsurf, int f, double u, double v, int nmer)
{
	getParams();

	// assume for now this is one of the points on the membrane neck.
#if 0
	// hack to start, get the site to align with the dynamin pdb.
	
	int best_f = -1;
	double best_uv[2];
	double best_chi2 = 1e10;	
	double r_target[3] = { 317.808973045822, 400.286776280324, 305.803188679245 };

	for( int i = 0; i < 10000; i ++)
	{

		theSurface->randomPointOnSurface( &f, &u, &v );			
	
		double rp[3], np[3];
		theSurface->evaluateRNRM( f, u, v, rp, np, rsurf );

		double dr[3] = { rp[0] - r_target[0], rp[1] - r_target[1], rp[2] - r_target[2] };

		double l = normalize(dr);

		if( l < best_chi2 )
		{
			best_f = f;
			best_uv[0] = u;
			best_uv[1] = v;
			best_chi2 = l;
		}
	}
	f = best_f;
	u = best_uv[0];
	v = best_uv[1];

	theSurface->find_near_spot( &f, &u, &v, r_target, rsurf ); 
#endif

	base_init();
	
	nmer_saved =nmer;

	nattach = 2 * nmer;
	nsites = nattach + aq_sites_per_nmer * nmer;

	alloc();

	double vdw_r = site_radius;

	for( int s = 0; s < nsites; s++ )
	{
		sigma[s] = vdw_r;
		mass[s] = monomer_MW;
	}

	for( int a = 0; a < nattach; a++ )
		sid[a] = theSurface->surface_id;

	bound = 1;

	double io_sign = 1;
	if( is_inside )
		io_sign *= -1;
	// we get the attachment point of one connection, then place the other site by attempting to wrap around the positive curvature direction.


	int init = 1;

	double dynamin_axis[3] = { 0,0,0}; // the current best estimate of the dynamin axis.

	double p_origin[3]={0,0,0};
	
	double **M;
	int mlow,mhigh;

	getM( &M, &mlow, &mhigh );

	for( int m = 0; m < nmer; m++ )
	{
		if( m == 0 )
		{
			fs[0] = f;
			puv[0] = u;			
			puv[1] = v;
		}

		placeSubunit( theSimulation, theSurface, rsurf, m, dynamin_axis );

		init = 0;		
	}

	memcpy( grad_fs, fs, sizeof(int) * nattach );
	memcpy( grad_puv, puv, sizeof(double) * 2*nattach );

	memset( PBC_ext, 0, sizeof(double) * 3*nsites );
	
	setrall(theSimulation);
	
	int nbonds = getNBonds();
	int theBonds[2*nbonds];

	putBonds(theBonds);

	int lock[nsites];
	memset( lock, 0, sizeof(int) * nsites);
	lock[0] = 1;

	for( int a = 0; a < nsites; a++ )
	{
		PBC_ext[3*a+0] = 0;
		PBC_ext[3*a+1] = 0;
		PBC_ext[3*a+2] = 0;
	}

	int done = 0;
	while( !done )
	{
		done = 1;
		for( int b = 0; b < nbonds; b++ )
		{
			int a1 = theBonds[2*b+0];
			int a2 = theBonds[2*b+1];
	
			if( lock[a2] && ! lock[a1] )
			{
				int t = a2;
				a2 = a1;
				a1 = t;
			}
	
			if( lock[a2] ) continue;
			if( !lock[a1] ) continue;
	
	
			double dr_r[3] = { rall[3*a2+0] - rall[3*a1+0], rall[3*a2+1] - rall[3*a1+1], rall[3*a2+2] - rall[3*a1+2] };
	
			double put[3];
	
			MinImage3D( dr_r, theSimulation->PBC_vec, put, theSimulation->alpha );
	
			if( a2 >= nattach )
			{
				rall[3*a2+0] = rall[3*a1+0] + dr_r[0];
				rall[3*a2+1] = rall[3*a1+1] + dr_r[1];
				rall[3*a2+2] = rall[3*a1+2] + dr_r[2];
	
				PBC_ext[3*a2+0] = 0;
				PBC_ext[3*a2+1] = 0;
				PBC_ext[3*a2+2] = 0;
			}	
			else
			{
				PBC_ext[3*a2+0] = put[0];
				PBC_ext[3*a2+1] = put[1];
				PBC_ext[3*a2+2] = put[2];
			}

			done = 0;	
			lock[a2] = 1;
		}
	}
}

// custom orient procedure.



void dynamin::bind( int f, double u, double v)
{
	bound = 1;
}

void dynamin::unbind( void )
{
	bound = 0;
}


void dynamin::loadParams( parameterBlock *block )
{
}

int inter_bonds_per_nmer = aq_sites_per_nmer-1;
int attach_per_nmer = 2;

int dynamin::getNBonds( void )
{
	
	int nsites_per_nmer = attach_per_nmer + aq_sites_per_nmer;
	
	int bonds_per_nmer = (aq_sites_per_nmer + attach_per_nmer - 1) + inter_bonds_per_nmer;

	int nbonds = bonds_per_nmer * nmer_saved - inter_bonds_per_nmer;

	return nbonds;
}

//#define MID_DIHE
#define SUPP_DIHE
	
int nangles_intra = attach_per_nmer + aq_sites_per_nmer - 2
#ifdef SUPP_DIHE
+ 2
#endif
;
int nangles_inter = inter_bonds_per_nmer; // one per bond 

int dynamin::getNAngles( void )
{
	int nangles = nangles_intra *nmer_saved + nangles_inter * (nmer_saved-1);


	return nangles;
}


void dynamin::putAngles( int *ang_list, double *angle_theta, double *angle_k )
{
	double use_angle_k = angle_force;


	int nangles = nangles_intra *nmer_saved + nangles_inter * (nmer_saved-1);

	// intra angles.
	
	int cnt = 0;
	for( int m = 0; m < nmer_saved; m++ )
	{
		int a0 = nattach + m * aq_sites_per_nmer;
		int s0 = attach_per_nmer * m;
// ang 1
		ang_list[3*cnt+0] = s0; 
		ang_list[3*cnt+1] = a0; 
		ang_list[3*cnt+2] = a0+1; 

		angle_theta[cnt] = 120.583149;
		angle_k[cnt] = use_angle_k;
	
		cnt++;
// ang 2	
		ang_list[3*cnt+0] = a0; 
		ang_list[3*cnt+1] = a0+1; 
		ang_list[3*cnt+2] = a0+2; 

		angle_theta[cnt] = 109.403969;
		angle_k[cnt] = use_angle_k;

		cnt++;
// ang 3	
		ang_list[3*cnt+0] = a0+1; 
		ang_list[3*cnt+1] = a0+2; 
		ang_list[3*cnt+2] = a0+3; 

		angle_theta[cnt] = 161.273414;
		angle_k[cnt] = use_angle_k;

		cnt++;
// ang 4	
		ang_list[3*cnt+0] = a0+2; 
		ang_list[3*cnt+1] = a0+3; 
		ang_list[3*cnt+2] = a0+4; 

		angle_theta[cnt] = 108.305316;
		angle_k[cnt] = use_angle_k;
		
		cnt++;
// ang 5
		ang_list[3*cnt+0] = a0+3; 
		ang_list[3*cnt+1] = a0+4; 
		ang_list[3*cnt+2] = s0+1; 

		angle_theta[cnt] = 108.305316;
		angle_k[cnt] = use_angle_k;
		
		cnt++;
#ifdef SUPP_DIHE
// ang S1
		ang_list[3*cnt+0] = a0; 
		ang_list[3*cnt+1] = a0+1; 
		ang_list[3*cnt+2] = a0+3; 

		angle_theta[cnt] = 100.0;
		angle_k[cnt] = use_angle_k;
		
		cnt++;

// ang S2
		ang_list[3*cnt+0] = a0+1; 
		ang_list[3*cnt+1] = a0+3; 
		ang_list[3*cnt+2] = a0+4; 

		angle_theta[cnt] = 100.0;
		angle_k[cnt] = use_angle_k;
		
		cnt++;


#endif


	}
	
	// inter angles.
	
	for( int m = 0; m < nmer_saved-1; m++ )
	{
		int a0 = nattach + m * aq_sites_per_nmer;
		int s0 = attach_per_nmer * m;
		
		int a1 = nattach + (m+1) * aq_sites_per_nmer;
		int s1 = attach_per_nmer * (m+1);
// ang 1
		ang_list[3*cnt+0] = a0+1; 
		ang_list[3*cnt+1] = a1; 
		ang_list[3*cnt+2] = a1+1; 

		angle_theta[cnt] = 69.990722;
		angle_k[cnt] = use_angle_k;
	
		cnt++;
// ang 2	
		ang_list[3*cnt+0] = a0+2; 
		ang_list[3*cnt+1] = a1+1; 
		ang_list[3*cnt+2] = a1+2; 

		angle_theta[cnt] = 73.935856;
		angle_k[cnt] = use_angle_k;

		cnt++;
// ang 3	
		ang_list[3*cnt+0] = a0+3; 
		ang_list[3*cnt+1] = a1+2; 
		ang_list[3*cnt+2] = a1+3; 

		angle_theta[cnt] = 54.940527; 
		angle_k[cnt] = use_angle_k;

		cnt++;
// ang 4	
		ang_list[3*cnt+0] = a0+4; 
		ang_list[3*cnt+1] = a1+3; 
		ang_list[3*cnt+2] = a1+4; 

		angle_theta[cnt] = 57.027801;
		angle_k[cnt] = use_angle_k;
		
		cnt++;
	}

	// convert to rads
	
	for( int a = 0; a < cnt; a++ )
		angle_theta[a] *= (M_PI/180.0);
}


int ndihe_inter = 4;

int ndihe_intra = 
2 + 

#ifdef MID_DIHE
2 +
#endif

#ifdef SUPP_DIHE
1
#else
0
#endif
;

int dynamin::getNDihe( void )
{
	int ndihe = ndihe_intra * nmer_saved + ndihe_inter * (nmer_saved-1);

	return ndihe;
}

void dynamin::putDihe( int *dihe_list, double *dihe_theta, double *dihe_k )
{
	double use_dihe_k = dihedral_force;

	int ndihe = getNDihe();

	// intra dihe.
	
	int cnt = 0;
	for( int m = 0; m < nmer_saved; m++ )
	{
		int a0 = nattach + m * aq_sites_per_nmer;
		int s0 = attach_per_nmer * m;
// dihe 1
		dihe_list[4*cnt+0] = s0; 
		dihe_list[4*cnt+1] = a0; 
		dihe_list[4*cnt+2] = a0+1; 
		dihe_list[4*cnt+3] = a0+2; 

		dihe_theta[cnt] = 26.080962*(M_PI/180.0);
		dihe_k[cnt] = use_dihe_k;

	
		cnt++;

#ifdef MID_DIHE
// dihe 2
		dihe_list[4*cnt+0] = a0; 
		dihe_list[4*cnt+1] = a0+1; 
		dihe_list[4*cnt+2] = a0+2; 
		dihe_list[4*cnt+3] = a0+3; 

		dihe_theta[cnt] = 4.219722*(M_PI/180.0);
		dihe_k[cnt] = use_dihe_k;
	
		cnt++;
// dihe 3
		dihe_list[4*cnt+0] = a0+1; 
		dihe_list[4*cnt+1] = a0+2; 
		dihe_list[4*cnt+2] = a0+3; 
		dihe_list[4*cnt+3] = a0+4; 

		dihe_theta[cnt] = 18.720648*(M_PI/180.0);
		dihe_k[cnt] = use_dihe_k;
	
		cnt++;
#endif

#ifdef SUPP_DIHE
		dihe_theta[cnt] = 22.032201 * (M_PI/180.0);
		dihe_list[4*cnt+0] = a0; 
		dihe_list[4*cnt+1] = a0+1; 

		dihe_list[4*cnt+2] = a0+3; 
		dihe_list[4*cnt+3] = a0+4; 
		dihe_k[cnt] = use_dihe_k;

		cnt++;
#endif
// dihe 4
		dihe_list[4*cnt+0] = a0+2; 
		dihe_list[4*cnt+1] = a0+3; 
		dihe_list[4*cnt+2] = a0+4; 
		dihe_list[4*cnt+3] = s0+1; 

		dihe_theta[cnt] = 11.862716*(M_PI/180.0);
		dihe_k[cnt] = use_dihe_k;
	
		cnt++;

	}
	
	// inter dihe.
	
	for( int m = 0; m < nmer_saved-1; m++ )
	{
		int a0 = nattach + m * aq_sites_per_nmer;
		int s0 = attach_per_nmer * m;
		
		int a1 = nattach + (m+1) * aq_sites_per_nmer;
		int s1 = attach_per_nmer * (m+1);
// dihe 1		
		dihe_list[4*cnt+0] = a0; 
		dihe_list[4*cnt+1] = a0+1; 
		dihe_list[4*cnt+2] = a1; 
		dihe_list[4*cnt+3] = a1+1; 

		dihe_theta[cnt] = 176.196137*(M_PI/180.0);
		dihe_k[cnt] = use_dihe_k;
	
		cnt++;

// dihe 2
		dihe_list[4*cnt+0] = a0+1;
		dihe_list[4*cnt+1] = a0+2; 
		dihe_list[4*cnt+2] = a1+1; 
		dihe_list[4*cnt+3] = a1+2; 

		dihe_theta[cnt] = 176.908217*(M_PI/180.0);
		dihe_k[cnt] = use_dihe_k;
	
		cnt++;

// dihe 3
		dihe_list[4*cnt+0] = a0+2;
		dihe_list[4*cnt+1] = a0+3; 
		dihe_list[4*cnt+2] = a1+2; 
		dihe_list[4*cnt+3] = a1+3; 

		dihe_theta[cnt] = 154.045015*(M_PI/180.0);
		dihe_k[cnt] = use_dihe_k;
	
		cnt++;

// dihe 3
		dihe_list[4*cnt+0] = a0+3;
		dihe_list[4*cnt+1] = a0+4; 
		dihe_list[4*cnt+2] = a1+3; 
		dihe_list[4*cnt+3] = a1+4; 

		dihe_theta[cnt] = 153.523730*(M_PI/180.0);
		dihe_k[cnt] = use_dihe_k;
	
		cnt++;
	}

	// convert to rads
	
//	for( int a = 0; a < cnt; a++ )
//		dihe_theta[a] *= (M_PI/180.0);
}

void dynamin::putBonds( int *bond_list, double *bond_r, double *bond_k_in)
{
	double bond_k = bond_force;
	int attach_per_nmer = 2;
	int nsites_per_nmer = attach_per_nmer + aq_sites_per_nmer;
	
	int bonds_per_nmer = (aq_sites_per_nmer + attach_per_nmer - 1) + aq_sites_per_nmer-1;

/*
double elbow_1_to_membrane = 45.0;
double elbow_angle	 = 120.0;
double elbow_to_elbow    = 70.0;
double plane_angle       = 21.094992;
double el_el_G_angle     = 70.596031;
double elbow_to_G	 = 84.113721;
double el_G_el_angle  = (161.273414);
double final_dihe	 = 11.862716;
*/
	for( int m = 0; m < nmer_saved; m++ )
	{
		// internal bonds.
		int a0 = nattach + m * aq_sites_per_nmer;
		int s0 = attach_per_nmer * m;
		int nb = 0;

		bond_list[(m*bonds_per_nmer)*2+nb*2+0] = s0;
		bond_list[(m*bonds_per_nmer)*2+nb*2+1] = a0;
		if( bond_r ) bond_r[m*bonds_per_nmer+nb] = elbow_1_to_membrane;
		if( bond_k_in ) bond_k_in[m*bonds_per_nmer+nb] = bond_k;
		nb++;
		
		bond_list[(m*bonds_per_nmer)*2+nb*2+0] = a0+aq_sites_per_nmer-1;
		bond_list[(m*bonds_per_nmer)*2+nb*2+1] = s0+1;
		if( bond_r ) bond_r[m*bonds_per_nmer+nb] = elbow_1_to_membrane;
		if( bond_k_in ) bond_k_in[m*bonds_per_nmer+nb] = bond_k;
		nb++;

		for( int aq = 0; aq < aq_sites_per_nmer-1; aq++ )
		{
			bond_list[(m*bonds_per_nmer)*2+nb*2+0] = a0+aq;
			bond_list[(m*bonds_per_nmer)*2+nb*2+1] = a0+aq+1;

			
			if( bond_r )  {
				switch( aq )
				{
					case 0:
					case 3: 
						bond_r[m*bonds_per_nmer+nb] = elbow_to_elbow;
						break;
					case 1:
					case 2: 
						bond_r[m*bonds_per_nmer+nb] = elbow_to_G;
						break;
					
				}
			}
			if( bond_k_in ) bond_k_in[m*bonds_per_nmer+nb] = bond_k;

			nb++;
		}

		// inter bonds.
		
		if( m < nmer_saved-1 )
		{
			int a0 = nattach + m * aq_sites_per_nmer;
			int a1 = nattach + (m+1) * aq_sites_per_nmer;

			for( int aq = 0; aq < aq_sites_per_nmer-1; aq++ )
			{
				bond_list[(m*bonds_per_nmer)*2+2*nb+0] = a0+aq+1;
				bond_list[(m*bonds_per_nmer)*2+2*nb+1] = a1+aq;
				if( bond_r )
				{
					switch( aq )
					{
						case 0:
							bond_r[m*bonds_per_nmer+nb] = 56.354136;
							break;
						case 1:
							bond_r[m*bonds_per_nmer+nb] = 47.392964;
							break;
						case 2:
							bond_r[m*bonds_per_nmer+nb] = 37.769028;
							break;
						case 3:
							bond_r[m*bonds_per_nmer+nb] = 38.011972;
							break;
					}
				}

				if( bond_k_in ) bond_k_in[m*bonds_per_nmer+nb] = bond_k;

				nb++;
			}
		}
	}

}

double dynamin::V( Simulation *theSimulation  )
{
	// normal potential? match 0-4 and 2-5 vectors to normal with angle potential

	double *alphas = theSimulation->alpha;

	double r[3*nsites];
	double n[3*nsites];

	int nbonds = getNBonds();
	int theBonds[2*nbonds];
	double bond_r_p[nbonds];
	double bond_k_p[nbonds];
	putBonds(theBonds, bond_r_p, bond_k_p );

	int nangles = getNAngles();
	int theAngles[3*nangles];
	double angle_th[nangles];
	double angle_k[nangles];
	putAngles(theAngles, angle_th, angle_k );
	
	int ndihe = getNDihe();
	int theDihe[4*ndihe];
	double dihe_th[ndihe];
	double dihe_k[ndihe];
	putDihe(theDihe, dihe_th, dihe_k );

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
		// zero nrm for aqueous
		memset(n+nattach*3 , 0, sizeof(double) * (nsites-nattach)*3 );
	}
	else
	{
		memcpy(r , rall, sizeof(double) * 3*nsites );
		memset(n, 0, sizeof(double) * 3*nsites ); 
	}


	double pot_nrm = 0, pot_bond =0, pot_angle = 0, pot_dihe = 0;

	for( int m = 0; m < nmer_saved; m++ )
	{
		int s0 = m * attach_per_nmer;
		int sx = m * attach_per_nmer+1;
		int a0 = nattach + m * aq_sites_per_nmer;
		int ax = nattach + m * aq_sites_per_nmer + aq_sites_per_nmer-1;

		// normal potential, first site:
		double dr[3] = { r[a0*3+0] - r[s0*3+0], r[a0*3+1] - r[s0*3+1], r[a0*3+2] - r[s0*3+2] };
		double ln = normalize(dr);
		double dp = dr[0] * n[s0*3+0] + dr[1] * n[s0*3+1] + dr[2] * n[s0*3+2];

		double dp_target = 1;	

		if( is_inside )
			dp_target *= -1;
		
		pot_nrm += 0.5 * k_nrm * (dp-dp_target)*(dp-dp_target);
		
		// normal potential, second site:

		dr[0] = r[3*ax+0] - r[3*(sx)+0];
		dr[1] = r[3*ax+1] - r[3*(sx)+1];
		dr[2] = r[3*ax+2] - r[3*(sx)+2];
		ln = normalize(dr);

		dp = dr[0] * n[sx*3+0] + dr[1] * n[sx*3+1] + dr[2] * n[sx*3+2];

		dp_target = 1;	

		if( is_inside )
			dp_target *= -1;
		
		pot_nrm += 0.5 * k_nrm * (dp-dp_target)*(dp-dp_target);
	}

	for( int b = 0; b < nbonds; b++ )
	{
		int a1 = theBonds[2*b+0];
		int a2 = theBonds[2*b+1];

		double dr[3] = { r[3*a1+0] - r[3*a2+0], r[3*a1+1] - r[3*a2+1], r[3*a1+2] - r[3*a2+2] };
		
		double r0 = bond_r_p[b];
		double k  = bond_k_p[b];

		double ln = normalize(dr);

		double l_pot = 0.5 * k * (ln - r0) * (ln-r0);
		pot_bond += l_pot;
	
	
	}

	// angles

	for( int a = 0; a < nangles; a++ )
	{
		int a1 = theAngles[3*a+0];
		int a2 = theAngles[3*a+1];
		int a3 = theAngles[3*a+2];

		int c1 = (a1-nattach) / aq_sites_per_nmer;
		int c2 = (a2-nattach) / aq_sites_per_nmer;
		int c3 = (a3-nattach) / aq_sites_per_nmer;


		double l_angle_k = angle_k[a];
		double phi0 = angle_th[a];
		double dr1[3] = { r[3*a1+0] - r[3*a2+0],
				 r[3*a1+1] - r[3*a2+1],
				 r[3*a1+2] - r[3*a2+2] };
		double dr2[3] = { r[3*a3+0] - r[3*a2+0],
				 r[3*a3+1] - r[3*a2+1],
				 r[3*a3+2] - r[3*a2+2] };
		double ln1 = normalize(dr1);
		double ln2 = normalize(dr2);

		double dp = dr1[0]*dr2[0]+dr1[1]*dr2[1]+dr1[2]*dr2[2];
		if( dp > 1 ) dp = 1 - 1e-13;
		double ac = acos(dp);

//#define DEBUG_DYNAMIN

#ifdef DEBUG_DYNAMIN
//		if( fabs(phi0-160*(M_PI/180.0)) < 5*(M_PI/180.0) )
			printf("same: %s ac: %lf phi0: %lf lpot: %lf\n", (c1 == c2 && c2 == c3 ? "yes" : "no" ), ac, phi0, 0.5 * l_angle_k * (ac - phi0) * (ac -phi0) ); 
#endif
		pot_angle += 0.5 * l_angle_k * (ac - phi0) * (ac -phi0);	
	
	}

	for( int d = 0; d < ndihe; d++ )
	{
//		if( d != dihe_loop_only ) continue;
/*
		if( dihe_loop_only < 0 ) 
		{
			if( d < - dihe_loop_only )
				continue;
		}
		else
		{
			if( d >= dihe_loop_only ) 
				continue;
		}
*/		int a1 = theDihe[4*d+0];
		int a2 = theDihe[4*d+1];
		int a3 = theDihe[4*d+2];
		int a4 = theDihe[4*d+3];
		
		int c1 = (a1-nattach) / aq_sites_per_nmer;
		int c2 = (a2-nattach) / aq_sites_per_nmer;
		int c3 = (a3-nattach) / aq_sites_per_nmer;
		int c4 = (a4-nattach) / aq_sites_per_nmer;
		double l_dihe_k = dihe_k[d];
		double th0 = dihe_th[d];

		
		double dr1A[3] = { r[3*a1+0] - r[3*a2+0],
				 r[3*a1+1] - r[3*a2+1],
				 r[3*a1+2] - r[3*a2+2] };
		double dr1B[3] = { r[3*a3+0] - r[3*a2+0],
				 r[3*a3+1] - r[3*a2+1],
				 r[3*a3+2] - r[3*a2+2] };
		
		double dr2A[3] = { r[3*a2+0] - r[3*a3+0],
				 r[3*a2+1] - r[3*a3+1],
				 r[3*a2+2] - r[3*a3+2] };
		double dr2B[3] = { r[3*a4+0] - r[3*a3+0],
				 r[3*a4+1] - r[3*a3+1],
				 r[3*a4+2] - r[3*a3+2] };
	
		double n1[3], n2[3];
	
		cross( dr1A,dr1B,n1);
		cross( dr2A,dr2B,n2);
		normalize(n1);
		normalize(n2);
		
		double dp = n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
	
		double x = dp;
		double m[3];
		double l = normalize(dr1B);
		
		cross( n1, dr1B, m );

		double y = dot( m, n2 );

		if( dp > 1 ) dp = 1 - 1e-13;
		double theta = acos(dp);
		double thetap = atan2( y, x );

		double dth = thetap-th0;
		while( dth > M_PI ) dth -= 2*M_PI;
		while( dth < -M_PI ) dth += 2*M_PI;

#ifdef DEBUG_DYNAMIN
		printf("dihe same: %s theta: %le th0: %le pot: %le\n",
			(c1==c2&&c2==c3&&c3==c4 ? "yes" : "no"), thetap, th0, l_dihe_k * (1-cos( thetap-th0)) );
#endif
		pot_dihe += l_dihe_k * (1-cos( thetap-th0));
//		pot_dihe += 0.5 * l_dihe_k * dth*dth;
	}




#ifdef DEBUG_DYNAMIN
	printf("POT NRM: %lf BOND: %lf ANGLE: %lf DIHE: %lf\n",
		pot_nrm, pot_bond, pot_angle, pot_dihe );
#endif
	double pot = pot_nrm + pot_bond + pot_angle + pot_dihe;

	if( ! (pot < 1 || pot > 0 )	 )
	{
		printf("dynamin nan.\n");
	}

	return pot;
}

// gets derivative of internal energy relative to position (surfacer_g) and the normal (surfacen_g).

double dynamin::grad( Simulation *theSimulation,  double *surfacer_g, double *surfacen_g )
{

	double *alphas = theSimulation->alpha;
	double r[3*nsites];
	double n[3*nsites];
	
	int nbonds = getNBonds();
	int theBonds[2*nbonds];
	double bond_r_p[nbonds];
	double bond_k_p[nbonds];
	putBonds(theBonds, bond_r_p, bond_k_p );
	
	int nangles = getNAngles();
	int theAngles[3*nangles];
	double angle_th[nangles];
	double angle_k[nangles];
	putAngles(theAngles, angle_th, angle_k );
	
	int ndihe = getNDihe();
	int theDihe[4*ndihe];
	double dihe_th[ndihe];
	double dihe_k[ndihe];
	putDihe(theDihe, dihe_th, dihe_k );
	
//	printf("dihe_loop_only: %d ndihe: %d\n", dihe_loop_only, ndihe );	

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
		// zero nrm for aqueous
		memset(n+nattach*3 , 0, sizeof(double) * (nsites-nattach)*3 );
	}

	double pot = 0;

	for( int m = 0; m < nmer_saved; m++ )
	{
		int s0 = m * attach_per_nmer;
		int sx = m * attach_per_nmer+1;
		int a0 = nattach + m * aq_sites_per_nmer;
		int ax = nattach + m * aq_sites_per_nmer + aq_sites_per_nmer-1;

		// normal potential, first site:
		double dr[3] = { r[a0*3+0] - r[s0*3+0], r[a0*3+1] - r[s0*3+1], r[a0*3+2] - r[s0*3+2] };
		double ln = normalize(dr);
		double dp = dr[0] * n[s0*3+0] + dr[1] * n[s0*3+1] + dr[2] * n[s0*3+2];

		double dp_target = 1;	

		if( is_inside )
			dp_target *= -1;
	
		double del = dp-dp_target;
	
		pot += 0.5 * k_nrm * del*del; 
		
		for( int c1 = 0; c1 < 3; c1++ )
		for( int c2 = 0; c2 < 3; c2++ )
		{
			if( c1 == c2 )
			{
				surfacer_g[3*a0+c1] += k_nrm*(del) * n[s0*3+c2] * (1/ln);
				surfacer_g[3*s0+c1] -= k_nrm*(del) * n[s0*3+c2] * (1/ln);
			}

			surfacer_g[3*a0+c1] += k_nrm*(del) * n[s0*3+c2]*(-(dr[c1]*dr[c2]/ln));		
			surfacer_g[3*s0+c1] -= k_nrm*(del) * n[s0*3+c2]*(-(dr[c1]*dr[c2]/ln));		
		}
		
		for( int c1 = 0; c1 < 3; c1++ )
			surfacen_g[s0*3+c1] += k_nrm*(del) * dr[c1];
		
		// normal potential, second site:

		dr[0] = r[3*ax+0] - r[3*(sx)+0];
		dr[1] = r[3*ax+1] - r[3*(sx)+1];
		dr[2] = r[3*ax+2] - r[3*(sx)+2];
		ln = normalize(dr);

		dp = dr[0] * n[sx*3+0] + dr[1] * n[sx*3+1] + dr[2] * n[sx*3+2];

		dp_target = 1;	

		if( is_inside )
			dp_target *= -1;
		del = dp - dp_target;
		pot += 0.5 * k_nrm * del*del;
		
		for( int c1 = 0; c1 < 3; c1++ )
		for( int c2 = 0; c2 < 3; c2++ )
		{
			if( c1 == c2 )
			{
				surfacer_g[3*ax+c1] += k_nrm*(del) * n[sx*3+c2] * (1/ln);
				surfacer_g[3*sx+c1] -= k_nrm*(del) * n[sx*3+c2] * (1/ln);
			}

			surfacer_g[3*ax+c1] += k_nrm*(del) * n[sx*3+c2]*(-(dr[c1]*dr[c2]/ln));		
			surfacer_g[3*sx+c1] -= k_nrm*(del) * n[sx*3+c2]*(-(dr[c1]*dr[c2]/ln));		
		}
	
	
		for( int c1 = 0; c1 < 3; c1++ )
			surfacen_g[sx*3+c1] += k_nrm*(del) * dr[c1];
	}

	for( int b = 0; b < nbonds; b++ )
	{
		int a1 = theBonds[2*b+0];
		int a2 = theBonds[2*b+1];

		double dr[3] = { r[3*a1+0] - r[3*a2+0], 
				 r[3*a1+1] - r[3*a2+1], 
				 r[3*a1+2] - r[3*a2+2] };
		
		double r0 = bond_r_p[b];
		double k  = bond_k_p[b];

		double ln = normalize(dr);

		pot += 0.5 * k * (ln - r0) * (ln-r0);	
		
		double f = k * (ln-r0); // dV / dx: 
		
		surfacer_g[3*a1+0] += dr[0] * f;
		surfacer_g[3*a1+1] += dr[1] * f;
		surfacer_g[3*a1+2] += dr[2] * f;

		surfacer_g[3*a2+0] -= dr[0] * f;
		surfacer_g[3*a2+1] -= dr[1] * f;
		surfacer_g[3*a2+2] -= dr[2] * f;
		
	}
	
	for( int a = 0; a < nangles; a++ )
	{
		int a1 = theAngles[3*a+0];
		int a2 = theAngles[3*a+1];
		int a3 = theAngles[3*a+2];

		double l_angle_k = angle_k[a];
		double phi0 = angle_th[a];

		double dr1[3] = { r[3*a1+0] - r[3*a2+0],
				 r[3*a1+1] - r[3*a2+1],
				 r[3*a1+2] - r[3*a2+2] };
		double dr2[3] = { r[3*a3+0] - r[3*a2+0],
				 r[3*a3+1] - r[3*a2+1],
				 r[3*a3+2] - r[3*a2+2] };
		double ln1 = normalize(dr1);
		double ln2 = normalize(dr2);

		double dp = dr1[0]*dr2[0]+dr1[1]*dr2[1]+dr1[2]*dr2[2];
		if( dp > 1 ) dp = 1 - 1e-13;
		double ac = acos(dp);

		pot += 0.5 * l_angle_k * (ac - phi0) * (ac -phi0);	
			
		double der_A = l_angle_k * (ac-phi0);
		double der_cos = -1.0/sqrt(1e-10+1.0-dp*dp); 
	
		// numerator derivative
	
		surfacer_g[3*a1+0] += der_A * der_cos * dr2[0] / ln1;
		surfacer_g[3*a1+1] += der_A * der_cos * dr2[1] / ln1;
		surfacer_g[3*a1+2] += der_A * der_cos * dr2[2] / ln1;
		
		surfacer_g[3*a2+0] += der_A * der_cos * (-dr2[0]) / ln1;
		surfacer_g[3*a2+1] += der_A * der_cos * (-dr2[1]) / ln1;
		surfacer_g[3*a2+2] += der_A * der_cos * (-dr2[2]) / ln1;
		
		surfacer_g[3*a2+0] += der_A * der_cos * (-dr1[0]) / ln2;
		surfacer_g[3*a2+1] += der_A * der_cos * (-dr1[1]) / ln2;
		surfacer_g[3*a2+2] += der_A * der_cos * (-dr1[2]) / ln2;
		
		surfacer_g[3*a3+0] += der_A * der_cos * (dr1[0]) / ln2;
		surfacer_g[3*a3+1] += der_A * der_cos * (dr1[1]) / ln2;
		surfacer_g[3*a3+2] += der_A * der_cos * (dr1[2]) / ln2;
	
		// denominator derivative.
		
		double ln1_2 = ln1*ln1;
		double ln2_2 = ln2*ln2;
	
		surfacer_g[3*a1+0] += -der_A * der_cos * dp * dr1[0] / ln1;
		surfacer_g[3*a1+1] += -der_A * der_cos * dp * dr1[1] / ln1;
		surfacer_g[3*a1+2] += -der_A * der_cos * dp * dr1[2] / ln1;
		
		surfacer_g[3*a2+0] += der_A * der_cos * dp * dr1[0] / ln1;
		surfacer_g[3*a2+1] += der_A * der_cos * dp * dr1[1] / ln1;
		surfacer_g[3*a2+2] += der_A * der_cos * dp * dr1[2] / ln1;
		
		surfacer_g[3*a2+0] += der_A * der_cos * dp * dr2[0] / ln2;
		surfacer_g[3*a2+1] += der_A * der_cos * dp * dr2[1] / ln2;
		surfacer_g[3*a2+2] += der_A * der_cos * dp * dr2[2] / ln2;
		
		surfacer_g[3*a3+0] += -der_A * der_cos * dp * dr2[0] / ln2;
		surfacer_g[3*a3+1] += -der_A * der_cos * dp * dr2[1] / ln2;
		surfacer_g[3*a3+2] += -der_A * der_cos * dp * dr2[2] / ln2;
	
	}

//#define ARC_COS_DIHE

#ifdef ARC_COS_DIHE	
	for( int d = 0; d < ndihe; d++ )
	{
		int a1 = theDihe[4*d+0];
		int a2 = theDihe[4*d+1];
		int a3 = theDihe[4*d+2];
		int a4 = theDihe[4*d+3];
		double l_dihe_k = dihe_k[d];
		double th0 = dihe_th[d];

		
		double dr1A[3] = { r[3*a1+0] - r[3*a2+0],
				 r[3*a1+1] - r[3*a2+1],
				 r[3*a1+2] - r[3*a2+2] };
		double dr1B[3] = { r[3*a3+0] - r[3*a2+0],
				 r[3*a3+1] - r[3*a2+1],
				 r[3*a3+2] - r[3*a2+2] };
		
		double dr2A[3] = { r[3*a2+0] - r[3*a3+0],
				 r[3*a2+1] - r[3*a3+1],
				 r[3*a2+2] - r[3*a3+2] };
		double dr2B[3] = { r[3*a4+0] - r[3*a3+0],
				 r[3*a4+1] - r[3*a3+1],
				 r[3*a4+2] - r[3*a3+2] };
	
		double n1[3], n2[3];
	
		cross( dr1A,dr1B,n1);
		cross( dr2A,dr2B,n2);
		normalize(n1);
		normalize(n2);
		double dp = n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
		if( dp > 1 ) dp = 1 - 1e-13;
		double thetao = acos(dp);
		double x = dp;
		double m[3];
		double l = normalize(dr1B);
		
		cross( n1, dr1B, m );

		double y = dot( m, n2 );

		if( dp > 1 ) dp = 1 - 1e-13;
		double theta = atan2( y, x );
		double the_sign = 1;

		if( theta < 0 )
			the_sign = -1;

		pot += 0.5 * l_dihe_k * (theta - th0 ) * (theta - th0 );

		// derivative of nrm wrt the three vectors.
		// slow: cartesian component of the normal
		// mid: which vector (r1,r2,r3)
		// fast: cartesian component of that vector
		double dnrm1[27];
		double dnrm2[27];

		// order here is mixed up, sorry.
		normal_cp_der( r+3*a2, r+3*a1, r+3*a3, dnrm1 );
		normal_cp_der( r+3*a3, r+3*a2, r+3*a4, dnrm2 );

		// everything prior to d(dp)/dr1
		double base_der = l_dihe_k * (theta - th0) * (-the_sign / sqrt(1e-10 + 1-dp*dp));
	
		// everything else is d(dp)/dx_i, etc.
		for( int pc = 0; pc < 3; pc++ )
		for( int nc = 0; nc < 3; nc++ )
		{
			surfacer_g[3*a1+pc] += base_der * n2[nc] * dnrm1[nc*9+1*3+pc];
			surfacer_g[3*a2+pc] += base_der * n2[nc] * dnrm1[nc*9+0*3+pc];
			surfacer_g[3*a3+pc] += base_der * n2[nc] * dnrm1[nc*9+2*3+pc];
			
			surfacer_g[3*a2+pc] += base_der * n1[nc] * dnrm2[nc*9+1*3+pc];
			surfacer_g[3*a3+pc] += base_der * n1[nc] * dnrm2[nc*9+0*3+pc];
			surfacer_g[3*a4+pc] += base_der * n1[nc] * dnrm2[nc*9+2*3+pc];
		}
	}
#else
	for( int d = 0; d < ndihe; d++ )
	{
//		if( d != dihe_loop_only ) continue;
/*	
		if( dihe_loop_only < 0 ) 
		{
			if( d < - dihe_loop_only )
				continue;
		}
		else
		{
			if( d >= dihe_loop_only ) 
				continue;
		}
*/
		int a1 = theDihe[4*d+0];
		int a2 = theDihe[4*d+1];
		int a3 = theDihe[4*d+2];
		int a4 = theDihe[4*d+3];
		double l_dihe_k = dihe_k[d];
		double th0 = dihe_th[d];

		
		double dr1A[3] = { r[3*a1+0] - r[3*a2+0],
				 r[3*a1+1] - r[3*a2+1],
				 r[3*a1+2] - r[3*a2+2] };
		double dr1B[3] = { r[3*a3+0] - r[3*a2+0],
				 r[3*a3+1] - r[3*a2+1],
				 r[3*a3+2] - r[3*a2+2] };
		
		double dr2A[3] = { r[3*a2+0] - r[3*a3+0],
				 r[3*a2+1] - r[3*a3+1],
				 r[3*a2+2] - r[3*a3+2] };
		double dr2B[3] = { r[3*a4+0] - r[3*a3+0],
				 r[3*a4+1] - r[3*a3+1],
				 r[3*a4+2] - r[3*a3+2] };
	
		double n1[3], n2[3];
	
		cross( dr1A,dr1B,n1);
		cross( dr2A,dr2B,n2);
		normalize(n1);
		normalize(n2);
		double dp = n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
		if( dp > 1 ) dp = 1 - 1e-13;
		double thetao = acos(dp);
		double x = dp;
		double m[3];
		double l = normalize(dr1B);
		
		cross( n1, dr1B, m );

		double y = dot( m, n2 );

		if( dp > 1 ) dp = 1 - 1e-13;
		double theta = atan2( y, x );
		double the_sign = 1;

		if( theta < 0 )
			the_sign = -1;

		double dth = theta-th0;
//		while( dth > M_PI ) dth -= 2*M_PI;
//		while( dth < -M_PI ) dth += 2*M_PI;
		//pot += 0.5 * l_dihe_k * dth*dth;
		
		double lpot = l_dihe_k * (1-cos(dth));			
		pot += l_dihe_k * (1-cos(dth));

#ifdef DEBUG_DYNAMIN
		printf("theta: %le th0: %le lpot: %le\n", theta, th0, lpot );
#endif
		double der_d_th[12];
		brute_force_dihedral( r+a1*3, r+a2*3, r+a3*3, r+a4*3, der_d_th ); 

		for( int c = 0; c < 3; c++ )
		{
#ifdef QUAD
			surfacer_g[3*a1+c] -= l_dihe_k * (theta-th0) * der_d_th[0+c];
			surfacer_g[3*a2+c] -= l_dihe_k * (theta-th0) * der_d_th[3+c];
			surfacer_g[3*a3+c] -= l_dihe_k * (theta-th0) * der_d_th[6+c];
			surfacer_g[3*a4+c] -= l_dihe_k * (theta-th0) * der_d_th[9+c];
#else
			surfacer_g[3*a1+c] -= l_dihe_k * sin(theta-th0) * der_d_th[0+c];
			surfacer_g[3*a2+c] -= l_dihe_k * sin(theta-th0) * der_d_th[3+c];
			surfacer_g[3*a3+c] -= l_dihe_k * sin(theta-th0) * der_d_th[6+c];
			surfacer_g[3*a4+c] -= l_dihe_k * sin(theta-th0) * der_d_th[9+c];
#endif
		}
	}
#endif


	double check = V( theSimulation );

	printf("GV: %le V: %le\n", pot, check );
	
	if( ! (pot < 1 || pot > 0 )	 )
	{
		printf("dynamin nan.\n");
	}

	return pot;
}

void dynamin::move_inside( void )
{
	pcomplex::move_inside();
}

void dynamin::move_outside( void )
{
	pcomplex::move_outside();
}

void dynamin::getDirName( char name[256] )
{
	sprintf(name,"dynamin");
}
void sdynamin::getDirName( char name[256] )
{
	sprintf(name,"sdynamin");
}


void dynamin::writeStructure( Simulation *theSimulation, surface_mask *upperSurfaceMask, surface_mask *lowerSurfaceMask, struct atom_rec **at_out, int *nat_out, char ***sequence_array, int *nseq_out, int **seq_at_array, char ***patches, struct ion_add **ions, int *nions, 
	aa_build_data *buildData, int *build_type )
{
	*build_type = BUILD_SEQUENCE;

	int nsegments = nmer_saved;

	char dirName[256];

	getDirName( dirName );

	struct atom_rec *at_fetch;
	int nat;
	int pool_code = pdbFetch( &at_fetch, &nat, dirName, "dynamin_martini", addToPool );
	
	const char *segment_names[] = { "PROA", "PROB" };

	int a_markers[nmer_saved*2];
	int tadd = 0;
	double com[nmer_saved*2][3];
	int ncom[nmer_saved];
	
	int a_start = *nat_out;
	
	*patches = (char **)malloc( sizeof(char *) * 2*nmer_saved );
	
	for( int m = 0; m < nmer_saved; m++ )
	for( int seg = 0; seg < 2; seg++,tadd++ )
	{
		patchFetch( (*patches)+tadd, dirName, segment_names[seg] );

		a_markers[tadd] = *nat_out;

		int ntries = 0;
		double random_noise = 0.01;

		int clash = 1;
		// addPeri... returns one if there's a clash.. keep trying, injecting more noise

		// dimer interface 178
		// elbow_1 	   388
		// elbow_2	   712	

		int res[3] = { 388, 712, 178 };
		int cmplx[3];

		if( seg == 0 )
		{
			res[0] = 388;	
			res[1] = 712;	
			res[2] = 178;	
			cmplx[0] = nattach + aq_sites_per_nmer * m + 0;	
			cmplx[1] = nattach + aq_sites_per_nmer * m + 1;	
			cmplx[2] = nattach + aq_sites_per_nmer * m + 2;	
		}
		else if( seg == 1 )
		{
			res[0] = 388;	
			res[1] = 712;	
			res[2] = 178;	

			cmplx[0] = nattach + aq_sites_per_nmer * m + 4;	
			cmplx[1] = nattach + aq_sites_per_nmer * m + 3;	
			cmplx[2] = nattach + aq_sites_per_nmer * m + 2;	
		}

		while( (clash = addAlignedProteinHelper( theSimulation, upperSurfaceMask, lowerSurfaceMask,
			at_out, nat_out,
			ions, nions,	
			pool_code,
			segment_names[seg],
			3,
			res,
			cmplx,	
			buildData, random_noise )) && ntries < 500)
		{
			// re-orient.
		
			if( ntries % 5 == 0 )
				random_noise *= 1.05;	
			ntries++;
		}

		if( clash )
		{
			printf("ADD FAILED. clash.\n");
			exit(1);
		}

		com[tadd][0] = 0;
		com[tadd][1] = 0;
		com[tadd][2] = 0;
		ncom[tadd]=0;

		for( int a = a_markers[tadd]; a < *nat_out; a++ )
		{
			com[tadd][0] += (*at_out)[a].x;
			com[tadd][1] += (*at_out)[a].y;
			com[tadd][2] += (*at_out)[a].z;
			ncom[tadd] += 1;
		}
	
		com[tadd][0] /= ncom[tadd];
		com[tadd][1] /= ncom[tadd];
		com[tadd][2] /= ncom[tadd];
	}

#ifdef UPDATE_THIS
	// fix centers-of-mass ?	

	double dr[3] = { com[1][0] - com[0][0], com[1][1] - com[0][1], com[1][2] - com[0][2] };
	double put[3];

	MinImage3D( dr, theSimulation->PBC_vec, put ); 

	for( int a = a_markers[1]; a < *nat_out; a++ )
	{
		(*at_out)[a].x += put[0] * theSimulation->PBC_vec[0][0] + put[1] * theSimulation->PBC_vec[1][0] + put[2] * theSimulation->PBC_vec[2][0];
		(*at_out)[a].y += put[0] * theSimulation->PBC_vec[0][1] + put[1] * theSimulation->PBC_vec[1][1] + put[2] * theSimulation->PBC_vec[2][1];
		(*at_out)[a].z += put[0] * theSimulation->PBC_vec[0][2] + put[1] * theSimulation->PBC_vec[1][2] + put[2] * theSimulation->PBC_vec[2][2];
	}

#endif

	struct atom_rec *at = *at_out;
	// make sure the domains are in the same PBC



	*sequence_array = (char **)malloc( sizeof(char*) * tadd );
	*seq_at_array = (int *)malloc( sizeof(int) * tadd );

	*nseq_out = tadd;


	for( int t = 0; t < tadd; t++ )
	{
		int seq_space = 10;
		int nseq = 0;

	
		char *seq = (char *)malloc( sizeof(char) * seq_space );
	
		int start = a_markers[t];
		
		seq[nseq] = threeToOne( at[start].resname );
		nseq++;
		int stop = *nat_out; 
		if( t < tadd-1 )
			stop = a_markers[t+1];
	
		(*seq_at_array)[t] = start;

		int pres = at[start].res;
		for( int a = start; a < stop; a++ )
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
	
				for( int ts = nseq; ts < nseq-1+doff; ts++ )
					seq[ts] = 'X';		
	
				nseq = nseq-1+doff;
	
				seq[nseq] = threeToOne( at[a].resname ); 
				nseq++;
			}
	
			pres = at[a].res;
		}	
	
		seq[nseq] = '\0';
	
		const char *missing_residues="GSG";
		int nmissing = 0;
		int ntodo = strlen(missing_residues);
	
		for( int ts = 0; ts < nseq; ts++ )
		{
			if( seq[ts] == 'X' )
			{
				if( nmissing >= ntodo )
				{
					printf("ERROR. found a missing residue but didn't have enough to insert.\n");
					exit(1);
				}
	
				seq[ts] = missing_residues[nmissing];
				nmissing++;
			}
		}
		(*sequence_array)[t] = seq;
		
	}
}

/* this function 
 * "gets" the coarse-grained coordinates of the observed protein from the molecular structure -- and its fit to the membrane surface.
 *
 * */


void dynamin::get( 
		Simulation *theSimulation, // this is the surface/PBC
		struct atom_rec *at,
		int syt7_start, int syt7_stop, int nat_tot )
{
#if 0
	int nsegments = 2;
	
	surface *theSurface;
	double *rsurf;
	theSimulation->fetch(sid[0],&theSurface,&rsurf);
	surface_record *sRec = theSimulation->fetch(sid[0]);

	int flip_sign = sRec->gather_flip;

	// how do we define the sites? we can tweak the definition here.
	
	// BRIDGE 263
	
	// Calcium binding ASPs: 166, 172, 225, 227, 233
 	//                         ARGs: 217
 	//                         LYS: 213

	// Calcium binding ASPs: 297, 303, 359, 365 
 	//		   ARGs: 347, 390, 392                         

	int npts = 6;
	// p is where we'll put the averages of the residues below.
	double p[npts*3]; // 
	memset( p, 0, sizeof(double) * 3 * npts );
	int nav[npts];
	memset( nav, 0, sizeof(int) * npts );
	// these are the residues we'll average over.
	int res_av[6][9] = {
		{166, 172, 225, 227, 233, -1},
		{217, 213, -1 },
		
		{ 297, 303, 359, 365, -1 },
		{ 347, 390, 392, -1 },
		
		{138,159,197,255,174,221,238,190},
		{399,270,289,332,320,307,353,369}
		//{139, 176, -1},
		//{289, 354, -1}
	};
	
	// for near pointing
	double **M;
	int mlow,mhigh;

	getM( &M, &mlow, &mhigh );

	// just for PBC wrapping:
	double alphas[3] = {1,1,1};
	int nentries[6] = { 5, 2, 2, 4, 8, 8 };


	for( int a = syt7_start; a < syt7_stop; a++ )
	{
		for( int px = 0; px < npts-2; px++ )
		{
			for( int rx = 0; rx < nentries[px]; rx++ )
			{
				if( at[a].res == res_av[px][rx] )
				{
					if( nav[px] == 0 )
					{
						p[3*px+0] = at[a].x;
						p[3*px+1] = at[a].y;
						p[3*px+2] = at[a].z;
					} 
					else
					{
						double dr[3] = { at[a].x - p[3*px+0]/nav[px],
								 at[a].y - p[3*px+1]/nav[px],
								 at[a].z - p[3*px+2]/nav[px] };

						theSimulation->wrapPBC(dr,alphas);

						p[3*px+0] += dr[0] + p[3*px+0]/nav[px];
						p[3*px+1] += dr[1] + p[3*px+1]/nav[px];
						p[3*px+2] += dr[2] + p[3*px+2]/nav[px];
					}

					nav[px] += 1;
				}	
			}
		}
	}	
	
	int min_res[2] = { 1, 265};
	int max_res[2] = { 266, 404 };
	int tp = 0;
	for( int px = 4; px <= 5; px++ )
	{
		for( int a = syt7_start; a < syt7_stop; a++ )
		{
			if( at[a].res >= min_res[tp] && at[a].res <= max_res[tp] )
			{
				if( nav[px] == 0 )
				{
					p[3*px+0] = at[a].x;
					p[3*px+1] = at[a].y;
					p[3*px+2] = at[a].z;
				} 
				else
				{
					double dr[3] = { at[a].x - p[3*px+0]/nav[px],
							 at[a].y - p[3*px+1]/nav[px],
							 at[a].z - p[3*px+2]/nav[px] };

					theSimulation->wrapPBC(dr,alphas);

					p[3*px+0] += dr[0] + p[3*px+0]/nav[px];
					p[3*px+1] += dr[1] + p[3*px+1]/nav[px];
					p[3*px+2] += dr[2] + p[3*px+2]/nav[px];
				}

				nav[px] += 1;
			}
		}
	}

	for( int px = 0; px < npts; px++ )
	{
		p[3*px+0] /= nav[px];
		p[3*px+1] /= nav[px];
		p[3*px+2] /= nav[px];
	}

	for( int px = 1; px < npts; px++ )
	{
		double dr[3] = { p[3*px+0] - p[0],
				 p[3*px+1] - p[1],
				 p[3*px+2] - p[2] };

		theSimulation->wrapPBC(dr,alphas);

		p[3*px+0] = p[0] + dr[0];
		p[3*px+1] = p[1] + dr[1];
		p[3*px+2] = p[2] + dr[2];
	}
	
	
	double nrm1[3],nrm2[3],nrm3[3],nrm4[4];
	double rpt1[3],rpt2[3],rpt3[3],rpt4[4];

	// find near points for the attachment sites.
	double dist;
	int f;
	double col_u, col_v;
	theSurface->nearPointOnBoxedSurface( p, &f, &col_u, &col_v, M, mlow, mhigh, &dist );					
	double rloop[3],nloop[3];
	theSurface->evaluateRNRM( f, col_u, col_v, rloop, nloop, rsurf );
	theSurface->evaluateRNRM( f, col_u, col_v, rpt1, nrm1, rsurf );

	if( flip_sign ) { nrm1[0] *= -1; nrm1[1] *= -1; nrm1[2] *= -1; }
	
	double c_vec_1_A[3];
	double c_vec_2_A[3];
	double c_val1_A, c_val2_A;
	double k_A;
	theSurface->c( f, col_u, col_v, rsurf, &k_A, c_vec_1_A, c_vec_2_A, &c_val1_A, &c_val2_A ); 
	if( flip_sign) { c_val1_A *= -1; c_val2_A *= -1; }

	for( int tp =0; tp < npts; tp++ )	
	{
		double dr[3];
		dr[0] = p[3*tp+0]-rloop[0];
		dr[1] = p[3*tp+1]-rloop[1];
		dr[2] = p[3*tp+2]-rloop[2];
		theSimulation->wrapPBC(dr,alphas);
		p[3*tp+0] = rloop[0] + dr[0];
		p[3*tp+1] = rloop[1] + dr[1];
		p[3*tp+2] = rloop[2] + dr[2];
	}

	puv[0] = col_u;
	puv[1] = col_v;
	fs[0] = f;
	theSurface->nearPointOnBoxedSurface( p+3, &f, &col_u, &col_v, M, mlow, mhigh, &dist );					

	theSurface->evaluateRNRM( f, col_u, col_v, rpt2, nrm2, rsurf );
	if( flip_sign ) { nrm2[0] *= -1; nrm2[1] *= -1; nrm2[2] *= -1; }
	

	puv[2] = col_u;
	puv[3] = col_v;
	fs[1] = f;
	theSurface->nearPointOnBoxedSurface( p+6, &f, &col_u, &col_v, M, mlow, mhigh, &dist );					
	theSurface->evaluateRNRM( f, col_u, col_v, rpt3, nrm3, rsurf );
	if( flip_sign ) { nrm3[0] *= -1; nrm3[1] *= -1; nrm3[2] *= -1; }
	double c_vec_1_B[3];
	double c_vec_2_B[3];
	double c_val1_B, c_val2_B;
	double k_B;
	theSurface->c( f, col_u, col_v, rsurf, &k_B, c_vec_1_B, c_vec_2_B, &c_val1_B, &c_val2_B ); 
	if( flip_sign) { c_val1_B *= -1; c_val2_B *= -1; }
	puv[4] = col_u;
	puv[5] = col_v;
	fs[2] = f;
	theSurface->nearPointOnBoxedSurface( p+9, &f, &col_u, &col_v, M, mlow, mhigh, &dist );					
	theSurface->evaluateRNRM( f, col_u, col_v, rpt4, nrm4, rsurf );
	if( flip_sign ) { nrm4[0] *= -1; nrm4[1] *= -1; nrm4[2] *= -1; }
	puv[6] = col_u;
	puv[7] = col_v;
	fs[3] = f;
	
	rall[4*3+0] = p[12];	
	rall[4*3+1] = p[13];	
	rall[4*3+2] = p[14];	
	
	rall[5*3+0] = p[15];	
	rall[5*3+1] = p[16];		
	rall[5*3+2] = p[17];		

	double rvecA[3] = { rall[4*3+0] - p[0],
			    rall[4*3+1] - p[1],
			    rall[4*3+2] - p[2] };
	normalize(rvecA);
	normalize(rvecB);

	double dpA = rvecA[0]*nrm1[0] + rvecA[1] * nrm1[1] + rvecA[2] * nrm1[2];
	double dpB = rvecB[0]*nrm3[0] + rvecB[1] * nrm3[1] + rvecB[2] * nrm3[2];
	

	// nrm points away from syt
	printf("SYT %d C2A_DP: %lf C2B_DP: %lf c_A: %le %le c_B: %le %le\n", my_id, dpA, dpB, c_val1_A, c_val2_A, c_val1_B, c_val2_B );
	
	int do_track = 1;

	if( do_track )
	{
		int track_com[] = { 392, 390, 382, 385, 347 };
	
		double com_to_track[3] = {0,0,0};
		int ntr_com=0;
	
		int ntr = sizeof(track_com)/sizeof(int);
	
		for( int a = syt7_start; a < syt7_stop; a++ )
		{
			for( int tx = 0; tx < ntr; tx++ )
			{
				if( track_com[tx] != at[a].res )
				 	continue;
	
				if( ntr_com == 0 )
				{
					double dr[3] = { at[a].x - p[5*3+0],
							 at[a].y - p[5*3+1],
							 at[a].z - p[5*3+2]};
	
					theSimulation->wrapPBC(dr,alphas);
	
					com_to_track[0] += p[5*3+0] + dr[0];
					com_to_track[1] += p[5*3+1] + dr[1];
					com_to_track[2] += p[5*3+2] + dr[2];
				} 
				else
				{
					double dr[3] = { at[a].x - com_to_track[0]/ntr_com,
							 at[a].y - com_to_track[1]/ntr_com,
							 at[a].z - com_to_track[2]/ntr_com};
	
					theSimulation->wrapPBC(dr,alphas);
	
					com_to_track[0] += dr[0] + com_to_track[0]/ntr_com;
					com_to_track[1] += dr[1] + com_to_track[1]/ntr_com;
					com_to_track[2] += dr[2] + com_to_track[2]/ntr_com;
				}
	
				ntr_com += 1;
			}
		}
	
		com_to_track[0] /= ntr_com;
		com_to_track[1] /= ntr_com;
		com_to_track[2] /= ntr_com;
		
		double dist;
		int f;
		double col_u, col_v;
		theSurface->nearPointOnBoxedSurface( com_to_track, &f, &col_u, &col_v, M, mlow, mhigh, &dist );					
		double ptr[3], nrmr[3];
		theSurface->evaluateRNRM( f, col_u, col_v, ptr, nrmr, rsurf );
		double vec1[3] = { com_to_track[0] - ptr[0], com_to_track[1] - ptr[1], com_to_track[2] - ptr[2] }; 
		double vec2[3] = { com_to_track[0] - p[5*3+0], com_to_track[1] - p[5*3+1], com_to_track[2] - p[5*3+2] }; 
		theSimulation->wrapPBC(vec1,alphas);
		normalize(vec1);
		normalize(vec2);
		double dp = vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
		printf("DEL %d dist %le dp: %le\n", my_id, dist, dp );
	}
#endif
}


char dynamin::getSiteCode( int p )
{
	if( p < nattach ) return 'N';
	return 'O';

}

/*int inter_bonds_per_nmer = aq_sites_per_nmer-1;
int attach_per_nmer = 2;
*/
int sdynamin::getNBonds( void )
{
	
	int nsites_per_nmer = attach_per_nmer + aq_sites_per_nmer;
	
	int bonds_per_nmer = (aq_sites_per_nmer + attach_per_nmer - 1) + inter_bonds_per_nmer;

	int nbonds = bonds_per_nmer * nmer_saved - inter_bonds_per_nmer;

	return nbonds;
}

/*
int nangles_intra = attach_per_nmer + aq_sites_per_nmer - 2
#ifdef SUPP_DIHE
+ 2
#endif
;
int nangles_inter = inter_bonds_per_nmer; // one per bond 
*/
int sdynamin::getNAngles( void )
{
	int nangles = nangles_intra *nmer_saved + nangles_inter * (nmer_saved-1);


	return nangles;
}


void sdynamin::putAngles( int *ang_list, double *angle_theta, double *angle_k )
{
	double use_angle_k = angle_force;


	int nangles = nangles_intra *nmer_saved + nangles_inter * (nmer_saved-1);

	// intra angles.
	
	int cnt = 0;
	for( int m = 0; m < nmer_saved; m++ )
	{
		int a0 = nattach + m * aq_sites_per_nmer;
		int s0 = attach_per_nmer * m;
// ang 1
		ang_list[3*cnt+0] = s0; 
		ang_list[3*cnt+1] = a0; 
		ang_list[3*cnt+2] = a0+1; 

		angle_theta[cnt] = 122.378814;
		angle_k[cnt] = use_angle_k;
	
		cnt++;
// ang 2	
		ang_list[3*cnt+0] = a0; 
		ang_list[3*cnt+1] = a0+1; 
		ang_list[3*cnt+2] = a0+2; 

		angle_theta[cnt] = 106.707580;
		angle_k[cnt] = use_angle_k;

		cnt++;
// ang 3	
		ang_list[3*cnt+0] = a0+1; 
		ang_list[3*cnt+1] = a0+2; 
		ang_list[3*cnt+2] = a0+3; 

		angle_theta[cnt] = 157.978543;
		angle_k[cnt] = use_angle_k;

		cnt++;
// ang 4	
		ang_list[3*cnt+0] = a0+2; 
		ang_list[3*cnt+1] = a0+3; 
		ang_list[3*cnt+2] = a0+4; 

		angle_theta[cnt] = 102.544231;
		angle_k[cnt] = use_angle_k;
		
		cnt++;
// ang 5
		ang_list[3*cnt+0] = a0+3; 
		ang_list[3*cnt+1] = a0+4; 
		ang_list[3*cnt+2] = s0+1; 

		angle_theta[cnt] = 108.416539;
		angle_k[cnt] = use_angle_k;
		
		cnt++;
#ifdef SUPP_DIHE
// ang S1
		ang_list[3*cnt+0] = a0; 
		ang_list[3*cnt+1] = a0+1; 
		ang_list[3*cnt+2] = a0+3; 

		angle_theta[cnt] = 100.0;
		angle_k[cnt] = use_angle_k;
		
		cnt++;

// ang S2
		ang_list[3*cnt+0] = a0+1; 
		ang_list[3*cnt+1] = a0+3; 
		ang_list[3*cnt+2] = a0+4; 

		angle_theta[cnt] = 100.0;
		angle_k[cnt] = use_angle_k;
		
		cnt++;


#endif


	}
	
	// inter angles.
	
	for( int m = 0; m < nmer_saved-1; m++ )
	{
		int a0 = nattach + m * aq_sites_per_nmer;
		int s0 = attach_per_nmer * m;
		
		int a1 = nattach + (m+1) * aq_sites_per_nmer;
		int s1 = attach_per_nmer * (m+1);
// ang 1
		ang_list[3*cnt+0] = a0+1; 
		ang_list[3*cnt+1] = a1; 
		ang_list[3*cnt+2] = a1+1; 

		angle_theta[cnt] = 66.595828;
		angle_k[cnt] = use_angle_k;
	
		cnt++;
// ang 2	
		ang_list[3*cnt+0] = a0+2; 
		ang_list[3*cnt+1] = a1+1; 
		ang_list[3*cnt+2] = a1+2; 

		angle_theta[cnt] = 68.575376;
		angle_k[cnt] = use_angle_k;

		cnt++;
// ang 3	
		ang_list[3*cnt+0] = a0+3; 
		ang_list[3*cnt+1] = a1+2; 
		ang_list[3*cnt+2] = a1+3; 

		angle_theta[cnt] = 54.117237; 
		angle_k[cnt] = use_angle_k;

		cnt++;
// ang 4	
		ang_list[3*cnt+0] = a0+4; 
		ang_list[3*cnt+1] = a1+3; 
		ang_list[3*cnt+2] = a1+4; 

		angle_theta[cnt] = 60.885045;
		angle_k[cnt] = use_angle_k;
		
		cnt++;
	}

	// convert to rads
	
	for( int a = 0; a < cnt; a++ )
		angle_theta[a] *= (M_PI/180.0);
}

/*

int ndihe_inter = 4;

int ndihe_intra = 
2 + 

#ifdef MID_DIHE
2 +
#endif

#ifdef SUPP_DIHE
1
#else
0
#endif
;
*/

int sdynamin::getNDihe( void )
{
	int ndihe = ndihe_intra * nmer_saved + ndihe_inter * (nmer_saved-1);

	return ndihe;
}

void sdynamin::putDihe( int *dihe_list, double *dihe_theta, double *dihe_k )
{
	double use_dihe_k = dihedral_force;

	int ndihe = getNDihe();

	// intra dihe.
	
	int cnt = 0;
	for( int m = 0; m < nmer_saved; m++ )
	{
		int a0 = nattach + m * aq_sites_per_nmer;
		int s0 = attach_per_nmer * m;
// dihe 1
		dihe_list[4*cnt+0] = s0; 
		dihe_list[4*cnt+1] = a0; 
		dihe_list[4*cnt+2] = a0+1; 
		dihe_list[4*cnt+3] = a0+2; 

		dihe_theta[cnt] = 33.294000*(M_PI/180.0);
		dihe_k[cnt] = use_dihe_k;

	
		cnt++;

#ifdef MID_DIHE
// dihe 2
		dihe_list[4*cnt+0] = a0; 
		dihe_list[4*cnt+1] = a0+1; 
		dihe_list[4*cnt+2] = a0+2; 
		dihe_list[4*cnt+3] = a0+3; 

		dihe_theta[cnt] = 5.590809*(M_PI/180.0);
		dihe_k[cnt] = use_dihe_k;
	
		cnt++;
// dihe 3
		dihe_list[4*cnt+0] = a0+1; 
		dihe_list[4*cnt+1] = a0+2; 
		dihe_list[4*cnt+2] = a0+3; 
		dihe_list[4*cnt+3] = a0+4; 

		dihe_theta[cnt] = 22.913731*(M_PI/180.0);
		dihe_k[cnt] = use_dihe_k;
	
		cnt++;
#endif

#ifdef SUPP_DIHE
		dihe_theta[cnt] = 22.913731 * (M_PI/180.0);
		dihe_list[4*cnt+0] = a0; 
		dihe_list[4*cnt+1] = a0+1; 

		dihe_list[4*cnt+2] = a0+3; 
		dihe_list[4*cnt+3] = a0+4; 
		dihe_k[cnt] = use_dihe_k;

		cnt++;
#endif
// dihe 4
		dihe_list[4*cnt+0] = a0+2; 
		dihe_list[4*cnt+1] = a0+3; 
		dihe_list[4*cnt+2] = a0+4; 
		dihe_list[4*cnt+3] = s0+1; 

		dihe_theta[cnt] = 7.612105*(M_PI/180.0);
		dihe_k[cnt] = use_dihe_k;
	
		cnt++;

	}
	
	// inter dihe.
	
	for( int m = 0; m < nmer_saved-1; m++ )
	{
		int a0 = nattach + m * aq_sites_per_nmer;
		int s0 = attach_per_nmer * m;
		
		int a1 = nattach + (m+1) * aq_sites_per_nmer;
		int s1 = attach_per_nmer * (m+1);
// dihe 1		
		dihe_list[4*cnt+0] = a0; 
		dihe_list[4*cnt+1] = a0+1; 
		dihe_list[4*cnt+2] = a1; 
		dihe_list[4*cnt+3] = a1+1; 

		dihe_theta[cnt] = 168.493706*(M_PI/180.0);
		dihe_k[cnt] = use_dihe_k;
	
		cnt++;

// dihe 2
		dihe_list[4*cnt+0] = a0+1;
		dihe_list[4*cnt+1] = a0+2; 
		dihe_list[4*cnt+2] = a1+1; 
		dihe_list[4*cnt+3] = a1+2; 

		dihe_theta[cnt] = 171.819312*(M_PI/180.0);
		dihe_k[cnt] = use_dihe_k;
	
		cnt++;

// dihe 3
		dihe_list[4*cnt+0] = a0+2;
		dihe_list[4*cnt+1] = a0+3; 
		dihe_list[4*cnt+2] = a1+2; 
		dihe_list[4*cnt+3] = a1+3; 

		dihe_theta[cnt] = 148.662895*(M_PI/180.0);
		dihe_k[cnt] = use_dihe_k;
	
		cnt++;

// dihe 3
		dihe_list[4*cnt+0] = a0+3;
		dihe_list[4*cnt+1] = a0+4; 
		dihe_list[4*cnt+2] = a1+3; 
		dihe_list[4*cnt+3] = a1+4; 

		dihe_theta[cnt] = 151.237723*(M_PI/180.0);
		dihe_k[cnt] = use_dihe_k;
	
		cnt++;
	}

	// convert to rads
	
//	for( int a = 0; a < cnt; a++ )
//		dihe_theta[a] *= (M_PI/180.0);
}

void sdynamin::putBonds( int *bond_list, double *bond_r, double *bond_k_in)
{
	double bond_k = bond_force;
	int attach_per_nmer = 2;
	int nsites_per_nmer = attach_per_nmer + aq_sites_per_nmer;
	
	int bonds_per_nmer = (aq_sites_per_nmer + attach_per_nmer - 1) + aq_sites_per_nmer-1;

/*
double elbow_1_to_membrane = 45.0;
double elbow_angle	 = 120.0;
double elbow_to_elbow    = 70.0;
double plane_angle       = 21.094992;
double el_el_G_angle     = 70.596031;
double elbow_to_G	 = 84.113721;
double el_G_el_angle  = (161.273414);
double final_dihe	 = 11.862716;
*/
	for( int m = 0; m < nmer_saved; m++ )
	{
		// internal bonds.
		int a0 = nattach + m * aq_sites_per_nmer;
		int s0 = attach_per_nmer * m;
		int nb = 0;

		bond_list[(m*bonds_per_nmer)*2+nb*2+0] = s0;
		bond_list[(m*bonds_per_nmer)*2+nb*2+1] = a0;
		if( bond_r ) bond_r[m*bonds_per_nmer+nb] = elbow_1_to_membrane;
		if( bond_k_in ) bond_k_in[m*bonds_per_nmer+nb] = bond_k;
		nb++;
		
		bond_list[(m*bonds_per_nmer)*2+nb*2+0] = a0+aq_sites_per_nmer-1;
		bond_list[(m*bonds_per_nmer)*2+nb*2+1] = s0+1;
		if( bond_r ) bond_r[m*bonds_per_nmer+nb] = elbow_1_to_membrane;
		if( bond_k_in ) bond_k_in[m*bonds_per_nmer+nb] = bond_k;
		nb++;

		for( int aq = 0; aq < aq_sites_per_nmer-1; aq++ )
		{
			bond_list[(m*bonds_per_nmer)*2+nb*2+0] = a0+aq;
			bond_list[(m*bonds_per_nmer)*2+nb*2+1] = a0+aq+1;

			
			if( bond_r )  {
				switch( aq )
				{
					case 0:
					case 3: 
						bond_r[m*bonds_per_nmer+nb] = elbow_to_elbow;
						break;
					case 1:
					case 2: 
						bond_r[m*bonds_per_nmer+nb] = elbow_to_G;
						break;
					
				}
			}
			if( bond_k_in ) bond_k_in[m*bonds_per_nmer+nb] = bond_k;

			nb++;
		}

		// inter bonds.
		
		if( m < nmer_saved-1 )
		{
			int a0 = nattach + m * aq_sites_per_nmer;
			int a1 = nattach + (m+1) * aq_sites_per_nmer;

			for( int aq = 0; aq < aq_sites_per_nmer-1; aq++ )
			{
				bond_list[(m*bonds_per_nmer)*2+2*nb+0] = a0+aq+1;
				bond_list[(m*bonds_per_nmer)*2+2*nb+1] = a1+aq;
				if( bond_r )
				{
					switch( aq )
					{
						case 0:
							bond_r[m*bonds_per_nmer+nb] = 51.656586;
							break;
						case 1:
							bond_r[m*bonds_per_nmer+nb] = 50.143111;
							break;
						case 2:
							bond_r[m*bonds_per_nmer+nb] = 38.694115;
							break;
						case 3:
							bond_r[m*bonds_per_nmer+nb] = 40.963655;
							break;
					}
				}

				if( bond_k_in ) bond_k_in[m*bonds_per_nmer+nb] = bond_k;

				nb++;
			}
		}
	}

}
