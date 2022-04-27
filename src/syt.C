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
#include "gather.h"
#include "io_mol_read.h"

#define SQRT_SAFETY (1e-7)

#define EPS_SMALL (1e-14)
#define WORKING

static double scale_k = 0.01;
	
static double attach_bond_k = 1.0 * scale_k;
static double attach_bond_2_k = 0.00001 * scale_k;
static double inter_bond_k = 1.0 * scale_k;
//static double dihedral_k = 1000.0 * scale_k;
static double dihedral_k = 100.0 * scale_k;
static double angle_k = 300.0 * scale_k;

// if one attachment is missing these will make a better structure... just for simulation building.
static double angle_2_k = 0.0 * scale_k;

#define OLD_NRM

#ifdef OLD_NRM
static double k_nrm = 1.0;
#else
static double k_nrm = 1.0;
#endif
static double syt7_phi0 = 0; // M_PI/3; //M_PI;

static int nbonds=6;
static int nangles=6;
static int ndihe=0;


static int bonds[6][2] = {
	{0,4},
	{1,4},
	{2,5},
	{3,5},
	{4,5},
	{0,2}
};	

static int angles[6][3] = {
	{0,4,1},
	{0,4,5},
	{1,4,5},
	{2,5,3},
	{2,5,4},
	{3,5,4}
};

static int dihedrals[2][4] = {
	{ 0,4,5,2 },
	{ 3,5,4,1 },
};

static double dihedral_target[2] = { M_PI/4, M_PI  }; 

static double monomer_MW = 10*26492.14; //amu 
static double C2_P_RADIUS = 25.0;
static double attach_p_radius = 5.0;
static double bond_length_attach_long = 45.0;
static double bond_length_attach_short = 35.0;
static double bond_length_inter = 40.0;
static double bond_length_inter_attach = 35.0;

static double bond_k[6][2] =
{
	{ bond_length_attach_long, attach_bond_k },	
	{ bond_length_attach_short, attach_bond_k },	
	{ bond_length_attach_long, attach_bond_k },	
	{ bond_length_attach_short, attach_bond_2_k },	
	{ bond_length_inter, inter_bond_k },
	{ bond_length_inter_attach, inter_bond_k }
};

static double all_angle_k[6] =
{
	  angle_k, 
	  angle_k, 
	  angle_k, 
	  angle_2_k, 
	  angle_k, 
	  angle_2_k 
};

void syt7::init( double *r )
{
	// aqueous initialization.

	// syt7 model has six sites
	

//	printf("Currently Syt7 can only be initialized on the membrane.\n");
//	exit(1);

	base_init();

	nsites = 6;
	nattach = 0;

	alloc();

	double vdw_r = 25.0;
	sigma[0] = vdw_r;
	sigma[1] = vdw_r;
	sigma[2] = vdw_r;
	sigma[3] = vdw_r;
	sigma[4] = vdw_r;
	sigma[5] = vdw_r;

	mass[0] = monomer_MW*100;
	mass[1] = monomer_MW*100;
	mass[2] = monomer_MW*100;
	mass[3] = monomer_MW*100;
	mass[4] = monomer_MW*100;
	mass[5] = monomer_MW*100;

	bound = 0;
	
	// initial geometry.

}

// initialize the BAR domain on the membrane.

void syt7::init( Simulation *theSimulation, surface *theSurface, double *rsurf, int f, double u, double v, int nmer )
{
	nmer_saved = nmer;
	// assume for now this is one of the points on the membrane neck.

	base_init();

	nsites = 6;
	nattach = 4;

	alloc();

	double vdw_r = 25.0;
	sigma[0] = vdw_r;
	sigma[1] = vdw_r;
	sigma[2] = vdw_r;
	sigma[3] = vdw_r;
	sigma[4] = vdw_r;
	sigma[5] = vdw_r;

	mass[0] = monomer_MW*100;
	mass[1] = monomer_MW*100;
	mass[2] = monomer_MW*100;
	mass[3] = monomer_MW*100;
	mass[4] = monomer_MW*100;
	mass[5] = monomer_MW*100;

	bound = 1;
	
	double rpt_attach1[3], nrm_attach1[3];
	theSurface->evaluateRNRM( f, u, v, rpt_attach1, nrm_attach1, rsurf );
	printf("attach1: %le %le %le\n", rpt_attach1[0], rpt_attach1[1], rpt_attach1[2] );
	double rp[3];
	double nrm[3];
	theSurface->evaluateRNRM( f, u, v, rp, nrm, rsurf);
	printf("rp: %le %le %le nrm: %le %le %le\n", rp[0], rp[1], rp[2], nrm[0], nrm[1], nrm[2] );
	double the_sign = rp[0]*nrm[0]+rp[1]*nrm[1]+rp[2]*nrm[2];
	int reverse = 1;
	if( the_sign < 0 )
		reverse = -1;
	nrm_attach1[0] *= reverse;
	nrm_attach1[1] *= reverse;
	nrm_attach1[2] *= reverse;
	nrm[0] *= reverse;
	nrm[1] *= reverse;
	nrm[2] *= reverse;
	sid[0] = theSurface->surface_id;
	sid[1] = theSurface->surface_id;
	sid[2] = theSurface->surface_id;
	sid[3] = theSurface->surface_id;
	
	//
	// This is the neck-attachment point of the first protein.
	// goes into site-position "0"

	puv[0] = u;
	puv[1] = v;
	fs[0] = f;

	// find the other points on the membrane.
	// assume we are on a saddle, in which case we will move along the positive curvature direction to place the second point.


	double dr_u[3];
	theSurface->ru( f, u, v, rsurf, dr_u );
	double drdu = normalize(dr_u);
	double dr_v[3];
	theSurface->rv( f, u, v, rsurf, dr_v );
	double drdv = normalize(dr_v);

	double c_vec_1[3];
	double c_vec_2[3];
	double c_val1, c_val2;
	double k;
	theSurface->c( f, u, v, rsurf, &k, c_vec_1, c_vec_2, &c_val1, &c_val2 ); 
	c_val1 *= reverse;
	c_val2 *= reverse;	


	//
	// Find the neck-attachment point of the second protein.
	// goes into site-position "2"

	double d_move_1[2] = { c_vec_1[0], c_vec_1[1] };
	double d_move_2[2] = { c_vec_2[0], c_vec_2[1] };

	if( c_val2 > c_val1 )
	{
		d_move_1[0] = c_vec_2[0];
		d_move_1[1] = c_vec_2[1];
		d_move_2[0] = c_vec_1[0];
		d_move_2[1] = c_vec_1[1];
	} 
	
	double d_move_expec[3] = {0,0,0};
	
	d_move_expec[0] += d_move_1[0] * dr_u[0]*drdu;
	d_move_expec[1] += d_move_1[0] * dr_u[1]*drdu;
	d_move_expec[2] += d_move_1[0] * dr_u[2]*drdu;
	
	d_move_expec[0] += d_move_1[1] * dr_v[0]*drdv;
	d_move_expec[1] += d_move_1[1] * dr_v[1]*drdv;
	d_move_expec[2] += d_move_1[1] * dr_v[2]*drdv;
	
	// save this to see if we need to reorient the other normal.
	double save_lat[3];
	memcpy( save_lat, d_move_expec, sizeof(double) * 3 );

	double len = normalize(d_move_expec);
	double r_move = bond_length_inter;

	d_move_1[0] *= r_move/len;
	d_move_1[1] *= r_move/len;

	int f_1 = f, nf = f;
	double uv1[2] = { u, v };
	double duv1[2] = { d_move_1[0], d_move_1[1] }; 
	
	do {
		f_1 = nf;
		nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf ); 
	
		double rp_t[3];
		double nrm_t[3];
		theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], rp_t, nrm_t, rsurf);
		printf("Moving to: %le %le %le\n", rp_t[0], rp_t[1], rp_t[2] );
	} while( nf != f_1 );
	
	// the neck-attachment point of the second protein.
	puv[4] = uv1[0];
	puv[5] = uv1[1];
	fs[2] = f_1;
	
	double rp_2[3];
	double nrm_2[3];
	theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], rp_2, nrm_2, rsurf);
	nrm_2[0] *= reverse;
	nrm_2[1] *= reverse;
	nrm_2[2] *= reverse;

	//
	// find the upper-leaflet attachment point of the first protein
	//	

	uv1[0] = u;
	uv1[1] = v;

	duv1[0] = d_move_2[0];
	duv1[1] = d_move_2[1];

	
	
	d_move_expec[0] = duv1[0] * dr_u[0]*drdu;
	d_move_expec[1] = duv1[0] * dr_u[1]*drdu;
	d_move_expec[2] = duv1[0] * dr_u[2]*drdu;
	
	d_move_expec[0] += duv1[1] * dr_v[0]*drdv;
	d_move_expec[1] += duv1[1] * dr_v[1]*drdv;
	d_move_expec[2] += duv1[1] * dr_v[2]*drdv;
	
	double test_v[3];

	cross(d_move_expec, save_lat, test_v );

	double dp = test_v[0] * nrm_attach1[0] + test_v[1] * nrm_attach1[1] + test_v[2] * nrm_attach1[2];

	if( dp > 0 )
	{
		duv1[0] *= -1;
		duv1[1] *= -1;
		d_move_expec[0] *= -1;
		d_move_expec[1] *= -1;
		d_move_expec[2] *= -1;
	}

	len = normalize(d_move_expec);
	r_move = bond_length_attach_short;
	// save the attachment displacement: we want the second-protein attachment displacement to be in the opposite direction (approximately).
	double attach_displacement_1[3] = { d_move_expec[0], d_move_expec[1], d_move_expec[2] };

	duv1[0] *= (r_move/len);
	duv1[1] *= (r_move/len);
	
	f_1 = f;
	nf = f;
	
	do {
		f_1 = nf;
		nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf ); 
	} while( nf != f_1 );
	
	// the upper-leaflet attachment point of the first protein.
	puv[2] = uv1[0];
	puv[3] = uv1[1];
	fs[1] = f_1;

	//
	// find the lower-leaflet attachment point of the second protein
	//	

	uv1[0] = puv[4];
	uv1[1] = puv[5];
	f_1 = fs[2];
	nf = fs[2];
	
	double dr_u_2[3];
	theSurface->ru( f_1, uv1[0], uv1[1], rsurf, dr_u_2 );
	double drdu2 = normalize(dr_u_2);
	double dr_v_2[3];
	theSurface->rv( f_1, uv1[0],uv1[1], rsurf, dr_v_2 );
	double drdv2 = normalize(dr_v_2);

	double c_vec_1_2[3];
	double c_vec_2_2[3];
	double c_val1_2, c_val2_2;
	double k_2;
	theSurface->c( f_1, uv1[0], uv1[1], rsurf, &k_2, c_vec_1_2, c_vec_2_2, &c_val1_2, &c_val2_2 ); 
	c_val1_2 *= reverse;
	c_val2_2 *= reverse;
	// choose the negative curvature direction.

	duv1[0] = c_vec_1_2[0];	
	duv1[1] = c_vec_1_2[1];	

	if( c_val1_2 > c_val2_2 )
	{
		duv1[0] = c_vec_2_2[0];	
		duv1[1] = c_vec_2_2[1];	
	}
	
	d_move_expec[0] = duv1[0] * dr_u_2[0]*drdu2;
	d_move_expec[1] = duv1[0] * dr_u_2[1]*drdu2;
	d_move_expec[2] = duv1[0] * dr_u_2[2]*drdu2;
	
	d_move_expec[0] += duv1[1] * dr_v_2[0]*drdv2;
	d_move_expec[1] += duv1[1] * dr_v_2[1]*drdv2;
	d_move_expec[2] += duv1[1] * dr_v_2[2]*drdv2;

	len = normalize(d_move_expec);
	double dpx =
		 d_move_expec[0] * attach_displacement_1[0] + 
		 d_move_expec[1] * attach_displacement_1[1] +
		 d_move_expec[2] * attach_displacement_1[2];

	if( dpx > 0 )
	{
		// wrong direction, send the other way.

		duv1[0] *= -1;
		duv1[1] *= -1;
	} 

	r_move = bond_length_attach_short;
	
	duv1[0] *= (r_move/len);
	duv1[1] *= (r_move/len);
	
	do {
		f_1 = nf;
		nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf ); 
	} while( nf != f_1 );

	puv[6] = uv1[0];
	puv[7] = uv1[1];
	fs[3] = f_1;

	// now for the ``aqueous'' sites.  they will be projected off the normal of the attachment sites.

	double io_sign = (is_inside ? 1 : -1 );

	rall[4*3+0] = rpt_attach1[0] + nrm_attach1[0] * bond_length_attach_long * io_sign;	
	rall[4*3+1] = rpt_attach1[1] + nrm_attach1[1] * bond_length_attach_long * io_sign;	
	rall[4*3+2] = rpt_attach1[2] + nrm_attach1[2] * bond_length_attach_long * io_sign;	
	
	double rpt_attach2[3], nrm_attach2[3];
	theSurface->evaluateRNRM( fs[2], puv[4], puv[5], rpt_attach2, nrm_attach2, rsurf );

	nrm_attach2[0] *= reverse;
	nrm_attach2[1] *= reverse;
	nrm_attach2[2] *= reverse;

	rall[5*3+0] = rpt_attach2[0] + nrm_attach2[0] * bond_length_attach_long * io_sign;	
	rall[5*3+1] = rpt_attach2[1] + nrm_attach2[1] * bond_length_attach_long * io_sign;	
	rall[5*3+2] = rpt_attach2[2] + nrm_attach2[2] * bond_length_attach_long * io_sign;	

	memcpy( grad_fs, fs, sizeof(int) * 4 );
	memcpy( grad_puv, puv, sizeof(double) * 4*2 );

	memset( PBC_ext, 0, sizeof(double) * 3*6 );
	
	setrall(theSimulation);
	
	// wrap them all to one atom

	for( int a = 1; a < 6; a++ )
	{
		double dr_r[3] = { rall[3*a+0] - rall[0], rall[3*a+1] - rall[1], rall[3*a+2] - rall[2] };

		double put[3];

		MinImage3D( dr_r, theSimulation->PBC_vec, put, theSimulation->alpha );

		PBC_ext[3*a+0] = put[0];
		PBC_ext[3*a+1] = put[1];
		PBC_ext[3*a+2] = put[2];
	}

}

// custom orient procedure.



void syt7::bind( int f, double u, double v)
{
	bound = 1;
}

void syt7::unbind( void )
{
	bound = 0;
}


void syt7::loadParams( parameterBlock *block )
{
}


int syt7::getNBonds( void )
{
	return 5;
}

void syt7::putBonds( int *bond_list, double *r, double *k)
{
	bond_list[0] = 0;
	bond_list[1] = 4;
	
	bond_list[2] = 1;
	bond_list[3] = 4;
	
	bond_list[4] = 2;
	bond_list[5] = 5;
	
	bond_list[6] = 3;
	bond_list[7] = 5;

	bond_list[8] = 4;
	bond_list[9] = 5;
	
	bond_list[10] = 0;
	bond_list[11] = 2;
}

double syt7::V( Simulation *theSimulation  )
{
	// normal potential? match 0-4 and 2-5 vectors to normal with angle potential

	double *alphas = theSimulation->alpha;

	double r[3*nsites];
	double n[3*nsites];

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
		
		memcpy(r+4*3 , rall+4*3, sizeof(double) * 3*2 );
		// zero nrm for aqueous
		memset(n+4*3 , 0, sizeof(double) * 3*2 );
	}
	else
	{
		memcpy(r , rall, sizeof(double) * 3*nsites );
		memset(n, 0, sizeof(double) * 3*nsites ); 
	}

	double pot = 0;

	// normal PMF potential for 0-4 and 2-5 vectors.
	double dr_05[3] = { r[5*3+0] - r[2*3+0], r[5*3+1] - r[2*3+1], r[5*3+2] - r[2*3+2] };
	double dr_04[3] = { r[4*3+0] - r[0*3+0], r[4*3+1] - r[0*3+1], r[4*3+2] - r[0*3+2] };
	double ln1 = normalize(dr_04);
	double ln2 = normalize(dr_05);
	double dp1 = dr_04[0] * n[0] + dr_04[1] * n[1] + dr_04[2] * n[2];
	double dp2 = dr_05[0] * n[2*3+0] + dr_05[1] * n[2*3+1] + dr_05[2] * n[2*3+2];

#ifdef OLD_NRM
	double dp_target = cos(syt7_phi0);
	if( is_inside )
		dp_target *= -1;
		
	pot += 0.5 * k_nrm * (dp1-dp_target)*(dp1-dp_target);
	pot += 0.5 * k_nrm * (dp2-dp_target)*(dp2-dp_target);

#else
	if( is_inside )
	{
		double phi1 = acos(dp1);
		double phi2 = acos(dp2);
	
		pot += 0.5 * k_nrm * (phi1 - syt7_phi0) * (phi1 - syt7_phi0);
		pot += 0.5 * k_nrm * (phi2 - syt7_phi0) * (phi2 - syt7_phi0);
	}
	else
	{
		double phi1 = acos(-dp1);
		double phi2 = acos(-dp2);
	
		pot += 0.5 * k_nrm * (phi1 - syt7_phi0) * (phi1 - syt7_phi0);
		pot += 0.5 * k_nrm * (phi2 - syt7_phi0) * (phi2 - syt7_phi0);
	}
#endif
	for( int b = 0; b < nbonds; b++ )
	{
		int a1 = bonds[b][0];
		int a2 = bonds[b][1];

		double r0 = bond_k[b][0];
		double k = bond_k[b][1];

		double dr[3] = { r[3*a1+0] - r[3*a2+0],
				 r[3*a1+1] - r[3*a2+1],
				 r[3*a1+2] - r[3*a2+2] };
		double ln = normalize(dr);

		pot += 0.5 * k * (ln - r0) * (ln-r0);	
	}
	
	for( int a = 0; a < nangles; a++ )
	{
		int a1 = angles[a][0];
		int a2 = angles[a][1];
		int a3 = angles[a][2];
		double l_angle_k = all_angle_k[a];

		double dr1[3] = { r[3*a1+0] - r[3*a2+0],
				 r[3*a1+1] - r[3*a2+1],
				 r[3*a1+2] - r[3*a2+2] };
		double dr2[3] = { r[3*a3+0] - r[3*a2+0],
				 r[3*a3+1] - r[3*a2+1],
				 r[3*a3+2] - r[3*a2+2] };
		double ln1 = normalize(dr1);
		double ln2 = normalize(dr2);

		double dp = dr1[0]*dr2[0]+dr1[1]*dr2[1]+dr1[2]*dr2[2];
		double nrm_dp = dp;
		pot += 0.5 * l_angle_k *nrm_dp*nrm_dp;	
	}

	for( int d = 0; d < ndihe; d++ )
	{
		int a1 = dihedrals[d][0];
		int a2 = dihedrals[d][1];
		int a3 = dihedrals[d][2];
		int a4 = dihedrals[d][3];
		
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
		double theta = acos(dp);

		pot += 0.5 * dihedral_k * (theta - dihedral_target[d] ) * (theta - dihedral_target[d] );
	}	

	return pot;
}

// gets derivative of internal energy relative to position (surfacer_g) and the normal (surfacen_g).

double syt7::grad( Simulation *theSimulation,  double *surfacer_g, double *surfacen_g )
{
	double *alphas = theSimulation->alpha;
	double r[3*nsites];
	double n[3*nsites];

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
		
		memcpy(r+4*3 , rall+4*3, sizeof(double) * 3*2 );
	}

	double pot = 0;
	
	// normal PMF potential for 0-4 and 2-5 vectors.
	double dr_05[3] = { r[5*3+0] - r[2*3+0], r[5*3+1] - r[2*3+1], r[5*3+2] - r[2*3+2] };
	double dr_04[3] = { r[4*3+0] - r[0*3+0], r[4*3+1] - r[0*3+1], r[4*3+2] - r[0*3+2] };
	double ln1 = normalize(dr_04);
	double ln2 = normalize(dr_05);
	double dp1 = dr_04[0] * n[0] + dr_04[1] * n[1] + dr_04[2] * n[2];
	double dp2 = dr_05[0] * n[2*3+0] + dr_05[1] * n[2*3+1] + dr_05[2] * n[2*3+2];


#ifdef OLD_NRM
	double dp_target = cos(syt7_phi0);
	if( is_inside )
		dp_target *= -1;
	double del1 = (dp1-dp_target);
	double del2 = (dp2-dp_target);

	pot += 0.5 * k_nrm * (del1)*(del1);
	pot += 0.5 * k_nrm * (del2)*(del2);

	for( int c1 = 0; c1 < 3; c1++ )
	for( int c2 = 0; c2 < 3; c2++ )
	{
		if( c1 == c2 )
		{
			surfacer_g[3*4+c1] += k_nrm*(del1) * n[0*3+c2] * (1/ln1);
			surfacer_g[3*0+c1] -= k_nrm*(del1) * n[0*3+c2] * (1/ln1);

			surfacer_g[3*5+c1] += k_nrm*(del2) * n[2*3+c2] * (1/ln2);
			surfacer_g[3*2+c1] -= k_nrm*(del2) * n[2*3+c2] * (1/ln2);
		}

		surfacer_g[3*4+c1] += k_nrm*(del1) * n[0*3+c2]*(-(dr_04[c1]*dr_04[c2]/ln1));		
		surfacer_g[3*0+c1] -= k_nrm*(del1) * n[0*3+c2]*(-(dr_04[c1]*dr_04[c2]/ln1));		
		surfacer_g[3*5+c1] += k_nrm*(del2) * n[2*3+c2]*(-(dr_05[c1]*dr_05[c2]/ln2));		
		surfacer_g[3*2+c1] -= k_nrm*(del2) * n[2*3+c2]*(-(dr_05[c1]*dr_05[c2]/ln2));		
	}
	
	for( int c1 = 0; c1 < 3; c1++ )
	{
		surfacen_g[0*3+c1] += k_nrm*(del1) * dr_04[c1];
		surfacen_g[2*3+c1] += k_nrm*(del2) * dr_05[c1];
	}
#else	
	double phi1,phi2;
	double fsign=1;

	if( is_inside )
	{
		phi1 = acos(dp1);
		phi2 = acos(dp2);
		fsign=1;
	}
	else
	{
		phi1 = acos(-dp1);
		phi2 = acos(-dp2);
		fsign=-1;
	}

	double del1 = (phi1-syt7_phi0);
	double del2 = (phi2-syt7_phi0);

	pot += 0.5 * k_nrm * del1*del1;
	pot += 0.5 * k_nrm * del2*del2;
	double phi_fac1 = fsign/sqrt(1-dp1*dp1);
	double phi_fac2	= fsign/sqrt(1-dp2*dp2);

	for( int c1 = 0; c1 < 3; c1++ )
	for( int c2 = 0; c2 < 3; c2++ )
	{
		if( c1 == c2 )
		{
			// 
			surfacer_g[3*4+c1] += phi_fac1*k_nrm*(del1) * n[0*3+c2] * (1/ln1);
			surfacer_g[3*0+c1] -= phi_fac1*k_nrm*(del1) * n[0*3+c2] * (1/ln1);

			surfacer_g[3*5+c1] += phi_fac2*k_nrm*(del2) * n[2*3+c2] * (1/ln2);
			surfacer_g[3*2+c1] -= phi_fac2*k_nrm*(del2) * n[2*3+c2] * (1/ln2);
		}

		surfacer_g[3*4+c1] += phi_fac1*k_nrm*(del1) * n[0*3+c2]*(-(dr_04[c1]*dr_04[c2]/ln1));		
		surfacer_g[3*0+c1] -= phi_fac1*k_nrm*(del1) * n[0*3+c2]*(-(dr_04[c1]*dr_04[c2]/ln1));		
		surfacer_g[3*5+c1] += phi_fac2*k_nrm*(del2) * n[2*3+c2]*(-(dr_05[c1]*dr_05[c2]/ln2));		
		surfacer_g[3*2+c1] -= phi_fac2*k_nrm*(del2) * n[2*3+c2]*(-(dr_05[c1]*dr_05[c2]/ln2));		
	}
	
	for( int c1 = 0; c1 < 3; c1++ )
	{
		surfacen_g[0*3+c1] += phi_fac2*k_nrm*(del1) * dr_04[c1];
		surfacen_g[2*3+c1] += phi_fac2*k_nrm*(del2) * dr_05[c1];
	}
#endif
	for( int b = 0; b < nbonds; b++ )
	{
		int a1 = bonds[b][0];
		int a2 = bonds[b][1];

		double r0 = bond_k[b][0];
		double k = bond_k[b][1];

		double dr[3] = { r[3*a1+0] - r[3*a2+0],
				 r[3*a1+1] - r[3*a2+1],
				 r[3*a1+2] - r[3*a2+2] };
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
		int a1 = angles[a][0];
		int a2 = angles[a][1];
		int a3 = angles[a][2];
		double l_angle_k = all_angle_k[a];

		double dr1[3] = { r[3*a1+0] - r[3*a2+0],
				 r[3*a1+1] - r[3*a2+1],
				 r[3*a1+2] - r[3*a2+2] };
		double dr2[3] = { r[3*a3+0] - r[3*a2+0],
				 r[3*a3+1] - r[3*a2+1],
				 r[3*a3+2] - r[3*a2+2] };
		double ln1 = normalize(dr1);
		double ln2 = normalize(dr2);

		double dp = dr1[0]*dr2[0]+dr1[1]*dr2[1]+dr1[2]*dr2[2];

		// numerator derivative

		surfacer_g[3*a1+0] += l_angle_k * dp * dr2[0] / ln1;
		surfacer_g[3*a1+1] += l_angle_k * dp * dr2[1] / ln1;
		surfacer_g[3*a1+2] += l_angle_k * dp * dr2[2] / ln1;
		
		surfacer_g[3*a2+0] += l_angle_k * dp * (-dr2[0]) / ln1;
		surfacer_g[3*a2+1] += l_angle_k * dp * (-dr2[1]) / ln1;
		surfacer_g[3*a2+2] += l_angle_k * dp * (-dr2[2]) / ln1;
		
		surfacer_g[3*a2+0] += l_angle_k * dp * (-dr1[0]) / ln2;
		surfacer_g[3*a2+1] += l_angle_k * dp * (-dr1[1]) / ln2;
		surfacer_g[3*a2+2] += l_angle_k * dp * (-dr1[2]) / ln2;
		
		surfacer_g[3*a3+0] += l_angle_k * dp * (dr1[0]) / ln2;
		surfacer_g[3*a3+1] += l_angle_k * dp * (dr1[1]) / ln2;
		surfacer_g[3*a3+2] += l_angle_k * dp * (dr1[2]) / ln2;

		// denominator derivative.
		
		double ln1_2 = ln1*ln1;
		double ln2_2 = ln2*ln2;

		surfacer_g[3*a1+0] += -l_angle_k * dp * dp * dr1[0] / ln1;
		surfacer_g[3*a1+1] += -l_angle_k * dp * dp * dr1[1] / ln1;
		surfacer_g[3*a1+2] += -l_angle_k * dp * dp * dr1[2] / ln1;
		
		surfacer_g[3*a2+0] += l_angle_k * dp * dp * dr1[0] / ln1;
		surfacer_g[3*a2+1] += l_angle_k * dp * dp * dr1[1] / ln1;
		surfacer_g[3*a2+2] += l_angle_k * dp * dp * dr1[2] / ln1;
		
		surfacer_g[3*a2+0] += l_angle_k * dp * dp * dr2[0] / ln2;
		surfacer_g[3*a2+1] += l_angle_k * dp * dp * dr2[1] / ln2;
		surfacer_g[3*a2+2] += l_angle_k * dp * dp * dr2[2] / ln2;
		
		surfacer_g[3*a3+0] += -l_angle_k * dp * dp * dr2[0] / ln2;
		surfacer_g[3*a3+1] += -l_angle_k * dp * dp * dr2[1] / ln2;
		surfacer_g[3*a3+2] += -l_angle_k * dp * dp * dr2[2] / ln2;

		pot += 0.5 *l_angle_k * dp*dp;	
	}

	for( int d = 0; d < ndihe; d++ )
	{
		int a1 = dihedrals[d][0];
		int a2 = dihedrals[d][1];
		int a3 = dihedrals[d][2];
		int a4 = dihedrals[d][3];
		double phi0 = dihedral_target[d];		

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
		double theta = acos(dp);

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
		double base_der = dihedral_k * (theta - phi0) * (-1 / sqrt(1-dp*dp));
	
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
	
		pot += 0.5 * dihedral_k * (theta - phi0) * (theta - phi0);
	}	

	return pot;
}

void syt7::move_inside( void )
{
	pcomplex::move_inside();
}

void syt7::move_outside( void )
{
	pcomplex::move_outside();
}


void syt7::writeStructure( Simulation *theSimulation, surface_mask *upperSurfaceMask, surface_mask *lowerSurfaceMask, struct atom_rec **at_out, int *nat_out, char ***sequence, int *nseq_out, int **seq_at_array, char ***patches, struct ion_add **ions, int *nions, 
	aa_build_data *buildData, int *build_type )
{
	*build_type = BUILD_SEQUENCE;

	int nsegments = 2;
	
	surface *theSurface;
	double *rsurf;
	theSimulation->fetch(sid[0],&theSurface,&rsurf);

	const char *domain_names[] = { "C2A", "C2B" };
	const char *segment_names[] = { "PROA", "PROC" };
	int segment_surface_sites[][2] = { {0, -1}, {2, -1} };
	//int segment_surface_sites[][2] = { {0, 1}, {2, 3} };
	int segment_aqueous_sites[][2] = { {4, 5}, {5, 4}  };
	int segment_residues[][3] = {    { 227, 222, 262},
					{ 359, 305, 267} };
	
	int a_start = *nat_out;
	
	int a_markers[2] = { 0,0};
	double com[2][3] = { {0,0,0}, {0,0,0} };
	int ncom[2] = {0,0};
	for( int seg = 0; seg < 2; seg++ )
	{
		a_markers[seg] = *nat_out;

		struct atom_rec *at;
		int nat;

		int pool_code = pdbFetch( &at, &nat, "syt7", domain_names[seg], addToPool );

		int ntries = 0;
		double random_noise = 0.01;

		int clash = 1;
		// addPeri... returns one if there's a clash.. keep trying, injecting more noise
		while( (clash = addPeripheralProteinHelper( theSimulation, upperSurfaceMask, lowerSurfaceMask,
			at_out, nat_out,
			ions, nions,	
			pool_code,
			segment_names[seg],
			segment_residues[seg],	
			segment_surface_sites[seg],	
			segment_aqueous_sites[seg],	
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
		

		for( int a = a_markers[seg]; a < *nat_out; a++ )
		{
			com[seg][0] += (*at_out)[a].x;
			com[seg][1] += (*at_out)[a].y;
			com[seg][2] += (*at_out)[a].z;
			ncom[seg] += 1;
		}
	
		com[seg][0] /= ncom[seg];
		com[seg][1] /= ncom[seg];
		com[seg][2] /= ncom[seg];
	}


	double dr[3] = { com[1][0] - com[0][0], com[1][1] - com[0][1], com[1][2] - com[0][2] };
	double put[3];

	MinImage3D( dr, theSimulation->PBC_vec, put ); 

	for( int a = a_markers[1]; a < *nat_out; a++ )
	{
		(*at_out)[a].x += put[0] * theSimulation->PBC_vec[0][0] + put[1] * theSimulation->PBC_vec[1][0] + put[2] * theSimulation->PBC_vec[2][0];
		(*at_out)[a].y += put[0] * theSimulation->PBC_vec[0][1] + put[1] * theSimulation->PBC_vec[1][1] + put[2] * theSimulation->PBC_vec[2][1];
		(*at_out)[a].z += put[0] * theSimulation->PBC_vec[0][2] + put[1] * theSimulation->PBC_vec[1][2] + put[2] * theSimulation->PBC_vec[2][2];
	}

/*	double del[3] = { 0,0,0};

	while( dr[0] + del[0] > La/2 ) del[0] -= La;
	while( dr[0] + del[0] <-La/2 ) del[0] += La;
	while( dr[1] + del[1] > Lb/2 ) del[1] -= Lb;
	while( dr[1] + del[1] <-Lb/2 ) del[1] += Lb;
	while( dr[2] + del[2] > Lc/2 ) del[2] -= Lc;
	while( dr[2] + del[2] <-Lc/2 ) del[2] += Lc;
*/	struct atom_rec *at = *at_out;

	// make sure the domains are in the same PBC


	int pres = at[a_start].res;
	
	int seq_space = 10;
	int nseq = 0;

	char *seq = (char *)malloc( sizeof(char) * seq_space );

	seq[nseq] = threeToOne( at[a_start].resname );
	nseq++;

	for( int a = a_start; a < *nat_out; a++ )
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

			for( int t = nseq; t < nseq-1+doff; t++ )
				seq[t] = 'X';		

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

	for( int t = 0; t < nseq; t++ )
	{
		if( seq[t] == 'X' )
		{
			if( nmissing >= ntodo )
			{
				printf("ERROR. found a missing residue but didn't have enough to insert.\n");
				exit(1);
			}

			seq[t] = missing_residues[nmissing];
			nmissing++;
		}
	}
	*sequence = (char **)malloc( sizeof(char *) * 1 );
	(*sequence)[0] = seq; 
	(*seq_at_array) = (int *)malloc( sizeof(int) * 1 );
	(*seq_at_array)[0] = 0;

	*patches = (char **)malloc( sizeof( char *) * 1 );
	(*patches)[0] = NULL;
}

/* this function 
 * "gets" the coarse-grained coordinates of the observed protein from the molecular structure -- and its fit to the membrane surface.
 *
 * */


void syt7::get( 
		Simulation *theSimulation, // this is the surface/PBC
		struct atom_rec *at,
		int syt7_start, int syt7_stop, int nat_tot )
{
	int nsegments = 2;
	
	surface *theSurface;
	double *rsurf;
	theSimulation->fetch(sid[0],&theSurface,&rsurf);
	surface_record *sRec = theSimulation->fetch(sid[0]);

	int flip_sign = sRec->gather_flip;
	double Ls[3] = {-1,-1,-1};

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
	double rvecB[3] = { rall[5*3+0] - p[6],
			    rall[5*3+1] - p[7],
			    rall[5*3+2] - p[8] };
	normalize(rvecA);
	normalize(rvecB);

	double dpA = rvecA[0]*nrm1[0] + rvecA[1] * nrm1[1] + rvecA[2] * nrm1[2];
	double dpB = rvecB[0]*nrm3[0] + rvecB[1] * nrm3[1] + rvecB[2] * nrm3[2];
	

	for( int c = 0; c < 3; c++ )
		Ls[c] = theSimulation->PBC_vec[c][c];
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

	// Orientation of the special binding domains.
	
	// first: is the domain bound to negatively-charged lipids?
	
	int C2_in_position[2] = { 0,0 };

	int strand_min_res[2] = { 179, 310 };
	int strand_max_res[2] = { 196, 329 }; 

	int pos_res[2][10] = {
			    { 183, 184, 186, 190, 192, 193, 194 }, 
			    { 313, 315, 316, 319, 320, 321 }
				};
	int sizes[2] = { 7, 6 };

	double avp[2][3] = { { 0,0,0},
			     { 0,0,0} };
	double navp[2] = {0,0};
	double I[2][9];
	int nbound[2] = {0,0};	
	int nbound_all[2] = { 0,0 };
	int nbound_ca[2] = { 0,0 };

	double domain_com[2][3] = { {0,0,0},
				    {0,0,0} };
	double ndom[2] = {0,0};

        int domain_start[2] = { 1, 255 };
        int domain_stop[2] = { 256, 500 }; // arbitrary end.

	int res_is_interacting[501];
	memset( res_is_interacting, 0, sizeof(int) * 501 );

	int nat2_space = 10;
	int nat2 = 0;
	int *a2_loop = (int *)malloc( sizeof(int) * nat2_space );

	double *depth = (double *)malloc( sizeof(double) * 500 );
	for( int r = 0; r < 500; r++)
		depth[r] = -1;

	int resv[2][2] = { {-1,-1},{-1,-1} };

	for( int a2 = syt7_start; a2 < syt7_stop ; a2++ )
	{
		if( !strcasecmp( at[a2].resname, "ASP") && !strcasecmp( at[a2].atname, "CA" ) )
		{
			if( at[a2].res == 303 )
				resv[1][1] = a2;
			else if( at[a2].res == 365 )
				resv[1][0] = a2;
			else if( at[a2].res == 233 )
				resv[0][1] = a2;
			else if( at[a2].res == 172 )
				resv[0][0] = a2;
		}
	}


	for( int a2 = 0; a2 < nat_tot; a2++ )
	{
		if( strncasecmp( at[a2].resname, "SAPI", 4) && strncasecmp( at[a2].resname, "PAPS", 4 ) )
			continue;
	
		if( at[a2].atname[0] != 'O' && at[a2].atname[0] != 'P' && strcasecmp( at[a2].atname, "C13") ) continue; // look for charge-bearing oxygens.
	
		if( nat2 == nat2_space )
		{
			nat2_space *= 2;
			a2_loop = (int *)realloc( a2_loop, sizeof(int) * nat2_space );
		}

		a2_loop[nat2] = a2;
		nat2++;	
	
	}

	for( int domain = 0; domain < 2; domain++ )
	{
		memset( I[domain], 0, sizeof(double) * 9 );
		
		for( int a = syt7_start; a < syt7_stop ; a++ )
		{
			if( at[a].res >= domain_start[domain] && at[a].res <= domain_stop[domain] )
			{
				domain_com[domain][0] += at[a].x;
				domain_com[domain][1] += at[a].y;
				domain_com[domain][2] += at[a].z;
				ndom[domain] += 1;	
			}

			if( at[a].res >= strand_min_res[domain] && at[a].res <= strand_max_res[domain] )
			{
				avp[domain][0] += at[a].x;
				avp[domain][1] += at[a].y;
				avp[domain][2] += at[a].z;
				navp[domain] += 1;
			}
		}

		avp[domain][0] /= navp[domain];
		avp[domain][1] /= navp[domain];
		avp[domain][2] /= navp[domain];

		for( int a = syt7_start; a < syt7_stop ; a++ )
		{
			if( at[a].res < domain_start[domain] || at[a].res > domain_stop[domain] )
				continue;

			int checkit = 0;

			if( at[a].res >= strand_min_res[domain] && at[a].res <= strand_max_res[domain] )
			{
				double dr[3] = {
					at[a].x - avp[domain][0],
					at[a].y - avp[domain][1],
					at[a].z - avp[domain][2] };
				
				for( int c1 = 0; c1 < 3; c1++ )
				for( int c2 = 0; c2 < 3; c2++ )
					I[domain][c1*3+c2] += dr[c1] * dr[c2];
			}

			for( int r = 0; r < sizes[domain]; r++ )
			{
				if( at[a].res == pos_res[domain][r]  )
					checkit = 1;
			}

			if( !strcasecmp( at[a].resname, "ASP") && at[a].atname[0] == 'O' )
			{
				for( int a2x = 0; a2x < nat2; a2x++ )
				{
					int a2 = a2_loop[a2x];

					if( strncasecmp( at[a2].resname, "SAPI", 4) && strncasecmp( at[a2].resname, "PAPS", 4 ) )
						continue;
	
					if( at[a2].atname[0] != 'O') continue; // look for charge-bearing oxygens.
	
					double dr[3] = { at[a].x - at[a2].x, at[a].y - at[a2].y, at[a].z - at[a2].z }; 
	
					theSimulation->wrapPBC(dr,alphas);
	
					double lr = normalize(dr);
	
					if( lr < 6.0 )
					{
//						res_is_interacting[at[a].res] = -1;
						nbound_ca[domain] += 1;
					}
				} 
					
			}

			if( !strcasecmp( at[a].resname, "ARG") || !strcasecmp( at[a].resname, "LYS") ) 
			{
				int arg_check = !strcasecmp( at[a].resname, "ARG" ) && ( !strcasecmp(at[a].atname, "CZ") );
				int lys_check = !strcasecmp( at[a].resname, "LYS" ) && (!strcasecmp( at[a].atname, "NZ") );

				if( lys_check || arg_check )
				{
					for( int a2x = 0; a2x < nat2; a2x++ )
					{
						int a2 = a2_loop[a2x];
						if( strncasecmp( at[a2].resname, "SAPI", 4) && strncasecmp( at[a2].resname, "PAPS", 4 ) )
							continue;
		
						if( at[a2].atname[0] != 'P' && (strcasecmp(at[a2].resname, "PAPS") || strcasecmp(at[a2].atname,"C13") )) continue; // look for charge-bearing phosphos, or carboxyls (PAPS).
		
						double dr[3] = { at[a].x - at[a2].x, at[a].y - at[a2].y, at[a].z - at[a2].z }; 
		
						theSimulation->wrapPBC(dr,alphas);
		
						double lr = normalize(dr);
		
						if( lr < 4.5 )
						{
							if( checkit )
							{
								nbound[domain] += 1;
								C2_in_position[domain] = 1; 
							}
		
							nbound_all[domain] += 1;	
							if( !strcasecmp(at[a2].resname, "PAPS") )
								res_is_interacting[at[a].res] = res_is_interacting[at[a].res] | (1<<0);
							if( !strncasecmp(at[a2].resname, "SAPI",4) )
								res_is_interacting[at[a].res] = res_is_interacting[at[a].res] | (1<<1);
						}
					}
				}
			} 
		}
	}

	for( int domain = 0; domain < 2; domain++ )
	{
		// diagonalize I.

		char jobz = 'V';
		char uplo = 'U';
		int order =3;
		double ev[3];
		int lwork = 100;	
		int info = 0;
		double work[100];
	
		dsyev(&jobz,&uplo,&order,I[domain],&order,ev,work,&lwork,&info);

		double vec[3];
		double maxev=-1;
		for( int iev = 0; iev < 3; iev++ )
		{
			if( ev[iev] > maxev )
			{
				maxev = ev[iev];
				vec[0] = I[domain][3*iev+0];
				vec[1] = I[domain][3*iev+1];
				vec[2] = I[domain][3*iev+2];
			}
		}
		normalize(vec);	

		double dist;
		int f;
		double col_u, col_v;
		theSurface->nearPointOnBoxedSurface( avp[domain], &f, &col_u, &col_v, M, mlow, mhigh, &dist );					

		double rp[3], np[3];
		theSurface->evaluateRNRM( f, col_u, col_v, rp, np, rsurf );

		double dr[3] = { rp[0] - avp[domain][0],
				 rp[1] - avp[domain][1],
				 rp[2] - avp[domain][2] };
		theSimulation->wrapPBC( dr, alphas );

		double c_vec_1[3];
		double c_vec_2[3];
		double c_val1, c_val2;
		double k;
		theSurface->c( f, col_u, col_v, rsurf, &k, c_vec_1, c_vec_2, &c_val1, &c_val2 ); 
	
		if( flip_sign) { c_val1 *= -1; c_val2 *= -1; }

		double del = normalize(dr);
	
		double drdu[3];
		theSurface->ru( f, col_u, col_v, rsurf, drdu);
		double drdv[3];
		theSurface->rv( f, col_u, col_v, rsurf, drdv );

		double c_vec_1_3[3] = { 
				c_vec_1[0]*drdu[0] + c_vec_1[1] * drdv[0],
				c_vec_1[0]*drdu[1] + c_vec_1[1] * drdv[1],
				c_vec_1[0]*drdu[2] + c_vec_1[1] * drdv[2]};
		double c_vec_2_3[3] = { 
				c_vec_2[0]*drdu[0] + c_vec_2[1] * drdv[0],
				c_vec_2[0]*drdu[1] + c_vec_2[1] * drdv[1],
				c_vec_2[0]*drdu[2] + c_vec_2[1] * drdv[2]};

		normalize(c_vec_1_3);
		normalize(c_vec_2_3);

		// orient in plane.

		double proj[3];

		double dp = dot( vec, np );

		proj[0] = vec[0] - dp * np[0];	
		proj[1] = vec[1] - dp * np[1];	
		proj[2] = vec[2] - dp * np[2];	

		normalize(proj);
	
		int n_interacting = 0, n_interacting_ca = 0;
		for( int r = domain_start[domain]; r <= domain_stop[domain]; r++ )		
		{	
			if( res_is_interacting[r] )
				n_interacting += 1;
			if( res_is_interacting[r] == -1 )
				n_interacting_ca += 1;	
		}
		printf("strand P%d C2%c distance: %lf nbound %d c1: %le c2: %le vec1: %lf %lf %lf vec2: %lf %lf %lf orientation: %lf %lf %lf phi: %lf ninteracting: %d interacting_ca: %d\n", my_id, (domain == 0 ? 'A' : 'B'), del, nbound[domain], c_val1, c_val2, c_vec_1_3[0], c_vec_1_3[1], c_vec_1_3[2], c_vec_2_3[0], c_vec_2_3[1], c_vec_2_3[2], proj[0], proj[1], proj[2], (180.0/M_PI) * acos(fabs(dp)), n_interacting, n_interacting_ca  );
			
		
	}
	
	free(a2_loop);

	double fp_center[3];
	double fp_r;
	double fp_h;

	domain_com[0][0] /= ndom[0];
	domain_com[0][1] /= ndom[0];
	domain_com[0][2] /= ndom[0];
	
	domain_com[1][0] /= ndom[1];
	domain_com[1][1] /= ndom[1];
	domain_com[1][2] /= ndom[1];

	double fp_thresh = 20.0;


	double p_CBL_1[3] = {  (at[resv[0][1]].x + at[resv[0][0]].x)/2,
			       (at[resv[0][1]].y + at[resv[0][0]].y)/2,
			       (at[resv[0][1]].z + at[resv[0][0]].z)/2 };
	double p_CBL_2[3] = {  (at[resv[1][1]].x + at[resv[1][0]].x)/2,
			       (at[resv[1][1]].y + at[resv[1][0]].y)/2,
			       (at[resv[1][1]].z + at[resv[1][0]].z)/2 };

	double dr_CBL_1[3] = { at[resv[0][1]].x - at[resv[0][0]].x,
			       at[resv[0][1]].y - at[resv[0][0]].y,
			       at[resv[0][1]].z - at[resv[0][0]].z };
	double dr_CBL_2[3] = { at[resv[1][1]].x - at[resv[1][0]].x,
			       at[resv[1][1]].y - at[resv[1][0]].y,
			       at[resv[1][1]].z - at[resv[1][0]].z };

	for( int c = 0; c < 3; c++ )
	{
		while( dr_CBL_1[c] < -Ls[c]/2 ) dr_CBL_1[c] += Ls[c]; 
		while( dr_CBL_1[c] >  Ls[c]/2 ) dr_CBL_1[c] -= Ls[c]; 
		while( dr_CBL_2[c] < -Ls[c]/2 ) dr_CBL_2[c] += Ls[c]; 
		while( dr_CBL_2[c] >  Ls[c]/2 ) dr_CBL_2[c] -= Ls[c]; 
	}

		

	if( theSimulation->getFusionPoreData( fp_center, &fp_r, &fp_h ) )
	{
		printf("FP cen: %lf %lf %lf\n", fp_center[0], fp_center[1], fp_center[2] );
		printf("FP r: %lf\n", fp_r );
		printf("FP h: %lf\n", fp_h );		

		double dr_A[3] = { domain_com[0][0] - fp_center[0], domain_com[0][1] - fp_center[1], domain_com[0][2] - fp_center[2] };
		double dr_B[3] = { domain_com[1][0] - fp_center[0], domain_com[1][1] - fp_center[1], domain_com[1][2] - fp_center[2] };

		double dist;	
		int f_com_A;
		double uv_A[2];
		theSurface->nearPointOnBoxedSurface( domain_com[0], &f_com_A, uv_A, uv_A+1, M, mlow, mhigh, &dist );		
	
		double c_vec_1[3];
		double c_vec_2[3];
		double c_val_A1, c_val_A2;
		double k;
		theSurface->c( f_com_A, uv_A[0], uv_A[1], rsurf, &k, c_vec_1, c_vec_2, &c_val_A1, &c_val_A2 ); 
		
		int f_com_B;
		double uv_B[2];
		theSurface->nearPointOnBoxedSurface( domain_com[1], &f_com_B, uv_B, uv_B+1, M, mlow, mhigh, &dist );		
	
		double c_val_B1, c_val_B2;
		theSurface->c( f_com_B, uv_B[0], uv_B[1], rsurf, &k, c_vec_1, c_vec_2, &c_val_B1, &c_val_B2 ); 
	
		double put[3];

		MinImage3D( dr_A, theSimulation->PBC_vec, put ); 
		MinImage3D( dr_B, theSimulation->PBC_vec, put ); 

		double r_A = normalize(dr_A);
		double r_B = normalize(dr_B);

		int is_FP[2] = {0,0};

		if( r_A < fp_r + fp_thresh )
			is_FP[0] = 1;
		if( r_B < fp_r + fp_thresh )
			is_FP[1] = 1;

		printf("CONTACTS %d C2A: r: %lf c: %lf %lf C2B: r: %lf c: %lf %lf res_depths:", my_id, r_A, c_val_A1, c_val_A2, r_B, c_val_B1, c_val_B2 );
		int max = -1;		
		for( int a = syt7_start; a < syt7_stop ; a++ )
		{
			if( !strcasecmp( at[a].atname, "CA") ) 
			{
				double rp[3] = { at[a].x, at[a].y, at[a].z };

				double dist;
				int f;
				double col_u, col_v;
				theSurface->nearPointOnBoxedSurface( rp, &f, &col_u, &col_v, M, mlow, mhigh, &dist );		

				double rpt[3], npt[3];

				theSurface->evaluateRNRM( f, col_u, col_v, rpt, npt, rsurf );

				double dr[3] = { rp[0] - rpt[0], rp[1] - rpt[1], rp[2] - rpt[2] };
	
				double put[3];
				MinImage3D( dr, theSimulation->PBC_vec, put );

				double d = normalize(dr);

				depth[at[a].res] = d; 
				max = at[a].res;
			}		
		}

		for( int r = 1; r <= max; r++ )
			printf(" %lf", depth[r] );
		printf("\n");
		
		printf("ARGLYS %d C2A: r: %lf c: %lf %lf C2B: r: %lf c: %lf %lf res_depths:", my_id, r_A, c_val_A1, c_val_A2, r_B, c_val_B1, c_val_B2 );
		
		for( int r = 1; r <= max; r++ )
			printf(" %d", res_is_interacting[r] );
		printf("\n");


		{
			double at_c_1,at_c_2;

			double dist;
			int f;
			double col_u, col_v;

			double c_vec_1A[2], c_vec_1B[2];		
			double c_vec_2A[2], c_vec_2B[2];		
			double c_val_1A, c_val_1B;
			double c_val_2A, c_val_2B;
			double k;
			theSurface->nearPointOnBoxedSurface( p_CBL_1, &f, &col_u, &col_v, M, mlow, mhigh, &dist );		
			theSurface->c( f, col_u, col_v, rsurf, &k, c_vec_1A, c_vec_1B, &c_val_1A, &c_val_1B ); 

			double dr_u[3];
			theSurface->ru( f, col_u, col_v, rsurf, dr_u );
//			double drdu = normalize(dr_u);
			double dr_v[3];
			theSurface->rv( f, col_u, col_v, rsurf, dr_v );
//			double drdv = normalize(dr_v);

			double RSP_CV_1A[3] = { 
				dr_u[0] * c_vec_1A[0] + dr_v[0] * c_vec_1A[1],
				dr_u[1] * c_vec_1A[0] + dr_v[1] * c_vec_1A[1],
				dr_u[2] * c_vec_1A[0] + dr_v[2] * c_vec_1A[1] };
			double RSP_CV_1B[3] = { 
				dr_u[0] * c_vec_1B[0] + dr_v[0] * c_vec_1B[1],
				dr_u[1] * c_vec_1B[0] + dr_v[1] * c_vec_1B[1],
				dr_u[2] * c_vec_1B[0] + dr_v[2] * c_vec_1B[1] };
					
			theSurface->nearPointOnBoxedSurface( p_CBL_2, &f, &col_u, &col_v, M, mlow, mhigh, &dist );		
			theSurface->ru( f, col_u, col_v, rsurf, dr_u );
			theSurface->rv( f, col_u, col_v, rsurf, dr_v );
//			drdu = normalize(dr_u);
//			drdv = normalize(dr_v);
			theSurface->c( f, col_u, col_v, rsurf, &k, c_vec_2A, c_vec_2B, &c_val_2A, &c_val_2B ); 

			double RSP_CV_2A[3] = { 
				dr_u[0] * c_vec_2A[0] + dr_v[0] * c_vec_2A[1],
				dr_u[1] * c_vec_2A[0] + dr_v[1] * c_vec_2A[1],
				dr_u[2] * c_vec_2A[0] + dr_v[2] * c_vec_2A[1] };
			double RSP_CV_2B[3] = { 
				dr_u[0] * c_vec_2B[0] + dr_v[0] * c_vec_2B[1],
				dr_u[1] * c_vec_2B[0] + dr_v[1] * c_vec_2B[1],
				dr_u[2] * c_vec_2B[0] + dr_v[2] * c_vec_2B[1] };

			normalize(RSP_CV_1A);
			normalize(RSP_CV_1B);
			normalize(RSP_CV_2A);
			normalize(RSP_CV_2B);

			// CBL
			printf("CBLs %d r1: %lf %lf %lf r2: %lf %lf %lf dr1: %lf %lf %lf dr2: %lf %lf %lf cvec1A: %lf %lf %lf cvec1B: %lf %lf %lf c1A: %lf c1B: %lf cvec2A: %lf %lf %lf cvec2B: %lf %lf %lf c2A: %lf c2B: %lf\n",
				my_id,
				p_CBL_1[0], p_CBL_1[1], p_CBL_1[2],
				p_CBL_2[0], p_CBL_2[1], p_CBL_2[2],
				dr_CBL_1[0], dr_CBL_1[1], dr_CBL_1[2],
				dr_CBL_2[0], dr_CBL_2[1], dr_CBL_2[2],
				RSP_CV_1A[0], RSP_CV_1A[1], RSP_CV_1A[2],
				RSP_CV_1B[0], RSP_CV_1B[1], RSP_CV_1B[2],
				c_val_1A, c_val_1B,
				RSP_CV_2A[0], RSP_CV_2A[1], RSP_CV_2A[2],
				RSP_CV_2B[0], RSP_CV_2B[1], RSP_CV_2B[2],
				c_val_2A, c_val_2B );
		}

	}	

	// angle data for creating a histogram.
	

	getLocalLipidComposition( theSimulation, at, syt7_start, syt7_stop, nat_tot );
	

	// CLASSIFICATION

	// * C2A near fusion pore
	
	// * C2B near fusion pore

	// * Both near fusion pore
	
	// Both in bulk region
}


char syt7::getSiteCode( int p )
{
	switch (p)
	{
		case 1:
		case 3:
			return 'N';

		default:
			return 'O';
	}
}



void syt7::getLocalLipidComposition( 
		Simulation *theSimulation, // this is the surface/PBC
		struct atom_rec *at,
		int syt7_start, int syt7_stop, int nat_tot )
{
	

	// for near pointing
	double **M;
	int mlow,mhigh;

	getM( &M, &mlow, &mhigh );
	int nsegments = 2;
	
	surface *theSurface;
	double *rsurf;
	theSimulation->fetch(sid[0],&theSurface,&rsurf);
	surface_record *sRec = theSimulation->fetch(sid[0]);

	int flip_sign = sRec->gather_flip;
	double Ls[3] = {-1,-1,-1};

	/*
 *		For each ``face'' of the protein that we're interested in, get the distance from the membrane and establish a local coordinate system based on our orientation. Then, output lipid centers-of-mass in the lipid coordinate system.
 *
 * 	the distance will be defined by a point-of-contact, and the orientation by the vector from the point-of-contact to the center of mass of the domain.
 *
 * 	*/

	// STEP 1: get domain COMs	
	
	double alphas[3] = {1,1,1};
	double domain_com[2][3] = { {0,0,0},
				    {0,0,0} };
	double ndom[2] = {0,0};

        int domain_start[2] = { 1, 255 };
        int domain_stop[2] = { 256, 500 }; // arbitrary end.
	
	for( int domain = 0; domain < 2; domain++ )
	{
		for( int a = syt7_start; a < syt7_stop ; a++ )
		{
			if( at[a].res >= domain_start[domain] && at[a].res <= domain_stop[domain] )
			{
				double cur_com[3] = {0,0,0};
		
				if( ndom[domain] > 0 )
				{
					cur_com[0] = domain_com[domain][0] / ndom[domain];
					cur_com[1] = domain_com[domain][1] / ndom[domain];
					cur_com[2] = domain_com[domain][2] / ndom[domain];
				}
		
				
				double dr[3] = { at[a].x - cur_com[0],
						 at[a].y - cur_com[1],
						 at[a].z - cur_com[2] };
					
				theSimulation->wrapPBC( dr, alphas );
				
				domain_com[domain][0] += cur_com[0] + dr[0];
				domain_com[domain][1] += cur_com[1] + dr[1];
				domain_com[domain][2] += cur_com[2] + dr[2];
				ndom[domain] += 1;	
			}
		}
	}
	
	domain_com[0][0] /= ndom[0];
	domain_com[0][1] /= ndom[0];
	domain_com[0][2] /= ndom[0];
	
	domain_com[1][0] /= ndom[1];
	domain_com[1][1] /= ndom[1];
	domain_com[1][2] /= ndom[1];



	// STEP 2: get all lipid COMs.

	// neutral surface atoms.

	int ns = io_nNSAtoms();
	struct atom_rec *at_lipid = io_getCurrentFrame();

	double *lcom = (double *)malloc( sizeof(double) * 3 * ns );
	for( int ax = 0; ax < ns; ax ++ )
	{
		int a = io_getNSAtom(ax);
			
		int astart=-1;
		int astop=-1;
		io_getNSAtomStartStop(ax,&astart,&astop); 

		lcom[3*ax+0]=0;
		lcom[3*ax+1]=0;
		lcom[3*ax+2]=0;
			
		for( int x = 0; x < astop-astart; x++ )
		{
			lcom[3*ax+0] += at[astart+x].x;
			lcom[3*ax+1] += at[astart+x].y;
			lcom[3*ax+2] += at[astart+x].z;
		}

		lcom[3*ax+0]/=(astop-astart);
		lcom[3*ax+1]/=(astop-astart);
		lcom[3*ax+2]/=(astop-astart);
	}
	
	// STEP 3: get points-of-contact.

	// number of points-of-contact:
	int nPofCon[2] = { 2, 3};
	// # of residues in PofCon:
	int nRes[2][3] = {
		{ 5, 3, 0}, // CBL, PBR_A
		{ 5, 5, 2}, /// CBL, PBR_B1, PBR_B2
	};

	int pocRes[2][3][6] =
	{
		{	// C2A
			{ 166, 172, 225, 227, 233 }, // Calcium-binding loop.
			{ 183, 184, 186 }, // PBR
		},

		{	// C2B
			{ 297, 303, 357, 359, 365 }, // calcium-binding loop.
			{ 313, 315, 316, 320, 321 }, // PBR1
			{ 390, 392 }, // PBR2
		}
	};
	
	int nfaces = 0;
	for( int domain = 0; domain < 2; domain++ )
		nfaces += nPofCon[domain];

	double pcom[2][3][3];	
	int f_ind_for_d[2][3];
	int domain_for_face[nfaces];
	int face_index[nfaces];
	double orientation[nfaces*3];
	double poc[nfaces*3];
	double near_point[nfaces*3];
	double tangent_plane[nfaces*9];		
	double dist[nfaces];

	int f_ind = 0;
	for( int domain = 0; domain < 2; domain++ )
	{
		for( int poc = 0; poc < nPofCon[domain]; poc++ )
		{
			domain_for_face[f_ind] = domain;
			face_index[f_ind] = poc;
			

			pcom[domain][poc][0] = 0;
			pcom[domain][poc][1] = 0;
			pcom[domain][poc][2] = 0;
			
			f_ind_for_d[domain][poc] = f_ind;
	
			f_ind++;
		}
	}

	for( int a = syt7_start; a < syt7_stop ; a++ )
	{
		int f_ind = 0;
		for( int domain = 0; domain < 2; domain++ )
		{
			for( int poc = 0; poc < nPofCon[domain]; poc++ )
			{
				for( int r = 0; r < nRes[domain][poc]; r++ )
				{
					if( at[a].res == pocRes[domain][poc][r] && !strcasecmp(at[a].atname, "CA") )
					{
						pcom[domain][poc][0] += at[a].x;
						pcom[domain][poc][1] += at[a].y;
						pcom[domain][poc][2] += at[a].z;
					} 
				}
				f_ind++;
			}
		}
	}

	f_ind=0;	
	for( int domain = 0; domain < 2; domain++ )
	{
		for( int tpoc = 0; tpoc < nPofCon[domain]; tpoc++ )
		{
			pcom[domain][tpoc][0] /= nRes[domain][tpoc];
			pcom[domain][tpoc][1] /= nRes[domain][tpoc];
			pcom[domain][tpoc][2] /= nRes[domain][tpoc];

			poc[f_ind*3+0] = pcom[domain][tpoc][0];	
			poc[f_ind*3+1] = pcom[domain][tpoc][1];	
			poc[f_ind*3+2] = pcom[domain][tpoc][2];	

			f_ind++;
		}
	}

	// Step 4: get each face's orientation and tangent plane coordinate system.
	
	int f_ref_CBL[2][3][2] =
	{
		{	{1,0}, {0,0} }, 	// domain 1 CBL, domain 0 CBL
		{	{0,0}, {1,0}, {1,0} }	// domain 0 CBL, domain 1 CBL, domain 1 CBL
	};
	
	f_ind=0;	
	for( int domain = 0; domain < 2; domain++ )
	{
		for( int tpoc = 0; tpoc < nPofCon[domain]; tpoc++ )
		{
			// pof - domain com
			orientation[f_ind*3+0] = poc[f_ind*3+0] - domain_com[domain][0]; 
			orientation[f_ind*3+1] = poc[f_ind*3+1] - domain_com[domain][1]; 
			orientation[f_ind*3+2] = poc[f_ind*3+2] - domain_com[domain][2]; 

			// normalized.
			normalize(orientation+f_ind*3);

			// near point on surface:
	
			// find near points for the attachment sites.
			double the_dist;
			int f;
			double col_u, col_v;
			double pt[3] = { poc[3*f_ind+0], poc[3*f_ind+1], poc[3*f_ind+2] };
			theSurface->nearPointOnBoxedSurface( pt, &f, &col_u, &col_v, M, mlow, mhigh, &the_dist );					
			double rpt[3], npt[3];
			theSurface->evaluateRNRM( f, col_u, col_v, rpt, npt, rsurf );

			double dr[3] = { 
				rpt[0] - poc[3*f_ind+0],
				rpt[1] - poc[3*f_ind+1],
				rpt[2] - poc[3*f_ind+2] };

			theSimulation->wrapPBC(dr,alphas);

			// near point , wrapped.
			near_point[3*f_ind+0] = poc[3*f_ind+0] + dr[0];
			near_point[3*f_ind+1] = poc[3*f_ind+1] + dr[1];
			near_point[3*f_ind+2] = poc[3*f_ind+2] + dr[2];

#define CBL_ORIENTATION

#ifdef CBL_ORIENTATION 
			int np_f = f_ind_for_d[f_ref_CBL[domain][tpoc][0]][f_ref_CBL[domain][tpoc][1]];
			double xvec[3] = { poc[3*f_ind+0] - poc[3*np_f+0],
					   poc[3*f_ind+1] - poc[3*np_f+1],
					   poc[3*f_ind+2] - poc[3*np_f+2] };
#else

			double xvec[3] = { orientation[f_ind*3+0], orientation[f_ind*3+1], orientation[3*f_ind+2] };
#endif
			double dp = dot( xvec, npt );

			xvec[0] -= dp * npt[0];
			xvec[1] -= dp * npt[1];
			xvec[2] -= dp * npt[2];

			normalize(xvec);

			double yvec[3];
			cross( xvec, npt, yvec );

			normalize(yvec );

			double check[3];

			cross( xvec, yvec, check );

			dp = dot( check, dr );
	
			// dr is near point - poc (points toward surface away from protein )
			// if x cross y dotted into dr is positive, the orientation is wrong (y points "down" if x to the right is positive.)

			if( dp > 0 )
			{
				yvec[0] *= -1;
				yvec[1] *= -1;
				yvec[2] *= -1;
			
			} 

			dp = dot( dr, npt );

			if( dp > 0 )
			{
				// make normal point to the protein.
				//
				npt[0] *= -1;
				npt[1] *= -1;
				npt[1] *= -1;
			}

			normalize(xvec);
			normalize(yvec);

			tangent_plane[9*f_ind+0] = xvec[0];
			tangent_plane[9*f_ind+1] = xvec[1];
			tangent_plane[9*f_ind+2] = xvec[2];

			tangent_plane[9*f_ind+3] = yvec[0];
			tangent_plane[9*f_ind+4] = yvec[1];
			tangent_plane[9*f_ind+5] = yvec[2];

			tangent_plane[9*f_ind+6] = npt[0];
			tangent_plane[9*f_ind+7] = npt[1];
			tangent_plane[9*f_ind+8] = npt[2];

			// orientation is POC-Domain_COM, 
			// dr is          pt-poc  
			double dp_dist_check = dot( orientation + 3*f_ind, dr );
	
			if( dp_dist_check < 0 )
				dist[f_ind] = -the_dist; 
			else
				dist[f_ind] = the_dist; 

			f_ind++;
		}
	}
	
	double max_d = 30.0; // 2*max_d by 2*max_d grid
	double max_dx = 60.0; // 2*max_d by 2*max_d grid

	// Step 5: print the positions of nearby lipids in the tangent-plane coordinate system.

	for( int f1 = 0; f1 < nfaces; f1++ )
	{
		int d1 = domain_for_face[f1];
		int p1 = face_index[f1];

		for( int f2 = 0; f2 < nfaces; f2++ )
		{
			int d2 = domain_for_face[f2];
			int p2 = face_index[f2];
			
			double dr[3] = { near_point[f2*3+0] - near_point[3*f1+0],
					 near_point[f2*3+1] - near_point[3*f1+1],
					 near_point[f2*3+2] - near_point[3*f1+2] };
			theSimulation->wrapPBC(dr,alphas);	
			
			double x_tangent_plane = dot( dr, tangent_plane+9*f1 );
			double y_tangent_plane = dot( dr, tangent_plane+9*f1 + 3 );
			double z_tangent_plane = dot( dr, tangent_plane+9*f1 + 6 ); 

			if( z_tangent_plane < 0 )
			{
				// opposite leaflet..
				continue;
			}
					
			if( z_tangent_plane > 30 )
			{
				// probably on the other side of the protein?
				continue;
			}

			if( fabs(x_tangent_plane) > max_dx ) continue;	
			if( fabs(y_tangent_plane) > max_d ) continue;	


			double dist_to_CBL1, dist_to_CBL2;

			int f_CBL1 = f_ind_for_d[0][0]; 
			int f_CBL2 = f_ind_for_d[1][0]; 
			double dr_check[3] = {  near_point[3*f1+0] - near_point[3*f_CBL1+0],
					        near_point[3*f1+1] - near_point[3*f_CBL1+1],
						near_point[3*f1+2] - near_point[3*f_CBL1+2] };
			theSimulation->wrapPBC( dr_check, alphas );
			
			double r_check1 = normalize(dr_check);
			dr_check[0] = near_point[3*f1+0] - near_point[3*f_CBL2+0];
			dr_check[1] = near_point[3*f1+1] - near_point[3*f_CBL2+1];
			dr_check[2] = near_point[3*f1+2] - near_point[3*f_CBL2+2];
			theSimulation->wrapPBC( dr_check, alphas );
			double r_check2 = normalize(dr_check);

			printf("SYT-SYT %d domain %d POC %d dist %lf domain2 %d POC %d dist2 %lf x %lf y %lf z %lf NP: %lf %lf %lf rCBL1: %lf rCBL2: %lf\n", my_id, d1, p1, dist[f1], d2, p2, dist[f2], x_tangent_plane, y_tangent_plane, z_tangent_plane, near_point[3*f1+0], near_point[3*f1+1], near_point[3*f1+2], r_check1, r_check2  ); 
		}
	}


	for( int l = 0; l < ns; l++ )
	{
		int a = io_getNSAtom(l);
				
		// find near points for the attachment sites.
		double the_dist;
		int l_f;
		double col_u, col_v;
		double l_pt[3] = { lcom[3*l+0], lcom[3*l+1], lcom[3*l+2] };
		theSurface->nearPointOnBoxedSurface( l_pt, &l_f, &col_u, &col_v, M, mlow, mhigh, &the_dist );					
		double l_rpt[3], l_npt[3];
		theSurface->evaluateRNRM( l_f, col_u, col_v, l_rpt, l_npt, rsurf );
			
		double c_vec_1[3];
		double c_vec_2[3];
		double c_val1, c_val2;
	
		if( 1 ) // ALEX CONTINUE CODING HERE
		{
			double k;
			theSurface->c( l_f, col_u, col_v, rsurf, &k, c_vec_1, c_vec_2, &c_val1, &c_val2 ); 

			if( flip_sign )
			{
				c_val1 *= -1;
				c_val2 *= -1;
			}
	
			double dr_u[3];
			theSurface->ru( l_f, col_u, col_v, rsurf, dr_u );
			double drdu = normalize(dr_u);
			double dr_v[3];
			theSurface->rv( l_f, col_u, col_v, rsurf, dr_v );
			double drdv = normalize(dr_v);

			double tvec1[3]={0,0,0};

			tvec1[0] = c_vec_1[0] * dr_u[0] + c_vec_1[1] * dr_v[0];
			tvec1[0] = c_vec_1[0] * dr_u[1] + c_vec_1[1] * dr_v[1];
			tvec1[0] = c_vec_1[0] * dr_u[2] + c_vec_1[1] * dr_v[2];
			
			double tvec2[3]={0,0,0};

			tvec2[0] = c_vec_2[0] * dr_u[0] + c_vec_2[1] * dr_v[0];
			tvec2[0] = c_vec_2[0] * dr_u[1] + c_vec_2[1] * dr_v[1];
			tvec2[0] = c_vec_2[0] * dr_u[2] + c_vec_2[1] * dr_v[2];

			memcpy( c_vec_1, tvec1, sizeof(double)*3);
			memcpy( c_vec_2, tvec2, sizeof(double)*3);
		}

		for( int f = 0; f < nfaces; f++ )
		{
			int domain = domain_for_face[f];
			int p = face_index[f];
			// lipid minus near point.
			double dr_l[3] = { lcom[3*l+0] - near_point[3*f+0],
					   lcom[3*l+1] - near_point[3*f+1],
					   lcom[3*l+2] - near_point[3*f+2] };
			theSimulation->wrapPBC(dr_l,alphas);	

			double x_tangent_plane = dot( dr_l, tangent_plane+9*f );
			double y_tangent_plane = dot( dr_l, tangent_plane+9*f + 3 );
			double z_tangent_plane = dot( dr_l, tangent_plane+9*f + 6 ); 

			double local_c1_x = dot( c_vec_1, tangent_plane+9*f );
			double local_c1_y = dot( c_vec_1, tangent_plane+9*f + 3 );
			double local_c1_z = dot( c_vec_1, tangent_plane+9*f + 6 );
			
			double local_c2_x = dot( c_vec_2, tangent_plane+9*f );
			double local_c2_y = dot( c_vec_2, tangent_plane+9*f + 3 );
			double local_c2_z = dot( c_vec_2, tangent_plane+9*f + 6 );

			if( z_tangent_plane < 0 )
			{
				// opposite leaflet..
				continue;
			}
					
			if( z_tangent_plane > 30 )
			{
				// probably on the other side of the protein?
				continue;
			}

			if( fabs(x_tangent_plane) > max_dx ) continue;	
			if( fabs(y_tangent_plane) > max_d ) continue;	


			double dist_to_CBL1, dist_to_CBL2;

			int f_CBL1 = f_ind_for_d[0][0]; 
			int f_CBL2 = f_ind_for_d[1][0]; 
			double dr_check[3] = {  near_point[3*f+0] - near_point[3*f_CBL1+0],
					        near_point[3*f+1] - near_point[3*f_CBL1+1],
						near_point[3*f+2] - near_point[3*f_CBL1+2] };
			theSimulation->wrapPBC( dr_check, alphas );
			
			double r_check1 = normalize(dr_check);
			dr_check[0] = near_point[3*f+0] - near_point[3*f_CBL2+0];
			dr_check[1] = near_point[3*f+1] - near_point[3*f_CBL2+1];
			dr_check[2] = near_point[3*f+2] - near_point[3*f_CBL2+2];
			theSimulation->wrapPBC( dr_check, alphas );
			double r_check2 = normalize(dr_check);

			

			printf("SYT7 %d domain %d POC %d dist %lf lipid %s segid %s res %d x %lf y %lf z %lf NP: %lf %lf %lf rCBL1: %lf rCBL2: %lf c1: %lf %lf %lf val %lf c2: %lf %lf %lf val %lf \n", my_id, domain, p, dist[f], at[a].resname, at[a].segid, at[a].res, x_tangent_plane, y_tangent_plane, z_tangent_plane, near_point[3*f+0], near_point[3*f+1], near_point[3*f+2], r_check1, r_check2,
		local_c1_x, local_c1_y, local_c1_z, c_val1, local_c2_x, local_c2_y, local_c2_z, c_val2 ); 
				 
		}
	}

	// Last: clean up.
	
	free(lcom);
}




