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

#define SQRT_SAFETY (1e-7)

#define EPS_SMALL (1e-14)
#define WORKING

static double bond_k = 0;//0.1;
static double angle_k = 1;
static double dihe_k = 0;//10;
#define OLD_NRM

#ifdef OLD_NRM
static double k_nrm = 1.0;
#else
static double k_nrm = 1.0;
#endif
static double syt7_phi0 = 0; // M_PI/3; //M_PI;


static double dynamin_stride = 6.35;
static double dynamin_rotor  = 23.68;
static double dynamin_R = 95; // the approximate radius of the tube we want.
static double dynamin_angle = 110; // the angle that one dynamin dimer takes up.

static double monomer_MW = 200000; //amu per dimer.
static double C2_P_RADIUS = 25.0;
static double attach_p_radius = 5.0;
static double inter_bond = 50.0;
static double bond_length_1 = 90.0;
static double bond_length_2 = 70.0;
static double bond_angle_1  = 120.0; // degs of course
static double bond_angle_2  = 140.0; // degs of course
static double dynamin_dihe	    = 3.0; // very small.

void dynamin::init( double *r, int nmer )
{
	// aqueous initialization.

	nmer_saved =nmer;

	printf("Currently dynamin can only be initialized on the membrane.\n");
	exit(1);

	base_init();

	nsites = 5 * nmer;
	nattach = 2 * nmer;

	alloc();

	double vdw_r = 25.0;

	for( int s = 0; s < nsites; s++ )
	{
		sigma[s] = vdw_r;
		mass[s] = monomer_MW;
	}

	bound = 0;
	
	// initial geometry.

}

// initialize the BAR domain on the membrane.

void dynamin::init( Simulation *theSimulation, surface *theSurface, double *rsurf, int f, double u, double v, int nmer)
{
	// assume for now this is one of the points on the membrane neck.

	base_init();
	
	nmer_saved =nmer;

	nsites = 5 * nmer;
	nattach = 2 * nmer;

	alloc();

	double vdw_r = 25.0;

	for( int s = 0; s < 5 * nmer; s++ )
	{
		sigma[s] = vdw_r;
		mass[s] = monomer_MW;
	}

	for( int a = 0; a < nattach; a++ )
		sid[a] = theSurface->surface_id;

	bound = 1;

	// we get the attachment point of one connection, then place the other site by attempting to wrap around the positive curvature direction.


	int init = 1;

	double dynamin_axis[3] = { 0,0,0}; // the current best estimate of the dynamin axis.

	for( int m = 0; m < nmer; m++ )
	{
		double r_1[3];
		double r_2[3];

		double nrm_1[3];
		double nrm_2[3];

		double neg_dir_in_1[3];
		double neg_dir_in_2[3];

		if( init )
		{
			fs[0] = f;
			puv[0] = u;			
			puv[1] = v;
		}
		else
		{
			// pick the next spot	


			int p_f = fs[2*(m-1)+0];
			double p_u = puv[4*(m-1)+0],
				p_v = puv[4*(m-1)+1];

			double dr_u[3];
			theSurface->ru( p_f, p_u, p_v, rsurf, dr_u );
			double dr_v[3];
			theSurface->rv( p_f, p_u, p_v, rsurf, dr_v );
	
			// direction that aligns with the current dynamin axis.
			

			double duv1[2]; 
			best_align( duv1, dr_u, dr_v, dynamin_axis );

			duv1[0] *= dynamin_stride;
			duv1[1] *= dynamin_stride;

			int f_1 = f, nf = f;
			double uv1[2] = { u, v };
		
			do {
				f_1 = nf;
				nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf ); 
			} while( nf != f_1 );
			
			// now for the angular move.	

			double rotor_stride = (2*M_PI*dynamin_R) * (dynamin_rotor / 360.0);

			double r_pt[3], npt[3];

			theSurface->evaluateRNRM(f_1, uv1[0], uv1[1], r_pt, npt, rsurf );

			double perp_dir[3];
			cross( npt, dynamin_axis, perp_dir );

			best_align( duv1, dr_u, dr_v, dynamin_axis );

			duv1[0] *= rotor_stride;
			duv1[1] *= rotor_stride;

			do {
				f_1 = nf;
				nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf ); 
			} while( nf != f_1 );
			  
			fs[2*m] = f_1;

			puv[4*m+0] = uv1[0];
			puv[4*m+1] = uv1[1];
		}
		
		theSurface->evaluateRNRM( fs[2*m+0], puv[4*m+0], puv[4*m+1], r_1, nrm_1, rsurf );	
			
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

		if( c_val1 > c_val2 )
		{
			memcpy( positive_dir_uv, c_vec_1, sizeof(double) * 2 );
			memcpy( negative_dir_uv, c_vec_2, sizeof(double) * 2 );
		}
		else
		{
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
		
		normalize(d_axis_expec_1);
			
		if( init )
		{
			dynamin_axis[0] = d_axis_expec_1[0];
			dynamin_axis[1] = d_axis_expec_1[1];
			dynamin_axis[2] = d_axis_expec_1[2];
		}
		double d_move[2] = { positive_dir_uv[0], positive_dir_uv[1] };

		double d_move_expec[3] = {0,0,0};
		
		d_move_expec[0] += d_move[0] * dr_u[0]*drdu;
		d_move_expec[1] += d_move[0] * dr_u[1]*drdu;
		d_move_expec[2] += d_move[0] * dr_u[2]*drdu;
		
		d_move_expec[0] += d_move[1] * dr_v[0]*drdv;
		d_move_expec[1] += d_move[1] * dr_v[1]*drdv;
		d_move_expec[2] += d_move[1] * dr_v[2]*drdv;
	
		// handedness rule.
	
		double nrm_expec[3];
		cross( dynamin_axis, d_move_expec, nrm_expec );

		double dp = dot( nrm_expec, nrm_1 );
	
		if( dp < 0 )
		{
			d_move[0] *= -1;
			d_move[1] *= -1;
		}
	
		double len = normalize(d_move_expec);
		double r_move = 2 * M_PI * dynamin_R * (dynamin_angle/360.0);
	
		d_move[0] *= r_move/len;
		d_move[1] *= r_move/len;
	
		int f_1 = f, nf = f;
		double uv1[2] = { u, v };
		double duv1[2] = { d_move[0], d_move[1] }; 
		
		do {
			f_1 = nf;
			nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf ); 
		
			double rp_t[3];
			double nrm_t[3];
			theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], rp_t, nrm_t, rsurf);
		} while( nf != f_1 );

		fs[2*m+1] = f_1;
		puv[4*m+2] = uv1[0];
		puv[4*m+3] = uv1[1];	
		
		theSurface->evaluateRNRM( fs[2*m+1], puv[4*m+2], puv[4*m+3], r_2, nrm_2, rsurf );	

		// this is the best estimate of the angular (+curvature) Cartesian direction according to the 2nd site.
		double use_pos_dir_2[3];

		double d_axis_expec_2[3] = {0,0,0};
		{
			double dr_u[3];
			theSurface->ru( f_1, uv1[0], uv1[1], rsurf, dr_u );
			double drdu = normalize(dr_u);
			double dr_v[3];
			theSurface->rv( f_1, uv1[0], uv1[1], rsurf, dr_v );
			double drdv = normalize(dr_v);
			
			// evaluate the curvature at the other attachment site to get a better (?) estimate of the dynamin axis.
	
			double c_vec_1[3];
			double c_vec_2[3];
			double c_val1, c_val2;
			double k;
			theSurface->c( f_1, uv1[0], uv1[1], rsurf, &k, c_vec_1, c_vec_2, &c_val1, &c_val2 ); 
		
			double negative_dir_uv[2];

			if( c_val1 > c_val2 )
				memcpy( negative_dir_uv, c_vec_2, sizeof(double) * 2 );
			else
				memcpy( negative_dir_uv, c_vec_1, sizeof(double) * 2 );
		
			double d_axis[2] = { negative_dir_uv[0], negative_dir_uv[1] };

			
			d_axis_expec_2[0] += d_axis[0] * dr_u[0]*drdu;
			d_axis_expec_2[1] += d_axis[0] * dr_u[1]*drdu;
			d_axis_expec_2[2] += d_axis[0] * dr_u[2]*drdu;
			
			d_axis_expec_2[0] += d_axis[1] * dr_v[0]*drdv;
			d_axis_expec_2[1] += d_axis[1] * dr_v[1]*drdv;
			d_axis_expec_2[2] += d_axis[1] * dr_v[2]*drdv;
	
			// handedness.	
			double nrm_expec[3];
			cross( dynamin_axis, d_move_expec, nrm_expec );

			double dp = dot( nrm_expec, nrm_2 );

			if( dp < 0 )
			{
				d_axis[0] *= -1;
				d_axis[1] *= -1;
				d_axis[2] *= -1;
			}

			normalize(d_axis_expec_2);

			cross( dynamin_axis, nrm_2, use_pos_dir_2 ); // shoot this back toward the first unit.

		 	normalize(use_pos_dir_2 );
		}

		// get the positions of the aqueous sites. put them in a plane perpendicular to the dynamin axis.

		double rotor_1 = 45.0;

		rall[3*nattach+m*3*3+0*3+0] = r_1[0] + nrm_1[0] * bond_length_1;	
		rall[3*nattach+m*3*3+0*3+1] = r_1[1] + nrm_1[1] * bond_length_1;	
		rall[3*nattach+m*3*3+0*3+2] = r_1[2] + nrm_1[2] * bond_length_1;	

		rotateArbitrary( rall + 3*nattach+m*3*3+0*3, dynamin_axis, r_1, 1, -rotor_1*(M_PI/180.0) ); 
		
		rall[3*nattach+m*3*3+2*3+0] = r_2[0] + nrm_2[0] * bond_length_1;
		rall[3*nattach+m*3*3+2*3+1] = r_2[1] + nrm_2[1] * bond_length_1;
		rall[3*nattach+m*3*3+2*3+2] = r_2[2] + nrm_2[2] * bond_length_1;	
		
		rotateArbitrary( rall + 3*nattach+m*3*3+2*3, dynamin_axis, r_2, 1, +rotor_1*(M_PI/180.0) ); 
	
		double shift_up = 20.0;
		double av_nrm[3] = { nrm_1[0]+nrm_2[0], nrm_1[1]+nrm_2[1], nrm_1[2]+nrm_2[2] };

		normalize(av_nrm);

		rall[3*nattach+m*3*3+1*3+0] = (rall[3*nattach+m*3*3+0*3+0]+rall[3*nattach+m*3*3+2*3+0])/2 + av_nrm[0] * shift_up;	
		rall[3*nattach+m*3*3+1*3+1] = (rall[3*nattach+m*3*3+0*3+1]+rall[3*nattach+m*3*3+2*3+1])/2 + av_nrm[1] * shift_up; 	
		rall[3*nattach+m*3*3+1*3+2] = (rall[3*nattach+m*3*3+0*3+2]+rall[3*nattach+m*3*3+2*3+2])/2 + av_nrm[2] * shift_up;	

		double dp1=1, dp2;
		dp1 = dot( d_axis_expec_1, dynamin_axis );
		dp2 = dot( d_axis_expec_2, dynamin_axis );
		
		if( dp1 < 0 )
		{
			d_axis_expec_1[0] *= -1;
			d_axis_expec_1[1] *= -1;
			d_axis_expec_1[2] *= -1;
		}

		if( dp2 < 0 )
		{
			d_axis_expec_2[0] *= -1;
			d_axis_expec_2[1] *= -1;
			d_axis_expec_2[2] *= -1;
		}
	
		// update the dynamin axis	

		double w = 0.75;

		dynamin_axis[0] = w * dynamin_axis[0] + (1-w) * (d_axis_expec_1[0] + d_axis_expec_2[0])/2;
		dynamin_axis[1] = w * dynamin_axis[1] + (1-w) * (d_axis_expec_1[1] + d_axis_expec_2[1])/2;
		dynamin_axis[2] = w * dynamin_axis[2] + (1-w) * (d_axis_expec_1[2] + d_axis_expec_2[2])/2;

		normalize( dynamin_axis );

		init = 0;
	}

	memcpy( grad_fs, fs, sizeof(int) * nattach );
	memcpy( grad_puv, puv, sizeof(double) * 2*nattach );

	memset( PBC_ext, 0, sizeof(double) * 3*nsites );
	
	setrall(theSimulation);
	
	// wrap them all to one atom

	for( int a = 1; a < nsites; a++ )
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


int dynamin::getNBonds( void )
{
	return nmer_saved * 4 + (nmer_saved-1) * 2;
}

void dynamin::putBonds( int *bond_list, double *bond_r, double *bond_k_in)
{
	int bonds_per_nmer = 4 + 2;
	int nsites_per_nmer = 5;

	int attach_per_nmer = 2;
	int aq_per_nmer = 3;

	for( int m = 0; m < nmer_saved; m++ )
	{
		// internal bonds.
		int a0 = nattach + m * aq_per_nmer;
		int s0 = attach_per_nmer * m;
		
		bond_list[(m*bonds_per_nmer)*2+0] = s0;
		bond_list[(m*bonds_per_nmer)*2+1] = a0;
		if( bond_r ) bond_r[m*bonds_per_nmer] = bond_length_1;
		if( bond_k_in ) bond_k_in[m*bonds_per_nmer] = bond_k;
		
		bond_list[(m*bonds_per_nmer)*2+2] = a0;
		bond_list[(m*bonds_per_nmer)*2+3] = a0+1;
		if( bond_r ) bond_r[m*bonds_per_nmer+1] = bond_length_2;
		if( bond_k_in ) bond_k_in[m*bonds_per_nmer+1] = bond_k;
		
		bond_list[(m*bonds_per_nmer)*2+4] = a0+1;
		bond_list[(m*bonds_per_nmer)*2+5] = a0+2;
		if( bond_r ) bond_r[m*bonds_per_nmer+2] = bond_length_2;
		if( bond_k_in ) bond_k_in[m*bonds_per_nmer+2] = bond_k;
		
		bond_list[(m*bonds_per_nmer)*2+6] = a0+2;
		bond_list[(m*bonds_per_nmer)*2+7] = s0+1;
		if( bond_r ) bond_r[m*bonds_per_nmer+3] = bond_length_1;
		if( bond_k_in ) bond_k_in[m*bonds_per_nmer+3] = bond_k;

		// inter bonds.
		
		if( m < nmer_saved-1 )
		{
			int a0 = nattach + m * aq_per_nmer;
			int a1 = nattach + (m+1) * aq_per_nmer;

			bond_list[(m*bonds_per_nmer)*2+8] = a0+0;
			bond_list[(m*bonds_per_nmer)*2+9] = a1+1;
			if( bond_r ) bond_r[m*bonds_per_nmer+4] = inter_bond;
			if( bond_k_in ) bond_k_in[m*bonds_per_nmer+4] = bond_k;
			
			bond_list[(m*bonds_per_nmer)*2+10] = a0+1;
			bond_list[(m*bonds_per_nmer)*2+11] = a1+2;
			if( bond_r ) bond_r[m*bonds_per_nmer+4] = inter_bond;
			if( bond_k_in ) bond_k_in[m*bonds_per_nmer+5] = bond_k;
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

	int aq_per_nmer = 3;

	double pot = 0;

	for( int b = 0; b < nbonds; b++ )
	{
		int a1 = theBonds[2*b+0];
		int a2 = theBonds[2*b+1];

		double dr[3] = { rall[3*a1+0] - rall[3*a2+0], rall[3*a1+1] - rall[3*a2+1], rall[3*a1+2] - rall[3*a2+2] };
		
		double r0 = bond_r_p[b];
		double k  = bond_k_p[b];

		double ln = normalize(dr);

		pot += 0.5 * k * (ln - r0) * (ln-r0);	
		
	}

	// angles

	int nattach_per_nmer = 2;

	for( int m = 0; m < nmer_saved; m++ )
	{
		int s0 = m*nattach_per_nmer;
		int a0 = nattach + m*aq_per_nmer;

		int angles[3][3] = { {s0,a0,a0+1}, {a0,a0+1,a0+2}, {a0+1,a0+2,s0+1} };
		double phi0[3] = { (M_PI/180.0)*bond_angle_1, (M_PI/180.0)*bond_angle_2, (M_PI/180.0)*bond_angle_1 };
		 
		for( int a = 0; a < 3; a++ )
		{
			int a1 = angles[a][0];
			int a2 = angles[a][1];
			int a3 = angles[a][2];
			
			double l_angle_k = angle_k;

			double dr1[3] = { r[3*a1+0] - r[3*a2+0],
					 r[3*a1+1] - r[3*a2+1],
					 r[3*a1+2] - r[3*a2+2] };
			double dr2[3] = { r[3*a3+0] - r[3*a2+0],
					 r[3*a3+1] - r[3*a2+1],
					 r[3*a3+2] - r[3*a2+2] };
			double ln1 = normalize(dr1);
			double ln2 = normalize(dr2);

			double dp = dr1[0]*dr2[0]+dr1[1]*dr2[1]+dr1[2]*dr2[2];
			double ac = acos(dp);

			pot += 0.5 * l_angle_k * (ac - phi0[a]) * (ac -phi0[a]);	
		}
	}

	for( int m = 0; m < nmer_saved; m++ )
	{
		int s0 = m*nattach_per_nmer;
		int a0 = nattach + m*aq_per_nmer;

		int dihedrals[2][4] = { {s0,a0,a0+1,a0+2}, {a0,a0+1,a0+2,s0+1} };
		double phi0[2] = { (M_PI/180.0)*dynamin_dihe, (M_PI/180.0)*dynamin_dihe };
		 
		for( int d = 0; d < 2; d++ )
		{
			double l_dihe_k = dihe_k;
			double dihedral_target = phi0[d];

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
	
			pot += 0.5 * l_dihe_k * (theta - dihedral_target ) * (theta - dihedral_target );
		}
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
	
	int aq_per_nmer = 3;


	for( int b = 0; b < nbonds; b++ )
	{
		int a1 = theBonds[2*b+0];
		int a2 = theBonds[2*b+1];

		double dr[3] = { rall[3*a1+0] - rall[3*a2+0], rall[3*a1+1] - rall[3*a2+1], rall[3*a1+2] - rall[3*a2+2] };
		
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

	// angles

	int nattach_per_nmer = 2;

	for( int m = 0; m < nmer_saved; m++ )
	{
		int s0 = m*nattach_per_nmer;
		int a0 = nattach + m*aq_per_nmer;

		int angles[3][3] = { {s0,a0,a0+1}, {a0,a0+1,a0+2}, {a0+1,a0+2,s0+1} };
		double phi0[3] = { (M_PI/180.0)*bond_angle_1, (M_PI/180.0)*bond_angle_2, (M_PI/180.0)*bond_angle_1 };
		 
		for( int a = 0; a < 3; a++ )
		{
			int a1 = angles[a][0];
			int a2 = angles[a][1];
			int a3 = angles[a][2];
			double l_angle_k = angle_k;

			double dr1[3] = { r[3*a1+0] - r[3*a2+0],
					 r[3*a1+1] - r[3*a2+1],
					 r[3*a1+2] - r[3*a2+2] };
			double dr2[3] = { r[3*a3+0] - r[3*a2+0],
					 r[3*a3+1] - r[3*a2+1],
					 r[3*a3+2] - r[3*a2+2] };
			double ln1 = normalize(dr1);
			double ln2 = normalize(dr2);

			double dp = dr1[0]*dr2[0]+dr1[1]*dr2[1]+dr1[2]*dr2[2];
			double ac = acos(dp);

			pot += 0.5 * l_angle_k * (ac - phi0[a]) * (ac -phi0[a]);	
	
			double der_A = l_angle_k * (ac-phi0[a]);
			double der_cos = -1.0/sqrt(1.0-dp*dp); 
	
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
	
	}

	for( int m = 0; m < nmer_saved; m++ )
	{
		int s0 = m*nattach_per_nmer;
		int a0 = nattach + m*aq_per_nmer;

		int dihedrals[2][4] = { {s0,a0,a0+1,a0+2}, {a0,a0+1,a0+2,s0+1} };
		double phi0[2] = { (M_PI/180.0)*dynamin_dihe, (M_PI/180.0)*dynamin_dihe };
		 
		for( int d = 0; d < 2; d++ )
		{
			double l_dihe_k = dihe_k;
			double dihedral_target = phi0[d];

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
	
			pot += 0.5 * l_dihe_k * (theta - dihedral_target ) * (theta - dihedral_target );


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
			double base_der = l_dihe_k * (theta - dihedral_target) * (-1 / sqrt(1-dp*dp));
		
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


void dynamin::writeStructure( Simulation *theSimulation, surface_mask *upperSurfaceMask, surface_mask *lowerSurfaceMask, struct atom_rec **at_out, int *nat_out, char **sequence, struct ion_add **ions, int *nions, 
	aa_build_data *buildData )
{
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
	*sequence = seq;
}

/* this function 
 * "gets" the coarse-grained coordinates of the observed protein from the molecular structure -- and its fit to the membrane surface.
 *
 * */


void dynamin::get( 
		Simulation *theSimulation, // this is the surface/PBC
		struct atom_rec *at,
		int syt7_start, int syt7_stop )
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
	double rvecB[3] = { rall[5*3+0] - p[6],
			    rall[5*3+1] - p[7],
			    rall[5*3+2] - p[8] };
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
	switch (p%5)
	{
		case 2:
		case 3:
		case 4:
			return 'N';

		default:
			return 'O';
	}
}



