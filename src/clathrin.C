#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "interp.h"
#include "input.h"
#include "simulation.h"
#include "mutil.h"
#include "dcd.h"
#include "alignSet.h"
#include <math.h>

static double *r_transform = NULL;
int ntransform = 0;
int do_rotor = 0;

void Simulation::setupClathrinFitter( struct parameterBlock *block )
{
	if( ! allSurfaces ) 
		return;

	if( !block->clathrinStructure )
	{
		printf("Clathrin fitting requires a clathrinStructure specified.\n");
		exit(1);
	}
	if( block->clathrin_rotor )
		do_rotor = 1;


	surface *mainSurface = allSurfaces->theSurface;

	// given a planar or spherical membrane, adds the potential for use by minimize.C	

	FILE *clathrinFile = fopen( block->clathrinStructure, "r" );

	if( !clathrinFile )
	{
		printf("Couldn't open clathrin file '%s'.\n", block->clathrinStructure );
		exit(1);
	}

	loadPSFfromPDB( clathrinFile );

	int nat = curNAtoms();

	struct atom_rec *at = (struct atom_rec *)malloc( sizeof(struct atom_rec ) * nat );

	rewind(clathrinFile);

	loadPDB(clathrinFile, at);
	
	for( int a = 0; a < nat; a++ )
	{
		at[a].x *= 10;
		at[a].y *= 10;
		at[a].z *= 10;
	}
	double com[3];
	double td[3][3];

	int nspace = 10;
	double *rpts = (double *)malloc( sizeof(double) * 3 * nspace );
	double *nrms = (double *)malloc( sizeof(double) * 3 * nspace );
	int npts = 0;

	double *oset = (double *)malloc( sizeof(double) * 3 * nat );

	for( int a = 0; a < nat; a++ )
	{
		oset[3*a+0] = at[a].x;
		oset[3*a+1] = at[a].y;
		oset[3*a+2] = at[a].z;
	
		if( !strcasecmp( at[a].atname, "COM") ) 
		{
			com[0] = at[a].x;	
			com[1] = at[a].y;	
			com[2] = at[a].z;	
		}
		else if( !strncasecmp( at[a].atname, "td", 2) && strlen(at[a].atname) > 2 )
		{
			int ind = at[a].atname[2] - '1';

			if( ind >= 0 && ind < 3 )
			{
				td[ind][0] = at[a].x;
				td[ind][1] = at[a].y;
				td[ind][2] = at[a].z;
			}
		}

		if( !strcasecmp( at[a].atname, "td3" ) )
		{


			double avp[3] = {0,0,0};
			for( int t =0; t < 3; t++ )
			{
				avp[0] += td[t][0]/3;
				avp[1] += td[t][1]/3;
				avp[2] += td[t][2]/3;
			}
			
			double dr[3] = { avp[0] - com[0], avp[1] - com[1], avp[2] - com[2] };
			double l = normalize(dr);

			double h = l + block->clathrin_h;
			
			if( npts == nspace )
			{
				nspace *= 2;
				rpts = (double *)realloc( rpts, sizeof(double)*3*nspace );
				nrms = (double *)realloc( nrms, sizeof(double)*3*nspace );
			}

			rpts[3*npts+0] =  com[0];
			rpts[3*npts+1] =  com[1];
			rpts[3*npts+2] =  com[2];

			nrms[3*npts+0] = dr[0];
			nrms[3*npts+1] = dr[1];
			nrms[3*npts+2] = dr[2];
			
			npts++;
				
		}
	}

	// find best fit center.
	
	double avt = 0;
	double rc[3]={0,0,0};

	for( int c = 0; c < 3; c++ )
	{
		double nrmrsum = 0;
		double nrmsum = 0;
		double rsum = 0;
		double num = npts;
		double nrm2 = 0;

		for( int p = 0; p < npts; p++ )
		{
			nrmrsum += nrms[3*p+c] * rpts[3*p+c];
			nrmsum += nrms[3*p+c];
			rsum += rpts[3*p+c];
			nrm2 += nrms[3*p+c]*nrms[3*p+c];			
		}

		double t = - (-nrmrsum * num + nrmsum * rsum)/(nrmsum*nrmsum-nrm2*num);
		double rcx = -(nrmrsum*nrmsum-nrm2*rsum)/(-nrmsum*nrmsum+nrm2*num);

		rc[c] = rcx;
		avt += t/3;
		printf("c: %d t: %lf rc: %lf\n", c, t, rcx );
	}
	
	printf("Center: %lf %lf %lf\n", rc[0], rc[1], rc[2] );		

	
	double av_norm[3] = { 0,0,0 };

	for( int p = 0; p < npts; p++ )
	{
		av_norm[0] += nrms[3*p+0];
		av_norm[1] += nrms[3*p+1];
		av_norm[2] += nrms[3*p+2];
	}

	normalize(av_norm);

	double crossp[3];
	double z_axis[3]={0,0,1};
	double ortho[3];
	cross( av_norm, z_axis, ortho );
	normalize(ortho);

	double h = block->clathrin_h;

	for( int p = 0; p < npts; p++ )
	{
		if(  block->recenter_clathrin )
		{
			rpts[3*p+0] -= rc[0];
			rpts[3*p+1] -= rc[1];
			rpts[3*p+2] -= rc[2];
		}
		rpts[3*p+0] += nrms[3*p+0] * h; 
		rpts[3*p+1] += nrms[3*p+1] * h; 
		rpts[3*p+2] += nrms[3*p+2] * h; 
	}

	if(  block->recenter_clathrin )
	{
		for( int p = 0; p < nat; p++ )
		{
			oset[3*p+0] -= rc[0];
			oset[3*p+1] -= rc[1];
			oset[3*p+2] -= rc[2];
		}
	}

	double origin[3]={0,0,0};

	if(  block->recenter_clathrin )
	{
		rotateArbitrary( rpts, ortho, origin, npts, acos(av_norm[2]) ); 
		rotateArbitrary( oset, ortho, origin, nat, acos(av_norm[2]) ); 
	}
	FILE *outputCheck = fopen("clathrin_check.xyz","w");

	for( int i = 0; i < nat; i++ )
	{
		at[i].x = oset[3*i+0];
		at[i].y = oset[3*i+1];
		at[i].z = oset[3*i+2];

		fprintf(outputCheck, "C %lf %lf %lf\n", oset[3*i+0], oset[3*i+1], oset[3*i+2] );
		//printATOM( outputCheck, 1+i, 1+i, at+i );
	}

	r_transform = (double *)malloc( sizeof(double) * 3 * nat );
	memcpy( r_transform, oset, sizeof(double) * 3 * nat ); 
	ntransform = nat;
	fclose(outputCheck);

	printf("rtransform0: %lf %lf %lf\n", r_transform[0], r_transform[1], r_transform[2] );

	for( int p = 0; p < npts; p++ )
	{	
		mainSurface->addFixedPoint( rpts+3*p, 0, block->clathrin_force_k );
/*		fixed_cut_point *new_point = (fixed_cut_point*)malloc( sizeof(fixed_cut_point) );
		new_point->rpt[0] = rpts[3*p+0]; 
		new_point->rpt[1] = rpts[3*p+1]; 
		new_point->rpt[2] = rpts[3*p+2]; 
		
		new_point->frc[0] = 0; 
		new_point->frc[1] = 0; 
		new_point->frc[2] = 0; 

		printf("rp: %lf %lf %lf\n", rpts[3*p+0], rpts[3*p+1], rpts[3*p+2] );
		new_point->k = block->clathrin_force_k;
		new_point->

		new_point->next = mainSurface->cutPoints;
		mainSurface->cutPoints = new_point; 
*/
	}
}

static double *rpos_0;
static double rcom_rot[3];

double rotate_pts( Simulation *theSimulation, double *eval_at, double dcom[3])
{
	surface *mainSurface = theSimulation->allSurfaces->theSurface;
	double *main_surface_r = theSimulation->allSurfaces->r;	
	double frot = 0;

	double rotor[3] = { eval_at[0], eval_at[1], eval_at[2] };

	int ind = 0;
	for( fixed_cut_point *thePt = mainSurface->cutPoints; thePt; thePt = thePt->next, ind++)
	{
		double dr[3] = { rpos_0[3*ind+0] - rcom_rot[0],
				 rpos_0[3*ind+1] - rcom_rot[1],
				 rpos_0[3*ind+2] - rcom_rot[2] };

		double dr2[3];

		dr2[0] = dr[0];
		dr2[1] = cos(rotor[0]) * dr[1] + sin(rotor[0]) * dr[2];
		dr2[2] = -sin(rotor[0]) * dr[1] + cos(rotor[0]) * dr[2];
		
		dr[0] = cos(rotor[1]) * dr2[0] + sin(rotor[1]) * dr2[2];
		dr[1] = dr2[1];
		dr[2] = -sin(rotor[1]) * dr2[0] + cos(rotor[1]) * dr2[2];
		
		dr2[0] = cos(rotor[2]) * dr[0] + sin(rotor[2]) * dr[1];
		dr2[1] = -sin(rotor[2]) * dr[0] + cos(rotor[2]) * dr[1];
		dr2[2] = dr[2];

		thePt->rpt[0] = rcom_rot[0] + dr2[0] + dcom[0];
		thePt->rpt[1] = rcom_rot[1] + dr2[1] + dcom[1];
		thePt->rpt[2] = rcom_rot[2] + dr2[2] + dcom[2];
		
	}

}

void set_saved_transform( Simulation *theSimulation, double rotor[3], double dcom[3] )
{
	for( int t = 0; t < ntransform; t++ )
	{
		double dr[3] = { r_transform[3*t+0] - rcom_rot[0],
				 r_transform[3*t+1] - rcom_rot[1],
				 r_transform[3*t+2] - rcom_rot[2] };
		
		double dr2[3];
	
		dr2[0] = dr[0];
		dr2[1] = cos(rotor[0]) * dr[1] + sin(rotor[0]) * dr[2];
		dr2[2] = -sin(rotor[0]) * dr[1] + cos(rotor[0]) * dr[2];
		
		dr[0] = cos(rotor[1]) * dr2[0] + sin(rotor[1]) * dr2[2];
		dr[1] = dr2[1];
		dr[2] = -sin(rotor[1]) * dr2[0] + cos(rotor[1]) * dr2[2];
		
		dr2[0] = cos(rotor[2]) * dr[0] + sin(rotor[2]) * dr[1];
		dr2[1] = -sin(rotor[2]) * dr[0] + cos(rotor[2]) * dr[1];
		dr2[2] = dr[2];
	
		r_transform[3*t+0] = rcom_rot[0] + dr2[0] + dcom[0];
		r_transform[3*t+1] = rcom_rot[1] + dr2[1] + dcom[1];
		r_transform[3*t+2] = rcom_rot[2] + dr2[2] + dcom[2];
	}
	

	
	
	
}

double rot_fenergy( Simulation *theSimulation, double rotor_g[6], double lambda, int apply_transform )
{
	surface *mainSurface = theSimulation->allSurfaces->theSurface;
	double *main_surface_r = theSimulation->allSurfaces->r;	
	double frot = 0;

	double rotor[3] = { -lambda * rotor_g[0], -lambda * rotor_g[1], -lambda * rotor_g[2] };
	double dcom[3] = { -lambda * rotor_g[3], -lambda * rotor_g[4], - lambda * rotor_g[5] };
	rotate_pts( theSimulation, rotor, dcom );


	if( apply_transform )
	{
		for( int t = 0; t < ntransform; t++ )
		{
			double dr[3] = { r_transform[3*t+0] - rcom_rot[0],
					 r_transform[3*t+1] - rcom_rot[1],
					 r_transform[3*t+2] - rcom_rot[2] };
			
			double dr2[3];
	
			dr2[0] = dr[0];
			dr2[1] = cos(rotor[0]) * dr[1] + sin(rotor[0]) * dr[2];
			dr2[2] = -sin(rotor[0]) * dr[1] + cos(rotor[0]) * dr[2];
			
			dr[0] = cos(rotor[1]) * dr2[0] + sin(rotor[1]) * dr2[2];
			dr[1] = dr2[1];
			dr[2] = -sin(rotor[1]) * dr2[0] + cos(rotor[1]) * dr2[2];
			
			dr2[0] = cos(rotor[2]) * dr[0] + sin(rotor[2]) * dr[1];
			dr2[1] = -sin(rotor[2]) * dr[0] + cos(rotor[2]) * dr[1];
			dr2[2] = dr[2];
	
			r_transform[3*t+0] = dcom[0] + rcom_rot[0] + dr2[0];
			r_transform[3*t+1] = dcom[1] + rcom_rot[1] + dr2[1];
			r_transform[3*t+2] = dcom[2] + rcom_rot[2] + dr2[2];
		}
	

		return 0;	
	}
	else
		return mainSurface->cutEnergy(main_surface_r);
}

void Simulation::clathrinGrad( double *g ) 
// g points to the x/y/z rotation variables
// 
{
	surface *mainSurface = allSurfaces->theSurface;

	// optimize the rotation of the clathrin surface with the membrane fixed.	

	g[0] = 0;
	g[1] = 0;
	g[2] = 0;

	// rho is the rotation around the x axis.
	// phi is the rotation around the y axis.
	// theta is the rotation around the z axis.

	double com[3] = {0,0,0};
	int npts = 0;
	for( fixed_cut_point *thePt = mainSurface->cutPoints; thePt; thePt = thePt->next )
	{
		com[0] += thePt->rpt[0];	
		com[1] += thePt->rpt[1];	
		com[2] += thePt->rpt[2];	

		npts += 1;
	}

	com[0] /= npts;
	com[1] /= npts;
	com[2] /= npts;

	rcom_rot[0] = com[0];
	rcom_rot[1] = com[1];
	rcom_rot[2] = com[2];

	rpos_0 = (double *)malloc( sizeof(double) * 3 * npts );
	
	int ind = 0;
	for( fixed_cut_point *thePt = mainSurface->cutPoints; thePt; thePt = thePt->next, ind++ )
	{
		double dr_d_angle[3] ={0,0,0};
		double rpt[3];	

		rpt[0] = thePt->rpt[0];
		rpt[1] = thePt->rpt[1];
		rpt[2] = thePt->rpt[2];

		// rotation about x axis.
		
		dr_d_angle[0] = 0;
		dr_d_angle[1] = (rpt[2] - com[2]);
		dr_d_angle[2] = -(rpt[1] - com[1]);

		g[0] += dr_d_angle[0] * thePt->frc[0];
		g[0] += dr_d_angle[1] * thePt->frc[1];
		g[0] += dr_d_angle[2] * thePt->frc[2];
		
		// rotation about y axis.
		
		dr_d_angle[0] = rpt[2] - com[2];
		dr_d_angle[1] = 0;
		dr_d_angle[2] = -(rpt[0] - com[0]);

		g[1] += dr_d_angle[0] * thePt->frc[0];
		g[1] += dr_d_angle[1] * thePt->frc[1];
		g[1] += dr_d_angle[2] * thePt->frc[2];
		
		// rotation about z axis.
		
		dr_d_angle[0] = rpt[1] - com[1];
		dr_d_angle[1] = -(rpt[0]-com[0]);
		dr_d_angle[2] = 0;

		g[2] += dr_d_angle[0] * thePt->frc[0];
		g[2] += dr_d_angle[1] * thePt->frc[1];
		g[2] += dr_d_angle[2] * thePt->frc[2];

		rpos_0[ind*3+0] = rpt[0];
		rpos_0[ind*3+1] = rpt[1];
		rpos_0[ind*3+2] = rpt[2];

		g[3] += thePt->frc[0]; 
		g[4] += thePt->frc[1]; 
		g[5] += thePt->frc[2]; 
	}

	if( !do_rotor )
	{
		g[0] = 0;
		g[1] = 0;
		g[2] = 0;
	}
}


void Simulation::OptClathrin( void )
{
	surface *mainSurface = allSurfaces->theSurface;

	// optimize the rotation of the clathrin surface with the membrane fixed.	

	double g[6] = {0,0,0,0,0,0}; // phi, theta, rho 
	// rho is the rotation around the x axis.
	// phi is the rotation around the y axis.
	// theta is the rotation around the z axis.
	clathrinGrad(g);

	// do one-dimensional bisection
	
	int opt_energy_mode = 0, write_transform = 1;

	double lambda = 0;
	double e0 = rot_fenergy(this, g, lambda, opt_energy_mode); 

	double lambda_low = 0;
	double elow = e0;
	double lambda_high = 0;
	double ehigh = e0;

	double lambda_max = 1e-10;

	printf("e0: %.14le\n", e0 );

	int done = 0;

	while( ! done )
	{
		double e1 = rot_fenergy( this, g, lambda_max, opt_energy_mode );
		printf("lambda: %.14le dE: %.14le\n", lambda_max, e1-e0 );

		if( e1 < elow )
		{
			elow = e1;
			lambda_low = lambda_max;
		}
		else
		{
			ehigh = e1;
			lambda_high = lambda_max;
			break;
		}

		lambda_max *= 1.25;
	}

	rot_fenergy( this, g, lambda_low, write_transform );

	free(rpos_0);
}

void Simulation::WriteClathrin(parameterBlock *block )
{
	char fileName[256];

	sprintf(fileName, "%s_clathrin.pdb", block->jobName );
	
	FILE *theFile = fopen(fileName,"w");


	surface *mainSurface = allSurfaces->theSurface;

	int npts = 0;
	for( fixed_cut_point *thePt = mainSurface->cutPoints; thePt; thePt = thePt->next )
		npts += 1;
	
	
	FILE *clathrinFile = fopen( block->clathrinStructure, "r" );

	if( !clathrinFile )
	{
		printf("Couldn't open clathrin file '%s'.\n", block->clathrinStructure );
		exit(1);
	}

	loadPSFfromPDB( clathrinFile );

	int nat = curNAtoms();

	struct atom_rec *at = (struct atom_rec *)malloc( sizeof(struct atom_rec ) * nat );

	rewind(clathrinFile);
	loadPDB(clathrinFile, at);

	fclose(clathrinFile);
	int ind = 0;
	
	for( int a = 0; a < nat; a++ )
	{
		at[a].x = r_transform[3*a+0]/10;
		at[a].y = r_transform[3*a+1]/10;
		at[a].z = r_transform[3*a+2]/10;

		printATOM( theFile, 1+a, at[a].res, at+a );
	}	

	fclose(theFile);

	sprintf(fileName, "%s_clathrin_A.pdb", block->jobName );
	
	theFile = fopen(fileName,"w");
	
	for( int a = 0; a < nat; a++ )
	{
		at[a].x = r_transform[3*a+0];
		at[a].y = r_transform[3*a+1];
		at[a].z = r_transform[3*a+2];

		printATOM( theFile, 1+a, at[a].res, at+a );
	}	

	fclose(theFile);
}

