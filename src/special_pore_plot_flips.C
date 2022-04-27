// This is a special purpose routine to mark the location of flip-flop events on a fusion pore.
//
//
//
//

// Else, draw only bilayer midplane
#define DRAW_LEAFLETS 	1
#define PLOT_CNT	1

#define __special_pore_plot_flipsC__
/*
 * Routines for gathering nearest-point curvature from a surface.
 * Initial purpose is to gather lipid's sampled curvature given a fit to the surface.
 * */
#include <stdio.h>
#include "simulation.h"
#include "interp.h"
#include "gather.h"
#include "M_matrix.h"
#include "meshCollisionLib.h"
#include "io_mol_read.h"
#include "input.h"
#include "util.h"
#include "mutil.h"
#include "pcomplex.h"
#include <math.h>
#include "maxc.h"

static const char *ps_header = 
"%%!PS-Adobe-3.0 EPSF-3.0\n"
"%%%%BoundingBox: 0 0 %lf %lf\n"
"%%%%Orientation: Portrait\n"
"%%%%Pages: 1\n"
"%%%%EndComments\n"
"%%%%Page: 1 1\n";

static const char *header2 =
"%%EndComments\n"
"\n"
"/arrowdict 14 dict def \n"
"arrowdict begin\n"
" /mtrx matrix def \n"
"end\n"
"\n"
"/arrow\n"
" { arrowdict begin\n"
"   /headlength exch def \n"
"   /halfheadthickness exch 2 div def \n"
"   /halfthickness exch 2 div def \n"
"   /tipy exch def /tipx exch def \n"
"   /taily exch def /tailx exch def \n"
"\n"
"   /dx tipx tailx sub def \n"
"   /dy tipy taily sub def \n"
"   /arrowlength dx dx mul dy dy mul add sqrt def \n"
"   /angle dy dx atan def \n"
"   /base arrowlength headlength sub def \n"
"   /savematrix mtrx currentmatrix def \n"
" \n"
"   tailx taily translate\n"
"   angle rotate\n"
"   0 halfthickness neg moveto\n"
"   base halfthickness neg lineto\n"
"   base halfheadthickness neg lineto\n"
"   arrowlength 0 lineto\n"
"   base halfheadthickness lineto\n"
"   base halfthickness lineto\n"
"   0 halfthickness lineto\n"
"   closepath\n"
"   savematrix setmatrix\n"
"  end\n"
"} def\n";

static double DPI = 300;
static double width_inches = 3.25;

void printAt( FILE *theFile, double x, double y, const char *string, double font_size);

void Simulation::plotFlipFlopEvents( parameterBlock *params )
{
	double center[3], r, h;
	getFusionPoreData( center, &r, &h ); 

	// this controls which way the water flow is drawn, and which way the shape is shown moving.
	int flow_and_force = params->special_force_plot;



	printf("Center: %lf %lf %lf\n", center[0], center[1], center[2] );

	if( !params->flipFile )
	{
		printf("printFlipFlopEvents called with flipFile input not set.\n");
		exit(1);
	}

	surface *theSurface = allSurfaces->theSurface;
	double *rsurf = allSurfaces->r;

	FILE *flipEvents = fopen(params->flipFile, "r" );

	if( !flipEvents )
	{
		printf("Couldn't open file '%s' of flip events.\n", params->flipFile );
		exit(1);
	}

	char *buffer = (char *)malloc( sizeof(char) * 100000 );

	char fileName[256 + strlen(params->jobName)];

	sprintf(fileName, "%s_flips.ps", params->jobName);	
	FILE *writePS = fopen(fileName,"w");

	double use_size[2] = { 200, 200 };


	double width = width_inches * DPI;
	double height = width_inches * DPI * use_size[1] / use_size[0];

	fprintf(writePS, ps_header, width, height );
	fprintf(writePS, "%s", header2 );
	fprintf(writePS, "%lf %lf translate\n", width/2, height/2 );
	fprintf(writePS, "%lf %lf scale\n", width / (use_size[0]), width / (use_size[1]) );
	fprintf(writePS, "%lf %lf translate\n", -center[1], -center[2] );

	
	int npts = 400;
	double cen[3];	
	int *f_ring = (int *)malloc( sizeof(int) * npts );
	double *uv_ring = (double *)malloc( sizeof(double) * npts * 2 );
	double *ring_pts = (double *)malloc( sizeof(double) * npts * 3 );
	int npts_out;
	// generate a trace of the pore, x = 0
	int do_reduce_convex = 2;
	theSurface->get_cut_points( 0, center[0], f_ring, uv_ring, ring_pts, npts, rsurf, &npts_out, cen, 0, do_reduce_convex  );

	
	double a_pt[3] = { ring_pts[0], ring_pts[1], ring_pts[2] };

	double translate[2] = {0,0};
	
	while( a_pt[1] - center[1] + translate[0] < -PBC_vec[1][1]/2 ) translate[0] += PBC_vec[1][1];
	while( a_pt[1] - center[1] + translate[0] >  PBC_vec[1][1]/2 ) translate[0] -= PBC_vec[1][1];
	while( a_pt[2] - center[2] + translate[1] < -PBC_vec[2][2]/2 ) translate[1] += PBC_vec[2][2];
	while( a_pt[2] - center[2] + translate[1] >  PBC_vec[2][2]/2 ) translate[1] -= PBC_vec[2][2];

	fprintf( writePS, "%lf %lf translate\n", translate[0], translate[1] );

	int CNT_f[3];
	double CNT_uv[3][2];
	int CNT_set[3] = {0,0};

	double approx_pt[3][2] = { { center[1] -70, center[2] + 50 }, { center[1]  -70, center[2] -50 }, {center[1]-15, center[2] } };
	double cnt_chi2[3] = { 1e10, 1e10, 1e10 };
	

//#define POINTS

#ifdef POINTS
	for( int i = 0; i < npts_out; i++ )
	{
		double rpt[3] ,npt[3];
		theSurface->evaluateRNRM( f_ring[i], uv_ring[2*i+0], uv_ring[2*i+1], rpt, npt, rsurf ); 

		fprintf(writePS, "newpath\n");
		fprintf(writePS, "%lf %lf %lf 0 360.0 arc\n", rpt[1], rpt[2], 3.0 );

		fprintf(writePS, "fill\n");
	}

#else
	for( int leaflet = 0; leaflet < 2; leaflet++ )
	{
#ifndef DRAW_LEAFLETS 
		if( leaflet == 1 ) 
			continue;
#endif
		int init = 1;
		fprintf(writePS, "newpath\n");
		double ppt[2];
		for( int i = 0; i <= npts_out; i++ )
		{
			int use_i = i;
			if( use_i >= npts_out )
				use_i -= npts_out;
	
			double rpt[3] ,npt[3];
			theSurface->evaluateRNRM( f_ring[use_i], uv_ring[2*use_i+0], uv_ring[2*use_i+1], rpt, npt, rsurf ); 

			for( int cnt = 0; cnt < 3; cnt++ )
			{
				double chi2 =  pow( rpt[1] - approx_pt[cnt][0], 2.0) + pow( rpt[2] - approx_pt[cnt][1],2.0);

				if( chi2 < cnt_chi2[cnt] )
				{
					printf("found %lf %lf near %lf %lf chi2 %le vs %le\n", rpt[1], rpt[2], approx_pt[cnt][0], approx_pt[cnt][1], chi2, cnt_chi2[cnt] );
					cnt_chi2[cnt] = chi2;
					CNT_f[cnt] = f_ring[use_i];
					CNT_uv[cnt][0] = uv_ring[2*use_i+0];
					CNT_uv[cnt][1] = uv_ring[2*use_i+1];
				} 
			}


#ifdef DRAW_LEAFLETS
			if( leaflet == 0 )
			{
				rpt[0] += npt[0] * 13.0;
				rpt[1] += npt[1] * 13.0;
				rpt[2] += npt[2] * 13.0;
			}
			else
			{
				rpt[0] -= npt[0] * 13.0;
				rpt[1] -= npt[1] * 13.0;
				rpt[2] -= npt[2] * 13.0;
			}
#endif
	
			if( init )
			{
				init=0;
				fprintf(writePS, "%lf %lf moveto\n", rpt[1], rpt[2] );
			}
			else
			{
				double dr[2] = { rpt[1] - ppt[0], rpt[2] - ppt[1] };
	
				if( dr[0] < PBC_vec[1][1] /2 && dr[0] > -PBC_vec[1][1]/2 &&
				    dr[1] < PBC_vec[2][2] /2 && dr[1] > -PBC_vec[2][2]/2 )
					fprintf(writePS, "%lf %lf lineto\n", rpt[1], rpt[2]);
				else
				{
					if( !init )
					{
						// this looks wrong but ppt is in 2D, center is in 3D.
						if( ppt[0] < center[1] )
							fprintf(writePS,"0.8 0.8 0.8 setrgbcolor\n");
						else
						fprintf(writePS,"0.0 0.0 0.0 setrgbcolor\n");
						fprintf(writePS, "stroke\n");
					}
					init = 1;
				}
			}
				
			ppt[0] = rpt[1];
			ppt[1] = rpt[2];
		}
		// this looks wrong but ppt is in 2D, center is in 3D.
		if( ppt[0] < center[1] )
			fprintf(writePS,"0.8 0.8 0.8 setrgbcolor\n");
		else
			fprintf(writePS,"0.0 0.0 0.0 setrgbcolor\n");
		if( ! init ) fprintf(writePS, "stroke\n");
	}
#endif
	double font_size = 5.0;

	double CNT_height = 30;
	double CNT_width  = 12;

	if( flow_and_force == -2 || flow_and_force == 2 )
	{
		for( int cnt = 0; cnt < 2; cnt++ )
		{ 
			double rpt[3] ,npt[3];
			theSurface->evaluateRNRM( CNT_f[cnt], CNT_uv[cnt][0], CNT_uv[cnt][1], rpt, npt, rsurf ); 
			
			double nrm_2d[3] = { npt[1], npt[2],0 };
			normalize(nrm_2d);
			double perp_2d[2] = { -nrm_2d[1], nrm_2d[0] };
	
			fprintf( writePS, "0.7 0.7 0.7 setrgbcolor\n");
			fprintf( writePS, "%lf %lf moveto\n", rpt[1] + nrm_2d[0] * CNT_height/2 + perp_2d[0] * CNT_width/2, 
							      rpt[2] + nrm_2d[1] * CNT_height/2 + perp_2d[1] * CNT_width/2  );
			fprintf( writePS, "%lf %lf lineto\n", rpt[1] + nrm_2d[0] * CNT_height/2 - perp_2d[0] * CNT_width/2, 
							      rpt[2] + nrm_2d[1] * CNT_height/2 - perp_2d[1] * CNT_width/2  );
			fprintf( writePS, "%lf %lf lineto\n", rpt[1] - nrm_2d[0] * CNT_height/2 - perp_2d[0] * CNT_width/2, 
							      rpt[2] - nrm_2d[1] * CNT_height/2 - perp_2d[1] * CNT_width/2  );
			fprintf( writePS, "%lf %lf lineto\n", rpt[1] - nrm_2d[0] * CNT_height/2 + perp_2d[0] * CNT_width/2, 
							      rpt[2] - nrm_2d[1] * CNT_height/2 + perp_2d[1] * CNT_width/2  );
			fprintf( writePS, "closepath\n fill\n");

			if( flow_and_force  == -2)
			{
				fprintf( writePS, "0.7 0.2 0.2 setrgbcolor\n");
				// draw flow going out
				
				fprintf( writePS, "%lf %lf %lf %lf %lf %lf %lf arrow\n fill\n",
					rpt[1] - nrm_2d[0] * 1.5 * CNT_height/2, 
					rpt[2] - nrm_2d[1] * 1.5 * CNT_height/2,
					rpt[1] + nrm_2d[0] * 1.5 * CNT_height/2, 
					rpt[2] + nrm_2d[1] * 1.5 * CNT_height/2,
						4., 8., 4. );
			}
			else if( flow_and_force == 2 )
			{
				fprintf( writePS, "0.2 0.2 0.7 setrgbcolor\n");
				// draw flow going out
				
				fprintf( writePS, "%lf %lf %lf %lf %lf %lf %lf arrow\n fill\n",
					rpt[1] + nrm_2d[0] * 1.5 * CNT_height/2, 
					rpt[2] + nrm_2d[1] * 1.5 * CNT_height/2,
					rpt[1] - nrm_2d[0] * 1.5 * CNT_height/2, 
					rpt[2] - nrm_2d[1] * 1.5 * CNT_height/2,
						4., 8., 4. );

			}
			
			if( cnt == 1 )
			{

				printAt( writePS,rpt[1] - nrm_2d[0] * 1.5 * CNT_height/2 + 10, 
					 rpt[2] - nrm_2d[1] * 1.5 * CNT_height/2 - 10.0, "Water flow across CNT", font_size );
			}
				
		}
		if( flow_and_force  == -2 )
		{
			double rpt[3] ,npt[3];
			theSurface->evaluateRNRM( CNT_f[2], CNT_uv[2][0], CNT_uv[2][1], rpt, npt, rsurf ); 
			double nrm_2d[3] = { npt[1], npt[2],0 };
			normalize(nrm_2d);
			
			
			fprintf( writePS, "%lf %lf %lf %lf %lf %lf %lf arrow\n fill\n",
				rpt[1] - nrm_2d[0] * 13.0, 
				rpt[2] - nrm_2d[1] * 13.0,
				rpt[1] - nrm_2d[0] * 25.0, 
				rpt[2] - nrm_2d[1] * 25.0,
					4., 8., 4. );

			printAt( writePS,rpt[1] - nrm_2d[0] * 35.0 , rpt[2] - nrm_2d[1] * 35.0 - 10.0, "Pore narrows", font_size );
		}
		else if( flow_and_force == 2 )
		{
			fprintf( writePS, "0.2 0.2 0.7 setrgbcolor\n");
			double rpt[3] ,npt[3];
			theSurface->evaluateRNRM( CNT_f[2], CNT_uv[2][0], CNT_uv[2][1], rpt, npt, rsurf ); 
			double nrm_2d[3] = { npt[1], npt[2],0 };
			normalize(nrm_2d);
			
			
			fprintf( writePS, "%lf %lf %lf %lf %lf %lf %lf arrow\n fill\n",
				rpt[1] + nrm_2d[0] * 13.0, 
				rpt[2] + nrm_2d[1] * 13.0,
				rpt[1] + nrm_2d[0] * 25.0, 
				rpt[2] + nrm_2d[1] * 25.0,
					4., 8., 4. );

			printAt( writePS,rpt[1] + nrm_2d[0] * 35 , rpt[2] + nrm_2d[1] * 35 - 10.0, "Pore expands", font_size );
		}
	

	}


	while( !feof(flipEvents) )
	{
		getLine( flipEvents, buffer );
		if( feof(flipEvents) ) break;

		char *segName = (char *)malloc( sizeof(char) * (strlen(buffer)+1) );
		char *resName = (char *)malloc( sizeof(char) * (strlen(buffer)+1) );
		char from, to;
		double z;
		int resid;
		int nr = sscanf(buffer, "%s %d %s %c %c z %lf cen %lf %lf %lf", resName, &resid, segName, &from, &to, &z ); 

		if( strlen(buffer) > 5 && nr != 6 )
		{
			printf("Error reading line '%s' of flip file.\n", buffer );
			exit(1);
		}

		if( nr == 6 )
		{
			// plot it ! 
	
			// find the near point, y > center[1]
	
			int near_pt = 0;
		
			double smalld = 1e10;
			for( int p = 0; p < npts_out; p++ )
			{
				double dz = z - ring_pts[3*p+2];

				while( dz   > PBC_vec[2][2]/2 ) dz -= PBC_vec[2][2];
				while( dz  < -PBC_vec[2][2]/2 ) dz += PBC_vec[2][2];
				double dz2 = dz*dz;
	
				if( dz2 < smalld && ring_pts[3*p+1] > center[1] )
				{
					smalld = dz2;
					near_pt = p;
				}
			}	
	
			double rpt[3] ,npt[3];
			theSurface->evaluateRNRM( f_ring[near_pt], uv_ring[2*near_pt+0], uv_ring[2*near_pt+1], rpt, npt, rsurf ); 
	
			double tail[2],head[2];
			
#ifndef DRAW_LEAFLETS
			if( from == 'u' )
			{
				head[0] = rpt[1];
				head[1] = rpt[2];
				tail[0] = rpt[1] + npt[1] * 5;		 
				tail[1] = rpt[2] + npt[2] * 5;		 
				fprintf(writePS, "0.0 0.0 0.5 setrgbcolor\n");
			}
			else
			{
				head[0] = rpt[1];
				head[1] = rpt[2];
				tail[0] = rpt[1] - npt[1] * 5;		 
				tail[1] = rpt[2] - npt[2] * 5;		 
				fprintf(writePS, "0.5 0.5 1.0 setrgbcolor\n");
			}
#else
			if( from == 'u' )
			{
				head[0] = rpt[1] + npt[1] * 5;
				head[1] = rpt[2] + npt[2] * 5;
				tail[0] = rpt[1] + npt[1] * 13;		 
				tail[1] = rpt[2] + npt[2] * 13;		 
				fprintf(writePS, "0.0 0.0 0.5 setrgbcolor\n");
			}
			else
			{
				head[0] = rpt[1] - npt[1] * 5;
				head[1] = rpt[2] - npt[2] * 5;
				tail[0] = rpt[1] - npt[1] * 13;		 
				tail[1] = rpt[2] - npt[2] * 13;		 
				fprintf(writePS, "0.5 0.5 1.0 setrgbcolor\n");
			}
#endif
				
			fprintf(writePS, "newpath\n%lf %lf %lf %lf %lf %lf %lf arrow\nfill\n", tail[0], tail[1], head[0], head[1], 2., 4., 2. );
		}

		free(resName);
		free(segName);
	}
	font_size = 7.0;			
	printAt( writePS,center[1] + 60, center[2] + 75, "Cholesterol flips", font_size );
	
	font_size = 10.0;
	if( flow_and_force < 0 )
		printAt( writePS,center[1] + 40, center[2] - 75, "Large pore", font_size );
	else
		printAt( writePS,center[1] + 40, center[2] - 75, "Small pore", font_size );

	fclose(writePS);
}
			
void printAt( FILE *theFile, double x, double y, const char *string, double font_size)
{
        double fs = font_size;

        fprintf(theFile,"/%s findfont\n", "Helvetica" );
        fprintf(theFile,"%lf scalefont\n", fs );  
        fprintf(theFile,"setfont\n");
        fprintf(theFile,"%lf %lf moveto\n", x, y - 0.35 * fs );
        fprintf(theFile,"(%s) stringwidth\n", string );
        fprintf(theFile,"exch -0.5 mul exch\n");    
        fprintf(theFile,"rmoveto\n");   
        fprintf(theFile,"(%s) show\n", string );
}
