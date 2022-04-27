#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "util.h"
#include "dcd.h"
#include "pdb.h"

int main( int argc, char **argv )
{
	if( argc != 2 )
	{
		printf("Syntax: getPoreInfo file.rho|file.pdb\n");
		exit(1);
	}

	double *rho;
	double *rho_xy;
	double *rho_z;
	double La,Lb,Lc;
	int nx, ny, nz;
	
	char *buffer = (char *)malloc( sizeof(char) * 100000 );

	if( !strcasecmp( argv[1]+strlen(argv[1])-3, "rho" ) )
	{
		FILE *rhoFile = fopen(argv[1], "r");
		if( !rhoFile )
		{
			printf("Cannot open density file '%s'.\n", argv[1] );
			exit(1);
		}

	
		fscanf(rhoFile, "%lf %lf %lf\n", &La, &Lb, &Lc );
		fscanf(rhoFile, "%d %d %d\n", &nx, &ny, &nz );

		rho = (double *)malloc( sizeof(double) * nx * ny * nz );
	
		for( int ix = 0; ix < nx; ix++ )
		for( int iy = 0; iy < ny; iy++ )
		{
			getLine( rhoFile, buffer );
	
			readNDoubles( buffer, rho+ix*ny*nz+iy*nz, nz );
		}
		fclose(rhoFile);
	}
	else if( !strcasecmp( argv[1]+strlen(argv[1])-3, "pdb" ) )
	{
		FILE *theFile = fopen(argv[1],"r");
		loadPSFfromPDB(theFile);
		rewind(theFile);
		int nat = curNAtoms();

		struct atom_rec *at = (struct atom_rec *)malloc( sizeof(struct atom_rec) * nat );
		
		loadPDB( theFile, at, nat );


		double alpha,beta,gamma;
		PBCD( &La, &Lb, &Lc, &alpha, &beta, &gamma );

		if( La < 0 )
		{
			printf("This requires the PBC cell set in the PDB.\n");
			exit(1);
		}
		
		nx = ceil(La/3);
		ny = ceil(Lb/3);
		nz = ceil(Lc/3);
		
		rho = (double *)malloc( sizeof(double) * nx * ny * nz );

		for( int a = 0; a < curNAtoms(); a++ )
		{
			if( strcasecmp( at[a].atname, "C21" ) ) continue;
			double f[3] = { at[a].x / La, at[a].y / Lb, at[a].z / Lc };
			for( int c = 0; c < 3; c++ )
			{
				while( f[c] < 0 ) f[c] += 1;
				while( f[c] >= 1 ) f[c] -= 1;
			}

			int bx = nx * f[0];
			int by = ny * f[1];
			int bz = nz * f[2];
			if( bx >= nx ) bx -= nx;
			if( by >= ny ) by -= ny;
			if( bz >= nz ) bz -= nz;

			rho[bx*ny*nz+by*nz+bz] += 1;
		}
		fclose(theFile);
	}
	
	int bin_range = 20.0 / (Lc/nz);
	
	rho_xy = (double *)malloc( sizeof(double) * nx * ny );
	rho_z = (double *)malloc( sizeof(double) * nz );
	double *rho_z_av = (double *)malloc( sizeof(double) * nz );

	for( int ix = 0; ix < nx; ix++ )
	for( int iy = 0; iy < ny; iy++ )
	{
		rho_xy[ix*ny+iy] = 0;
		for( int iz = 0; iz < nz; iz++ )
			 rho_xy[ix*ny+iy] += rho[ix*ny*nz+iy*nz+iz];
	}

	for( int iz = 0; iz < nz; iz++ )
	{
		for( int ix = 0; ix < nx; ix++ )
		for( int iy = 0; iy < ny; iy++ )
			rho_z[iz] += rho[ix*ny*nz+iy*nz+iz];
#ifdef DEBUG_PRINT
		printf("%d %lf\n", iz, rho_z[iz] );
#endif
	}
	for( int iz = 0; iz < nz; iz++ )
	{
		for( int dz = -bin_range; dz <= bin_range; dz++ )
		{
			int z = iz + dz;

			while( z < 0 ) z += nz;
			while( z >= nz ) z -= nz;

			rho_z_av[iz] += rho_z[z];
		}
//		printf("%d %lf\n", iz, rho_z_av[iz] );
	}

//	memcpy( rho_z, rho_z_av, sizeof(double) * nz );

	int i_max_z = 0;
	double max_z = rho_z[0];
	int i_min_z = 0;
	double min_z = rho_z[0];
	for( int iz = 0; iz  <nz; iz++ )
	{
		if( rho_z[iz] > max_z )
		{
			max_z = rho_z[iz];
			i_max_z = iz;
		}
		if( rho_z[iz] < min_z )
		{
			min_z = rho_z[iz];
			i_min_z = iz;
		}
	}


	int leaflet_mode = 0;

	int trigger = 0;
	int i_leaflets[4] = {-1,-1,-1,-1};
	double leaflets[4]={-1,-1,-1,-1};
	double trigger_value = 0;
	int ileaflet = 0;
	double midplane_thresh = 0.5;
	int ntries= 0;
	while( ileaflet != 4 && ntries < 10 )
	{
		ileaflet=0;
		for( int iz = 0; iz < nz; iz++ )
		{
			int use_z = iz + i_min_z;
			if( use_z >= nz )
				use_z -= nz;	
			// if we find a big peak count it until it comes down again.
			if( !trigger && rho_z[use_z] > midplane_thresh * max_z ) 
			{
				trigger = 1;
				trigger_value = rho_z[use_z];
				leaflets[ileaflet] = use_z;
				i_leaflets[ileaflet] = use_z;
			}
	
			if( trigger && rho_z[use_z] > trigger_value )
			{
				trigger_value = rho_z[use_z];
				leaflets[ileaflet] = use_z;
				i_leaflets[ileaflet] = use_z;
			}
	
			if( trigger && rho_z[use_z] < midplane_thresh * trigger_value )
			{
				trigger = 0;
				ileaflet++;
			}
		
			if( ileaflet == 4 ) break;
		}
#ifdef DEBUG_PRINT
		printf("Midplane_thresh: %lf ileaflets: %d %d %d %d\n", midplane_thresh, i_leaflets[0], i_leaflets[1], i_leaflets[2], i_leaflets[3] );
#endif
		midplane_thresh += 0.025;
		ntries++;
	}

	if( ileaflet != 4 )
	{
		ileaflet = 0;
		leaflet_mode = 1;

		// look for two indicators.

		trigger = 0;
		for( int iz = 0; iz < nz; iz++ )
		{
			int use_z = iz + i_min_z;
			if( use_z >= nz )
			use_z -= nz;	
			// if we find a big peak count it until it comes down again.
			if( !trigger && rho_z[use_z] > 0.5 * max_z ) 
			{
				trigger = 1;
				trigger_value = rho_z[use_z];
				leaflets[ileaflet] = use_z;
				i_leaflets[ileaflet] = use_z;
			}
	
			if( trigger && rho_z[use_z] > trigger_value )
			{
				trigger_value = rho_z[use_z];
				leaflets[ileaflet] = use_z;
				i_leaflets[ileaflet] = use_z;
			}
	
			if( trigger && (rho_z[use_z] < 0.5 * trigger_value || use_z - leaflets[ileaflet] > 2*bin_range) )
			{
				trigger = 0;
				ileaflet++;
			}
		
			if( ileaflet == 4 ) break;
		}

		int leaflet_1 = i_leaflets[0];
		int leaflet_2 = i_leaflets[1];

		int dbin = nz * 15.0 / Lc; 

		i_leaflets[0] = leaflet_1 - dbin;
		while( i_leaflets[0] < 0 ) i_leaflets[0] += nz;
		i_leaflets[1] = leaflet_1 + dbin;
		while( i_leaflets[1] >= nz ) i_leaflets[1] -= nz;
		
		i_leaflets[2] = leaflet_2 - dbin;
		while( i_leaflets[2] < 0 ) i_leaflets[2] += nz;
		i_leaflets[3] = leaflet_2 + dbin;
		while( i_leaflets[3] >= nz ) i_leaflets[3] -= nz;
		
		leaflets[0] = i_leaflets[0];	
		leaflets[1] = i_leaflets[1];	
		leaflets[2] = i_leaflets[2];	
		leaflets[3] = i_leaflets[3];	
	}

#ifdef DEBUG_PRINT
	printf("i_leaflets %d %d ; %d %d\n", i_leaflets[0], i_leaflets[1], i_leaflets[2], i_leaflets[3] );
#endif
	int swap_upper_lower = 0; // we may have picked the wrong assignment of bilayers (need the pore between them).


	int midp1 = ((leaflets[0]+leaflets[1]) + (leaflets[2]+leaflets[3]))/4;
	int midp2 = ((leaflets[0]+leaflets[1]+nz+nz) + (leaflets[2]+leaflets[3]))/4;

	while( midp1 < 0 ) midp1 += nz;
	while( midp1 >= nz) midp1 -= nz;
	
	while( midp2 < 0 ) midp2 += nz;
	while( midp2 >= nz) midp2 -= nz;

//	printf("midp1: %d midp2: %d rhozm1: %le rhozm2: %le\n", midp1, midp2, rho_z[midp1], rho_z[midp2] );

	if( rho_z[midp1] < rho_z[midp2] && leaflets[0] < leaflets[2] )
		swap_upper_lower =1;
	else if( rho_z[midp1] > rho_z[midp2] && leaflets[0] > leaflets[2] )
		swap_upper_lower = 1;

	if( swap_upper_lower )
	{
//		printf("Swapping.\n");
		int t[2] = { leaflets[2], leaflets[3] };
		leaflets[2] = leaflets[0];
		leaflets[3] = leaflets[1];
		leaflets[0] = t[0];
		leaflets[1] = t[1];
	}

	double lz[4] =
	{
		Lc * (leaflets[0]+0.5) / nz,
		Lc * (leaflets[1]+0.5) / nz,
		Lc * (leaflets[2]+0.5) / nz,
		Lc * (leaflets[3]+0.5) / nz
	};

	if( lz[2] < lz[0] )
	{
		lz[2] += Lc;
		lz[3] += Lc;
	}

	
	//printf("leaflets: %lf %lf %lf %lf\n", leaflets[0], leaflets[1], leaflets[2], leaflets[3] );

#ifdef WRAP_STANDARD
	for( int dz = 1; dz < 4; dz++ )	
	{
		while( lz[dz] - lz[dz-1] < -Lc/2 )
		{
			for( int ddz = dz; ddz < 4; ddz++ )
				lz[ddz] += Lc;
		}
		
		while( lz[dz] - lz[dz-1] > Lc/2 )
		{
			for( int ddz = dz; ddz < 4; ddz++ )
				lz[ddz] -= Lc;
		}
	}	
#else
		
	while( lz[1] - lz[0] < -Lc/2 )
	{
		lz[1] += Lc;
	}
	while( lz[1] - lz[0] >= Lc/2 )
	{
		lz[1] -= Lc;
	}
	
	while( lz[3] - lz[2] < -Lc/2 )
	{
		lz[3] += Lc;
	}
	while( lz[3] - lz[2] >= Lc/2 )
	{
		lz[3] -= Lc;
	}

/*
	// wrap such that we put water outside.

	// loop "down".

	int down_lim = i_leaflets[3];

	while( down_lim < i_leaflets[0] )
		down_lim += nz; 
	while( down_lim > i_leaflets[0] )
		down_lim -= nz; 

	double sum_down = 0;
	for( int lz = i_leaflets[0]; lz >= down_lim; lz-- )
	{
		int use_lz = lz;
		while( use_lz < 0 )
			use_lz += nz;
		sum_down += rho_z[lz];
	} 
	
	int up_lim = i_leaflets[2];

	while( up_lim > i_leaflets[2] )
		up_lim += nz; 
	while( up_lim > i_leaflets[1] )
		up_lim -= nz; 

	double sum_up = 0;
	for( int lz = i_leaflets[1]; lz <= up_lim; lz++ )
	{
		int use_lz = lz;
		while( use_lz >= nz)
			use_lz -= nz;
		sum_up += rho_z[lz];
	}
*/

	double sum_up=0,sum_down=0;

	int loop_high = i_leaflets[0];
	int loop_low  = i_leaflets[3];
	
	while( loop_high > loop_low )
		loop_high -= nz; // reset
	while( loop_high < loop_low )
		loop_high += nz; // exactly one PBC
	
	for( int ilx = loop_high; ilx >= loop_low; ilx-- )
	{
		int il = ilx;

		while( il < 0 ) il += nz;
		while( il >= nz ) il -= nz;

		sum_up += rho_z[il];
	}
	loop_low = i_leaflets[1];
	loop_high = i_leaflets[2];
	
	while( loop_high > loop_low )
		loop_high -= nz; // reset
	while( loop_high < loop_low )
		loop_high += nz; // exactly one PBC
	
	for( int ilx = loop_low; ilx <= loop_high; ilx++ )
	{
		int il = ilx;

		while( il < 0 ) il += nz;
		while( il >= nz ) il -= nz;

		sum_down += rho_z[il];
	}
//	printf("sum_up: %lf sum_down: %lf\n", sum_up, sum_down );

	// 2/3 is the upper
	// 0/1 is the lower, it's z value must be below.

//	printf("leaflets: %lf %lf, %lf %lf\n", leaflets[0], leaflets[1], leaflets[2], leaflets[3] );
	while( leaflets[1] < leaflets[0] )
		leaflets[1] += nz; 
	while( leaflets[3] < leaflets[2] )
		leaflets[3] += nz; 

	
	while( leaflets[2] < leaflets[1] )
	{
		leaflets[2] += nz;
		leaflets[3] += nz;
	}
	while( leaflets[2] - leaflets[1] >= nz )
	{
		leaflets[2] -= nz;
		leaflets[3] -= nz;
	}

	while( leaflets[0] >= nz )
	{
		leaflets[0] -= nz;
		leaflets[1] -= nz;
		leaflets[2] -= nz;
		leaflets[3] -= nz;
	}


	lz[0] = leaflets[0] * Lc / nz; 
	lz[1] = leaflets[1] * Lc / nz; 
	lz[2] = leaflets[2] * Lc / nz; 
	lz[3] = leaflets[3] * Lc / nz; 

//	printf("LZ: lower: %lf %lf upper: %lf %lf\n", lz[0], lz[1], lz[2], lz[3] );

	midp1 = ((leaflets[0]+leaflets[1]) + (leaflets[2]+leaflets[3]))/4;
	midp2 = ((leaflets[0]+leaflets[1]+nz+nz) + (leaflets[2]+leaflets[3]))/4;


	while( midp1 < 0 ) midp1 += nz;
	while( midp1 >= nz) midp1 -= nz;
	
	while( midp2 < 0 ) midp2 += nz;
	while( midp2 >= nz) midp2 -= nz;
	
//	printf("midp1: %d val %lf\n", midp1, rho_z[midp1] ); 
//	printf("midp2: %d val %lf\n", midp2, rho_z[midp2] ); 

	if( rho_z[midp1] < rho_z[midp2] )
	{
		double t[2] = {lz[2],lz[3]};
		lz[2] = lz[0];
		lz[3] = lz[1];
		lz[0] = t[0];
		lz[1] = t[1];
	}

		while( lz[2] < lz[1] )
		{
			lz[2] += Lc;
			lz[3] += Lc;
		}
		while( lz[2] - lz[1] >= Lc )
		{
			lz[2] -= Lc;
			lz[3] -= Lc;
		}



//	printf("Upper: %lf Lower %lf\n", lz[2], lz[0] );

#endif
//	printf("lz2: %le %le %le %le\n", lz[0], lz[1], lz[2], lz[3] );

	if( lz[3] < lz[1] )
	{
		while( lz[3] < lz[1] )
		{
			lz[3] += nz;
			lz[2] += nz;
		}
	}

	double thickness = fabs((lz[1]+lz[0])/2-(lz[3]+lz[2])/2);
	
	double wrapto = 0;

#ifdef WRAP_STANDARD
	double best_chi2 = 1e10;
	int nbins = nz;
	
	for( int zb = 0; zb < nz; zb++ )
	{
	        double zv = Lc * (zb+0.5) / (double)nz;
//		printf("zv: %lf rho: %lf\n", zv, rho_z[zb] );
	
	         int zlow  = zb- nz/2;
	         int zhigh = zlow + nz;
	
	         double lchi2 = 0;
	         for( int iz = zlow; iz < zhigh; iz++ )
	         {
	                 double dz = Lc * (iz+0.5) / nz - zv;
	
	                 int iiz = iz;
	                 while( iiz < 0 ) iiz += nz;
	                 while( iiz >= nz ) iiz -= nz;
	
	                 lchi2 += rho_z[iiz] * (dz) * (dz);
	         }
	
	         if( lchi2 < best_chi2 )
	         {
	                 best_chi2 = lchi2;
	                 wrapto = zv;
	         }
	}
#else
	wrapto = (lz[1]+lz[0]+lz[3]+lz[2])/4;

#endif

	while( wrapto < 0 )
		wrapto += Lc;
	while( wrapto >= Lc )
		wrapto -= Lc;
//	printf("wrapto: %lf\n",wrapto );
	
	double com[3] = {0,0,0};
	double ncom[3] = {0,0,0};

	com[2] = wrapto;
	ncom[2] = 1;

	// find the xy max
	
	double max_val = rho_xy[0];
	double min_val = rho_xy[0];
	double av_val  = 0;
	int min_xy = 0;
	int min_x = 0;
	int min_y = 0;

	for( int ix = 0; ix < nx; ix++ )
	for( int iy = 0; iy < ny; iy++ )
	{
		if( rho_xy[ix*ny+iy] > max_val )
			max_val = rho_xy[ix*ny+iy];
		if( rho_xy[ix*ny+iy] < min_val )
		{
			min_val = rho_xy[ix*ny+iy];
			min_x = ix;
			min_y = iy;
		}
		av_val += rho_xy[ix*ny+iy];
	}		


	com[0] = (min_x+0.5) * La / nx;
	com[1] = (min_y+0.5) * Lb / ny;
	ncom[0] += 1;
	ncom[1] += 1;

	av_val /= nx*ny;
	
	int *contig = (int *)malloc( sizeof(int) * nx * ny );
	memset( contig, 0, sizeof(int) * nx *ny );
	int *list = (int *)malloc( sizeof(int) * nx * ny );
	int nlist = 0;
	list[nlist++] = min_x*ny+min_y;
	contig[min_x*ny+min_y] = 1;

	double thresh = min_val + (max_val-min_val) * 0.025; 

	int done = 0;
	

	while( !done )
	{
		done = 1;
		
		for( int il = 0; il < nlist; il++ )
		{
			int t_x = list[il] / ny;
			int t_y = list[il] % ny;

			for( int dx = -1; dx <= 1; dx++ )
			for( int dy = -1; dy <= 1; dy++ )
			{
				int n_x = t_x+dx;
				int n_y = t_y+dy;

				if( n_x < 0 ) n_x += nx;
				if( n_x >= nx) n_x -= nx;
				if( n_y < 0 ) n_y += ny;
				if( n_y >= ny ) n_y -= ny;

				if( contig[n_x*ny+n_y] ) continue;

				if( rho_xy[n_x*ny+n_y] > thresh ) continue;

				list[nlist] = n_x*ny+n_y;
				nlist++;

				contig[n_x*ny+n_y] = 1;

				done = 0;

				int del[2] = { n_x - min_x, n_y - min_y };

				if( del[0] < -nx/2 ) del[0] += nx;
				if( del[0] >  nx/2 ) del[0] -= nx;
				if( del[1] < -ny/2 ) del[1] += ny;
				if( del[1] >  ny/2 ) del[1] -= ny;

				com[0] += (min_x + del[0])*La/nx;
				com[1] += (min_y + del[1])*Lb/ny;

				ncom[0] += 1;
				ncom[1] += 1;
			}
		}
	}

	double thresh_area = nlist * (La/nx) * (Lb/ny);
	
	double r_thresh = sqrt(thresh_area/M_PI);

	printf("xcom: %lf ycom: %lf zcom: %lf rthresh %lf thickness: %lf\n", com[0]/ncom[0], com[1]/ncom[1], com[2] / ncom[2], r_thresh, thickness);

	return -1;
}
