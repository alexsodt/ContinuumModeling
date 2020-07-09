#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "util.h"
#include "mutil.h"

// approximate positions of the junctions.
double basis_centers[8][3] = 
{
	{  0.12, 0.63, 0.85 },
	{  0.6,	 0.05, 0.39 },
	{  0.08, 0.38, 0.61 },
	{  0.28, 0.39, 0.37 },
	{  0.34, 0.66, 0.11 },
	{  0.65, 0.89, 0.12 },
	{  0.84, 0.10, 0.61 },
	{  0.83, 0.91, 0.89 }
};

int main( int argc, char **argv )
{
	int xy_dupe = 1;
	int z_dupe = 1;

	if( argc > 1 )
		xy_dupe = atoi(argv[1]);
	if( argc > 2 )
		z_dupe = atoi(argv[2]);

	int nmain = 8;
	nmain *= (1+2*xy_dupe) * (1+2*xy_dupe) * (1+2*z_dupe);
	
	int main_bonds[nmain][3];
	int nmain_bonds[nmain];
	memset( nmain_bonds, 0, sizeof(int) * nmain );

	double orientations[nmain][3];

	double b_cut = 0.5;

	double Lx = (1+xy_dupe);
	double Ly = (1+xy_dupe);
	double Lz = 1 + (1+z_dupe);
	double Ls[3] = {Lx,Ly,Lz};
	double *all_centers = (double *)malloc( sizeof(double) * nmain * 3 );

	nmain = 0;

	for( int dx = 0; dx <= xy_dupe; dx++ )
	for( int dy = 0; dy <= xy_dupe; dy++ )
	for( int dz = 0; dz <= z_dupe; dz++ )
	{
		for( int b = 0; b < 8; b++ )
		{
			all_centers[nmain*3+0] = basis_centers[b][0] + dx;
			all_centers[nmain*3+1] = basis_centers[b][1] + dy;
			all_centers[nmain*3+2] = basis_centers[b][2] + dz;
			nmain++;
		}
	}

	for( int b1 = 0; b1 < nmain; b1++ )
	{
		for( int b2 = 0; b2 < nmain; b2++ )
		{
			if( b2 == b1 ) continue;
	
			double dr[3] = { 
				all_centers[b1*3+0] - all_centers[3*b2+0],
				all_centers[b1*3+1] - all_centers[3*b2+1],
				all_centers[b1*3+2] - all_centers[3*b2+2] };
			while( dr[0] < -Lx/2 ) dr[0] += Lx;
			while( dr[0] >  Lx/2 ) dr[0] -= Lx;
			while( dr[1] < -Ly/2 ) dr[1] += Ly;
			while( dr[1] >  Ly/2 ) dr[1] -= Ly;
			while( dr[2] < -Lz/2 ) dr[2] += Lz;
			while( dr[2] >  Lz/2 ) dr[2] -= Lz;
		
			double ln = normalize(dr);
	
			if( ln < b_cut )
			{
				if( nmain_bonds[b1] >= 3 )
				{
					printf("ERR. nmain_bonds: %d\n", nmain_bonds[b1]+1);
					exit(1);
				}
	
				main_bonds[b1][nmain_bonds[b1]] = b2;
				nmain_bonds[b1]++;
			}
		}

	}



	double dr = 0.25;


	double *all_r = (double *)malloc( sizeof(double) * nmain * 6  * 3 );
	int npts = 0;
	
	int faces[nmain][3][4];
		
	int foffs[3][4] = { 
			    { 2, 0, 3, 5 },  // face to zero bond
			    { 0, 1, 4, 3 }, // face to one bondd
			    { 1, 2, 5, 4 }  // could be end-cap.
				};
		
	int sub_bonds[6][3] =
	{
		{ 1, 2, 3 },
		{ 0, 2, 4 },
		{ 0, 1, 5 },
		{ 0, 4, 5 },
		{ 1, 3, 5 },
		{ 2, 3, 4 }
	};

	for( int b1 = 0; b1 < nmain; b1++ )
	{
		int b2 = main_bonds[b1][0];
		int b3 = main_bonds[b1][1];
		double dr1[3] = { 
				all_centers[3*b2+0] - all_centers[3*b1+0],
				all_centers[3*b2+1] - all_centers[3*b1+1],
				all_centers[3*b2+2] - all_centers[3*b1+2] };
			while( dr1[0] < -Lx/2 ) dr1[0] += Lx;
			while( dr1[0] >  Lx/2 ) dr1[0] -= Lx;
			while( dr1[1] < -Ly/2 ) dr1[1] += Ly;
			while( dr1[1] >  Ly/2 ) dr1[1] -= Ly;
			while( dr1[2] < -Lz/2 ) dr1[2] += Lz;
			while( dr1[2] >  Lz/2 ) dr1[2] -= Lz;
		

		double dr2[3] = { 
				all_centers[3*b3+0] - all_centers[3*b1+0],
				all_centers[3*b3+1] - all_centers[3*b1+1],
				all_centers[3*b3+2] - all_centers[3*b1+2] };
			while( dr2[0] < -Lx/2 ) dr2[0] += Lx;
			while( dr2[0] >  Lx/2 ) dr2[0] -= Lx;
			while( dr2[1] < -Ly/2 ) dr2[1] += Ly;
			while( dr2[1] >  Ly/2 ) dr2[1] -= Ly;
			while( dr2[2] < -Lz/2 ) dr2[2] += Lz;
			while( dr2[2] >  Lz/2 ) dr2[2] -= Lz;

		normalize(dr1);
		normalize(dr2);

		cross( dr1, dr2, orientations[b1] );

		normalize( orientations[b1] );				

		// make a prism.

		double ortho[3];

		cross( orientations[b1], dr1, ortho );

		double dp_ortho = ortho[0] * dr2[0] + ortho[1] * dr2[1] + ortho[2] * dr2[2];

		if( dp_ortho < 0 )
		{
			ortho[0]*=-1;
			ortho[1]*=-1;
			ortho[2]*=-1;
		}

		normalize(ortho);

		for( int dz = -1; dz <= 1; dz +=2 )
		for( int ith = 0; ith < 3; ith++ )
		{
			double l = (dr/2)/sin(60*M_PI/180.0);

			double phi = (2*M_PI/6) + ith * (2*M_PI/3);
			all_r[3*npts+0] = all_centers[3*b1+0] - orientations[b1][0] * dr/2 * dz + cos(phi) * dr1[0]*l + sin(phi) * ortho[0]*l;
			all_r[3*npts+1] = all_centers[3*b1+1] - orientations[b1][1] * dr/2 * dz + cos(phi) * dr1[1]*l + sin(phi) * ortho[1]*l;
			all_r[3*npts+2] = all_centers[3*b1+2] - orientations[b1][2] * dr/2 * dz + cos(phi) * dr1[2]*l + sin(phi) * ortho[2]*l;
			npts++;
		}
	}	

	int max_bonds = 20;
	int *bonds = (int *)malloc( sizeof(int) * npts * max_bonds );
	int *nbonds = (int *)malloc( sizeof(int) * npts );
	memset( nbonds, 0, sizeof(int) * npts );

	int n_sub_pts = 6;
	for( int b = 0; b < nmain; b++ )
	{

		int b_off = b * n_sub_pts;
		for( int sub = 0; sub < n_sub_pts; sub++ )
		for( int t = 0; t < 3; t++ )
		{
			int p = b_off + sub;

			bonds[p*max_bonds+nbonds[p]] = b_off + sub_bonds[sub][t];
			nbonds[p]++;
		}
	}

	int end_cap[nmain][3];
	double fcens[nmain][3][3];


	for( int b = 0; b < nmain; b++ )
	{	
		int b_off = b * n_sub_pts;

		for( int f = 0; f < 3; f++ )
		{
			fcens[b][f][0] = 0;
			fcens[b][f][1] = 0;
			fcens[b][f][2] = 0;

			for( int x = 0; x < 4; x++ )
			{
				faces[b][f][x] = b_off + foffs[f][x]; 

				fcens[b][f][0] += all_r[3*faces[b][f][x]+0]/4;
				fcens[b][f][1] += all_r[3*faces[b][f][x]+1]/4;
				fcens[b][f][2] += all_r[3*faces[b][f][x]+2]/4;
			}
		}

		end_cap[b][0]=0;
		end_cap[b][1]=0;
		end_cap[b][2]=0;
		if( nmain_bonds[b] == 2 )
			end_cap[b][2] = 1; 
	}

	// connect the faces.

	int done[nmain][3];
	for( int b = 0; b < nmain; b++ )
	for( int f = 0; f < 3; f++ )
	{
		done[b][f] = 0;
		if( end_cap[b][f] )
			done[b][f] = 1;
	}
	int all_done = 0;
	int one_more_loop = 0;

	while( !all_done )
	{
		int got_one = 0;
		int best_b1 = -1;
		int best_f1 = -1;
		int o_best_b2 = -1;
		int o_best_f2 = -1;
		double best_or = 1e10;
		for( int b1 = 0; b1 < nmain; b1++ )
		{
	//		if( b1 != 0 ) continue;
			for( int f1 = 0; f1 < 3; f1++ )
			{
	//			if( f1 != 2 ) continue;
				if( done[b1][f1] ) continue;
				
				int best_b = -1;	
				int best_f = -1;
	
				double small_r = 1e10;

				for( int xb = 0; xb < nmain_bonds[b1]; xb++ )
				{
					int b2 = main_bonds[b1][xb];
	
					for( int f2 = 0; f2 < 3; f2++ )
					{
						if( done[b2][f2] ) continue;
						double dr[3] = { fcens[b1][f1][0] - fcens[b2][f2][0],
								 fcens[b1][f1][1] - fcens[b2][f2][1],
								 fcens[b1][f1][2] - fcens[b2][f2][2] };
						for( int c = 0; c < 3; c++ )
						{
							while( dr[c] < -Ls[c]/2 ) dr[c] += Ls[c];
							while( dr[c] >  Ls[c]/2 ) dr[c] -= Ls[c];
						}
	
						double l = normalize(dr);
	
						if( l < small_r )
						{
							best_b = b2;
							best_f = f2;
							small_r =l;
						}
					}
				}	

				if( small_r < best_or )
				{
					got_one = 1;
					best_or = small_r;
					best_b1 = b1;
					best_f1 = f1;
					o_best_b2 = best_b;
					o_best_f2 = best_f;
				}
			}
		}
			
		if( got_one  )
		{	
			int b1 = best_b1;
			int f1 = best_f1;
			int rev_cycle = 0;
			int f2 = o_best_f2;
			int b2 = o_best_b2;
			int face_1[4] = { faces[b1][f1][0], faces[b1][f1][1], faces[b1][f1][2], faces[b1][f1][3] };
			int face_2[4] = { faces[b2][f2][0], faces[b2][f2][1], faces[b2][f2][2], faces[b2][f2][3] }; 

			double small_r = 1e10;	
			int t_start = -1;
			for( int s = 0; s < 4; s++ )
			{
				double dr[3] = { 
					all_r[3*face_1[0]+0] - all_r[3*face_2[s]+0],
					all_r[3*face_1[0]+1] - all_r[3*face_2[s]+1],
					all_r[3*face_1[0]+2] - all_r[3*face_2[s]+2] };
	
				for( int c = 0; c < 3; c++ )
				{
					while( dr[c] < -Ls[c]/2 ) dr[c] += Ls[c];
					while( dr[c] >  Ls[c]/2 ) dr[c] -= Ls[c];
				}
	
				double l = normalize(dr);
	
				if( l < small_r )
				{
					small_r = l;
					t_start = s;
				}
			}
			// the vector separating.
	
			double vec[3] = { fcens[b2][f2][0] - fcens[b1][f1][0],
					  fcens[b2][f2][1] - fcens[b1][f1][1],
					  fcens[b2][f2][2] - fcens[b1][f1][2] };
			
			for( int c = 0; c < 3; c++ )
			{
				while( vec[c] < -Ls[c]/2 ) vec[c] += Ls[c];
				while( vec[c] >  Ls[c]/2 ) vec[c] -= Ls[c];
			}
	
			double cross1[3], cross2[3];
			double dr1_A[3], dr1_B[3];
			double dr2_A[3], dr2_B[3];
	
			dr1_A[0] = all_r[3*face_1[2]+0] - all_r[3*face_1[1]+0];
			dr1_A[1] = all_r[3*face_1[2]+1] - all_r[3*face_1[1]+1];
			dr1_A[2] = all_r[3*face_1[2]+2] - all_r[3*face_1[1]+2];
			
			dr1_B[0] = all_r[3*face_1[0]+0] - all_r[3*face_1[1]+0];
			dr1_B[1] = all_r[3*face_1[0]+1] - all_r[3*face_1[1]+1];
			dr1_B[2] = all_r[3*face_1[0]+2] - all_r[3*face_1[1]+2];
	
			dr2_A[0] = all_r[3*face_2[2]+0] - all_r[3*face_2[1]+0];
			dr2_A[1] = all_r[3*face_2[2]+1] - all_r[3*face_2[1]+1];
			dr2_A[2] = all_r[3*face_2[2]+2] - all_r[3*face_2[1]+2];
			
			dr2_B[0] = all_r[3*face_2[0]+0] - all_r[3*face_2[1]+0];
			dr2_B[1] = all_r[3*face_2[0]+1] - all_r[3*face_2[1]+1];
			dr2_B[2] = all_r[3*face_2[0]+2] - all_r[3*face_2[1]+2];
	
			cross( dr1_A, dr1_B, cross1 );
			cross( dr2_A, dr2_B, cross2 );
			normalize(cross1);
			normalize(cross2);
			double dp = cross1[0] * cross2[0] + cross1[1] * cross2[1] + cross1[2] * cross2[2];
	
			int incr = 1;
			if( dp < 0 )
				incr = -1;
	
			int ind = t_start;
			for( int s = 0; s < 4; s++, ind += incr )
			{
				if( ind < 0 ) ind = 3;
				if( ind >= 4 ) ind = 0;
	
				face_2[s] = faces[b2][f2][ind];
			}
	
			for( int x = 0; x < 4; x++ )
			{
				int xp1 = x+1;
				if( xp1 >= 4 ) xp1 = 0;
				bonds[face_1[x]*max_bonds+nbonds[face_1[x]]] = face_2[x];
				nbonds[face_1[x]]+=1;
				bonds[face_1[x]*max_bonds+nbonds[face_1[x]]] = face_2[xp1];
				nbonds[face_1[x]]+=1;
				bonds[face_2[x]*max_bonds+nbonds[face_2[x]]] = face_1[x];
				nbonds[face_2[x]]+=1;
				bonds[face_2[xp1]*max_bonds+nbonds[face_2[xp1]]] = face_1[x];
				nbonds[face_2[xp1]]+=1;
			}
	
			done[b1][f1]=1;
			done[b2][f2]=1;
		}
		else if( !one_more_loop )
		{	
			one_more_loop = 1;

		}
		else break;
	}


	int *all_bonds = (int *)malloc( sizeof(int) * npts * max_bonds * 2 );
	int all_nbonds = 0;

	for( int b = 0; b < npts; b++ )
	{
		for( int xb = 0; xb < nbonds[b]; xb++ )
		{
			int b2 = bonds[b*max_bonds+xb];

			if( b2 > b )
			{
				double dr[3] = { all_r[3*b+0] - all_r[3*b2+0], all_r[3*b+1] - all_r[3*b2+1], all_r[3*b+2]-all_r[3*b2+2]};
				double ln = normalize(dr);
				if( ln < 0.6 )
				{
					all_bonds[2*all_nbonds+0] = b;
					all_bonds[2*all_nbonds+1] = b2;
					all_nbonds++;
				}
			}
		}
	}

	// handle caps.

	for( int b = 0; b < nmain; b++ )
	{
		for( int f = 0; f < 3; f++ )
		{
			if( end_cap[b][f] )
			{
				int valence[2] = 
				{
					nbonds[faces[b][f][0]] + nbonds[faces[b][f][2]],	
					nbonds[faces[b][f][1]] + nbonds[faces[b][f][3]]	
				};

				if( valence[0] <= valence[1] )
				{
					int v1 = faces[b][f][0];
					int v2 = faces[b][f][2];
					bonds[v1*max_bonds+nbonds[v1]] = v2;
					nbonds[v1]++;
					bonds[v2*max_bonds+nbonds[v2]] = v1;
					nbonds[v2]++;
				}
				else
				{
					int v1 = faces[b][f][1];
					int v2 = faces[b][f][3];
					bonds[v1*max_bonds+nbonds[v1]] = v2;
					nbonds[v1]++;
					bonds[v2*max_bonds+nbonds[v2]] = v1;
					nbonds[v2]++;
				}
			}
		}
	}

#if 1
	FILE *qiiPSF = fopen("gyroid.psf","w");
	writePSF(qiiPSF, npts, NULL, all_bonds, all_nbonds );
	fclose(qiiPSF);

	FILE *qiiXYZ = fopen("gyroid.xyz","w");
	fprintf(qiiXYZ,"%d\n", npts );
	fprintf(qiiXYZ,"ok\n");
	for( int p = 0; p  < npts; p++ )
		fprintf(qiiXYZ,"C %lf %lf %lf\n", all_r[3*p+0], all_r[3*p+1], all_r[3*p+2] );
#endif

	double L = 1;
	printf("3d saved surface\n");

	printf("%lf 0.0 0.0\n", Lx);
	printf("0.0 %lf 0.0\n", Ly);
	printf("0.0 0.0 %lf\n", Lz);

	for( int p = 0; p < npts; p++ )
	{
		printf("%d %lf %lf %lf %d", p, all_r[p*3+0], all_r[3*p+1], all_r[3*p+2], nbonds[p] ); 

		for( int b = 0; b < nbonds[p]; b++ )
			printf(" %d", bonds[p*max_bonds+b] );
		printf("\n");
	}
}




