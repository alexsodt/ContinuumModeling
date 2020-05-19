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
	int main_bonds[8][3];
	int nmain_bonds[8];
	memset( nmain_bonds, 0, sizeof(int) * 8 );

	double orientations[8][3];

	double b_cut = 0.5;

	for( int b1 = 0; b1 < 8; b1++ )
	{
		for( int b2 = 0; b2 < 8; b2++ )
		{
			if( b2 == b1 ) continue;
	
			double dr[3] = { 
				basis_centers[b1][0] - basis_centers[b2][0],
				basis_centers[b1][1] - basis_centers[b2][1],
				basis_centers[b1][2] - basis_centers[b2][2] };
			while( dr[0] < -0.5 ) dr[0] += 1;
			while( dr[0] >  0.5 ) dr[0] -= 1;
			while( dr[1] < -0.5 ) dr[1] += 1;
			while( dr[1] >  0.5 ) dr[1] -= 1;
			while( dr[2] < -0.5 ) dr[2] += 1;
			while( dr[2] >  0.5 ) dr[2] -= 1;
		
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

	double *all_r = (double *)malloc( sizeof(double) * 8 * 6  * 3 );
	int npts = 0;

	double dr = 0.25;

	for( int b1 = 0; b1 < 8; b1++ )
	{
		int b2 = main_bonds[b1][0];
		int b3 = main_bonds[b1][1];
		double dr1[3] = { 
				basis_centers[b2][0] - basis_centers[b1][0],
				basis_centers[b2][1] - basis_centers[b1][1],
				basis_centers[b2][2] - basis_centers[b1][2] };
			while( dr1[0] < -0.5 ) dr1[0] += 1;
			while( dr1[0] >  0.5 ) dr1[0] -= 1;
			while( dr1[1] < -0.5 ) dr1[1] += 1;
			while( dr1[1] >  0.5 ) dr1[1] -= 1;
			while( dr1[2] < -0.5 ) dr1[2] += 1;
			while( dr1[2] >  0.5 ) dr1[2] -= 1;
		

		double dr2[3] = { 
				basis_centers[b3][0] - basis_centers[b1][0],
				basis_centers[b3][1] - basis_centers[b1][1],
				basis_centers[b3][2] - basis_centers[b1][2] };
			while( dr2[0] < -0.5 ) dr2[0] += 1;
			while( dr2[0] >  0.5 ) dr2[0] -= 1;
			while( dr2[1] < -0.5 ) dr2[1] += 1;
			while( dr2[1] >  0.5 ) dr2[1] -= 1;
			while( dr2[2] < -0.5 ) dr2[2] += 1;
			while( dr2[2] >  0.5 ) dr2[2] -= 1;

		normalize(dr1);
		normalize(dr2);

		cross( dr1, dr2, orientations[b1] );

		normalize( orientations[b1] );				

		// make a prism.

		double ortho[3];

		cross( orientations[b1], dr1, ortho );

		normalize(ortho);

		for( int dz = -1; dz <= 1; dz +=2 )
		for( int ith = 0; ith < 3; ith++ )
		{
			double l = (dr/2)/sin(60*M_PI/180.0);

			double phi = (2*M_PI/6) + ith * (2*M_PI/3);
			all_r[3*npts+0] = basis_centers[b1][0] - orientations[b1][0] * dr/2 * dz + cos(phi) * dr1[0]*l + sin(phi) * ortho[0]*l;
			all_r[3*npts+1] = basis_centers[b1][1] - orientations[b1][1] * dr/2 * dz + cos(phi) * dr1[1]*l + sin(phi) * ortho[1]*l;
			all_r[3*npts+2] = basis_centers[b1][2] - orientations[b1][2] * dr/2 * dz + cos(phi) * dr1[2]*l + sin(phi) * ortho[2]*l;
			npts++;
		}
	}	

	int max_bonds = 20;
	int *bonds = (int *)malloc( sizeof(int) * npts * max_bonds );
	int *nbonds = (int *)malloc( sizeof(int) * npts );
	memset( nbonds, 0, sizeof(int) * npts );

	int n_sub_pts = 6;
	for( int b = 0; b < 8; b++ )
	{
		int sub_bonds[6][3] =
		{
			{ 1, 2, 3 },
			{ 0, 2, 4 },
			{ 0, 1, 5 },
			{ 0, 4, 5 },
			{ 1, 3, 5 },
			{ 2, 3, 4 }
		};

		int b_off = b * n_sub_pts;
		for( int sub = 0; sub < n_sub_pts; sub++ )
		for( int t = 0; t < 3; t++ )
		{
			int p = b_off + sub;

			bonds[p*max_bonds+nbonds[p]] = b_off + sub_bonds[sub][t];
			nbonds[p]++;
		}
	}

	int faces[8][3][4];
	double fcens[8][3][3];

	for( int b = 0; b < 8; b++ )
	{	
		int b_off = b * n_sub_pts;
		int foffs[3][4] = { { 0, 1, 4, 3 },
				    { 1, 2, 5, 4 },
				    { 2, 0, 3, 5 } };

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
	}

	// connect the faces.

	int done[8][3];
	for( int b = 0; b < 8; b++ )
	for( int f = 0; f < 3; f++ )
		done[b][f] = 0;

	int all_done = 0;

	while( !all_done )
	{
		int got_one = 0;
		int best_b1 = -1;
		int best_f1 = -1;
		int o_best_b2 = -1;
		int o_best_f2 = -1;
		double best_or = 1e10;
		for( int b1 = 0; b1 < 8; b1++ )
		{
	//		if( b1 != 0 ) continue;
			for( int f1 = 0; f1 < 3; f1++ )
			{
	//			if( f1 != 2 ) continue;
				if( done[b1][f1] ) continue;
				
				int best_b = -1;	
				int best_f = -1;
	
				double small_r = 1e10;
				for( int xb = 0; xb < 3; xb++ )
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
							while( dr[c] < -0.5 ) dr[c] += 1;
							while( dr[c] >  0.5 ) dr[c] -= 1;
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
					while( dr[c] < -0.5 ) dr[c] += 1;
					while( dr[c] >  0.5 ) dr[c] -= 1;
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
				while( vec[c] < -0.5 ) vec[c] += 1;
				while( vec[c] >  0.5 ) vec[c] -= 1;
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
		else
			break;
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

#if 0
	FILE *qiiPSF = fopen("gyroid.psf","w");
	writePSF(qiiPSF, npts, NULL, all_bonds, all_nbonds );

	printf("%d\n", npts );
	printf("ok\n");
	for( int p = 0; p  < npts; p++ )
		printf("C %lf %lf %lf\n", all_r[3*p+0], all_r[3*p+1], all_r[3*p+2] );
#endif

	double L = 1;
	printf("3d saved surface\n");

	printf("%lf 0.0 0.0\n", L);
	printf("0.0 %lf 0.0\n", L);
	printf("0.0 0.0 %lf\n", L);

	for( int p = 0; p < npts; p++ )
	{
		printf("%d %lf %lf %lf %d", p, all_r[p*3+0], all_r[3*p+1], all_r[3*p+2], nbonds[p] ); 

		for( int b = 0; b < nbonds[p]; b++ )
			printf(" %d", bonds[p*max_bonds+b] );
		printf("\n");
	}
}




