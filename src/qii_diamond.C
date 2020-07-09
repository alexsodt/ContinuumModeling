#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
const char *ops[] = 
{
"x,y,z",
"z+0.25,y+0.25,-x+0.25",
"y+0.25,x+0.25,-z+0.25",
"x+0.25,z+0.25,-y+0.25",
"z+0.25,x+0.25,-y+0.25",
"y+0.25,z+0.25,-x+0.25",
"x+0.25,y+0.25,-z+0.25",
"z+0.25,-y+0.25,x+0.25",
"y+0.25,-x+0.25,z+0.25",
"x+0.25,-z+0.25,y+0.25",
"z+0.25,-x+0.25,y+0.25",
"y+0.25,-z+0.25,x+0.25",
"x+0.25,-y+0.25,z+0.25",
"-z+0.25,y+0.25,x+0.25",
"-y+0.25,x+0.25,z+0.25",
"-x+0.25,z+0.25,y+0.25",
"-z+0.25,x+0.25,y+0.25",
"-y+0.25,z+0.25,x+0.25",
"-x+0.25,y+0.25,z+0.25",
"-z+0.25,-y+0.25,-x+0.25",
"-y+0.25,-x+0.25,-z+0.25",
"-x+0.25,-z+0.25,-y+0.25",
"-z+0.25,-x+0.25,-y+0.25",
"-y+0.25,-z+0.25,-x+0.25",
"-x+0.25,-y+0.25,-z+0.25",
"-z,-y,x",
"-y,-x,z",
"-x,-z,y",
"-z,-x,y",
"-y,-z,x",
"-x,-y,z",
"-z,y,-x",
"-y,x,-z",
"-x,z,-y",
"-z,x,-y",
"-y,z,-x",
"-x,y,-z",
"z,-y,-x",
"y,-x,-z",
"x,-z,-y",
"z,-x,-y",
"y,-z,-x",
"x,-y,-z",
"z,y,x",
"y,x,z",
"x,z,y",
"z,x,y",
"y,z,x",
"z+0.25,y-0.25,-x-0.25",
"y+0.25,x-0.25,-z-0.25",
"x+0.25,z-0.25,-y-0.25",
"z+0.25,x-0.25,-y-0.25",
"y+0.25,z-0.25,-x-0.25",
"x+0.25,y-0.25,-z-0.25",
"z+0.25,-y-0.25,x-0.25",
"y+0.25,-x-0.25,z-0.25",
"x+0.25,-z-0.25,y-0.25",
"z+0.25,-x-0.25,y-0.25",
"y+0.25,-z-0.25,x-0.25",
"x+0.25,-y-0.25,z-0.25",
"-z+0.25,y-0.25,x-0.25",
"-y+0.25,x-0.25,z-0.25",
"-x+0.25,z-0.25,y-0.25",
"-z+0.25,x-0.25,y-0.25",
"-y+0.25,z-0.25,x-0.25",
"-x+0.25,y-0.25,z-0.25",
"-z+0.25,-y-0.25,-x-0.25",
"-y+0.25,-x-0.25,-z-0.25",
"-x+0.25,-z-0.25,-y-0.25",
"-z+0.25,-x-0.25,-y-0.25",
"-y+0.25,-z-0.25,-x-0.25",
"-x+0.25,-y-0.25,-z-0.25",
"-z,-y+0.5,x+0.5",
"-y,-x+0.5,z+0.5",
"-x,-z+0.5,y+0.5",
"-z,-x+0.5,y+0.5",
"-y,-z+0.5,x+0.5",
"-x,-y+0.5,z+0.5",
"-z,y+0.5,-x+0.5",
"-y,x+0.5,-z+0.5",
"-x,z+0.5,-y+0.5",
"-z,x+0.5,-y+0.5",
"-y,z+0.5,-x+0.5",
"-x,y+0.5,-z+0.5",
"z,-y+0.5,-x+0.5",
"y,-x+0.5,-z+0.5",
"x,-z+0.5,-y+0.5",
"z,-x+0.5,-y+0.5",
"y,-z+0.5,-x+0.5",
"x,-y+0.5,-z+0.5",
"z,y+0.5,x+0.5",
"y,x+0.5,z+0.5",
"x,z+0.5,y+0.5",
"z,x+0.5,y+0.5",
"y,z+0.5,x+0.5",
"x,y+0.5,z+0.5",
"z-0.25,y+0.25,-x-0.25",
"y-0.25,x+0.25,-z-0.25",
"x-0.25,z+0.25,-y-0.25",
"z-0.25,x+0.25,-y-0.25",
"y-0.25,z+0.25,-x-0.25",
"x-0.25,y+0.25,-z-0.25",
"z-0.25,-y+0.25,x-0.25",
"y-0.25,-x+0.25,z-0.25",
"x-0.25,-z+0.25,y-0.25",
"z-0.25,-x+0.25,y-0.25",
"y-0.25,-z+0.25,x-0.25",
"x-0.25,-y+0.25,z-0.25",
"-z-0.25,y+0.25,x-0.25",
"-y-0.25,x+0.25,z-0.25",
"-x-0.25,z+0.25,y-0.25",
"-z-0.25,x+0.25,y-0.25",
"-y-0.25,z+0.25,x-0.25",
"-x-0.25,y+0.25,z-0.25",
"-z-0.25,-y+0.25,-x-0.25",
"-y-0.25,-x+0.25,-z-0.25",
"-x-0.25,-z+0.25,-y-0.25",
"-z-0.25,-x+0.25,-y-0.25",
"-y-0.25,-z+0.25,-x-0.25",
"-x-0.25,-y+0.25,-z-0.25",
"-z+0.5,-y,x+0.5",
"-y+0.5,-x,z+0.5",
"-x+0.5,-z,y+0.5",
"-z+0.5,-x,y+0.5",
"-y+0.5,-z,x+0.5",
"-x+0.5,-y,z+0.5",
"-z+0.5,y,-x+0.5",
"-y+0.5,x,-z+0.5",
"-x+0.5,z,-y+0.5",
"-z+0.5,x,-y+0.5",
"-y+0.5,z,-x+0.5",
"-x+0.5,y,-z+0.5",
"z+0.5,-y,-x+0.5",
"y+0.5,-x,-z+0.5",
"x+0.5,-z,-y+0.5",
"z+0.5,-x,-y+0.5",
"y+0.5,-z,-x+0.5",
"x+0.5,-y,-z+0.5",
"z+0.5,y,x+0.5",
"y+0.5,x,z+0.5",
"x+0.5,z,y+0.5",
"z+0.5,x,y+0.5",
"y+0.5,z,x+0.5",
"x+0.5,y,z+0.5",
"z-0.25,y-0.25,-x+0.25",
"y-0.25,x-0.25,-z+0.25",
"x-0.25,z-0.25,-y+0.25",
"z-0.25,x-0.25,-y+0.25",
"y-0.25,z-0.25,-x+0.25",
"x-0.25,y-0.25,-z+0.25",
"z-0.25,-y-0.25,x+0.25",
"y-0.25,-x-0.25,z+0.25",
"x-0.25,-z-0.25,y+0.25",
"z-0.25,-x-0.25,y+0.25",
"y-0.25,-z-0.25,x+0.25",
"x-0.25,-y-0.25,z+0.25",
"-z-0.25,y-0.25,x+0.25",
"-y-0.25,x-0.25,z+0.25",
"-x-0.25,z-0.25,y+0.25",
"-z-0.25,x-0.25,y+0.25",
"-y-0.25,z-0.25,x+0.25",
"-x-0.25,y-0.25,z+0.25",
"-z-0.25,-y-0.25,-x+0.25",
"-y-0.25,-x-0.25,-z+0.25",
"-x-0.25,-z-0.25,-y+0.25",
"-z-0.25,-x-0.25,-y+0.25",
"-y-0.25,-z-0.25,-x+0.25",
"-x-0.25,-y-0.25,-z+0.25",
"-z+0.5,-y+0.5,x",
"-y+0.5,-x+0.5,z",
"-x+0.5,-z+0.5,y",
"-z+0.5,-x+0.5,y",
"-y+0.5,-z+0.5,x",
"-x+0.5,-y+0.5,z",
"-z+0.5,y+0.5,-x",
"-y+0.5,x+0.5,-z",
"-x+0.5,z+0.5,-y",
"-z+0.5,x+0.5,-y",
"-y+0.5,z+0.5,-x",
"-x+0.5,y+0.5,-z",
"z+0.5,-y+0.5,-x",
"y+0.5,-x+0.5,-z",
"x+0.5,-z+0.5,-y",
"z+0.5,-x+0.5,-y",
"y+0.5,-z+0.5,-x",
"x+0.5,-y+0.5,-z",
"z+0.5,y+0.5,x",
"y+0.5,x+0.5,z",
"x+0.5,z+0.5,y",
"z+0.5,x+0.5,y",
"y+0.5,z+0.5,x",
"x+0.5,y+0.5,z"};

double op( const char *the_op, double *src, double *put )
{
	put[0] = 0;	
	put[1] = 0;	
	put[2] = 0;	

	const char *read = the_op;

	for( int c = 0; c < 3; c++ )
	{
		double the_sign = 1;
		if( *read == '-' )
		{
			the_sign = -1;
			read+=1;
		}
		
		switch( *read )
		{
			case 'x':
				put[c] = the_sign * src[0];
				break;

			case 'y':
				put[c] = the_sign * src[1];
				break;

			case 'z':
				put[c] = the_sign * src[2];
				break;
		}

		read += 1;

		if( *read != ',' && *read )
		{
			double shift = atof( read );

			put[c] += shift;

			while( *read && *read != ',' )
			{
				read += 1;
			}
		}

		if( *read ) read += 1;
	}
}

int main( int argc, char **argv )
{
	//creates an idealized Pn3m topology mesh (that must be subdivided).
	// something like icosahedron.

	double bases[6] = { 0.0,0.0,0.0,
			   0.25,0.25,0.25};
	int nbases = 2;
	// each diamond center is a tetrahedron.

	int nspace = 10;
	int nr = 0;
	double *all_r = (double *)malloc( sizeof(double) * 3 * nspace );	
	int nops = sizeof(ops)/sizeof(const char*);
	for( int b = 0; b < nbases; b++ )
	{
		for( int top = 0; top < nops; top++ )
		{
			double put[3];

			op( ops[top], bases+3*b, put );
			int is_dupe = 0;
			for( int x = 0; x < nr; x++ )
			{
				double dr[3] = { all_r[3*x+0] - put[0],
						 all_r[3*x+1] - put[1],
						 all_r[3*x+2] - put[2] };
				while( dr[0] < -0.5 ) dr[0] += 1;
				while( dr[1] < -0.5 ) dr[1] += 1;
				while( dr[2] < -0.5 ) dr[2] += 1;
				while( dr[0] >  0.5 ) dr[0] -= 1;
				while( dr[1] >  0.5 ) dr[1] -= 1;
				while( dr[2] >  0.5 ) dr[2] -= 1;

				if( fabs(dr[0]) > 1e-6 || fabs(dr[1]) > 1e-6 || fabs(dr[2])>1e-6 )
					continue;

				is_dupe = 1;
				break;
			}

			if( !is_dupe )
			{	
				if( nspace == nr )
				{
					nspace *= 2;
					all_r = (double *)realloc( all_r, sizeof(double) * 3 * nspace );
				} 

				all_r[3*nr+0] = put[0];
				all_r[3*nr+1] = put[1];
				all_r[3*nr+2] = put[2];

				nr++;
			}
		}
	}

	// get all bonds.
	
	int MAX_BONDS = 20;
	int *nbonds = (int *)malloc( sizeof(int) * nr );
	memset( nbonds, 0, sizeof(int) * nr );
	int *bonds = (int *)malloc( sizeof(int) * MAX_BONDS * nr );

	double dbas[3] = { bases[3]-bases[0], bases[4]-bases[1], bases[5]-bases[2] };

	double rbase = sqrt(dbas[0]*dbas[0]+dbas[1]*dbas[1]+dbas[2]*dbas[2]);
	double rcut = 1.05 * rbase;

	for( int t1 = 0; t1 < nr; t1++ )
	{
		for( int t2 = 0; t2 < nr; t2++ )
		{
			if( t1 == t2 ) continue;

			double dr[3] = { 
				all_r[3*t1+0] - all_r[3*t2+0],
				all_r[3*t1+1] - all_r[3*t2+1],
				all_r[3*t1+2] - all_r[3*t2+2] };
				
			while( dr[0] < -0.5 ) dr[0] += 1;
			while( dr[1] < -0.5 ) dr[1] += 1;
			while( dr[2] < -0.5 ) dr[2] += 1;
			while( dr[0] >  0.5 ) dr[0] -= 1;
			while( dr[1] >  0.5 ) dr[1] -= 1;
			while( dr[2] >  0.5 ) dr[2] -= 1;
			double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

			if( r < rcut )
			{
				bonds[t1*MAX_BONDS+nbonds[t1]] = t2;
				nbonds[t1] += 1;
			}
		}
	} 	

#if 0
	double dscale = 3.567;
	printf("%d\n", nr );
	printf("Diamond.\n");
	for( int r = 0; r < nr; r++ )
	{
		printf("C %lf %lf %lf", dscale * all_r[3*r+0], dscale * all_r[3*r+1], dscale * all_r[3*r+2] );

		printf(" # %d", nbonds[r] );
		for( int x = 0; x < nbonds[r]; x++ )	
			printf(" %d", bonds[r*MAX_BONDS+x] );
		printf("\n");
	}
#endif

	// each diamond center will be a tetrahedron.
	// 1,1,1 xx -1,1,-1 xx -1,-1,1 xx 1,-1,-1	

	struct tetrahedron
	{	
		double pts[4][3];
		int neighbors[4];
	};
	
	tetrahedron *tets = (struct tetrahedron *)malloc( sizeof( struct tetrahedron) * nr );

	double tet_size = 0.25;// fraction of C-C bond

#if 0
	printf("%d\n", 4 * nr );
	printf("diamond tetrahedrae.\n");
#endif
	for( int r = 0; r < nr; r++ )
	{
		double cen[3];

		cen[0] = all_r[3*r+0];	
		cen[1] = all_r[3*r+1];	
		cen[2] = all_r[3*r+2];	

		double tbonds[4*3];

		for( int b = 0; b < nbonds[r]; b++ )
		{
			tets[r].neighbors[b] = bonds[r*MAX_BONDS+b];
	
			tbonds[b*3+0] = all_r[3*bonds[r*MAX_BONDS+b]+0];		
			tbonds[b*3+1] = all_r[3*bonds[r*MAX_BONDS+b]+1];		
			tbonds[b*3+2] = all_r[3*bonds[r*MAX_BONDS+b]+2];		
		
			double dr[3] = { 
				tbonds[b*3+0] - cen[0],
				tbonds[b*3+1] - cen[1],
				tbonds[b*3+2] - cen[2] };
			
			while( dr[0] < -0.5 ) dr[0] += 1;
			while( dr[1] < -0.5 ) dr[1] += 1;
			while( dr[2] < -0.5 ) dr[2] += 1;
			while( dr[0] >  0.5 ) dr[0] -= 1;
			while( dr[1] >  0.5 ) dr[1] -= 1;
			while( dr[2] >  0.5 ) dr[2] -= 1;
		
			tets[r].pts[b][0] = cen[0] - dr[0] * tet_size; 
			tets[r].pts[b][1] = cen[1] - dr[1] * tet_size; 
			tets[r].pts[b][2] = cen[2] - dr[2] * tet_size; 

#if 0
			printf("C %lf %lf %lf\n", 
				tets[r].pts[b][0],
				tets[r].pts[b][1],
				tets[r].pts[b][2] );
#endif
		}
	}

	// each point is bonded to the other tetrahedral points (3) plus 3 * 2 == 9 neighbors?
	
	printf("3d saved surface\n");

	printf("1.0 0.0 0.0\n");
	printf("0.0 1.0 0.0\n");
	printf("0.0 0.0 1.0\n");

	int ntri_space =10;
	int ntri_excluded_space =10;
	int ntri = 0;
	int ntri_excluded = 0;

	int *all_tri = (int *)malloc( sizeof(int) * 3 * ntri_space );
	int *all_tri_excluded = (int *)malloc( sizeof(int) * 3 * ntri_excluded_space );

	for( int t = 0; t < nr; t++ )
	{
		// get the excluded triangle

		for( int b1 = 0; b1 < 4; b1++ )
		for( int b2 = b1+1; b2 < 4; b2++ )
		for( int b3 = b2+1; b3 < 4; b3++ )
		{
			if( ntri_excluded_space == ntri_excluded )
			{
				ntri_excluded_space *= 2;
				all_tri_excluded = (int *)realloc( all_tri_excluded, sizeof(int) * 3 * ntri_excluded_space );
			}	

			all_tri_excluded[3*ntri_excluded+0] = t*4+b1;
			all_tri_excluded[3*ntri_excluded+1] = t*4+b2; 
			all_tri_excluded[3*ntri_excluded+2] = t*4+b3;  
			ntri_excluded++;
		}
		for( int b = 0; b < 4; b++ )
		{
			double myp[3] = { tets[t].pts[b][0], tets[t].pts[b][1], tets[t].pts[b][2] };
			printf("%d %lf %lf %lf", t*4+b, 
				tets[t].pts[b][0],	
				tets[t].pts[b][1],	
				tets[t].pts[b][2] );	
			int local_bonds[9];
			int nlocal=0;

			for( int z = 0; z < 4; z++ )
			{
				if( z != b ) { local_bonds[nlocal] = t*4+z; nlocal++; }
			}

			double all_neighbors[3*16];
			int ids[16];
			int tet_links[16];
			double dist[16];
			for( int t2x = 0; t2x < 4; t2x++ )
			{
				int t2 = tets[t].neighbors[t2x];

				for( int b2 = 0; b2 < 4; b2++ )
				{
					all_neighbors[t2x*4*3+b2*3+0] = tets[t2].pts[b2][0];
					all_neighbors[t2x*4*3+b2*3+1] = tets[t2].pts[b2][1];
					all_neighbors[t2x*4*3+b2*3+2] = tets[t2].pts[b2][2];
					ids[t2x*4+b2] = t2*4+b2;
					tet_links[t2x*4+b2] = t2;
				}	 
			}	

			int sorter[16];
			for( int n = 0; n < 16; n++ )
			{
				double dr[3] = { 
					all_neighbors[n*3+0] - myp[0], 
					all_neighbors[n*3+1] - myp[1], 
					all_neighbors[n*3+2] - myp[2] }; 

				while( dr[0] < -0.5 ) dr[0] += 1;
				while( dr[1] < -0.5 ) dr[1] += 1;
				while( dr[2] < -0.5 ) dr[2] += 1;
				while( dr[0] >  0.5 ) dr[0] -= 1;
				while( dr[1] >  0.5 ) dr[1] -= 1;
				while( dr[2] >  0.5 ) dr[2] -= 1;
				
				double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

				dist[n] = r;		 
				sorter[n] = n;
			}
			int done = 0;

			while( !done )
			{
				done = 1;
				for( int nx = 0; nx < 15; nx++ )
				{
					if( dist[sorter[nx+1]] < dist[sorter[nx]] )
					{
						int t = sorter[nx];
						sorter[nx] = sorter[nx+1];
						sorter[nx+1] = t;
						done = 0;
					}
				}
			}

			int add_bonds = 6;

			for( int tt = 0; tt < 6; tt++ )
			{
				local_bonds[nlocal] = ids[sorter[tt]];
				nlocal++;
			}
			printf(" %d", nlocal );
			for( int tt = 0; tt < nlocal; tt++ )
				printf(" %d", local_bonds[tt] );
			printf("\n");
			

			
			for( int p = 0; p < add_bonds; p++ )
			for( int p2 = p+1; p2 < add_bonds; p2++ )
			{
				if( tet_links[sorter[p]] == tet_links[sorter[p2]] )
				{
					if( ntri_space == ntri )
					{
						ntri_space *= 2;
						all_tri = (int *)realloc( all_tri, sizeof(int) * 3 * ntri_space );
					}	

					all_tri[3*ntri+0] = t*4+b;
					all_tri[3*ntri+1] = ids[sorter[p]]; 
					all_tri[3*ntri+2] = ids[sorter[p2]]; 
					ntri++;
				}
			}
			
//			for( int t = 0; t < 16; t++ )
//				printf("%d %lf\n", t, dist[sorter[t]] );
	
		}

	}	
	printf("ntri %d\n", ntri );
	for( int t = 0; t < ntri; t++ )
		printf("%d %d %d\n", all_tri[3*t+0], all_tri[3*t+1], all_tri[3*t+2] );
	printf("excl_tri %d\n", ntri_excluded );
	for( int t = 0; t < ntri_excluded; t++ )
		printf("%d %d %d\n", all_tri_excluded[3*t+0], all_tri_excluded[3*t+1], all_tri_excluded[3*t+2] );
}




