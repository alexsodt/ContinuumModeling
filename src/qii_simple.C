#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
const char *ops[] = 
{
"x,y,z",
"x+0.5,y,z",
"x,y+0.5,z",
"x,y,z+0.5",
"x+0.5,y+0.5,z",
"x,y+0.5,z+0.5",
"x+0.5,y,z+0.5",
"x+0.5,y+0.5,z+0.5"
};

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
	int nuniq = 2;

	double basep[24] = {0,0,0,
			   0,0,0.5,
			   0,0.5,0,
			   0,0.5,0.5,
			   0.5,0,0,
			   0.5,0.0,0.5,
			   0.5,0.5,0,
			   0.5,0.5,0.5 };
	int nbase = 8;

	double dr = 0.1;
	int pts_per_base = 8;
	
	double *all_r = (double *)malloc( sizeof(double) * nbase * pts_per_base * 3 );

	for( int b = 0; b < nbase; b++ )
	{
		int sp = 0;

		for( int dx = -1; dx <= 1; dx += 2)
		for( int dy = -1; dy <= 1; dy += 2)
		for( int dz = -1; dz <= 1; dz += 2, sp++)
		{
			all_r[b*pts_per_base*3+sp*3+0] = basep[b*3+0] + dx * dr;
			all_r[b*pts_per_base*3+sp*3+1] = basep[b*3+1] + dy * dr;
			all_r[b*pts_per_base*3+sp*3+2] = basep[b*3+2] + dz * dr;
		}
	}

	int total_pts = nbase * pts_per_base;
	int max_bonds = 20;	

	int *nbonds = (int *)malloc( sizeof(int) * total_pts  ); 
	int *bonds = (int *)malloc( sizeof(int) * total_pts * max_bonds ); 

	memset( nbonds, 0, sizeof(int) * total_pts );

	// internal bonds
	
	for( int b1 = 0; b1 < nbase; b1++ )
	{
		for( int subx = 0; subx < 2; subx++ )
		for( int suby = 0; suby < 2; suby++ )
		for( int subz = 0; subz < 2; subz++ )
		{
			// three bonds
			
			int subp1[3] = { 1-subx, suby, subz };
			int subp2[3] = { subx, 1-suby, subz };
			int subp3[3] = { subx, suby, 1-subz };

			int p1 = b1*pts_per_base+subx*4+suby*2+subz;
			int p2 = b1*pts_per_base+subp1[0]*4+subp1[1]*2+subp1[2];
			int p3 = b1*pts_per_base+subp2[0]*4+subp2[1]*2+subp2[2];
			int p4 = b1*pts_per_base+subp3[0]*4+subp3[1]*2+subp3[2];
			
			bonds[p1*max_bonds+nbonds[p1]] = p2; nbonds[p1]++;
			bonds[p1*max_bonds+nbonds[p1]] = p3; nbonds[p1]++;
			bonds[p1*max_bonds+nbonds[p1]] = p4; nbonds[p1]++;
		}	
	}

	for( int b1 = 0; b1 < nbase; b1++ )
	{
		int to_bonds[6][3] = {
			{0,0,-1},
			{0,0,1},
			{0,-1,0},
			{0,1,0},
			{-1,0,0},
			{1,0,0} };
		
		int b_prof[3] = { b1/4, (b1/2)%2, (b1%2) };	

		for( int bond  = 0; bond < 6; bond++ )
		{
			int b_prof2[3] = { b_prof[0] + to_bonds[bond][0],
					   b_prof[1] + to_bonds[bond][1],
					   b_prof[2] + to_bonds[bond][2] };

			if( to_bonds[bond][0] != 0 ) b_prof2[0] = 1-b_prof[0];
			if( to_bonds[bond][1] != 0 ) b_prof2[1] = 1-b_prof[1];
			if( to_bonds[bond][2] != 0 ) b_prof2[2] = 1-b_prof[2];

			int b2 = b_prof2[0] * 4 + b_prof2[1] * 2 + b_prof2[2];

			for( int subx = 0; subx < 2; subx++ )
			for( int suby = 0; suby < 2; suby++ )
			for( int subz = 0; subz < 2; subz++ )
			{
				// is this sub-point bonded ?

				if( to_bonds[bond][0] == 1 && subx == 0 ) continue;
				if( to_bonds[bond][1] == 1 && suby == 0 ) continue;
				if( to_bonds[bond][2] == 1 && subz == 0 ) continue;
				
				if( to_bonds[bond][0] == -1 && subx == 1 ) continue;
				if( to_bonds[bond][1] == -1 && suby == 1 ) continue;
				if( to_bonds[bond][2] == -1 && subz == 1 ) continue;

				int sub_p1 = subx*4+suby*2+subz;

				// bonded to which point?
				
				int subx2 = subx;				
				int suby2 = suby;				
				int subz2 = subz;				
	
				// the direct point, straight across

				if( to_bonds[bond][0] == 1 ) subx2=0;
				if( to_bonds[bond][1] == 1 ) suby2=0;
				if( to_bonds[bond][2] == 1 ) subz2=0;
				
				if( to_bonds[bond][0] == -1  ) subx2 = 1;
				if( to_bonds[bond][1] == -1  ) suby2 = 1;
				if( to_bonds[bond][2] == -1  ) subz2 = 1;

				int sub_p2 = subx2*4+suby2*2+subz2;

				int p1 = b1*pts_per_base+sub_p1;
				int p2 = b2*pts_per_base+sub_p2;


				bonds[p1*max_bonds+nbonds[p1]] = p2;
				nbonds[p1] += 1;
				
				// the right-hand-rule point, diagonal across
				
				int do_x = 0, do_y = 0, do_z = 0;

				if( to_bonds[bond][1] == 0 && to_bonds[bond][2] == 0 ) do_x=1;
				if( to_bonds[bond][0] == 0 && to_bonds[bond][2] == 0 ) do_y=1;
				if( to_bonds[bond][0] == 0 && to_bonds[bond][1] == 0 ) do_z=1;

				int do_neg = 0;
				int map[2]={0,1};
				if( to_bonds[bond][0] + to_bonds[bond][1] + to_bonds[bond][2] < 0 )
				{
					do_neg = 1;
					map[0]=1;
					map[1]=0;
				}

				if( do_x )
				{
					if( !do_neg )
					{
						if( suby == 0 && subz == 0 ) suby2=1; // across in x, y pos
						if( suby == 1 && subz == 0 ) subz2=1; // across in x, y pos
						if( suby == 0 && subz == 1 ) subz2=0; // across in x, y pos
						if( suby == 1 && subz == 1 ) suby2=0; // across in x, y pos
					}
					else
					{
						if( suby == 0 && subz == 0 ) subz2=1; // across in x, y pos
						if( suby == 1 && subz == 0 ) suby2=0; // across in x, y pos
						if( suby == 0 && subz == 1 ) suby2=1; // across in x, y pos
						if( suby == 1 && subz == 1 ) subz2=0; // across in x, y pos
					}
				}
				if( do_y )
				{
					if( !do_neg )
					{
						if( subx == 0 && subz == 0 ) subx2=1; // across in x, y pos
						if( subx == 1 && subz == 0 ) subz2=1; // across in x, y pos
						if( subx == 0 && subz == 1 ) subz2=0; // across in x, y pos
						if( subx == 1 && subz == 1 ) subx2=0; // across in x, y pos
					}
					else
					{
						if( subx == 0 && subz == 0 ) subz2=1; // across in x, y pos
						if( subx == 1 && subz == 0 ) subx2=0; // across in x, y pos
						if( subx == 0 && subz == 1 ) subx2=1; // across in x, y pos
						if( subx == 1 && subz == 1 ) subz2=0; // across in x, y pos
					}
#if 0
					if( !do_neg )
					{
						if( subx == 0 && subz == 0 ) subx2=1; // across in x, y pos
						if( subx == 1 && subz == 0 ) subz2=1; // across in x, y pos
						if( subx == 0 && subz == 1 ) subz2=0; // across in x, y pos
						if( subx == 1 && subz == 1 ) subx2=0; // across in x, y pos
					}
					else
					{
						if( subx == 0 && subz == 0 ) subz2=1; // across in x, y pos
						if( subx == 1 && subz == 0 ) subx2=0; // across in x, y pos
						if( subx == 0 && subz == 1 ) subx2=2; // across in x, y pos
						if( subx == 1 && subz == 1 ) subz2=0; // across in x, y pos
					}
#endif
				}
				if( do_z )
				{
					if( !do_neg )
					{
						if( subx == 0 && suby == 0 ) subx2=1; // across in x, y pos
						if( subx == 1 && suby == 0 ) suby2=1; // across in x, y pos
						if( subx == 0 && suby == 1 ) suby2=0; // across in x, y pos
						if( subx == 1 && suby == 1 ) subx2=0; // across in x, y pos
					}
					else
					{
						if( subx == 0 && suby == 0 ) suby2=1; // across in x, y pos
						if( subx == 1 && suby == 0 ) subx2=0; // across in x, y pos
						if( subx == 0 && suby == 1 ) subx2=1; // across in x, y pos
						if( subx == 1 && suby == 1 ) suby2=0; // across in x, y pos
					}
				}
		

				sub_p2 = subx2*4+suby2*2+subz2;

				p2 = b2*pts_per_base+sub_p2;

				bonds[p1*max_bonds+nbonds[p1]] = p2;
				nbonds[p1] += 1;
			}			
		}	
	}
	
	printf("3d saved surface\n");

	printf("1.0 0.0 0.0\n");
	printf("0.0 1.0 0.0\n");
	printf("0.0 0.0 1.0\n");

	for( int p = 0; p < nbase*pts_per_base; p++ )
	{
		printf("%d %lf %lf %lf %d", p, all_r[p*3+0], all_r[3*p+1], all_r[3*p+2], nbonds[p] ); 

		for( int b = 0; b < nbonds[p]; b++ )
			printf(" %d", bonds[p*max_bonds+b] );
		printf("\n");
	}

}




