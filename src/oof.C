#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "util.h"
#include "mutil.h"

const char *ops[]  = {
"x,y,z",
"x,-y,-z+0.5",
"-x+0.5,y,-z",
"-x,-y+0.5,z",
"z,x,y",
"-z,-x+0.5,y",
"z,-x,-y+0.5",
"-z+0.5,x,-y",
"y,z,x",
"-y+0.5,z,-x",
"-y,-z+0.5,x",
"y,-z,-x+0.5",
"+x+0.25,-z+0.25,+y+0.75",
"+x+0.25,+z+0.75,-y+0.75",
"-x+0.25,-z+0.25,-y+0.25",
"-x+0.75,+z+0.25,+y+0.75",
"+z+0.25,+y+0.75,-x+0.75",
"-z+0.25,+y+0.75,+x+0.25",
"-z+0.25,-y+0.25,-x+0.25",
"+z+0.75,-y+0.75,+x+0.25",
"-y+0.75,+x+0.25,+z+0.75",
"+y+0.75,-x+0.75,+z+0.25",
"-y+0.25,-x+0.25,-z+0.25",
"+y+0.75,+x+0.25,-z+0.25",
"-x,-y,-z",
"-x,y,z+0.5",
"+x+0.5,-y,z",
"x,y+0.5,-z",
"-z,-x,-y",
"z,x+0.5,-y",
"-z,x,y+0.5",
"z+0.5,-x,y",
"-y,-z,-x",
"y+0.5,-z,x",
"y,z+0.5,-x",
"-y,z,x+0.5",
"-x+0.25,+z+0.25,-y+0.75",
"-x+0.25,-z+0.75,+y+0.75",
"+x+0.25,+z+0.25,+y+0.25",
"+x+0.75,-z+0.25,-y+0.75",
"-z+0.25,-y+0.75,+x+0.75",
"+z+0.25,-y+0.75,-x+0.25",
"+z+0.25,+y+0.25,+x+0.25",
"-z+0.75,+y+0.75,-x+0.25",
"+y+0.75,-x+0.25,-z+0.75",
"-y+0.75,+x+0.75,-z+0.25",
"+y+0.25,+x+0.25,+z+0.25",
"-y+0.75,-x+0.25,+z+0.25",
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
		if( *read == '+' )
		{
			the_sign = 1;
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
	double dx = 0.025, dy=0.0125, dz = 0.03;
	double move_away = 0.2;
	double offset[3] = { 0.2, 0 , 0.0 };
	int ndel = 4;
	double basep[(1+ndel*3)*3];

	basep[0] = 0.843;
	basep[1] = 0.054;
	basep[2] = 0.758;

	double p1[3] = { 1.167, 0.394, 0.65 };
	double p2[3] = { 0.608, 0.1, 0.33 };
	double p3[3] = { 0.87, -0.1, 0.93 };

	double dr_axis[3][3] = { { p1[0] - basep[0], p1[1] - basep[1], p1[2] - basep[2] },
			         { p2[0] - basep[0], p2[1] - basep[1], p2[2] - basep[2] },
			         { p3[0] - basep[0], p3[1] - basep[1], p3[2] - basep[2] }};

	basep[0] = 0.75;
	basep[1] = 0;
	basep[2] = 0.75;

	double dr = 0.05;


	
	for( int axis = 0; axis < 3; axis++ )
	{
		for( int del = 0; del < ndel; del++ )
		{
			basep[(1+axis*ndel+del)*3+0] = basep[0] + dr_axis[axis][0] * (1+del)*dr;
			basep[(1+axis*ndel+del)*3+1] = basep[1] + dr_axis[axis][1] * (1+del)*dr;
			basep[(1+axis*ndel+del)*3+2] = basep[2] + dr_axis[axis][2] * (1+del)*dr;
		}
	}

	int nbase = 1+ndel*3;
	// points to 0,1,1
	// 	     1,0,1
	//           1,1,0
	int nspace = 10 + nbase;
	int  npts = nbase;		 
	double *allr = (double *)malloc( sizeof(double) * nspace * 3 );
	int *pt_src = (int *)malloc( sizeof(int) * nspace );

	memset( pt_src, 0, sizeof(int) * nspace );

	int nops = sizeof(ops) / sizeof(const char *);

	memcpy( allr, basep, sizeof(double) * 3 * nbase );

	for( int  b = 0; b < nbase; b++ )
		pt_src[b] = b;
	int done = 0;

	while( !done )
	{ 
		done = 1;
		for( int b = 0; b < npts; b++ )
		{
			for(  int o  = 0; o  < nops; o++ )
			{
				double *rp = allr+3*b;
				double put[3];
				op( ops[o], rp, put );
	
				for( int c = 0; c < 3; c++ )
				{
					while(put[c]>1) put[c]-=1;
					while(put[c]<0) put[c]+=1;
				}
	
				int gotit = 0;
				for( int x = 0; x < npts; x++ )
				{
					double dr[3] = { allr[3*x+0] - put[0], allr[3*x+1] - put[1], allr[3*x+2] - put[2] };		
					for( int c = 0; c < 3; c++ )
					{
						while(dr[c]>0.5) dr[c]-=1;
						while(dr[c]<-0.5) dr[c]+=1;
					}
	
					double rn = normalize(dr);
	
					if( rn <  1e-5 )
					{
						gotit=1;
					}		
				} 
				if( !gotit )
				{
					done = 0;
					if( npts == nspace )
					{
						nspace *= 2;
						allr = (double *)realloc( allr, sizeof(double) * nspace * 3 );
						pt_src = (int *)realloc( pt_src, sizeof(int) * nspace  );
					}
	
					allr[3*npts+0] = put[0];
					allr[3*npts+1] = put[1];
					allr[3*npts+2] = put[2];
	
					pt_src[npts] = pt_src[b];
	
					npts++;
				}
			}
		}
	}
	printf("%d\n", npts );
	printf("great!\n");

	const char *atoms="CNOFHPS";
	for( int x = 0; x < npts; x++ )
		printf("%c %lf %lf %lf\n", atoms[pt_src[x]%7], allr[3*x+0], allr[3*x+1], allr[3*x+2] );
}

/*{
	int nuniq = 2;

	double dsep = 0.3;

	double s = sin(30.0*M_PI/180.0);
	double cc = cos(30.0*M_PI/180.0);

	double L = 2*dsep*cc;

	// four "unique" points.
	double basep[12] = {0,0,0,
		           0,dsep,0, //
			   0,dsep + dsep*s,dsep*cc,
			   dsep*cc,-dsep*s, 0};		
	int main_bonds[4][3] =
	{
		{1, 3, 3 },
		{0, 2, 2 },
		{1, 1, 0 },
		{0, 0, 0 } 	
	};
		
	double pbc_sep[6][3][3] =
	{
		{ {0,0,0}, {0,0,0}, {-1,0,0} },
		{ {0,0,0}, {0,0,0}, {0,0,-1} },
		{ {0,0,0}, {0,0,-1}, { 
	}

	int pts_per_base = 6;
	int nbase = 4;

	double orientations[4][3] = 
	{
		{ 0, 0, 1 },
		{ 1, 0, 0 },
		{ 1, 0, 0 },
		{ 0, 0, 1 }
	};

	double rotors[4] =
	{	// really just for my benefit..
		0, // flat face along y
		M_PI, // flat face along -y,
		0, // flat face along y
		M_PI, // flat face along -y	
	};
	
	double dr = 0.15;
	
	double *all_r = (double *)malloc( sizeof(double) * nbase * pts_per_base * 3 );

	int total_pts = nbase * pts_per_base;
	int max_bonds = 20;	

	int *nbonds = (int *)malloc( sizeof(int) * total_pts  ); 
	int *bonds = (int *)malloc( sizeof(int) * total_pts * max_bonds ); 

	memset( nbonds, 0, sizeof(int) * total_pts );

	for( int b = 0; b < nbase; b++ )
	{
		int sp = 0;

		double *cen = basep  + b*3;
 
		int b_off = b * pts_per_base*3;

		double *axis = orientations[b];
		double face_dir[3] = { 0, cos(rotors[b]), 0 };
		double off_dir[3];
		cross( axis, face_dir, off_dir );

		// base of the prism.
		for( int c = 0; c < 3; c++ )
			all_r[b_off+ 0*3+c] = cen[c] - axis[c]*dr/2 + face_dir[c] * (s/cc) * dr/2 + off_dir[c] * dr/2;  
		for( int c = 0; c < 3; c++ )
			all_r[b_off+ 1*3+c] = cen[c] - axis[c]*dr/2 + face_dir[c] * (s/cc) * dr/2 - off_dir[c] * dr/2;  
		for( int c = 0; c < 3; c++ )
			all_r[b_off+ 2*3+c] = cen[c] - axis[c]*dr/2 - face_dir[c] * (1/cc) * dr/2;  
		
		// apex of the prism
		for( int c = 0; c < 3; c++ )
			all_r[b_off+ 3*3+c] = cen[c] + axis[c]*dr/2 + face_dir[c] * (s/cc) * dr/2 + off_dir[c] * dr/2;  
		for( int c = 0; c < 3; c++ )
			all_r[b_off+ 4*3+c] = cen[c] + axis[c]*dr/2 + face_dir[c] * (s/cc) * dr/2 - off_dir[c] * dr/2;  
		for( int c = 0; c < 3; c++ )
			all_r[b_off+ 5*3+c] = cen[c] + axis[c]*dr/2 - face_dir[c] * (1/cc) * dr/2;  
		
	
		int bonds_base[6][3] = {
			{ 1, 2, 3 },
			{ 0, 2, 4 },
			{ 0, 1, 5 },
			{ 4, 5, 0 },
			{ 3, 5, 1 },
			{ 3, 4, 2 }
		}; 


		for( int p = 0; p < 6; p++ )
		{
			int p1 = b * pts_per_base+p;
			for( int xb = 0; xb < 3; xb++ )
			{
				int p2 = b * pts_per_base+bonds_base[p][xb];

				bonds[p1*max_bonds+nbonds[p1]] = p2;
				nbonds[p1] += 1;
			}
		}
	}

	printf("%d\n", total_pts );
	printf("comment\n");
	
	double scale = 20;
	for( int p = 0; p < total_pts; p++ )
		printf("C %lf %lf %lf\n", all_r[3*p+0]*scale, all_r[3*p+1]*scale, all_r[3*p+2]*scale );





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
#endif
	
	printf("3d saved surface\n");

	printf("%lf 0.0 0.0\n", L);
	printf("0.0 %lf 0.0\n", L);
	printf("0.0 0.0 %lf\n", L);

	for( int p = 0; p < nbase*pts_per_base; p++ )
	{
		printf("%d %lf %lf %lf %d", p, all_r[p*3+0], all_r[3*p+1], all_r[3*p+2], nbonds[p] ); 

		for( int b = 0; b < nbonds[p]; b++ )
			printf(" %d", bonds[p*max_bonds+b] );
		printf("\n");
	}
}*/




