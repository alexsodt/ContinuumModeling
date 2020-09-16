#include "face_mask.h"
#include "mutil.h"


void surface_mask::build( surface *theSurface, double *rsurf, double target_area )
{
	int nt = theSurface->nt;

	reg_for_f = (int *)malloc( sizeof(int) * nt );

	int *regions_for_face = reg_for_f;
	int *regions_for_tri = (int *)malloc( sizeof(int) * nt );
	
	double cur_area,area0;
	theSurface->area(rsurf,-1, &cur_area,&area0);

	int nregions = area0 / target_area;

	nregions = nt/4;

	theSurface->getRegions(regions_for_tri, nregions );

	for( int t = 0;  t < nt; t++ )
		regions_for_face[theSurface->theTriangles[t].f] = regions_for_tri[t];

	nreg = nregions;
	nregSpace = nreg;
	
	regions = (surface_region *)malloc( sizeof(surface_region) * nregSpace  );

	for( int r = 0; r < nregions; r++ )
	{
		regions[r].unset = 1;
		regions[r].pool_code = -1;
		regions[r].r_code = r;
		regions[r].com[0] = 0;
		regions[r].com[1] = 0;
		regions[r].com[2] = 0;
		regions[r].is_aligned = 0;
		regions[r].align_x_on[0] = 0;
		regions[r].align_x_on[1] = 0;
		regions[r].flipped = 0;
		regions[r].mod_region = 0;

		int nt_for_r = 0;

		for( int t = 0; t < nt; t++ )
			if( regions_for_face[t] == r )
				nt_for_r++;

		regions[r].f_list = (int *)malloc( sizeof(int) * nt_for_r );
		regions[r].nf = nt_for_r;
		regions[r].nfSpace = nt_for_r;
		nt_for_r = 0;		

		for( int t = 0; t < nt; t++ )
			if( regions_for_face[t] == r )
			{
				regions[r].f_list[nt_for_r] = t; 
				nt_for_r++;
			}
	}
}

void surface_mask::applyPoolToAllRegions( int pool_code )
{
	for( int r = 0; r < nreg; r++ )
	{
		regions[r].pool_code = pool_code;
		regions[r].unset = 0;
	}	
}

void surface_mask::modifyMaskWithPoolAtPoint( surface *theSurface, double *rsurf, double PBC[3][3],
							int pool,
						   int f_in, double u, double v,
							double radius, // in Angstroms.
							double src_xy[2], 
							double src_orientation[2], int flipped )
{
	// move these regions to the beginning so that they are done first.

	// use a different pool structure at this face and surrounding faces.
	// moreover, the src_xy is matched to the u/v point here.
	// the src_orientation is also set.

	// find all the nearby triangles (within a radius) form a new region with this pool structure.	

	double rpt[3], nrm[3];

	theSurface->evaluateRNRM( f_in, u, v, rpt, nrm, rsurf );
	
	if( nreg <= nregSpace )
	{
		nregSpace *= 2;
		regions = (struct surface_region *)realloc( regions, sizeof(struct surface_region) * nregSpace );
	}

	// move all the regions down the list.

	for( int r = nreg; r > 0; r-- )
		regions[r] = regions[r-1];

	for( int t = 0; t < theSurface->nt; t++ )
		reg_for_f[t] = reg_for_f[t]+1; 

	int r = 0;
	
	regions[r].unset = 0;
	regions[r].pool_code = pool;
	regions[r].r_code = r;
	regions[r].com[0] = src_xy[0];
	regions[r].com[1] = src_xy[1];
	regions[r].com[2] = 0;
	regions[r].is_aligned = 1;
	regions[r].align_x_on[0] = src_orientation[0];
	regions[r].align_x_on[1] = src_orientation[1];
	regions[r].f_list = (int *)malloc( sizeof(int) * theSurface->nt );
	regions[r].f_align_on = f_in;
	regions[r].u_align_on = u;
	regions[r].v_align_on = v;
	regions[r].nf = 0;
	regions[r].nfSpace = theSurface->nt;
	regions[r].flipped = flipped;
	regions[r].mod_region = 1;

	for( int f = 0; f < theSurface->nt; f++ )
	{
		double rpt2[3], nrm2[3];

		theSurface->evaluateRNRM( f, 1.0/3.0, 1.0/3.0, rpt2, nrm2, rsurf );
	

		double put[3];
		double alphas[3]={1,1,1};

		double dr[3] = { rpt2[0] - rpt[0], rpt2[1] - rpt[1], rpt2[2] - rpt[2] };

		MinImage3D( dr, PBC, put, alphas );

		double rl = normalize(dr);	

		if( rl < radius || f == f_in )
		{
			// add this face to the list.
		
			// remove from wherever it is.
			
			int from_r = reg_for_f[f];
			int take_index = -1;

			for( int t = 0; t < regions[from_r].nf; t++ )
			{
				if( regions[from_r].f_list[t] == f )
					take_index = t;
			} 
		
			if( take_index == -1 )
			{
				printf("Logical error establishing pool.\n");
				exit(1);
			}

			regions[from_r].f_list[take_index] = regions[from_r].f_list[regions[from_r].nf-1];
			regions[from_r].nf--;

			regions[r].f_list[regions[r].nf] = f;
			regions[r].nf++;				

			reg_for_f[f] = r;			
		}
	}

	regions[r].f_list = (int *)realloc( regions[r].f_list, sizeof(int) * regions[r].nf );
	regions[r].nfSpace = regions[r].nf;



	nreg++;
}

int surface_mask::getPool( int r )
{
	if( r < 0 || r >= nreg )
	{
		printf("Invalid region requested.\n");
		exit(1);	
	}

	return regions[r].pool_code;
}

const struct surface_region *surface_mask::getRegion( int r )
{
	if( r < 0 || r >= nreg )
	{
		printf("Invalid region requested.\n");
		exit(1);	
	}

	return regions+r;
}




