#ifndef __facemaskh__
#define __facemaskh__

#include "interp.h"

/*
 * This class contains a list of every ``region'' of the surface and its list of triangles
 * Data is associated with those regions that can link to, for example, build the system from a planar patch.
 * 
 *
 *
 * */

struct surface_region 
{
	int unset; // initial value if un-assigned.
	int r_code; // the index into the regions array below for this struct.
	int pool_code; // the code for which to load this lipid patch's molecular structure.
	int mod_region;
	double com[3];
	int f_align_on;
	double u_align_on;
	double v_align_on;
	int is_aligned; // should we align the base to this?
	double align_x_on[2]; // the u,v direction on which to align the x-vector of the base simulation.
	int *f_list; // list of faces for which we use this.
	int nf;
	int nfSpace;
	int flipped; // this means, do we need to flip over the pooled structure to get the right lipids.
		     // we might want the lower leaflet of the pool structure to go on the upper leaflet of the real one.
};

struct surface_mask
{
	int nreg;
	int nregSpace;
	struct surface_region *regions;
	int *reg_for_f;	// which region (offset) for each t.

	// build the region breakdown.
	void build( surface *theSurface, double *rsurf, double target_area );
	void clear( void ); // free memory.
	void applyPoolToAllRegions( int pool_code );
	int getPool(int region);
	const struct surface_region *getRegion( int region );
	void modifyMaskWithPoolAtPoint( surface *theSurface, double *rsurf, double PBC[3][3],
							int pool,
						   int f, double u, double v,
							double radius, // in Angstroms.
							double src_xy[2], 
							double src_orientation[2], int flipped );
};

#endif
