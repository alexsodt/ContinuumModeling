#ifndef __rd_kaylah__
#define __rd_kaylah__

#include "interp.h"
#include "input.h"

typedef struct
{
	int id;
	double prev_norm;
	double prev_sep;
	int info;
} RD_tracked_info;

typedef struct
{
        int ntracked;
        int ntracked_space;
        RD_tracked_info *tracked_info;
} RD_tracked;

struct RD
{
	RD_tracked **tracked;	

	double k_on;
	double k_off;
	double binding_radius;
	
	void init(int ncomplex);
	void get_tracked(surface *theSurface, double *rsurf, pcomplex **allComplexes, int ncomplex);
};

#endif
