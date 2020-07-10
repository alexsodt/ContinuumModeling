#ifndef __maxch__
#define __maxch__

#include "interp.h"

void max_c( surface *theSurface, int *f, double *u, double *v, int niter, double *rsurf, int mode );

#define MAX_C_GAUSS 		0
#define MAX_C_POS		1
#define MAX_C_NEG		2
#define MAX_C_TOT_POS		3
#define MAX_C_TOT_NEG		4
#define MAX_C_MIN		5


#endif
