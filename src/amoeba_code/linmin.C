#define NRANSI
#include "nrutil.h"
#include <stdio.h>
#define TOL 2.0e-2

int ncom;
double *pcom,*xicom,(*nrfunc)(double []);

void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []))
{
	double brent(double ax, double bx, double cx,
		double (*f)(double), double tol, double *xmin);
	double f1dim(double x);
	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
		double *fc, double (*func)(double));
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;
	pcom=vector(1,n);
	xicom=vector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=0.001;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_vector(xicom,1,n);
	free_vector(pcom,1,n);
}

void linmap(double p[], double xi[], int n, double *fret, double (*func)(double []))
{
	double brent(double ax, double bx, double cx,
		double (*f)(double), double tol, double *xmin);
	double f1dim(double x);
	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
		double *fc, double (*func)(double));
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;
	printf("linmap\n");
	ncom=n;
	pcom=vector(1,n);
	xicom=vector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=0.001;

	for( double del = -100; del <= 100; del += 0.1 )
		printf("map %lf %lf\n", xx*del, f1dim( xx*del ) );	

	printf("end linmap\n");
	free_vector(xicom,1,n);
	free_vector(pcom,1,n);
}
#undef TOL
#undef NRANSI
