#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "interp.h"

surface *getAuxSphere( int nsubdiv )
{
	surface *newSurface = (surface *)malloc( sizeof(surface) );

	newSurface->clear();

	double R = 1.0;
	double phi = (sqrt(5.0) -1)/2.0;
	double a = 1.0 / sqrt(3.0);
	double b = a / phi;
	double c = a * phi;
	
	double vertices[12*3];

	vertices[0] = 0;
	vertices[1] = 0;
	vertices[2] = R;
	vertices[3+0] = 0;
	vertices[3+1] = 0;
	vertices[3+2] = -R;

	double phi1 = atan(0.5);
	double theta = 2*M_PI/5.0;
	double phase = 2*M_PI/10.0;

	for( int i = 0; i < 5; i++ )
	{
		vertices[2*3+i*3+0] = R * sin( i * theta ) * cos( phi1 ); 
		vertices[2*3+i*3+1] = R * cos( i * theta ) * cos( phi1 ); 
		vertices[2*3+i*3+2] = R * sin( phi1 ); 
	}
	for( int i = 0; i < 5; i++ )
	{
		vertices[7*3+i*3+0] = R * sin( i * theta +phase) * cos(phi1); 
		vertices[7*3+i*3+1] = R * cos( i * theta +phase ) * cos(phi1); 
		vertices[7*3+i*3+2] = -R * sin( phi1 ); 
	}

//	printf("hi\n");
//	for( int i = 0; i < 12; i++ )
//		printf("C %lf %lf %lf\n", vertices[3*i+0], vertices[3*i+1], vertices[3*i+2] );

	int neighbors[12][5] = 
	{
		{ 2, 3, 4, 5, 6 },
		{ 7, 8, 9, 10, 11 },
	
		{ 0, 6, 3, 7, 11 }, 
		{ 0, 2, 4, 8, 7 },
		{ 0, 3, 5, 9, 8 },
		{ 0, 4, 6, 10, 9 },
		{ 0, 5, 2, 11, 10 },
		
		{ 1, 8, 11, 2, 3 },
		{ 1, 9,  7, 3, 4 },
		{ 1,10,  8, 4, 5 },
		{ 1,11,  9, 5, 6 },
		{ 1, 7, 10, 6, 2 }
	};

	

	for( int v = 0; v < 12; v++ )
	{
		printf("%d %lf %lf %lf 5",
			v, vertices[3*v+0], vertices[3*v+1], vertices[3*v+2] );
		for( int p = 0; p < 5; p++ )
			printf(" %d", neighbors[v][p] ); 
		printf("\n");
	}

	newSurface->nv = 12;
	newSurface->theVertices = (vertex *)malloc( sizeof(vertex) * newSurface->nv );

	vertex *theVertices = newSurface->theVertices;

	int nv = newSurface->nv;
	
	for( int i = 0; i < 12; i++ )
	{
		theVertices[i].protein_pt = 0;
		theVertices[i].id = i;
		theVertices[i].c0 = 0;
		theVertices[i].valence = 5;
		theVertices[i].r[0] = vertices[3*i+0]; 
		theVertices[i].r[1] = vertices[3*i+1]; 
		theVertices[i].r[2] = vertices[3*i+2]; 

		for( int e = 0; e < 5; e++ )
		{
			theVertices[i].edges[e] = neighbors[i][e]; 
			theVertices[i].edgeIndices[e] = -1;
		}
	}

	newSurface->finishLoad( 60, NULL, 0, NULL );

	for( int sd = 0; sd < nsubdiv; sd++ )
	{
		surface *temp = (surface *)malloc( sizeof(surface) );
		surface *old = newSurface;
		temp->subdivideSurface( newSurface );
		
		newSurface = temp;
	
		old->destroy();
		free(old);	

	}

	newSurface->generatePlan();

	return newSurface;
}
