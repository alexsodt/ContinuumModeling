#include <stdlib.h>
#include <string.h>
#include "mutil.h"
#include "aa_build_util.h"
#include "GJK.h"
#include "pdb.h"
#include "globals.h"

void aa_build_data::init( void )
{
	global_nspace = 10;
	global_cycles = (int **)malloc( sizeof(int*) * global_nspace );
	global_cycle_len = (int *)malloc( sizeof(int) * global_nspace );	
	global_ncycles = 0;

	global_atom_space = 10;
	global_bonds_tot = 0;

	// space for the bonds themselves:
	global_bond_space = 10;
	global_bonds = (int *)malloc( sizeof(int) * global_bond_space );

	// atoms need links to bonds.
	global_bond_offsets = (int *)malloc( sizeof(int) * global_atom_space );
	global_nbonds = (int *)malloc( sizeof(int) * global_atom_space );
	memset( global_nbonds, 0, sizeof(int) * global_atom_space );

	nplaced = 0;
	nplaced_pcut = 0;
	nplacedSpace = 100;

	placed_atoms = (double *)malloc( sizeof(double) * 3 * nplacedSpace );	
}

int aa_build_data::curPlace( void )
{
	return nplaced;
}

int aa_build_data::addAtom( double *r )
{
	if( nplaced == nplacedSpace )
	{
		nplacedSpace *= 2;
		placed_atoms = (double *)realloc( placed_atoms, sizeof(double) * nplacedSpace * 3 );		
	}

	placed_atoms[3*nplaced+0] = r[0];
	placed_atoms[3*nplaced+1] = r[1];
	placed_atoms[3*nplaced+2] = r[2];

	boxit( placed_atoms+3*nplaced, nplaced, theBoxes, PBC_vec, nx, ny, nz );

	nplaced++;

	return nplaced-1;
}

// no requirement here for continuity, instead there is a map.

int aa_build_data::checkMappedBonds( double *coords, int *map, int nmapped, int *bonds, int *bond_offsets, int *nbonds )
{
	int clash = 0;
	for( int m = 0; m < nmapped; m++ )
	{
		if( map[m] < 0 ) continue;

		for( int ax = 0; ax < nbonds[m]; ax++ )
		{
			int m2 = bonds[bond_offsets[m]+ax];
	
			if( map[m2] < 0 ) continue;
		
			double *r1 = coords + 3 * map[m];
			double *r2 = coords + 3 * map[m2];

			if( bondClash( r1, r2 ) )
				clash = 1;					
		}
	}	
	return clash;
}

int aa_build_data::checkMappedCycles( double *coords, int *map, int nmapped, int **cycles, int *cycle_lens, int ncycles )
{
	int clash = 0;
	for( int c = 0; c < ncycles && !clash; c++ )
	{
		int relev = 1;

		for( int x = 0; x < cycle_lens[c]; x++ )
		{
			if( map[cycles[c][x]] < 0 )
				relev = 0;
		}

		if( relev )
		{
			int the_cycle[cycle_lens[c]];
			double cycler[3*cycle_lens[c]];

			for( int i = 0; i < cycle_lens[c]; i++ )
				the_cycle[i] = map[cycles[c][i]];

			if( cycleClash( coords, 0, the_cycle, cycle_lens[c] ) )
				clash = 1;
		}
	}
	return clash;
} 

void aa_build_data::addMappedBonds( int offset, int *map, int nmapped, int *bonds, int *bond_offsets, int *nbonds )
{
	int atoms_to_add = 0;
	int bonds_to_add = 0;

	for( int m = 0; m < nmapped; m++ )
	{
		if( map[m] >= 0 )
			atoms_to_add++;
	}

	for( int m = 0; m < nmapped; m++ )
	{
		if( map[m] < 0 ) continue;

		for( int ax = 0; ax < nbonds[m]; ax++ )
		{
			int m2 = bonds[bond_offsets[m]+ax];
	
			if( map[m2] < 0 ) continue;
		
			bonds_to_add++;	
		}
	}
	if( offset + atoms_to_add >= global_atom_space )
	{
		int marker = global_atom_space;
		global_atom_space = offset+atoms_to_add + 1024;

		global_bond_offsets = (int *)realloc( global_bond_offsets, sizeof(int) * global_atom_space );	
		global_nbonds       = (int *)realloc( global_nbonds, sizeof(int) * global_atom_space );	

		memset( global_nbonds+marker, 0, sizeof(int) * (global_atom_space-marker) );
	} 

	if( global_bonds_tot + bonds_to_add >= global_bond_space )
	{
		global_bond_space = global_bonds_tot + bonds_to_add + 1024;

		global_bonds = (int *)realloc( global_bonds, sizeof(int) * global_bond_space );
	}

	for( int m = 0; m < nmapped; m++ )
	{
		if( map[m] < 0 ) continue;

		global_bond_offsets[offset+map[m]] = global_bonds_tot;
		global_nbonds[offset+map[m]] = 0;

		for( int ax = 0; ax < nbonds[m]; ax++ )
		{
			int m2 = bonds[bond_offsets[m]+ax];
	
			if( map[m2] < 0 ) continue;
				
			global_bonds[global_bonds_tot] = offset + map[m2];

			global_bonds_tot++;
			global_nbonds[offset+map[m]] += 1;	
		}
	}	
}

// used for adding a contiguous residue or segment.

void aa_build_data::addBondsInRun( int offset, int a_start, int a_stop, int *bonds, int *bond_offsets, int *nbonds )
{
	int atoms_to_add = a_stop-a_start;
	int bonds_to_add = 0;

	for( int a = a_start; a < a_stop; a++ )
	{
		for( int ax = 0; ax < nbonds[a]; ax++ )
		{
			int m2 = bonds[bond_offsets[a]+ax];
	
			if( m2 < a_start || m2 >= a_stop ) continue;
	
			bonds_to_add++;	
		}
	}

	if( offset + atoms_to_add >= global_atom_space )
	{
		int marker = global_atom_space;

		global_atom_space = offset+atoms_to_add + 1024;

		global_bond_offsets = (int *)realloc( global_bond_offsets, sizeof(int) * global_atom_space );	
		global_nbonds       = (int *)realloc( global_nbonds, sizeof(int) * global_atom_space );	
		
		memset( global_nbonds+marker, 0, sizeof(int) * (global_atom_space-marker) );
	} 

	if( global_bonds_tot + bonds_to_add >= global_bond_space )
	{
		global_bond_space = global_bonds_tot + bonds_to_add + 1024;

		global_bonds = (int *)realloc( global_bonds, sizeof(int) * global_bond_space );
	}

	for( int a = a_start; a < a_stop; a++ )
	{
		global_bond_offsets[offset+a-a_start] = global_bonds_tot;
		global_nbonds[offset+a-a_start] = 0;

		for( int ax = 0; ax < nbonds[a]; ax++ )
		{
			int a2 = bonds[bond_offsets[a]+ax];
	
			if( a2 < a_start || a2 >= a_stop ) continue; 
				
			global_bonds[global_bonds_tot] = offset + a2-a_start;

			global_bonds_tot++;
			global_nbonds[offset+a-a_start] += 1;	
		}
	}	
}


void aa_build_data::addCyclesInRun( int offset, double *coords, int a_start, int a_stop, int **cycles, int *cycle_lens, int ncycles )
{
	int ncycles_to_add = 0;

	for( int c = 0; c < ncycles; c++ )
	{
		int relev = 1;

		for( int x = 0; x < cycle_lens[c]; x++ )
		{
			if( cycles[c][x] < a_start || cycles[c][x] >= a_stop )
				relev = 0;
		}

		if( relev ) ncycles_to_add++;
	}

	if( global_ncycles + ncycles_to_add >= global_nspace )
	{
		global_nspace *= 2;
		global_nspace += ncycles_to_add;
	
		global_cycles = (int **)realloc( global_cycles, sizeof(int *) * global_nspace );
		global_cycle_len = (int *)realloc( global_cycle_len, sizeof(int) * global_nspace );
	}

	for( int c = 0; c < ncycles; c++ )
	{
		int relev = 1;

		for( int x = 0; x < cycle_lens[c]; x++ )
		{
			if( cycles[c][x] < a_start || cycles[c][x] >= a_stop )
				relev = 0;
		}

		if( relev )
		{
			global_cycles[global_ncycles] = (int *)malloc( sizeof(int) * cycle_lens[c] );
			global_cycle_len[global_ncycles] = cycle_lens[c];

			for( int t = 0; t < cycle_lens[c]; t++ )
				global_cycles[global_ncycles][t] = offset + cycles[c][t] - a_start;
			
			// box the cycle.
	
			double box_com[3] = {0,0,0};
			for( int t = 0; t < cycle_lens[c]; t++ )
			{
				box_com[0] += coords[3*(cycles[c][t]-a_start)+0];
				box_com[1] += coords[3*(cycles[c][t]-a_start)+1];
				box_com[2] += coords[3*(cycles[c][t]-a_start)+2];
			}

			box_com[0] /= cycle_lens[c];
			box_com[1] /= cycle_lens[c];
			box_com[2] /= cycle_lens[c];


			boxit( box_com, global_ncycles, cycleBoxes, PBC_vec, nx_c, ny_c, nz_c );
			
			global_ncycles++;
		}
	}

}

void aa_build_data::addMappedCycles( int offset, double *coords, int *map, int nmap, int **cycles, int *cycle_lens, int ncycles )
{
	int ncycles_to_add = 0;

	for( int c = 0; c < ncycles; c++ )
	{
		int relev = 1;

		for( int x = 0; x < cycle_lens[c]; x++ )
		{
			if( map[cycles[c][x]] < 0 )
				relev = 0;
		}

		if( relev ) ncycles_to_add++;
	}

	if( global_ncycles + ncycles_to_add >= global_nspace)
	{
		global_nspace *= 2;
		global_nspace += ncycles_to_add;
	
		global_cycles = (int **)realloc( global_cycles, sizeof(int *) * global_nspace );
		global_cycle_len = (int *)realloc( global_cycle_len, sizeof(int) * global_nspace );
	}

	for( int c = 0; c < ncycles; c++ )
	{
		int relev = 1;

		for( int x = 0; x < cycle_lens[c]; x++ )
		{
			if( map[cycles[c][x]] < 0 )
				relev = 0;
		}

		if( relev )
		{
			global_cycles[global_ncycles] = (int *)malloc( sizeof(int) * cycle_lens[c] );
			global_cycle_len[global_ncycles] = cycle_lens[c];

			for( int t = 0; t < cycle_lens[c]; t++ )
				global_cycles[global_ncycles][t] = offset + map[cycles[c][t]];
			
			// box the cycle.
	
			double box_com[3] = {0,0,0};
			for( int t = 0; t < cycle_lens[c]; t++ )
			{
				box_com[0] += coords[3*(map[cycles[c][t]])+0];
				box_com[1] += coords[3*(map[cycles[c][t]])+1];
				box_com[2] += coords[3*(map[cycles[c][t]])+2];
			}

			box_com[0] /= cycle_lens[c];
			box_com[1] /= cycle_lens[c];
			box_com[2] /= cycle_lens[c];

			boxit( box_com, global_ncycles, cycleBoxes, PBC_vec, nx_c, ny_c, nz_c );
			
			global_ncycles++;
		}
	}

}

void aa_build_data::addSpecialCycle( int *cycle, int len, double *coords )
{
	if( global_ncycles >= global_nspace)
	{
		global_nspace *= 2;
	
		global_cycles = (int **)realloc( global_cycles, sizeof(int *) * global_nspace );
		global_cycle_len = (int *)realloc( global_cycle_len, sizeof(int) * global_nspace );
	}

	global_cycles[global_ncycles] = (int *)malloc( sizeof(int) * len );
	global_cycle_len[global_ncycles] = len;

	for( int t = 0; t < len; t++ )
		global_cycles[global_ncycles][t] = cycle[t];
			
	// box the cycle.
	
	double box_com[3] = {0,0,0};
	for( int t = 0; t < len; t++ )
	{
		box_com[0] += coords[3*t+0];
		box_com[1] += coords[3*t+1];
		box_com[2] += coords[3*t+2];
	}

	box_com[0] /= len;
	box_com[1] /= len;
	box_com[2] /= len;
			
	boxit( box_com, global_ncycles, cycleBoxes, PBC_vec, nx_c, ny_c, nz_c );
	
	global_ncycles++;

}

int aa_build_data::bondClash( double *r1, double *r2 )
{
	int clash = 0;

	// modify this to use boxing data I guess.

	double bond_com[3] = { (r1[0]+r2[0])/2, (r1[1]+r2[1])/2, (r1[2]+r2[2])/2 };

	while( bond_com[0] < 0 ) bond_com[0] += PBC_vec[0][0];
	while( bond_com[1] < 0 ) bond_com[1] += PBC_vec[1][1];
	while( bond_com[2] < 0 ) bond_com[2] += PBC_vec[2][2];
	while( bond_com[0] > PBC_vec[0][0] ) bond_com[0] -= PBC_vec[0][0];
	while( bond_com[1] > PBC_vec[1][1] ) bond_com[1] -= PBC_vec[1][1];
	while( bond_com[2] > PBC_vec[2][2] ) bond_com[2] -= PBC_vec[2][2];
	
	int bx = bond_com[0] * nx_c / PBC_vec[0][0];
	int by = bond_com[1] * ny_c / PBC_vec[1][1];
	int bz = bond_com[2] * nz_c / PBC_vec[2][2];
			
	for( int dx = -1; dx <= 1 && !clash; dx++ )
	for( int dy = -1; dy <= 1 && !clash; dy++ )
	for( int dz = -1; dz <= 1 && !clash; dz++ )
	{
		int n_b_x = bx + dx;
		int n_b_y = by + dy;
		int n_b_z = bz + dz;
	
		if( n_b_x >= nx_c ) n_b_x -= nx_c;
		if( n_b_x < 0 ) n_b_x += nx_c;
		if( n_b_y >= ny_c ) n_b_y -= ny_c;
		if( n_b_y < 0 ) n_b_y += ny_c;
		if( n_b_z >= nz_c ) n_b_z -= nz_c;
		if( n_b_z < 0 ) n_b_z += nz_c;

		int nb = n_b_x*ny_c*nz_c+n_b_y*nz_c+n_b_z;

		for( int px = 0; px < cycleBoxes[nb].np && !clash; px++ )
		{
			int c = cycleBoxes[nb].plist[px];

			double ring_com[3] = {0,0,0};
	
			for( int t = 0; t < global_cycle_len[c]; t++ )
			{
				int loff = global_cycles[c][t];
				ring_com[0] += (placed_atoms)[3*loff+0];
				ring_com[1] += (placed_atoms)[3*loff+1];
				ring_com[2] += (placed_atoms)[3*loff+2];
			}
	
			ring_com[0] /= global_cycle_len[c];
			ring_com[1] /= global_cycle_len[c];
			ring_com[2] /= global_cycle_len[c];
	
			double dr[3];
	
			dr[0] = r1[0] - ring_com[0];	
			dr[1] = r1[1] - ring_com[1];	
			dr[2] = r1[2] - ring_com[2];	
				
			double shift[3] = {0,0,0};
	
			while( dr[0] + shift[0] < -PBC_vec[0][0]/2 ) shift[0] += PBC_vec[0][0]; 
			while( dr[1] + shift[1] < -PBC_vec[1][1]/2 ) shift[1] += PBC_vec[1][1]; 
			while( dr[2] + shift[2] < -PBC_vec[2][2]/2 ) shift[2] += PBC_vec[2][2]; 
			while( dr[0] + shift[0] >  PBC_vec[0][0]/2 ) shift[0] -= PBC_vec[0][0]; 
			while( dr[1] + shift[1] >  PBC_vec[1][1]/2 ) shift[1] -= PBC_vec[1][1]; 
			while( dr[2] + shift[2] >  PBC_vec[2][2]/2 ) shift[2] -= PBC_vec[2][2];
	
			dr[0] += shift[0];
			dr[1] += shift[1];
			dr[2] += shift[2];
	
			double r = normalize(dr);
				
			if( r < 10.0 )
			{
				double r1_s[3] = { r1[0] + shift[0],
						   r1[1] + shift[1],
						   r1[2] + shift[2] };
	
				double r2_s[3] = { r2[0] + shift[0],
						   r2[1] + shift[1],
						   r2[2] + shift[2] };
	
				double convex_set[3*global_cycle_len[c]];
			
				for( int t = 0; t < global_cycle_len[c] && !clash; t++ )
				{
					int loff1 = global_cycles[c][t];
	
					convex_set[3*t+0] = (placed_atoms)[loff1*3+0];
					convex_set[3*t+1] = (placed_atoms)[loff1*3+1];
					convex_set[3*t+2] = (placed_atoms)[loff1*3+2];
					if( t > 0 )
					{
						double dr[3] = { convex_set[3*t+0] - convex_set[3*(t-1)+0],
									convex_set[3*t+1] - convex_set[3*(t-1)+1],
									convex_set[3*t+2] - convex_set[3*(t-1)+2] };
						while( dr[0]  < -PBC_vec[0][0]/2 ) dr[0] += PBC_vec[0][0]; 
						while( dr[1]  < -PBC_vec[1][1]/2 ) dr[1] += PBC_vec[1][1]; 
						while( dr[2]  < -PBC_vec[2][2]/2 ) dr[2] += PBC_vec[2][2]; 
						while( dr[0]  >  PBC_vec[0][0]/2 ) dr[0] -= PBC_vec[0][0]; 
						while( dr[1]  >  PBC_vec[1][1]/2 ) dr[1] -= PBC_vec[1][1]; 
						while( dr[2]  >  PBC_vec[2][2]/2 ) dr[2] -= PBC_vec[2][2];

						convex_set[3*t+0] = convex_set[3*(t-1)+0] + dr[0];
						convex_set[3*t+1] = convex_set[3*(t-1)+1] + dr[1];
						convex_set[3*t+2] = convex_set[3*(t-1)+2] + dr[2];
					}
				}
	
				if( box_GJK( convex_set, global_cycle_len[c], r1_s, r2_s, 0.75 ) )
					clash = 2;
			}
		}
	}

	return clash;
}

int aa_build_data::cycleClash( double *coords, int a_start, int *cycle, int len )
{
	// loop over all the bonds to see what this cycle clashes with
	double ring_com[3] = {0,0,0};

	for( int t = 0; t < len; t++ )
	{
		int loff = cycle[t] - a_start;

		ring_com[0] += coords[3*loff+0];
		ring_com[1] += coords[3*loff+1];
		ring_com[2] += coords[3*loff+2];
	}

	ring_com[0] /= len;
	ring_com[1] /= len;
	ring_com[2] /= len;

	double ring_com_pbc[3] = {ring_com[0], ring_com[1], ring_com[2] };
	
	while( ring_com_pbc[0] < 0 ) ring_com_pbc[0] += PBC_vec[0][0];
	while( ring_com_pbc[1] < 0 ) ring_com_pbc[1] += PBC_vec[1][1];
	while( ring_com_pbc[2] < 0 ) ring_com_pbc[2] += PBC_vec[2][2];
	while( ring_com_pbc[0] >= PBC_vec[0][0] ) ring_com_pbc[0] -= PBC_vec[0][0];
	while( ring_com_pbc[1] >= PBC_vec[1][1] ) ring_com_pbc[1] -= PBC_vec[1][1];
	while( ring_com_pbc[2] >= PBC_vec[2][2] ) ring_com_pbc[2] -= PBC_vec[2][2];

	int clash = 0;
			
	int bx = ring_com_pbc[0] * nx / PBC_vec[0][0];
	int by = ring_com_pbc[1] * ny / PBC_vec[1][1];
	int bz = ring_com_pbc[2] * nz / PBC_vec[2][2];
			
	for( int dx = -1; dx <= 1 && !clash; dx++ )
	for( int dy = -1; dy <= 1 && !clash; dy++ )
	for( int dz = -1; dz <= 1 && !clash; dz++ )
	{
		int n_b_x = bx + dx;
		int n_b_y = by + dy;
		int n_b_z = bz + dz;
	
		if( n_b_x >= nx ) n_b_x -= nx;
		if( n_b_x < 0 ) n_b_x += nx;
		if( n_b_y >= ny ) n_b_y -= ny;
		if( n_b_y < 0 ) n_b_y += ny;
		if( n_b_z >= nz ) n_b_z -= nz;
		if( n_b_z < 0 ) n_b_z += nz;

		int nb = n_b_x*ny*nz+n_b_y*nz+n_b_z;

		for( int px = 0; px < theBoxes[nb].np && !clash; px++ )
		{
			int p = theBoxes[nb].plist[px];

			double dr[3] = { 
				placed_atoms[3*p+0] - ring_com[0], 
				placed_atoms[3*p+1] - ring_com[1],
				placed_atoms[3*p+2] - ring_com[2] };

			double shift[3] = {0,0,0};
			while( dr[0] + shift[0] < -PBC_vec[0][0]/2 ) shift[0] += PBC_vec[0][0]; 
			while( dr[1] + shift[1] < -PBC_vec[1][1]/2 ) shift[1] += PBC_vec[1][1]; 
			while( dr[2] + shift[2] < -PBC_vec[2][2]/2 ) shift[2] += PBC_vec[2][2]; 
			while( dr[0] + shift[0] >  PBC_vec[0][0]/2 ) shift[0] -= PBC_vec[0][0]; 
			while( dr[1] + shift[1] >  PBC_vec[1][1]/2 ) shift[1] -= PBC_vec[1][1]; 
			while( dr[2] + shift[2] >  PBC_vec[2][2]/2 ) shift[2] -= PBC_vec[2][2];

			dr[0] += shift[0];
			dr[1] += shift[1];
			dr[2] += shift[2];
			double dist = normalize(dr);

			if( dist < 7.0 )
			{
				double convex_set[3*len];
			
				for( int t = 0; t < len && !clash; t++ )
				{
					int loff1 = cycle[t] - a_start;

					convex_set[3*t+0] = coords[loff1*3+0];
					convex_set[3*t+1] = coords[loff1*3+1];
					convex_set[3*t+2] = coords[loff1*3+2];
				}
					
				for( int bx = 0; bx < global_nbonds[p] && !clash; bx++ )
				{
					int p2 = global_bonds[global_bond_offsets[p]+bx];

					double r1[3] = { (placed_atoms)[3*p+0], (placed_atoms)[3*p+1],(placed_atoms)[3*p+2]};
					double r2[3] = { (placed_atoms)[3*p2+0], (placed_atoms)[3*p2+1],(placed_atoms)[3*p2+2]};

					r1[0] += shift[0];
					r1[1] += shift[1];
					r1[2] += shift[2];
					
					r2[0] += shift[0];
					r2[1] += shift[1];
					r2[2] += shift[2];

					if( box_GJK( convex_set, len, r1, r2, 0.75 ) )
						clash = 1;
				}
			}
		}
	}

	return clash;
}

int aa_build_data::nclash_aa( double *coords, int lc, int is_mod, double cutoff, int *halve )
{
	int nclash = 0;
	
	int do_halved[3] = {0,0,0};

	if( halve && halve[0] ) do_halved[0]=1;
	if( halve && halve[1] ) do_halved[1]=1;
	if( halve && halve[2] ) do_halved[2]=1;

	for( int tdx = 0; tdx <= do_halved[0]; tdx++ )
	for( int tdy = 0; tdy <= do_halved[1]; tdy++ )
	for( int tdz = 0; tdz <= do_halved[2]; tdz++ )
	{
		for( int t = 0; t < lc; t++ )
		{
			double tr[3] = { coords[3*t+0] + tdx * PBC_vec[0][0]/2, coords[3*t+1] + tdy*PBC_vec[1][1]/2, coords[3*t+2] + tdz * PBC_vec[2][2]/2 };
	
			while( tr[0] < 0 ) tr[0] += PBC_vec[0][0];
			while( tr[1] < 0 ) tr[1] += PBC_vec[1][1];
			while( tr[2] < 0 ) tr[2] += PBC_vec[2][2];
			while( tr[0] >= PBC_vec[0][0] ) tr[0] -= PBC_vec[0][0];
			while( tr[1] >= PBC_vec[1][1] ) tr[1] -= PBC_vec[1][1];
			while( tr[2] >= PBC_vec[2][2] ) tr[2] -= PBC_vec[2][2];
	
			int bx = tr[0] * nx / PBC_vec[0][0];
			int by = tr[1] * ny / PBC_vec[1][1];
			int bz = tr[2] * nz / PBC_vec[2][2];
			
			if( bx >= nx ) bx -= nx;
			if( bx < 0 ) bx += nx;
			if( by >= ny ) by -= ny;
			if( by < 0 ) by += ny;
			if( bz >= nz ) bz -= nz;
			if( bz < 0 ) bz += nz;
		
			int b = bx*ny*nz+by*nz+bz;
			int bad_res = 0;
		
			for( int dx = -1; dx <= 1; dx++ )
			for( int dy = -1; dy <= 1; dy++ )
			for( int dz = -1; dz <= 1; dz++ )
			{
				int n_b_x = bx + dx;
				int n_b_y = by + dy;
				int n_b_z = bz + dz;
			
				if( n_b_x >= nx ) n_b_x -= nx;
				if( n_b_x < 0 ) n_b_x += nx;
				if( n_b_y >= ny ) n_b_y -= ny;
				if( n_b_y < 0 ) n_b_y += ny;
				if( n_b_z >= nz ) n_b_z -= nz;
				if( n_b_z < 0 ) n_b_z += nz;
	
				int nb = n_b_x*ny*nz+n_b_y*nz+n_b_z;
	
	
				for( int px = 0; px < theBoxes[nb].np; px++ )
				{
					int p = theBoxes[nb].plist[px];
	
					if( p < nplaced_pcut && is_mod ) continue;
					if( p >= nplaced_pcut && (tdx > 0 || tdy > 0 || tdz > 0) ) continue;	

					double dr[3] = { 
						placed_atoms[3*p+0] - tr[0], 
						placed_atoms[3*p+1] - tr[1],
						placed_atoms[3*p+2] - tr[2] };
	
					while( dr[0] < -PBC_vec[0][0]/2 ) dr[0] += PBC_vec[0][0]; 
					while( dr[1] < -PBC_vec[1][1]/2 ) dr[1] += PBC_vec[1][1]; 
					while( dr[2] < -PBC_vec[2][2]/2 ) dr[2] += PBC_vec[2][2]; 
					while( dr[0] >  PBC_vec[0][0]/2 ) dr[0] -= PBC_vec[0][0]; 
					while( dr[1] >  PBC_vec[1][1]/2 ) dr[1] -= PBC_vec[1][1]; 
					while( dr[2] >  PBC_vec[2][2]/2 ) dr[2] -= PBC_vec[2][2];
	
					double dist = normalize(dr);

	
					if( dist < cutoff )
					{
						nclash++; 
					}
				}
			}
		}
	}

	return nclash;
}
	
void aa_build_data::setupBoxing( double PBC_in[3][3], int nx_in, int ny_in, int nz_in )
{
	for( int tx = 0; tx < 3; tx++ )
	for( int ty = 0; ty < 3; ty++ )
		PBC_vec[tx][ty] = PBC_in[tx][ty];

	nx = nx_in;
	ny = ny_in;
	nz = nz_in;

	nx_c = nx/3; if( nx_c < 1 ) nx_c = 1;
	ny_c = ny/4; if( ny_c < 1 ) ny_c = 1;
	nz_c = nz/4; if( nz_c < 1 ) nz_c = 1;

	theBoxes = (caa_box *)malloc( sizeof(caa_box) * nx * ny * nz );
	
	for( int b = 0; b < nx*ny*nz; b++ )
	{
		theBoxes[b].np = 0;
		theBoxes[b].npSpace = 2;
		theBoxes[b].plist = (int *)malloc( sizeof(int) * theBoxes[b].npSpace );
	}
	
	cycleBoxes = (caa_box *)malloc( sizeof(caa_box) * nx_c * ny_c * nz_c );
	
	for( int b = 0; b < nx_c*ny_c*nz_c; b++ )
	{
		cycleBoxes[b].np = 0;
		cycleBoxes[b].npSpace = 2;
		cycleBoxes[b].plist = (int *)malloc( sizeof(int) * cycleBoxes[b].npSpace );
	}
		
}


void boxit( double *r_in, int index, caa_box *theBoxes, double PBC_vec[3][3], int nx, int ny, int nz ) 
{ 
	double r[3] = { r_in[0], r_in[1], r_in[2] };
	
	while( r[0] <  0 ) r[0] += PBC_vec[0][0];
	while( r[1] <  0 ) r[1] += PBC_vec[1][1];
	while( r[2] <  0 ) r[2] += PBC_vec[2][2];
	while( r[0] >  PBC_vec[0][0] ) r[0] -= PBC_vec[0][0];
	while( r[1] >  PBC_vec[1][1] ) r[1] -= PBC_vec[1][1];
	while( r[2] >  PBC_vec[2][2] ) r[2] -= PBC_vec[2][2];

	int bx = r[0] * nx / PBC_vec[0][0];
	int by = r[1] * ny / PBC_vec[1][1];
	int bz = r[2] * nz / PBC_vec[2][2];

	if( bx >= nx ) bx -= nx;
	if( bx < 0 ) bx += nx;
	if( by >= ny ) by -= ny;
	if( by < 0 ) by += ny;
	if( bz >= nz ) bz -= nz;
	if( bz < 0 ) bz += nz;

	int b = bx*ny*nz+by*nz+bz;

	if( theBoxes[b].np == theBoxes[b].npSpace )
	{
		theBoxes[b].npSpace *= 2;
		theBoxes[b].plist = (int * )realloc( theBoxes[b].plist, sizeof(int) * theBoxes[b].npSpace );
	}

	theBoxes[b].plist[theBoxes[b].np] = index;
	theBoxes[b].np += 1;
}

void aa_build_data::addCrossedBonds( int a_start, int a_stop ) // notice I did not call this fn hot crossed bonds
{
	// add crossed bonds as cycles.
	// if three bonds cross, and it happens, it can be unbreakable.
	// loop over the atoms we added
	for( int xl = a_start; xl < a_stop; xl++ )
	{	
		int bx = (placed_atoms)[3*xl+0] * nx / PBC_vec[0][0];
		int by = (placed_atoms)[3*xl+1] * ny / PBC_vec[1][1];
		int bz = (placed_atoms)[3*xl+2] * nz / PBC_vec[2][2];
				
		double r1A[3] = { (placed_atoms)[3*xl+0], (placed_atoms)[3*xl+1],(placed_atoms)[3*xl+2]};
				
		// loop over nearby boxes
		for( int dx = -1; dx <= 1 ; dx++ )
		for( int dy = -1; dy <= 1 ; dy++ )
		for( int dz = -1; dz <= 1 ; dz++ )
		{
			int n_b_x = bx + dx;
			int n_b_y = by + dy;
			int n_b_z = bz + dz;
		
			if( n_b_x >= nx ) n_b_x -= nx;
			if( n_b_x < 0 ) n_b_x += nx;
			if( n_b_y >= ny ) n_b_y -= ny;
			if( n_b_y < 0 ) n_b_y += ny;
			if( n_b_z >= nz ) n_b_z -= nz;
			if( n_b_z < 0 ) n_b_z += nz;
		
			int nb = n_b_x*ny*nz+n_b_y*nz+n_b_z;
			// loop over the atoms in the box.	
			for( int px = 0; px < theBoxes[nb].np; px++ )
			{
				int p = theBoxes[nb].plist[px];

				// ignore atoms of our own molecule: intramolecular bonds should all be fine on input.
				if( p >= a_start && p < a_stop ) continue;
	
				double dr[3] = { 
					(placed_atoms)[3*p+0] - r1A[0], 
					(placed_atoms)[3*p+1] - r1A[1],
					(placed_atoms)[3*p+2] - r1A[2] };
		
				double shift[3] = {0,0,0};
				while( dr[0] + shift[0] < -PBC_vec[0][0]/2 ) shift[0] += PBC_vec[0][0]; 
				while( dr[1] + shift[1] < -PBC_vec[1][1]/2 ) shift[1] += PBC_vec[1][1]; 
				while( dr[2] + shift[2] < -PBC_vec[2][2]/2 ) shift[2] += PBC_vec[2][2]; 
				while( dr[0] + shift[0] >  PBC_vec[0][0]/2 ) shift[0] -= PBC_vec[0][0]; 
				while( dr[1] + shift[1] >  PBC_vec[1][1]/2 ) shift[1] -= PBC_vec[1][1]; 
				while( dr[2] + shift[2] >  PBC_vec[2][2]/2 ) shift[2] -= PBC_vec[2][2];
	
				dr[0] += shift[0];
				dr[1] += shift[1];
				dr[2] += shift[2];
				double dist = normalize(dr);
		
			
				// if atom is close enough it's worth checking closely.
				if( dist < 7.0 )
				{
					double r2A[3] = { (placed_atoms)[p*3+0], (placed_atoms)[p*3+1], (placed_atoms)[p*3+2] };

					r2A[0] += shift[0];
					r2A[1] += shift[1];
					r2A[2] += shift[2];

					for( int bx = 0; bx < global_nbonds[xl]; bx++ )
					for( int px = 0; px < global_nbonds[p]; px++ )
					{
						int b2 = global_bonds[global_bond_offsets[xl]+bx];
						int p2 = global_bonds[global_bond_offsets[p]+px];

						if( b2 < xl || p2 < p )
							continue;
	
						// no shift for r1B since it's in r1A's molecule.
						double r1B[3] = { (placed_atoms)[3*b2+0], (placed_atoms)[3*b2+1],(placed_atoms)[3*b2+2]};
						double r2B[3] = { (placed_atoms)[3*p2+0], (placed_atoms)[3*p2+1],(placed_atoms)[3*p2+2]};
	
						// PBC shift determined from first atom of other molecule.
						r2B[0] += shift[0];
						r2B[1] += shift[1];
						r2B[2] += shift[2];
			
						double t1, t2;
						double dr = segmentSegmentDist( r1A, r1B, r2A, r2B, &t1, &t2 );
		
						if( dr < 1.0 && t1 > -0.2 && t1 < 1.2 && t2 > -0.2 && t2 < 1.2 )
						{
							// add a cycle		
							int special_cycle[4] = { xl, p, b2, p2 };		

							double c_cycle[12] =
							{
								(placed_atoms)[3*xl+0], (placed_atoms)[3*xl+1], (placed_atoms)[3*xl+2],
								(placed_atoms)[3*p+0],  (placed_atoms)[3*p+1],  (placed_atoms)[3*p+2],
								(placed_atoms)[3*b2+0], (placed_atoms)[3*b2+1], (placed_atoms)[3*b2+2],
								(placed_atoms)[3*p2+0], (placed_atoms)[3*p2+1], (placed_atoms)[3*p2+2]
							};

							addSpecialCycle( special_cycle, 4, c_cycle );
						}
					}
				}
			}
		}
	}
}

void autoBonds( struct atom_rec *at, int nat, int **bonds_out, int *nbonds_out)
{
	double bond_cutoff = 2.0;

	if( doMartini() )
		bond_cutoff = 6.0;

	int nbondsSpace = 10;
	int nbonds = 0;

	int *bonds = (int *)malloc( sizeof(int) * 2 * nbondsSpace );
	
	for( int a1 = 0; a1 < nat; a1++ )
	{
		for( int a2 = a1+1; a2 < nat; a2++ )
		{
			double dr[3] = { at[a2].x - at[a1].x,
					 at[a2].y - at[a1].y,
					 at[a2].z - at[a1].z };

			double r = normalize(dr);

			if( r < bond_cutoff )
			{
				if( nbonds == nbondsSpace )
				{
					nbondsSpace *= 2;
					bonds = (int *)realloc( bonds, sizeof(int) * 2 * nbondsSpace );
				}
				bonds[2*nbonds+0] = a1;
				bonds[2*nbonds+1] = a2;
				nbonds++;
			}
		}						
	}
	*bonds_out = bonds;
	*nbonds_out = nbonds;
}
