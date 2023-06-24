/***************************************************************************
 *   Copyright (C) $2009$ by $Amirreza Rastegari$                          *
 *   $arstgri@gmail.com$                                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


// #define NTAB 32
// #define diamDev1 10.0
// #define diamDev2 10.0
// #define D1 15.0
// #define D2 10.0

typedef enum type_definition {one,two} SPHERE_TYPE;

typedef struct position_record
{
	double x;
	double y;
	double z;
} POSITION;

typedef struct extent_definition
{
	double Xxtnt;
	double Yxtnt;
	double Zxtnt;
} EXTENT;

typedef struct sphere_poition_record 
{
	int num;
	int mask;
	SPHERE_TYPE type;
	double diameter;
	POSITION position;
	struct sphere_poition_record *next;
} SPHERE;

typedef struct neighbour_list
{
	int number;
	SPHERE *address;
	struct neighbour_list *next;
} NE;

typedef struct neighbour_record
{
	int size;
	NE *next;
} NEIGHBOUR;

typedef struct neighbour_address
{
	NE *location;
	int size;
} NA;

typedef struct sphere_properties
{
	SPHERE *sph;
// 	NEIGHBOUR *nelist;
	int numsph;
	unsigned long int iter;
} SP;

typedef struct glpk_data_prep_input
{
	int iindice;
	int jindice;
	double value;
	struct glpk_data_prep_input *next;
} GL;

typedef struct glpk_temp_constraint_input
{
	double L;
	struct glpk_temp_constraint_input *next;
} LLIST;

typedef struct dataoutput_for_glpk
{
	int a;
	int q;
	int *ia;
	int *ja;
	double *ar;
	double *l;
} OU;

typedef struct moving_particle_list
{
	int number;
	struct moving_particle_list *next;
} MOVING;

typedef struct data_read_and_write
{
	int numsph;
	unsigned long int iter;
	SPHERE *sp;
	NEIGHBOUR *adj;
} DRWS;

typedef struct moving_particle_record
{
	int *k;
	int e;
} MVPR;

NE *find_place(NE *alpha, NE *temp);

void analysis(SPHERE *s, int NumberOfSpheres, unsigned long int iter, double height, double diam);

NA ncd( EXTENT t, int i, int j, int k,NEIGHBOUR *ne, int snum);

SP drop(SPHERE *s, double J, double diam, double height, double Sxyz, int NumberOfSpheres, EXTENT ex, NEIGHBOUR *ne, double g, double diamDev1, double diamDev2, double D1, double D2);

void output_vtk(SPHERE *s, int NumberOfSpheres);
void output_vtkb(SPHERE *s, int NumberOfSpheres);
void output_data(SPHERE *s);

SP monte_carlo(double J, SPHERE *s, NEIGHBOUR *ne, double g, double height, double diam, double Sxyz, EXTENT ex, int NumberOfSpheres, int qtop, int ptop, double diamDev1, double diamDev2, double D1, double D2, double Dmin, double temperature, double denuminator, unsigned long int iter);

SP new_monte_carlo(double J, SPHERE *s, NEIGHBOUR *ne, double g, double height, double diam, double Sxyz, EXTENT ex, int NumberOfSpheres, int qtop, int ptop, double diamDev1, double diamDev2, double D1, double D2, double Dmin, double temperature);

SP test(SPHERE *s, NEIGHBOUR *ne,int NumberOfSpheres, EXTENT ex);

void chain_free(NE *nlist);

OU data_quoe_prep(int *ia, int *ja, double *ar, double *l, SPHERE *s, NEIGHBOUR *ne,int NumberOfSpheres, EXTENT ex,double Sxyz, double gap, double height, double diam);

MOVING find_jammed(int *ia, int *ja, double *ar, double *l, double *load, int q, int a, double gap, int NumberOfSpheres, double deltaMax, double deltaMin);

MVPR unique(MOVING mov);

void fsave_data(SPHERE *s, EXTENT ex, double Sxyz, int NumberOfSpheres, unsigned long int iter);

DRWS fload_data(NEIGHBOUR *adjacent, double diam);

void input_data(double *J, double *diamDev1, double *diamDev2, double *D1, double *D2, double *diam, double *height, double *g, double *gap, double *deltamax, double *deltamin, int *q, int *a, double *temperature, double *denuminator);

SP monte_carlo_complementary(SPHERE *s, NEIGHBOUR *ne, double g, double height, double diam, double Sxyz, EXTENT ex, int NumberOfSpheres, double he, double Dmin, double t);

void output_moving_vtk(SPHERE *s, int NumberOfSpheres, int e, int *k);
