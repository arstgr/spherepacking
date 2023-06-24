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
#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>
#include <math.h>
#include "packing.h"

MOVING find_jammed(int *ia, int *ja, double *ar, double *l, double *load, int q, int a, double gap, int NumberOfSpheres, double deltaMax, double deltaMin)
{
	glp_prob *lp;
	int i,j;
	double z,x,y;
	MOVING *temp, *ptr, *eptr, *tptr, movingout;
	FILE *myfile;
	SPHERE *s;
	DRWS uv;
	
	lp = glp_create_prob();
	glp_set_prob_name(lp, "jamming");
	glp_set_obj_dir(lp, GLP_MAX);
	glp_add_rows(lp, a-1);
	
	for (i=1;i<a;i++)
	{
		if (l[i]>deltaMin)
			glp_set_row_bnds(lp, i, GLP_DB, 0.0, l[i]);
		else
			glp_set_row_bnds(lp, i, GLP_FX, 0.0, 0.0);
	}
	
	
	glp_add_cols(lp, 3*NumberOfSpheres);
	
	
// 	for (i=1;i<=3*NumberOfSpheres;i++)
// 	{
// 		if (load[i-1]>=0)
// 			glp_set_col_bnds(lp, i, GLP_DB, 0.0, deltaMax);
// 		else
// 			glp_set_col_bnds(lp, i, GLP_DB, -deltaMax, 0.0);
// 	}
	
	for (i=1;i<=3*NumberOfSpheres;i++)
		glp_set_col_bnds(lp, i, GLP_DB, -deltaMax, deltaMax);
	
	for (i=1;i<=3*NumberOfSpheres;i++)
		glp_set_obj_coef(lp, i, load[i-1]);
	
	glp_load_matrix(lp, q-1, ja, ia, ar);
	
	glp_simplex(lp, NULL);

	z = glp_get_obj_val(lp);
	printf("maximum possible total displacement=%lf\n",z);
	free(load);
	load=NULL;
	
	j=0;
	ptr=temp=NULL;
	for (i=1;i<=NumberOfSpheres;i++)
	{
		x=glp_get_col_prim(lp, 3*i-2);
		y=glp_get_col_prim(lp, 3*i-1);
		z=glp_get_col_prim(lp, 3*i);
//		if ((fabs(z)>gap)||(fabs(x)>gap)||(fabs(y)>gap))
//		if ((fabs(z)>gap))
		if (sqrt(z*z+x*x+y*y)>gap)
		{
			temp=(MOVING *)malloc(sizeof(MOVING));
			temp->number=i;
			temp->next=NULL;
			
			if (j==0)
				ptr=temp;
			else
				eptr->next=temp;
			eptr=temp;
			j++;
		}	
	}
	myfile=fopen("glpk-solution-report","w");
	for (i=1;i<=NumberOfSpheres;i++)
		fprintf(myfile,"%d %lf %lf %lf\n",i,glp_get_col_prim(lp, 3*i-2),glp_get_col_prim(lp, 3*i-1),glp_get_col_prim(lp, 3*i));
	fclose(myfile);
	
	
	glp_delete_prob(lp);
	movingout.number=j;
	movingout.next=ptr;
	return movingout;
}
