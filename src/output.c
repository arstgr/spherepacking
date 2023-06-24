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
#include <string.h>
#include "packing.h"

void output_vtk(SPHERE *s, int NumberOfSpheres)
{
	int i,j,k;
	FILE *myfile;
	SPHERE *temp;
	
	myfile=fopen("sphere.vtk","w");
	fprintf(myfile, "# vtk DataFile Version 3.0 \n");
	fprintf(myfile, "vtk output \n");
	fprintf(myfile, "ASCII \n");
	fprintf(myfile, "DATASET POLYDATA \n");
	fprintf(myfile, "POINTS  %d float \n", NumberOfSpheres);
	
	
	for (i=0, temp=s;i<NumberOfSpheres;temp=temp->next, i++)
		fprintf(myfile, "%lf %lf %lf\n", temp->position.x, temp->position.z, temp->position.y);
	fprintf(myfile, "POINT_DATA %d\n", NumberOfSpheres);
	fprintf(myfile, "SCALARS diameters float \n");
	fprintf(myfile, "LOOKUP_TABLE myTable \n");
	
	for (temp=s;temp!=NULL;temp=temp->next)
		fprintf(myfile,"%lf\n",temp->diameter);
	
	fprintf(myfile, "\n\n\n");
	fclose(myfile);
}

void output_moving_vtk(SPHERE *s, int NumberOfSpheres, int e, int *k)
{
	int i,j;
	FILE *myfile;
	SPHERE *temp;

	myfile=fopen("moving-sphere.vtk","w");
	fprintf(myfile, "# vtk DataFile Version 3.0 \n");
	fprintf(myfile, "vtk output \n");
	fprintf(myfile, "ASCII \n");
	fprintf(myfile, "DATASET POLYDATA \n");
	fprintf(myfile, "POINTS  %d float \n", e);

	for (i=0;i<e;i++)
		for (temp=s;temp!=NULL;temp=temp->next)
			if (temp->num==(k[i]-1))
				fprintf(myfile, "%lf %lf %lf\n", temp->position.x, temp->position.z, temp->position.y);

	fprintf(myfile, "POINT_DATA %d\n", e);
	fprintf(myfile, "SCALARS diameters float \n");
	fprintf(myfile, "LOOKUP_TABLE myTable \n");
	
	for (i=0;i<e;i++)
		for (temp=s;temp!=NULL;temp=temp->next)
			if (temp->num==(k[i]-1))
				fprintf(myfile,"%lf\n",temp->diameter);
	
	fprintf(myfile, "\n\n\n");
	fclose(myfile);
}

void output_vtkb(SPHERE *s, int NumberOfSpheres)
{
	int i,j,k;
	FILE *myfile;
	SPHERE *temp;
	
	myfile=fopen("sphere.vtk","wb");
	fprintf(myfile, "# vtk DataFile Version 3.0 \n");
	fprintf(myfile, "vtk output \n");
	fprintf(myfile, "BINARY \n");
	fprintf(myfile, "DATASET POLYDATA \n");
	fprintf(myfile, "POINTS  %d float \n", NumberOfSpheres);
	
	
	for (i=0, temp=s;i<NumberOfSpheres;temp=temp->next, i++)
	{
		fwrite(&temp->position.x, sizeof(double),1,myfile);
		fwrite(&temp->position.y, sizeof(double),1,myfile);
		fwrite(&temp->position.z, sizeof(double),1,myfile);
	}
	
	fprintf(myfile, "POINT_DATA %d\n", NumberOfSpheres);
	fprintf(myfile, "SCALARS diameters float \n");
	fprintf(myfile, "LOOKUP_TABLE myTable \n");
	
	for (temp=s;temp!=NULL;temp=temp->next)
		fwrite(&temp->diameter,sizeof(double),1,myfile);
	
	fclose(myfile);
}

void output_data(SPHERE *s)
{
	int i,j;
	FILE *myfile,*myotherfile;
	SPHERE *temp;
	
	myfile=fopen("C.packing","w");
	myotherfile=fopen("d.packing","w");
	for (temp=s;temp!=NULL;temp=temp->next)
	{
		fprintf(myfile,"%lf %lf %lf\n",temp->position.x,temp->position.y,temp->position.z);
		fprintf(myotherfile,"%lf\n",temp->diameter);
	}
	fclose(myfile);
	fclose(myotherfile);
}

	

	
