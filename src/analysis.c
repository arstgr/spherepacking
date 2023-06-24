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
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "packing.h"

void analysis(SPHERE *s, int NumberOfSpheres, unsigned long int iter, double height, double diam)
{
	double PF=0.,U=0.0,G=0.;
	SPHERE *gamma,*temp;
	FILE *my,*my1,*my2;
	
	my=fopen("analysis","a");
	my1=fopen("analysis-PF","a");
	my2=fopen("analysis-U","a");
	
	for (gamma=s; gamma != NULL; gamma=gamma->next)
		PF += pow((gamma->diameter), 3)*0.5*0.5*0.5;
	
	for (gamma=s; gamma != NULL; gamma=gamma->next)
		U += /*pow((gamma->diameter), 3)*0.5*0.5*0.5*/(gamma->position.z);
	
	U = U/NumberOfSpheres;

	PF *= (4./3.)/(height*0.25*diam*diam);
	
	
	U *= (4./3.)*M_PI*9.81;
	
	for (gamma=s; gamma != NULL; gamma=gamma->next)
		for (temp=s; temp!=NULL; temp=temp->next)
			G+=sqrt(pow((temp->position.x-gamma->position.x),2)+pow((temp->position.y-gamma->position.y),2)+pow((temp->position.z-gamma->position.z),2));
	
	fprintf(my,"%ld %.6f %.6f %.6f\n",iter,1.-PF,U,G);
	fprintf(my1,"%ld %.6f\n",iter,1.-PF);
	fprintf(my2,"%ld %.6f\n",iter,U);
	fclose(my);
	fclose(my1);
	fclose(my2);
}

	
