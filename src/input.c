/***************************************************************************
 *   $Amirreza Rastegari$                                                  *
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
#include <errno.h>
#include <stdlib.h>
#include "packing.h"

void input_data(double *J, double *diamDev1, double *diamDev2, double *D1, double *D2, double *diam, double *height, double *g, double *gap, double *deltamax, double *deltamin, int *q, int *a, double *temperature, double *denuminator)
{
	int i,j,k;
	char c;
	FILE *myfile;
	
	myfile=fopen("info.packing","r");
	fscanf(myfile,"%*[^=]");
	fscanf(myfile,"%1c",&c);
	fscanf(myfile," %lf",J);
	if ((*J<0.0)||(*J>100.0))
		error(EXIT_FAILURE,0,"Error: J must be between 0-100 ");
	
	fscanf(myfile,"%*[^=]");
	fscanf(myfile,"%1c",&c);
	fscanf(myfile," %lf",diamDev1);
	
	fscanf(myfile,"%*[^=]");
	fscanf(myfile,"%1c",&c);
	fscanf(myfile," %lf",diamDev2);
	
	fscanf(myfile,"%*[^=]");
	fscanf(myfile,"%1c",&c);
	fscanf(myfile," %lf",D1);
	
	fscanf(myfile,"%*[^=]");
	fscanf(myfile,"%1c",&c);
	fscanf(myfile," %lf",D2);
	
	fscanf(myfile,"%*[^=]");
	fscanf(myfile,"%1c",&c);
	fscanf(myfile," %lf",diam);
	
	fscanf(myfile,"%*[^=]");
	fscanf(myfile,"%1c",&c);
	fscanf(myfile," %lf",height);
	
	fscanf(myfile,"%*[^=]");
	fscanf(myfile,"%1c",&c);
	fscanf(myfile," %lf",g);
	
	fscanf(myfile,"%*[^=]");
	fscanf(myfile,"%1c",&c);
	fscanf(myfile," %lf",gap);
	
	fscanf(myfile,"%*[^=]");
	fscanf(myfile,"%1c",&c);
	fscanf(myfile," %lf",deltamin);
	
	fscanf(myfile,"%*[^=]");
	fscanf(myfile,"%1c",&c);
	fscanf(myfile," %lf",deltamax);
	
	fscanf(myfile,"%*[^=]");
	fscanf(myfile,"%1c",&c);
	fscanf(myfile," %d",a);
	
	fscanf(myfile,"%*[^=]");
	fscanf(myfile,"%1c",&c);
	fscanf(myfile," %d",q);

	fscanf(myfile,"%*[^=]");
	fscanf(myfile,"%1c",&c);
	fscanf(myfile," %lf",temperature);

	fscanf(myfile,"%*[^=]");
	fscanf(myfile,"%1c",&c);
	fscanf(myfile," %lf",denuminator);
	
	fclose(myfile);
}
