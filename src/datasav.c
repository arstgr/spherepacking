/***************************************************************************
 *   Amirreza Rastegari$                                                   *
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
#include <errno.h>
#include <math.h>
#include "packing.h"

void fsave_data(SPHERE *s, EXTENT ex, double Sxyz, int NumberOfSpheres, unsigned long int iter)
{
	FILE *myfile;
	int i,j;
	SPHERE *temp;
	NE *tmp;
	
	myfile=fopen("sphere.packing","wb");
	fwrite(&NumberOfSpheres,sizeof(int),1,myfile);
	fwrite(&Sxyz,sizeof(double),1,myfile);
	fwrite(&iter,sizeof(unsigned long int),1,myfile);
	fwrite(&ex,sizeof(EXTENT),1,myfile);
	for (temp=s; temp!=NULL; temp=temp->next)
		fwrite(temp,sizeof(SPHERE),1,myfile);
	fclose(myfile);
}

DRWS fload_data(NEIGHBOUR *adjacent, double diam)
{
	SPHERE *s, *temp,*etemp;
	NE *tmp,*etmp, *zeta;
	int i,j,k,numsph;
	int ii,jj,kk,mm;
	FILE *myfile;
	DRWS u;
	EXTENT ex;
	double Sxyz;
	unsigned long int iter;
	
	myfile=fopen("sphere.packing","rb");
	if (myfile==NULL)
		error(EXIT_FAILURE, ENOENT ,"Can not restart the program: No prevoius results were found");
	
	fread(&i,sizeof(int),1,myfile);
	numsph=i;
	fread(&Sxyz,sizeof(double),1,myfile);
	fread(&iter,sizeof(unsigned long int),1,myfile);
	fread(&ex,sizeof(EXTENT),1,myfile);
	
// 	adjacent=(NEIGHBOUR *)calloc((int)(ex.Xxtnt*ex.Yxtnt*ex.Zxtnt),sizeof(NEIGHBOUR));
// 	if (adjacent == NULL)
// 		error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
// 	
// 	for (i=0;i<(int)(ex.Xxtnt*ex.Yxtnt*ex.Zxtnt);i++)
// 		adjacent[i].next=NULL;
	
	for (j=0;j<numsph;j++)
	{
		temp=(SPHERE *)malloc(sizeof(SPHERE));
		if(temp==NULL)
			error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
		fread(temp,sizeof(SPHERE),1,myfile);
		temp->next=NULL;
		if(j==0)
			s=temp;
		else
			etemp->next=temp;
		etemp=temp;
		
		ii=(int)floor((temp->position.x+0.5*diam)/Sxyz);
		jj=(int)floor((temp->position.y+0.5*diam)/Sxyz);
		kk=(int)floor(temp->position.z/Sxyz);
		mm=ii+jj*(int)ex.Xxtnt+kk*(int)ex.Xxtnt*(int)ex.Yxtnt;
		
		tmp=(NE *)malloc(sizeof(NE));
		if (tmp==NULL)
			error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
		
		tmp->address=temp;
		tmp->number=j;
		tmp->next=NULL;
		
		zeta=find_place( adjacent[mm].next, tmp);
		if (zeta != NULL)
		{
			tmp->next=zeta->next;
			zeta->next=tmp;
		}
		else
		{
			tmp->next=adjacent[mm].next;
			adjacent[mm].next=tmp;
		}
		adjacent[mm].size++;
		
	}
	fclose(myfile);

	u.numsph=numsph;
	u.sp=s;
	u.iter=iter;
// 	u.adj=adjacent;
	return u;
}
