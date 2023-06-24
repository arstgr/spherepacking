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
#include <math.h>
#include <errno.h>
#include "packing.h"


OU data_quoe_prep(int *ia, int *ja, double *ar, double *l, SPHERE *s, NEIGHBOUR *ne,int NumberOfSpheres, EXTENT ex,double Sxyz, double gap, double height, double diam)
{
	int i,j,k,a,q,nn;
	NE *nlist, *tmp;
	SPHERE *temp;
	NA output;
	double den,p,x,y;
	GL *ptr, *eptr, *tptr;
	LLIST *lptr, *leptr, *ltptr;
	OU u;
	
	tmp=nlist=NULL;
	q=a=1;
	
	ptr=eptr=tptr=NULL;
	lptr=leptr=ltptr=NULL;
	
	for (temp=s;temp != NULL; temp=temp->next)
	{
		i=(int)floor((temp->position.x+0.5*diam)/Sxyz);
		j=(int)floor((temp->position.y+0.5*diam)/Sxyz);
		k=(int)floor(temp->position.z/Sxyz);
		
		chain_free(nlist);
		output=ncd( ex, i, j, k, ne, temp->num);
		nn=output.size;
		nlist=output.location;
		
		for (tmp=nlist; tmp != NULL; tmp=tmp->next)
		{
			den=sqrt(pow(tmp->address->position.x - temp->position.x,2)+pow(tmp->address->position.y - temp->position.y,2)+pow(tmp->address->position.z - temp->position.z,2));
			
			if (den < gap*0.5*(tmp->address->diameter + temp->diameter))
			{
				if (tmp->number > temp->num)
				{
					tptr=(GL *)malloc(sizeof(GL));
					if (tptr == NULL)
						error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
					
					tptr->iindice=3*temp->num+1;
					tptr->jindice=a;
					tptr->value=(tmp->address->position.x - temp->position.x)/den;		
					tptr->next=NULL;
					q++;
					
					if (a==1)
						ptr=eptr=tptr;
					else
						eptr->next=tptr;
					eptr=tptr;
					
					tptr=(GL *)malloc(sizeof(GL));
					if (tptr == NULL)
						error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
					
					tptr->iindice=3*temp->num+2;
					tptr->jindice=a;
					tptr->value=(tmp->address->position.y - temp->position.y)/den;
					tptr->next=NULL;
					q++;
					
					eptr->next=tptr;
					eptr=tptr;
					
					tptr=(GL *)malloc(sizeof(GL));
					if (tptr == NULL)
						error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
					
					tptr->iindice=3*temp->num+3;
					tptr->jindice=a;
					tptr->value=(tmp->address->position.z - temp->position.z)/den;
					tptr->next=NULL;
					q++;
					
					eptr->next=tptr;
					eptr=tptr;
					
					tptr=(GL *)malloc(sizeof(GL));
					if (tptr == NULL)
						error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
					
					tptr->iindice=3*tmp->number+1;
					tptr->jindice=a;
					tptr->value=-(tmp->address->position.x - temp->position.x)/den;
					tptr->next=NULL;
					q++;
					
					eptr->next=tptr;
					eptr=tptr;
					
					tptr=(GL *)malloc(sizeof(GL));
					if (tptr == NULL)
						error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
					
					tptr->iindice=3*tmp->number+2;
					tptr->jindice=a;
					tptr->value=-(tmp->address->position.y - temp->position.y)/den;
					tptr->next=NULL;
					q++;
					
					eptr->next=tptr;
					eptr=tptr;
					
					tptr=(GL *)malloc(sizeof(GL));
					if (tptr == NULL)
						error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
					
					tptr->iindice=3*tmp->number+3;
					tptr->jindice=a;
					tptr->value=-(tmp->address->position.z - temp->position.z)/den;
					tptr->next=NULL;
					q++;
					
					eptr->next=tptr;
					eptr=tptr;
					
					ltptr=(LLIST *)malloc(sizeof(LLIST));
					if (ltptr == NULL)
						error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
					
					ltptr->L=den-0.5*(temp->diameter+tmp->address->diameter);
					ltptr->next=NULL;
					
					if (a==1)
						lptr=leptr=ltptr;
					else
						leptr->next=ltptr;
					leptr=ltptr;
					
					a++;
				}
			}
		}
		if (sqrt(pow(temp->position.x,2) + pow(temp->position.y,2)) > (0.5*(diam+temp->diameter)+gap))
		{
			p=atan(temp->position.y/temp->position.x);
			x=0.5*diam*cos(p);
			y=0.5*diam*sin(p);
			den=sqrt(pow(x-temp->position.x,2)+pow(y-temp->position.y,2));
			
			tptr=(GL *)malloc(sizeof(GL));
			if (tptr == NULL)
				error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
			
			tptr->iindice=3*temp->num+1;
			tptr->jindice=a;
			tptr->value=(x- temp->position.x)/den;
			tptr->next=NULL;
			q++;
			
			if (a==1)
				ptr=eptr=tptr;
			else
				eptr->next=tptr;
			eptr=tptr;
			
			tptr=(GL *)malloc(sizeof(GL));
			if (tptr == NULL)
				error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
				
			tptr->iindice=3*temp->num+2;
			tptr->jindice=a;
			tptr->value=(y-temp->position.y)/den;
			tptr->next=NULL;
			q++;
				
			eptr->next=tptr;
			eptr=tptr;
			
			tptr=(GL *)malloc(sizeof(GL));
			if (tptr == NULL)
				error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
			
			tptr->iindice=3*temp->num+3;
			tptr->jindice=a;
			tptr->value=0.0;
			tptr->next=NULL;
			q++;
			
			eptr->next=tptr;
			eptr=tptr;
			
			ltptr=(LLIST *)malloc(sizeof(LLIST));
			if (ltptr == NULL)
				error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
					
			ltptr->L=den-0.5*temp->diameter;
			ltptr->next=NULL;
					
			if (a==1)
				lptr=leptr=ltptr;
			else
				leptr->next=ltptr;
			leptr=ltptr;
			
			a++;
		}
		if (temp->position.z > (height-0.5*temp->diameter - gap))
		{
			tptr=(GL *)malloc(sizeof(GL));
			if (tptr == NULL)
				error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
			
			tptr->iindice=3*temp->num+1;
			tptr->jindice=a;
			tptr->value=0.0;
			tptr->next=NULL;
			q++;
			
			if (a==1)
				ptr=eptr=tptr;
			else
				eptr->next=tptr;
			eptr=tptr;
			
			tptr=(GL *)malloc(sizeof(GL));
			if (tptr == NULL)
				error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
				
			tptr->iindice=3*temp->num+2;
			tptr->jindice=a;
			tptr->value=0.0;
			tptr->next=NULL;
			q++;
				
			eptr->next=tptr;
			eptr=tptr;
			
			tptr=(GL *)malloc(sizeof(GL));
			if (tptr == NULL)
				error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
			
			tptr->iindice=3*temp->num+3;
			tptr->jindice=a;
			tptr->value=1.0;
			tptr->next=NULL;
			q++;
			
			eptr->next=tptr;
			eptr=tptr;
			
			ltptr=(LLIST *)malloc(sizeof(LLIST));
			if (ltptr == NULL)
				error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
					
			ltptr->L=height - temp->position.z - 0.5*temp->diameter;
			ltptr->next=NULL;
					
			if (a==1)
				lptr=leptr=ltptr;
			else
				leptr->next=ltptr;
			leptr=ltptr;
			
			a++;
		}
		if (temp->position.z <= 0.5*temp->diameter+gap)
		{
			tptr=(GL *)malloc(sizeof(GL));
			if (tptr == NULL)
				error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
			
			tptr->iindice=3*temp->num+1;
			tptr->jindice=a;
			tptr->value=0.0;
			tptr->next=NULL;
			q++;
			
			if (a==1)
				ptr=eptr=tptr;
			else
				eptr->next=tptr;
			eptr=tptr;
			
			tptr=(GL *)malloc(sizeof(GL));
			if (tptr == NULL)
				error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
				
			tptr->iindice=3*temp->num+2;
			tptr->jindice=a;
			tptr->value=0.0;
			tptr->next=NULL;
			q++;
				
			eptr->next=tptr;
			eptr=tptr;
			
			tptr=(GL *)malloc(sizeof(GL));
			if (tptr == NULL)
				error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
			
			tptr->iindice=3*temp->num+3;
			tptr->jindice=a;
			tptr->value=-1.0;
			tptr->next=NULL;
			q++;
			
			eptr->next=tptr;
			eptr=tptr;
			
			ltptr=(LLIST *)malloc(sizeof(LLIST));
			if (ltptr == NULL)
				error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
					
			ltptr->L=temp->position.z - 0.5*temp->diameter;
			ltptr->next=NULL;
					
			if (a==1)
				lptr=leptr=ltptr;
			else
				leptr->next=ltptr;
			leptr=ltptr;
			
			a++;
		}
	}
	sp_free(s);
	for (i=0;i<(int)(ex.Xxtnt*ex.Yxtnt*ex.Zxtnt);i++)
		chain_free(ne[i].next);
	
	ia=(int *)calloc(q,sizeof(int));
	if (ia == NULL)
		error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
	ja=(int *)calloc(q,sizeof(int));
	if (ja == NULL)
		error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
	ar=(double *)calloc(q,sizeof(double));
	if (ar == NULL)
		error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
	
	eptr=NULL;
	
	for (tptr = ptr, i=1; tptr != NULL; tptr=tptr->next, i++)
	{
		ia[i]=tptr->iindice;
		ja[i]=tptr->jindice;
		ar[i]=tptr->value;
		
		free(eptr);
		eptr = tptr;
	}
	
	l=(double *)calloc(a,sizeof(double));
	if (l == NULL)
		error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
	
	leptr=NULL;
	
	for (ltptr=lptr, i=1; ltptr != NULL; ltptr=ltptr->next, i++)
	{
		l[i]=ltptr->L;
		free(leptr);
		leptr=ltptr;
	}
	
	u.a=a;
	u.q=q;
	u.ia=ia;
	u.ja=ja;
	u.ar=ar;
	u.l=l;
	
	return u;
}

	
