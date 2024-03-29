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
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include "packing.h"

SP monte_carlo_complementary(SPHERE *s, NEIGHBOUR *ne, double g, double height, double diam, double Sxyz, EXTENT ex, int NumberOfSpheres, double he, double Dmin, double t)
{
	int h,p,q,a,nn,mm,m,n,mn;
	SPHERE *gamma;
	double r[3],P[3], deltaE;
	int i,j,k,ii,jj,kk;
	NA output;
	NE *nlist, *temp, *ptrbef, *ptr;
	SP out, spheres;
	double randomnumber;
	double kappa[3]={Dmin,0.5*Dmin,0.25*Dmin},alpha;
	
	for (p=0;p<3;p++)
	{
		alpha=kappa[p];
		for (q=0;q<1100;q++)
		{
			for (gamma=s;gamma != NULL; gamma=gamma->next)
			{
				if (gamma->mask)
				{
					if (gamma->position.z >= he)
					{
						r[2]=2.0*drand48()-1.0;
						P[2]=gamma->position.z+alpha*r[2];
/*						r[0]=alpha*drand48();
						r[1]=2.0*M_PI*drand48();
						r[2]=M_PI*drand48();
						P[2]=gamma->position.z+r[0]*cos(r[2])*/;
						
						if ((P[2]<=(height-0.5*(gamma->diameter)))&&(P[2]>=0.5*(gamma->diameter)))
						{
							r[0]=2.0*drand48()-1.0;
							r[1]=2.0*drand48()-1.0;
							P[0]=gamma->position.x+alpha*r[0];
							P[1]=gamma->position.y+alpha*r[1];
// 								P[0]=gamma->position.x+r[0]*sin(r[2])*cos(r[1]);
// 								P[1]=gamma->position.y+r[0]*sin(r[2])*sin(r[1]);
							
							if (sqrt(pow(P[0],2)+pow(P[1],2))<=0.5*(diam-(gamma->diameter)))
							{
								i=(int)floor((P[0]+0.5*diam)/Sxyz);
								j=(int)floor((P[1]+0.5*diam)/Sxyz);
								k=(int)floor(P[2]/Sxyz);
								m=i+j*(int)ex.Xxtnt+k*(int)ex.Xxtnt*(int)ex.Yxtnt;
							
								output=ncd( ex, i, j, k, ne, gamma->num);
								nn=output.size;
								nlist=output.location;
								mm=1;
								for (temp=nlist; mm&&(temp != NULL); temp=temp->next)
								{
									if (sqrt(pow((temp->address->position.x-P[0]),2)+pow((temp->address->position.y-P[1]),2)+pow((temp->address->position.z-P[2]),2))<0.5*(temp->address->diameter+gamma->diameter))
									{
										mm=0;
									}
								}
								
								chain_free(nlist);
								if (mm != 0)
								{
									deltaE=P[2]-(gamma->position.z);
									randomnumber=drand48();
									if ((deltaE < 0.01)||(exp(-deltaE/t)<=randomnumber))
									{
										ii=(int)floor((gamma->position.x+0.5*diam)/Sxyz);
										jj=(int)floor((gamma->position.y+0.5*diam)/Sxyz);
										kk=(int)floor((gamma->position.z)/Sxyz);
										mn=ii+jj*(int)ex.Xxtnt+kk*(int)ex.Xxtnt*(int)ex.Yxtnt;
										
										if (mn != m)
										{
										ptr=ptrbef=NULL;
										for (nlist=ne[mn].next; nlist != NULL; nlist=nlist->next)
										{
										if (nlist->number < gamma->num)
										{
										ptrbef=nlist;
										}
										else if (nlist->number == gamma->num)
										{
										ptr=nlist;
										break;
										}
										}
										if (ptrbef != NULL)
										ptrbef->next=ptr->next;
										else
										ne[mn].next=ptr->next;
										
										ptrbef=NULL;
										ptrbef=find_place( ne[m].next, ptr);
										if (ptrbef != NULL)
										{
										ptr->next=ptrbef->next;
										ptrbef->next=ptr;
										}
										else
										{
										ptr->next=ne[m].next;
										ne[m].next=ptr;
										}
										}
										
										gamma->position.x=P[0];
										gamma->position.y=P[1];
										gamma->position.z=P[2];
									}
								}
							}
					       }
					}
				}
			}
		}
	}
	out.sph=s;
	out.numsph=NumberOfSpheres;
	return out;
}
