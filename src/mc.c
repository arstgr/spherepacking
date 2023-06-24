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

#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>

#include "packing.h"

SP monte_carlo(double J, SPHERE *s, NEIGHBOUR *ne, double g, double height, double diam, double Sxyz, EXTENT ex, int NumberOfSpheres, int qtop, int ptop, double diamDev1, double diamDev2, double D1, double D2, double Dmin, double temperature, double denuminator, unsigned long int iter)
{
	double kappa[8]={20.0,12.0,7.0,5.0,3.0,2.0,1.0,1.0},alpha;
	double deltaEmax=0.0, tt, ti[4],t;
	int h,p,q,a,nn,mm,m,n,mn;
	SPHERE *gamma;
	double r[3],P[3], deltaE;
	int i,j,k,ii,jj,kk;
	NA output;
	NE *nlist, *temp, *ptrbef, *ptr;
	SP out, spheres;
	double randomnumber;
	
	for (gamma=s; gamma != NULL; gamma=gamma->next)
		deltaEmax += /*(gamma->type)*/pow((gamma->diameter), 3)*(gamma->position.z);
	
	deltaEmax *= 9.81*M_PI*4./3.;
	
	tt=-deltaEmax/(temperature*log(0.99));
	
	for (i=0;i<4;i++)
		ti[i]=denuminator*tt/pow(denuminator,i);
	
	printf("Inner Loop=%d\n",ptop);
	printf("Outer Loop=%d\n",qtop);
	for (a=0;a<4;a++)
	{
		t=ti[a];
		for (h=0;h<8;h++)
		{
			alpha=kappa[h];
			for (p=0;p<qtop;p++)
			{
				for (q=0;q<ptop;q++)
				{
					for (gamma=s;gamma != NULL; gamma=gamma->next)
					{
						if (gamma->mask)
						{
							r[2]=2.0*drand48()-1.0;
							P[2]=gamma->position.z+alpha*r[2];
/*							r[0]=alpha*drand48();
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
//				out=monte_carlo_complementary(s, ne, g, height, diam, Sxyz, ex, NumberOfSpheres, height-2.5*Dmin, Dmin, t*4.);
				analysis(s, NumberOfSpheres, iter, height, diam);
				iter++;
				for (i=0;i<1000;i++)
				{
					spheres= drop( s, J, diam, height, Sxyz, NumberOfSpheres, ex, ne, g, diamDev1, diamDev2, D1, D2);
					s=spheres.sph;
					NumberOfSpheres=spheres.numsph;
				}
			}
		}
	}
	out=monte_carlo_complementary(s, ne, g, height, diam, Sxyz, ex, NumberOfSpheres, height-2.5*Dmin, Dmin, t*4.);
	out=monte_carlo_complementary(s, ne, g, height, diam, Sxyz, ex, NumberOfSpheres, height-1.5*Dmin, Dmin, t*4.);
	out=monte_carlo_complementary(s, ne, g, height, diam, Sxyz, ex, NumberOfSpheres, height-1.*Dmin, Dmin, t*4.);
// 	out.sph=s;
// 	out.numsph=NumberOfSpheres;
	out.iter=iter;
	return out;	
}
