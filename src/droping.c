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
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "packing.h"

SP drop(SPHERE *s, double J, double diam, double height, double Sxyz, int NumberOfSpheres, EXTENT ex, NEIGHBOUR *ne, double g, double diamDev1, double diamDev2, double D1, double D2)
{
	SPHERE temp, *beta, *gamma;
	NE *nlist, *tmp, *chai, *zeta;
	int i,j,k,kk,nn,mm;
	int u=1,h;
	double t, test=drand48();
	const double rad=0.5*diam;
	double theta, radious;
	double coef;
	NA output;
	SP spheres;
	double distance,delta;
	while (u)
	{
		t=100.0*drand48();
		if (t<J)
		{
			temp.diameter = diamDev1*drand48()+D1;
			temp.type = one;
		}
		else
		{
			temp.diameter = diamDev2*drand48()+D2;
			temp.type = two;
		}
		h=1;
		while (h)
		{
			radious = (rad-0.5*temp.diameter)*drand48();
			theta = 2.0*M_PI*drand48();
			temp.position.x = radious*cos(theta);
			temp.position.y = radious*sin(theta);
// 			temp.position.x = diam*(drand48()-0.5);
// 			temp.position.y = diam*(drand48()-0.5);	
			
			if ((sqrt((temp.position.x*temp.position.x)+(temp.position.y*temp.position.y))+0.5*temp.diameter)<=rad)
				h=0;
		}

		temp.num=NumberOfSpheres;
		temp.mask=1;
		
		i=(int)floor((temp.position.x+0.5*diam)/Sxyz);
		j=(int)floor((temp.position.y+0.5*diam)/Sxyz);
		k=(int)floor((height-0.0001)/Sxyz);

		nlist=NULL;
		output=ncd( ex, i, j, k, ne, -1);
		nn=output.size;
		nlist=output.location;
		
		coef=1.0;
		for (temp.position.z=height;coef*(temp.position.z) > 0.5*temp.diameter; temp.position.z -= 0.001*(height+1.0-temp.diameter))
		{
			kk=(int)floor(temp.position.z/Sxyz);
			if (kk<k)
			{
				k=kk;
				chain_free(nlist);
				nlist=NULL;
				output=ncd( ex, i, j, k, ne, -1);
				nn=output.size;
				nlist=output.location;
			}
			if (nn)
			{
				for (mm=0,tmp=nlist;mm<nn;tmp=tmp->next,mm++)
				{
					distance=sqrt(pow(tmp->address->position.x-temp.position.x,2)+pow(tmp->address->position.y-temp.position.y,2)+pow(tmp->address->position.z-temp.position.z,2));
					delta=0.5*(temp.diameter+tmp->address->diameter);
					if (distance<(1.+g)*delta)
					{
						coef=0.0;
						mm=nn;
					}
				}
			}
		}
		k=(int)floor(temp.position.z/Sxyz);
		chain_free(nlist);
		if ((temp.position.z<=(height-0.5*temp.diameter))/*&&(distance>=delta)*/)
		{
			beta=(SPHERE *)malloc(sizeof(SPHERE));
			if (beta==NULL)
				error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
			*beta=temp;
			beta->next = NULL;
			if (NumberOfSpheres==0)
				s=beta;
			else
			{
				for (h=0,gamma=s;h<(NumberOfSpheres-1);h++)
					gamma=gamma->next;
				gamma->next=beta;
			}
			mm=i+j*(int)ex.Xxtnt+k*(int)ex.Xxtnt*(int)ex.Yxtnt;
			chai=(NE *)malloc(sizeof(NE));
			if (chai==NULL)
				error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
			chai->address=beta;
			chai->number=NumberOfSpheres;
			chai->next=NULL;
			
			zeta=find_place( ne[mm].next, chai);
			
			if (zeta != NULL)
			{
				chai->next=zeta->next;
				zeta->next=chai;
			}
			else
			{
				chai->next=ne[mm].next;
				ne[mm].next=chai;
			}
			NumberOfSpheres++;
			ne[mm].size++;
		}
		else
		{
			u=0;	
		}
	}
	spheres.numsph=NumberOfSpheres;
	spheres.sph=s;
	return spheres;
}

main(int argc, char *argv[])
{
	double Dmax=drand48()*drand48();
	double Dmin;
	double Sxyz,xx,yy,zz,diamDev1,diamDev2,D1,D2;
	EXTENT xtnt;
	double g=20.0,gap=40.0;
	double diam=90.0,height=200.0;
	NEIGHBOUR *adjacent;
	SPHERE *s,*gamma;
	double J=40.0;
	int NumberOfSpheres=0;
	int i,j,aaa,qqq;
	SP spheres;
	NE *nlist;
	double *ar, *l;
	int *ia, *ja, *k,e, *cs;
	OU u;
	double deltaMin=1.,deltaMax=60.0;
	double *load, temperature;
	MOVING movingout;
	double deltaEmax=0.,denuminator;
	DRWS uv;
	MVPR mv;
	unsigned long int iter;

	input_data(&J, &diamDev1, &diamDev2, &D1, &D2, &diam, &height, &g, &gap, &deltaMax, &deltaMin, &qqq, &aaa, &temperature, &denuminator);
	iter=0;
	if ((D1*D2) != 0.0)
	{
		Dmax=((diamDev1+D1)>=(diamDev2+D2) ? (diamDev1+D1) : (diamDev2+D2));
		Dmin=((D2>=D1) ? D1 : D2);
	}
	else
	{
		if (D1 != 0.0)
		{
			Dmax=D1+diamDev1;
			Dmin=D1;
		}
		else
		{
			Dmax=D2+diamDev2;
			Dmin=D2;
		}
	}
	if ((D1==0.0)&&(D2==0.0))
		error(EXIT_FAILURE, 0 ,"Bad Initial Conditions: Problem is illposed !");
	if ((J==0.0)&&(D2==0.0))
		error(EXIT_FAILURE, 0 ,"Bad Initial Conditions: Problem is illposed !");
	if ((J==100.0)&&(D1==0.0))
		error(EXIT_FAILURE, 0 ,"Bad Initial Conditions: Problem is illposed !");
	if (Dmin > diam)
		error(EXIT_FAILURE, 0 ,"Bad Initial Conditions: Problem is illposed !");
	if (Dmin > height)
		error(EXIT_FAILURE, 0 ,"Bad Initial Conditions: Problem is illposed !");
// 	if ((D1*D2 == 0.)&&((J != 0.)||(J != 100.)))
// 		error(EXIT_FAILURE, 0 ,"Bad Initial Conditions: Problem is illposed !");
	
	printf("J=%lf  diamDev1=%lf   diamDev2=%lf   D1=%lf  D2=%lf   diam=%lf   height=%lf   g=%lf  gap=%lf   deltaMax=%lf  deltaMin=%lf   qqq=%d  aaa=%d    Temperature=%lf  M Factor=%lf\n",J,diamDev1,diamDev2,D1,D2,diam,height,g,gap,deltaMax,deltaMin,qqq,aaa,temperature, denuminator);
	
	if (Dmin == height)
		height=height + 5;
	if (((height-Dmin)<=3)&&(height-Dmin)>=0)
		height=height + 4;
// 	Sxyz=ceil(0.5*(Dmin+Dmax)+0.0001);
	Sxyz=ceil(Dmax);
	printf("Sxyz=%lf\n",Sxyz);
	xtnt.Xxtnt=xtnt.Yxtnt=ceil(diam/Sxyz);

	xtnt.Zxtnt=ceil(height/Sxyz);
	printf("Xbox=%f  Ybox=%f   Zbox=%f\n",xtnt.Xxtnt,xtnt.Yxtnt,xtnt.Zxtnt);
	adjacent=(NEIGHBOUR *)calloc((int)(xtnt.Xxtnt*xtnt.Yxtnt*xtnt.Zxtnt),sizeof(NEIGHBOUR));
	if (adjacent == NULL)
		error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");

	for (i=0;i<(int)(xtnt.Xxtnt*xtnt.Yxtnt*xtnt.Zxtnt);i++)
		adjacent[i].next=NULL;

	printf("sizeof box=%d\n",(int)(xtnt.Xxtnt*xtnt.Yxtnt*xtnt.Zxtnt));
	if (argc == 1)
	{
		remove("analysis");
		remove("analysis-PF");
		remove("analysis-U");
		printf("Dropping Started\n");
		for (i=0;i<1000;i++)
		{
			spheres= drop( s, J, diam, height, Sxyz, NumberOfSpheres, xtnt, adjacent, 0.0001, diamDev1, diamDev2, D1, D2);
			s=spheres.sph;
			NumberOfSpheres=spheres.numsph;
		}
	}
	else if ((*argv[1] == '-')&&(*(argv[1]+1) == 'r'))
	{
		uv=fload_data(adjacent, diam);
		s=uv.sp;
		NumberOfSpheres=uv.numsph;
		iter=uv.iter;
		printf("Iteration %ld\n",iter);
	}
	else
		error(EXIT_FAILURE, -1 ,"Bad input options");
	printf("NumberOfSpheres=%d\n",NumberOfSpheres);
	if (NumberOfSpheres==0)
		error(EXIT_FAILURE, -1 ,"Bad Initial Conditions: Problem is illposed !");
	
	
	printf("Monte Carlo Iterations Have been Started\n");
	spheres=monte_carlo( J, s, adjacent, 0.0001, height, diam, Sxyz, xtnt, NumberOfSpheres, qqq, aaa,diamDev1, diamDev2, D1, D2, Dmin, temperature, denuminator,iter);
	s=spheres.sph;
	NumberOfSpheres=spheres.numsph;
	printf("NumberOfSpheres=%d\n",NumberOfSpheres);
	iter=spheres.iter;
	fsave_data(s, xtnt, Sxyz, NumberOfSpheres,iter);
	
// 	printf("New Monte Carlo Iterations Have been Started\n");
// 	spheres=new_monte_carlo( J, s, adjacent, 0.0001, height, diam, Sxyz, xtnt, NumberOfSpheres, qqq, aaa,diamDev1, diamDev2, D1, D2, Dmin, temperature);
// 	s=spheres.sph;
// 	NumberOfSpheres=spheres.numsph;
// 	printf("NumberOfSpheres=%d\n",NumberOfSpheres);

	for (gamma=s; gamma != NULL; gamma=gamma->next)
		deltaEmax += pow((gamma->diameter), 3)*0.5*0.5*0.5;

	deltaEmax *= (4./3.)/(height*0.25*diam*diam);
	
	printf("Packing Factor of this system = %lf\n",deltaEmax);
	deltaEmax=0.;
	for (gamma=s; gamma != NULL; gamma=gamma->next)
		deltaEmax += pow((gamma->diameter), 3)*0.5*0.5*0.5*(gamma->position.z);
	
	deltaEmax *= (4./3.)*M_PI*9.81;
	printf("Potential Energy of this system = %lf\n",deltaEmax);

	output_data(s);
	output_vtk( s, NumberOfSpheres);
// 	g=10.0;
	
	u=data_quoe_prep(ia, ja, ar, l, s, &adjacent[0], NumberOfSpheres, xtnt, Sxyz, g, height, diam);
	printf("q=%d\n",u.q);
// 	g=50.0;
	
	load=(double *)calloc(3*NumberOfSpheres,sizeof(double));
	
// 	for (i=1;i<=3*NumberOfSpheres;i++)
// 		load[i-1]=200.0*(drand48()-0.5);
	
	for (i=1;i<=NumberOfSpheres;i++)
		load[3*i-1]=-200.0;
 
	movingout=find_jammed(u.ia, u.ja, u.ar, u.l, load, u.q, u.a, gap, NumberOfSpheres, deltaMax, deltaMin);
	free(adjacent);
	mv=unique(movingout);
	printf("Number of Moving Particles=%d\n",mv.e);
	
	adjacent=(NEIGHBOUR *)calloc((int)(xtnt.Xxtnt*xtnt.Yxtnt*xtnt.Zxtnt),sizeof(NEIGHBOUR));
	if (adjacent == NULL)
		error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");

	for (i=0;i<(int)(xtnt.Xxtnt*xtnt.Yxtnt*xtnt.Zxtnt);i++)
		adjacent[i].next=NULL;
	uv=fload_data(adjacent, diam);
	s=uv.sp;
	NumberOfSpheres=uv.numsph;
	printf("NumberOfSpheres=%d\n",NumberOfSpheres);
	
	output_moving_vtk(s, NumberOfSpheres, mv.e, mv.k);
	
	free(mv.k);
	
	return 0;
}
