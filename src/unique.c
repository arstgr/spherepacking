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
#include "packing.h"

int compare_ints( const void *a, const void *b)
{
	const int *ia =(const int *)a;
	const int *ib =(const int *)b;
	return ( *ia > *ib )-( *ia < *ib );
}

MVPR unique(MOVING mov)
{
	int i,*k;
	MOVING *temp,*tmp;
	int d,e,f,g;
	FILE *myfile;
	MVPR mv;
	
	k=(int *)calloc(mov.number,sizeof(int));
	
	for (i=0,temp=mov.next;temp!=NULL;i++)
	{
		tmp=temp;
		k[i]=temp->number;
		temp=temp->next;
		free(tmp);
	}	
	
	d=mov.number;
	qsort(k,d,sizeof(int),compare_ints);
	g=k[d-1]+10;
	f=k[0];
	e=d;
	for (i=1;i<d;i++)
	{
		if (k[i]==f)
		{
			k[i]=g;
			e--;
		}
		else
			f=k[i];
	}
	qsort(k,d,sizeof(int),compare_ints);
	
	myfile=fopen("k.packing","w");
	if (d)
	{
		for (i=0;i<e;i++)
			fprintf(myfile,"%d\n",k[i]);
		fclose(myfile);
	}
	else
	{
		fprintf(myfile,"%d\n",-1);
		fclose(myfile);
	}
	mv.k=k;
	mv.e=e;
	return mv;
}
