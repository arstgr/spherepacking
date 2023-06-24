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
#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include "packing.h"

SP test(SPHERE *s, NEIGHBOUR *ne,int NumberOfSpheres, EXTENT ex)
{
	int i=0,j,mn;
	SP out;
	SPHERE *gamma;
	NE *nlist, *temp;
	
	for (gamma=s;gamma != NULL; gamma=gamma->next)
	{
		printf("num=%d    diam=%lf\n",gamma->num, gamma->diameter);
		i++;
	}
	printf("NumberOfSpheres=%d         i=%d\n",NumberOfSpheres,i);
	
	mn=(int)(ex.Xxtnt*ex.Yxtnt*ex.Zxtnt);
	printf("mn=%d\n",mn);
	for (j=0;j<mn;j++)
	{
		for (nlist=ne[j].next; nlist != NULL; nlist=nlist->next)
			printf("j=%d     number=%d    diam=%lf    num=%d\n",j,nlist->number, nlist->address->diameter, nlist->address->num);
	}
	
	out.sph=s;
	out.numsph=NumberOfSpheres;
	return out;
}
