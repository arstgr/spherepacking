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

#include  <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "packing.h"

NA ncd( EXTENT t, int i, int j, int k,NEIGHBOUR *ne, int snum)
		/* The outputs are nlist[0..n-1] and n: the sorted list of neighbours and the size of this list respectively */
{
	int istencil[3], jstencil[3], kstencil[3];
	int u,v,w,m,uu,vv,ww,a,n=0;
	NE *temp, *tmp, *pos, *alpha;
	NA output;
	

	uu=0;
	for (u=-1;u<2;u++)
	{
		if (((i+u)<(int)t.Xxtnt)&&((i+u)>=0))
		{
			istencil[uu]=i+u;
			uu++;
		}
	}

	vv=0;
	for (v=-1;v<2;v++)
	{
		if (((j+v)<(int)t.Yxtnt)&&((j+v)>=0))
		{
			jstencil[vv]=j+v;
			vv++;
		}
	}

	ww=0;
	for (w=-1;w<2;w++)
	{
		if (((k+w)<(int)t.Zxtnt)&&((k+w)>=0))
		{
			kstencil[ww]=k+w;
			ww++;
		}
	}
	alpha=NULL;
	n=0;
	for (w=0;w<ww;w++)
	{
		for (v=0;v<vv;v++)
		{
			for (u=0;u<uu;u++)
			{
				m=(int)(istencil[u]+t.Xxtnt*jstencil[v]+t.Yxtnt*t.Xxtnt*kstencil[w]);
				for (temp=ne[m].next;temp!=NULL;temp=temp->next)
				{
					if (temp->number != snum)
					{
						tmp=(NE *)malloc(sizeof(NE));
						if (tmp==NULL)
							error(EXIT_FAILURE, ENOMEM ,"Insufficient Memory");
							
						*tmp=*temp;
						tmp->next=NULL;
						if (n==0)
						{
							alpha=tmp;
							n++;
						}
						else
						{
							pos=find_place(alpha,tmp); //pos is smaller than tmp
							if (pos != NULL)
							{
								tmp->next=pos->next;
								pos->next=tmp;
								n++;
							}
							else
							{
								tmp->next=alpha;
								alpha=tmp;
								n++;
							}
						}
					}
				}
			}
		}
	}

// 	pos=tmp=temp=NULL;
// 	for (temp=alpha; temp != NULL; temp=temp->next)
// 	{
// 		if (temp->number < snum)
// 			tmp=temp;
// 		else if (temp->number == snum)
// 			pos=temp;
// 		     else
// 			     break;
// 	}
// 	if ((pos != NULL)&&(tmp!=NULL))
// 	{
// 		tmp->next=pos->next;
// 		free(pos);
// 		n--;
// 	}
	output.location=alpha;
	output.size=n;
	return output;
}

