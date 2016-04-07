/************************************************************************************
FronTier is a set of libraries that implements differnt types of Front Traking algorithms.
Front Tracking is a numerical method for the solution of partial differential equations 
whose solutions have discontinuities.  


Copyright (C) 1999 by The University at Stony Brook. 
 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

******************************************************************************/


/*
*			sesspln.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include <cdecs.h>


/*
*				splcomp():
*
*	This routine  computes the spline function for given points
*	x, y where x is the independent variable and y the dependent
*	using the "not a knot" condition. slp is the output and conains
*	the slope of the cubic spline for each i, 0 <= i < n.
*	The routine is based on C. de Boor's book, A Practical
*	Guide to Splines.
*	Note that imput is assumed to be in log varaiables.
*/

EXPORT	int	splcomp(
	double	*x,
	double	*y,
	double	*slp,
	int	n,
	double	*c,
	double	*d,
	double	*e,
	double	*dx,
	double	*dy)
{
	int	i, info;

	for (i = 0; i < n; i++)
	{
	    slp[i] = 0.0;
	    c[i] = 0.0;
	    d[i] = 0.0;
	    e[i] = 0.0;
	}
	if (n == 2)
	{
	    slp[1] = slp[0] = (y[0] - y[1]) / (x[0] - x[1]);
	    return -1;
        }
	for (i = 0; i < n-1; i++)
	{
	    dy[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
	    dx[i] = x[i + 1] - x[i];
	}
	d[0] = dx[1];
	e[0] = x[2] - x[0];
	slp[0] = ((dx[0]+2.0*(x[2]-x[0]))*dx[1]*dy[0] + dx[0]*dx[0]*dy[1]) /
		(x[2] - x[0]);
	if (n > 3)
	{
	    d[n-1] = dx[n-3];
	    c[n-1] = x[n-1] - x[n-3];
	    slp[n-1] = (dx[n-2]*dx[n-2]*dy[n-3] +
			(2.0*(x[n-1] - x[n-3]) + dx[n-2])*dx[n-3]*dy[n-2]) /
			(x[n-1] - x[n-3]);
	}
	if (n == 3)
	{
	    d[n-1] = 1.0;
	    c[n-1] = 1.0;
	    slp[n-1] = 2.0*dy[n-1];
	}
	for (i = 1; i < n-1; i++)
	{
	    slp[i] = 3.0*(dx[i] * dy[i-1] + dx[i-1]*dy[i]);
	    d[i] = 2.0*(dx[i-1] + dx[i]);
	    c[i] = dx[i];
	    e[i] = dx[i-1];
	}
	gtsl(&n,c,d,e,slp,&info);

	return info-1;
}		/*end splcomp*/

/*
*				spline():
*
*	To compute y = f(x) where f is a spline function with x
*	as the independent variable, s1 and s2 are slopes associated
*	wtih the data points y1,y2,x1,and x2. Output is y.
*/

EXPORT	void	spline(
	double	x1,
	double	x2,
	double	y1,
	double	y2,
	double	s1,
	double	s2,
	double	x,
	double	*y,
	double	*dy)
{
	double c1, c2, c3, c4, xp;

	if (y1 == 0.0 && y2 == 0.0)
	{
	    *y = 0.0;
	    *dy = 0.0;
	    return;
	}

	c1 = y1;
	c2 = s1;
	c4 = s1 + s2;
	c4 -= (y1 - y2) / (x1 - x2) * 2;
	c4 /= (x1 - x2) * (x1 - x2);
	c3 = (y1 - y2) / (x1 - x2);
	c3 -= s1;
	c3 /= x2 - x1;
	c3 -= c4 * (x2 - x1);
	xp = x - x1;
	*y = c1 + c2*xp + c3*xp*xp + c4*xp*xp*xp;
	*dy = c2 + 2.0*c3*xp + 3.0*c4*xp*xp;
}		/*end spline*/

/*
*			splcomp2():
*
*	This routine  computes the spline function for given points
*	x, y where x is the independent variable and y the dependent
*	using the "complete" condition. slp is the output and conains
*	the slope of the cubic spline for each i, 0 <= i < n.
*	Note that imput is assumed to be in log varaiables.
*/

EXPORT	int	splcomp2(
	double	*x,
	double	*y,
	double	*slp,
	int	n,
	double	*c,
	double	*d,
	double	*e,
	double	*dx,
	double	*dy,
	double	a,
	double	b)
{
	int	i, info;

	for (i = 0; i < n; i++)
	{
	    slp[i] = 0.0;
	    c[i] = 0.0;
	    d[i] = 0.0;
	    e[i] = 0.0;
	}
	for (i = 0; i < n-1; i++)
	{
	    dy[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
	    dx[i] = x[i + 1] - x[i];
	}
	d[0] = 1.0;
	e[0] = 0.0;
	slp[0] = a;
	d[n-1] = 1.0;
	c[n-1] = 0.0;
	slp[n-1] = b;

	for (i = 1; i < n-1; i++)
	{
	    slp[i] = 3.0*(dx[i] * dy[i-1] + dx[i-1]*dy[i]);
	    d[i] = 2.0*(dx[i-1] + dx[i]);
	    c[i] = dx[i];
	    e[i] = dx[i-1];
	}

	gtsl(&n,c,d,e,slp,&info);
	return info - 1;
}		/* end splcomp2*/

