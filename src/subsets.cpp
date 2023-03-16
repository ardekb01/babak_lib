/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \
*  Copyright 2003 by Babak A. Ardekani, Ph.D.                *
*  ALL RIGHTS RESERVED.  No part of this program may be      *
*  used, transferred, or modified by any means without       *
*  prior written permission.  Making copies of any part      *
*  of this program for any purpose is a violation of         *
*  copyright laws.                                           *
\ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define _subsets


#include <stdlib.h>
#include <math.h>
#include <stdio.h>

// calculates the binomial coefficient nCm for 'small' n amd m
int binomialCoeff(int n,int m)
{
	double a=1.0;
	double b=1.0;

	if(m==n) return(1);
	if(m>n || n<0 || m<0) return(0);

	// calculates n!/m! 
	for(int i=m+1; i<=n; i++) a *= i;

	// calculates (n-m)!
	for(int i=1; i<=(n-m); i++) b *= i;

	return( (int)rint(a/b));
}

// This function returns an array of nCm integers. In each integer, exactly m bits are 1.
int *subsets(int n, int m)
{
	int maxsubsetindex; // maximum subset index
	int count;
	int shiftedx;
	int c,C;
	int *subsetlist;

	if(m>n)
	{
		printf("\n\nError: m must be less than or equal to n.\n\n");
		return(NULL);
	}

	if( n >= (int)(8*sizeof(int)) )
	{
		printf("\n\nError: n must be less than %lu.\n\n",8*sizeof(int));
		return(NULL);
	}

	maxsubsetindex = (int)pow(2.0,1.0*n)-1;

	C =  binomialCoeff(n,m);
	subsetlist = (int *)calloc(C,sizeof(int));

	c = 0;
	for(int x=0; x<=maxsubsetindex; x++)
	{
		count=0;
		for(int i=0; i<n; i++)
		{
			shiftedx = x >> i;
			count += (shiftedx & 1);
		}

		if(count==m) subsetlist[c++]=x;
	}

	return(subsetlist);
}
