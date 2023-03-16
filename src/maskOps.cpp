#define _maskOps

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <babak_lib.h>

//////////////////////////////////////////////////////////////////////////////////////////////////
char *find_foreground_mask(short *im, int nv, int nb, int nclass, int niter, short *thresh)
{
	char *mask;
	double *hist;
	double *mean;
	double *var;
	double *p;
	double *fit;
	short *label;
	//FILE *fp;
	int low,high,v;
	int bw;
	int T;
	int nbv;

	mask=(char *)calloc(nv,sizeof(char));

	nbv=0;
	for(int i=0; i<nv; i++) if(im[i]>0) nbv++;

	mean = (double *)calloc(nclass,sizeof(double));
	var = (double *)calloc(nclass,sizeof(double));
	p = (double *)calloc(nclass,sizeof(double));
	fit = (double *)calloc(nb,sizeof(double));
	label = (short *)calloc(nb,sizeof(short));
	if(mean==NULL || var==NULL || p==NULL || fit==NULL || label==NULL)
	{
		printf("\nMemory allocation problem in 'find_foreground_mask'\n");
		return(NULL);
	}

	{
		short *im_mskd;
		int k=0;

		im_mskd = (short *)calloc(nbv,sizeof(short));

		if(im_mskd==NULL) 
		{
			printf("\nMemory allocation problem for 'im_mskd', aborting ...\n");
			exit(0);
		}

		for(int i=0; i<nv; i++) if(im[i]>0) im_mskd[k++]=im[i];

		setLowHigh(im_mskd, nbv, &low, &high, 0.01);
		hist=findHistogram(im_mskd, nbv, nb, low, high, &bw);

		free(im_mskd);
	}

	for(int b=0; b<nb; b++) hist[b]/=nbv;

	EMFIT1d(hist, fit, label, nb, mean, var, p, nclass, niter);

	T=0;
	for(int b=1; b<nb-1; b++)
	if( fit[b-1]>fit[b] && fit[b+1]>fit[b] ) 
	{
		T=b;
		break;
	}
			
	*thresh=0;
	for(int i=0;i<nv;i++) 
	if(im[i]>0)
	{
		v = im[i];

		if(v<=low) v=0;
		else if(v>=high) v=high-low;
		else v-=low;

		v = (int)(v/bw);

		if(v==T && im[i] > *thresh ) *thresh=im[i];
	}

	for(int i=0;i<nv;i++) if(im[i] <= *thresh) mask[i]=0; else mask[i]=1;

	free(mean);
	free(var);
	free(p);
	free(fit);
	free(label);
	free(hist);

	return(mask);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////

// mask values are 0 and non-zero (don't have to be 0 or 1)
short *mask_image(short *im, char *mask, int nv, int *nbv_out)
{
	int nbv;
	short *im_mskd;
	int k;

	nbv=0;
	for(int i=0; i<nv; i++) if(mask[i]!=0) nbv++;

	im_mskd = (short *)calloc(nbv,sizeof(short));

	k=0;
	for(int i=0; i<nv; i++) if(mask[i]!=0) im_mskd[k++]=im[i];

	*nbv_out = nbv;

	return(im_mskd);
}
//////////////////////////////////////////////////////////////////////////////////////////////////

void cc_filter(char *mask, int nx, int ny, int nz)
{
	int ncc, maxsize;
	int nv;

	nv = nx*ny*nz;

	max_Connected_Component(mask, nx, ny, nz, &ncc, &maxsize);

	for(int i=0; i<nv; i++) if(mask[i]==1) mask[i]=0; else mask[i]=1;

	max_Connected_Component(mask, nx, ny, nz, &ncc, &maxsize);

	for(int i=0; i<nv; i++) if(mask[i]==1) mask[i]=0; else mask[i]=1;
}
