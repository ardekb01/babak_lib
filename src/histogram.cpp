#define _histogram

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <minmax.h>

double *findHistogram(short *im, int nv, int nb, int low, int high, int *bw_return);
double *findHistogram(short *im1, short *im2, int nv, int nb1, int nb2, int *bw1_r, int *bw2_r, int low1, int high1, int low2, int high2);

void trimExtremes(short *image, short *msk, int nv, float percent)
{
   short low, high;
   short min=0, max=0;
   int *histogram;
   int hsize;			/* histogram size */
   int b;
   int i;
   int nbv;

   int nmax;
   int n;

   nbv = minmax(image, msk, nv, min, max);

   hsize = max-min+1;

   histogram=(int *)calloc(hsize,sizeof(int));

   for(int i=0; i<nv; i++)
   {
      if(msk[i]>0)
      {
         b = image[i] - min;

         if(b>=0 && b<hsize) histogram[ b ]++;
      }
   }

   nmax = (int)( percent * nbv/100.0);

   n=0;
   for(i=0;i<hsize;i++)
   {
      n += histogram[i];
      if(n>nmax) break;
   }
   low=i+min; 

   n=0;
   for(i=0;i<hsize;i++)
   {
      n += histogram[hsize-1-i];
      if(n>nmax) break;
   }
   high=hsize-1-i+min; 

   for(int v=0; v<nv; v++)
   {
      if(msk[v]>0)
      {
         if(image[v]<low) 
         {
            image[v]=low;
         }
         else if(image[v]>high) 
         {
            image[v]=high;
         }
      }
   }
 
   free(histogram);
}

// This function creates a histogram with 'nb' bins, i.e., from 0 to nb-1.
// It returns the histogram as array hist[0], hist[1], ..., hist[nb-1].
// Image voxels with values in the interval [low,high] are 
// mapped to interval [0,nb-1]. All voxels with values less than 'low'
// are mapped to h[0] and all voxels with values greater than 'high' are
// mapped to h[nb-1].
double *findHistogram(short *im, int nv, int nb, int low, int high, int *bw_return)
{
	int	bw;		// histogram bin width
	int v;
	double *hist;
	int b;

	hist = (double *)calloc(nb, sizeof(double) );
	if(hist==NULL) 
	{
		printf("\nMemory allocation problem for 'hist'\n");
		return(NULL);
	}

	for(int i=0;i<nb;i++) hist[i]=0.0;

	bw = (int)ceil( (high-low)/(nb-1.0) );
	if(bw<=0) bw=1;

	*bw_return = bw;

	// im[i] is mapped to histogram bin (int)(v/bw) where:
	// v=0 if im[i]<=low 
	// v=high-low if im[i]>=high
	// v=im[i]-low if low<im[i]<high 
	// note that bw >= (high-low)/(nb-1.0), or (high-low)/bw <= (nb-1.0)
	// therefore when v=high-low we have v/bw <= nb-1.0
	for(int i=0;i<nv;i++)
	{
		v = im[i];

		if(v<=low) v=0;
		else if(v>=high) v=high-low;
		else v-=low;

		b = (int)(v/bw);

		if(b>=0 && b<nb) hist[b]++;
	};

	return(hist);
}

// This function creates a histogram with 'nb' bins, i.e., from 0 to nb-1.
// It returns the histogram as array hist[0], hist[1], ..., hist[nb-1].
// Image voxels with values in the interval [low,high] are 
// mapped to interval [0,nb-1]. All voxels with values less than 'low'
// are mapped to h[0] and all voxels with values greater than 'high' are
// mapped to h[nb-1].
double *findHistogram(short *im1, short *im2, int nv, int nb1, int nb2, int *bw1_r, int *bw2_r, int low1, int high1, int low2, int high2)
{
	int	bw1, bw2;
	int v1, v2;
	int b1, b2;
	double *hist;

	hist = (double *)calloc(nb2*nb1, sizeof(double) );
	if(hist==NULL) return(NULL);

	for(int i=0;i<nb2*nb1;i++) hist[i]=0.0;

	bw1 = (int)ceil( (high1-low1)/(nb1-1.0) );
	if(bw1<=0) bw1=1;
	*bw1_r=bw1;

	bw2 = (int)ceil( (high2-low2)/(nb2-1.0) );
	if(bw2<=0) bw2=1;
	*bw2_r=bw2;

	// im[i] is mapped to histogram bin (int)(v/bw) where:
	// v=0 if im[i]<=low 
	// v=high-low if im[i]>=high
	// v=im[i]-low if low<im[i]<high 
	// note that bw >= (high-low)/(nb-1.0), or (high-low)/bw <= (nb-1.0)
	// therefore when v=high-low we have v/bw <= nb-1.0
	for(int i=0;i<nv;i++)
	{
		v1 = im1[i];
		if(v1<=low1) v1=0;
		else if(v1>=high1) v1=high1-low1;
		else v1 -= low1;
		b1 = (int)(v1/bw1);

		v2 = im2[i];
		if(v2<=low2) v2=0;
		else if(v2>=high2) v2=high2-low2;
		else v2 -= low2;
		b2 = (int)(v2/bw2);

		
		hist[nb1*b2 + b1]++;
	};

	return(hist);
}

int otsu(double *h, int nb)
{
	double sum=0.0;
	int maxIndx=-1;
	double muT=0.0, mu=0.0;
	double sigma2=0.0, maxsigma2=0.0;
	double omega=0.0;
	
	sum = 0.0;
	for(int i=0; i<nb; i++) sum+=h[i];

	if(sum<=0.0) return(-1);

	muT=0.0;
	for(int i=0; i<nb; i++) 
	{
		h[i]/=sum;
		muT += i*h[i];
	}

	for(int k=0; k<nb; k++)
	{
		sigma2=mu=omega=0.0;

		for(int i=0; i<=k; i++)
		{
			mu += i*h[i];
			omega += h[i];
		}

		if(omega==0.0 || omega==1.0) sigma2=0.0;
		else sigma2 = (muT*omega-mu)*(muT*omega-mu) / ( omega * (1.0-omega) );

		if(sigma2>=maxsigma2) { maxIndx=k; maxsigma2=sigma2; }

		// printf("k=%d sigma2=%lf mu=%lf omega=%lf\n", k, sigma2, mu, omega);
	}

	return(maxIndx);
}

/*
double *findHistogram(short *im, int nv, int *nb, short *min, short *max)
{
   double *hist;

   minmax(im, nv, min, max);

   *nb = (*max) - (*min) + 1;

   hist = (double *)calloc( *nb, sizeof(double) );

   if(hist==NULL) 
   {
      printf("\nMemory allocation problem for 'hist'\n");
      return(NULL);
   }

   for(int i=0;i<*nb;i++) hist[i]=0.0;

   for(int i=0;i<nv;i++)
   {
      hist[im[i] - (*min)] += 1.0;
   };

	return(hist);
}
*/

double *findHistogram(short *im, int nv, int *nb, short &min, short &max)
{
   double *hist;

   minmax(im, nv, min, max);

   *nb = max - min + 1;

   hist = (double *)calloc( *nb, sizeof(double) );

   if(hist==NULL) 
   {
      printf("\nMemory allocation problem for 'hist'\n");
      return(NULL);
   }

   for(int i=0;i<*nb;i++) hist[i]=0.0;

   for(int i=0;i<nv;i++)
   {
      hist[im[i] - min] += 1.0;
   };

	return(hist);
}
