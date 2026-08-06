#ifndef _stats_h

#include <math.h>

// Computes SSD
template <class TYPE> double SSD(TYPE *x, int n)
{
   double sum=0.0;
   double sum2=0.0;

   if(n<=0) return(0.0);

   for(int i=0; i<n; i++) 
   {
      sum += x[i];
      sum2 += (x[i]*x[i]);
   }

   return(sum2 - sum*sum/n);
}

// Computes the sample mean of a set of n observations {x_1,x_2,...,x_n} from a given distribution. 
template <class TYPE> double sample_mean(TYPE *x, int n)
{
   double mean=0.0;

   if(n<=0) return(0.0);

   for(int i=0; i<n; i++) mean += x[i];

   mean /= n;

   return(mean);
}

// Computes the unbiased sample variance of a set of n observations {x_1,x_2,...,x_n} from a given distribution.
template <class TYPE> double sample_variance(TYPE *x, int n, double &mean)
{
   double sum_of_sq=0.0;
   double sum=0.0;
   double var;

   if(n<2) 
   {
      mean = 0.0;
      return(0.0);
   }

   for(int i=0; i<n; i++)
   {
      sum_of_sq += x[i]*x[i];
      sum += x[i];
   }

   mean = sum/n;

   var = ( sum_of_sq - sum*sum/n )/(n-1.0); 

   return(var);
}

// Compute the dot product between n-dimensional vectors a and b
template <class TYPE> TYPE dot(TYPE *a, TYPE *b, int n)
{
   TYPE dotproduct=0.0;

   if(n<0) return(dotproduct);

   for(int i=0; i<n; i++)
   {
      dotproduct += a[i]*b[i];
   }

   return(dotproduct);
}

// Compute the dot product between n-dimensional vectors a and b
template <class TYPE> TYPE dot(TYPE *a, TYPE *b, TYPE *w, int n)
{
   TYPE dotproduct=0.0;

   if(n<0) return(dotproduct);

   for(int i=0; i<n; i++)
   {
      dotproduct += a[i]*b[i]*w[i];
   }

   return(dotproduct);
}

// Implements Eq. (2.2.1) of J. Cohen & P. Cohen (2nd ed.)
template <class TYPE> double standardDeviation(TYPE *x, int n)
{
	double sx,sxx; 
	double sd;

	if(n<=0) return(0.0);

	sx=sxx=0.0;
	for(int i=0; i<n; i++)
	{
		sx += x[i];
		sxx += x[i]*x[i];
	}

	sd = (sxx - sx*sx/n)/n;

	if(sd > 0.0) sd = sqrt(sd); else sd=0.0;
	
	return (sd);
}

// Computes the Pearson correlation coefficient between x and y.
// Implements Eq. (2.3.2) of J. Cohen & P. Cohen (2nd ed.)
template <class TYPE1, class TYPE2> float pearsonCorrelation(TYPE1 *x, TYPE2 *y, int n)
{
	double sx,sy,sxx,syy,sxy; 

	double Sxx, Sxy, Syy;

	double dum=0.0;

	sx=sy=sxx=syy=sxy=0.0;
	for(int i=0; i<n; i++)
	{
		sx += x[i];
		sxx += x[i]*x[i];

		sy += y[i];
		syy += y[i]*y[i];

		sxy += x[i]*y[i];
	}

	Sxx = n*sxx - sx*sx;
	Syy = n*syy - sy*sy;
	Sxy = n*sxy - sx*sy;

	if(Sxx*Syy > 0.0) dum = sqrt(Sxx*Syy);
	
	if(dum != 0.0)
		return ( (float)(Sxy/dum) );
	else
		return( 0.0 );
}

// Computes the Pearson correlation coefficient between x and y.
// Implements Eq. (2.3.2) of J. Cohen & P. Cohen (2nd ed.)
template <class TYPE1, class TYPE2> float pearsonCorrelation(TYPE1 *x, TYPE2 *y, short *msk, int n)
{
   int count;
   double sx,sy,sxx,syy,sxy; 
   double Sxx, Sxy, Syy;
   double dum=0.0;

	sx=sy=sxx=syy=sxy=0.0;
    count=0;
	for(int i=0; i<n; i++)
    if(msk[i] != 0)
	{
        count++;

		sx += x[i];
		sxx += x[i]*x[i];

		sy += y[i];
		syy += y[i]*y[i];

		sxy += x[i]*y[i];
	}

    if(count == 0) return(0.0);

	Sxx = count*sxx - sx*sx;
	Syy = count*syy - sy*sy;
	Sxy = count*sxy - sx*sy;

	if(Sxx*Syy > 0.0) dum = sqrt(Sxx*Syy);
	
	if(dum != 0.0)
		return ( (float)(Sxy/dum) );
	else
		return( 0.0 );
}

template <class TYPE> double independent_samples_t(TYPE *x1, int n1, TYPE *x2, int n2, int &df, double &meandiff)
{
   double var1, var2;
   double mean1, mean2;
   double var, sd;
   double t;

   df = n1+n2-2;

   if(n1<2 || n2<2) return(0.0);

   var1 = sample_variance(x1,n1,mean1);
   var2 = sample_variance(x2,n2,mean2);

   var = ( ( (n1-1)*var1 + (n2-1)*var2 ) / (n1 + n2 - 2.0) ) * (1.0/n1 + 1.0/n2);
   sd = sqrt(var);

   if( sd==0.0) return(0.0);

   meandiff = mean1 - mean2;

   t = (mean1 - mean2)/sd;

   return(t);
}

template <class TYPE> double paired_samples_t(TYPE *x1, TYPE *x2, int n, int &df, double &meandiff)
{
   double var;
   double mean;
   double sd;
   double t;

   df = n-1;

   if(n<2) return(0.0);

   for(int i=0;i<n;i++) x1[i] -= x2[i];

   var = sample_variance(x1,n,mean)/n;

   sd = sqrt(var);

   if( sd==0.0) return(0.0);

   meandiff = mean;

   t = mean/sd;

   return(t);
}

template <class TYPE> double one_sample_t(TYPE *x, int n, int &df, double &mean)
{
   double var;
   double sd;
   double t;

   df = n-1;

   if(n<2) return(0.0);

   var = sample_variance(x,n,mean)/n;

   sd = sqrt(var);

   if( sd==0.0) return(0.0);

   t = mean/sd;

   return(t);
}

#define _stats_h

#endif
