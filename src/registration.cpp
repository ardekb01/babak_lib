#define _registration

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <babak_lib.h>

#define NBIN 256
#define MCC     2048 /* maximum number of allowed 2D connected components */
#define MINSIZE  101 /* minimum size of a 2D connected components */

///////////////////////////////////////////////////////////////////////////////////////////////
// babak_lib global variables

float (*interpolator)(float x, float y, float z, short *array, int nx, int ny, int nz, int np);
float P[6];
struct im_params IP;
///////////////////////////////////////////////////////////////////////////////////////////////

float *xjit, *yjit, *zjit;

//float alpha_max=-10000.0;
//float alpha_min=10000.0;

static float PCOM[6];
static float XICOM[6];
static int ncom;
static float global_min;
static float Pmin[6];
static float v1,v2,v3,v4;
static float w1,w2;
static float (*obj_fnc)(short *KMI, float *P, struct im_params *IP);
static int hist[NBIN];
static int low,high;        /* range of pixel considered in histogram computation */


static void  findInterval(short *KMI, float *ax,float *bx,float *cx, float *fa, float *fb, float *fc, struct im_params *IP);
static float f1dim(short *KMI, float x, struct im_params *IP);
static float optimize(short *KMI, float ax,float bx,float cx,float tol, float *xmin, struct im_params *IP);
static float minimize1D(short *KMI, float *X, int ndim, struct im_params *IP);
static float Gradient_Descent(short *KMI, int	ndim, float ftol, struct im_params *IP);
float *transformation(float x, float y, float z, float ax, float ay, 
float az, float sx, float sy, float sz, int rX, int rY, int rZ, char *code);
unsigned char linearInterpolatorUC(float x, float y, float z, unsigned char *array, int nx, int ny, int nz, int np);
static float costFunction1(short *KMI, float *P, struct im_params *IP);
static float newCostFunction(short *KMI, float *P, struct im_params *IP);
//static float costFunction2(short *KMI, float *P, struct im_params *IP);
void scale_short_minmax(short *imagein, unsigned char **imageout, int np, int min,int max);
void label_3d_cc(short *KMI,unsigned short label,int i,int j, int k,int *size,short CC, struct im_params *IP);
static void resetCC(short *KMI,int i,int j,int k,unsigned char CC, struct im_params *IP);
int label_CCI(short *KMI, int size_thresh,struct im_params * IP, int nvoxels);
short *KMcluster(short *ccImage, short *im_in, int nclass, int maxiter, int thresh, int Low, int High, int nv);

static short PL2d(short *im, int nc, int i);
static short PR2d(short *im, int nc,  int i);
static short PA2d(short *im, int nc, int i);
static short PB2d(short *im, int nc, int np, int i);

static void label_2d_cc(short *im,int nc,int np,int i,short oldlabel, short newlabel, int *size);
static void delete_small_ccs(short *im, int nc, int nr, int nz);
void cca(short *im, int nx, int ny, int nz);
int findThresholdLevel(short *image_in, int nv);

static short label;
static int size[MCC];

static int nclass=8;
static int maxiter=500;

float *findTransMatrix(short *trg, int Tnx, int Tny, int Tnz, float Tdx, float Tdy, float Tdz, 
short *obj, int Onx, int Ony, int Onz, float Odx, float Ody, float Odz)
{
   char transcode[5]="ZXYT";
   int thresh;
   short *ccImage;
   short *im_out;
   int size_thresh=100;
   int OL,OH;
   float *T;

   int Tnp,Tnv,Onv;
	
   for(int i=0; i<6; i++) P[i]=0.0;

   IP.FROM = 1;
   IP.TO = Tnz;

   IP.nx2=Tnx;
   IP.ny2=Tny;
   IP.nz2=Tnz;
   IP.dx2=Tdx;
   IP.dy2=Tdy;
   IP.dz2=Tdz;

   IP.nv2=IP.nx2*IP.ny2*IP.nz2;
   IP.np2=IP.nx2*IP.ny2;

   IP.xc2=IP.dx2*(IP.nx2-1)/2.0; /* +---+---+ */
   IP.yc2=IP.dy2*(IP.ny2-1)/2.0;
   IP.zc2=IP.dz2*(IP.nz2-1)/2.0;
	
   IP.sf = 1;

   IP.nx1=Onx;
   IP.ny1=Ony;
   IP.nz1=Onz;
   IP.dx1=Odx;
   IP.dy1=Ody;
   IP.dz1=Odz;

   IP.nv1=IP.nx1*IP.ny1*IP.nz1;
   IP.np1=IP.nx1*IP.ny1;

   IP.xc1=IP.dx1*(IP.nx1-1)/2.0; /* +---+---+ */
   IP.yc1=IP.dy1*(IP.ny1-1)/2.0;
   IP.zc1=IP.dz1*(IP.nz1-1)/2.0;

   Tnv = Tnx*Tny*Tnz;
   Tnp = Tnx*Tny;
   Onv = Onx*Ony*Onz;

   thresh=findThresholdLevel(trg, Tnv);
   printf("Threshold level = %d\n",thresh);

   ccImage=(short *)calloc(Tnv,sizeof(short));

   /* thresholding */
   for(int i=0;i<Tnv;i++)
   if( trg[i] <= thresh)
      ccImage[i]=0;
   else
      ccImage[i]=1;

   cca(ccImage, Tnx, Tny, Tnz);

   printf("Number of classes = %d\n",nclass);
   printf("Maximum number of iterations allowed = %d\n",maxiter);

   printf("\nK-means clustering ...\n");
   im_out=KMcluster(ccImage, trg, nclass, maxiter, thresh, low, high , Tnv);

   free(ccImage);
   ccImage=im_out;

   // make sure edges are zero
   for(int k=0; k<Tnz; k++)
   {
      for(int i=0; i<Tnx; i++) im_out[k*Tnp+0*Tnx+i]=0;
      for(int i=0; i<Tnx; i++) im_out[k*Tnp+(Tny-1)*Tnx+i]=0;
      for(int j=0; j<Tny; j++) im_out[k*Tnp+j*Tnx+0]=0;
      for(int j=0; j<Tny; j++) im_out[k*Tnp+j*Tnx+(Tnx-1)]=0;
   }

   fprintf(stdout,"\nCluster CC analysis ...\n");
   IP.CCI=(unsigned *)calloc(Tnv,sizeof(unsigned int));
   IP.NCC=label_CCI(ccImage,size_thresh,&IP,Tnv);

   fprintf(stdout,"\nCost function minimization ...\n");
   IP.Size=(int *)calloc(IP.NCC,sizeof(int));
   IP.S=(float *)calloc(IP.NCC,sizeof(float));
   IP.SS=(float *)calloc(IP.NCC,sizeof(float));

	IP.data1=(unsigned char *)malloc(IP.nv1);
	setLowHigh(obj, Onv, &OL, &OH);
	scale_short_minmax( obj, &(IP.data1), IP.nv1, OL, OH);
	IP.max1=255;

	obj_fnc=costFunction1;

	// (void)Gradient_Descent(ccImage,6,0.001,&IP);  4/10/03
	(void)Gradient_Descent(ccImage,6,0.00001,&IP);

	free(ccImage);
	free(IP.data1);
	free(IP.Size);
	free(IP.S);
	free(IP.SS);
	free(IP.CCI);

   T=transformation(P[0],P[1],P[2],P[3],P[4],P[5], 1.,1.,1.,0,0,0,transcode);

	return(T);
}

static void findHistogram(short *im, int nv)
{
	int j;			/* histogram index */
	int	delta;		/* histogram bin width */
	short v;         /* voxel value */

	delta = (int)( (high-low)/(1.0*NBIN) ) + 1;

	for(j=0;j<NBIN;j++)
		hist[j]=0;

	for(int i=0; i<nv; i++)
	{
		v=im[i];

		/* It is important to not let v=high because if v=high then j=NBIN and
		array h will be out of index. */
		if( v!=0 && v>=low && v<=high)
		{
			j = (int)( (v-low)/(1.0*delta) );
			hist[j]++;
		}
   }
}

// revised version
int loadTransformation( char *filename, float *T)
{
   FILE *fp;
   char line[1000];
   char s[11];
   int i;

   if(T==NULL)
   {
      printf("Warning: Null pointer passed as argument to loadTransformation()\n");
      return(1);
   }

   if( filename==NULL || filename[0]=='\0')
   {
      for(int i=0; i<16; i++) T[i]=0.0;
      T[0]=T[5]=T[10]=T[15]=1.0;

      return(1);
   }

   strcpy(s,"0123456789");

   fp=fopen(filename,"r");

   if(fp==NULL)
   {
      printf("Warning: Could not open %s passed as argument to loadTransformation()\n", filename);
      return(1);
   }

   i=0;
   while( fgets(line,1000,fp) != NULL  && i != 4 ) 
   {
      if(line[0]!='#' && line[0]!='\n' && strpbrk(line,s)!=NULL )
      {
         sscanf(line,"%f  %f  %f  %f\n",T+4*i+0,T+4*i+1,T+4*i+2,T+4*i+3);
         i++;
      }
   }
   fclose(fp);

   if(i!=4)
   {
      printf("Warning: Could not read %s passed as argument to loadTransformation()\n", filename);
      return(1);
   }

   return(0);
}

int findThresholdLevel(short *image_in, int nv)
{
	int j;
	double histmax1, histmax2; 
	int histmin;
	int delta;		/* histogram bin width */
	int thresh;
	int jmax1, jmax2,jmin;

	
	printf("\nAutomatic Thresholding ...\n");

	// from the image header, find the lower and upper limits of the image grey levels
	setLowHigh(image_in, nv, &low, &high);

	findHistogram(image_in,nv);

	histmax1=0.0;
	jmax1=1;
	for(j=1;j<NBIN-1;j++)
	if( hist[j]>=hist[j-1] && hist[j]>=hist[j+1] && hist[j]>histmax1)
	{
		histmax1=hist[j];
		jmax1 = j;
	}

	histmax2=0.0;
	jmax2=1;
	for(j=1;j<NBIN-1;j++)
	if( hist[j]>=hist[j-1] && hist[j]>=hist[j+1] && hist[j]>histmax2 && j!=jmax1)
	{
		histmax2=hist[j];
		jmax2 = j;
	}

	histmin=hist[jmax1];
	jmin=jmax1;
	for(j=jmax1;j<=jmax2;j++)
	if( hist[j]<histmin )
	{
		jmin=j;	
		histmin=hist[j];
	}

	delta = (int)( (high-low)/(1.0*NBIN) ) + 1;

	thresh = (jmin+1)*delta + low - 1;

	return(thresh);
}

void cca(short *im, int nx, int ny, int nz)
{
	int i,j,k;

	int np,nv;
	int maxsize;
	int clss;
	int *codebook;	/* label-class equivalence code book */
	int *class_size;

	short *dumptr;
	short *dumptr2;
	short p1,p2;

	unsigned char *equiv_table;

	np=nx*ny;
	nv=np*nz;

	delete_small_ccs(im, nx, ny, nz);

	/* allocate memory for the equivalence table */
	equiv_table=(unsigned char *)calloc( label*label, sizeof(unsigned char) );

	/* set the diagonal enteries of the equivalence table to 1 */
   for(i=0;i<label;i++) equiv_table[i*label + i]=1;

	/* set the enteries of the equivalence table */
	for(k=0;k<nz-1;k++)
	{
      dumptr=im+np*k;

      for(i=0;i<np;i++)
   	if( (p1=dumptr[i]) && (p2=dumptr[i+np]) )
         equiv_table[p1*label + p2]=equiv_table[p2*label + p1]=1;
	}

	/* find the transitive enclosure of the equivalence table */
   for(j=0;j<label;j++)
   {
      for(i=0;i<label;i++)
      if( equiv_table[i*label + j] )
      {
         for(k=0;k<label;k++)
         if(equiv_table[i*label + k] || equiv_table[j*label + k])
            equiv_table[i*label + k] = 1;
      }
   }

	codebook=(int *)calloc( label, sizeof(int) );
	class_size=(int *)calloc( label, sizeof(int) );
	clss=1;
   for(i=2;i<label;i++)
	{
   	for(j=2;j<label;j++)
   	if( equiv_table[i*label + j] )
   	{
			codebook[j]=clss;
			class_size[clss] += size[j];

      	for(k=i;k<label;k++)
         	equiv_table[k*label + j]=0;
   	}

		if(class_size[clss]>0)
			clss++;
	}

	maxsize=class_size[1];
	k=1;
	for(i=2;i<clss;i++)
	{
		if(class_size[i]>maxsize)
		{
			maxsize=class_size[i];
			k=i;
		}
	}

   for(i=0;i<nv;i++)
   {
      if( im[i] ) 
      {
         if( codebook[ im[i] ]==k )
         {
            im[i] = 1;
         }
         else
         {
            im[i] = 0;
         }
      }
   }

   for(k=0;k<nz;k++)
   {
      dumptr=im+np*k;

		for(j=0;j<ny;j++)
		{
			dumptr2 = dumptr+j*nx;

			for(i=0;i<nx;i++)
			if( dumptr2[i] == 0 )
			{
				dumptr2[i] = 2;
			}
			else
			{
				/* if(i!=0) dumptr[j*nx + i - 1] = 0; */
				break;
			}

         for(i=nx-1;i>=0;i--)
         if( dumptr2[i] == 0 )
         {
            dumptr2[i] = 2;
         }
         else
         {
            /* if(i!=0) dumptr[j*nx + i - 1] = 0; */
            break;
         }
		}

      for(i=0;i<nx;i++)
      {
         for(j=0;j<ny;j++)
         if( dumptr[j*nx + i] == 2 )
			{
				continue;
			}
         else if( dumptr[j*nx + i] == 0 )
         {
            dumptr[j*nx + i] = 2;
         }
         else
         {
            /* if(i!=0) dumptr[j*nx + i - 1] = 0; */
            break;
         }

         for(j=ny-1;j>=0;j--)
         if( dumptr[j*nx + i] == 2 )
         {
            continue;
         }
         else if( dumptr[j*nx + i] == 0 )
         {
            dumptr[j*nx + i] = 2;
         }
         else
         {
            break;
         }
      }

		for(i=nx;i<(np-nx);i++)
		if( dumptr[i]==0 && ( dumptr[i-1]==2 || dumptr[i+1]==2 || 
		dumptr[i-nx]==2 || dumptr[i+nx]==2 ) )
		{
			size[0]=0;
      	label_2d_cc(dumptr,nx,np,i,0,2,&size[0]);
		}

		for(i=0;i<np;i++)
			if(dumptr[i]==2) dumptr[i]=0;
			else dumptr[i]=1;
   }

	free(equiv_table);
	free(codebook);
	free(class_size);
}

static void delete_small_ccs(short *im, int nc, int nr, int nz)
{
	int i,k;

	int np;
	short *dumptr;

	np=nc*nr;

	size[0]=size[1]=0;
	label=2;
	for(k=0;k<nz;k++)
	{
		dumptr=im+np*k;
		for(i=0;i<np;i++)
		if(dumptr[i]==1)
		{
			size[label]=0;
	   		label_2d_cc(dumptr,nc,np,i,1,label,&size[label]);

			if(size[label] < MINSIZE)
			{
	   		label_2d_cc(dumptr,nc,np,i,label,0,&size[label]);
			}
			else
			{
				label++;
			}

			if(label>MCC) break;
		}
	}

	// printf("Number of labels = %d\n",label); 
}

static void label_2d_cc(short *im,int nc,int np,int i,short oldlabel, short newlabel, int *size)
{
	(*size)++;

	im[i]=newlabel;

	if( PL2d(im,nc,i)==oldlabel )
	   	label_2d_cc(im,nc,np,i-1,oldlabel,newlabel,size);
   if( PR2d(im,nc,i)==oldlabel )
	   	label_2d_cc(im,nc,np,i+1,oldlabel,newlabel,size);
   if( PA2d(im,nc,i)==oldlabel )
	   	label_2d_cc(im,nc,np,i-nc,oldlabel,newlabel,size);
   if( PB2d(im,nc,np,i)==oldlabel )
	   	label_2d_cc(im,nc,np,i+nc,oldlabel,newlabel,size);
}

static short PL2d(short *im, int nc, int i)
{
   if( i%nc == 0 )
      return( -1 );
   else
      return( im[i-1] );
}

static short PR2d(short *im, int nc,  int i)
{
   if( (i + 1)%nc == 0 )
      return( -1 );
   else
      return( im[i+1] );
}

static short PA2d(short *im, int nc, int i)
{
   if( (i - nc) < 0 )
      return( -1 );
   else
      return( im[i-nc] );
}

static short PB2d(short *im, int nc, int np, int i)
{
   if( (i + nc) >= np  )
      return( -1 );
   else
      return( im[i+nc] );
}

short *KMcluster(short *ccImage, short *im_in, int nclass, int maxiter, int , int Low, int High, int nv)
{
   int i,n;
   int h; // h=hist[n] as n varies 

   int converge; // converge is set to 1 when the algorithm has converged.

   int delta;		// histogram bin width
   int k;
   int v;
   int   iter; // iter=current Kmeans iteration.
   short *im_out;
	
   char    db[NBIN];
   double *mean;
   double *oldmean;
   double *PR;

   int clss=0;   // the class with minimum distance from a given pixel value 

   float d;    // square distance between a given pixel value and a class mean 

   float dmin; // minimum square distance between a given pixel value and a class mean 

   float dmax; // dmax=maximum possible square distance between a given pixel value and class mean 

   mean = (double *)calloc(nclass,sizeof(double));
   oldmean = (double *)calloc(nclass,sizeof(double));
   PR = (double *)calloc(nclass,sizeof(double));
   im_out = (short *)calloc(nv,sizeof(short));

   // from the image header, find the lower and upper limits of the image grey levels 
   low  = Low;
   high = High;

   findHistogram(im_in,nv);

   delta = (int)( (high-low)/(1.0*NBIN) ) + 1;

   for(i=0;i<NBIN;i++)
      db[i]=0;

   // initialize class means 
   for(i=0;i<nclass;i++)
      mean[i]=(i+1)*NBIN*1.0/(nclass+1.0);

   dmax=256*256;

   iter=1;

	do {
      if(!(iter%50)) printf("\n\titeration %d ...",iter);

      /* initialize npts and olmean */
      for (i = 0 ; i < nclass; i++) {
         PR[i]=0.0;
         oldmean[i]=mean[i];
         mean[i]=0.0;
      }

		for(n=0;n<NBIN;n++)
      if( (h=hist[n])!=0) 
		{
            dmin=dmax;
            for(i=0;i<nclass;i++) 
				{
               d=oldmean[i]-n;
               d *= d;
               if(d<dmin) 
					{
                  dmin=d;
                  clss=i;
               }
            }

            db[n]=clss;

            PR[clss] += h;

            mean[clss] += h*n;
		}

      // Check for convergence. Set converge=1 if the means
      // didn't chage from the previous iteration. 
		converge = 1;
		for (i = 0 ; i < nclass; i++)
		if(PR[i]!=0.0) 
		{
			mean[i]/=PR[i];
			if (mean[i] != oldmean[i])
				converge = 0;
		}

      iter++;
   } while(!converge && iter<=maxiter);

   for(i=0;i<nv;i++)
   if( ccImage[i] )
   {
      v=im_in[i];
      if(v<low) v=low;
      if(v>high) v=high;

      k = (int)( (v-low)/(1.0*delta) );
      im_out[i]=db[k]+1;
   }

   free(mean);
   free(oldmean);
   free(PR);
   return(im_out);
}

int label_CCI(short *KMI, int size_thresh,struct im_params * IP, int nvoxels)
{
   int i,j,k;
   unsigned short label;
   int size, maxsize; 
   int NCC; 	 // number of connected components in CCI 
   short CC; 

   label=0;
   NCC=0;
   IP->NP=0;
   maxsize=0;

   // make sure edges are zero
   for(int kk=0; kk<IP->nz2; kk++)
   {
      for(int ii=0; ii<IP->nx2; ii++) KMI[kk*IP->np2+0*IP->nx2+ii]=0;
      for(int ii=0; ii<IP->nx2; ii++) KMI[kk*IP->np2+(IP->ny2-1)*IP->nx2+ii]=0;
      for(int jj=0; jj<IP->ny2; jj++) KMI[kk*IP->np2+jj*IP->nx2+0]=0;
      for(int jj=0; jj<IP->ny2; jj++) KMI[kk*IP->np2+jj*IP->nx2+(IP->nx2-1)]=0;
   }

   for(k=0;k<(IP->TO-IP->FROM+1);k++) 
   for(j=0;j<IP->ny2;j++)
   for(i=0;i<IP->nx2;i++) 
   if( (CC=KMI[k*IP->np2+j*IP->nx2+i]) && !IP->CCI[k*IP->np2+j*IP->nx2+i] ) 
   {  
      label++;
      size=0;
      label_3d_cc(KMI,label,i,j,k,&size,CC,IP);
      if(size<=size_thresh) 
      {
         resetCC(KMI,i,j,k,CC,IP);
         label--;
      } 
      else 
      {
         NCC++;
         IP->NP += size;
         if(size>maxsize) 
         {
            maxsize=size;
         }
      }
   }

// printf("NCC=%d NP=%d\n",NCC,IP->NP);

   j=0;
	for(i=0;i<nvoxels;i++)
	if(KMI[i]) 
	{
			if( (j%IP->sf)!=0 ) 
				KMI[i]=0;
			j++;
	}

	return(NCC);
}

static void resetCC(short *KMI,int i,int j,int k,unsigned char CC, struct im_params *IP)
{
	KMI[k*IP->np2+j*IP->nx2+i]=0;

	if(KMI[k*IP->np2+j*IP->nx2+i-1]==CC)
      resetCC(KMI,i-1,j,k,CC,IP);
	if(KMI[k*IP->np2+j*IP->nx2+i+1]==CC)
      resetCC(KMI,i+1,j,k,CC,IP);
   if(KMI[k*IP->np2+(j+1)*IP->nx2+i]==CC)
      resetCC(KMI,i,j+1,k,CC,IP);
	if(KMI[k*IP->np2+(j-1)*IP->nx2+i]==CC)
      resetCC(KMI,i,j-1,k,CC,IP);
	if(k>0 && KMI[(k-1)*IP->np2+j*IP->nx2+i]==CC)
      resetCC(KMI,i,j,k-1,CC,IP);
	if(k<IP->TO-IP->FROM && KMI[(k+1)*IP->np2+j*IP->nx2+i]==CC)
      resetCC(KMI,i,j,k+1,CC,IP);
}

void label_3d_cc(short *KMI,unsigned short label,int i,int j, int k,int *size,short CC, struct im_params *IP)
{
   (*size)++;

   IP->CCI[k*IP->np2+j*IP->nx2+i]=label;

   if(KMI[k*IP->np2+j*IP->nx2+i-1]==CC && !IP->CCI[k*IP->np2+j*IP->nx2+i-1])
      label_3d_cc(KMI,label,i-1,j,k,size,CC,IP);
   if(KMI[k*IP->np2+j*IP->nx2+i+1]==CC && !IP->CCI[k*IP->np2+j*IP->nx2+i+1])
		label_3d_cc(KMI,label,i+1,j,k,size,CC,IP);
   if(KMI[k*IP->np2+(j+1)*IP->nx2+i]==CC && !IP->CCI[k*IP->np2+(j+1)*IP->nx2+i])
		label_3d_cc(KMI,label,i,j+1,k,size,CC,IP);
   if(KMI[k*IP->np2+(j-1)*IP->nx2+i]==CC && !IP->CCI[k*IP->np2+(j-1)*IP->nx2+i])
		label_3d_cc(KMI,label,i,j-1,k,size,CC,IP);
   if(k>0 && KMI[(k-1)*IP->np2+j*IP->nx2+i]==CC && !IP->CCI[(k-1)*IP->np2+j*IP->nx2+i])
		label_3d_cc(KMI,label,i,j,k-1,size,CC,IP);
   if(k<IP->TO-IP->FROM && KMI[(k+1)*IP->np2+j*IP->nx2+i]==CC && !IP->CCI[(k+1)*IP->np2+j*IP->nx2+i])
		label_3d_cc(KMI,label,i,j,k+1,size,CC,IP);
}

void scale_short_minmax(short *imagein, unsigned char **imageout, int np, int min,int max)
{
   float scale;
   short value;
   int tmp;

   if(min>max) 
   {
      tmp=min; min=max; max=tmp;
   }

   scale=254.99/(max-min);

   for(int i=0;i<np;i++)
   {
      value=imagein[i];

      if(value < min)
         (*imageout)[i]=0;
      else if(value > max)
         (*imageout)[i]=254;
      else
         (*imageout)[i]=(unsigned char)((value-min)*scale);
   }
}

static float costFunction1(short *KMI, float *P, struct im_params *IP)
{
   char transcode[5]="ZXYT";

	float Ax,Bx;
	float Ay,By;
	float Az,Bz;

	int i,j,k;
	int q;	

	float  	F;
	float  	*T;
	float  	*invT;

	static float FMAX=0.0;

	int	npixels;

	int jj;

	float   x,y,z;	 
	float   xx,yy,zz;
	unsigned char VAL1; 

	for(j=0;j<IP->NCC;j++) 
	{
		IP->S[j]=IP->SS[j]=0.0;
		IP->Size[j]=0;
	}

	invT=transformation(P[0],P[1],P[2],P[3],P[4],P[5], 1.,1.,1.,0,0,0,transcode);
	T = inv4(invT); free(invT);

	T[0] /= IP->dx1;
	T[1] /= IP->dx1;
	T[2] /= IP->dx1;
	T[3] /= IP->dx1;
	T[3] += IP->xc1/IP->dx1;

	T[4] /= IP->dy1;
	T[5] /= IP->dy1;
	T[6] /= IP->dy1;
	T[7] /= IP->dy1;
	T[7] += IP->yc1/IP->dy1;

	T[8] /= IP->dz1;
	T[9] /= IP->dz1;
	T[10] /= IP->dz1;
	T[11] /= IP->dz1;
	T[11] += IP->zc1/IP->dz1;

	npixels=0;
	q=0;
	for(k=IP->FROM-1;k<IP->TO;k++) 
	{
		zz=k*IP->dz2-IP->zc2;
		Bx=T[2]*zz+T[3];
		By=T[6]*zz+T[7];
		Bz=T[10]*zz+T[11];
		for(j=0;j<IP->ny2;j++) 
		{
			yy=j*IP->dy2-IP->yc2;
			Ax=T[1]*yy+Bx;
			Ay=T[5]*yy+By;
			Az=T[9]*yy+Bz;

			for(i=0;i<IP->nx2;i++) 
			{
		   		if(KMI[q]) 
				{
            		xx=i*IP->dx2-IP->xc2;

            		x=T[0]*xx+Ax;
	    	      	y=T[4]*xx+Ay;
	    	      	z=T[8]*xx+Az;

		        	VAL1=linearInterpolatorUC(x,y,z,IP->data1,IP->nx1,IP->ny1,IP->nz1,IP->np1);

		        	if(VAL1>0) 
					{
		         		npixels++;
			    		if(VAL1 > IP->max1) VAL1=IP->max1;
			    		jj=IP->CCI[q]-1;
			    		IP->S[jj] += VAL1;
			    		IP->SS[jj] += (float)(VAL1)*VAL1;
			    		IP->Size[jj]++;
		       		}
				}
		    	q++;
			}
		}
	}

	free(T);

	F=0.0;
	for(j=0;j<IP->NCC;j++) 
		if(IP->Size[j])
			F +=  ( IP->SS[j] - IP->S[j]*IP->S[j]/IP->Size[j] );

	if(F>FMAX)
		FMAX=F;

	if(npixels<IP->NP/(IP->sf*20)) 
	{
		//printf("\t%10.3f",FMAX);
		fflush(NULL);
  		return(FMAX);
	} 
	else 
	{
		//printf("\t%10.3f",F/npixels);
		fflush(NULL);
		return(F/npixels);
	}
}

static float newCostFunction(short *KMI, float *P, struct im_params *IP)
{
   char transcode[5]="ZXYT";

	float alpha;
	double mu, var;

	float Ax,Bx;
	float Ay,By;
	float Az,Bz;

	int i,j,k;
	int q;	

	float  	F;
	float  	*T;
	float  	*invT;

	static float FMAX=0.0;

	int	npixels;

	int jj;

	float   x,y,z;	 
	float   xx,yy,zz;
	float VAL1; 

	for(j=0;j<IP->NCC;j++) 
	{
		IP->S[j]=IP->SS[j]=0.0;
		IP->Size[j]=0;
	}

	invT=transformation(P[0],P[1],P[2],P[3],P[4],P[5], 1.,1.,1.,0,0,0, transcode);
	T = inv4(invT); free(invT);

	T[0] /= IP->dx1;
	T[1] /= IP->dx1;
	T[2] /= IP->dx1;
	T[3] /= IP->dx1;
	T[3] += IP->xc1/IP->dx1;

	T[4] /= IP->dy1;
	T[5] /= IP->dy1;
	T[6] /= IP->dy1;
	T[7] /= IP->dy1;
	T[7] += IP->yc1/IP->dy1;

	T[8] /= IP->dz1;
	T[9] /= IP->dz1;
	T[10] /= IP->dz1;
	T[11] /= IP->dz1;
	T[11] += IP->zc1/IP->dz1;

	npixels=0;
	q=0;
	for(k=IP->FROM-1;k<IP->TO;k++) 
	{
		zz=(k + zjit[q])*IP->dz2-IP->zc2;
		Bx=T[2]*zz+T[3];
		By=T[6]*zz+T[7];
		Bz=T[10]*zz+T[11];
		for(j=0;j<IP->ny2;j++) 
		{
			yy=(j + yjit[q])*IP->dy2-IP->yc2;
			Ax=T[1]*yy+Bx;
			Ay=T[5]*yy+By;
			Az=T[9]*yy+Bz;

			for(i=0;i<IP->nx2;i++) 
			{
		   		if(KMI[q]) 
				{
            		xx=(i + xjit[q])*IP->dx2-IP->xc2;

            		x=T[0]*xx+Ax;
	    	      	y=T[4]*xx+Ay;
	    	      	z=T[8]*xx+Az;

					VAL1 = (float)linearInterpolator(x, y, z, IP->data1, IP->nx1, IP->ny1, IP->nz1, IP->np1, &alpha);
					// VAL1 = (float)PNN(x, y, z, IP->data1, IP->nx1, IP->ny1, IP->nz1);
					// VAL1 = (float)nearestNeighbor(x, y, z, IP->data1, IP->nx1, IP->ny1, IP->nz1, IP->np1);

		        	if(VAL1>0.0) 
					{
		         		npixels++;
			    		if(VAL1 > IP->max1) VAL1=IP->max1;
			    		jj=IP->CCI[q]-1;

			    		IP->S[jj] += VAL1;
			    		IP->SS[jj] += VAL1*VAL1;
			    		IP->Size[jj]++;
		       		}
				}
		    	q++;
			}
		}
	}

	free(T);

	F=0.0;
	for(j=0;j<IP->NCC;j++) 
	if(IP->Size[j]>0)
	{
		mu = IP->S[j]/IP->Size[j];
		var =  IP->SS[j]/IP->Size[j] - mu*mu;

		if(var>0.0) F +=  (float)( IP->Size[j]*log(sqrt(var)) );
	}

	if(F>FMAX) FMAX=F;

	if(npixels<IP->NP/(IP->sf*20)) 
	{
		//printf("\t%10.3f",FMAX);
		fflush(NULL);
  		return(FMAX);
	} 
	else 
	{
		//printf("\t%10.3f",F/npixels);
		fflush(NULL);
		return(F/npixels);
	}
}

/*
static float costFunction2(short *KMI, float *P, struct im_params *IP)
{
   char transcode[5]="ZXYT";

	double sd;
	float alpha;
	float mu, E;
	double *sum1, *sum2, *sum3;
	float Ax,Bx;
	float Ay,By;
	float Az,Bz;

	int i,j,k;
	int q;	

	float  	F;
	float  	*T;
	float  	*invT;

	static float FMAX=0.0;

	int	npixels;

	int jj;

	float   x,y,z;	 
	float   xx,yy,zz;
	unsigned char VAL1; 

	sum1 = (double *)calloc( IP->NCC, sizeof(double));
	sum2 = (double *)calloc( IP->NCC, sizeof(double));
	sum3 = (double *)calloc( IP->NCC, sizeof(double));

	for(j=0;j<IP->NCC;j++) 
	{
		IP->S[j]=IP->SS[j]=0.0;
		IP->Size[j]=0;
		sum1[j]=sum2[j]=sum3[j]=0.0;
	}

	invT=transformation(P[0],P[1],P[2],P[3],P[4],P[5], 1.,1.,1.,0,0,0, transcode);
	T = inv4(invT); free(invT);

	T[0] /= IP->dx1;
	T[1] /= IP->dx1;
	T[2] /= IP->dx1;
	T[3] /= IP->dx1;
	T[3] += IP->xc1/IP->dx1;

	T[4] /= IP->dy1;
	T[5] /= IP->dy1;
	T[6] /= IP->dy1;
	T[7] /= IP->dy1;
	T[7] += IP->yc1/IP->dy1;

	T[8] /= IP->dz1;
	T[9] /= IP->dz1;
	T[10] /= IP->dz1;
	T[11] /= IP->dz1;
	T[11] += IP->zc1/IP->dz1;

	npixels=0;
	q=0;
	for(k=IP->FROM-1;k<IP->TO;k++) 
	{
		zz=k*IP->dz2-IP->zc2;
		Bx=T[2]*zz+T[3];
		By=T[6]*zz+T[7];
		Bz=T[10]*zz+T[11];
		for(j=0;j<IP->ny2;j++) 
		{
			yy=j*IP->dy2-IP->yc2;
			Ax=T[1]*yy+Bx;
			Ay=T[5]*yy+By;
			Az=T[9]*yy+Bz;

			for(i=0;i<IP->nx2;i++) 
			{
		   		if(KMI[q]) 
				{
            				xx=i*IP->dx2-IP->xc2;

            				x=T[0]*xx+Ax;
	    	      			y=T[4]*xx+Ay;
	    	      			z=T[8]*xx+Az;

		        		VAL1=linearInterpolator(x,y,z,IP->data1,IP->nx1,IP->ny1,IP->nz1,IP->np1,&alpha);

		        		if(VAL1>0) 
					{
		         			npixels++;
			    			if(VAL1 > IP->max1) VAL1=IP->max1;
			    			jj=IP->CCI[q]-1;

			    			IP->S[jj] += VAL1;
			    			IP->SS[jj] += (float)(VAL1)*VAL1;
			    			IP->Size[jj]++;

						if(alpha>0.0) 
						{
							sum1[jj] += VAL1*VAL1/alpha;
							sum2[jj] += VAL1/alpha;
							sum3[jj] += 1.0/alpha;
						}		
		       			}
				}
		    		q++;
			}
		}
	}

	free(T);

	F=0.0;
	for(j=0;j<IP->NCC;j++) 
	if(IP->Size[j]>1 && sum3[j]>0.0)
	{
		mu = sum2[j]/sum3[j];
		sd = sqrt( (sum1[j]-mu*sum2[j])/IP->Size[j] );

		E = (float)log(sd);

		F += (IP->Size[j]*E);
	}

	free(sum1);
	free(sum2);
	free(sum3);

	if(F>FMAX) FMAX=F;

	if(npixels<IP->NP/(IP->sf*20)) 
	{
		//printf("\t%10.3f",FMAX);
		fflush(NULL);
  		return(FMAX);
	} 
	else 
	{
		//printf("\t%10.3f",F/npixels);
		fflush(NULL);
		return(F/npixels);
	}
}
*/

unsigned char linearInterpolatorUC(float x, float y, float z, unsigned char *array, int nx, int ny, int nz, int np)
{
	int     i,j,k,n;
	float   u,uu;
	
	i=(int)(x);
	j=(int)(y);
	k=(int)(z);

	if(i<0 || i>(nx-2) || j<0 || j>(ny-2) )
	{
		return(0);
	}

	if( k>=0 && k<(nz-1) )
	{
		u = x - i; if(u<0.0) u=0.0;
		uu = 1.0-u;

		n=k*np + j*nx +i;
		v1 = (float)(array[n])*uu + (float)(array[n+1])*u;
		v2 = (float)(array[n+nx])*uu + (float)(array[n+nx+1])*u;
		v3 = (float)(array[n+np])*uu + (float)(array[n+np+1])*u;
		v4 = (float)(array[n+np+nx])*uu + (float)(array[n+np+nx+1])*u;

		u = y - j; if(u<0.0) u=0.0;
		uu = 1.0-u;
		w1 = v1*uu + v2*u;
		w2 = v3*uu + v4*u;

		u = z - k; if(u<0.0) u=0.0;
		return( (unsigned char)( w1*(1.0-u) + w2*u  + 0.5) );
	}

	if( k==(nz-1) )
	{
		u = x - i; if(u<0.0) u=0.0;
		uu = 1.0-u;

		n=k*np + j*nx +i;
		v1 = (float)(array[n])*uu + (float)(array[n+1])*u;

		n=k*np + (j+1)*nx +i;
		v2 = (float)(array[n])*uu + (float)(array[n+1])*u;

		u = y - j; if(u<0.0) u=0.0;
		return( (unsigned char)( v1*(1.0-u) + v2*u  + 0.5) );
	}

	return(0);
}

/* nine parameter version */
/* x,y,z translations in mm */ 
/* ax,ay,az rotations in degrees */
/* sx,sy,sz scaling parameters */
/* rX,rY,rZ reflection parameters */
/* code determines order of transformations */
float *transformation(float x, float y, float z, float ax, float ay, 
float az, float sx, float sy, float sz, int rX, int rY, int rZ, char *code)
{
	int i;
	float *Tx;    /* 4x4 matrix of rotaion about x-axis */
	float *Ty;    /* 4x4 matrix of rotaion about y-axis */
	float *Tz;    /* 4x4 matrix of rotaion about z-axis */
	float *Tt;    /* 4x4 matrix of translation */
	float *Ts;    /* 4x4 matrix of scaling */
	float *Tr;    /* 4x4 matrix of reflection */
	float *T;     /* 4x4 overall transformation matrix */
	double pi;

	pi=(float)(4.0*atan(1.0));

	Tx=(float *)calloc(16,sizeof(float));
	Ty=(float *)calloc(16,sizeof(float));
	Tz=(float *)calloc(16,sizeof(float));
	Tt=(float *)calloc(16,sizeof(float));
	Ts=(float *)calloc(16,sizeof(float));
	Tr=(float *)calloc(16,sizeof(float));
	T= (float *)calloc(16,sizeof(float));

   /* set T equal to the identity matrix first */
	T[0]=T[5]=T[10]=T[15]=1.0;

   /* set matrix of translation */
	Tt[0]=Tt[5]=Tt[10]=Tt[15]=1.0;
	Tt[3]=x;
	Tt[7]=y;
	Tt[11]=z;

   /* set matrix of scaling */
	Ts[0]=sx;
	Ts[5]=sy;
	Ts[10]=sz;
	Ts[15]=1.0;

   /* set matrix of scaling */
	if(rX) Tr[0]=-1.0; else Tr[0]=1.0;
	if(rY) Tr[5]=-1.0; else Tr[5]=1.0;
	if(rZ) Tr[10]=-1.0; else Tr[10]=1.0;
	Tr[15]=1.0;

   /* set matrix of rotaion about x-axis */
	Tx[0]=Tx[15]=1.0;
	Tx[5]=Tx[10]=(float)cos(ax*pi/180.0);
	Tx[6]=-(float)sin(ax*pi/180.0);
	Tx[9]=(float)sin(ax*pi/180.0);

   /* set matrix of rotaion about y-axis */
   Ty[5]=Ty[15]=1.0;
	Ty[0]=Ty[10]=(float)cos(ay*pi/180.0);
	Ty[8]=-(float)sin(ay*pi/180.0);
	Ty[2]=(float)sin(ay*pi/180.0);

   /* set matrix of rotaion about z-axis */
	Tz[10]=Tz[15]=1.0;
	Tz[0]=Tz[5]=(float)cos(az*pi/180.0);
	Tz[1]=-(float)sin(az*pi/180.0);
	Tz[4]=(float)sin(az*pi/180.0);

	i=0;
	while( code[i] != '\0' ) {
		switch(code[i++])	{
			case 'X': {
 				/* set T=Tx*T */
				multi(Tx,4,4,T,4,4,T);
				break;
			}
			case 'Y': {
 				/* set T=Ty*T */
				multi(Ty,4,4,T,4,4,T);
				break;
			}
			case 'Z': {
 				/* set T=Tz*T */
				multi(Tz,4,4,T,4,4,T);
				break;
			}
			case 'T': {
 				/* set T=Tt*T */
				multi(Tt,4,4,T,4,4,T);
				break;
			}
			case 'S': {
 				/* set T=Ts*T */
				multi(Ts,4,4,T,4,4,T);
				break;
 			}
			case 'R': {
 				/* set T=Tr*T */
				multi(Tr,4,4,T,4,4,T);
				break;
 			}
		}
	}

	free(Tx);
	free(Ty);
	free(Tz);
	free(Tt);
	free(Ts);
	free(Tr);

	return(T);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Coverts angular rotations (ax,ay,az) and translations (x,y,z) to a rigid-body transformation matrix (T)
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void set_transformation(float x, float y, float z, float ax, float ay, float az, const char *code, float *T)
{
   int i;
   float Tx[16];    /* 4x4 matrix of rotaion about x-axis */
   float Ty[16];    /* 4x4 matrix of rotaion about y-axis */
   float Tz[16];    /* 4x4 matrix of rotaion about z-axis */
   float Tt[16];    /* 4x4 matrix of translation */
   double pi;

   pi=(float)(4.0*atan(1.0));

   /* set T equal to the identity matrix first */
   set_to_I(T,4);

   /* set matrix of translation */
   set_to_I(Tt,4);
   Tt[3]=x;
   Tt[7]=y;
   Tt[11]=z;

   /* set matrix of rotaion about x-axis */
   set_to_I(Tx,4);
   Tx[5]=Tx[10]=(float)cos(ax*pi/180.0);
   Tx[6]=-(float)sin(ax*pi/180.0);
   Tx[9]=(float)sin(ax*pi/180.0);

   /* set matrix of rotaion about y-axis */
   set_to_I(Ty,4);
   Ty[0]=Ty[10]=(float)cos(ay*pi/180.0);
   Ty[8]=-(float)sin(ay*pi/180.0);
   Ty[2]=(float)sin(ay*pi/180.0);

   /* set matrix of rotaion about z-axis */
   set_to_I(Tz,4);
   Tz[0]=Tz[5]=(float)cos(az*pi/180.0);
   Tz[1]=-(float)sin(az*pi/180.0);
   Tz[4]=(float)sin(az*pi/180.0);

   i=0;
   while( code[i] != '\0' ) {
      switch(code[i++])	{
         case 'X': {
            /* set T=Tx*T */
            multi(Tx,4,4,T,4,4,T);
            break;
         }
         case 'Y': {
            /* set T=Ty*T */
            multi(Ty,4,4,T,4,4,T);
            break;
         }
         case 'Z': {
            /* set T=Tz*T */
            multi(Tz,4,4,T,4,4,T);
            break;
         }
         case 'T': {
            /* set T=Tt*T */
            multi(Tt,4,4,T,4,4,T);
            break;
         }
      }
   }
}
/////////////////////////////////////////////////////////////////////////////////////////////////

static float Gradient_Descent(short *KMI, int	ndim, float ftol, struct im_params *IP)
{
	float 	PT[6],XIT[6];
	int	iter=0;  /* number of iterations taken */
	int i,j;
	float	fret,fp;
	float   sum;

	fp=(*obj_fnc)(KMI,P,IP); 

	global_min=fp;

	do {
		iter++;

		for(j=0;j<ndim;j++)
			PT[j]=P[j];

		printf("\nIteration %d: ",iter);

		PT[0]=P[0]+0.25 ;
		XIT[0]=fp-(*obj_fnc)(KMI,PT,IP);

		if( (fp-XIT[0]) < global_min )
		{
			global_min = (fp-XIT[0]);
			for(i=0;i<6;i++) Pmin[i]=P[i];
		}

		PT[0]=P[0];
		PT[1]=P[1]+0.25 ;
		XIT[1]=fp-(*obj_fnc)(KMI,PT,IP);

		if( (fp-XIT[1]) < global_min )
		{
			global_min = (fp-XIT[1]);
			for(i=0;i<6;i++) Pmin[i]=P[i];
		}

		PT[1]=P[1];
		PT[2]=P[2]+0.25 ;
		XIT[2]=fp-(*obj_fnc)(KMI,PT,IP);

		if( (fp-XIT[2]) < global_min )
		{
			global_min = (fp-XIT[2]);
			for(i=0;i<6;i++) Pmin[i]=P[i];
		}

		PT[2]=P[2];
		PT[3]=P[3]+0.25 ;
		XIT[3]=fp-(*obj_fnc)(KMI,PT,IP);

		if( (fp-XIT[3]) < global_min )
		{
			global_min = (fp-XIT[3]);
			for(i=0;i<6;i++) Pmin[i]=P[i];
		}

		PT[3]=P[3];
		PT[4]=P[4]+0.25 ;
		XIT[4]=fp-(*obj_fnc)(KMI,PT,IP);

		if( (fp-XIT[4]) < global_min )
		{
			global_min = (fp-XIT[4]);
			for(i=0;i<6;i++) Pmin[i]=P[i];
		}

		PT[4]=P[4];
		PT[5]=P[5]+0.25 ;
		XIT[5]=fp-(*obj_fnc)(KMI,PT,IP);

		if( (fp-XIT[5]) < global_min )
		{
			global_min = (fp-XIT[5]);
			for(i=0;i<6;i++) Pmin[i]=P[i];
		}

		sum=0.0;
		for(i=0;i<ndim;i++)
			sum+= (XIT[i]*XIT[i]);

		sum=(float)sqrt( (double)sum);				
		//printf("Search Direction: "); 
		for(i=0;i<ndim;i++) {
			XIT[i]/=sum;
			//printf("   %5.2f",XIT[i]);
		}
		//printf("\n");

		fret=minimize1D(KMI,XIT,ndim,IP);

		if( fret < global_min )
		{
			global_min = fret;
			for(i=0;i<6;i++) Pmin[i]=P[i];
		}
		else 
		{
			fret = global_min;
			for(i=0;i<6;i++) P[i]=Pmin[i];
		}

/*
		printf("Minimum Point: \n");
		printf("   %5.2f mm",P[0]);
		printf("   %5.2f mm",P[1]);
		printf("   %5.2f mm",P[2]);
		printf("   %5.2f deg",P[3]);
		printf("   %5.2f deg",P[4]);
		printf("   %5.2f deg",P[5]);
*/

		printf("Minimum = %f\n",fret);

		// write(registerer.fd[1],Pmin,sizeof(float)*6);
		// kill(getppid(),SIGUSR1);

		if(2.0*fabs((double)fp-fret)<= ftol*(fabs((double)fp)+fabs((double)fret)))
		break;
		
		fp=fret;

	} while(1);

	return(fret);
}

static float minimize1D(short *KMI, float *X, int ndim, struct im_params *IP)
{
	float ax,xx,bx,fa,fx,fb;
	float fret;
	float tol,xmin;
	int i;

	/* tol=0.0001; */
	// tol=0.1;  4/10/03
	tol=0.001;

	ax=0.0;
	xx=2.0;

	ncom=ndim;
	for(i=0;i<ndim;i++) {
		PCOM[i]=P[i];
		XICOM[i]=X[i];
	}

	findInterval(KMI,&ax,&xx,&bx,&fa,&fx,&fb,IP);

	fret=optimize(KMI,ax,xx,bx,tol,&xmin,IP);
	
	for(i=0;i<ndim;i++) { 
		X[i]*=xmin;
		P[i]+=X[i];
	}

	return(fret);
}

static float optimize(short *KMI, float ax,float bx,float cx,float tol, float *xmin, struct im_params *IP)
{
	float cgold=0.381966;
	float zeps=1.0e-10;
	int itmax=100;
	float a,b,v,w,x,u,fu,fx,fv,fw,e=0.0,xm,d,r,q,p,etemp;
	int i;
	float tol1,tol2;

	if(ax>cx)  {
		b=ax; a=cx;
	} else {
		a=ax; b=cx;
	}

	w=x=v=bx;
	e=0.0;
	fv=fw=fx=f1dim(KMI,x,IP);

	for(i=0;i<itmax;i++) {
		xm=0.5*(a+b);
		tol1=tol*(float)fabs((double)x)+zeps;
		tol2=2. * tol1;

		if( fabs((double)x-xm) <= (tol2-0.5*(b-a)))
			break;

		if( fabs((double)e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;	
			q=2.0*(q-r);
			if(q>0.0) p=-p;
			q=(float)fabs((double)q);
			etemp=e;
			e=d;
			if( fabs((double)p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x) )
				goto ONE;
			d=p/q;
			u=x+d;
			if( (u-a)<tol2 || (b-u)<tol2) 
				d= ( (xm-x)>0.0 ) ? tol1 : -tol1;
			goto TWO;
		}

		ONE:
		if( x >= xm)
			e=a-x;
		else
			e=b-x;
		d=cgold*e;

		TWO:
		if(fabs((double)d) >= tol1)
			u=x+d;
		else
			u= x + ( (d>0) ? tol1 : -tol1 );

		fu=f1dim(KMI,u,IP);
		if(fu<=fx) {
			if(u>=x)
				a=x;
			else
				b=x;
			v=w;
			fv=fw;
			w=x;
			fw=fx;
			x=u;
			fx=fu;
		} else {
			if(u<x)
				a=u;
			else
				b=u;
			if(fu<=fw || w==x) {
				v=w;
				fv=fw;
				w=u;
				fw=fu;
			} else if(fu<=fv || v==x || v==w) {
				v=u;
				fv=fu;
			}
		}
	}
	if(i==itmax-1)
		printf("\nExceeded maximum iteratons");

	*xmin=x;
	return(fx);
}

static float f1dim(short *KMI, float x, struct im_params *IP)
{
	int i;
	float   XT[6];

	for(i=0;i<ncom;i++) 
		XT[i]=PCOM[i]+x*XICOM[i];
	
	return((*obj_fnc)(KMI,XT,IP)); 
}

static void  findInterval(short *KMI, float *ax,float *bx,float *cx, float *fa, float *fb, float *fc, struct im_params *IP)
{
	float glimit=100.0;

	float dum;
	float r,q,u,ulim;
	float fu;

	*fa=f1dim(KMI,*ax,IP);
	*fb=f1dim(KMI,*bx,IP);

	if( *fb > *fa) {
		dum=*ax;
		*ax=*bx;
		*bx=dum;
		dum=*fb;
		*fb=*fa;
		*fa=dum;
	}

	*cx= *bx + 1.618034*(*bx - *ax);
	*fc=f1dim(KMI,*cx,IP);

	while(*fb >= *fc) {	
		r = (*bx - *ax) * (*fb - *fc);	
		q = (*bx - *cx) * (*fb - *fa);	
		dum=(float)fabs((double)q-r);
		dum= (dum>1.0e-20) ? dum : 1.0e-20;
		dum= ( (q-r)>=0 ) ? dum : -dum;
		u = *bx - ((*bx-*cx)*q-(*bx-*ax)*r)/(2.*dum);
		ulim = *bx + glimit*(*cx - *bx);
		if( (*bx - u)*(u - *cx) > 0.0 ) {
			fu=f1dim(KMI,u,IP);
			if(fu < *fc) {
				*ax=*bx;
				*fa=*fb;
				*bx=u;
				*fb=fu;
				continue;
			} else if( fu > *fb) {
				*cx=u;
				*fc=fu;
				continue;
			}
			u=*cx+1.618034*(*cx-*bx);
			fu=f1dim(KMI,u,IP);
		} else if( (*cx - u)*(u-ulim) > 0.0 ) {
			fu=f1dim(KMI,u,IP);
			if(fu < *fc) {
				*bx=*cx;
				*cx=u;
				u=*cx + 1.618034*(*cx - *bx);
				*fb=*fc;
				*fc=fu;
				fu=f1dim(KMI,u,IP);
			}
		} else if( (u-ulim)*(ulim- *cx) >= 0.0 ) {
			u = ulim;
			fu=f1dim(KMI,u,IP);
		} else {
			u = *cx + 1.618034*(*cx - *bx);
			fu=f1dim(KMI,u,IP);
		}
		*ax=*bx;
		*bx=*cx;
		*cx=u;
		*fa=*fb;
		*fb=*fc;
		*fc=fu;
	}
}

void testCostFunc1(short *trg, int Tnx, int Tny, int Tnz, float Tdx, float Tdy, float Tdz, 
short *obj, int Onx, int Ony, int Onz, float Odx, float Ody, float Odz)
{
	FILE *fp;
	float cost;
	int thresh;
	short *ccImage;
	short *im_out;
	int size_thresh=100;
	int OL,OH;
	//float *T;

	int Tnp,Tnv,Onv;
	
	for(int i=0; i<6; i++) P[i]=0.0;

	IP.FROM = 1;
	IP.TO = Tnz;

	IP.nx2=Tnx;
	IP.ny2=Tny;
	IP.nz2=Tnz;
	IP.dx2=Tdx;
	IP.dy2=Tdy;
	IP.dz2=Tdz;

	IP.nv2=IP.nx2*IP.ny2*IP.nz2;
	IP.np2=IP.nx2*IP.ny2;

	IP.xc2=IP.dx2*(IP.nx2-1)/2.0; /* +---+---+ */
	IP.yc2=IP.dy2*(IP.ny2-1)/2.0;
	IP.zc2=IP.dz2*(IP.nz2-1)/2.0;
	
	IP.sf = 1;

	IP.nx1=Onx;
	IP.ny1=Ony;
	IP.nz1=Onz;
	IP.dx1=Odx;
	IP.dy1=Ody;
	IP.dz1=Odz;

	IP.nv1=IP.nx1*IP.ny1*IP.nz1;
	IP.np1=IP.nx1*IP.ny1;

	IP.xc1=IP.dx1*(IP.nx1-1)/2.0; /* +---+---+ */
	IP.yc1=IP.dy1*(IP.ny1-1)/2.0;
	IP.zc1=IP.dz1*(IP.nz1-1)/2.0;

	Tnv = Tnx*Tny*Tnz;
	Tnp = Tnx*Tny;
	Onv = Onx*Ony*Onz;

	thresh=findThresholdLevel(trg, Tnv);
	printf("Threshold level = %d\n",thresh);

	ccImage=(short *)calloc(Tnv,sizeof(short));

	/* thresholding */
	for(int i=0;i<Tnv;i++)
	if( trg[i] <= thresh)
		ccImage[i]=0;
	else
		ccImage[i]=1;

	cca(ccImage, Tnx, Tny, Tnz);

	printf("Number of classes = %d\n",nclass);
	printf("Maximum number of iterations allowed = %d\n",maxiter);

	printf("\nK-means clustering ...\n");
	im_out=KMcluster(ccImage, trg, nclass, maxiter, thresh, low, high , Tnv);

	free(ccImage);
	ccImage=im_out;

	// make sure edges are zero
	for(int k=0; k<Tnz; k++)
	{
		for(int i=0; i<Tnx; i++) im_out[k*Tnp+0*Tnx+i]=0;
		for(int i=0; i<Tnx; i++) im_out[k*Tnp+(Tny-1)*Tnx+i]=0;
		for(int j=0; j<Tny; j++) im_out[k*Tnp+j*Tnx+0]=0;
		for(int j=0; j<Tny; j++) im_out[k*Tnp+j*Tnx+(Tnx-1)]=0;
	}

	fprintf(stdout,"\nCluster CC analysis ...\n");
	IP.CCI=(unsigned *)calloc(Tnv,sizeof(unsigned int));
	IP.NCC=label_CCI(ccImage,size_thresh,&IP,Tnv);

	fprintf(stdout,"\nCost function minimization ...\n");
	IP.Size=(int *)calloc(IP.NCC,sizeof(int));
	IP.S=(float *)calloc(IP.NCC,sizeof(float));
	IP.SS=(float *)calloc(IP.NCC,sizeof(float));

	IP.data1=(unsigned char *)malloc(IP.nv1);
	setLowHigh(obj, Onv, &OL, &OH);
	scale_short_minmax( obj, &(IP.data1), IP.nv1, OL, OH);
	IP.max1=255;

	initializeRandomNumberGenerator();

	xjit=(float *)calloc(Tnv,sizeof(float));
	yjit=(float *)calloc(Tnv,sizeof(float));
	zjit=(float *)calloc(Tnv,sizeof(float));

/*
*/
	for(int i=0; i<Tnv; i++)
	{
		xjit[i] = (float) (drand48() - 0.5);
		yjit[i] = (float) (drand48() - 0.5);
		zjit[i] = (float) (drand48() - 0.5);
	}

	fp = fopen("tt","w");
    if(fp==NULL) file_open_error("tt");

	for(int i=-40; i<=40; i++)
	// for(int i=-240; i<=240; i++)
	{

		P[0]=P[1]=P[2]=P[3]=P[4]=P[5]=0.0;

		// P[0]=i*1.0/100.0;
		P[5]=i*1.0/20.0;
		cost=newCostFunction(ccImage,P,&IP); 
	
		// fprintf(fp,"%f %f\n",i*1.0/100.0, cost);	
		fprintf(fp,"%f %f\n",i*1.0/20.0, cost);	
	}
	fclose(fp);

//	printf("\n%f %f\n",alpha_min, alpha_max);

	free(ccImage);
	free(IP.data1);
	free(IP.Size);
	free(IP.S);
	free(IP.SS);
	free(IP.CCI);
	free(xjit); free(yjit); free(zjit);
}
