#include <stdlib.h>
//#include <malloc.h>  
#include <math.h>
#include <strings.h>
#include <string.h>		// required by strlen()
#include <stdio.h>
#include <sys/types.h>  //      required by time(), stat()
#include <time.h>       //      required by time()
#include <sys/stat.h>   //      required by stat() 
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include "babak_lib.h"

#define YES 1
#define NO 0

char prefix[512];
int opt;

static struct option options[] =
{
	{"-d", 1, 'd'},
	{"-c", 1, 'c'},
	{"-o", 1, 'o'},
	{0, 0, 0}
};

int opt_d=NO;
int opt_c=NO;
int opt_o=NO;

void print_help_and_exit()
{
	printf("\n\nUsage: ccstats -d <image list> -c <cc image> -o <output prefix>\n\n");
	exit(0);
}

int numberOfRows(const char *filename, int nc)
{
	FILE *fp;
	int n=0;
	char s[512];

	if(nc<=0) return(0);

	fp = fopen(filename,"r");
	if(fp == NULL) file_open_error(filename);

	while( fscanf(fp, "%s", s) != EOF ) n++;

	fclose(fp);

	printf("\nNumber of images = %d\n",n/nc);

	if( (n/nc) <= 1)
	{
		printf("\n\nNumber of rows of the data matrix must be greater than 1, aborting ...\n\n");
		exit(0);
	}

	return(n/nc);
}

char **read_idata(const char *dataFile, const char *dataTypeCode, const char *dataMaskCode, int nr, int nc, int pi)
{
	char **idata=NULL;
	FILE *fp;
	char s[512];
	int count=0;

	if( (nr*pi) > 0)
	{
		// memory allocation
		idata = (char **)calloc(nr*pi,sizeof(char *));
    	if(idata==NULL) memory_allocation_error("idata");
		for(int i=0; i<nr*pi; i++) 
		{
			idata[i] = (char *)calloc(512,sizeof(char));
			if(idata[i]==NULL) memory_allocation_error("idata[]");
		}
	}
	else	return(NULL);

	fp = fopen(dataFile,"r");
	if(fp == NULL) file_open_error(dataFile);

	count=0;
	for(int i=0; i<nr; i++)
	for(int j=0; j<nc; j++)
	{
		if( dataTypeCode[j]=='i' && dataMaskCode[j]=='1' )
			fscanf(fp, "%s", idata[count++]);
		else
			fscanf(fp, "%s", s);
	}

	fclose(fp);

	printf("\nImage files:\n");
	for(int i=0; i<nr; i++) 
	{
		for(int j=0; j<pi; j++) 
			printf("%s\t",idata[i*pi + j]);
		printf("\n");
	}

	return(idata);
}

void checkDimension_vancova(int N, char **imagefile, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz)
{
   nifti_1_header hdr;
   short dataType;

   if(N==0) return;

   hdr = read_NIFTI_hdr(imagefile[0]);
   *nx = hdr.dim[1];
   *ny = hdr.dim[2];
   *nz = hdr.dim[3];
   *dx = hdr.pixdim[1];
   *dy = hdr.pixdim[2];
   *dz = hdr.pixdim[3];
   dataType = hdr.datatype;

   if(dataType != 4)
   {
      printf("\n\ncheckDimension_vancova(): This program cannot handle NIFTI image data type %d, aborting ...\n",dataType);
      exit(0);
   }

   if(opt_v)
   {
      printf("\nImage matrix size (nx,ny,nz) = %d x %d x %d",*nx,*ny,*nz);
      printf("\nImage voxel size (dx,dy,dz)= %f x %f x %f\n",*dx,*dy,*dz);
   }

   for(int i=1; i<N; i++)
   {
      hdr = read_NIFTI_hdr(imagefile[i]);

      if( *nx != hdr.dim[1] ||  *ny != hdr.dim[2] ||  *nz != hdr.dim[3]) 
      {
            printf("\n\nImage %d: %s",i+1,imagefile[i]);
            printf("\n\tMatrix size = %d x %d x %d", hdr.dim[1], hdr.dim[2], hdr.dim[3]);
            printf("\n\nAll input images must be of size: %d x %d x %d\n\n",*nx,*ny,*nz);
            exit(0);
      }

      if(dataType != 4)
      {
            printf("\n\ncheckDimension_vancova(): This program cannot handle NIFTI image data type %d, aborting ...\n",dataType);
            exit(0);
      }
   }
}

// N is the number of images
void save_cluster_avg(float *spm,int nx,int ny,int nz, int N, char **imlist)
{
   nifti_1_header hdr;
   short *im;
   char ccfile[1024];
   char *temp_im;
   int np, nv; 
   int nvox; // total number of voxels in all CC's
   int *L;
   int *S;
   int ncc; // number of CC's
   int CCsize;
   FILE *fp;

   float *ccavg;
   short *ccmax;
   short *ccmin;
   float *ccsd;

   float min,max,avg;
   float iavg, javg, kavg;
   int imax=0,jmax=0,kmax=0;
   int imin=0,jmin=0,kmin=0;

   nv=nx*ny*nz;
   np=nx*ny;

   //////////////////////////////////////////////////////////////////
   temp_im = (char *)calloc(nv,1);
   nvox=0;
   for(int i=0; i<nv; i++) 
   if(spm[i]!=0.0) 
   {
      temp_im[i]=1;  
      nvox++;
   } else temp_im[i]=0;

   Connected_Component_location(temp_im, nx, ny, nz, &ncc, &S, &L);
   free(temp_im);
   //////////////////////////////////////////////////////////////////

   ccavg = (float *)calloc(N*ncc, sizeof(float));  // (N x ncc) matrix
   ccmax = (short *)calloc(N*ncc, sizeof(short));  // (N x ncc) matrix
   ccmin = (short *)calloc(N*ncc, sizeof(short));  // (N x ncc) matrix
   ccsd = (float *)calloc(N*ncc, sizeof(float));   // (N x ncc) matrix
   for(int i=0; i<N*ncc; i++) 
   {
      ccavg[i] = 0.0;
      ccmax[i] = 0;
      ccmin[i] = 0;
      ccsd[i] = 0.0;
   }

   for(int i=0; i<N; i++) 		// i is image index
   {
      im = (short *)read_nifti_image(imlist[i], &hdr);

      for(int j=0; j<ncc; j++)	// j is CC index
      {
         // determine size of the connected component from array S
         if( j != (ncc-1) ) CCsize=S[j+1]-S[j]; else CCsize=nvox-S[j];

         for(int k=S[j]; k<S[j]+CCsize; k++) 
         {
            if(k==S[j]) // set both ccmin and ccmax to the value of the first voxel
            {
               ccmin[i*ncc + j] = ccmax[i*ncc + j] = im[ L[k] ]; 
            }
            else if( im[ L[k] ] > ccmax[i*ncc + j] )
            {
               ccmax[i*ncc + j] = im[ L[k] ];
            }
            else if( im[ L[k] ] < ccmin[i*ncc + j] )
            {
               ccmin[i*ncc + j] = im[ L[k] ];
            }

            ccavg[i*ncc + j] += im[ L[k] ];
         }

         if(CCsize>0) ccavg[i*ncc + j] /= CCsize;

         for(int k=S[j]; k<S[j]+CCsize; k++) 
         {
            ccsd[i*ncc + j] += (im[ L[k] ] - ccavg[i*ncc + j])*(im[ L[k] ] - ccavg[i*ncc + j]);
         }
     
         if(CCsize == 1) 
         {
            ccsd[i*ncc + j] = sqrtf( ccsd[i*ncc + j] );
         } 
         else 
         {
            ccsd[i*ncc + j] = sqrtf( ccsd[i*ncc + j]/(CCsize-1.0) );
         }
      }

      free(im);
   }



	for(int j=0; j<ncc; j++)	// j is CC index
	{
        sprintf(ccfile,"%s_cc_stats%d.csv",prefix,j);
        fp = fopen(ccfile,"w");

		// determine size of the connected component from array S
		if( j != (ncc-1) ) CCsize=S[j+1]-S[j]; else CCsize=nvox-S[j];

		min = max = spm[ L[ S[j]] ];
		avg = 0.0;
		iavg=javg=kavg=0.0;

		for(int k=S[j]; k<S[j]+CCsize; k++)
		{
			if(spm[L[k]] <= min ) 
			{
				min = spm[L[k]];
				imin = (L[k]%np)%nx;
				jmin = (L[k]%np)/nx;
				kmin = (L[k]/np);
			}

			if(spm[L[k]] >= max ) 
			{
				max = spm[L[k]];
				imax = (L[k]%np)%nx;
				jmax = (L[k]%np)/nx;
				kmax = (L[k]/np);
			}

			avg += spm[L[k]];

			iavg += (L[k]%np)%nx;
			javg += (L[k]%np)/nx;
			kavg += (L[k]/np);
		}
		if(CCsize != 0) 
		{
			avg /= CCsize;
			iavg /= CCsize;
			javg /= CCsize;
			kavg /= CCsize;
		}
		
      printf("Cluster number %d:\n", j+1);
      printf("size=%d\n",CCsize);
      //printf("minimum statistical parameter value=%f\n",min);
      //printf("maximum statistical parameter value=%f\n",max);
      printf("average statistical parameter value=%f\n",avg);
      //printf("centroid (i,j,k)=(%.2f,%.2f,%.2f)\n",iavg,javg,kavg);
      printf("peak (i,j,k)=(%d,%d,%d)\n",imax,jmax,kmax);
      printf("trough (i,j,k)=(%d,%d,%d)\n",imin,jmin,kmin);

      fprintf(fp,"Cluster number %d:\n", j+1);
      fprintf(fp,"size=%d\n",CCsize);
      //fprintf(fp,"minimum statistical parameter value=%f\n",min);
      //fprintf(fp,"maximum statistical parameter value=%f\n",max);
      fprintf(fp,"average statistical parameter value=%f\n",avg);
      //fprintf(fp,"centroid (i,j,k)=(%.2f,%.2f,%.2f)\n",iavg,javg,kavg);
      fprintf(fp,"peak (i,j,k)=(%d,%d,%d)\n",imax,jmax,kmax);
      fprintf(fp,"trough (i,j,k)=(%d,%d,%d)\n",imin,jmin,kmin);

      //fprintf(fp,"Cluster avg,min,max,SD:\n");
      for(int i=0; i<N; i++) 
      {
         fprintf(fp,"%s,%f\n",imlist[i],ccavg[i*ncc+j]);
         //fprintf(fp,"%s,%f,%d,%d,%f\n",imlist[i],ccavg[i*ncc+j],ccmin[i*ncc+j],ccmax[i*ncc+j],ccsd[i*ncc+j]);
      }
      fclose(fp);
   }


   free(ccavg);
   free(ccmax);
   free(ccmin);
   free(ccsd);
   free(L);
   free(S);
}

int main(int argc, char **argv)
{
   nifti_1_header hdr;

	char dataFile[512];
	char ccFile[512];
	int n;	// number of images
	char **imlist=NULL;
	int nx,ny,nz,nv;
	float dx,dy,dz;
	float *ccim;

	while( (opt=getoption(argc, argv, options)) != -1)
	{
		switch (opt) {
			case 'd':
				sprintf(dataFile,"%s",optarg);
				opt_d=YES;
				break;
			case 'c':
				sprintf(ccFile,"%s",optarg);
				opt_c=YES;
				break;
			case 'o':
				sprintf(prefix,"%s",optarg);
				opt_o=YES;
				break;
			case '?':
				print_help_and_exit();
		}
	}

	if(!opt_d || !opt_c) print_help_and_exit();
	if(!opt_o) sprintf(prefix,"cc");

	printf("\nImage data file = %s\n",dataFile);
	printf("\nConnected components image = %s\n",ccFile);

    n = numberOfRows(dataFile,1);

	imlist = read_idata(dataFile, "i", "1", n, 1, 1);

   checkDimension_vancova(n, imlist, &nx, &ny, &nz, &dx, &dy, &dz);
   nv = nx*ny*nz;

   // ccim=(float *)read_analyze_image(ccFile, &nx0, &ny0, &nz0, &dx, &dy, &dz, &type, 0);
   {
       char *tmp;

       tmp = (char *)read_nifti_image(ccFile, &hdr);

	   if(nx!=hdr.dim[1] || ny!=hdr.dim[2]|| nz!=hdr.dim[3])
	   {
	   	  printf("\nIncompatible matrix dimensions between the input images and the clusters image, aborting ....\n");
		  exit(0);
	   }

	   if(hdr.datatype == 16)
       {
          ccim = (float *)tmp;
       }
       else if(hdr.datatype==4 || hdr.datatype==512 )
       {
          short *tmps;    
          tmps = (short *)tmp;
          ccim = (float *)calloc(nx*ny*nz, sizeof(float));
          for(int i=0; i<nx*ny*nz; i++) ccim[i]=tmps[i];
          free(tmp);
       }
       else
	   {
	   	   printf("\nThe clusters image is not of type 4, 16 or 512, aborting ...\n");
		   exit(0);
	   }
    }

	save_cluster_avg(ccim, nx, ny, nz, n, imlist);

   free(ccim);    
}
