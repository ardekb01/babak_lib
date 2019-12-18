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
#include "spm_analyze.h"
#include <f2c.h>
#include "babak_lib.h"
#include "niftiimage.h"
#include <minmax.h>
#include <stats.h>

#define YES 1
#define NO 0

extern "C" int cdft_(integer *which,doublereal *P,doublereal *Q, doublereal *t, doublereal *df,integer *status,integer *bound);
extern "C" int cdff_(integer *which,doublereal *P,doublereal *Q, doublereal *F, doublereal *dfn, doublereal *dfd,integer *status,integer *bound);

FILE *logFilePtr;

int opt;

static struct option options[] =
{
	{"-i", 1, 'i'}, //new
	{"-d", 1, 'd'},
	{"-p", 1, 'p'},
	{"-rw", 1, 'r'},

	{"-nozero", 0, 'n'},

	{"-v", 0, 'v'},
	{"-verbose", 0, 'v'},

	{"-data", 1, 'd'},

	{"-dataType", 1, 'T'},
	{"-dataMask", 1, 'M'},

   {"-m", 1, 'm'},
   {"-mask", 1, 'm'},

	{"-output", 1, 'o'},
	{"-o", 1, 'o'},

	{"-c", 1, 'c'},
	{"-contrast", 1, 'c'},

	{"-FWHM", 1, 'F'},

	{"-x", 1, 'x'},
	{"-y", 1, 'y'},
	{"-z", 1, 'z'},
	{0, 0, 0}
};

int opt_nozero=NO;
int opt_d=NO;
int opt_dataType=NO;
int opt_dataMask=NO;
int opt_m=NO;
int opt_o=NO;
int opt_c=NO;

int opt_x=NO;
int opt_y=NO;
int opt_z=NO;
int opt_FWHM=NO;

void print_help_and_exit()
{
	printf("\n\nUsage: vancova [-v/-verbose] [-dataMask <dataMaskCode>] [-FWHM <FWHM>]\n"
	"[-o/-output <prefix>] [-m/-mask mask] [-x <column>] [-y <row>] [-z <slice>] [-nozero]\n"
	"-d/-data <dataFile> -dataType <dataTypeCode> -c/-contrast <contrastFile>\n\n");

	printf("This program performs voxelwise analysis of covariance on image data.\n\n");

	printf("Required arguments:\n\n");

	printf("-d or -data <dataFile>\n"
	"\tSpecifies a text file, <dataFile>, containing a data table used in the ANCOVA analyses.\n"
	"\tColumns of this table can be either numberical or list image files. However, the first\n"
	"\tcolumn is reserved for the dependent variable, and must contain images.  All images\n"
	"\tare expected to be in NIFTI format and of type `short', be spatially registered, and\n"
	"\thave equal matrix size and number of slices.\n\n"
	"\tThe program runs a voxelwise ANCOVA, where at each voxel the image columns are replaced\n"
	"\tby the image values at the voxel being analyzed. The advantage of this program over\n"
	"\tsimilar programs (e.g., SPM) is that it allows the linear model to change from one voxel\n"
	"\tto another.\n\n");

	printf("-dataType <dataTypeCode>\n"
	"\t<dataTypeCode> is a string exclusively of characters `i' or `n' that specifies the type of data\n"
	"\t(numerical or image list) contained in the <dataFile>, where `i' signifies image list columns\n"
	"\tand `n' signifies numerical columns. For example, code \"innnni\" would indicate that <dataFile>\n"
	"\thas 6 columns.  The first and 6th columns are image lists, and columns 2-5 are numerical.\n\n");

	printf("-c or -contrast <contrastFile>\n"
	"\tSpecifies a text file containing the contrast being tested.  The numbers in <contrastFile>\n"
	"\tshould correspond to the columns in the <dataFile>.  That is, one number fo reach column.\n"
	"\tWhen a zero is specified in the contrast file, it means that the corresponding independent\n"
	"\tvariable is controlled for or covaried out in the ANCOVA analyses.\n\n");

	printf("Optional arguments:\n\n");

	printf("-dataMask <dataMaskCode>\n"
	"\t<dataMaskCode> is a string exclusively of charaters `0' or `1'. The length of the string must\n"
	"\tequal the number of columns in the <dataFile>.  If a 0 is specified in a given position in\n"
	"\t<dataMaskCode>, the corresponding column is ignored and not used in the ANCOVA analyses.\n"
	"\tThe default value for this argument is \"11...1\", that is, all codes are given the value\n"
	"\tof 1. This means that by default, all columns will be used in the analyses.  For example,\n"
	"\ta code \"101101\" would instruct the program to ignore columns 2 and 5 of the <dataFile> in\n"
	"\tthe ANCOVA analyses.\n\n");

	printf("-FWHM <FWHM>\n"
	"\tSpecifies the full width at half maximum (in millimeters) of a Gaussian filter that is\n"
	"\tapplied to all input images before analyses are performed.  Default FWHM=0.0\n\n");

	printf("-nozero\n"
	"\tWhen this option is selected, voxels where the dependent variable takes on a value of zero\n"
	"\tare not analyzed, that is, the t-value is set to zero at those voxels.\n\n");

	printf("-m or -mask <mask>\n"
	"\tThis option specifies a NIFTI image of type `short'.  Voxels where the <mask> image\n"
	"\tvalue is zero are not used in the ANCOVA analyses.  If no mask is specified, all\n"
	"\timage voxels are analyzed.\n\n");

	printf("-o or -output <prefix>\n"
	"\tThis option can be used to specify a prefix for the output images.  The default value is\n"
	"\t<prefix>=V.\n\n");

	printf("-v or -verbose\n"
	"\tRuns the program in verbose mode.\n\n");

	printf("-x <column> -y <row> -z <slice>\n"
	"\tThese options are intended for testing the software.  When specified, a table is created\n"
	"\tand written to the log file.  This table can be ported into other statistical programs\n"
	"\t(e.g., SPSS) in order to validate the results of the current program at the specified voxel.\n\n");

	printf("Outputs:\n"
	"\tThe program outputs alog file <prefix>.log, and two NIFTI images: <prefix>_t.nii, and\n"
	"\t<prefix>_df.nii.  <prefix>_t is of type `float' and stores the t-values computed from ANOCOVA\n"
	"\tanalyses and the specified contrast. The <prefix>_df is of type `char' and stores the voxelwise\n"
	"\tdegrees of freedom of each ANCOVA analysis.\n\n");

	printf("Related programs:\n"
	"\tthreshold_tmap ccstats\n\n");
	exit(0);
}

// Ensures that the dataTypeCode consists entirely of 'i' and 'n' characters.
// Determines nc, number of columns of the data file.
int checkDataTypeCode(const char *dataTypeCode)
{
	int nc=0;

	if(opt_v) printf("dataTypeCode = \"%s\"\n",dataTypeCode);
	fprintf(logFilePtr,"dataTypeCode = \"%s\"\n",dataTypeCode);

	nc = strlen(dataTypeCode);

	if(nc<=0) 
	{
		printf("\n\ncheckDataTypeCode(): Illegal dataTypeCode: \"%s\", aborting ...\n\n",dataTypeCode);
		exit(0);
	}

	if(dataTypeCode[0] != 'i')
	{
		printf("\n\ncheckDataTypeCode(): The first element of dataTypeCode must be 'i'");
		printf("\ncheckDataTypeCode(): Illegal dataTypeCode: \"%s\", aborting ...\n\n",dataTypeCode);
		exit(0);
	}

	for(int i=0; i<nc; i++)
	if( dataTypeCode[i] != 'i' && dataTypeCode[i] != 'n' )
	{
		printf("\n\ncheckDataTypeCode(): Illegal dataTypeCode: \"%s\", aborting ...\n\n",dataTypeCode);
		exit(0);
	}

	if(opt_v) printf("Number of variables (nc) = %d\n",nc);
	fprintf(logFilePtr,"Number of variables (nc) = %d\n",nc);

	return(nc);
}

void checkDataMaskCode(char *dataMaskCode, int nc)
{
	if(!opt_dataMask)
	{
		for(int i=0; i<nc; i++) dataMaskCode[i]='1';
		dataMaskCode[nc]='\0';
	}

	if(opt_v) printf("dataMaskCode = \"%s\"\n",dataMaskCode);
	fprintf(logFilePtr,"dataMaskCode = \"%s\"\n",dataMaskCode);

	if( (unsigned int)nc != strlen(dataMaskCode) )
	{
		printf("\n\ncheckDataMaskCode(): dataTypeCode and dataMaskCode have different dimensions, aborting ...\n\n");
		exit(0);
	}

	if(dataMaskCode[0] != '1')
	{
		printf("\n\ncheckDataMaskCode(): The first element of dataMaskCode must be '1'");
		printf("\ncheckDataMaskCode(): Illegal dataMaskCode: \"%s\", aborting ...\n\n",dataMaskCode);
		exit(0);
	}

	for(int i=0; i<nc; i++)
	if( dataMaskCode[i] != '0' && dataMaskCode[i] != '1' )
	{
		printf("\n\ncheckDataMaskCode(): Illegal dataMaskCode: \"%s\", aborting ...\n\n",dataMaskCode);
		exit(0);
	}
}

void numberOfUnmaskedVariables(const char *dataTypeCode, const char *dataMaskCode, int *pi, int *pn, int nc)
{
	*pn = *pi = 0;

	for(int i=0; i<nc; i++)
	if( dataMaskCode[i] == '1')
	{
		if(dataTypeCode[i] == 'i') (*pi)++;
		else if(dataTypeCode[i] == 'n') (*pn)++;
	}

	if(opt_v) printf("Number of unmasked image-type variables  = %d\n",*pi);
	if(opt_v) printf("Number of unmasked numerical variables  = %d\n",*pn);
	fprintf(logFilePtr,"Number of unmasked image-type variables  = %d\n",*pi);
	fprintf(logFilePtr,"Number of unmasked numerical variables  = %d\n",*pn);

	if(*pi<=0)
	{
		printf("\n\nNumber of image-type variables must be greater than 0, aborting ...\n\n");
		exit(0);
	}

	if((*pi + *pn)<=1)
	{
		printf("\n\nNumber of unmasked variables must be greater than 1, aborting ...\n\n");
		exit(0);
	}
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

	if(opt_v) printf("Number of rows of the data matrix = %d\n",n/nc);
	fprintf(logFilePtr,"Number of rows of the data matrix = %d\n",n/nc);

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

	if(opt_v)
	{
		printf("\nImage-type data matrix (idata):\n");
		for(int i=0; i<nr; i++) 
		{
			for(int j=0; j<pi; j++) 
				printf("%s\t",idata[i*pi + j]);
			printf("\n");
		}
	}

	fprintf(logFilePtr,"\nImage-type data matrix (idata):\n");
	for(int i=0; i<nr; i++) 
	{
		for(int j=0; j<pi; j++) 
			fprintf(logFilePtr,"%s\t",idata[i*pi + j]);
		fprintf(logFilePtr,"\n");
	}

	return(idata);
}

double *read_ndata(const char *dataFile, const char *dataTypeCode, const char *dataMaskCode, int nr, int nc, int pn)
{
	double *ndata=NULL;
	FILE *fp;
	char s[512];
	int count=0;

	if( (nr*pn) > 0)
	{
		// memory allocation
		ndata = (double *)calloc(nr*pn,sizeof(double));
    	if(ndata==NULL) memory_allocation_error("ndata");
	}
	else	return(NULL);

	fp = fopen(dataFile,"r");
	if(fp == NULL) file_open_error(dataFile);

	count=0;
	for(int i=0; i<nr; i++)
	for(int j=0; j<nc; j++)
	{
		if( dataTypeCode[j]=='n' && dataMaskCode[j]=='1' )
			fscanf(fp, "%lf", &ndata[count++]);
		else
			fscanf(fp, "%s", s);
	}

	fclose(fp);

	if(opt_v)
	{
		printf("\nNumerical data matrix (ndata):\n");
		for(int i=0; i<nr; i++) 
		{
			for(int j=0; j<pn; j++) 
				printf("%lf\t",ndata[i*pn + j]);
			printf("\n");
		}
	}

	fprintf(logFilePtr,"\nNumerical data matrix (ndata):\n");
	for(int i=0; i<nr; i++) 
	{
		for(int j=0; j<pn; j++) 
			fprintf(logFilePtr,"%lf\t",ndata[i*pn + j]);
		fprintf(logFilePtr,"\n");
	}

	return(ndata);
}

void checkDimension_vancova(int N, char **imagefile, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz)
{
   NIFTIIMAGE im0;;

   if(N==0) return;

   im0.readheader(imagefile[0]);

   *nx = im0.nx();
   *ny = im0.ny();
   *nz = im0.nz();
   *dx = im0.dx();
   *dy = im0.dy();
   *dz = im0.dz();

   if(im0.datatype() != 4)
   {
      printf("\n\ncheckDimension_vancova(): This program cannot handle NIFTI image data type %d, aborting ...\n",im0.datatype());
      exit(0);
   }

   if(opt_v)
   {
      printf("\nImage matrix size (nx,ny,nz) = %d x %d x %d",*nx,*ny,*nz);
      printf("\nImage voxel size (dx,dy,dz)= %f x %f x %f\n",*dx,*dy,*dz);
   }
   fprintf(logFilePtr,"\nImage matrix size (nx,ny,nz) = %d x %d x %d",*nx,*ny,*nz);
   fprintf(logFilePtr,"\nImage voxel size (dx,dy,dz)= %f x %f x %f\n",*dx,*dy,*dz);

   for(int i=1; i<N; i++)
   {
      NIFTIIMAGE im;;

      im.readheader(imagefile[i]);

      if( *nx != im.nx() ||  *ny != im.ny() ||  *nz != im.nz() ) 
      {
            printf("\n\nImage %d: %s",i+1,imagefile[i]);
            printf("\n\tMatrix size = %d x %d x %d", im.nx(), im.ny(), im.nz());
            printf("\n\nAll input images must be of size: %d x %d x %d\n\n",*nx,*ny,*nz);
            exit(0);
      }

      if(im.datatype() != 4)
      {
            printf("\n\ncheckDimension_vancova(): This program cannot handle NIFTI image data type %d, aborting ...\n",im.datatype());
            exit(0);
      }
   }
}

short *readMask(char *maskFile, int nx, int ny, int nz, int *nbv)
{
   NIFTIIMAGE im;
   short *imdata;
   short *mask;
   int nv;

   nv = nx*ny*nz;

   mask = (short *)calloc(nv,sizeof(short));
   for(int i=0; i<nv; i++) mask[i]=1;

   if(opt_m)
   {
      if(opt_v) printf("\nMask image = %s\n",maskFile);
      fprintf(logFilePtr,"\nMask image = %s\n",maskFile);

      // read the mask
      im.read(maskFile);
      imdata = (short *)im.getdata();

      if(opt_v)
      {
         printf("Mask matrix size = %d x %d x %d\n", im.nx(), im.ny(), im.nz());
         printf("Mask voxel size = %f x %f x %f\n", im.dx(), im.dy(), im.dz());
      }

      fprintf(logFilePtr,"Mask matrix size = %d x %d x %d\n",im.nx(), im.ny(), im.nz());
      fprintf(logFilePtr,"Mask voxel size = %f x %f x %f\n", im.dx(), im.dy(), im.dz());

      if( im.nv() == nv )
      {
         for(int i=0; i<nv; i++) mask[i] = imdata[i];
      }
      else
      {
         printf("readmask(): Incompatible mask dimensions, reverting to default mask ...\n");
         fprintf(logFilePtr,"readmask(): Incompatible mask dimensions, reverting to default mask ...\n");
      }
   }

   *nbv = 0;
   for(int i=0; i<nv; i++) 
   if(mask[i]!=0) 
   {
      mask[i]=1;
      (*nbv)++;
   }

   return(mask);
}

char **maskTransposeSave(char **idata, short *mask, int nr, int pi, int nbv, float FWHM)
{
   char tempFilename[512];
   char **maskedImageFile=NULL;
   short *masked_image=NULL;
   struct stat fileinfo;   // structure into which information is placed about the file

   if(pi<0) return(NULL);
   
   if(nbv<=0) return(NULL);
   
   if(FWHM<0.0) FWHM=0.0;

   // allocate memory for masked images.
   masked_image = (short *)calloc(nbv,sizeof(short));
   if(masked_image==NULL) memory_allocation_error("masked_image");

   // allocat memory for maskedImageFile;
   maskedImageFile= (char **)calloc(pi,sizeof(char *));
   if(maskedImageFile==NULL) memory_allocation_error("maskedImageFile");

   for(int j=0; j<pi; j++) 
   {
      // memory allocation
      maskedImageFile[j] = (char *)calloc(512,sizeof(char));
      if(maskedImageFile[j]==NULL) memory_allocation_error("maskedImageFile[]");

      get_temp_filename(maskedImageFile[j]);
      get_temp_filename(tempFilename);

      for(int i=0; i<nr; i++) 
         mask_and_save_nii(idata[i*pi+j], tempFilename, mask, masked_image, nbv, FWHM);

		// obtain information about the masked image series file.
		if( stat(tempFilename, &fileinfo) == -1 )
		{
			printf("\n\nstat() failure. Cannot read %s, aborting ...\n\n",tempFilename);
			remove(tempFilename);
			exit(0);
		}

		if( fileinfo.st_size != (unsigned int)(sizeof(short)*nr*nbv) )
		{
			printf("\n\nmaskTransposeSave(): File %s does not have the expected size, aborting ...\n\n",tempFilename);
			remove(tempFilename);
			exit(0);
		}

		read_transpose_save(tempFilename, maskedImageFile[j], nr, 0);

		// obtain information about the masked image series file.
		if( stat(maskedImageFile[j], &fileinfo) == -1 )
		{
			printf("\n\nstat() failure. Cannot read %s, aborting ...\n\n",maskedImageFile[j]);
			remove(maskedImageFile[j]);
			remove(tempFilename);
			exit(0);
		}

		if( fileinfo.st_size != (unsigned int)(sizeof(short)*nr*nbv) )
		{
			printf("\n\nmaskTransposeSave(): File %s does not have the expected size, aborting ...\n\n",maskedImageFile[j]);
			remove(maskedImageFile[j]);
			remove(tempFilename);
			exit(0);
		}

		remove(tempFilename);
	}

   if(masked_image!=NULL) free(masked_image);

   return(maskedImageFile);
}

void set_X_matrix(double *ndata, double *X, int nr, int pn, int p, int nc, const char *dataTypeCode, const char *dataMaskCode)
{
   int X_j=0, ndata_j=0;

   // X_j goes from 0 to p
   // ndata_j goes from 0 to pn
   X_j=ndata_j=0;
   for(int j=1; j<nc; j++)		// j starts from one since the 1st column is the dependent variable
   {
      if(dataTypeCode[j]=='i' && dataMaskCode[j]=='1') X_j++;	// skip image-type data column

      if(dataTypeCode[j]=='n' && dataMaskCode[j]=='1')		// record numerical data column
      {
         for(int i=0; i<nr; i++) X[i*p + X_j]=ndata[i*pn + ndata_j];
         X_j++;
         ndata_j++;
      }
   }
}

double *readContrast(const char *contrastFile, int nc, int p, const char *dataMaskCode, int *partial_corr_flag, int *F_flag)
{
	double *contrast=NULL;
	double *c=NULL;
	int count;
	FILE *fp;

	if(nc<=0 || p<=0) return(NULL);

	contrast = (double *)calloc(nc, sizeof(double));
    if(contrast==NULL) memory_allocation_error("readContrast(): contrast");

	c = (double *)calloc(p, sizeof(double));
    if(c==NULL) memory_allocation_error("readContrast(): c");

	fp = fopen(contrastFile,"r");
	if(fp == NULL) file_open_error(contrastFile);

	for(int i=0; i<nc; i++) fscanf(fp, "%lf", contrast+i); 

	fclose(fp);

	count = 0;
	for(int i=1; i<nc; i++)
	if(dataMaskCode[i]=='1') c[count++]=contrast[i];

	if(opt_v) 
	{
		printf("contrast = ");
		for(int i=0; i<p; i++) printf("%5.2lf\t",c[i]);
		printf("\n");
	}
	fprintf(logFilePtr,"contrast = ");
	for(int i=0; i<p; i++) fprintf(logFilePtr,"%5.2lf\t",c[i]);
	fprintf(logFilePtr,"\n");

	if(contrast!=NULL) free(contrast);

	// determines if the contrast indicates a partial correlation analysis
	*partial_corr_flag=0;
	for(int i=0; i<p; i++) if(c[i]!=0.0) (*partial_corr_flag)++;
	if(*partial_corr_flag==1)
	{
		for(int i=0; i<p; i++) 
		{ 
			if(c[i]!=0.0 && c[i]!=1.0) 
				*partial_corr_flag=0;
		}
	}
	else
	{
		*partial_corr_flag=0;
	}

	// determines if the contrast indicates an F-test 
	*F_flag=0;
	for(int i=0; i<p; i++) if(c[i]!=0.0) (*F_flag)++;
	if(*F_flag>1)
	{
		for(int i=0; i<p; i++) 
		{ 
			if(c[i]!=0.0 && c[i]!=1.0) 
				*F_flag=0;
		}
	}
	else
	{
		*F_flag=0;
	}

	return(c);
}

double compute_t_value(double *y, double *X, double *c, int n, int p, double *df, float *r2)
{
	int rank;
	float *G;
	double *beta;
	double *Xty;
	double ctbeta;
	double ctGc, var, t;
	double yty, ytPy;
	double mu;

	mu=removeVectorMean(y, n);
	removeVectorMean(X, n, p);
	scaleAbsToOne(X, n, p);

	G = (float *)calloc(p*p, sizeof(float));
	Xty = (double *)calloc(p,sizeof(double));
	beta = (double *)calloc(p,sizeof(double));

	rank = ginverse(X, n, p, G);

	mat_trans_mat(X, n, p, y, 1, Xty);
	multi(G,p,p,Xty,p,1,beta);

	ctGc = xtAx(G, c, p);		// computes ct * G * c
	ctbeta = dot(beta,c,p);		// computes ct * beta
	ytPy = dot(Xty,beta,p);		// computes yt * X * G * Xt * y
	yty = dot(y,y,n);			// computes yt * y

	*df = n-1.0-rank;

	if( (*df) > 0.0 )
		var = (yty-ytPy)/(*df);
	else
		var = 0.0;

	if(var*ctGc > 0.0)
		t = ctbeta / sqrt(var*ctGc);
	else
		t = 0.0;

	{
		double E1;
		double E;

		if(ctGc!=0.0)
			E1=ctbeta*ctbeta/ctGc;
		else
			E1=0.0;

		E=yty-ytPy+E1;

		if(E!=0.0)
			*r2=(float)(E1/E);
		else
			*r2=0.0;
	}

/***
{
	double *tt;
	double sum=0.0;
	tt = (double *)calloc(n, sizeof(double));

	for(int i=0; i<n; i++) tt[i] = mu + beta[0]*X[i*p + 0] + beta[1]*X[i*p + 1];

	sum=0.0;
	for(int i=0; i<28; i++) sum += tt[i];
	printf("\n%lf\n",sum/28.0);

	sum=0.0;
	for(int i=28; i<n; i++) sum += tt[i];
	printf("\n%lf\n",sum/36.0);
}
**/

	free(G);
	free(Xty);
	free(beta);

	return(t);
}

// Model: y = X*beta + e = X1*beta1 + X2*beta2 + e
double compute_F_value(double *y, double *X, double *c, int n, int p, double *dfn, double *dfd, float *r2)
{
	int rank, rank1, rank2;
	float *G=NULL, *G1=NULL, *G2=NULL;
	double *X1=NULL, *X2=NULL;
	double *beta=NULL, *beta1=NULL, *beta2=NULL;
	double *Xty;
	double evar, var1, F;
	double yty, ytPy;
	double mu;
	double *P2X1=NULL, *X2tX1=NULL, *G2X2tX1=NULL;
	double *X1ty=NULL, *X2ty=NULL;
	

	int p1; 	// number of columns of matrix X1
	int p2; 	// number of columns of matrix X2
	int X1_column, X2_column;	// column indices

	
	//////////////////////////////////////////////////////////////////
	mu=removeVectorMean(y, n);
	removeVectorMean(X, n, p);
	scaleAbsToOne(X, n, p);
	yty = dot(y,y,n);			// computes yt * y

	G = (float *)calloc(p*p, sizeof(float));
	Xty = (double *)calloc(p,sizeof(double));
	beta = (double *)calloc(p,sizeof(double));

	rank = ginverse(X, n, p, G);
	mat_trans_mat(X, n, p, y, 1, Xty);
	multi(G,p,p,Xty,p,1,beta);

	ytPy = dot(Xty,beta,p);		// computes yt * X * G * Xt * y

	evar = (yty-ytPy);
	*dfd = n-1.0-rank;
	//////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////
	// This part of the code separates the matrix X into X1 and X2.
	// X1 is the subspace of interest.
	// X2 is the covariate subspace.
	//////////////////////////////////////////////////////////////////

	// determine p1
	p1=0;
	for(int i=0; i<p; i++) if(c[i]==1.0) p1++;
	if(p1==0)
	{	
		*dfd=*dfn=0.0;
		return(0.0);
	}

	X1 = (double *)calloc(n*p1, sizeof(double));

	X1_column=0;
	for(int j=0; j<p; j++) 
	{
		if(c[j]==1.0) 
		{
			// copies the jth column of X into the X1_column of X1
			for(int i=0; i<n; i++) X1[i*p1 + X1_column] = X[i*p + j];

			X1_column++;

		}
	}

	p2 = p - p1;

	if(p2!=0)
	{
		X2 = (double *)calloc(n*p2, sizeof(double));

		X2_column=0;
		for(int j=0; j<p; j++) 
		{
			if(c[j]!=1.0) 
			{
				// copies the jth column of X into the X2_column of X2
				for(int i=0; i<n; i++) X2[i*p2 + X2_column] = X[i*p + j];

				X2_column++;
			}
		}
	}
	//////////////////////////////////////////////////////////////////

	if(p2 != 0)
	{
		G2 = (float *)calloc(p2*p2, sizeof(float));
		rank2 = ginverse(X2, n, p2, G2);

		X2ty= (double *)calloc(p2, sizeof(double));
		mat_trans_mat(X2, n, p2, y, 1, X2ty);

		beta2 = (double *)calloc(p2,sizeof(double));
		multi(G2,p2,p2,X2ty,p2,1,beta2);

		X2tX1 = (double *)calloc(p2*p1, sizeof(double));
		mat_trans_mat(X2, n, p2, X1, p1, X2tX1);
	
		G2X2tX1 = (double *)calloc(p2*p1, sizeof(double));
		multi(G2,p2,p2,X2tX1,p2,p1,G2X2tX1);

		P2X1= (double *)calloc(n*p1, sizeof(double));
		multi(X2,n,p2,G2X2tX1,p2,p1,P2X1);

		for(int i=0; i<n*p1; i++) X1[i] -= P2X1[i];
	}

	G1 = (float *)calloc(p1*p1, sizeof(float));
	rank1 = ginverse(X1, n, p1, G1);

	X1ty= (double *)calloc(p1, sizeof(double));
	mat_trans_mat(X1, n, p1, y, 1, X1ty);

	beta1 = (double *)calloc(p1,sizeof(double));
	multi(G1,p1,p1,X1ty,p1,1,beta1);

	var1 = dot(beta1,X1ty,p1);
	*dfn = rank1;

	if( (evar*rank1) > 0.0)
		F = var1*(*dfd) / (evar*rank1);
	else
		F = 0.0;

	if((evar+var1)!=0.0)
		*r2=(float)(var1/(evar+var1));
	else
		*r2=0.0;

	free(G);
	free(Xty);
	free(beta);

	free(G1);
	free(X1);
	free(X1ty);
	free(beta1);

	if(p2!=0) 
	{
		free(G2);
		free(X2);
		free(X2ty);
		free(beta2);

		free(X2tX1);
		free(G2X2tX1);
		free(P2X1);
	}

	return(F);
}

int main(int argc, char **argv)
{
   float dum;
   float *motion_params;
   float *rw;
   char outputfile[1024]="";
   char imagefile[1024]=""; // 4D NIFTI image
   char paramfile[1024]=""; 
   char rwfile[1024]=""; 
   char prefix[512];
   nifti_1_header hdr;
   int nx,ny,nz,nv;
   float dx,dy,dz;
   SHORTIM im4D;
   char dataFile[512]="";
   FILE *fp;

   int nbv;
   NIFTIIMAGE im0;
   int partial_corr_flag=0;
   int F_flag=0;

   FILE **maskedImageFilePtr;

   float FWHM=0.0; // FWHM of the smoothing filter specified using the -FWHM <FWHM>. Default is set to 0 here.

   unsigned char *dfmap, *df1map, *df2map;
	char s[512];
	char **maskedImageFile=NULL;
   char **idata=NULL;
	char contrastFile[512];
	char maskFile[512];
	char logFile[512];
	char dataMaskCode[128];
	char dataTypeCode[128];		// a string composed of characters 'i' and 'n' only, where 'i' denotes
								// image type variables and 'n' denotes numerical type variables
   double *ndata=NULL;
	double *X, *y;
	double *c;
	double t, df, dfn, dfd, F;
   float *spm_t, *spm_r2, *spm_F;

	int nc=0;	// number of columns of the data file
	int nr=0; 	// number of rows of the data file
	int pi=0; 	// number of unmasked image-type independent variables + 1. 
				// The additinoal 1 is for the dependent variable, which is always image-type
	int pn=0; 	// number of unmasked number-type independent variables 
	int p=0; 	// number of unmasked independent variables

	int X_j=0, idata_j=0;
	int xi=0, yi=0, zi=0, vi=0;
	
	short *mask=NULL;
	short *sdum=NULL;

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) {
         case 'r':
            sprintf(rwfile,"%s",optarg);
            break;
         case 'i':
            sprintf(imagefile,"%s",optarg);
            break;
         case 'd':
            sprintf(dataFile,"%s",optarg);
            break;
         case 'p':
            sprintf(paramfile,"%s",optarg);
            break;
			case 'x':
				xi = atoi(optarg);
				opt_x=YES;
				break;
			case 'y':
				yi = atoi(optarg);
				opt_y=YES;
				break;
			case 'z':
				zi = atoi(optarg);
				opt_z=YES;
				break;
			case 'F':
				FWHM = atof(optarg);
				break;
			case 'n':
				opt_nozero=YES;
				break;
			case 'v':
				opt_v=YES;
				break;
			case 'c':
				sprintf(contrastFile,"%s",optarg);
				opt_c=YES;
				break;
			case 'm':
				sprintf(maskFile,"%s",optarg);
				opt_m=YES;
				break;
			case 'T':
				sprintf(dataTypeCode,"%s",optarg);
				opt_dataType=YES;
				break;
			case 'M':
				sprintf(dataMaskCode,"%s",optarg);
				opt_dataMask=YES;
				break;
			case 'o':
				sprintf(prefix,"%s",optarg);
				opt_o=YES;
				break;
			case '?':
				print_help_and_exit();
		}
	}

   ///////////////////////////////////////////////////////////////////////////////////
   // read and split the 4D image
   
   if( imagefile[0] == '\0' )
   {
      printf("Please specify a 4D NIFTI image file using -i <filename> option.\n");
      exit(0);
   }

   if(opt_v)
   {
      printf("4D NIFTI image file: %s\n",imagefile);
   }

   if( niftiFilename(prefix, imagefile) == 0 ) exit(0);
   if(opt_v)
   {
      printf("Image prefix: %s\n",prefix);
   }

   im4D.v = (short *)read_nifti_image(imagefile, &hdr);

   if(hdr.datatype != 4)
   {
      printf("Sorry, this program cannot handle NIFTI image data type %d, aborting ...\n",hdr.datatype);
      exit(0);
   }

   set_dim(im4D, hdr);

   if(opt_v)
   {
      printf("Matrix: %d x %d x %d x %d\n",im4D.nx, im4D.ny, im4D.nz, im4D.nt);
      printf("Voxel: %f x %f x %f\n",im4D.dx, im4D.dy, im4D.dz);
   }

   
   hdr.dim[0]=3;
   hdr.dim[4]=1;
   for(int i=0; i<im4D.nt; i++)
   {
      sprintf(outputfile,"%s%04d.nii",prefix,i);
      save_nifti_image(outputfile, im4D.v + i*im4D.nv, &hdr);
   }
   ///////////////////////////////////////////////////////////////////////////////////
exit(0);

   ///////////////////////////////////////////////////////////////////////////////////
   
   if( paramfile[0] == '\0' )
   {
      printf("Please specify the motion parameters file using -p <filename> option.\n");
      exit(0);
   }

   if(opt_v)
   {
      printf("Motion parameters file: %s\n",paramfile);
   }

   motion_params = (float *)calloc(im4D.nt*6,sizeof(float));
   fp = fopen(paramfile,"r");
   for(int i=0; i<im4D.nt; i++)
   {
      fscanf(fp,"%f",&dum);

      for(int j=0; j<6; j++)
      {
         fscanf(fp,"%f",&motion_params[i*6 + j]);
      }

      fscanf(fp,"%f",&dum);
      fscanf(fp,"%f",&dum);
   }
   fclose(fp);

   ///////////////////////////////////////////////////////////////////////////////////
   
   ///////////////////////////////////////////////////////////////////////////////////
   
   if( rwfile[0] == '\0' )
   {
      printf("Please specify the reference waveform file using -rw <filename> option.\n");
      exit(0);
   }

   if(opt_v)
   {
      printf("Reference waveform file: %s\n",rwfile);
   }

   rw = (float *)calloc(im4D.nt*4,sizeof(float));
   fp = fopen(rwfile,"r");
   for(int i=0; i<im4D.nt; i++)
   {
      for(int j=0; j<4; j++)
      {
         fscanf(fp,"%f",&rw[i*4 + j]);
      }
   }
   fclose(fp);

   ///////////////////////////////////////////////////////////////////////////////////
   
   if( dataFile[0] == '\0' )
   {
      printf("Please specify an output data file using -d <filename> option.\n");
      exit(0);
   }

   if(opt_v)
   {
      printf("Output data file: %s\n",dataFile);
   }

   fp = fopen(dataFile,"w");
   for(int i=0; i<im4D.nt; i++)
   {
      sprintf(outputfile,"%s%04d.nii",prefix,i);
      fprintf(fp,"%s\t",outputfile);

      for(int j=0; j<6; j++)
      {
         fprintf(fp,"%7.4f\t",motion_params[i*6 + j]);
      }

      for(int j=0; j<4; j++)
      {
         fprintf(fp,"%8.6f\t",rw[i*4 + j]);
      }

      fprintf(fp,"\n");
   }
   fclose(fp);
   exit(0);
   

   if(!opt_d || !opt_dataType || !opt_c)
      print_help_and_exit();

   // if no prefix is specified at the command line, set the default prefix to "V"
   if(!opt_o) sprintf(prefix,"%s","V");

   // name the log file <prefix>.log
   sprintf(logFile,"%s.log",prefix);

   // open the log file for writting
   logFilePtr = fopen(logFile,"w");
   if(logFilePtr == NULL) file_open_error(logFile);

   // log information
   if(opt_v) printf("Output file prefix = %s\n",prefix);
   fprintf(logFilePtr,"Output file prefix = %s\n",prefix);

   // log information
   if(opt_v)   printf("Smoothing filter FWHM = %7.3f mm\n",FWHM);
   fprintf(logFilePtr,"Smoothing filter FWHM = %7.3f mm\n",FWHM);

   // ensure dataTypeCode consists entirely of i's and n's 
   nc = checkDataTypeCode(dataTypeCode);

   // ensure the size of dataMaskCode matches the size of dataTypeCode
   // ensure dataMaskCode consists entirely of 1's and 0's 
   // if no dataMaskCode is specified, set the default to "11...1" (all 1's)
   checkDataMaskCode(dataMaskCode,nc);

   // pi is the number of image-type variables (both the dep. variable and the indep. variables)
   // pn is the number of number-type variables both the dep. variable and the indep. variables)
   numberOfUnmaskedVariables(dataTypeCode, dataMaskCode, &pi, &pn, nc);
   p = pi+pn-1;	// p is the total number of unmasked independent variables

   // must be cased after nc and p are known
   // c is a double array of size p
   c = readContrast(contrastFile, nc, p, dataMaskCode, &partial_corr_flag, &F_flag);

   // nr is the number of rows of the dataFile
   // nc MUST equal to the number of columns of the dataFile, this will not be checked by the program
   nr = numberOfRows(dataFile,nc);

   if(pi>0) idata = read_idata(dataFile, dataTypeCode, dataMaskCode, nr, nc, pi);

   if(pn>0) ndata = read_ndata(dataFile, dataTypeCode, dataMaskCode, nr, nc, pn);

   // ensure all input image data have the same matrix size and data type 4
   checkDimension_vancova(nr*pi, idata, &nx, &ny, &nz, &dx, &dy, &dz);
   nv = nx*ny*nz;

   // read mask image
   // if no mask image is specified or the mask image matrix size does not match the
   // input images' matrix size, create a default mask where all voxels are 1.
   mask = readMask(maskFile, nx, ny, nz, &nbv);

   // if nozero option is selected, adjust mask to remove voxels where the dependent
   // variable assumes a value of 0 at any of it's nc elements.
   if(opt_nozero)
   for(int i=0; i<nr; i++)
   {
      NIFTIIMAGE im;
      short *imdata;

      im.read(idata[pi*i]);
      imdata = (short *)im.getdata();

      for(int v=0; v<nv; v++) 
      if(mask[v]!=0 && imdata[v]==0)
      {
         mask[v]=0; 
         nbv--;
      }
   }

   if(opt_v) printf("Number of mask voxels = %d\n",nbv);
   fprintf(logFilePtr,"Number of mask voxels = %d\n",nbv);

   if(nbv<=0) 
   {
      printf("Error: number of mask voxels must be positive, aborting ...\n");
      exit(0);
   }

   if(nr>0) 
   {	
      y = (double *)calloc(nr, sizeof(double));
      if(y==NULL) memory_allocation_error("y");

      sdum = (short *)calloc(nr, sizeof(short));
      if(sdum==NULL) memory_allocation_error("sdum");
   }

   // X is (nr by p) matrix (stores all unmasked independent variables)
   if(p>0 && nr>0) X = (double *)calloc(nr*p, sizeof(double));
   if(X==NULL) memory_allocation_error("X");

   set_X_matrix(ndata, X, nr, pn, p, nc, dataTypeCode, dataMaskCode);

   if(pi>0) maskedImageFile=maskTransposeSave(idata, mask, nr, pi, nbv, FWHM);

   if(pi>0) maskedImageFilePtr = (FILE **)calloc(pi,sizeof(FILE *));

   // open all maskdImageFile's
   for(int i=0; i<pi; i++)
   {
      maskedImageFilePtr[i]=fopen(maskedImageFile[i],"r");
      if(maskedImageFile[i] == NULL) file_open_error(maskedImageFile[i]);
   }

   // vi=voxel index
   if(opt_x && opt_y && opt_z) 
   {
      vi = zi*(nx*ny) + yi*nx + xi;
      if(vi<0 || vi>=nv)
      {
         printf("Warning: specified voxel (%d,%d,%d) is out of image matrix range.\n",xi,yi,zi);
         vi = -1;
      }
   } else vi = -1;

   spm_r2 = (float *)calloc(nv, sizeof(float));
   if(spm_r2==NULL) memory_allocation_error("spm_r2");

   if(F_flag==0)
   {
      spm_t = (float *)calloc(nv, sizeof(float));
      if(spm_t==NULL) memory_allocation_error("spm_t");

      dfmap = (unsigned char *)calloc(nv, sizeof(unsigned char));
      if(dfmap==NULL) memory_allocation_error("dfmap");
   }
   else
   {
      spm_F = (float *)calloc(nv, sizeof(float));
      if(spm_F==NULL) memory_allocation_error("spm_F");

      df1map = (unsigned char *)calloc(nv, sizeof(unsigned char));
      if(df1map==NULL) memory_allocation_error("df1map");

      df2map = (unsigned char *)calloc(nv, sizeof(unsigned char));
      if(df2map==NULL) memory_allocation_error("df2map");
   }

	for(int vox=0; vox<nv; vox++)
	if(mask[vox]!=0)
	{
		// since matrix X is mean-corrected in subsequent processing
		// and we don't want to save the mean corrected matrix at voxel 'vi'
		if(vi == vox )  set_X_matrix(ndata, X, nr, pn, p, nc, dataTypeCode, dataMaskCode);

		fread(sdum, sizeof(short), nr, maskedImageFilePtr[0]);
		for(int i=0; i<nr; i++) y[i]=sdum[i];

		X_j=idata_j=0;
		for(int j=1; j<nc; j++)		// j starts from one since the 1st column is the dependent variable
		{
			if(dataTypeCode[j]=='n' && dataMaskCode[j]=='1') X_j++;	// skip numerical data column

			if(dataTypeCode[j]=='i' && dataMaskCode[j]=='1')		// record image-type data column
			{
				fread(sdum, sizeof(short), nr, maskedImageFilePtr[idata_j + 1]);
				for(int i=0; i<nr; i++) X[i*p + X_j]=sdum[i];
				X_j++;
				idata_j++;
			}
		}

		if(vi == vox )  
		{
			sprintf(s,"y vector at voxel (%d,%d,%d):",xi,yi,zi);
			printMatrix(y,nr,1,s, logFilePtr);
			if(opt_v) printMatrix(y,nr,1,s, NULL);

			sprintf(s,"X matrix at voxel (%d,%d,%d):",xi,yi,zi);
			printMatrix(X,nr,p,s, logFilePtr);
			if(opt_v) printMatrix(X,nr,p,s, NULL);
		}

		if(F_flag==0)
		{
			spm_t[vox] = (float) compute_t_value(y, X, c, nr, p, &df, &spm_r2[vox]);

			dfmap[vox] = (unsigned char)df;

			if(vi == vox )  
			{
				t=spm_t[vi];
				if(opt_v) 
				{
					printf("\nt-value at (%d,%d,%d) = %f with %d d.f.\n",xi,yi,zi,spm_t[vi],dfmap[vi]);
					if(partial_corr_flag) printf("Coefficient of Determination at (%d,%d,%d) = %f\n",xi,yi,zi,spm_r2[vi]);
				}
				fprintf(logFilePtr,"\nt-value at (%d,%d,%d) = %f with %d d.f.\n",xi,yi,zi,spm_t[vi],dfmap[vi]);
				if(partial_corr_flag) fprintf(logFilePtr,"Coefficient of Determination at (%d,%d,%d) = %f\n",xi,yi,zi,spm_r2[vi]);
			}
		}
		else
		{
			spm_F[vox] = (float) compute_F_value(y, X, c, nr, p, &dfn, &dfd, &spm_r2[vox]);

			df1map[vox] = (unsigned char)dfn;
			df2map[vox] = (unsigned char)dfd;

			if(vi == vox )  
			{
				F=spm_F[vi];

				if(opt_v) 
				{
					printf("\nF-value at (%d,%d,%d) = %f with (%d,%d) d.f.\n",xi,yi,zi,spm_F[vi],df1map[vi],df2map[vi]);
					printf("Coefficient of Determination at (%d,%d,%d) = %f\n",xi,yi,zi,spm_r2[vi]);
				}
				fprintf(logFilePtr,"\nF-value at (%d,%d,%d) = %f with (%d,%d) d.f.\n",xi,yi,zi,spm_F[vi],df1map[vi],df2map[vi]);
				fprintf(logFilePtr,"Coefficient of Determination at (%d,%d,%d) = %f\n",xi,yi,zi,spm_r2[vi]);
			}
		}
	}

   // get an example header
   im0.readheader(idata[0]);
   hdr = im0.getheader();

   hdr.scl_slope = 1.0;
   hdr.scl_inter = 0.0;

   if(partial_corr_flag || F_flag)
   {
      hdr.datatype=NIFTI_TYPE_INT32;
      hdr.bitpix = 32;
      sprintf(s,"%s_r2.nii",prefix);
      save_nifti_image(s, spm_r2, &hdr);
   }

   if(F_flag==0)
   {
      hdr.datatype=NIFTI_TYPE_INT32;
      hdr.bitpix = 32;
      sprintf(s,"%s_t.nii",prefix);
      save_nifti_image(s, spm_t, &hdr);

      hdr.datatype=NIFTI_TYPE_UINT8;
      hdr.bitpix = 8;
      sprintf(s,"%s_df.nii",prefix);
      save_nifti_image(s, dfmap, &hdr);
   }
   else
   {
      hdr.datatype=NIFTI_TYPE_INT32;
      hdr.bitpix = 32;
      sprintf(s,"%s_F.nii",prefix);
      save_nifti_image(s, spm_F, &hdr);

      hdr.datatype=NIFTI_TYPE_UINT8;
      hdr.bitpix = 8;
      sprintf(s,"%s_df1.nii",prefix);
      save_nifti_image(s, df1map, &hdr);
   
      sprintf(s,"%s_df2.nii",prefix);
      save_nifti_image(s, df2map, &hdr);
   }

	/////////////////////////////////////////////////////////////////////////////////////
	if(F_flag==0)
	{
		FILE *fp;
		int nbv=0, N;
		float mint, maxt;
		float tval;
		double sum = 0.0;;
		int count;
		double df=1;

		count = 0;
		for(int i=0; i<nv; i++)
		if(spm_t[i]!=0.0)
		{
			sum += dfmap[i];
			count++;
		}

		if(count!=0) df = sum/count;
		df = nearbyint(df);

		for(int i=0; i<nv; i++)
		{
			if(mask[i]!=0) nbv++;
		}
		//printf("\nnbv=%d\n",nbv);

		minmax(spm_t, nv, mint, maxt);
		//printf("mint=%f, maxt=%f\n",mint,maxt);

		sprintf(s,"%s_FDR.tbl",prefix);
		fp=fopen(s,"w");

		fprintf(fp,"df=%lf\n\n",df);
		fprintf(fp,"t-value threshold\tp-value threshold\tFDR\n");

		for(tval=0.0; tval<=maxt; tval+=0.1)
		{
			N=0;
			for(int i=0; i<nv; i++)
			if(mask[i]!=0 && (spm_t[i]>tval || -spm_t[i]>tval)) N++;
		}

		fclose(fp);
	}
    /////////////////////////////////////////////////////////////////////////////////////


	// close files
	for(int i=0; i<pi; i++)
	if(maskedImageFilePtr[i]!=NULL) fclose(maskedImageFilePtr[i]);
	if(maskedImageFilePtr!=NULL) free(maskedImageFilePtr);

	// remove these files before freeing the memory of maskedImageFile
	for(int i=0; i<pi; i++) remove(maskedImageFile[i]);

	// free memory
	if(maskedImageFile!=NULL)
	{
		for(int i=0; i<pi; i++) 
			if(maskedImageFile[i]!=NULL) free(maskedImageFile[i]);
    	free(maskedImageFile);
	}

	// free memory
	if(idata!=NULL)
	{
		for(int i=0; i<nr*pi; i++) 
			if(idata[i]!=NULL) free(idata[i]);
    	free(idata);
	}

	// free memory
	if(ndata!=NULL) free(ndata);
	if(mask!=NULL) free(mask);
	if(y!=NULL) free(y);
	if(X!=NULL) free(X);
	if(sdum!=NULL) free(sdum);
	if(c!=NULL) free(c);

	if( spm_r2 != NULL ) free(spm_r2);

	if(F_flag==0)
	{
		if( spm_t  != NULL ) free(spm_t);
		if( dfmap  != NULL ) free(dfmap);
	}
	else
	{
		if( spm_F  != NULL ) free(spm_F);
		if( df1map  != NULL ) free(df1map);
		if( df2map  != NULL ) free(df2map);
	}

	// close log file
	fclose(logFilePtr);
}
