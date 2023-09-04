#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <babak_lib.h>

#define NO 0
#define YES 1
#define LIN 1
#define NEARN 2

//static float v1,v2,v3,v4;
//static float w1,w2,w3,w4;

//////////////////////////////////////////////////////////////////////////////////////////////////
int opt;

static struct option options[] =
{
  {"-h",  0, 'h'},
  {"-help",  0, 'h'},
  {"-i",  1, 'i'},
  {"-t",  1, 't'},
  {"-o",  1, 'o'},
  {"-v",  0, 'v'},
  {"-T", 1, 'T'},
  {"-nx", 1, '1'},
  {"-ny", 1, '2'},
  {"-nz", 1, '3'},
  {"-dx", 1, '4'},
  {"-dy", 1, '5'},
  {"-dz", 1, '6'},
  {"-cubicspline", 0, '7'},
  {"-nn", 0, '8'},
  {0, 0, 0}
};

int opt_T=NO;
int opt_nx=NO;
int opt_ny=NO;
int opt_nz=NO;
int opt_dx=NO;
int opt_dy=NO;
int opt_dz=NO;
int opt_o=NO;
int opt_cubicspline=NO;
int opt_nn=NO;
//////////////////////////////////////////////////////////////////////////////////////////////////

void print_help_and_exit()
{
  printf("\nUsage: reslice [-v] [-cubicspline -nn] [-nx <nx> -ny <ny> -nz <nz>]\n"
  "[-dx <dx> -dy <dy> -dz <dz>] [-t <targetImageFile>]\n" 
  "-i <inputImageFile> -o <outputImageFile> -T <transformation matrix>\n\n"
  "Required flag:\n"
  "-i <inputImageFile> the input image to be resliced\n\n"
  "-o <outputImageFile> name for the output resliced image\n\n"
  "-T <transformation matrix> Affine transformation used for reslicing\n\n"
  
  "Optional flag:\n"
  "-cubicspline applies the cubic spline interpolation method\n\n"
  "-nn applies the nearest neighbor interpolation method\n\n"
  "-nx <nx> -ny <ny> -nz <nz> specifies the output matrix size\n\n"
  "-dx <dx> -dy <dy> -dz <dz> specifies the output voxel size\n\n" 
  "-t <targetImageFile> used to specify the output matrix and voxel sizes\n\n"); 
  exit(0);
}

int main(int argc, char **argv)
{
  nifti_1_header iphdr;
  nifti_1_header ophdr;
  char inputImageFile[1024]="";
  char outputImageFile[1024]="";
  char targetImageFile[1024]="";
  char transFile[1024]="";
  int nx=0,ny=0,nz=0;
  int nx2=0,ny2=0,nz2=0;
  float dx=0.0,dy=0.0,dz=0.0;
  float dx2=0.0,dy2=0.0,dz2=0.0;
  float T[16],*invT;
  short *im_in,*im_out;

  while ((opt = getoption(argc, argv, options)) != -1 )
  {
		switch (opt) {
			case 'v':
				opt_v=YES;
				break;
			case 'T':
				sprintf(transFile,"%s",optarg);
				opt_T=YES;
				break;
			case '1':
				nx2=atoi(optarg);
				opt_nx=YES;
				break;
			case '2':
				ny2=atoi(optarg);
				opt_ny=YES;
				break;
			case '3':
				nz2=atoi(optarg);
				opt_nz=YES;
				break;
			case '4':
				dx2=atof(optarg);
				opt_dx=YES;
				break;
			case '5':
				dy2=atof(optarg);
				opt_dy=YES;
				break;
			case '6':
				dz2=atof(optarg);
				opt_dz=YES;
				break;
			case 'i':
				sprintf(inputImageFile,"%s",optarg);
				break;
			case 't':
				sprintf(targetImageFile,"%s",optarg);
				break;
			case 'o':
				sprintf(outputImageFile,"%s",optarg);
				opt_o=YES;
				break;
			case '7':
				opt_cubicspline=YES;
				break;
			case '8':
				opt_nn=YES;
				break;
			case 'h':
				print_help_and_exit();
		}
  }

//////////////////////////////////////////////////////////////////////////////////////////////////
//receive input image
//////////////////////////////////////////////////////////////////////////////////////////////////
  if(inputImageFile[0]=='\0')
  {
    printf("Please specify the image to be resliced using the \"-i <inputImageFile>\" flag.\n");
    exit(0);
  }

  if(opt_v)
  {
    printf("Input image = %s\n",inputImageFile);
  }

  im_in = (short *)read_nifti_image(inputImageFile, &iphdr);
  if(im_in == NULL)
  {
    printf("Error reading %s, aborting ...\n", inputImageFile);
    exit(0);
  }

  nx=iphdr.dim[1]; ny=iphdr.dim[2]; nz=iphdr.dim[3];
  dx=iphdr.pixdim[1]; dy=iphdr.pixdim[2]; dz=iphdr.pixdim[3]; 

  // to deal with the sometimes -tive voxel dimensions in SPM/FSL data sets
  if(dx<0.0) dx *= -1.0; if(dy<0.0) dy *= -1.0; if(dz<0.0) dz *= -1.0;

  if(opt_v)
  {
    printf("Matrix size = %d x %d x %d\n",nx,ny,nz);
    printf("Voxel size = %f x %f x %f\n",dx,dy,dz);
  }
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
//Receive output image information
//////////////////////////////////////////////////////////////////////////////////////////////////
  if(outputImageFile[0]=='\0')
  {
    printf("Please specify a name for the resliced image using \"-o <outputImageFile>\" flag.\n");
    exit(0);
  }

  if(opt_v)
  {
    printf("Output image = %s\n",outputImageFile);
  }

  if(!opt_nx || nx2<=0) nx2=nx;
  if(!opt_ny || ny2<=0) ny2=ny;
  if(!opt_nz || nz2<=0) nz2=nz;
  if(!opt_dx || dx2<=0) dx2=dx;
  if(!opt_dy || dy2<=0) dy2=dy;
  if(!opt_dz || dz2<=0) dz2=dz;

  if(targetImageFile[0]!='\0')
  {
    ophdr = read_NIFTI_hdr(targetImageFile);

    nx2=ophdr.dim[1]; ny2=ophdr.dim[2]; nz2=ophdr.dim[3];
    dx2=ophdr.pixdim[1]; dy2=ophdr.pixdim[2]; dz2=ophdr.pixdim[3]; 

    // to deal with the sometimes -tive voxel dimensions in SPM/FSL data sets
    if(dx2<0.0) dx2 *= -1.0; if(dy2<0.0) dy2 *= -1.0; if(dz2<0.0) dz2 *= -1.0;
  }

  if(opt_v)
  {
    printf("Matrix size = %d x %d x %d\n",nx2,ny2,nz2);
    printf("Voxel size = %f x %f x %f\n",dx2,dy2,dz2);
  }
//////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
//Receive transformation matrix
//////////////////////////////////////////////////////////////////////////////////////////////////
  if(transFile[0]=='\0')
  {
    printf("Please specify the transformation matrix using the \"-T <transformation matrix>\" flag.\n");
    exit(0);
  }

  if(opt_v)
  {
    printf("Transformation matrix = %s\n",transFile);
  }

  // ensure the input transformatin file exists, is readable, and has the expected size
  if ( checkFileExistence(transFile)==0 )
  {
    printf("Error: File %s does not exist! Aborting ...\n",transFile);
    exit(0);
  }

  if ( checkFileReadOK(transFile)==0 )
  {
    printf("Error: Read permission for %s denied! Aborting ...\n",transFile);
    exit(0);
  }

  loadTransformation( transFile, T);

  if(opt_v)
  {
    printf("%9.6f  %9.6f  %9.6f  %9.6f\n",T[0],T[1],T[2],T[3]);
    printf("%9.6f  %9.6f  %9.6f  %9.6f\n",T[4],T[5],T[6],T[7]);
    printf("%9.6f  %9.6f  %9.6f  %9.6f\n",T[8],T[9],T[10],T[11]);
    printf("%9.6f  %9.6f  %9.6f  %9.6f\n",T[12],T[13],T[14],T[15]);
  }

  invT=inv4(T);
//////////////////////////////////////////////////////////////////////////////////////////////////

  if(opt_v)
  {
    if(opt_cubicspline)
      printf("Applying the cubic spline interpolation method ...\n");
    else if(opt_nn)
      printf("Applying the nearest neighbor interpolation method ...\n");
    else
      printf("Applying the trilinear interpolation method ...\n");
  }

  if(opt_cubicspline)
    im_out=resliceImageCubicSpline(im_in, nx, ny, nz, dx, dy, dz, nx2, ny2, nz2, dx2, dy2, dz2, invT);
  else if(opt_nn)
    im_out=resliceImage(im_in, nx, ny, nz, dx, dy, dz, nx2, ny2, nz2, dx2, dy2, dz2, invT, NEARN);
  else
  {
    im_out=resliceImage(im_in, nx, ny, nz, dx, dy, dz, nx2, ny2, nz2, dx2, dy2, dz2, invT, LIN);
  }

  if(targetImageFile[0]=='\0')
  {
    ophdr = read_NIFTI_hdr(inputImageFile);
    ophdr.dim[0]=3; ophdr.dim[1]=nx2; ophdr.dim[2]=ny2; ophdr.dim[3]=nz2;
    ophdr.pixdim[1]=dx2; ophdr.pixdim[2]=dy2; ophdr.pixdim[3]=dz2; 
  }

  save_nifti_image(outputImageFile, im_out, &ophdr);
}
