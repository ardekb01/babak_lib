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
        {"-i",  1, 'i'},
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
int opt_i=NO;
int opt_o=NO;
int opt_cubicspline=NO;
int opt_nn=NO;
//////////////////////////////////////////////////////////////////////////////////////////////////

void print_help_and_exit()
{
	printf("Usage: reslice [-v] [-cubicspline -nn] -T <transformation> [-nx <nx> -ny <ny> -nz <nz>]\n"
	"[-dx <dx> -dy <dy> -dz <dz>] -i <inputImageFile> -o <output filename>\n"
	"-cubicspline: applies the cubic spline interpolation method.\n"
	"-nn: applies the nearest neighbor interpolation method.\n");
	exit(0);
}

int main(int argc, char **argv)
{
	char inputImageFile[1024];
	char transFile[1024];
	char outputfile[1024];
	int nx,ny,nz;
	int nx2=0,ny2=0,nz2=0;
	float dx,dy,dz;
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
				opt_i=YES;
				break;
			case 'o':
				sprintf(outputfile,"%s",optarg);
				opt_o=YES;
				break;
			case '7':
				opt_cubicspline=YES;
				break;
			case '8':
				opt_nn=YES;
				break;
			case '?':
				print_help_and_exit();
		}
	}

	if(!opt_i || !opt_o || !opt_T)
		print_help_and_exit();

   if(opt_v)
   {
      printf("Transformation matrix file = %s\n",transFile);
      printf("Input image file = %s\n",inputImageFile);
      printf("Output file = %s\n",outputfile);
   }

	if(opt_cubicspline)
		printf("Applying the cubic spline interpolation method.\n");
	else if(opt_nn)
		printf("Applying the nearest neighbor interpolation method.\n");
	else
		printf("Applying the trilinear interpolation method.\n");

	//////////////////////////////////////////////////////////////////////////////////////////////////
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
	//////////////////////////////////////////////////////////////////////////////////////////////////

	loadTransformation( transFile, T);

	if(opt_v)
	{
		printf("\nTransformation parameters:\n");
		printf("%6.3f  %6.3f  %6.3f  %6.3f\n",T[0],T[1],T[2],T[3]);
		printf("%6.3f  %6.3f  %6.3f  %6.3f\n",T[4],T[5],T[6],T[7]);
		printf("%6.3f  %6.3f  %6.3f  %6.3f\n",T[8],T[9],T[10],T[11]);
		printf("%6.3f  %6.3f  %6.3f  %6.3f\n",T[12],T[13],T[14],T[15]);
		printf("\n");
	}

	invT=inv4(T);

   {
      DIM dim;

      im_in = readNiftiImage(inputImageFile, &dim, 1);
      nx = dim.nx; ny = dim.ny; nz = dim.nz;
      dx = dim.dx; dy = dim.dy; dz = dim.dz;
   }

	// to deal with the sometimes -tive voxel dimensions in SPM/FSL data sets
	if(dx<0.0) dx *= -1.0; if(dy<0.0) dy *= -1.0; if(dz<0.0) dz *= -1.0;

	if(!opt_nx) nx2=nx;
	if(!opt_ny) ny2=ny;
	if(!opt_nz) nz2=nz;
	if(!opt_dx) dx2=dx;
	if(!opt_dy) dy2=dy;
	if(!opt_dz) dz2=dz;

	if(opt_v)
	{
		printf("\nResliced image dimensions:\n");
		printf("Matrix size = %d x %d x %d\n",nx2,ny2,nz2);
		printf("Voxel size = %f x %f x %f\n",dx2,dy2,dz2);
	}

   if(opt_cubicspline)
      im_out=resliceImageCubicSpline(im_in, nx, ny, nz, dx, dy, dz, nx2, ny2, nz2, dx2, dy2, dz2, invT);
   else if(opt_nn)
      im_out=resliceImage(im_in, nx, ny, nz, dx, dy, dz, nx2, ny2, nz2, dx2, dy2, dz2, invT, NEARN);
   else
   {
      im_out=resliceImage(im_in, nx, ny, nz, dx, dy, dz, nx2, ny2, nz2, dx2, dy2, dz2, invT, LIN);
   }

   {
      nifti_1_header hdr;
      hdr = read_NIFTI_hdr(inputImageFile);
      hdr.dim[0]=3; hdr.dim[1]=nx2; hdr.dim[2]=ny2; hdr.dim[3]=nz2;
      hdr.pixdim[1]=dx2; hdr.pixdim[2]=dy2; hdr.pixdim[3]=dz2; 
      save_nifti_image(outputfile, im_out, &hdr);
   }
}
