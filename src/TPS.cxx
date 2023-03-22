#define _TPS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <volume.h>
#include <ctype.h>
#include <nifti1_io.h>
#include <niftiimage.h>
#include <babak_lib.h>
#include <minmax.h>
#include <interpolator.h>

#define YES 1
#define NO 0

int opt;

static struct option options[] =
{
   {"-nx",1,'x'},
   {"-ny",1,'y'},
   {"-nz",1,'z'},
   {"-dx",1,'X'},
   {"-dy",1,'Y'},
   {"-dz",1,'Z'},
   {"-version",0,'V'},
   {"-v",0,'v'},
   {"-lmfg", 1, '1'},
   {"-lmxy", 1, '2'},
   {"-i",1,'i'},
   {"-T", 1, 'T'},
   {"-output-orient",1,'u'},
   {"-oo",1,'u'},
   {0,0,0}
/*
   {"-center-AC", 0, 'M'},
   {"-nn",0,'n'}, // does nearest neighbor interpolation
   {"-standard", 0, 'S'},
   {"-s", 0, 'S'},
   {"-no-tilt-correction", 0, 'R'},
   {"-noppm",0,'N'},
   {"-nopng",0,'g'},
   {"-notxt",0,'t'},
   {"-rvsps",1,'0'},
   {"-rac",1,'1'},
   {"-rpc",1,'2'},
   {"-input-orient",1,'O'}, // secret option
   {"-help",0,'h'},
*/
};

//////////////////////////////////////////////////////////////////////////////////////////////////

void print_help_and_exit()
{
  printf(
  "\nUsage: TSP [options] -i <input NIFTI image> -lmxy <(x,y) landmarks> -lmfg <(f,g) landmarks>\n\n"

  "-i <input NIFTI image>\n"
  "\t3D T1W MRI brain volume in NIFTI format of type short or unsigned short\n\n"

  "Options:\n\n"

  "-v\n"
  "\tEnables verbose mode\n\n"

  "-lm <landmarks-file>\n"
  "\tA text file containing manually determined (i, j, k) coordinates\n"
  "\tof the AC, PC and VSPS, respectively. When this file is supplied,\n"
  "\tautomatic detection of these landmarks is suppressed. This is useful\n" 
  "\tin cases when automatic landmark detection fails.\n\n"

  "-no-tilt-correction\n"
  "\tDoes not tilt-correct the output, but the SFORM and QFORM are set\n"
  "\tcorrectly in the output volume header. This is useful for applications\n"
  "\tthat would like to use acpcdetect as a preprocessing tilt-correction\n"
  "\tstep without applying interpolation at this stage.\n\n"

  "-center-AC\n"
  "\tPlaces the output volume's FOV center at AC\n\n"

  "-standard\n"
  "\tTilt-correction is performed without using the Orion landmarks.\n"
  "\tTilt-correction is done using the AC, PC and MSP only. This is the\n"
  "\tmethod used in version 1.0 of acpcdetect.  In the current version, the\n"
  "\t8 orion landmarks are also used to stabilize the standardization of the\n"
  "\torientation. Using this option, therefore, reverts back to the method\n"
  "\tof version 1.0 without using the additional Orion landmarks.\n\n"

  "-output-orient <orientation-code>\n"
  "\tSpecifies the orientation of the output volume (default: RAS).\n"
  "\tIn ART, orientation codes are 3-letter codes consisting of 6 letters:\n"
  "\tA, P, I, S, L, R.  There are 48 possible combinations. For example\n"
  "\tPIL for Posterior-Inferior-Left or RAS for Right-Anterior-Superior.\n\n"

  "-nx <int>\n"
  "\tNumber of voxels in i direction (the fastest varying index) of the\n"
  "\toutput volume. The default value is determined from the input volume.\n\n"

  "-ny <int>\n"
  "\tNumber of voxels in j direction (the 2nd fastest varying index) of the\n"
  "\toutput volume. The default value is determined from the input volume.\n\n"

  "-nz <int>\n"
  "\tNumber of voxels in k direction (the slowest varying index) of the\n"
  "\toutput volume. The default value is determined from the input volume.\n\n"

  "-dx <float>\n"
  "\tVoxel dimension of the output volume in i direction. The default\n"
  "\tvalue is determined from the input volume.\n\n"

  "-dy <float>\n"
  "\tVoxel dimension of the output volume in j direction. The default\n"
  "\tvalue is determined from the input volume.\n\n"

  "-dz <float>\n"
  "\tVoxel dimension of the output volume in k direction. The default\n"
  "\tvalue is determined from the input volume.\n\n"

  "-version\n"
  "\tPrints software version\n\n"

  "-help\n"
  "\tPrints help information\n\n"

  "-noppm\n"
  "\tPrevents outputting *.ppm images\n\n"

  "-nopng\n"
  "\tPrevents outputting *.png images\n\n"

  "-notxt\n"
  "\tPrevents outputting *.txt files\n\n"

  "-rvsps <r>\n"
  "\tSearch radius for VSPS (default = 50 mm)\n\n"

  "-rac <r>\n"
  "\tSearch radius for AC (default = 15 mm)\n\n"

  "-rpc <r>\n"
  "\tSearch radius for PC (default = 15 mm)\n\n"

  "-nn\n"
  "\tUses the nearest neighbor interpolation for tilt-correction.\n\n"

  "Outputs:\n\n"
  "<output-volume>.nii\n"
  "\tWhere the output volume is saved. The default filename is\n"
  "\t<input-volume>_<output-orientation-code> (default\n"
  "\t<output-orientation-code> is RAS). This volume will be the tilt-corrected\n"
  "\tversion of the input volume. However, if the -no-tilt-correction option is\n"
  "\tselected, the output volume will not be resliced (only reoriented). The\n"
  "\ttilt-correction information, however, are still written in the QFORM and\n"
  "\tSFORM entries of the image header as well as in the *.mrx and *.mat files\n"
  "\t(described below).\n\n"

  "<input-volume>.mrx\n"
  "\tTransformation matrix for tilt-correction in ART format\n\n"

  "<input-volume>_FSL.mat\n"
  "\tTransformation matrix for tilt-correction in FSL format\n\n"

  "<input-volume>_ACPC_sagittal.ppm\n"
  "\tSagittal view of the detected AC/PC locations in\n"
  "\tPPM format (output suppressed by -noppm option)\n\n"

  "<input-volume>_ACPC_sagittal.png\n"
  "\tSagittal view of the detected AC/PC locations in\n"
  "\tPNG format (output suppressed by -nopng option)\n\n"

  "<input-volume>_ACPC_axial.ppm\n"
  "\tAxial view of the detected AC/PC locations in PPM\n"
  "\tformat (output suppressed by -noppm option)\n\n"

  "<input-volume>_ACPC_axial.png\n"
  "\tAxial view of the detected AC/PC locations in PNG\n"
  "\tformat (output suppressed by -nopng option)\n\n"

  "<input-volume>_orion.ppm\n"
  "\tMid-sagittal view of the detected Orion landmarks in\n"
  "\tPPM format (output suppressed by -noppm option)\n\n"

  "<input-volume>_orion.png\n"
  "\tMid-sagittal view of the detected Orion landmarks in\n"
  "\tPNG format (output suppressed by -nopng option)\n\n"

  "<input-volume>_ACPC.txt\n"
  "\tStores the detected AC, PC and VSPS (i, j, k) coordinates and the\n"
  "\testimated mid-sagittal plane (output suppressed by -notxt option)\n\n"

  "<input-volume>_orion.txt\n"
  "\tStores (i, j, k) coordinates of the 8 detected Orion\n"
  "\tlandmarks (output suppressed by -notxt option)\n\n"
  );

  exit(0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
void invert_symmetric_matrix(float4 *X, int p)
{
   float4 *D;
   float4 *Ut;
   float4 *DUt;
   float4 *Xinv;

   // ensure memory was allocated for X and G
   if(X==NULL)
   {
      printf("\n\ninvert_symmetric_matrix(): NULL input matrix encountered.\n\n");
      exit(-1);
   }

   // make sure p is positive
   if(p<=0)
   {
      printf("\n\ninvert_symmetric_matrix(): Non-positive matrix dimension encountered.\n\n");
      exit(-1);
   }

   // Diagnolized X, such that X = U diag(D) U'
   // This function returns U' in place of X
   D=diagATA_float(X, p, 'L');
   Ut=X;

   //printMatrix(Ut,p,p,"Ut",NULL);

   if(D==NULL)
   {
      printf("\n\nninvert_symmetric_matrix(): Could not diagnolize X.\n\n");
      exit(-1);
   }

   DUt=(float4 *)calloc(p*p,sizeof(float4));
   if(DUt==NULL) 
   {
      printf("\n\ninvert_symmetric_matrix(): Memory allocation problem.\n\n");
      exit(-1);
   }

   mat_trans_mat(Ut, p, p, Ut, p, DUt);
   //printMatrix(DUt,p,p,"UUt",NULL);

   // compute diag(1/D) * U'
   for(int i=0; i<p; i++)
   {
      for(int j=0; j<p; j++) DUt[i*p + j]=Ut[i*p +j]/D[i];
   }

   Xinv=(float4 *)calloc(p*p,sizeof(float4));
   if(Xinv==NULL) 
   {
      printf("\n\ninvert_symmetric_matrix(): Memory allocation problem.\n\n");
      exit(-1);
   }

   // Compute Xinv = U * diag(1/D) * U'
   mat_trans_mat(Ut, p, p, DUt, p, Xinv);
   //printMatrix(Xinv,p,p,"Xinv",NULL);

   for(int i=0; i<p*p; i++) X[i]=Xinv[i];

   free(DUt);
   free(D);
   free(Xinv);
}

float4 tpsU(float4 x, float4 y)
{
   double d;
   float4 result;

   d = x*x + y*y;

   if( d > 0.0 )
      result = (float4)(d*log(d));
   else
      result=0.0;

   return(result);
}

float4 *compute_K_inverse(float4 *x, float4 *y, int n)
{
   int d;
   float4 *K;
   float4 *Kinv;
   
   d = n+3;

   K = (float4 *)calloc(d*d, sizeof(float4));
   Kinv = (float4 *)calloc(d*d, sizeof(float4));

   for(int i=0; i<n; i++)
   {
      for(int j=0; j<=i; j++)
      {
         K[i*d + j] = tpsU(x[i]-x[j],y[i]-y[j]);
      }
   }

   { 
      int i1=n;
      int i2=n+1;
      int i3=n+2;
      for(int j=0; j<n; j++)
      {
         K[i1*d + j] = x[j];
         K[i2*d + j] = y[j];
         K[i3*d + j] = 1.0;
      }
   }

   for(int j=0; j<d; j++)
   {
      for(int i=0; i<j; i++)
      {
         K[i*d + j] =  K[j*d + i];
      }
   }

   // copy K in Kinv
   for(int i=0; i<d*d; i++) Kinv[i]=K[i];

   invert_symmetric_matrix(Kinv, d);

   free(K);
   return(Kinv);
}

void read_landmarks(char *filename, int &n, float * &X, float * &Y)
{  
   FILE *fp;
   float dum;
   
   n = 0;
   
   fp=fopen(filename,"r");
   while( fscanf(fp,"%f", &dum) != EOF ) n++;
   fclose(fp);
   
   n /= 3;

   X = (float *)calloc(n+3, sizeof(float) );
   Y = (float *)calloc(n+3, sizeof(float) );

   fp=fopen(filename,"r");
   for(int i=0; i<n; i++) 
   {
      fscanf(fp,"%f", X+i);
      fscanf(fp,"%f", Y+i);
      fscanf(fp,"%f", &dum);
   }
   fclose(fp);

   X[n]=X[n+1]=X[n+2]=0.0;
   Y[n]=Y[n+1]=Y[n+2]=0.0;
}

void tpsTransform(float *X, float *Y, float &x, float &y, float *wx, float *wy, int n, float *k)
{
  for(int i=0; i<n; i++) k[i] = tpsU(X[i]-x,Y[i]-y);
  k[n] = x;
  k[n+1] = y;
  k[n+2] = 1;

  x=0.0;
  for(int i=0; i<n+3; i++) x += k[i]*wx[i];

  y=0.0;
  for(int i=0; i<n+3; i++) y += k[i]*wy[i];
}

int main(int argc, char **argv)
{
  int interpolation_method=LIN;
  float *invT=NULL;
  float PIL2RAS[16];
  float PIL2OPORIENT[16];
  float OPORIENT2PIL[16];
  float Tijk2xyz[16];

  // It is very important to have these initializations.
  int nx=0, ny=0, nz=0;
  float dx=0.0, dy=0.0, dz=0.0;

  char landmarks_fg[DEFAULT_STRING_LENGTH]="";
  float *F=NULL, *G=NULL;  // Arrays that store the (f,g) set of landmarks

  char landmarks_xy[DEFAULT_STRING_LENGTH]="";
  float *X=NULL, *Y=NULL;  // Arrays that store the (x,y) set of landmarks

  int dum;
  int nlm=0;  // number of landmarks 

  float *Kinv;
  float *wx, *wy, *urow;

  char ipimagepath[DEFAULT_STRING_LENGTH]="";
  char ipimagename[DEFAULT_STRING_LENGTH]="";
  char ipimagedir[DEFAULT_STRING_LENGTH]="";  // important to initialize to ""
  short *ipimage;
  nifti_1_header iphdr; // 348 bytes
  char iporient[4]="";
  DIM ipdim;  // input (ip) dimension structures (dim)

  char opimagepath[DEFAULT_STRING_LENGTH]="";
  short *opimage;
  nifti_1_header ophdr; // 348 bytes
  char oporient[4]="PIL"; // default output orientation is PIL
  DIM opdim; 
  opdim.dx = opdim.dy = opdim.dz = 0.0;
  opdim.nx = opdim.ny = opdim.nz = 0;

  char transformation_matrix[DEFAULT_STRING_LENGTH]="";
  float T[16];

  int nPA=0, nLR=0, nIS=0;
  float dPA=0.0, dLR=0.0, dIS=0.0;

  while ((opt = getoption(argc, argv, options)) != -1 )
  {
      switch (opt) 
      {
         case 'u':
            sprintf(oporient,"%s",optarg);
            oporient[0]=(char)toupper((int)oporient[0]);
            oporient[1]=(char)toupper((int)oporient[1]);
            oporient[2]=(char)toupper((int)oporient[2]);
            break;
         case 'x':
            nx = atoi(optarg);
            break;
         case 'y':
            ny = atoi(optarg);
            break;
         case 'z':
            nz = atoi(optarg);
            break;
         case 'X':
            dx = atof(optarg);
            break;
         case 'Y':
            dy = atof(optarg);
            break;
         case 'Z':
            dz = atof(optarg);
            break;
         case 'v':
            opt_v=YES;
            break;
         case '1':
            sprintf(landmarks_fg,"%s",optarg);
            break;
         case '2':
            sprintf(landmarks_xy,"%s",optarg);
            break;
         case 'i':
            sprintf(ipimagepath,"%s",optarg);
            break;
         case 'T':
            sprintf(transformation_matrix,"%s",optarg);
            break;
         case 'V':
            printf("TPS v1.0, March 20, 2023 release\n");
            exit(0);
      }
  }

  /////////////////////////////////////////////////////////////////////////////////

  if( landmarks_fg[0]=='\0' )
  {
    printf("Please specify the (f,g) set of landmarks using: -lmfg <filename>\n");
    exit(0);
  }
  if(opt_v) printf("(f,g) set of landmarks: %s\n",landmarks_fg);

  read_landmarks(landmarks_fg, dum, F, G);
  if(opt_v) 
  {
    printf("Number of landmarks = %d\n",dum);
    for(int i=0; i<dum; i++)
       printf("(f[%d],g[%d]) = (%f, %f)\n",i,i,F[i],G[i]);
    printf("\n");
  }

  /////////////////////////////////////////////////////////////////////////////////

  if( landmarks_xy[0]=='\0' )
  {
    printf("Please specify the (x,y) set of landmarks using: -lmxy <filename>\n");
    exit(0);
  }
  if(opt_v) printf("(x,y) set of landmarks: %s\n",landmarks_xy);

  read_landmarks(landmarks_xy, nlm, X, Y);
  if( nlm != dum)
  {
    printf("Number (x,y) landmarks must equal the number of (f,g) landmarks, aborting ...\n");
    exit(0);
  }

  if(opt_v) 
  {
    printf("Number of landmarks = %d\n",nlm);
    for(int i=0; i<nlm; i++)
       printf("(x[%d],y[%d]) = (%f, %f)\n",i,i,X[i],Y[i]);
    printf("\n");
  }

  /////////////////////////////////////////////////////////////////////////////////

  if( ipimagepath[0]=='\0' )
  {
    printf("Please specify an input image using: -i <inputimage.nii>\n");
    exit(0);
  }
  if(opt_v) printf("Input image: %s\n",ipimagepath);

  getNiftiImageOrientation(ipimagepath, iporient);
  if(opt_v) printf("Orientation: %s\n",iporient);

  ipimage = (short  *)read_nifti_image(ipimagepath, &iphdr);
  if(ipimage==NULL)
  {
    printf("Error reading %s, aborting ...\n", ipimagepath);
    exit(1);
  }

  if( iphdr.datatype != DT_SIGNED_SHORT && iphdr.datatype != DT_UINT16)
  {
    printf("This program can only handle short and unsigned short data types, aborting ...\n");
    exit(0);
  }

  set_dim(ipdim, iphdr);

  if(opt_v) 
  {
    printf("Matrix size: %d x %d x %d\n",ipdim.nx,ipdim.ny,ipdim.nz);
    printf("Voxel size: %6.4f x %6.4f x %6.4f\n",ipdim.dx,ipdim.dy,ipdim.dz);
    printf("scl_slope = %f\n",iphdr.scl_slope);
    printf("scl_inter = %f\n",iphdr.scl_inter);
    printf("\n");
  }

  // set nPA, nIS, nLR, dPA, dIS, dLR depending on iporient
  if(iporient[0]=='P' || iporient[0]=='A') { nPA=ipdim.nx; dPA=ipdim.dx; }
  if(iporient[0]=='I' || iporient[0]=='S') { nIS=ipdim.nx; dIS=ipdim.dx; }
  if(iporient[0]=='L' || iporient[0]=='R') { nLR=ipdim.nx; dLR=ipdim.dx; }
  if(iporient[1]=='P' || iporient[1]=='A') { nPA=ipdim.ny; dPA=ipdim.dy; }
  if(iporient[1]=='I' || iporient[1]=='S') { nIS=ipdim.ny; dIS=ipdim.dy; }
  if(iporient[1]=='L' || iporient[1]=='R') { nLR=ipdim.ny; dLR=ipdim.dy; }
  if(iporient[2]=='P' || iporient[2]=='A') { nPA=ipdim.nz; dPA=ipdim.dz; }
  if(iporient[2]=='I' || iporient[2]=='S') { nIS=ipdim.nz; dIS=ipdim.dz; }
  if(iporient[2]=='L' || iporient[2]=='R') { nLR=ipdim.nz; dLR=ipdim.dz; }

  // determine input image filename without the .nii suffix
  if( niftiFilename(ipimagename, ipimagepath)==0 ) { exit(1); }

  // determine input image directory
  getDirectoryName(ipimagepath, ipimagedir);

  /////////////////////////////////////////////////////////////////////////////////

  // This block of code sets: oporient, opdim, opimagepath, PIL2OPORIENT
  // OPORIENT2PIL, Tijk2xyz, ophdr

  // if ouput orientation is specified using -oo option, make sure it's valid
  if(isOrientationCodeValid(oporient)==0 )
  {
    printf("%s is not a valid orientation code, aborting ...\n",oporient);
    exit(0);
  }

  // setting opdim
  if(oporient[0]=='P' || oporient[0]=='A') { opdim.nx=nPA; opdim.dx=dPA; }
  if(oporient[0]=='I' || oporient[0]=='S') { opdim.nx=nIS; opdim.dx=dIS; }
  if(oporient[0]=='L' || oporient[0]=='R') { opdim.nx=nLR; opdim.dx=dLR; }
  if(oporient[1]=='P' || oporient[1]=='A') { opdim.ny=nPA; opdim.dy=dPA; }
  if(oporient[1]=='I' || oporient[1]=='S') { opdim.ny=nIS; opdim.dy=dIS; }
  if(oporient[1]=='L' || oporient[1]=='R') { opdim.ny=nLR; opdim.dy=dLR; }
  if(oporient[2]=='P' || oporient[2]=='A') { opdim.nz=nPA; opdim.dz=dPA; }
  if(oporient[2]=='I' || oporient[2]=='S') { opdim.nz=nIS; opdim.dz=dIS; }
  if(oporient[2]=='L' || oporient[2]=='R') { opdim.nz=nLR; opdim.dz=dLR; }
  if(nx > 0) opdim.nx=nx; 
  if(ny > 0) opdim.ny=ny; 
  if(nz > 0) opdim.nz=nz;
  if(dx > 0.0) opdim.dx=dx; 
  if(dy > 0.0) opdim.dy=dy; 
  if(dz > 0.0) opdim.dz=dz;
  opdim.nt=1; 
  opdim.dt=0.0; 
  opdim.np=opdim.nx*opdim.ny; 
  opdim.nv=opdim.np*opdim.nz; 

  sprintf(opimagepath,"%s/%s_%s.nii",ipimagedir,ipimagename,oporient);

  if(opt_v) 
  {
    printf("Output image: %s\n",opimagepath);
    printf("Orientation: %s\n",oporient);
    printf("Matrix size: %d x %d x %d\n",opdim.nx,opdim.ny,opdim.nz);
    printf("Voxel size: %6.4f x %6.4f x %6.4f\n",opdim.dx,opdim.dy,opdim.dz);
    printf("scl_slope = %f\n",1.0);
    printf("scl_inter = %f\n",0.0);
    printf("\n");
  }

  // PIL2OPORIENT takes points from PIL space to oporient space
  inversePILtransform(oporient, PIL2OPORIENT);

  // OPORIENT2PIL takes points from oporient space to PIL space
  PILtransform(oporient, OPORIENT2PIL);

  // (i,j,k) -> (x,y,z) in oporient
  ijk2xyz(Tijk2xyz, opdim);

  ophdr = iphdr;
  ophdr.pixdim[1]=opdim.dx; 
  ophdr.pixdim[2]=opdim.dy; 
  ophdr.pixdim[3]=opdim.dz;
  ophdr.dim[0] = 4;
  ophdr.dim[1]=opdim.nx; 
  ophdr.dim[2]=opdim.ny; 
  ophdr.dim[3]=opdim.nz;
  ophdr.dim[4] = 1;
  ophdr.scl_slope = 1.0;
  ophdr.scl_inter = 0.0;
  sprintf(ophdr.descrip,"Created by ART acpcdetect");

  // setup ophdr sform and qform
  float Ttmp[16];
  inversePILtransform("RAS", PIL2RAS); // PIL2RAS takes (x,y,z) points from PIL to RAS space
  multi(OPORIENT2PIL, 4, 4,  Tijk2xyz, 4,  4, Ttmp);
  multi(PIL2RAS, 4, 4,  Ttmp, 4,  4, Ttmp);
  update_qsform(ophdr, Ttmp);

  //////////////////////////////////////////////////////////////////////////////

  T[0]= 1.0;  T[1]= 0.0;  T[2]= 0.0;  T[3]= 0.0;
  T[4]= 0.0;  T[5]= 1.0;  T[6]= 0.0;  T[7]= 0.0;
  T[8]= 0.0;  T[9]= 0.0;  T[10]=1.0;  T[11]=0.0;
  T[12]=0.0;  T[13]=0.0;  T[14]=0.0;  T[15]=1.0;

  if( transformation_matrix[0]!='\0' )
  {
     if(opt_v) printf("Transformation matrix: %s\n",transformation_matrix);
     loadTransformation(transformation_matrix, T);
  }
  else
  {
     if(opt_v) printf("Transformation matrix = Identity\n");
  }

  if(opt_v)
  {
    printf("%6.3f  %6.3f  %6.3f  %6.3f\n",T[0],T[1],T[2],T[3]);
    printf("%6.3f  %6.3f  %6.3f  %6.3f\n",T[4],T[5],T[6],T[7]);
    printf("%6.3f  %6.3f  %6.3f  %6.3f\n",T[8],T[9],T[10],T[11]);
    printf("%6.3f  %6.3f  %6.3f  %6.3f\n",T[12],T[13],T[14],T[15]);
    printf("\n");
  }

  invT = inv4(T);

//  opimage = resliceImage(ipimage,ipdim,opdim,invT,LIN);

  /////////////////////////////////////////////////////////////////////////////////

  Kinv = compute_K_inverse(X, Y, nlm);
  wx=(float *)calloc(nlm+3, sizeof(float));
  wy=(float *)calloc(nlm+3, sizeof(float));
  urow=(float *)calloc(nlm+3, sizeof(float));
  multi(Kinv,nlm+3,nlm+3,F,nlm+3,1,wx);
  multi(Kinv,nlm+3,nlm+3,G,nlm+3,1,wy);

  {
    int nx1, ny1, nz1;
    int nx2, ny2, nz2;
    float dx1, dy1, dz1;
    float dx2, dy2, dz2;

    float  x,y,z;   
    float  xx,yy,zz;   
    int q;
    int np1;
    float xc1, yc1, zc1;
    float xc2, yc2, zc2;
    float *beta=NULL, del;
    float *c=NULL;

/*
for(int i=0; i<nlm; i++)
{
    xx = X[i]; yy=Y[i];
    printf("(f,g) = %f %f\n",F[i],G[i]);
    printf("(x,y) = %f %f\n",xx,yy);
    tpsTransform(X, Y, xx, yy, wx, wy, nlm, urow);
    printf("(x,y) = %f %f\n",xx,yy);
}
*/

    nx1 = iphdr.dim[1]; ny1 = iphdr.dim[2]; nz1 = iphdr.dim[3];
    nx2 = ophdr.dim[1]; ny2 = ophdr.dim[2]; nz2 = ophdr.dim[3];
    dx1 = iphdr.pixdim[1]; dy1 = iphdr.pixdim[2]; dz1 = iphdr.pixdim[3];
    dx2 = ophdr.pixdim[1]; dy2 = ophdr.pixdim[2]; dz2 = ophdr.pixdim[3];

    //printf("Target voxel size: %f %f %f\n",dx2,dy2,dz2);
    //printf("Target matrix size: %d %d %d\n",nx2,ny2,nz2);

    //printf("Subject voxel size: %f %f %f\n",dx1,dy1,dz1);
    //printf("Subject matrix size: %d %d %d\n",nx1,ny1,nz1);

    if(interpolation_method == CUBICSPLINE)
    {
      beta=computeBeta(&del);
      c = (float *)calloc(nx1*ny1*nz1, sizeof(float));
      cubicSplineAnalysis(ipimage, c, nx1, ny1, nz1);
    }

    np1=nx1*ny1;

    opimage=(short *)calloc(nx2*ny2*nz2,sizeof(short));
    if( opimage == NULL )
    {
      memory_allocation_error("opimage");
    }

    xc1=dx1*(nx1-1)/2.0;     /* +---+---+ */
    yc1=dy1*(ny1-1)/2.0;
    zc1=dz1*(nz1-1)/2.0;

    xc2=dx2*(nx2-1)/2.0;     /* +---+---+ */
    yc2=dy2*(ny2-1)/2.0;
    zc2=dz2*(nz2-1)/2.0;

    q=0;
    for(int k=0;k<nz2;k++) 
    {
      for(int j=0;j<ny2;j++) 
      {
         for(int i=0;i<nx2;i++) 
         {
            xx = i*dx2 - xc2;
            yy = j*dy2 - yc2;
            zz = k*dz2 - zc2;
            tpsTransform(X, Y, xx, yy, wx, wy, nlm, urow);

            x = ( invT[0]*xx +invT[1]*yy +invT[2]*zz  +invT[3]   + xc1 )/dx1;
            y = ( invT[4]*xx +invT[5]*yy +invT[6]*zz  +invT[7]   + yc1 )/dy1;
            z = ( invT[8]*xx +invT[9]*yy +invT[10]*zz +invT[11]  + zc1 )/dz1;

            if(interpolation_method == LIN )
            {
              opimage[q++]=(short)(linearInterpolator(x,y,z,ipimage,nx1,ny1,nz1,np1)+0.5);
            }
            else if(interpolation_method == NEARN)
            {
	      opimage[q++]=(short)(nearestNeighbor(x,y,z,ipimage,nx1,ny1,nz1,np1)+0.5);
            }
            else if(interpolation_method == CUBICSPLINE)
            {
               opimage[q++] = (short)(cubicSplineSynthesis(c, nx1, ny1, nz1, x, y, z, beta, del)+0.5);
            }
         }
      }
    }

    if(interpolation_method == CUBICSPLINE)
    {
      free(beta);
      free(c);
    }
  }

  save_nifti_image(opimagepath, opimage, &ophdr);

  free(F), free(G);
  free(X), free(Y);
  free(wx), free(wy), free(urow);
  delete ipimage;
  delete opimage;
  free(invT);
}
