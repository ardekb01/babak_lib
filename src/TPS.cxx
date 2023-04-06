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
   {"-lmtarget", 1, '1'},
   {"-lmsubject", 1, '2'},
   {"-i",1,'i'},
   {"-T", 1, 'T'},
   {"-output-orient",1,'u'},
   {"-oo",1,'u'},
   {"-help",0,'h'},
   {"-h",0,'h'},
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
*/
};

//////////////////////////////////////////////////////////////////////////////////////////////////

void print_help_and_exit()
{
  printf(
  "\nUsage: TSP [options] -i <input NIFTI image> -lmsubject <file> -lmtarget <file>\n\n"

  "-lmtarget <file>\n"
  "\tSpecifies a text file containing the \"target\" set of landmarks. These are fixed\n"
  "\tlandmarks to which the \"subject\" landmarks are mapped. These landmarks are\n"
  "\tassumed to be in PIL space and are usually obtained using the acpcdetect program.\n\n"

  "-lmsubject <file>\n"
  "\tSpecifies a text file containing the \"subject\" set of landmarks. A TPS deformation\n"
  "\tis found that maps these landmarks to the \"target\" landmarks. These landmarks are\n"
  "\tassumed to be in PIL space and are usually obtained using the acpcdetect program.\n\n"

  "-i <input NIFTI image>\n"
  "\t3D T1W MRI brain volume in NIFTI format of type short or unsigned short.\n"
  "\tA TPS deformation will be applied to this image so that the \"subject\"\n"
  "\tlandmarks are precisely mapped to the \"target\" landmarks.\n\n"

  "Options:\n\n"

  "-v\n"
  "\tEnables verbose mode\n\n"

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
void invert_symmetric_matrix(float *X, int p)
{
   float *D;
   float *Ut;
   float *DUt;
   float *Xinv;

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

   //printMatrix(D,p,1,"D",NULL);

   if(D==NULL)
   {
      printf("\n\nninvert_symmetric_matrix(): Could not diagnolize X.\n\n");
      exit(-1);
   }

   DUt=(float *)calloc(p*p,sizeof(float));
   if(DUt==NULL) 
   {
      printf("\n\ninvert_symmetric_matrix(): Memory allocation problem.\n\n");
      exit(-1);
   }

   // compute diag(1/D) * U'
   for(int i=0; i<p; i++)
   {
      for(int j=0; j<p; j++) DUt[i*p + j]=Ut[i*p +j]/D[i];
   }

   //mat_trans_mat(Ut, p, p, DUt, p, DUt);
   //printMatrix(DUt,p,p,"UDUt",NULL);

   Xinv=(float *)calloc(p*p,sizeof(float));
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

float tpsU(float x, float y)
{
   double d;
   float result;

   d = x*x + y*y;

   if( d > 0.0 )
      result = (float)(d*log(d));
   else
      result=0.0;

   return(result);
}

float *compute_L_inverse(float *x, float *y, int n)
{
   float *L;
   float *Li;
   
   int d;
   d = n+3;
   L = (float *)calloc(d*d, sizeof(float));
   Li = (float *)calloc(d*d, sizeof(float));

   ////////////////////////////////////////////////
   // set K
   ////////////////////////////////////////////////
   float *K;  // nxn

   K = (float *)calloc(n*n, sizeof(float));

   for(int i=0; i<n; i++)
   {
      for(int j=0; j<=i; j++)
      {
         K[i*n + j] = tpsU(x[i]-x[j],y[i]-y[j]);
      }
   }

   for(int j=0; j<n; j++)
   {
      for(int i=0; i<j; i++)
      {
         K[i*n + j] =  K[j*n + i];
      }
   }

   //printMatrix(K,n,n,"K:",NULL);  // for testing only
   ////////////////////////////////////////////////

   ////////////////////////////////////////////////
   // compute Ki=inverse(K)
   ////////////////////////////////////////////////
   float *Ki; // inverse(K)  nxn

   Ki = (float *)calloc(n*n, sizeof(float));

   for(int i=0; i<n*n; i++) 
   {
      Ki[i]=K[i]; // copy K in Ki temporarily
   }

   invert_symmetric_matrix(Ki, n);

   //multi(Ki,n,n,K,n,n,Ki);  // for testing only
   //printMatrix(Ki,n,n,"K*invserse(K):",NULL);  // for testing only
   ////////////////////////////////////////////////

   ////////////////////////////////////////////////
   // set P and Pt=transpose(P)
   ////////////////////////////////////////////////
   float *P;  // nx3
   float *Pt; // transpose(P) 3xn

   P = (float *)calloc(n*3, sizeof(float));

   for(int i=0; i<n; i++)
   {
      P[3*i] = 1.0;
      P[3*i + 1] = x[i];
      P[3*i + 2] = y[i];
   }

   Pt = trans(P, n, 3);

   //printMatrix(P,n,3,"P:",NULL);  // for testing only
   //printMatrix(Pt,3,n,"Pt:",NULL);  // for testing only
   ////////////////////////////////////////////////

   //////////////////////////////////////////////////////
   // compute B and Bt=transpose(B) where B=inverse(K)*P
   //////////////////////////////////////////////////////
   float *B; // B=inverse(K)*P
   float *Bt; // transpose(B)

   B = (float *)calloc(n*3, sizeof(float));

   multi(Ki,n,n,P,n,3,B); // compute B=inverse(K)*P 
   Bt = trans(B, n, 3);

   //printMatrix(B,n,3,"B:",NULL);  // for testing only
   //printMatrix(Bt,3,n,"Bt:",NULL);  // for testing only
   //////////////////////////////////////////////////////

   //////////////////////////////////////////////////////
   // compute A and Ai=inverse(A) 
   // where A=transpose(P)*inverse(K)*P
   //////////////////////////////////////////////////////
   float *A;
   float *Ai;

   A = (float *)calloc(3*3, sizeof(float));

   multi(Pt,3,n,B,n,3,A); // compute PtKiP=transpose(P)inverse(K)*P 
   Ai = inv3(A);

   //printMatrix(A,3,3,"A:",NULL);  // for testing only
   //printMatrix(Ai,3,3,"Ai:",NULL);  // for testing only
   //multi(Ai,3,3,A,3,3,Ai);  // for testing only
   //printMatrix(Ai,3,3,"A*invserse(A):",NULL);  // for testing only
   //////////////////////////////////////////////////////

   //////////////////////////////////////////////////////
   // compute BAi = B*inverse(A)
   // and AiBt = inverse(A)*transpose(B)
   //////////////////////////////////////////////////////
   float *BAi;
   float *AiBt;

   BAi = (float *)calloc(n*3, sizeof(float));
   multi(B,n,3,Ai,3,3,BAi);
   AiBt = trans(BAi, n, 3);

   //printMatrix(BAi,n,3,"BAi:",NULL);  // for testing only
   //printMatrix(AiBt,3,n,"AiBt:",NULL);  // for testing only
   //////////////////////////////////////////////////////

   //////////////////////////////////////////////////////
   // compute C
   //////////////////////////////////////////////////////
   float *C;
   C = (float *)calloc(n*n, sizeof(float));

   multi(BAi,n,3,Bt,3,n,C);

   //printMatrix(C,n,n,"C:",NULL);  // for testing only
   //////////////////////////////////////////////////////

   for(int i=0; i<n; i++)
   for(int j=0; j<n; j++)
   {
      L[i*d + j] = K[i*n + j];

      Li[i*d + j] = Ki[i*n + j] - C[i*n + j];
   }

   for(int i=0; i<n; i++)
   {
      L[i*d + n] = P[i*3];
      L[i*d + n + 1] = P[i*3 + 1];
      L[i*d + n + 2] = P[i*3 + 2];

      Li[i*d + n] = BAi[i*3];
      Li[i*d + n + 1] = BAi[i*3 + 1];
      Li[i*d + n + 2] = BAi[i*3 + 2];
   }

   for(int j=0; j<n; j++)
   {
      L[n*d + j] = Pt[j];
      L[(n+1)*d + j] = Pt[n + j];
      L[(n+2)*d + j] = Pt[2*n + j];

      Li[n*d + j] = AiBt[j];
      Li[(n+1)*d + j] = AiBt[n + j];
      Li[(n+2)*d + j] = AiBt[2*n + j];
   }

   for(int i=0; i<3; i++)
   for(int j=0; j<3; j++)
   {
      Li[(i+n)*d + j+n] = -Ai[i*3 + j];
   }

   //printMatrix(L,d,d,"L:",NULL);  // for testing only
   //printMatrix(Li,d,d,"Li:",NULL);  // for testing only
   //multi(L,d,d,Li,d,d,Li);
   //printMatrix(Li,d,d,"I:",NULL);  // for testing only

   free(K); free(Ki);
   free(P); free(Pt);
   free(B); free(Bt);
   free(A); free(Ai);
   free(BAi); free(AiBt);
   free(C);
   free(L);

   return(Li);
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
  k[n] = 1;
  k[n+1] = x;
  k[n+2] = y;

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

  char target_landmarks[DEFAULT_STRING_LENGTH]="";
  float *F=NULL, *G=NULL;  // Arrays that store the (f,g) set of landmarks

  char subject_landmarks[DEFAULT_STRING_LENGTH]="";
  float *X=NULL, *Y=NULL;  // Arrays that store the (x,y) set of landmarks

  int dum;
  int nlm=0;  // number of landmarks 

  float *Linv;
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
         case 'h':
            print_help_and_exit();
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
            sprintf(target_landmarks,"%s",optarg);
            break;
         case '2':
            sprintf(subject_landmarks,"%s",optarg);
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

  if( target_landmarks[0]=='\0' )
  {
    printf("Please specify the target set of landmarks using: -lmtarget <file>\n");
    exit(0);
  }
  if(opt_v) printf("Target set of landmarks: %s\n",target_landmarks);

  read_landmarks(target_landmarks, dum, F, G);
  if(opt_v) 
  {
    printf("Number of landmarks = %d\n",dum);
    for(int i=0; i<dum; i++)
       printf("(%f, %f)\n",F[i],G[i]);
    printf("\n");
  }

  /////////////////////////////////////////////////////////////////////////////////

  if( subject_landmarks[0]=='\0' )
  {
    printf("Please specify the subject set of landmarks using: -lmsubject <file>\n");
    exit(0);
  }
  if(opt_v) printf("Subject set of landmarks: %s\n",subject_landmarks);

  read_landmarks(subject_landmarks, nlm, X, Y);
  if( nlm != dum)
  {
    printf("Number subject landmarks must equal the number of target landmarks, aborting ...\n");
    exit(0);
  }

  if(opt_v) 
  {
    printf("Number of landmarks = %d\n",nlm);
    for(int i=0; i<nlm; i++)
       printf("(%f, %f)\n",X[i],Y[i]);
    printf("\n");
  }

/*
  //test segment for observing wx and wy
  Linv = compute_L_inverse(X, Y, nlm);
  wx=(float *)calloc(nlm+3, sizeof(float));
  wy=(float *)calloc(nlm+3, sizeof(float));
  urow=(float *)calloc(nlm+3, sizeof(float));
  multi(Linv,nlm+3,nlm+3,F,nlm+3,1,wx);
  multi(Linv,nlm+3,nlm+3,G,nlm+3,1,wy);
  for(int i=0; i<nlm+3; i++) printf("%f, %f\n",wx[i],wy[i]);
exit(0);
*/
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

  sprintf(opimagepath,"%s/%s_TPS.nii",ipimagedir,ipimagename);

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

  Linv = compute_L_inverse(X, Y, nlm);
  wx=(float *)calloc(nlm+3, sizeof(float));
  wy=(float *)calloc(nlm+3, sizeof(float));
  urow=(float *)calloc(nlm+3, sizeof(float));
  multi(Linv,nlm+3,nlm+3,F,nlm+3,1,wx);
  multi(Linv,nlm+3,nlm+3,G,nlm+3,1,wy);

  {
    int nxi, nyi, nzi;
    int nxo, nyo, nzo;
    float dxi, dyi, dzi;
    float dxo, dyo, dzo;

    float  x,y,z;   
    float  xx,yy,zz;   
    int q;
    int npi;
    float xci, yci, zci;
    float xco, yco, zco;
    float *beta=NULL, del;
    float *c=NULL;

//test code to make sure that the subject landmarks are indeed precisely
//mapped to the target landmarks
for(int i=0; i<nlm; i++)
{
    xx = X[i]; yy=Y[i];
    printf("Target = %f %f\n",F[i],G[i]);
    printf("Subject = %f %f\n",xx,yy);
    tpsTransform(X, Y, xx, yy, wx, wy, nlm, urow);
    printf("(x,y) = %f %f\n",xx,yy);
}
//exit(0);

    nxi = iphdr.dim[1]; nyi = iphdr.dim[2]; nzi = iphdr.dim[3];
    nxo = ophdr.dim[1]; nyo = ophdr.dim[2]; nzo = ophdr.dim[3];
    dxi = iphdr.pixdim[1]; dyi = iphdr.pixdim[2]; dzi = iphdr.pixdim[3];
    dxo = ophdr.pixdim[1]; dyo = ophdr.pixdim[2]; dzo = ophdr.pixdim[3];

    if(interpolation_method == CUBICSPLINE)
    {
      beta=computeBeta(&del);
      c = (float *)calloc(nxi*nyi*nzi, sizeof(float));
      cubicSplineAnalysis(ipimage, c, nxi, nyi, nzi);
    }

    npi=nxi*nyi;

    opimage=(short *)calloc(nxo*nyo*nzo,sizeof(short));
    if( opimage == NULL )
    {
      memory_allocation_error("opimage");
    }

    xci=dxi*(nxi-1)/2.0;     /* +---+---+ */
    yci=dyi*(nyi-1)/2.0;
    zci=dzi*(nzi-1)/2.0;

    xco=dxo*(nxo-1)/2.0;     /* +---+---+ */
    yco=dyo*(nyo-1)/2.0;
    zco=dzo*(nzo-1)/2.0;

    q=0;
    for(int k=0;k<nzo;k++) 
    {
      for(int j=0;j<nyo;j++) 
      {
         for(int i=0;i<nxo;i++) 
         {
            xx = i*dxo - xco;
            yy = j*dyo - yco;
            zz = k*dzo - zco;

            // here xx and yy are updated using TPS
            tpsTransform(X, Y, xx, yy, wx, wy, nlm, urow);

            x = ( invT[0]*xx +invT[1]*yy +invT[2]*zz  +invT[3]   + xci )/dxi;
            y = ( invT[4]*xx +invT[5]*yy +invT[6]*zz  +invT[7]   + yci )/dyi;
            z = ( invT[8]*xx +invT[9]*yy +invT[10]*zz +invT[11]  + zci )/dzi;

            if(interpolation_method == LIN )
            {
              opimage[q++]=(short)(linearInterpolator(x,y,z,ipimage,nxi,nyi,nzi,npi)+0.5);
            }
            else if(interpolation_method == NEARN)
            {
	      opimage[q++]=(short)(nearestNeighbor(x,y,z,ipimage,nxi,nyi,nzi,npi)+0.5);
            }
            else if(interpolation_method == CUBICSPLINE)
            {
               opimage[q++] = (short)(cubicSplineSynthesis(c, nxi, nyi, nzi, x, y, z, beta, del)+0.5);
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
