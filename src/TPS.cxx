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
#include "bbk_linear_algebra.h"

#define YES 1
#define NO 0

int opt;

static struct CmdOption options[] =
{
   {"-nx",1,'x'},
   {"-ny",1,'y'},
   {"-nz",1,'z'},
   {"-dx",1,'X'},
   {"-dy",1,'Y'},
   {"-dz",1,'Z'},
   {"-version",0,'V'},
   {"-v",0,'v'},

   {"-lmatlas", 1, '1'},
   {"-lma", 1, '1'},
   {"-la", 1, '1'},

   {"-lmsubject", 1, '2'},
   {"-lms", 1, '2'},
   {"-ls", 1, '2'},

   {"-i",1,'i'},
   {"-a",1,'i'},

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
  "\nUsage: TSP [options] -i <input NIFTI image> -lmsubject <file> -lmatlas <file>\n\n"

  "-lmatlas <file>\n"
  "\tSpecifies a text file containing the \"atlas\" set of landmarks. These are fixed\n"
  "\tlandmarks to which the \"subject\" landmarks are mapped. These landmarks are\n"
  "\tassumed to be in PIL space and are usually obtained using the acpcdetect program.\n\n"

  "-lmsubject <file>\n"
  "\tSpecifies a text file containing the \"subject\" set of landmarks. A TPS deformation\n"
  "\tis found that maps these landmarks to the \"atlas\" landmarks. These landmarks are\n"
  "\tassumed to be in PIL space and are usually obtained using the acpcdetect program.\n\n"

  "-i <input NIFTI image>\n"
  "\t3D T1W MRI brain volume in NIFTI format of type short or unsigned short.\n"
  "\tA TPS deformation will be applied to this image so that the \"subject\"\n"
  "\tlandmarks are precisely mapped to the \"atlas\" landmarks.\n\n"

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

//float *compute_L_inverse(float *x, float *y, int n)
float *compute_L_inverse(float *x, float *y, size_t n)
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

   for(int i=0; i<(int)n; i++)
   {
      for(int j=0; j<=i; j++)
      {
         K[i*n + j] = tpsU(x[i]-x[j],y[i]-y[j]);
      }
   }

   for(int j=0; j<(int)n; j++)
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

   for(int i=0; i<(int)(n*n); i++) 
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

   for(int i=0; i<(int)n; i++)
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
   Ai = inv3x3(A);

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

   for(int i=0; i<(int)n; i++)
   for(int j=0; j<(int)n; j++)
   {
      L[i*d + j] = K[i*n + j];

      Li[i*d + j] = Ki[i*n + j] - C[i*n + j];
   }

   for(int i=0; i<(int)n; i++)
   {
      L[i*d + n] = P[i*3];
      L[i*d + n + 1] = P[i*3 + 1];
      L[i*d + n + 2] = P[i*3 + 2];

      Li[i*d + n] = BAi[i*3];
      Li[i*d + n + 1] = BAi[i*3 + 1];
      Li[i*d + n + 2] = BAi[i*3 + 2];
   }

   for(int j=0; j<(int)n; j++)
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
   if(fp==NULL) { printf("Error: Cound not open %s, aborting ...\n",filename); exit(0); }
   while( fscanf(fp,"%f", &dum) != EOF ) n++;
   fclose(fp);
   
   n /= 3;

   X = (float *)calloc(n+3, sizeof(float) );
   Y = (float *)calloc(n+3, sizeof(float) );

   fp=fopen(filename,"r");
   for(int i=0; i<n; i++) 
   {
      if( fscanf(fp,"%f", X+i) != 1) { 
         fprintf(stderr, "Error: Failed to read X from file.\n");
      }

      if( fscanf(fp,"%f", Y+i) != 1) { 
         fprintf(stderr, "Error: Failed to read X from file.\n");
      }

      if( fscanf(fp,"%f", &dum) != 1) { 
         fprintf(stderr, "Error: Failed to read dum from file.\n");
      }
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
  char atlas_landmarks[DEFAULT_STRING_LENGTH]="";
  float *F=NULL, *G=NULL;  // Arrays that store the (f,g) set of landmarks

  char subject_landmarks[DEFAULT_STRING_LENGTH]="";
  float *X=NULL, *Y=NULL;  // Arrays that store the (x,y) set of landmarks

  int dum;
  int nlm=0;  // number of landmarks 

  char ipimagepath[DEFAULT_STRING_LENGTH]="";

  char oporient[4]="PIL"; // default output orientation is PIL
  DIM opdim; 
  opdim.dx = opdim.dy = opdim.dz = 0.0;
  opdim.nx = opdim.ny = opdim.nz = 0;

  char transformation_matrix[DEFAULT_STRING_LENGTH]="";

  while ((opt = getoption(argc, argv, options)) != -1 )
  {
      switch (opt) 
      {
         case 'h':
            print_help_and_exit();
            break;
         case 'u':
            snprintf(oporient,sizeof(oporient),"%s",optArg);
            oporient[0]=(char)toupper((int)oporient[0]);
            oporient[1]=(char)toupper((int)oporient[1]);
            oporient[2]=(char)toupper((int)oporient[2]);
            break;
         case 'v':
            opt_v=YES;
            break;
         case '1':
            snprintf(atlas_landmarks,sizeof(atlas_landmarks),"%s",optArg);
            break;
         case '2':
            snprintf(subject_landmarks,sizeof(subject_landmarks),"%s",optArg);
            break;
         case 'i':
            snprintf(ipimagepath,sizeof(ipimagepath),"%s",optArg);
            break;
         case 'T':
            snprintf(transformation_matrix,sizeof(transformation_matrix),"%s",optArg);
            break;
         case 'V':
            printf("TPS v1.0, March 20, 2023 release\n");
            exit(0);
      }
  }

  /////////////////////////////////////////////////////////////////////////////////

  if( atlas_landmarks[0]=='\0' )
  {
    printf("Please specify the atlas landmarks using: -lmatlas <file>\n");
    exit(0);
  }
  if(opt_v) printf("Atlas landmarks: %s\n",atlas_landmarks);

  read_landmarks(atlas_landmarks, dum, F, G);
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
    printf("Please specify the subject landmarks using: -lmsubject <file>\n");
    exit(0);
  }
  if(opt_v) printf("Subject landmarks: %s\n",subject_landmarks);

  read_landmarks(subject_landmarks, nlm, X, Y);
  if( nlm != dum)
  {
    printf("Number of subject landmarks is not equal the number of atlas landmarks, aborting ...\n");
    exit(0);
  }

  if(opt_v) 
  {
    printf("Number of landmarks = %d\n",nlm);
    for(int i=0; i<nlm; i++)
       printf("(%f, %f)\n",X[i],Y[i]);
    printf("\n");
  }
}
