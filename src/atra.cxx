
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>  //      required by time(), stat()
#include <time.h>       //      required by time()
#include <sys/stat.h>   //      required by stat() 
#include <unistd.h>
#include <spm_analyze.h>
#include <babak_lib.h>
#include <sph.h>
#include <landmarks.h>
#include <minmax.h>
#include <ctype.h>

#define YES 1
#define NO 0

#define DEFAULT_SEARCH_RADIUS 3
#define DEFAULT_PATCH_RADIUS 7
#define MAX_RADIUS 100
#define MAXITER 25  // default maxiter

/////////////////////////////////////
// Global variables required by ATRA
int PILcloudthreshold=50;
int del=15;
int patch_radius=DEFAULT_PATCH_RADIUS;
int search_radius=DEFAULT_SEARCH_RADIUS;
int maxiter=MAXITER;
int opt_nn=NO;
/////////////////////////////////////

int opt;

static struct option options[] =
{
   {"-nn",0,'e'},
   {"-orient",1,'O'},
   {"-o",1,'o'},
   {"-nx",1,'x'},
   {"-ny",1,'y'},
   {"-nz",1,'z'},
   {"-dx",1,'X'},
   {"-dy",1,'Y'},
   {"-dz",1,'Z'},
   {"-noppm",0,'N'},
   {"-nopng",0,'g'},
   {"-notxt",0,'n'},
   {"-v", 0, 'v'},
   {"-imlist", 1, 'i'},
   {"-i", 1, 'i'},
   {"-maxiter", 1, 'm'},
   {"-del", 1, 'd'},
   {"-threshold", 1, 't'},
   {"-r", 1, 'r'},
   {"-R", 1, 'R'},
   {"-h", 0, 'h'},
   {"-help", 0, 'h'},
   {0, 0, 0}
};

void print_help_and_exit()
{
   printf("\nUsage: atra [-v -nx <n> -ny <n> -nz <n> -dx <f> -dy <f> -dz <f> -orient <code> -o <OutputPrefix> -nopng -noppm -notxt -nn] -i <volume list>\n\n"
   "Required:\n"
   "\t-i <volume list>: A text file containing the list of 3D T1W volumes to be registered.\n"
   "\tThe input volumes are required to be in NIFTI-1 format of type 'short int'.\n\n"
   "Options:\n"
   "\t-nopng Prevents outputting *.png images\n\n" 
   "\t-noppm Prevents outputting *.ppm images\n\n" 
   "\t-notxt Prevents outputting *ACPC.txt\n\n" 
   "\t-nn Uses the Nearest Neighbor inperpolation for output images\n\n"
   "\t-v Enables verbose mode\n\n" 
   "\t-nx <int> -ny <int> -nz <int>: Output matrix size (default: 255x255x189)\n\n"
   "\t-dx <float> -dy <float> -dz <float>: Output voxel size (default: 1.0x1.0x1.0 mm^3)\n\n"
   "\t-orient <code>: Output orientation (default <code> is PIL (Posterior-Inferior-Left))\n"
   "\t<code> is a string of 3 letters with 48 possibilities, e.g., RAS (Right-Anterior-Superior).\n\n"
   "\t-o <OutputPrefix>: Prefix used for saving some output files (default=atra).\n\n"
   );

   exit(0);
}

// returns 1 if all images have the same dimensions nx, ny, and nz, 0 otherwise
int checkDimension_avgImage(int N, char **imagefile, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz)
{
   nifti_1_header hdr;
   short dataType;

   if(N==0) return(1);

   printf("Image %d: %s\n",1,imagefile[0]);

   hdr = read_NIFTI_hdr(imagefile[0]);
   *nx = hdr.dim[1];
   *ny = hdr.dim[2];
   *nz = hdr.dim[3];
   *dx = hdr.pixdim[1];
   *dy = hdr.pixdim[2];
   *dz = hdr.pixdim[3];
   dataType = hdr.datatype;

   for(int i=1; i<N; i++)
   {
      printf("Image %d: %s\n",i+1,imagefile[i]);
      hdr = read_NIFTI_hdr(imagefile[i]);

      if( *nx != hdr.dim[1] ||  *ny != hdr.dim[2] ||  *nz != hdr.dim[3]) 
      {
            return(0);
      }
   }

   return(1);
}

void compute_landmark_cm(int nim, int lm_in[][3], int lmcm[], char loomsk[])
{
   float x=0.0, y=0.0, z=0.0;
   int n=0;

   for(int i=0; i<nim; i++)
   if( loomsk[i]==1 )
   {
      x += lm_in[i][0];
      y += lm_in[i][1];
      z += lm_in[i][2];
      n++;
   }

   if( n>0 )
   {
      lmcm[0] = (int)(x/n + 0.5);
      lmcm[1] = (int)(y/n + 0.5);
      lmcm[2] = (int)(z/n + 0.5);
   }
}

// The average of rows of A is returned in sph.v
void extract_patches(int nim, SPH &sph, SHORTIM im[], int lm[][3], float *A)
{
   for(int i=0; i<nim; i++)
   {
      sph.set(im[i], lm[i][0], lm[i][1], lm[i][2]); 
      sph.zeromean();
      sph.normalize();
      sph.get(A+i*sph.n);
   }

   avgRow(A, nim, sph.n, sph.v);
}

// The average of rows of A is returned in sph.v
void extract_patches(int nim, SPH &sph, SHORTIM im[], int lm_in[][3], float *A, char loomsk[])
{
   int nrow=0;

   for(int i=0; i<nim; i++)
   if( loomsk[i]==1 )
   {
      sph.set(im[i], lm_in[i][0], lm_in[i][1], lm_in[i][2]); 
      sph.zeromean();
      sph.normalize();
      sph.get(A + nrow*sph.n);
      nrow++;
   }

   avgRow(A, nrow, sph.n, sph.v);
}

int lm_distance(int lm1[][3], int lm2[][3], int n)
{
   int d=0;

   for(int i=0; i<n; i++)
   {
      if( lm1[i][0] > lm2[i][0] ) d += (lm1[i][0]-lm2[i][0]);
      else d += (lm2[i][0]-lm1[i][0]);

      if( lm1[i][1] > lm2[i][1] ) d += (lm1[i][1]-lm2[i][1]);
      else d += (lm2[i][1]-lm1[i][1]);

      if( lm1[i][2] > lm2[i][2] ) d += (lm1[i][2]-lm2[i][2]);
      else d += (lm2[i][2]-lm1[i][2]);
   }

   return(d);
}

void lm_copy(int lm1[][3], int lm2[][3], int n)
{
   for(int i=0; i<n; i++)
   {
      lm2[i][0]=lm1[i][0];
      lm2[i][1]=lm1[i][1];
      lm2[i][2]=lm1[i][2];
   }
}

// nim: number of images in the training set
// refsph: spherical patch holder
int seek_lm(int nim, SPH &refsph, SHORTIM im[], int lm_in[][3], float *A, float *ssd, char loomsk[], int lmcm[],
SPH &searchsph, SPH &testsph, int lm_out[][3])
{
   int d=1;  // compiler wants this initialized
   int search_center[3];

   search_center[0]=lm_in[0][0];
   search_center[1]=lm_in[0][1];
   search_center[2]=lm_in[0][2];

   for(int i=0; i<nim; i++) loomsk[i]=1;

   for(int iter=0; iter<maxiter; iter++)
   {
      float dum=0.0; 
      extract_patches(nim, refsph, im, lm_in, A);
      ssdRow(A, nim, refsph.n, refsph.v, ssd);
      for(int i=0; i<refsph.n; i++) dum += ssd[i];
      ssd[iter]=dum*refsph.n;
      if(iter>=2 && ssd[iter]==ssd[iter-2]) break;

      for(int s=0; s<nim; s++)
      {
         loomsk[s]=0; // leaves s out

         extract_patches(nim, refsph, im, lm_in, A, loomsk);
         detect_lm(searchsph, testsph, im[s], search_center, refsph, lm_out[s]);

         loomsk[s]=1; // add s in 
      }

      d=lm_distance(lm_in, lm_out, nim);
      lm_copy(lm_out, lm_in, nim);
      compute_landmark_cm(nim, lm_in, lmcm, loomsk);
      
      if(d==0) break;
   }

   return(d);
}

void atra(const char *imagelistfile, DIM output_dim, const char *outputOrientationCode, const char *outputPrefix)
{
   char orient[4]="";
   //PIL2OUT maps points from PIL orientation to the output 
   //orientation specified by outputOrientationCode
   float PIL2OUT[16];
   float *ssd;
   int nim=0; // number of volumes listed as input in imagelistfile
   int *seedi, *seedj, *seedk;
   char temporaryFilename[DEFAULT_STRING_LENGTH]; // place holder for temporary filenames
   short *PILbraincloud;
   SHORTIM PILim[MAXIM];
   SHORTIM im[MAXIM];
   DIM PILbraincloud_dim;
   DIM input_dim;
   nifti_1_header PILbraincloud_hdr; 
   nifti_1_header tmp_hdr; 

   FILE *fp;
   char **imagefile; // the nim image files
   char **landmarksfile; // the (potential) nim landmarks files
   char **imagefileprefix;
   char **imagedir;

   float **TPIL; // nim-array of 4x4 transformation matrices which transform each of
                 // the nim images into the standardized PIL system
   float **TGPA; // nim-array of 4x4 transformation matrices obtained from Generalized Procrustes Analysis 
   nifti_1_header imhdr;
   float *invT;
   int nseeds=0;

   SPH refsph(patch_radius);
   SPH testsph(patch_radius);
   SPH searchsph(search_radius);
   int *alignedLM;
   float **Pt;
   float **P;
   int nlm;
   int nlm_nonconvergent;
   int current_nlm_aligned;
   int final_nlm_aligned=0;
   float *A; // (NxM) matrix where N=n; M=sph.n
   char *sameflag;
   int lm_in[MAXIM][3]; // (i,j,k) coordinates of the input landmarks
   int lm_out[MAXIM][3];
   int cbd; //city block distance 
   char loomsk[MAXIM];
   int lmcm[3]; // landmarks center of mass
   float TOUT[16], TOUT_FSL[16];
   float OUT2PIL[16];
   float PIL2RAS[16];
   float OUT2RAS[16];
   float T_ijk2xyz[16];
   float *scalefactor;
   float *avgOutputImage;
   short *avgOutputImageShort;

   //PIL2OUT maps points from PIL orientation to the output 
   //orientation specified by outputOrientationCode
   inversePILtransform(outputOrientationCode, PIL2OUT);

   opt_qform=YES;
   opt_sform=YES;

   // find nim
   fp=fopen(imagelistfile,"r");
   if(fp==NULL) file_open_error(imagelistfile);

   nim=0;
   while(fscanf(fp,"%s",temporaryFilename) != EOF ) 
   {
      if( not_magical_nifti(temporaryFilename,0)==0 ) nim++;
   }
   fclose(fp);

   if(opt_v) printf("Number of volumes to be registered = %d\n", nim);

   if(nim==0) 
   {
      printf("Number of input volumes = %d, I have nothing to do!\n", nim);
      exit(0);
   }

   // allocate memory for various arrays
   // All 5 variables freed before atra() returns
   imagefile = (char **)calloc(nim, sizeof(char *));
   landmarksfile = (char **)calloc(nim, sizeof(char *));
   imagefileprefix = (char **)calloc(nim, sizeof(char *));
   imagedir = (char **)calloc(nim, sizeof(char *));
   TPIL = (float **)calloc(nim, sizeof(float *));
   TGPA = (float **)calloc(nim, sizeof(float *));

   for(int i=0; i<nim; i++)
   {
      // freed before atra() returns
      imagefile[i] = (char *)calloc(DEFAULT_STRING_LENGTH, sizeof(char));
      imagefile[i][0]='\0';

      // freed before atra() returns
      landmarksfile[i] = (char *)calloc(DEFAULT_STRING_LENGTH, sizeof(char));
      landmarksfile[i][0]='\0';

      // freed before atra() returns
      imagefileprefix[i] = (char *)calloc(DEFAULT_STRING_LENGTH, sizeof(char));
      imagefileprefix[i][0]='\0';

      // freed before atra() returns
      imagedir[i] = (char *)calloc(DEFAULT_STRING_LENGTH, sizeof(char));
      imagedir[i][0]='\0';

      // both variables freed before atra() returns
      TPIL[i] = (float *)calloc(16, sizeof(float));
      TGPA[i] = (float *)calloc(16, sizeof(float));
   }

   ////////////////////////////////////////////////////////////////////////////////////////////
   // fill imagefile, landmarksfile and imagefileprefix arrays
   ////////////////////////////////////////////////////////////////////////////////////////////
   fp=fopen(imagelistfile,"r");
   if(fp==NULL) file_open_error(imagelistfile);
   nim=0;
   while(fscanf(fp,"%s",temporaryFilename) != EOF ) 
   {
      if( not_magical_nifti(temporaryFilename,0)==0 ) 
      {
         strcpy(imagefile[nim], temporaryFilename);
         if( niftiFilename(imagefileprefix[nim], imagefile[nim])==0 ) { exit(0); }
         getDirectoryName(imagefile[nim], imagedir[nim]);
         nim++;
      }
      else
      {
         // ensure nim>0 in case they put a lm file first
         if(nim>0) strcpy(landmarksfile[nim-1], temporaryFilename);
      }
   }
   fclose(fp);

   if(opt_v) 
   {
      printf("Input volume list:\n");
      for(int i=0; i<nim; i++) 
      {
         printf("Volume %d: %s\n",i+1,imagefile[i]);
         if(landmarksfile[i][0]!='\0') printf("Volume %d landmarks: %s\n",i+1,landmarksfile[i]);
      }
   }
   ////////////////////////////////////////////////////////////////////////////////////////////

   /////////////////////////////////////////////////////////////////////////////////////////////
   // read PILbraincloud.nii from the $ARTHOME directory
   /////////////////////////////////////////////////////////////////////////////////////////////
   snprintf(temporaryFilename,sizeof(temporaryFilename),"%s/PILbrain.nii",ARTHOME);

   // freed before atra() returns
   PILbraincloud = (short  *)read_nifti_image(temporaryFilename, &PILbraincloud_hdr);

   set_dim(PILbraincloud_dim, PILbraincloud_hdr);

   if(PILbraincloud==NULL)
   {
         printf("Error reading %s, aborting ...\n", temporaryFilename);
         exit(1);
   }
   /////////////////////////////////////////////////////////////////////////////////////////////

   /////////////////////////////////////////////////////////////////////////////////////////////
   // check output_dim
   /////////////////////////////////////////////////////////////////////////////////////////////
   if( output_dim.dx <= 0.0) output_dim.dx=PILbraincloud_dim.dx;
   if( output_dim.dy <= 0.0) output_dim.dy=PILbraincloud_dim.dy;
   if( output_dim.dz <= 0.0) output_dim.dz=PILbraincloud_dim.dz;
   if( output_dim.nx <= 0 || output_dim.nx > 1024) output_dim.nx=PILbraincloud_dim.nx;
   if( output_dim.ny <= 0 || output_dim.ny > 1024) output_dim.ny=PILbraincloud_dim.ny;
   if( output_dim.nz <= 0 || output_dim.nz > 1024) output_dim.nz=PILbraincloud_dim.nz;
   output_dim.np = output_dim.nx*output_dim.ny;
   output_dim.nv = output_dim.np*output_dim.nz;
   if(opt_v)
   {
      printf("Output matrix size: %d x %d x %d\n", output_dim.nx, output_dim.ny, output_dim.nz);
      printf("Output voxel size: %7.5f x %7.5f x %7.5f mm^3\n", output_dim.dx, output_dim.dy, output_dim.dz);
   }

   /////////////////////////////////////////////////////////////////////////////////////////////
   // find TPIL
   /////////////////////////////////////////////////////////////////////////////////////////////
   if(opt_v) printf("PIL transformation:\n");
   for(int i=0; i<nim; i++)
   {
      if(opt_v) printf("Processing %s ...\n",imagefile[i]);
      orient[0]='\0';
      new_PIL_transform(imagefile[i],landmarksfile[i],orient, TPIL[i],1);
   }
   /////////////////////////////////////////////////////////////////////////////////////////////

   for(int i=0; i<nim; i++)
   {
      set_dim(PILim[i], PILbraincloud_hdr);

      //freed before atra() returns
      im[i].v = (short *)read_nifti_image(imagefile[i], &imhdr);
      set_dim(im[i], imhdr);

      invT = inv4(TPIL[i]);
      //freed before atra() returns
      PILim[i].v = resliceImage(im[i].v, im[i].nx, im[i].ny, im[i].nz, im[i].dx, im[i].dy, im[i].dz,
      PILim[i].nx, PILim[i].ny, PILim[i].nz, PILim[i].dx, PILim[i].dy, PILim[i].dz, invT, LIN);
      free(invT);
   }

   if(nim==1)
   {
      short *tmp;

      tmp_hdr = read_NIFTI_hdr(imagefile[0]);
      set_dim(input_dim, tmp_hdr);

      multi(PIL2OUT,4,4,TPIL[0],4,4,TOUT);
      invT = inv4(TOUT);
      tmp = resliceImage(im[0].v,input_dim,output_dim,invT,LIN);
      free(invT);

      set_dim(tmp_hdr, output_dim);
      tmp_hdr.magic[0]='n'; tmp_hdr.magic[1]='+'; tmp_hdr.magic[2]='1';
      strcpy(tmp_hdr.descrip,"Created by ART ATRA program.");
      snprintf(temporaryFilename,sizeof(temporaryFilename),"%s/%s_%s.nii",imagedir[0],imagefileprefix[0],outputOrientationCode);
      save_nifti_image(temporaryFilename, tmp, &tmp_hdr);

      PILtransform(outputOrientationCode, OUT2PIL);
      inversePILtransform("RAS", PIL2RAS);
      ijk2xyz(T_ijk2xyz, output_dim.nx, output_dim.ny, output_dim.nz, output_dim.dx, output_dim.dy, output_dim.dz);
      multi(PIL2RAS, 4, 4,  OUT2PIL, 4,  4, OUT2RAS);
      multi(OUT2RAS, 4, 4,  T_ijk2xyz, 4,  4, OUT2RAS);
      update_qsform(temporaryFilename, OUT2RAS);

      // Yes! saving the same image again under a different name! 
      snprintf(temporaryFilename,sizeof(temporaryFilename),"%s_avg_%s.nii",outputPrefix, outputOrientationCode);
      save_nifti_image(temporaryFilename, tmp, &tmp_hdr);
      update_qsform(temporaryFilename, OUT2RAS);

      free(tmp);

      snprintf(temporaryFilename,sizeof(temporaryFilename),"%s/%s_PIL.mrx",imagedir[0],imagefileprefix[0]);
      remove(temporaryFilename);

      multi(PIL2OUT,4,4,TPIL[0],4,4,TOUT);
      snprintf(temporaryFilename,sizeof(temporaryFilename),"%s/%s.mrx",imagedir[0],imagefileprefix[0]);
      fp=fopen(temporaryFilename,"w");
      if(fp==NULL) file_open_error(temporaryFilename);
      printMatrix(TOUT,4,4,"",fp);
      fclose(fp);

      art_to_fsl(TOUT, TOUT_FSL, input_dim, output_dim);
      snprintf(temporaryFilename,sizeof(temporaryFilename),"%s/%s_FSL.mat",imagedir[0],imagefileprefix[0]);
      fp=fopen(temporaryFilename,"w");
      if(fp==NULL) file_open_error(temporaryFilename);
      printMatrix(TOUT_FSL,4,4,"",fp);
      fclose(fp);

      snprintf(temporaryFilename,sizeof(temporaryFilename),"%s.txt",outputPrefix);
      fp = fopen(temporaryFilename,"w");
      if(fp==NULL) file_open_error(temporaryFilename);
      fprintf(fp,"%d\n",nim);
      fprintf(fp,"%s %f\n",imagefile[0],1.0);
      snprintf(temporaryFilename,sizeof(temporaryFilename),"%s/%s.mrx",imagedir[0],imagefileprefix[0]);
      fprintf(fp,"%s\n",temporaryFilename);
      fclose(fp);

      free(imagefile[0]);
      free(imagefile);
      free(landmarksfile[0]);
      free(landmarksfile);
      free(imagefileprefix[0]);
      free(imagefileprefix);
      free(TPIL[0]);
      free(TPIL);
      free(TGPA[0]);
      free(TGPA);
      free(im[0].v);
      free(PILim[0].v);
      free(PILbraincloud);

      return;
   }

   //////////////////////////////////////////////////////////
   // Compute the number of seeds
   nseeds=0;
   for(int i=0; i<PILbraincloud_dim.nx; i += del)
   for(int j=0; j<PILbraincloud_dim.ny; j += del)
   for(int k=14; k<PILbraincloud_dim.nz; k += del)
   if ( PILbraincloud[k*PILbraincloud_dim.np + j*PILbraincloud_dim.nx + i] > PILcloudthreshold )
   {
      nseeds++;
   }

   // all 5 variables freed before atra() returns
   seedi = (int *)calloc(nseeds, sizeof(int)); 
   seedj = (int *)calloc(nseeds, sizeof(int)); 
   seedk = (int *)calloc(nseeds, sizeof(int)); 
   sameflag = (char *)calloc(nseeds, sizeof(char)); 
   alignedLM = (int *)calloc(nseeds*3, sizeof(int));

   nseeds=0;
   for(int i=0; i<PILbraincloud_dim.nx; i += del)
   for(int j=0; j<PILbraincloud_dim.ny; j += del)
   for(int k=14; k<PILbraincloud_dim.nz; k += del)
   if ( PILbraincloud[k*PILbraincloud_dim.np + j*PILbraincloud_dim.nx + i] > PILcloudthreshold )
   {
      seedi[nseeds]=i;
      seedj[nseeds]=j;
      seedk[nseeds]=k;
      nseeds++;
   }
   //////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////
   // LOOC landmark identification
   //////////////////////////////////////////////////////////
   // All 4 variables free before atra() returns
   ssd = (float *)calloc(refsph.n,sizeof(float));
   A = (float *)calloc(nim*refsph.n,sizeof(float));
   Pt = (float **)calloc(nim, sizeof(float *) );
   P = (float **)calloc(nim, sizeof(float *) );

   // the Pt[i] are actually nlm*3 matrices but we allocate more rows (nseeds rows)
   // since at this time we only know that nlm<=nseeds
   // Pt[i] are freed before atra() returns
   for(int i=0;i<nim;i++) Pt[i] = (float *)calloc(nseeds*3, sizeof(float) );

   if(opt_v) printf("LOOC landmark identification ...\n");
   if(opt_v && maxiter!=MAXITER) printf("Maximum number of allowed iterations per seed = %d\n",maxiter);
   if(opt_v && del!=20) printf("Del = %d\n",del);

   for(int iteration=0; iteration<maxiter; iteration++)
   {
      nlm=0;
      nlm_nonconvergent=0;
      current_nlm_aligned=0;

      for(int s=0; s<nseeds; s++)
      {
         sameflag[s]=0;

         // initialize the landmark location for all images to be the same
         for(int i=0; i<nim; i++)
         {
            lm_in[i][0]=seedi[s];
            lm_in[i][1]=seedj[s];
            lm_in[i][2]=seedk[s];
         }
   
         cbd = seek_lm(nim, refsph, PILim, lm_in, A, ssd, loomsk, lmcm, searchsph, testsph, lm_out);
   
         if(cbd>0) { nlm_nonconvergent++; continue; }
   
         /////////////////////////////////////////////////////////////////////////////////////////////////////
         // check to see if the LOOC converged to the same location in all images
         /////////////////////////////////////////////////////////////////////////////////////////////////////
         sameflag[s]=1;
         for(int i=1; i<nim; i++)
         {
            if(lm_out[i][0] != lm_out[0][0] || lm_out[i][1] != lm_out[0][1] || lm_out[i][2] != lm_out[0][2] )
            { sameflag[s]=0; break; }
         }
   
         if(sameflag[s]) { current_nlm_aligned++;}
         /////////////////////////////////////////////////////////////////////////////////////////////////////

         for(int i=0; i<nim; i++)
         {
            Pt[i][nlm*3 + 0] = lm_out[i][0]; 
            Pt[i][nlm*3 + 1] = lm_out[i][1]; 
            Pt[i][nlm*3 + 2] = lm_out[i][2]; 
         }
   
         nlm++;
      }

      if( current_nlm_aligned <= final_nlm_aligned && iteration>0) break;
      else  //final_nlm_aligned does not get updated, nor is alignedLM
      {
         int incrementflg;

         final_nlm_aligned=0;
         for(int lm=0; lm<nlm; lm++)
         {
            incrementflg=1;
            for(int i=1; i<nim; i++)
            {
               if(Pt[i][lm*3+0] != Pt[0][lm*3+0] || Pt[i][lm*3+1] != Pt[0][lm*3+1] || Pt[i][lm*3+2] != Pt[0][lm*3+2] ) 
               { incrementflg=0; break; }
            }

            if(incrementflg)
            {
               alignedLM[3*final_nlm_aligned+0]=Pt[0][lm*3+0];
               alignedLM[3*final_nlm_aligned+1]=Pt[0][lm*3+1];
               alignedLM[3*final_nlm_aligned+2]=Pt[0][lm*3+2];
               final_nlm_aligned++;
            }
         }
      }

      if(opt_v) printf("Iteration %d ...\n", iteration+1);
      //if(opt_v) printf("Number of seeds = %d\n", nseeds);
      //if(opt_v) printf("Number of LOOC landmarks identified = %d\n",nlm);
      if(opt_v) printf("N* = %d\n",current_nlm_aligned);

      if(nlm<6)
      {
         printf("Warning: insufficient LOOC detected. Skipping the generalized Procrustes anlaysis ...\n");
         break;
      }
      else
      {
         float *Q;
         float dum[3];
         float *Ptmp;
         float *Qtmp;
   
         if(opt_v) printf("Generalized Procrustes anlaysis ...\n");
   
         for(int i=0; i<nim; i++)
         {
            P[i] = (float *)calloc(3*nlm, sizeof(float) );
            transpose_matrix(Pt[i], nlm,  3, P[i]);
            convert_to_xyz(P[i], nlm, PILim[i]);
         }
   
         Q = (float *)calloc(3*nlm, sizeof(float) );
         Qtmp = (float *)calloc(3*nlm, sizeof(float) );
         Ptmp = (float *)calloc(3*nlm, sizeof(float) );
   
         // compute initial Q
         for(int k=0; k<3*nlm; k++)
         {
               Q[k]=0.0;
               for(int i=0; i<nim; i++) Q[k] += P[i][k];
               Q[k]/=nim;
         }

         for(int GPiter=0; GPiter<100; GPiter++)
         {
            for(int i=0; i<nim; i++) 
            {
               for(int k=0; k<3*nlm; k++) 
               {  
                  Ptmp[k]=P[i][k]; 
                  Qtmp[k]=Q[k]; 
               }

               Procrustes(Qtmp, nlm, Ptmp, TGPA[i]);
            }
   
            // update Q: find a new Q as the average of transformed P[i]
            for(int k=0; k<nlm; k++) 
            {
               Q[k]=Q[k+nlm]=Q[k+nlm*2]=0.0;
               for(int i=0; i<nim; i++) 
               {
                  dum[0] = P[i][k];
                  dum[1] = P[i][nlm + k];
                  dum[2] = P[i][2*nlm + k];
   
                  Q[k]         += (dum[0]*TGPA[i][0] + dum[1]*TGPA[i][1] + dum[2]*TGPA[i][2]  + TGPA[i][3]);
                  Q[nlm + k]   += (dum[0]*TGPA[i][4] + dum[1]*TGPA[i][5] + dum[2]*TGPA[i][6]  + TGPA[i][7]);
                  Q[2*nlm + k] += (dum[0]*TGPA[i][8] + dum[1]*TGPA[i][9] + dum[2]*TGPA[i][10] + TGPA[i][11]);
               }
               Q[k]       /= nim;
               Q[k+nlm]   /= nim;
               Q[k+2*nlm] /= nim;
            }
         }

         free(Q);
         free(Ptmp);
         free(Qtmp);
         for(int i=0; i<nim; i++)
         {
            free(P[i]); 
         }
      }
   
      // update TPIL 
      for(int i=0; i<nim; i++)
      {
         multi(TGPA[i],4,4,TPIL[i],4,4,TPIL[i]);
      }

      // update PILim
      for(int i=0; i<nim; i++)
      {
         free(PILim[i].v);
         invT = inv4(TPIL[i]);
         PILim[i].v = resliceImage(im[i].v, im[i].nx, im[i].ny, im[i].nz, im[i].dx, im[i].dy, im[i].dz,
         PILim[i].nx, PILim[i].ny, PILim[i].nz, PILim[i].dx, PILim[i].dy, PILim[i].dz, invT, LIN);
         free(invT);
      }
   }
   // end of iterations
   //////////////////////////////////////////////////////////

   float *mean; 
   int ii,jj,kk,np,ny;

   snprintf(temporaryFilename,sizeof(temporaryFilename),"%s.txt",outputPrefix);
   fp = fopen(temporaryFilename,"w");
   if(fp==NULL) file_open_error(temporaryFilename);
   fprintf(fp,"%d\n",nim);
   mean = (float *)calloc(nim, sizeof(float));
   scalefactor = (float *)calloc(nim, sizeof(float));

   for(int i=0; i<nim; i++) mean[i]=0.0;

   for(int lm=0; lm<final_nlm_aligned; lm++)
   {
      ii=alignedLM[3*lm+0];
      jj=alignedLM[3*lm+1];
      kk=alignedLM[3*lm+2];
      for(int i=0; i<nim; i++)
      {
         np = PILim[i].np;
         ny = PILim[i].ny;
         mean[i] += PILim[i].v[kk*np + jj*ny + ii];
      }
   }

   for(int i=0; i<nim; i++) 
   {
      mean[i]/=final_nlm_aligned;
   }

   float scale=0.0;
   for(int i=0; i<nim; i++) scale += mean[i];
   scale /= nim;
   
   for(int i=0; i<nim; i++) 
   {
      scalefactor[i] = scale/mean[i];
      fprintf(fp,"%s %f\n",imagefile[i],scalefactor[i]);
      snprintf(temporaryFilename,sizeof(temporaryFilename),"%s/%s.mrx",imagedir[i],imagefileprefix[i]);
      fprintf(fp,"%s\n",temporaryFilename);
   }
   fclose(fp);

   free(mean);

   //////////////////////////////////////////////////////////
   // save outputs
   //////////////////////////////////////////////////////////
   avgOutputImage = (float *)calloc(output_dim.nv, sizeof(float));
   avgOutputImageShort = (short *)calloc(output_dim.nv, sizeof(short));

   for(int i=0; i<nim; i++)
   {
      short *tmp;

      tmp_hdr = read_NIFTI_hdr(imagefile[i]);
      set_dim(input_dim, tmp_hdr);

      multi(PIL2OUT,4,4,TPIL[i],4,4,TOUT);
      invT = inv4(TOUT);
      if(opt_nn)
         tmp = resliceImage(im[i].v,input_dim,output_dim,invT,NEARN);
      else
         tmp = resliceImage(im[i].v,input_dim,output_dim,invT,LIN);
      free(invT);

      set_dim(tmp_hdr, output_dim);
      tmp_hdr.magic[0]='n'; tmp_hdr.magic[1]='+'; tmp_hdr.magic[2]='1';
      strcpy(tmp_hdr.descrip,"Created by ART ATRA program.");
      snprintf(temporaryFilename,sizeof(temporaryFilename),"%s/%s_%s.nii",imagedir[i],imagefileprefix[i],outputOrientationCode);
      save_nifti_image(temporaryFilename, tmp, &tmp_hdr);

      PILtransform(outputOrientationCode, OUT2PIL);
      inversePILtransform("RAS", PIL2RAS);
      ijk2xyz(T_ijk2xyz, output_dim.nx, output_dim.ny, output_dim.nz, output_dim.dx, output_dim.dy, output_dim.dz);
      multi(PIL2RAS, 4, 4,  OUT2PIL, 4,  4, OUT2RAS);
      multi(OUT2RAS, 4, 4,  T_ijk2xyz, 4,  4, OUT2RAS);
      update_qsform(temporaryFilename, OUT2RAS);

      for(int v=0; v<output_dim.nv; v++) avgOutputImage[v] += tmp[v]*scalefactor[i];

      free(tmp);

      // save transformation matrix
      snprintf(temporaryFilename,sizeof(temporaryFilename),"%s/%s_PIL.mrx",imagedir[i], imagefileprefix[i]);
      remove(temporaryFilename);

      multi(PIL2OUT,4,4,TPIL[i],4,4,TOUT);
      snprintf(temporaryFilename,sizeof(temporaryFilename),"%s/%s.mrx",imagedir[i],imagefileprefix[i]);
      fp=fopen(temporaryFilename,"w");
      if(fp==NULL) file_open_error(temporaryFilename);
      printMatrix(TOUT,4,4,"",fp);
      fclose(fp);

      art_to_fsl(TOUT, TOUT_FSL, input_dim, output_dim);
      snprintf(temporaryFilename,sizeof(temporaryFilename),"%s/%s_FSL.mat",imagedir[i],imagefileprefix[i]);
      fp=fopen(temporaryFilename,"w");
      if(fp==NULL) file_open_error(temporaryFilename);
      printMatrix(TOUT_FSL,4,4,"",fp);
      fclose(fp);
   }
   //////////////////////////////////////////////////////////

   for(int v=0; v<output_dim.nv; v++) avgOutputImageShort[v] = (short)(avgOutputImage[v]/nim + 0.5);

   snprintf(temporaryFilename,sizeof(temporaryFilename),"%s_avg_%s.nii",outputPrefix, outputOrientationCode);
   save_nifti_image(temporaryFilename, avgOutputImageShort, &tmp_hdr);
   update_qsform(temporaryFilename, OUT2RAS);

   free(avgOutputImage);
   free(avgOutputImageShort);

   // free memory
   for(int i=0; i<nim; i++)
   {
      free(imagefile[i]);
      free(landmarksfile[i]);
      free(imagefileprefix[i]);
      free(imagedir[i]);
      free(TPIL[i]);
      free(TGPA[i]);
      free(PILim[i].v);
      free(im[i].v);
      free(Pt[i]);
   }
   free(imagefile);
   free(landmarksfile);
   free(imagefileprefix);
   free(imagedir);
   free(TPIL);
   free(TGPA);
   free(PILbraincloud);
   free(seedi);
   free(seedj);
   free(seedk);
   free(sameflag);
   free(alignedLM);
   free(ssd);
   free(A);
   free(Pt);
   free(P);
   free(scalefactor);

   return;
}

int main(int argc, char **argv)
{
   time_t time_start, time_end;

   time(&time_start);

   char outputOrientationCode[4]="PIL";

   DIM output_dim;
   output_dim.nx=output_dim.ny=output_dim.nz=0;
   output_dim.dx=output_dim.dy=output_dim.dz=0.0;

   // filename where the list of images to be registered are input 
   char imagelistfile[DEFAULT_STRING_LENGTH]=""; 

   char outputPrefix[DEFAULT_STRING_LENGTH]="atra";

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) 
      {
         case 'e':
            opt_nn=YES;
            break;
         case 'X':
            output_dim.dx = atof(optarg);
            break;
         case 'Y':
            output_dim.dy = atof(optarg);
            break;
         case 'Z':
            output_dim.dz = atof(optarg);
            break;
         case 'x':
            output_dim.nx = atoi(optarg);
            break;
         case 'y':
            output_dim.ny = atoi(optarg);
            break;
         case 'd':
            del = atoi(optarg);
            break;
         case 'z':
            output_dim.nz = atoi(optarg);
            break;
         case 'm':
            maxiter = atoi(optarg);
            if(maxiter<10 || maxiter>100) maxiter=MAXITER;
            break;
         case 'o':
            strcpy(outputPrefix,optarg);
            break;
         case 'i':
            strcpy(imagelistfile,optarg);
            break;
         case 'v':
            opt_v=YES;
            break;
         case 'N':
            opt_ppm = NO;
            break;
         case 'g':
            opt_png = NO;
            break;
         case 'n':
            opt_txt = NO;
            break;
         case 'O':
            strcpy(outputOrientationCode,optarg);
            break;
         case 'r':
            patch_radius = atoi(optarg);
            if(patch_radius<=0 || patch_radius>MAX_RADIUS) patch_radius=DEFAULT_PATCH_RADIUS;
            break;
         case 'R':
            search_radius = atoi(optarg);
            if(search_radius<=0 || search_radius>MAX_RADIUS) search_radius=DEFAULT_SEARCH_RADIUS;
            break;
         case 't':
            PILcloudthreshold = atoi(optarg);
            if(PILcloudthreshold<0 || PILcloudthreshold>100) PILcloudthreshold=50;
            break;
         case 'h':
            print_help_and_exit();
            break;
         case '?':
            print_help_and_exit();
      }
   }

   for(int c=0; c<4; c++)
      outputOrientationCode[c]=toupper(outputOrientationCode[c]);

   getARTHOME();

   if(argc==1) print_help_and_exit();

   if(imagelistfile[0]=='\0')
   {
      printf("Please specify an image list using the -imlist argument\n");
      exit(0);
   }

   if( isOrientationCodeValid(outputOrientationCode)==0 )
   {
      printf("Sorry, %s is not a valid orientation code.\n",outputOrientationCode);
      exit(0);
   }

   if(opt_v) printf("Post-registration volume orientation: %s\n",outputOrientationCode);

   atra(imagelistfile, output_dim, outputOrientationCode, outputPrefix);

#if 0
   ///////////////////////////////////////////////////////////////////////////////////////////////////
   fp=fopen("secret","w");
   srand48(time(NULL));
   double random_number;
   float random_mat[16];
   random_number=drand48();
   fprintf(fp,"random_number 1 = %lf\n",random_number);
   set_to_I(random_mat,4);
   if(random_number<=(1./6.)) 
   {
      random_mat[3]=0.5;
   }
   else if(random_number<=(2./6.)) 
   {
      random_mat[7]=0.5;
   }
   else if(random_number<=(3./6.)) 
   {
      random_mat[11]=0.5;
   }
   else if(random_number<=(4./6.)) 
   {
      rotate(random_mat, .5*3.14159/180., 1.0, 0.0, 0.0);
   }
   else if(random_number<=(5./6.)) 
   {
      rotate(random_mat, .5*3.14159/180., 0.0, 1.0, 0.0);
   }
   else
   {
      rotate(random_mat, .5*3.14159/180., 0.0, 0.0, 1.0);
   }

   random_number=drand48();
   fprintf(fp,"random_number 2 = %lf\n",random_number);
   if(random_number<=0.25)
   {
      multi(random_mat,4,4,TPIL[2],4,4,TPIL[2]);
      free(PILim[2].v);
      invT = inv4(TPIL[2]);
      PILim[2].v = resliceImage(im[2].v, im[2].nx, im[2].ny, im[2].nz, im[2].dx, im[2].dy, im[2].dz,
      PILim[2].nx, PILim[2].ny, PILim[2].nz, PILim[2].dx, PILim[2].dy, PILim[2].dz, invT, LIN);
      free(invT);
      fprintf(fp,"F0\n");
   }
   else if(random_number<=0.50)
   {
      multi(random_mat,4,4,TPIL[3],4,4,TPIL[3]);
      free(PILim[3].v);
      invT = inv4(TPIL[3]);
      PILim[3].v = resliceImage(im[3].v, im[3].nx, im[3].ny, im[3].nz, im[3].dx, im[3].dy, im[3].dz,
      PILim[3].nx, PILim[3].ny, PILim[3].nz, PILim[3].dx, PILim[3].dy, PILim[3].dz, invT, LIN);
      free(invT);
      fprintf(fp,"F1\n");
   }
   else if(random_number<=0.75)
   {
      multi(random_mat,4,4,TPIL[1],4,4,TPIL[1]);
      free(PILim[1].v);
      invT = inv4(TPIL[1]);
      PILim[1].v = resliceImage(im[1].v, im[1].nx, im[1].ny, im[1].nz, im[1].dx, im[1].dy, im[1].dz,
      PILim[1].nx, PILim[1].ny, PILim[1].nz, PILim[1].dx, PILim[1].dy, PILim[1].dz, invT, LIN);
      free(invT);
      fprintf(fp,"B1\n");
   }
   else
   {
      multi(random_mat,4,4,TPIL[0],4,4,TPIL[0]);
      free(PILim[0].v);
      invT = inv4(TPIL[0]);
      PILim[0].v = resliceImage(im[0].v, im[0].nx, im[0].ny, im[0].nz, im[0].dx, im[0].dy, im[0].dz,
      PILim[0].nx, PILim[0].ny, PILim[0].nz, PILim[0].dx, PILim[0].dy, PILim[0].dz, invT, LIN);
      free(invT);
      fprintf(fp,"B0\n");
   }
   fclose(fp);
   ///////////////////////////////////////////////////////////////////////////////////////////////////
#endif

   time(&time_end);
   if(opt_v) printf("Execution time = %ld minutes and %ld seconds = %ld seconds\n",(time_end-time_start)/60,
   (time_end-time_start)%60, (time_end-time_start));

   return 0;
}
