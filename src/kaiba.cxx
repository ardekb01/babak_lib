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
#include <sph.h>
#include <smooth.h>
#include <landmarks.h>
#include <interpolator.h>
#include <minmax.h>

#define NBIN 256

#define YES 1
#define NO 0
#define MAXNCLASS 15

#ifndef MAXFRAC
#define MXFRAC 0.4
#endif

#ifndef ALPHA_PARAM
#define ALPHA_PARAM 0.0025
#endif

#ifndef BETA_PARAM
#define BETA_PARAM 0.2
// value determined systematically using OASIS_longitudinal data
//#define BETA_PARAM 0.16   
#endif

// Won't be re-defined if this variable is defined at compilation time
#ifndef MAXITER
#define MAXITER 20
#endif

// PIL brain cloud threshold for definining a brain mask
#ifndef CLOUD_THRESH
#define CLOUD_THRESH 90 
#endif

#ifndef TOLERANCE
#define TOLERANCE 1.e-7
#endif

int opt;

/////////////////////////////////////////////////////////////////////////
// Global variables

char opt_flip=YES;
char flipped; // Takes YES or NO
float alpha_param;

/////////////////////////////////////////////////////////////////////////

static struct option options[] =
{
   {"-v",0,'v'},
   {"-F",0,'F'},
   {"-noflip",0,'F'},
   {"-version",0,'V'},
   {"-V",0,'V'},
   {"-nopng",0,'g'},
   {"-p",1,'p'},   // output prefix
   {"-o",1,'p'},   // output prefix
   {"-b",1,'b'},   // baseline image
   {"-i",1,'b'},   // baseline image
   {"-lm",1,'l'},  // landmark  file
   {"-alpha",1,'a'}, 
   {"-a",1,'a'}, 
   {"-noppm",0,'N'}, 
   {"-notxt",0,'t'}, 
   {0,0,0}
};

void print_help_and_exit()
{
   printf("\nUsage: kaiba [options] -o <output prefix> -i <T1W NIFTI>.nii or -i <image list>.txt\n"
   "\nRequired arguments:\n"
   "   -o <prefix>: Prefix used for naming output files\n"
   "   -i <T1W NIFTI>.nii: T1-weighted 3D structural MRI (must be NIFTI format of type short)\n"
   "   -i <image list>.txt: Image list output from ATRA program\n"
   "\nOptions:\n"
   "   -v : Enables verbose mode\n"
   "   -V or -version : Prints program version\n"
   "   -notxt: Prevents outputtign *ACPC.txt\n"
   "   -noppm : Prevents outputtign *.ppm and *.png images\n"
   "   -nopng : Prevents output images in PNG format (still outputs PPM)\n"
   "   -lm <filename>: Manually specifies AC/PC/RP landmarks for <T1W NIFIT>.nii\n"
//   "   -a or -alpha <alpha>: Specifies to alpha parameter\n"
   "\n");

   exit(0);
}

// this is worked out in my technical note notebook
// compute 1/sqrt(2*var) * (integral from -inf to x) of exp[ (x-mu)^2/(2*var) ]
float8 normalCDF(float8 x, float8 mu, float8 var)
{
   if( var == 0.0 ) 
   {
      printf("Warning: Zero variance passed to normalCDF()\n");
      return(0.0);
   }

   return( 0.5 + 0.5*erf ( (x-mu)/sqrt(2.0*var) ) );
}


//////////////////////////////////////////////////////////////////////////////////////////////////
float8 ssd_cost_function(float4 *T, DIM dimb, DIM dimf, float4 *sclbim, float4 *sclfim, int2 *bmsk, int2 *fmsk)
{
   int kmin_f=0; 
   int kmax_f=dimf.nz-1;
   int jmin_f=0; 
   int jmax_f=dimf.ny-1;
   int imin_f=0; 
   int imax_f=dimf.nx-1;
   int kmin_b=0;
   int kmax_b=dimb.nz-1;
   int jmin_b=0;
   int jmax_b=dimb.ny-1;
   int imin_b=0;
   int imax_b=dimb.nx-1;

   float4 *invT;
   float4 dif;
   float8 cost=0.0;
   int v, slice_offset, offset;
   float4 Tmod[16]; //modified T
   float4 invTmod[16]; //modified invT

   float4 psub0, psub1, psub2;
   float4 nxsub2, nysub2, nzsub2; 

   float4 ptrg0, ptrg1, ptrg2;
   float4 nxtrg2, nytrg2, nztrg2;

   invT = inv4(T);

   nxsub2 = (dimf.nx-1)/2.0;
   nysub2 = (dimf.ny-1)/2.0;
   nzsub2 = (dimf.nz-1)/2.0;

   nxtrg2 = (dimb.nx-1)/2.0;
   nytrg2 = (dimb.ny-1)/2.0;
   nztrg2 = (dimb.nz-1)/2.0;

   ////////////////////////////////////////
   Tmod[0] = T[0]*dimf.dx/dimb.dx;
   Tmod[1] = T[1]*dimf.dy/dimb.dx;
   Tmod[2] = T[2]*dimf.dz/dimb.dx;
   Tmod[3] = T[3]/dimb.dx + nxtrg2;

   Tmod[4] = T[4]*dimf.dx/dimb.dy;
   Tmod[5] = T[5]*dimf.dy/dimb.dy;
   Tmod[6] = T[6]*dimf.dz/dimb.dy;
   Tmod[7] = T[7]/dimb.dy + nytrg2;

   Tmod[8] = T[8]*dimf.dx/dimb.dz;
   Tmod[9] = T[9]*dimf.dy/dimb.dz;
   Tmod[10] = T[10]*dimf.dz/dimb.dz;
   Tmod[11] = T[11]/dimb.dz + nztrg2;
   ////////////////////////////////////////
   invTmod[0] = invT[0]*dimb.dx/dimf.dx;
   invTmod[1] = invT[1]*dimb.dy/dimf.dx;
   invTmod[2] = invT[2]*dimb.dz/dimf.dx;
   invTmod[3] = invT[3]/dimf.dx + nxsub2;

   invTmod[4] = invT[4]*dimb.dx/dimf.dy;
   invTmod[5] = invT[5]*dimb.dy/dimf.dy;
   invTmod[6] = invT[6]*dimb.dz/dimf.dy;
   invTmod[7] = invT[7]/dimf.dy + nysub2;

   invTmod[8] = invT[8]*dimb.dx/dimf.dz;
   invTmod[9] = invT[9]*dimb.dy/dimf.dz;
   invTmod[10] = invT[10]*dimb.dz/dimf.dz;
   invTmod[11] = invT[11]/dimf.dz + nzsub2;
   ////////////////////////////////////////
   
   float4 t2, t6, t10;
   float4 t1, t5, t9;

   for(int k=kmin_f; k<=kmax_f; k++)
   {
      slice_offset = k*dimf.np;
      psub2 = (k-nzsub2);
      t2  = Tmod[2]*psub2  + Tmod[3];
      t6  = Tmod[6]*psub2  + Tmod[7];
      t10 = Tmod[10]*psub2 + Tmod[11];
      for(int j=jmin_f; j<=jmax_f; j++)
      {
         offset = slice_offset + j*dimf.nx;
         psub1 = (j-nysub2);
         t1 = Tmod[1]*psub1;
         t5 = Tmod[5]*psub1;
         t9 = Tmod[9]*psub1;
         for(int i=imin_f; i<=imax_f; i++)
         {
            v = offset + i;

            if( fmsk[v]>0)
            {
               psub0 = (i-nxsub2);

               ptrg0 = Tmod[0]*psub0 + t1 + t2;
               ptrg1 = Tmod[4]*psub0 + t5 + t6;
               ptrg2 = Tmod[8]*psub0 + t9 + t10;

               dif = sclfim[v] - linearInterpolator(ptrg0, ptrg1, ptrg2, sclbim, dimb.nx, dimb.ny, dimb.nz, dimb.np);
               cost += (dif*dif);
            }
         }
      }
   }

   for(int k=kmin_b; k<=kmax_b; k++)
   {
      slice_offset = k*dimb.np;
      ptrg2 = (k-nztrg2);
      t2  = invTmod[2]*ptrg2  + invTmod[3];
      t6  = invTmod[6]*ptrg2  + invTmod[7];
      t10 = invTmod[10]*ptrg2 + invTmod[11];
      for(int j=jmin_b; j<=jmax_b; j++)
      {
         offset = slice_offset + j*dimb.nx;
         ptrg1 = (j-nytrg2);
         t1 = invTmod[1]*ptrg1;
         t5 = invTmod[5]*ptrg1;
         t9 = invTmod[9]*ptrg1;
         for(int i=imin_b; i<=imax_b; i++)
         {
            v = offset + i;

            if( bmsk[v]>0)
            {
               ptrg0 = (i-nxtrg2);

               psub0 = invTmod[0]*ptrg0 + t1 + t2;
               psub1 = invTmod[4]*ptrg0 + t5 + t6;
               psub2 = invTmod[8]*ptrg0 + t9 + t10;

               dif = sclbim[v] - linearInterpolator(psub0, psub1, psub2, sclfim, dimf.nx, dimf.ny, dimf.nz, dimf.np);
               cost += (dif*dif);
            }
         }
      }
   }

   free(invT);
   return(cost);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

float8 ncc_cost_function(float4 *T, DIM dimb, DIM dimf, float4 *sclbim, float4 *sclfim, int2 *bmsk, int2 *fmsk)
{
   int kmin_f=0; 
   int kmax_f=dimf.nz-1;
   int jmin_f=0; 
   int jmax_f=dimf.ny-1;
   int imin_f=0; 
   int imax_f=dimf.nx-1;
   int kmin_b=0;
   int kmax_b=dimb.nz-1;
   int jmin_b=0;
   int jmax_b=dimb.ny-1;
   int imin_b=0;
   int imax_b=dimb.nx-1;

   int n=0; // number of voxels that take part in the NCC calculation
   float8 sum1, sum2, sum11, sum22, sum12;
   float4 subject_image_value, target_image_value;

   float4 *invT;
   float8 cost=0.0;
   int v, slice_offset, offset;
   float4 Tmod[16]; //modified T
   float4 invTmod[16]; //modified invT

   float4 psub0, psub1, psub2;
   float4 nxsub2, nysub2, nzsub2; 

   float4 ptrg0, ptrg1, ptrg2;
   float4 nxtrg2, nytrg2, nztrg2;

   invT = inv4(T);

   nxsub2 = (dimf.nx-1)/2.0;
   nysub2 = (dimf.ny-1)/2.0;
   nzsub2 = (dimf.nz-1)/2.0;

   nxtrg2 = (dimb.nx-1)/2.0;
   nytrg2 = (dimb.ny-1)/2.0;
   nztrg2 = (dimb.nz-1)/2.0;

   ////////////////////////////////////////
   Tmod[0] = T[0]*dimf.dx/dimb.dx;
   Tmod[1] = T[1]*dimf.dy/dimb.dx;
   Tmod[2] = T[2]*dimf.dz/dimb.dx;
   Tmod[3] = T[3]/dimb.dx + nxtrg2;

   Tmod[4] = T[4]*dimf.dx/dimb.dy;
   Tmod[5] = T[5]*dimf.dy/dimb.dy;
   Tmod[6] = T[6]*dimf.dz/dimb.dy;
   Tmod[7] = T[7]/dimb.dy + nytrg2;

   Tmod[8] = T[8]*dimf.dx/dimb.dz;
   Tmod[9] = T[9]*dimf.dy/dimb.dz;
   Tmod[10] = T[10]*dimf.dz/dimb.dz;
   Tmod[11] = T[11]/dimb.dz + nztrg2;
   ////////////////////////////////////////
   invTmod[0] = invT[0]*dimb.dx/dimf.dx;
   invTmod[1] = invT[1]*dimb.dy/dimf.dx;
   invTmod[2] = invT[2]*dimb.dz/dimf.dx;
   invTmod[3] = invT[3]/dimf.dx + nxsub2;

   invTmod[4] = invT[4]*dimb.dx/dimf.dy;
   invTmod[5] = invT[5]*dimb.dy/dimf.dy;
   invTmod[6] = invT[6]*dimb.dz/dimf.dy;
   invTmod[7] = invT[7]/dimf.dy + nysub2;

   invTmod[8] = invT[8]*dimb.dx/dimf.dz;
   invTmod[9] = invT[9]*dimb.dy/dimf.dz;
   invTmod[10] = invT[10]*dimb.dz/dimf.dz;
   invTmod[11] = invT[11]/dimf.dz + nzsub2;
   ////////////////////////////////////////
   
   float4 t2, t6, t10;
   float4 t1, t5, t9;

   // initialize sums to zero
   sum1=sum2=sum11=sum22=sum12=0.0;

   for(int k=kmin_f; k<=kmax_f; k++)
   {
      slice_offset = k*dimf.np;
      psub2 = (k-nzsub2);
      t2  = Tmod[2]*psub2  + Tmod[3];
      t6  = Tmod[6]*psub2  + Tmod[7];
      t10 = Tmod[10]*psub2 + Tmod[11];
      for(int j=jmin_f; j<=jmax_f; j++)
      {
         offset = slice_offset + j*dimf.nx;
         psub1 = (j-nysub2);
         t1 = Tmod[1]*psub1;
         t5 = Tmod[5]*psub1;
         t9 = Tmod[9]*psub1;
         for(int i=imin_f; i<=imax_f; i++)
         {
            v = offset + i;

            if( fmsk[v]>0)
            {
               n++;

               psub0 = (i-nxsub2);

               ptrg0 = Tmod[0]*psub0 + t1 + t2;
               ptrg1 = Tmod[4]*psub0 + t5 + t6;
               ptrg2 = Tmod[8]*psub0 + t9 + t10;

               subject_image_value = sclfim[v];
               target_image_value = linearInterpolator(ptrg0, ptrg1, ptrg2, sclbim, dimb.nx, dimb.ny, dimb.nz, dimb.np);
               sum1 += subject_image_value;
               sum2 += target_image_value;
               sum12 += (subject_image_value*target_image_value);
               sum11 += (subject_image_value*subject_image_value);
               sum22 += (target_image_value*target_image_value);
            }
         }
      }
   }

   for(int k=kmin_b; k<=kmax_b; k++)
   {
      slice_offset = k*dimb.np;
      ptrg2 = (k-nztrg2);
      t2  = invTmod[2]*ptrg2  + invTmod[3];
      t6  = invTmod[6]*ptrg2  + invTmod[7];
      t10 = invTmod[10]*ptrg2 + invTmod[11];
      for(int j=jmin_b; j<=jmax_b; j++)
      {
         offset = slice_offset + j*dimb.nx;
         ptrg1 = (j-nytrg2);
         t1 = invTmod[1]*ptrg1;
         t5 = invTmod[5]*ptrg1;
         t9 = invTmod[9]*ptrg1;
         for(int i=imin_b; i<=imax_b; i++)
         {
            v = offset + i;

            if( bmsk[v]>0)
            {
               n++;

               ptrg0 = (i-nxtrg2);

               psub0 = invTmod[0]*ptrg0 + t1 + t2;
               psub1 = invTmod[4]*ptrg0 + t5 + t6;
               psub2 = invTmod[8]*ptrg0 + t9 + t10;

               subject_image_value = linearInterpolator(psub0, psub1, psub2, sclfim, dimf.nx, dimf.ny, dimf.nz, dimf.np);
               target_image_value = sclbim[v];

               sum1 += subject_image_value;
               sum2 += target_image_value;
               sum12 += (subject_image_value*target_image_value);
               sum11 += (subject_image_value*subject_image_value);
               sum22 += (target_image_value*target_image_value);
            }
         }
      }
   }

   free(invT);

   if( n > 0 )
   {
//      cost = (sum12 - sum1*sum2/n);
      cost = -(sum12 - sum1*sum2/n);  // since it is a minimization problem
      cost /= sqrt( sum11 - sum1*sum1/n ); 
      cost /= sqrt( sum22 - sum2*sum2/n ); 
   }

   return(cost);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

// This function computes and returns a 4x4 transformation matrix T.
// Applying this transformation to a point p=(x',y',z') in the image
// coordinates system (ICS) yields the coordinates (x,y,z) of the same point with
// respect to the magnet coordinates system.  That is: (x,y,z)=T(x',y',z').
// In ART, the origin of the ICS is the center of the image
// volume.  The x-axis in the ICS is pointing to the right of the image
// (which is not necessarily the right of the subject!).  The y-axis in the
// ICS is pointing down.  The z-axis is determined by the right-hand-rule.
void transformation_to_magnet_coordinates(nifti_1_header hdr, float4 *T)
{
   float4 rowvec[3]={0.,0.,0.};
   float4 columnvec[3]={0.,0.,0.};
   float4 normalvec[3]={0.,0.,0.};
   float4 centervec[3]={0.,0.,0.};

   //printf("sizeof_hdr = %d\n", hdr.sizeof_hdr);
   //printf("number of dimensions = %d\n", hdr.dim[0]);
   //printf("matrix size = %d x %d x %d\n", hdr.dim[1], hdr.dim[2], hdr.dim[3]);
   //printf("voxel size = %8.6f x %8.6f x %8.6f\n", hdr.pixdim[1], hdr.pixdim[2], hdr.pixdim[3]);
   //printf("datatype = %d\n",hdr.datatype);
   //printf("vox_offset = %d\n", (int)hdr.vox_offset);
   //printf("qform_code = %d\n", hdr.qform_code);
   //printf("sform_code = %d\n", hdr.sform_code);
   //printf("magic code = %s\n", hdr.magic);
   //printf("srow_x: %f %f %f %f\n", hdr.srow_x[0],hdr.srow_x[1],hdr.srow_x[2],hdr.srow_x[3]);
   //printf("srow_y: %f %f %f %f\n", hdr.srow_y[0],hdr.srow_y[1],hdr.srow_y[2],hdr.srow_y[3]);
   //printf("srow_z: %f %f %f %f\n", hdr.srow_z[0],hdr.srow_z[1],hdr.srow_z[2],hdr.srow_z[3]);
   //printf("quatern_b = %f\n", hdr.quatern_b);
   //printf("quatern_c = %f\n", hdr.quatern_c);
   //printf("quatern_d = %f\n", hdr.quatern_d);
   //printf("qoffset_x = %f\n",hdr.qoffset_x);
   //printf("qoffset_y = %f\n",hdr.qoffset_y);
   //printf("qoffset_z = %f\n",hdr.qoffset_z);

   int nx,ny,nz;
   float4 dx,dy,dz;

   nx = hdr.dim[1];
   ny = hdr.dim[2];
   nz = hdr.dim[3];

   dx=hdr.pixdim[1];
   dy=hdr.pixdim[2];
   dz=hdr.pixdim[3];

   if(hdr.qform_code>0 )
   {
      float4 dum;
      mat44 R;
   
      R=nifti_quatern_to_mat44( hdr.quatern_b, hdr.quatern_c, hdr.quatern_d,
      hdr.qoffset_x, hdr.qoffset_y, hdr.qoffset_z, dx, dy, dz, hdr.pixdim[0]);

      dum = R.m[0][0]*R.m[0][0] + R.m[1][0]*R.m[1][0] + R.m[2][0]*R.m[2][0];
      dum = sqrtf(dum);

      if(dum != 0.0)
      {
         // note that the -tive signs converts from NIFTI's RAS to ART's LAI
         rowvec[0]=-R.m[0][0]/dum;
         rowvec[1]=R.m[1][0]/dum;
         rowvec[2]=-R.m[2][0]/dum;
      }
      else
      {
         rowvec[0]=1.0;
         rowvec[1]=0.0;
         rowvec[2]=0.0;
      }

      dum = R.m[0][1]*R.m[0][1] + R.m[1][1]*R.m[1][1] + R.m[2][1]*R.m[2][1];
      dum = sqrtf(dum);

      if(dum != 0.0)
      {
         // note that the -tive signs converts from NIFTI's RAS to ART's LAI
         columnvec[0]=-R.m[0][1]/dum;
         columnvec[1]=R.m[1][1]/dum;
         columnvec[2]=-R.m[2][1]/dum;
      }
      else
      {
         columnvec[0]=0.0;
         columnvec[1]=1.0;
         columnvec[2]=0.0;
      }

      dum = R.m[0][2]*R.m[0][2] + R.m[1][2]*R.m[1][2] + R.m[2][2]*R.m[2][2];
      dum = sqrtf(dum);

      if(dum != 0.0)
      {
         // note that the -tive signs converts from NIFTI's RAS to ART's LAI
         normalvec[0]=-R.m[0][2]/dum;
         normalvec[1]=R.m[1][2]/dum;
         normalvec[2]=-R.m[2][2]/dum;
      }
      else
      {
         normalvec[0]=0.0;
         normalvec[1]=0.0;
         normalvec[2]=1.0;
      }

      // note that the -tive signs converts from NIFTI's RAS to ART's LAI
      centervec[0] = -R.m[0][3] +
      rowvec[0] * dx*(nx-1.0)/2.0 +
      columnvec[0] * dy*(ny-1.0)/2.0;

      centervec[1] = R.m[1][3] +
      rowvec[1] * dx*(nx-1.0)/2.0 +
      columnvec[1] * dy*(ny-1.0)/2.0;

      centervec[2] = -R.m[2][3] +
      rowvec[2] * dx*(nx-1.0)/2.0 +
      columnvec[2] * dy*(ny-1.0)/2.0; 
   }
   else if(hdr.sform_code>0)
   {
      float4 dum;

      dum = hdr.srow_x[0]*hdr.srow_x[0] + hdr.srow_y[0]*hdr.srow_y[0] + hdr.srow_z[0]*hdr.srow_z[0];
      dum = sqrtf(dum);

      if(dum != 0.0)
      {
         // note that the -tive signs converts from NIFTI's RAS to ART's LAI
         rowvec[0]=-hdr.srow_x[0]/dum;
         rowvec[1]=hdr.srow_y[0]/dum;
         rowvec[2]=-hdr.srow_z[0]/dum;
      }
      else
      {
         rowvec[0]=1.0;
         rowvec[1]=0.0;
         rowvec[2]=0.0;
      }

      dum= hdr.srow_x[1]*hdr.srow_x[1] + hdr.srow_y[1]*hdr.srow_y[1] + hdr.srow_z[1]*hdr.srow_z[1];
      dum = sqrtf(dum);

      if(dum != 0.0)
      {
         // note that the -tive signs converts from NIFTI's RAS to ART's LAI
         columnvec[0]=-hdr.srow_x[1]/dum;
         columnvec[1]=hdr.srow_y[1]/dum;
         columnvec[2]=-hdr.srow_z[1]/dum;
      }
      else
      {
         columnvec[0]=0.0;
         columnvec[1]=1.0;
         columnvec[2]=0.0;
      }

      dum= hdr.srow_x[2]*hdr.srow_x[2] + hdr.srow_y[2]*hdr.srow_y[2] + hdr.srow_z[2]*hdr.srow_z[2];
      dum = sqrtf(dum);

      if(dum != 0.0)
      {
         // note that the -tive signs converts from NIFTI's RAS to ART's LAI
         normalvec[0]=-hdr.srow_x[2]/dum;
         normalvec[1]=hdr.srow_y[2]/dum;
         normalvec[2]=-hdr.srow_z[2]/dum;
      }
      else
      {
         normalvec[0]=0.0;
         normalvec[1]=0.0;
         normalvec[2]=1.0;
      }

      // note that the -tive signs converts from NIFTI's RAS to ART's LAI
      centervec[0] = -hdr.srow_x[3] +
      rowvec[0] * dx*(nx-1.0)/2.0 +
      columnvec[0] * dy*(ny-1.0)/2.0;

      centervec[1] = hdr.srow_y[3] +
      rowvec[1] * dx*(nx-1.0)/2.0 +
      columnvec[1] * dy*(ny-1.0)/2.0;

      centervec[2] = -hdr.srow_z[3] +
      rowvec[2] * dx*(nx-1.0)/2.0 +
      columnvec[2] * dy*(ny-1.0)/2.0; 
   }
   else
   {
      printf("\n**WARNING**: NIFTI header did not contain image orientation information.\n\n");
   }

   //printf("Center Vector = (%7.5lf %7.5lf %7.5lf)\n", centervec[0],centervec[1], centervec[2]);
   //printf("Normal Vector = (%7.5lf %7.5lf %7.5lf)\n", normalvec[0],normalvec[1], normalvec[2]);
   //printf("Row Vector = (%7.5lf %7.5lf %7.5lf)\n", rowvec[0],rowvec[1], rowvec[2]);
   //printf("Column Vector = (%7.5lf %7.5lf %7.5lf)\n", columnvec[0],columnvec[1], columnvec[2]);

   for(int i=0; i<16; i++) T[i]=0.0;

   // The fifteenth element (i.e., element 4,4) is always 1.0.
   T[15]=1.0;
   
   T[0] = rowvec[0];
   T[4] = rowvec[1];
   T[8] = rowvec[2];

   T[1] = columnvec[0];
   T[5] = columnvec[1];
   T[9] = columnvec[2];

   T[2] = normalvec[0];
   T[6] = normalvec[1];
   T[10]= normalvec[2];

   T[3] = centervec[0] + (dz*(nz-1.0)/2.0)*normalvec[0];
   T[7] = centervec[1] + (dz*(nz-1.0)/2.0)*normalvec[1];
   T[11]= centervec[2] + (dz*(nz-1.0)/2.0)*normalvec[2];

   return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////

void atra_to_fsl(float4 *Matra, float4 *Mfsl, DIM dimf, DIM trg_dim)
{
   float4 Tsub[16], Ttrg[16];
   float4 *inv_Tsub;

   Tsub[0]=1.0;  Tsub[1]=0.0;  Tsub[2]=0.0;  Tsub[3]=(dimf.nx-1.0)*dimf.dx/2.0;
   Tsub[4]=0.0;  Tsub[5]=1.0;  Tsub[6]=0.0;  Tsub[7]=(dimf.ny-1.0)*dimf.dy/2.0;
   Tsub[8]=0.0;  Tsub[9]=0.0;  Tsub[10]=1.0; Tsub[11]=(dimf.nz-1.0)*dimf.dz/2.0;
   Tsub[12]=0.0; Tsub[13]=0.0; Tsub[14]=0.0; Tsub[15]=1.0;

   Ttrg[0]=1.0;  Ttrg[1]=0.0;  Ttrg[2]=0.0;  Ttrg[3]=(trg_dim.nx-1.0)*trg_dim.dx/2.0;
   Ttrg[4]=0.0;  Ttrg[5]=1.0;  Ttrg[6]=0.0;  Ttrg[7]=(trg_dim.ny-1.0)*trg_dim.dy/2.0;
   Ttrg[8]=0.0;  Ttrg[9]=0.0;  Ttrg[10]=1.0; Ttrg[11]=(trg_dim.nz-1.0)*trg_dim.dz/2.0;
   Ttrg[12]=0.0; Ttrg[13]=0.0; Ttrg[14]=0.0; Ttrg[15]=1.0;

   inv_Tsub = inv4(Tsub);

   multi(Ttrg,4,4,Matra,4,4,Mfsl);
   multi(Mfsl,4,4,inv_Tsub,4,4, Mfsl);

   free(inv_Tsub);
}

/////////////////////////////////////////////////
// find sqrt(T) and inverse_sqrt(T)
/////////////////////////////////////////////////
void sqrt_matrix(float4 *T, float4 *sqrtT, float4 *invsqrtT)
{
   float4 w[3], v[3];
   float4 theta;

   SE3_to_se3(T, w, v, theta);

   theta /= 2.0;
   se3_to_SE3(sqrtT, w, v, theta);

   theta *= -1.0;
   se3_to_SE3(invsqrtT, w, v, theta);

   return;
}

// bfile: baseline image filename
// ffile: follow-up image filename
void symmetric_registration(SHORTIM &aimpil, const char *bfile, const char *ffile, const char *blmfile,const char *flmfile, int verbose)
{
   char orient[4]="";
   int2 *fmsk, *bmsk;
   float4 *sclfim, *sclbim;
   int2 *PILbraincloud;
   DIM PILbraincloud_dim;
   nifti_1_header PILbraincloud_hdr; 

   char filename[DEFAULT_STRING_LENGTH]="";  // a generic filename for reading/writing stuff

   DIM dimf; // follow-up image dimensions structure
   DIM dimb; // baseline image dimensions structure
   char bprefix[DEFAULT_STRING_LENGTH]=""; //baseline image prefix
   char fprefix[DEFAULT_STRING_LENGTH]=""; //follow-up image prefix
   float4 T[16]; // The unknown transformation matrix that takes points from the follow-up to baseline space 
   float4 Tf[16]; // The unknown transformation matrix that takes points from the follow-up to mid PIL space 
   float4 Tb[16]; // The unknown transformation matrix that takes points from the baseline to mid PIL space 
   float4 Tinter[16]; // Transforms points from the follow-up PIL to baseline PIL spaces
   float4 *invT;  // inverse of T
   float4 sqrtTinter[16];
   float4 invsqrtTinter[16];

   /////////////////////////////////////////////////////////////////////////////////////////////
   // read PILbraincloud.nii from the $ARTHOME directory
   /////////////////////////////////////////////////////////////////////////////////////////////
   snprintf(filename,sizeof(filename),"%s/PILbrain.nii",ARTHOME);

   PILbraincloud = (int2 *)read_nifti_image(filename, &PILbraincloud_hdr);

   if(PILbraincloud==NULL)
   {
      printf("Error reading %s, aborting ...\n", filename);
      exit(1);
   }

   if(verbose)
   {
      printf("PIL brain cloud: %s\n",filename);
      printf("PIL brain cloud threshold level = %d%%\n",CLOUD_THRESH);
      printf("PIL brain cloud matrix size = %d x %d x %d (voxels)\n", 
      PILbraincloud_hdr.dim[1], PILbraincloud_hdr.dim[2], PILbraincloud_hdr.dim[3]);
      printf("PIL brain cloud voxel size = %8.6f x %8.6f x %8.6f (mm^3)\n", 
      PILbraincloud_hdr.pixdim[1], PILbraincloud_hdr.pixdim[2], PILbraincloud_hdr.pixdim[3]);
   }

   set_dim(PILbraincloud_dim, PILbraincloud_hdr);
   /////////////////////////////////////////////////////////////////////////////////////////////

   if(verbose)
   {
      printf("Starting unbiased symmetric registration ...\n");
      printf("ARTHOME: %s\n",ARTHOME);
      printf("Maximum number of iterations = %d\n", MAXITER);
      printf("Baseline image: %s\n",bfile);
      printf("Follow-up image: %s\n",ffile);
   }

   // Note: niftiFilename does a few extra checks to ensure that the file has either
   // .hdr or .nii extension, the magic field in the header is set correctly, 
   // the file can be opened and a header can be read.
   if( niftiFilename(bprefix, bfile)==0 )
   {
      exit(0);
   }

   if(verbose)
   {
      printf("Baseline image prefix: %s\n",bprefix);
   }

   if( niftiFilename(fprefix, ffile)==0 )
   {
      exit(0);
   }

   if(verbose)
   {
      printf("Follow-up image prefix: %s\n",fprefix);
   }
   //////////////////////////////////////////////////////////////////////////////////

   float4 bTPIL[16]; // takes the baseline image to standard PIL orientation 
   float4 fTPIL[16]; // takes the follow-up image to standard PIL orientation 
   float4 *ibTPIL; // inverse of bTPIL
   float4 *ifTPIL; // inverse of fTPIL

   if(verbose) printf("Computing baseline image PIL transformation ...\n");
   orient[0]='\0';
   new_PIL_transform(bfile, blmfile, orient, bTPIL, 0);

   ibTPIL= inv4(bTPIL);

   if(verbose) printf("Computing follow-up image PIL transformation ...\n");
   orient[0]='\0';
   new_PIL_transform(ffile, flmfile, orient, fTPIL, 0);

   ifTPIL= inv4(fTPIL);

   ///////////////////////////////////////////////////////////////////////////////////////////////
   // Read baseline and follow-up images
   ///////////////////////////////////////////////////////////////////////////////////////////////
   int2 *bim; // baseline image
   int2 *fim; // follow-up image
   nifti_1_header bhdr;  // baseline image NIFTI header
   nifti_1_header fhdr;  // follow-up image NIFTI header

   bim = (int2 *)read_nifti_image(bfile, &bhdr);

   if(bim==NULL)
   {
      printf("Error reading %s, aborting ...\n", bfile);
      exit(1);
   }

   set_dim(dimb, bhdr);

   fim = (int2 *)read_nifti_image(ffile, &fhdr);

   if(fim==NULL)
   {
         printf("Error reading %s, aborting ...\n", ffile);
         exit(1);
   }

   set_dim(dimf, fhdr);
   ///////////////////////////////////////////////////////////////////////////////////////////////

   ///////////////////////////////////////////////////////////////////////////////////////////////
   // determine subject and target masks
   ///////////////////////////////////////////////////////////////////////////////////////////////
   {
      float4 Tdum[16];

      for(int i=0; i<16; i++) Tdum[i]=bTPIL[i];
      bmsk = resliceImage(PILbraincloud, PILbraincloud_dim, dimb, Tdum, LIN);
      for(int v=0; v<dimb.nv; v++) if(bmsk[v]<CLOUD_THRESH) bmsk[v]=0;
      //save_nifti_image("bmsk.nii", bmsk, &bhdr);

      for(int i=0; i<16; i++) Tdum[i]=fTPIL[i];
      fmsk = resliceImage(PILbraincloud, PILbraincloud_dim, dimf, Tdum, LIN);
      for(int v=0; v<dimf.nv; v++) if(fmsk[v]<CLOUD_THRESH) fmsk[v]=0;
      //save_nifti_image("fmsk.nii", fmsk, &fhdr);
      
      delete PILbraincloud;
   }
   ///////////////////////////////////////////////////////////////////////////////////////////////
   
   ///////////////////////////////////////////////////////////////////////////////////////////////
   {
      float4 bscale;
      float4 fscale;

      trimExtremes(bim, bmsk, dimb.nv, 0.05);
      trimExtremes(fim, fmsk, dimf.nv, 0.05);

      bscale=imageMean(bim, bmsk, dimb.nv);
      fscale=imageMean(fim, fmsk, dimf.nv);

      sclfim = (float4 *)calloc(dimf.nv, sizeof(float4));
      sclbim = (float4 *)calloc(dimb.nv, sizeof(float4));

      for(int v=0; v<dimf.nv; v++) sclfim[v] = fim[v]/fscale;
      for(int v=0; v<dimb.nv; v++) sclbim[v] = bim[v]/bscale;
   }

   {
      float8 relative_change;
      float8 mincost, oldmincost, cost;
      float8 (*cost_function)(float4 *T, DIM dimb, DIM dimf, float4 *sclbim, float4 *sclfim,int2 *bmsk, int2 *fmsk);
      float4 P[6];
      float4 Pmin[6];
      float4 stepsize[6]={0.25, 0.25, 0.25, 0.1, 0.1, 0.1};  // stepsize used in optimization
      //float4 iP[6]={3.0, 3.0, 3.0, 1.5, 1.5, 1.5}; // interval used in optimization
      // New interval makes it twice as fast with same resutls
      float4 iP[6]={1.0, 1.0, 1.0, 1.0, 1.0, 1.0}; // interval used in optimization 

      cost_function=ssd_cost_function; // used for T1 to T1 registration
//      cost_function=ncc_cost_function; 
      set_to_I(Tinter,4);

      // initially assume Tinter=Identity matrix
      multi(ibTPIL, 4, 4,  fTPIL, 4,  4, T);
      oldmincost = mincost = cost_function(T, dimb, dimf, sclbim, sclfim,bmsk,fmsk);

      for(int j=0; j<6; j++) P[j]=Pmin[j]=0.0;

      // automatically sets the step size for the first three variables
      //if( dimf.dx<dimb.dx) stepsize[0]=dimf.dx/5.0;
      //else stepsize[0]=dimb.dx/5.0;
      //if( dimf.dy<dimb.dy) stepsize[1]=dimf.dy/5.0;
      //else stepsize[1]=dimb.dy/5.0;
      //if( dimf.dz<dimb.dz) stepsize[2]=dimf.dz/5.0;
      //else stepsize[2]=dimb.dz/5.0;
      //printf("Step sizes = %f %f %f %f %f %f\n",stepsize[0],stepsize[1],stepsize[2],stepsize[3],stepsize[4],stepsize[5]); 

      if(verbose)
      {
         printf("Tolerance = %3.1e\n",TOLERANCE);
         printf("Initial cost = %f\n", mincost);
      }

      for(int iter=1; iter<=MAXITER; iter++)
      {
         if(verbose)
         {
            printf("Iteration %d ...\n",iter);
         }

         for(int i=0; i<6; i++)
         {
            for(P[i] = Pmin[i]-iP[i]; P[i]<=Pmin[i]+iP[i]; P[i]+=stepsize[i] )
            {
               set_transformation(P[0], P[1], P[2], P[3], P[4], P[5], "ZXYT", Tinter);
               multi(Tinter, 4, 4,  fTPIL, 4,  4, T);
               multi(ibTPIL, 4, 4,  T, 4,  4, T);
   
               cost = cost_function(T, dimb, dimf, sclbim, sclfim,bmsk,fmsk);
   
               if( cost < mincost )
               {
                  Pmin[i]=P[i];
                  mincost = cost;
               }
            }
            P[i]=Pmin[i];

            if(verbose)
            {
               printf("P0=%f P1=%f P2=%f P3=%f P4=%f P5=%f\n", Pmin[0], Pmin[1], Pmin[2], Pmin[3], Pmin[4], Pmin[5]);
            }
         }

         if(oldmincost != 0.0)
         {
           relative_change=(oldmincost-mincost)/fabs(oldmincost);
         }
         else
         {
           relative_change=0.0;
         }

         if(verbose)
         {
            printf("Cost = %f\n", mincost);
            printf("Relative change = %3.1e x 100%%\n", relative_change );
         }
   
         if( oldmincost==0.0 || relative_change <= TOLERANCE )
            break;
         else
            oldmincost = mincost;
      }

      set_transformation(P[0], P[1], P[2], P[3], P[4], P[5], "ZXYT", Tinter);

      if( Tinter[0]!=1.0 || Tinter[1]!=0.0 || Tinter[2]!=0.0 || Tinter[3]!=0.0 ||
      Tinter[4]!=0.0 || Tinter[5]!=1.0 || Tinter[6]!=0.0 || Tinter[7]!=0.0 ||
      Tinter[8]!=0.0 || Tinter[9]!=0.0 || Tinter[10]!=1.0 || Tinter[11]!=0.0 ||
      Tinter[12]!=0.0 || Tinter[13]!=0.0 || Tinter[14]!=0.0 || Tinter[15]!=1.0)
      {
         // Tinter does not equal identity matrix
         sqrt_matrix(Tinter, sqrtTinter, invsqrtTinter);
      }
      {
         // Tinter equals identity matrix
         for(int i=0; i<16; i++) sqrtTinter[i]=invsqrtTinter[i]=Tinter[i];
      }
      multi(sqrtTinter, 4, 4,  fTPIL, 4,  4, Tf);
      multi(invsqrtTinter, 4, 4,  bTPIL, 4,  4, Tb);
   }
   
   /////////////////////////////////////////////////
   // save transformation matrices
   /////////////////////////////////////////////////
   {
      FILE *fp;

      ////////////////////////////////////////////////////////////////////////////////////////////
      snprintf(filename,sizeof(filename),"%s_PIL.mrx",fprefix);
      fp = fopen(filename,"w");
      if(fp != NULL)
      {
         fprintf(fp,"# %s to midpoint rigid-body registration matrix computed by KAIBA",ffile);
         printMatrix(Tf, 4, 4, "", fp);
         fclose(fp);
      }
      else
      {
         printf("Warning: cound not write to %s\n", filename);
      }
      ////////////////////////////////////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////////////////////////////////////
      snprintf(filename,sizeof(filename),"%s_PIL.mrx",bprefix);
      fp = fopen(filename,"w");
      if(fp != NULL)
      {
         fprintf(fp,"# %s to midpoint rigid-body registration matrix computed by KAIBA",bfile);
         printMatrix(Tb, 4, 4, "", fp);
         fclose(fp);
      }
      else
      {
         printf("Warning: cound not write to %s\n", filename);
      }
      ////////////////////////////////////////////////////////////////////////////////////////////
   }
   /////////////////////////////////////////////////

   /////////////////////////////////////////////////
   // save registred images
   /////////////////////////////////////////////////
   {
      SHORTIM bimpil; // baseline image after transformation to standard PIL space
      SHORTIM fimpil; // follow-up image after transformation to standard PIL space

      invT = inv4(Tf);
      fimpil.v= resliceImage(fim, dimf, PILbraincloud_dim, invT, LIN);
      set_dim(fimpil, PILbraincloud_dim);
      strcpy(PILbraincloud_hdr.descrip,"Created by ART's KAIBA module");
      snprintf(filename,sizeof(filename),"%s_PIL.nii",fprefix);
      save_nifti_image(filename, fimpil.v, &PILbraincloud_hdr);
      free(invT);

      invT = inv4(Tb);
      bimpil.v = resliceImage(bim, dimb, PILbraincloud_dim, invT, LIN);
      set_dim(bimpil, PILbraincloud_dim);
      strcpy(PILbraincloud_hdr.descrip,"Created by ART's KAIBA module");
      snprintf(filename,sizeof(filename),"%s_PIL.nii",bprefix);
      save_nifti_image(filename, bimpil.v, &PILbraincloud_hdr);
      free(invT);

      set_dim(aimpil, PILbraincloud_dim);
      aimpil.v = (int2 *)calloc(aimpil.nv, sizeof(int2));
      for(int i=0; i<aimpil.nv; i++) 
         aimpil.v[i] = (int2)( (bimpil.v[i] + fimpil.v[i])/2.0 + 0.5 );

      delete bimpil.v;
      delete fimpil.v;
   }
   /////////////////////////////////////////////////

   delete sclbim;
   delete sclfim;
   delete bmsk;
   delete fmsk;
}

///////////////////////////////////////////////////////////////////////////////////////////////

void compute_lm_transformation(char *lmfile, SHORTIM im, float4 *A)
{
   FILE *fp;
   int NLM;
   int r;
   int R;
   float4 *LM; // 4xNLM matrix
   float4 *CM; // 4xNLM matrix
   int cm[3]; // landmarks center of mass
   int lm[3];

   fp=fopen(lmfile, "r");

   if(fp==NULL) 
   {
      printf("Could not find %s, aborting ...\n",lmfile);
      exit(0);
   }

   fread(&NLM, sizeof(int), 1, fp);
   fread(&r, sizeof(int), 1, fp);
   fread(&R, sizeof(int), 1, fp);
   SPH searchsph(R);
   SPH testsph(r);
   SPH refsph(r);
   LM = (float4 *)calloc(4*NLM, sizeof(float4));
   CM = (float4 *)calloc(4*NLM, sizeof(float4));

   for(int n=0; n<NLM; n++)
   {
      fread(&cm[0], sizeof(int), 1, fp);
      fread(&cm[1], sizeof(int), 1, fp);
      fread(&cm[2], sizeof(int), 1, fp);
      fread(refsph.v, sizeof(float4), refsph.n, fp);

      CM[0*NLM + n]=(cm[0] - (im.nx-1)/2.0)*im.dx; 
      CM[1*NLM + n]=(cm[1] - (im.ny-1)/2.0)*im.dy;
      CM[2*NLM + n]=(cm[2] - (im.nz-1)/2.0)*im.dz;
      CM[3*NLM + n]=1;

      detect_lm(searchsph, testsph, im, cm, refsph, lm);

      LM[0*NLM + n]=(lm[0] - (im.nx-1)/2.0)*im.dx; 
      LM[1*NLM + n]=(lm[1] - (im.ny-1)/2.0)*im.dy;
      LM[2*NLM + n]=(lm[2] - (im.nz-1)/2.0)*im.dz;
      LM[3*NLM + n]=1;
   }

   fclose(fp);

   leastSquaresAffineTrans(LM, CM, NLM, A);

   free(LM);
   free(CM);
}

///////////////////////////////////////////////////////////////////////////////////////////////

void find_roi(nifti_1_header *subimhdr, SHORTIM pilim, float4 pilT[],const char *side, const char *prefix, float *Tlm)
{
   DIM subdim;

   set_dim(subdim, subimhdr);

   // these are stored in $ARTHOME/<side>.nii file as hdr.dim[5,6,7]
   int imin;
   int jmin;
   int kmin;

   int2 *stndrd_roi;

   char filename[DEFAULT_STRING_LENGTH];

   // hcim is a pre-processed (transformed) version of subim (the original input image in native orientation)
   // readied for hippocampus segmentation.
   SHORTIM hcim; 
  
   // hcT is an affine transformation from subim to hcim
   float4 hcT[16];

   // hcim matrix and voxel dimensions are set to a starndard size
   set_dim(hcim, pilim);

   stndrd_roi = (int2 *)calloc(hcim.nv, sizeof(int2));

   multi(Tlm,4,4, pilT, 4,4, hcT);

   //snprintf(filename,sizeof(filename),"%s_%s.mrx",prefix,side);
   //fp=fopen(filename,"w");
   //printMatrix(hcT,4,4,"",fp);
   //fclose(fp);

   ////////////////////////////////////////////////////////////////////////////////////////////

   SHORTIM msk; 
   nifti_1_header mskhdr;
   
   snprintf(filename,sizeof(filename),"%s/%s.nii",ARTHOME,side);
   msk.v = (int2 *)read_nifti_image(filename, &mskhdr);

   if(msk.v==NULL) exit(0);

   msk.nx = mskhdr.dim[1];
   msk.ny = mskhdr.dim[2];
   msk.nz = mskhdr.dim[3];
   msk.np = msk.nx*msk.ny;
   msk.nv = msk.np*msk.nz;
   msk.dx = mskhdr.pixdim[1];
   msk.dy = mskhdr.pixdim[2];
   msk.dz = mskhdr.pixdim[3];
   imin = mskhdr.dim[5];
   jmin = mskhdr.dim[6];
   kmin = mskhdr.dim[7];
   //number of atlases: (mskhdr.dim[4]-1)/2;

   ////////////////////////////////////////////////////////////////////////////////////////////
   {
      int2 *ntv_spc_roi;
      float4 T[16];

      // if image was flipped (PIR), then 'lhc3' finds the RHROI
      // and 'rhc3' finds the LHROI
 /*
      if( side[0]=='r')
      {
         if(!flipped)
            snprintf(filename,sizeof(filename),"%s_RHROI1.nii",prefix);
         else
            snprintf(filename,sizeof(filename),"%s_LHROI2.nii",prefix);
      }
      if( side[0]=='l')
      {
         if(!flipped)
            snprintf(filename,sizeof(filename),"%s_LHROI1.nii",prefix);
         else
            snprintf(filename,sizeof(filename),"%s_RHROI2.nii",prefix);
      }
*/
      snprintf(filename,sizeof(filename),"%s",prefix);

      for(int n=0; n<hcim.nv; n++) stndrd_roi[n] = 0;

      for(int k=0; k<msk.nz; k++)
      for(int j=0; j<msk.ny; j++)
      for(int i=0; i<msk.nx; i++)
      {
         stndrd_roi[(k+kmin)*hcim.np + (j+jmin)*hcim.nx + (i+imin)] = msk.v[k*msk.np + j*msk.nx + i];
      }

      for(int i=0; i<16; i++) T[i]=hcT[i];

      ntv_spc_roi = resliceImage(stndrd_roi, hcim.nx, hcim.ny, hcim.nz, hcim.dx, hcim.dy, hcim.dz,
      subdim.nx, subdim.ny, subdim.nz, subdim.dx, subdim.dy, subdim.dz, T, LIN);

      save_nifti_image(filename, ntv_spc_roi, subimhdr);

      free(ntv_spc_roi); free(stndrd_roi);
   }

   return;
}

///////////////////////////////////////////////////////////////////////////////////////////////

void compute_hi(char *imfile, char *roifile, float4 &parenchymasize, int &voisize)
{
   int binw; // histogram bin width
   int gm_pk_srch_strt;
   int2 *roi;
   int2 *im;
   nifti_1_header hdr;
   int nx, ny, nz, np, nv;
   float4 dx, dy, dz;
   int I_alpha;
   int nbin;
   char roifileprefix[DEFAULT_STRING_LENGTH]=""; //baseline image prefix
   char roifiledir[DEFAULT_STRING_LENGTH]=""; //baseline image prefix
   char filename[DEFAULT_STRING_LENGTH]=""; //baseline image prefix

   if( niftiFilename(roifileprefix, roifile)==0 ) exit(1);
   getDirectoryName(roifile, roifiledir);

   if(opt_v)
   {
      printf("Computing HPF ...\n");
      printf("Image file: %s\n", imfile);
      printf("ROI file: %s\n", roifile);
   }

   roi = (int2 *)read_nifti_image(roifile, &hdr);
   nx = hdr.dim[1];
   ny = hdr.dim[2];
   nz = hdr.dim[3];
   dx = hdr.pixdim[1];
   dy = hdr.pixdim[2];
   dz = hdr.pixdim[3];
   nv = nx*ny*nz;
   np = nx*ny;

   im = (int2 *)read_nifti_image(imfile, &hdr);
   setMX(im, roi, nv, I_alpha, alpha_param);

   if(opt_v) printf("I_alpha = %d\n",I_alpha);

   for(int i=0; i<nv; i++)
      if( im[i] > I_alpha ) roi[i]=0;

   voisize = 0;
   for(int i=0; i<nv; i++) 
   {
      if( roi[i]>0 ) voisize++;
   }

//   if(opt_v) printf("ROI size = %d voxels\n", voisize);

   /////////////////////////////////////////////////////////////
   int im_min=0, im_max=0;
   float8 *hist;
   float8 *fit;
   float8 mean[MAXNCLASS+1];
   float8 var[MAXNCLASS+1];
   float8 p[MAXNCLASS+1];
   int2 *label;
   int I_gm=0;
   int I_csf;
   
   // initialize min and max variables
   for(int i=0; i<nv; i++)
   {
      if( roi[i] > 0 )
      {
         im_min=im_max=im[i];
         break;
      }
   }

   // find im_min and im_max amongst the core voxels
   for(int i=0; i<nv; i++)
   {
      if( roi[i] > 0 )
      {
         if( im[i]<im_min ) im_min=im[i];
         else if( im[i]>im_max ) im_max=im[i];
      }
   }

   if(opt_v) printf("Image intensity range within ROI: min=%d max=%d\n",im_min,im_max);

   if(im_max<NBIN) {nbin=im_max+1; binw=1; }
   else { nbin = NBIN; binw = im_max/nbin + 1; }

   hist = (float8 *)calloc(nbin, sizeof(float8));
   fit = (float8 *)calloc(nbin, sizeof(float8));
   label = (int2 *)calloc(nbin, sizeof(int2));

   // initialize hist to 0 (to be sure)
   for(int i=0; i<nbin; i++) hist[i]=0.0;

   // set hist of the core voxels
   for(int i=0; i<nv; i++)
   {
      if( roi[i] > 0)
      {
         hist[im[i]/binw]++;
      }
   }

   while( hist[nbin-1] == 0 ) nbin--;
   if(opt_v) printf("Number of histogram bins = %d\n",nbin);
   if(opt_v) printf("Histogram bin width = %d\n",binw);

   // normalize hist
   for(int i=0; i<nbin; i++) hist[i]/=voisize;

   gm_pk_srch_strt = (MXFRAC*I_alpha/binw);
   if( gm_pk_srch_strt < 0) gm_pk_srch_strt=0;

   int nclass=5;

   EMFIT1d(hist, fit, label, nbin, mean, var, p, nclass, 1000);

   // find I_gm
   float8 hmax=0.0;
   for(int i=gm_pk_srch_strt; i<nbin; i++)
   {
      if( fit[i] > hmax ) 
      { 
         hmax=fit[i]; 
         I_gm=i; 
      }
   }

#if 0
   // implementation of minimum error threshold selection method
   // Do not remove
   
   // find GM class
   int gmclass;
   {
      float8 del;
      del = fabs(I_gm - mean[0]);
      gmclass=0;

      for(int i=1; i<nclass; i++) 
      {
         if( fabs(I_gm - mean[i]) < del ) 
         {
            del = fabs(I_gm - mean[i]);
            gmclass=i;
         }
      }
   }

   {
      float8 total_error;
      float8 min_total_error=1.0;
      int x_at_min_total_error;

      for(int x=0; x<I_gm; x++)
      {
         total_error=0.0;
         for(int i=0; i<gmclass; i++)
         {
            if(p[i]>0.0) total_error += p[i]*(1.0 - normalCDF(x,mean[i],var[i]) );
         }
         for(int i=gmclass; i<nclass; i++)
         {
            if(p[i]>0.0) total_error += p[i]*normalCDF(x,mean[i],var[i]);
         }
         if( total_error < min_total_error )
         {
            min_total_error = total_error;
            x_at_min_total_error = x;
         }
      }

      I_csf = x_at_min_total_error + im_min;
   }
#endif

   // Al's GOLDEN method for finding the I_csf
   I_csf = (int)(I_gm - BETA_PARAM*I_alpha/binw + 0.5);

   // New method that only relies on I_gm
   // The jury is still out whether this is better than the above method
   //I_csf =(int)(0.67*I_gm + 0.5);

   //float K;
   //K = (BETA_PARAM*I_alpha/binw)/I_gm;
   /////////////////////////////////////////
   
   if(opt_v)
   {
      printf("A histogram peak detected at: %d\n",I_gm);
      printf("Image intensity threshold: %d\n",I_csf);
      //printf("K=%f\n",K);
   }

   {
      int n=nbin;
      float *data1, *data2;
      int *bin;

      data1 = (float *)calloc(n, sizeof(float));
      data2 = (float *)calloc(n, sizeof(float));
      bin = (int *)calloc(n, sizeof(int));

      for(int i=0; i<n; i++)
      {
         data1[i] = (float)hist[i];
         data2[i] = (float)fit[i];
         bin[i] = i;
      }

      snprintf(filename,sizeof(filename),"%s/%s_hist",roifiledir,roifileprefix);
      hist1D_plot(filename, n, bin, data1, data2, I_gm, I_csf);

      free(data1);
      free(data2);
   }

   /////////////////////////////////////////////////////////////
   
   float8 csfvol=0.0;

   for(int i=0; i<I_csf; i++)
   {
      csfvol += hist[i];
   }

   csfvol += 0.5*hist[I_csf];

   free(hist); free(fit); free(label);
   /////////////////////////////////////////////////////////////

   if(opt_v)
   {
      printf("HPF = %6.4f\n",1.0-csfvol);
   }

   parenchymasize = (float)((1.0 - csfvol)*voisize);
}

///////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   char orient[4]="";
   float Tleft0[16], Tleft1[16], Tright0[16], Tright1[16];
   float4 *invT;
   SHORTIM aimpil; // average of input images after transformation to standard PIL space
   SHORTIM im[MAXIM];
   DIM im_dim[MAXIM];
   nifti_1_header im_hdr[MAXIM];
   DIM PILbraincloud_dim;
   float TPIL[MAXIM][16];
   //float4 RHI0, RHI1, LHI0, LHI1;
   float RHI, LHI, BHI, HIasymm;
   int Rvoisize0[MAXIM], Rvoisize1[MAXIM], Lvoisize0[MAXIM], Lvoisize1[MAXIM];
   float Rparenchymasize0[MAXIM], Rparenchymasize1[MAXIM], Lparenchymasize0[MAXIM], Lparenchymasize1[MAXIM];

   int nim=0; // number of input images
   char **imagefile; // the nim input image files
   char **mrxfile; // the nim input image files
   char **imagefileprefix; 
   char **imagedir; 
   float *scalefactor;

   short *tmp;

   FILE *fp;
   char filename[DEFAULT_STRING_LENGTH]="";  // a generic filename for reading/writing stuff

   char opprefix[DEFAULT_STRING_LENGTH]=""; // prefix used for reading/writing output files

   char roifile[DEFAULT_STRING_LENGTH]="";

   char lmfile[DEFAULT_STRING_LENGTH]="";

   char bfile[DEFAULT_STRING_LENGTH]=""; // baseline image filename

   if(argc==1) print_help_and_exit();

   alpha_param = ALPHA_PARAM;

   while ((opt = getoption(argc, argv, options)) != -1 )
   {
      switch (opt) 
      {
         case 'V':
            printf("KAIBA Version 3.1 released Dec. 18, 2019.\n");
            printf("Author: Babak A. Ardekani, Ph.D.\n");
            exit(0);
         case 't':
            opt_txt = NO;
            break;
         case 'N':
            opt_ppm=NO;
            break;
         case 'F':
            opt_flip=NO;
            break;
         case 'v':
            opt_v=YES;
            break;
         case 'g':
            opt_png=NO;
            break;
         case 'p':
            strcpy(opprefix,optarg);
            break;
         case 'l':
            strcpy(lmfile,optarg);
            break;
         case 'b':
            strcpy(bfile,optarg);
            break;
         case 'a':
            alpha_param = atof(optarg); 
            break;
         case '?':
            print_help_and_exit();
      }
   }

   getARTHOME();
   if(opt_v) printf("ARTHOME = %s\n", ARTHOME);

   /////////////////////////////////////////////////////////////////////////////////////
   // input image file names and transformation file names
   /////////////////////////////////////////////////////////////////////////////////////
   // Ensure that a baseline image has been specified at the command line.
   if( bfile[0]=='\0' )
   {
      printf("Please specify an image or an image list using -i argument. As in:\n");
      printf("$ kaiba -i <image>.nii\nor\n");
      printf("$ kaiba -i <imagelist>\n");
      exit(0);
   }

   if( not_magical_nifti(bfile,0)==0 )  // a single image was specified using -i <image>.nii
   {
      nim=1;
   }
   else // an image list was specified using -i <imagelistfile>
   {
      fp=fopen(bfile,"r");
      if(fp==NULL) file_open_error(bfile);
      fscanf(fp,"%d",&nim);
      fclose(fp);
   }

   if(opt_v) printf("Number of input images = %d\n\n", nim);
   if(nim<1) exit(0);

   imagefile = (char **)calloc(nim, sizeof(char *));
   mrxfile = (char **)calloc(nim, sizeof(char *));
   imagefileprefix = (char **)calloc(nim, sizeof(char *));
   imagedir = (char **)calloc(nim, sizeof(char *));
   scalefactor = (float *)calloc(nim, sizeof(float));

   for(int i=0; i<nim; i++)
   {
      imagefile[i] = (char *)calloc(DEFAULT_STRING_LENGTH, sizeof(char));
      imagefile[i][0]='\0';

      mrxfile[i] = (char *)calloc(DEFAULT_STRING_LENGTH, sizeof(char));
      mrxfile[i][0]='\0';

      imagefileprefix[i] = (char *)calloc(DEFAULT_STRING_LENGTH, sizeof(char));
      imagefileprefix[i][0]='\0';

      imagedir[i] = (char *)calloc(DEFAULT_STRING_LENGTH, sizeof(char));
      imagedir[i][0]='\0';
   }

   if( not_magical_nifti(bfile,0)==0 )  // a single image was specified using -i <image>.nii
   {
      strcpy(imagefile[0], bfile);
      if( niftiFilename(imagefileprefix[0], imagefile[0])==0 ) { exit(0); }
      scalefactor[0]=1.0;
      getDirectoryName(imagefile[0], imagedir[0]);

      if(opt_v) printf("Input image: %s\n", imagefile[0]);
   }
   else
   {
      //////////////////////////////////////////////////////////////////////////////////
      // fill imagefile mrxfile arrays
      //////////////////////////////////////////////////////////////////////////////////
      fp=fopen(bfile,"r");
      if(fp==NULL) file_open_error(bfile);
      fscanf(fp,"%d",&nim);
      for(int i=0; i<nim; i++)
      {
         fscanf(fp,"%s",imagefile[i]);
         fscanf(fp,"%f",&scalefactor[i]);
         fscanf(fp,"%s",mrxfile[i]);
         if( niftiFilename(imagefileprefix[i], imagefile[i])==0 ) { exit(0); }
         getDirectoryName(imagefile[i], imagedir[i]);
      }
      fclose(fp);

      if(opt_v)
      {
         for(int i=0; i<nim; i++)
         {
            printf("Input image %d: %s\n",i+1, imagefile[i]);
            printf("Transformation matrix: %s\n",mrxfile[i]);
            printf("Scale factor: %f\n\n",scalefactor[i]);
         }
      }
   }
   /////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////

printf("Point1\n");
   // Ensure that an output prefix has been specified at the command line.
   if( opprefix[0]=='\0' )
   {
      printf("Please specify an output prefix using argument: -o <prefix>\n");
      exit(0);
   }

printf("Point2\n");
   /////////////////////////////////////////////////////////////////////////////////////
   // read input images
   /////////////////////////////////////////////////////////////////////////////////////
   for(int i=0; i<nim; i++)
   {
      im[i].v = (int2 *)read_nifti_image(imagefile[i], &im_hdr[i]);
      if(im[i].v==NULL)
      {
         printf("Error reading %s, aborting ...\n", imagefile[i]);
         exit(1);
      }
      set_dim(im_dim[i], im_hdr[i]);
   }
   /////////////////////////////////////////////////////////////////////////////////////
printf("Point3\n");

   /////////////////////////////////////////////////////////////////////////////////////////////
   // read PILbraincloud.nii from the $ARTHOME directory
   // The only reason this is done is to read dimensions of PILbrain.nii
   /////////////////////////////////////////////////////////////////////////////////////////////
   int2 *PILbraincloud;
   nifti_1_header PILbraincloud_hdr; 

   strcpy(filename,""); strcat(filename,ARTHOME); strcat(filename,"/PILbrain.nii");

   PILbraincloud = (int2 *)read_nifti_image(filename, &PILbraincloud_hdr);

   if(PILbraincloud==NULL)
   {
         printf("Error reading %s, aborting ...\n", filename);
         exit(1);
   }

   set_dim(PILbraincloud_dim, PILbraincloud_hdr);
   set_dim(aimpil, PILbraincloud_dim); 
   delete PILbraincloud;
   /////////////////////////////////////////////////////////////////////////////////////////////
printf("Point4\n");

   /////////////////////////////////////////////////////////////////////////////////////
   // average input images after PIL transformation and store in aimpil
   /////////////////////////////////////////////////////////////////////////////////////
   //if(nim>1)
   if( not_magical_nifti(bfile,0)==1 )  // uses the output from ATRA
   {
      float *sum;
      short *tmp;

      aimpil.v = (short *)calloc(aimpil.nv, sizeof(short));
      sum = (float *)calloc(aimpil.nv, sizeof(float));
      for(int v=0; v<aimpil.nv; v++) sum[v]=0.0;

      for(int i=0; i<nim; i++)
      {
         loadTransformation(mrxfile[i], TPIL[i]);
         invT = inv4(TPIL[i]);
         tmp = resliceImage(im[i].v, im_dim[i], PILbraincloud_dim, invT, LIN);

         for(int v=0; v<aimpil.nv; v++) sum[v] += tmp[v]*scalefactor[i];

         free(invT);
         free(tmp);
      }

      for(int v=0; v<aimpil.nv; v++) aimpil.v[v] = (short)(sum[v]/nim + 0.5);
      free(sum);

      strcpy(filename,""); strcat(filename,opprefix); strcat(filename,"_avg_PIL.nii");
      save_nifti_image(filename, aimpil.v, &PILbraincloud_hdr);
   } 
   else if (nim==1) 
   {
      if(opt_v) printf("Computing PIL transformation ...\n");
      orient[0]='\0';
      new_PIL_transform(imagefile[0],lmfile, orient, TPIL[0], 0);

      invT = inv4(TPIL[0]);
      aimpil.v = resliceImage(im[0].v, im_dim[0], PILbraincloud_dim, invT, LIN);
      free(invT);
   } 
   /////////////////////////////////////////////////////////////////////////////////////
printf("Point5\n");
   
   /////////////////////////////////////////////////////////////////////////////////////
   
   strcpy(filename,""); strcat(filename,opprefix); strcat(filename,".csv");

   fp = fopen(filename,"a");
   if(fp==NULL) file_open_error(filename);
   fprintf(fp,"Volume,Hemisphere,HPF\n");

   if(opt_v) printf("\nLandmark detection ...\n");
   strcpy(filename,""); strcat(filename,ARTHOME); strcat(filename,"/rhc3.mdl");
   compute_lm_transformation(filename, aimpil, Tright0);
   strcpy(filename,""); strcat(filename,ARTHOME); strcat(filename,"/lhc3.mdl");
   compute_lm_transformation(filename, aimpil, Tleft0);

   flipped=NO;
   for(int i=0; i<nim; i++)
   {
      //find_roi(&im_hdr[i], aimpil, TPIL[i], "rhc3", imagefileprefix[i], Tright0);
      snprintf(roifile,sizeof(roifile),"%s/%s_RHROI1.nii",imagedir[i],imagefileprefix[i]);
      find_roi(&im_hdr[i], aimpil, TPIL[i], "rhc3", roifile, Tright0);
      compute_hi(imagefile[i], roifile, Rparenchymasize0[i], Rvoisize0[i]);

      //find_roi(&im_hdr[i], aimpil, TPIL[i], "lhc3", imagefileprefix[i], Tleft0);
      snprintf(roifile,sizeof(roifile),"%s/%s_LHROI1.nii",imagedir[i],imagefileprefix[i]);
      find_roi(&im_hdr[i], aimpil, TPIL[i], "lhc3", roifile, Tleft0);
      compute_hi(imagefile[i], roifile, Lparenchymasize0[i], Lvoisize0[i]);
   }

   if(opt_flip) 
   { 
      flipped=YES;
      float I[16];

      set_to_I(I,4);
      I[10]=-1.0;
      tmp = resliceImage(aimpil.v, PILbraincloud_dim, PILbraincloud_dim, I, LIN);
      free(aimpil.v);
      aimpil.v=tmp;

      if(opt_v) printf("\nLandmark detection ...\n");
      snprintf(filename,sizeof(filename),"%s/%s.mdl",ARTHOME,"lhc3");
      compute_lm_transformation(filename, aimpil, Tright1);
      snprintf(filename,sizeof(filename),"%s/%s.mdl",ARTHOME,"rhc3");
      compute_lm_transformation(filename, aimpil, Tleft1);

      for(int i=0; i<nim; i++)
      {
         TPIL[i][8]*=-1.0; TPIL[i][9]*=-1.0; TPIL[i][10]*=-1.0; TPIL[i][11]*=-1.0; 

         //find_roi(&im_hdr[i], aimpil, TPIL[i], "lhc3", imagefileprefix[i],Tright1);
         snprintf(roifile,sizeof(roifile),"%s/%s_RHROI2.nii",imagedir[i],imagefileprefix[i]);
         find_roi(&im_hdr[i], aimpil, TPIL[i], "lhc3", roifile,Tright1);
         compute_hi(imagefile[i], roifile, Rparenchymasize1[i], Rvoisize1[i]);

         //find_roi(&im_hdr[i], aimpil, TPIL[i], "rhc3", imagefileprefix[i],Tleft1);
         snprintf(roifile,sizeof(roifile),"%s/%s_LHROI2.nii",imagedir[i],imagefileprefix[i]);
         find_roi(&im_hdr[i], aimpil, TPIL[i], "rhc3", roifile,Tleft1);
         compute_hi(imagefile[i], roifile, Lparenchymasize1[i], Lvoisize1[i]);

         // important: restore TPIL
         TPIL[i][8]*=-1.0; TPIL[i][9]*=-1.0; TPIL[i][10]*=-1.0; TPIL[i][11]*=-1.0; 
      }
   }
   else
   {
      for(int i=0; i<nim; i++)
      {
         Rvoisize1[i]=Rvoisize0[i]; Lvoisize1[i]=Lvoisize0[i];
         Rparenchymasize1[i]=Rparenchymasize0[i]; Lparenchymasize1[i]=Lparenchymasize0[i];
      }
   }

   for(int i=0; i<nim; i++) 
   {
      //RHI1 = Rparenchymasize1[i]/Rvoisize1[i];
      //RHI0 = Rparenchymasize0[i]/Rvoisize0[i];
      //fprintf(fp,"%s,Right,%f*\n",imagefile[i],(RHI1+RHI0)/2.0);
      RHI = (Rparenchymasize1[i]+Rparenchymasize0[i])/(Rvoisize1[i]+Rvoisize0[i]);
      fprintf(fp,"%s,Right,%f\n",imagefile[i],RHI);
   }

   for(int i=0; i<nim; i++) 
   {
      //LHI1 = Lparenchymasize1[i]/Lvoisize1[i];
      //LHI0 = Lparenchymasize0[i]/Lvoisize0[i];
      //fprintf(fp,"%s,Left,%f*\n",imagefile[i],(LHI1+LHI0)/2.0);
      LHI = (Lparenchymasize1[i]+Lparenchymasize0[i])/(Lvoisize1[i]+Lvoisize0[i]);
      fprintf(fp,"%s,Left,%f\n",imagefile[i],LHI);
   }

   for(int i=0; i<nim; i++) 
   {
      BHI = (Lparenchymasize1[i]+Lparenchymasize0[i]+Rparenchymasize1[i]+Rparenchymasize0[i]);
      BHI /= (Rvoisize1[i]+Rvoisize0[i]+Lvoisize1[i]+Lvoisize0[i]);
      fprintf(fp,"%s,Bilateral,%f\n",imagefile[i],BHI);
   }

   for(int i=0; i<nim; i++) 
   {
      RHI = (Rparenchymasize1[i]+Rparenchymasize0[i])/(Rvoisize1[i]+Rvoisize0[i]);
      LHI = (Lparenchymasize1[i]+Lparenchymasize0[i])/(Lvoisize1[i]+Lvoisize0[i]);
      BHI = (Lparenchymasize1[i]+Lparenchymasize0[i]+Rparenchymasize1[i]+Rparenchymasize0[i]);
      BHI /= (Rvoisize1[i]+Rvoisize0[i]+Lvoisize1[i]+Lvoisize0[i]);
      HIasymm = 100.0*(RHI-LHI)/BHI;
      fprintf(fp,"%s,Asymm,%f\n",imagefile[i],HIasymm);
   }

   delete aimpil.v;

   fclose(fp);

#if 0
#endif
  for(int i=0; i<nim; i++)
  {
    short *roi1, *roi2, *roi, *roiPIL;
    nifti_1_header hdr;
    DIM roidim;
    DIM opdim;
    float cm[3];
    double C[6];  
    double L[3];
    double UT[9];
    float T[16];
    float tmpT[16];
    float alpha;
    short *im, *imHC;

    // set dimensions of the output hippocampus image
    set_dim(opdim, 59, 29, 45, 1.0, 1.0, 1.0);

    // read roi1 
    snprintf(roifile,sizeof(roifile),"%s/%s_LHROI1.nii",imagedir[i],imagefileprefix[i]);
    roi1 = (int2 *)read_nifti_image(roifile, &hdr);

    // read roi2 
    snprintf(roifile,sizeof(roifile),"%s/%s_LHROI2.nii",imagedir[i],imagefileprefix[i]);
    roi2 = (int2 *)read_nifti_image(roifile, &hdr);
   
    set_dim(roidim, hdr);

    // combine roi1 and roi2
    roi = roi1;
    for(int v=0; v<roidim.nv; v++) { roi[v] = roi1[v] + roi2[v]; }
 
    // transform compined roi to PIL space
    invT = inv4(TPIL[i]);
    roiPIL = resliceImage(roi, roidim, PILbraincloud_dim, invT, LIN);
    free(invT);

    // find the center of mass (cm) and covrariance matrix (C) of 
    // the roi in PIL space, then find it's first principal axis
    compute_cov(roiPIL, PILbraincloud_dim, cm, C);
    s3eigenval(C, L);
    s3eigenvec(C, L, UT);
    // ensure that the first principal axis points posteriorly
    if(UT[0]<0.0) { UT[0]*=-1.0; UT[1]*=-1; UT[2]*=-1; }

    // modify TPIL and store in T so that the hippocampus center
    // of mass will be the center of the FOV
    tmpT[0]=1.0;  tmpT[1]=0.0;  tmpT[2]=0.0;  tmpT[3]=-cm[0];
    tmpT[4]=0.0;  tmpT[5]=1.0;  tmpT[6]=0.0;  tmpT[7]=-cm[1];
    tmpT[8]=0.0;  tmpT[9]=0.0;  tmpT[10]=1.0; tmpT[11]=-cm[2];
    tmpT[12]=0.0; tmpT[13]=0.0; tmpT[14]=0.0; tmpT[15]=1.0;
    multi(tmpT, 4, 4, TPIL[i], 4, 4, T);

    // furture modify T so that the hippocampus first principal
    // axis would be aligned with the x-axis
    if(UT[0]>1.0) UT[0]=1.0; // just in case to prevent acos from getting into trouble
    alpha = (float)acos((double)UT[0]);
    rotate(tmpT, alpha, 0, UT[2], -UT[1]);
    multi(tmpT, 4, 4, T, 4, 4, T);

    set_dim(PILbraincloud_hdr,opdim);

    free(roiPIL);
    invT = inv4(T);
    roiPIL = resliceImage(roi, roidim, opdim, invT, LIN);
    free(invT);

    snprintf(roifile,sizeof(roifile),"%s/%s_LHC_ROI.nii",imagedir[i],imagefileprefix[i]);
    save_nifti_image(roifile, roiPIL, &PILbraincloud_hdr);

    im = (int2 *)read_nifti_image(imagefile[i], &hdr);
    invT = inv4(T);
    imHC = resliceImage(im, roidim, opdim, invT, LIN);
    free(invT);

    snprintf(roifile,sizeof(roifile),"%s/%s_LHC.nii",imagedir[i],imagefileprefix[i]);
    save_nifti_image(roifile, imHC, &PILbraincloud_hdr);

    free(im); free(imHC);
    free(roi1); free(roi2); free(roiPIL);
    /////////////////////////////////////////////////////////////////////////
    // Processing the right side
    
    // read roi1 
    snprintf(roifile,sizeof(roifile),"%s/%s_RHROI1.nii",imagedir[i],imagefileprefix[i]);
    roi1 = (int2 *)read_nifti_image(roifile, &hdr);

    // read roi2 
    snprintf(roifile,sizeof(roifile),"%s/%s_RHROI2.nii",imagedir[i],imagefileprefix[i]);
    roi2 = (int2 *)read_nifti_image(roifile, &hdr);
   
    set_dim(roidim, hdr);

    // combine roi1 and roi2
    roi = roi1;
    for(int v=0; v<roidim.nv; v++) { roi[v] = roi1[v] + roi2[v]; }
 
    // transform compined roi to PIL space
    invT = inv4(TPIL[i]);
    roiPIL = resliceImage(roi, roidim, PILbraincloud_dim, invT, LIN);
    free(invT);

    // find the center of mass (cm) and covrariance matrix (C) of 
    // the roi in PIL space, then find it's first principal axis
    compute_cov(roiPIL, PILbraincloud_dim, cm, C);
    s3eigenval(C, L);
    s3eigenvec(C, L, UT);
    // ensure that the first principal axis points posteriorly
    if(UT[0]<0.0) { UT[0]*=-1.0; UT[1]*=-1; UT[2]*=-1; }

    // modify TPIL and store in T so that the hippocampus center
    // of mass will be the center of the FOV
    tmpT[0]=1.0;  tmpT[1]=0.0;  tmpT[2]=0.0;  tmpT[3]=-cm[0];
    tmpT[4]=0.0;  tmpT[5]=1.0;  tmpT[6]=0.0;  tmpT[7]=-cm[1];
    tmpT[8]=0.0;  tmpT[9]=0.0;  tmpT[10]=1.0; tmpT[11]=-cm[2];
    tmpT[12]=0.0; tmpT[13]=0.0; tmpT[14]=0.0; tmpT[15]=1.0;
    multi(tmpT, 4, 4, TPIL[i], 4, 4, T);

    // furture modify T so that the hippocampus first principal
    // axis would be aligned with the x-axis
    if(UT[0]>1.0) UT[0]=1.0; // just in case to prevent acos from getting into trouble
    alpha = (float)acos((double)UT[0]);
    rotate(tmpT, alpha, 0, UT[2], -UT[1]);
    multi(tmpT, 4, 4, T, 4, 4, T);

    set_dim(PILbraincloud_hdr,opdim);

    // coverts T to PIR
    T[8]*=-1.0; T[9]*=-1.0; T[10]*=-1.0; T[11]*=-1.0; 

    // converts the header from PIL to PIR
    PILbraincloud_hdr.pixdim[0] *= -1.0;

    free(roiPIL);
    invT = inv4(T);
    roiPIL = resliceImage(roi, roidim, opdim, invT, LIN);
    free(invT);

    snprintf(roifile,sizeof(roifile),"%s/%s_RHC_ROI.nii",imagedir[i],imagefileprefix[i]);
    save_nifti_image(roifile, roiPIL, &PILbraincloud_hdr);

    im = (int2 *)read_nifti_image(imagefile[i], &hdr);
    invT = inv4(T);
    imHC = resliceImage(im, roidim, opdim, invT, LIN);
    free(invT);

    snprintf(roifile,sizeof(roifile),"%s/%s_RHC.nii",imagedir[i],imagefileprefix[i]);
    save_nifti_image(roifile, imHC, &PILbraincloud_hdr);

    // reverts the header to PIL 
    PILbraincloud_hdr.pixdim[0] *= -1.0;

    free(im); free(imHC);
    free(roi1); free(roi2); free(roiPIL);
  }
}
