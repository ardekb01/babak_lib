#include <stdlib.h>
#include <nifti1.h>

extern void PILtransform(const char *orientCode, float *orientMat);
extern void inversePILtransform(const char *orientCode, float *orientMat);
extern void multi(float *A,int iA,int jA,float *B,int iB,int jB,float *C);
extern void multi(double *A,int iA,int jA,double *B,int iB,int jB,double *C);
extern void multi(float *A,int iA,int jA, double *B,int iB,int jB,double *C);
extern void multi(double *A,int iA,int jA, float *B,int iB,int jB,float *C);
extern void memory_allocation_error(const char *s);
extern void getNiftiImageOrientation(const char *filename, char *orientation);
extern void getNiftiImageOrientation(nifti_1_header hdr, char *orientation);
extern void update_qsform(nifti_1_header &hdr, const char *orientationcode);

// T =
// T[0] T[1] T[2]  0
// T[4] T[5] T[6]  0
// T[8] T[9] T[10] 0
// 0    0    0     1
//
// where T[0], T[1], T[2], T[4], T[5], T[6], T[8], T[9], T[10] can be either 1 or -1
// Also, the first three columns (rows) can only have one non-zero (1 or -1) value

short *reorientVolume(short *input_image, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1, float *T,
int &nx2, int &ny2, int &nz2, float &dx2, float &dy2, float &dz2)
{
   // set nx2, ny2, nz2, dx2, dy2, dz2
   if     ( T[0] != 0.0 ) 	{ nx2 = nx1; dx2 = dx1; }
   else if( T[4] != 0.0 ) 	{ ny2 = nx1; dy2 = dx1; }
   else if( T[8] != 0.0 ) 	{ nz2 = nx1; dz2 = dx1; }

   if     ( T[1] != 0.0 ) 	{ nx2 = ny1; dx2 = dy1; }
   else if( T[5] != 0.0 ) 	{ ny2 = ny1; dy2 = dy1; }
   else if( T[9] != 0.0 ) 	{ nz2 = ny1; dz2 = dy1; }

   if     ( T[2] != 0.0 ) 	{ nx2 = nz1; dx2 = dz1; }
   else if( T[6] != 0.0 ) 	{ ny2 = nz1; dy2 = dz1; }
   else if( T[10]!= 0.0 )   { nz2 = nz1; dz2 = dz1; }

   // allocate memory for the output image
   short *output_image;
   output_image  = (short *)calloc(nx2*ny2*nz2, sizeof(short));
   if(output_image==NULL)
   {
      memory_allocation_error("output_image");
   }

   int xc,yc,zc;
 
   if ( T[0] < 0.0 || T[4] < 0.0 || T[8] < 0.0 )
      xc=nx1-1;
   else
      xc=0;

   if ( T[1] < 0.0 || T[5] < 0.0 || T[9] < 0.0 )
      yc=ny1-1;
   else
      yc=0;

   if ( T[2] < 0.0 || T[6] < 0.0 || T[10] < 0.0 )
      zc=nz1-1;
   else
      zc=0;

   int q=0;
   int np1 = nx1 * ny1;
   int i1_a, j1_a, k1_a; 
   int i1_b, j1_b, k1_b; 
   int i1, j1, k1; 
   for(int k2=0; k2<nz2; k2++)
   {
      i1_a = (int)T[8] * k2 + xc;
      j1_a = (int)T[9] * k2 + yc;
      k1_a = (int)T[10]* k2 + zc;

      for(int j2=0; j2<ny2; j2++)
      {
         i1_b = (int)T[4] * j2 + i1_a;
         j1_b = (int)T[5] * j2 + j1_a;
         k1_b = (int)T[6] * j2 + k1_a;

         for(int i2=0; i2<nx2; i2++)
         {
            i1 = (int)T[0] * i2 + i1_b;
            j1 = (int)T[1] * i2 + j1_b;
            k1 = (int)T[2] * i2 + k1_b;

            output_image[q++]= input_image[k1*np1 + j1*nx1 + i1];
         }
      }
   }

   return(output_image);
}

short *reorientVolume(short *input_image, nifti_1_header oldhdr, const char *neworient, nifti_1_header &newhdr,
float *T_oldorient_to_neworient)
{
   // we will make some changes to newhdr later in this function
   newhdr = oldhdr;

   char oldorient[4];
   getNiftiImageOrientation(oldhdr, oldorient);

   float T_oldorient_to_PIL[16];
   float T_PIL_to_neworient[16];
   PILtransform(oldorient, T_oldorient_to_PIL);
   inversePILtransform(neworient, T_PIL_to_neworient);
   multi(T_PIL_to_neworient,4,4, T_oldorient_to_PIL,4,4,T_oldorient_to_neworient);

   int nx1=0, ny1=0, nz1=0;
   int nx2=0, ny2=0, nz2=0;

   float dx1=0.0, dy1=0.0, dz1=0.0;
   float dx2=0.0, dy2=0.0, dz2=0.0;

   nx1=oldhdr.dim[1];
   ny1=oldhdr.dim[2];
   nz1=oldhdr.dim[3];
   dx1=oldhdr.pixdim[1];
   dy1=oldhdr.pixdim[2];
   dz1=oldhdr.pixdim[3];

   // set nx2, ny2, nz2, dx2, dy2, dz2
   if     ( T_oldorient_to_neworient[0] != 0.0 ) 	{ nx2 = nx1; dx2 = dx1; }
   else if( T_oldorient_to_neworient[4] != 0.0 ) 	{ ny2 = nx1; dy2 = dx1; }
   else if( T_oldorient_to_neworient[8] != 0.0 ) 	{ nz2 = nx1; dz2 = dx1; }

   if     ( T_oldorient_to_neworient[1] != 0.0 ) 	{ nx2 = ny1; dx2 = dy1; }
   else if( T_oldorient_to_neworient[5] != 0.0 ) 	{ ny2 = ny1; dy2 = dy1; }
   else if( T_oldorient_to_neworient[9] != 0.0 ) 	{ nz2 = ny1; dz2 = dy1; }

   if     ( T_oldorient_to_neworient[2] != 0.0 ) 	{ nx2 = nz1; dx2 = dz1; }
   else if( T_oldorient_to_neworient[6] != 0.0 ) 	{ ny2 = nz1; dy2 = dz1; }
   else if( T_oldorient_to_neworient[10]!= 0.0 )   { nz2 = nz1; dz2 = dz1; }

   // allocate memory for the output image
   short *output_image;
   output_image  = (short *)calloc(nx2*ny2*nz2, sizeof(short));
   if(output_image==NULL)
   {
      memory_allocation_error("output_image");
   }

   int xc,yc,zc;
 
   if ( T_oldorient_to_neworient[0] < 0.0 || T_oldorient_to_neworient[4] < 0.0 || T_oldorient_to_neworient[8] < 0.0 )
      xc=nx1-1;
   else
      xc=0;

   if ( T_oldorient_to_neworient[1] < 0.0 || T_oldorient_to_neworient[5] < 0.0 || T_oldorient_to_neworient[9] < 0.0 )
      yc=ny1-1;
   else
      yc=0;

   if ( T_oldorient_to_neworient[2] < 0.0 || T_oldorient_to_neworient[6] < 0.0 || T_oldorient_to_neworient[10] < 0.0 )
      zc=nz1-1;
   else
      zc=0;

   int q=0;
   int np1 = nx1 * ny1;
   int i1_a, j1_a, k1_a; 
   int i1_b, j1_b, k1_b; 
   int i1, j1, k1; 
   for(int k2=0; k2<nz2; k2++)
   {
      i1_a = (int)T_oldorient_to_neworient[8] * k2 + xc;
      j1_a = (int)T_oldorient_to_neworient[9] * k2 + yc;
      k1_a = (int)T_oldorient_to_neworient[10]* k2 + zc;

      for(int j2=0; j2<ny2; j2++)
      {
         i1_b = (int)T_oldorient_to_neworient[4] * j2 + i1_a;
         j1_b = (int)T_oldorient_to_neworient[5] * j2 + j1_a;
         k1_b = (int)T_oldorient_to_neworient[6] * j2 + k1_a;

         for(int i2=0; i2<nx2; i2++)
         {
            i1 = (int)T_oldorient_to_neworient[0] * i2 + i1_b;
            j1 = (int)T_oldorient_to_neworient[1] * i2 + j1_b;
            k1 = (int)T_oldorient_to_neworient[2] * i2 + k1_b;

            output_image[q++]= input_image[k1*np1 + j1*nx1 + i1];
         }
      }
   }

   newhdr.dim[1]=nx2;
   newhdr.dim[2]=ny2;
   newhdr.dim[3]=nz2;
   newhdr.pixdim[1]=dx2;
   newhdr.pixdim[2]=dy2;
   newhdr.pixdim[3]=dz2;

   update_qsform(newhdr, neworient);

   return(output_image);
}
