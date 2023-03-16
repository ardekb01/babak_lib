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

#define _update_qsform

int opt_sform=NO;
int opt_qform=NO;

void update_qsform(nifti_1_header &hdr, float *matrix)
{
   //FILE *fp;
   mat44 R;

   hdr.sform_code = NIFTI_XFORM_ALIGNED_ANAT;
   hdr.srow_x[0]=matrix[0]; hdr.srow_x[1]=matrix[1]; hdr.srow_x[2]=matrix[2]; hdr.srow_x[3]=matrix[3];
   hdr.srow_y[0]=matrix[4]; hdr.srow_y[1]=matrix[5]; hdr.srow_y[2]=matrix[6]; hdr.srow_y[3]=matrix[7];
   hdr.srow_z[0]=matrix[8]; hdr.srow_z[1]=matrix[9]; hdr.srow_z[2]=matrix[10]; hdr.srow_z[3]=matrix[11];

   hdr.qform_code = NIFTI_XFORM_ALIGNED_ANAT;
   R.m[0][0]=matrix[0];  R.m[0][1]=matrix[1];  R.m[0][2]=matrix[2];  R.m[0][3]=matrix[3];
   R.m[1][0]=matrix[4];  R.m[1][1]=matrix[5];  R.m[1][2]=matrix[6];  R.m[1][3]=matrix[7];
   R.m[2][0]=matrix[8];  R.m[2][1]=matrix[9];  R.m[2][2]=matrix[10]; R.m[2][3]=matrix[11];
   R.m[3][0]=matrix[12]; R.m[3][1]=matrix[13]; R.m[3][2]=matrix[14]; R.m[3][3]=matrix[15];

   nifti_mat44_to_quatern( R,  &(hdr.quatern_b), &(hdr.quatern_c), &(hdr.quatern_d),
   &(hdr.qoffset_x), &(hdr.qoffset_y), &(hdr.qoffset_z), 
   &(hdr.pixdim[1]), &(hdr.pixdim[2]), &(hdr.pixdim[3]), &(hdr.pixdim[0]));
}

void update_qsform(const char *imagefilename , float *matrix)
{
   FILE *fp;
   nifti_1_header hdr; // 348 bytes
   nifti1_extender ext; // 4 bytes
   char *extension;
   int extension_size=0;
   char *data;
   int data_size=0;
   char swapflg=0;
   mat44 R;

   fp = fopen(imagefilename,"r");
   if(fp==NULL) file_open_error(imagefilename);
   fread(&hdr, sizeof(nifti_1_header), 1, fp);

   if(hdr.dim[0]<1 || hdr.dim[0]>7)
   {
      swapniftiheader(&hdr);
      swapflg=1;
   }

   if( opt_sform)
   {
      hdr.sform_code = NIFTI_XFORM_ALIGNED_ANAT;
      hdr.srow_x[0]=matrix[0]; hdr.srow_x[1]=matrix[1]; hdr.srow_x[2]=matrix[2]; hdr.srow_x[3]=matrix[3];
      hdr.srow_y[0]=matrix[4]; hdr.srow_y[1]=matrix[5]; hdr.srow_y[2]=matrix[6]; hdr.srow_y[3]=matrix[7];
      hdr.srow_z[0]=matrix[8]; hdr.srow_z[1]=matrix[9]; hdr.srow_z[2]=matrix[10]; hdr.srow_z[3]=matrix[11];
   }

   if( opt_qform)
   {
      hdr.qform_code = NIFTI_XFORM_ALIGNED_ANAT;
      R.m[0][0]=matrix[0];  R.m[0][1]=matrix[1];  R.m[0][2]=matrix[2];  R.m[0][3]=matrix[3];
      R.m[1][0]=matrix[4];  R.m[1][1]=matrix[5];  R.m[1][2]=matrix[6];  R.m[1][3]=matrix[7];
      R.m[2][0]=matrix[8];  R.m[2][1]=matrix[9];  R.m[2][2]=matrix[10]; R.m[2][3]=matrix[11];
      R.m[3][0]=matrix[12]; R.m[3][1]=matrix[13]; R.m[3][2]=matrix[14]; R.m[3][3]=matrix[15];

      nifti_mat44_to_quatern( R,  &(hdr.quatern_b), &(hdr.quatern_c), &(hdr.quatern_d),
      &(hdr.qoffset_x), &(hdr.qoffset_y), &(hdr.qoffset_z), 
      &(hdr.pixdim[1]), &(hdr.pixdim[2]), &(hdr.pixdim[3]), &(hdr.pixdim[0]));
   }

   if( hdr.magic[0]=='n' && hdr.magic[1]=='+' && hdr.magic[2]=='1' )
   {
      fread(&ext, sizeof(nifti1_extender), 1, fp);

      extension_size = (int)(hdr.vox_offset)-352;

      if( extension_size > 0 )
      {
         extension = (char *)calloc(extension_size, 1);
         fread(extension, 1, extension_size, fp);
      }

      data_size = 1;
      for(int i=1; i<=hdr.dim[0]; i++)
      {
         data_size *= hdr.dim[i];
      }
      data_size *= (hdr.bitpix/8);

      if( data_size > 0 )
      {
         data = (char *)calloc(data_size, 1);
         fread(data, 1, data_size, fp);
      }
   }

   fclose(fp);

   fp = fopen(imagefilename,"w");

   if(fp==NULL)
   {
      printf("\nWarning: Could not write to %s\n", imagefilename);
      return;
   }

   if(swapflg)
   {
      swapniftiheader(&hdr);
   }

   fwrite(&hdr, sizeof(nifti_1_header), 1, fp);

   if( hdr.magic[0]=='n' && hdr.magic[1]=='+' && hdr.magic[2]=='1' )
   {
      fwrite(&ext, sizeof(nifti1_extender), 1, fp);
      if( extension_size > 0 )
      {
         fwrite(extension, 1, extension_size, fp);
         delete extension;
      }
      if( data_size > 0 )
      {
         fwrite(data, 1, data_size, fp);
         delete data;
      }
   }
   fclose(fp);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// imagefile1 = source
// imagefile2 = target 

void update_qsform( const char *imagefile1, const char *imagefile2)
{
   FILE *fp;
   nifti_1_header hdr; // 348 bytes
   //nifti1_extender ext; // 4 bytes
   //char *extension;
   //int extension_size=0;
   //char *data;
   //int data_size=0;

   mat44 R;
   float q1[16], q2[16], s1[16], s2[16];
   float *invq1;

   ///////////////////////////////////////////////////////////////////////////
   fp = fopen(imagefile1,"r");
   if(fp==NULL) file_open_error(imagefile1);

   fread(&hdr, sizeof(nifti_1_header), 1, fp);
   if(hdr.dim[0]<1 || hdr.dim[0]>7)
   {
      swapniftiheader(&hdr);
   }

   s1[0]=hdr.srow_x[0]; s1[1]=hdr.srow_x[1]; s1[2]=hdr.srow_x[2]; s1[3]=hdr.srow_x[3];
   s1[4]=hdr.srow_y[0]; s1[5]=hdr.srow_y[1]; s1[6]=hdr.srow_y[2]; s1[7]=hdr.srow_y[3];
   s1[8]=hdr.srow_z[0]; s1[9]=hdr.srow_z[1]; s1[10]=hdr.srow_z[2]; s1[11]=hdr.srow_z[3];
   s1[12]=0.0; s1[13]=0.0; s1[14]=0.0; s1[15]=1.0;

   R=nifti_quatern_to_mat44( hdr.quatern_b, hdr.quatern_c, hdr.quatern_d,
   hdr.qoffset_x, hdr.qoffset_y, hdr.qoffset_z, hdr.pixdim[1], hdr.pixdim[2], hdr.pixdim[3], hdr.pixdim[0]);

   q1[0]=R.m[0][0]; q1[1]=R.m[0][1]; q1[2]=R.m[0][2]; q1[3]=R.m[0][3];
   q1[4]=R.m[1][0]; q1[5]=R.m[1][1]; q1[6]=R.m[1][2]; q1[7]=R.m[1][3];
   q1[8]=R.m[2][0]; q1[9]=R.m[2][1]; q1[10]=R.m[2][2]; q1[11]=R.m[2][3];
   q1[12]=0.0; q1[13]=0.0; q1[14]=0.0; q1[15]=1.0;

   fclose(fp);

   ///////////////////////////////////////////////////////////////////////////

   fp = fopen(imagefile2,"r");
   if(fp==NULL) file_open_error(imagefile2);

   fread(&hdr, sizeof(nifti_1_header), 1, fp);
   if(hdr.dim[0]<1 || hdr.dim[0]>7)
   {
      swapniftiheader(&hdr);
   }

   R=nifti_quatern_to_mat44( hdr.quatern_b, hdr.quatern_c, hdr.quatern_d,
   hdr.qoffset_x, hdr.qoffset_y, hdr.qoffset_z, hdr.pixdim[1], hdr.pixdim[2], hdr.pixdim[3], hdr.pixdim[0]);

   q2[0]=R.m[0][0]; q2[1]=R.m[0][1]; q2[2]=R.m[0][2]; q2[3]=R.m[0][3];
   q2[4]=R.m[1][0]; q2[5]=R.m[1][1]; q2[6]=R.m[1][2]; q2[7]=R.m[1][3];
   q2[8]=R.m[2][0]; q2[9]=R.m[2][1]; q2[10]=R.m[2][2]; q2[11]=R.m[2][3];
   q2[12]=0.0; q2[13]=0.0; q2[14]=0.0; q2[15]=1.0;

   fclose(fp);
   ///////////////////////////////////////////////////////////////////////////

   invq1 = inv4(q1);

   multi(invq1,4,4,q2,4,4,s2);
   multi(s1,4,4,s2,4,4,s2);

   update_qsform( imagefile2, s2);

   free(invq1);
}
