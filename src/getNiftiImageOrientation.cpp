///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// Copyright (C) 2024 Babak A. Ardekani, PhD - All Rights Reserved.  //
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <nifti1.h>
#include <nifti1_io.h>
#include <babak_lib.h>

void getNiftiImageOrientation(const char *filename, char *orientation)
{
   FILE *fp;
   nifti_1_header hdr;
   int swapflg=0;

   orientation[0]='\0';

   // ensure that the specified image has either a .hdr or a .nii extension
   if( !checkNiftiFileExtension(filename) )
   {
      errorMessage("The image filename must have a `.hdr' or `.nii' extension.");
   }

   fp = fopen(filename,"r");

   if(fp==NULL)
   {
      file_open_error(filename);
   }

   if( fread(&hdr, sizeof(nifti_1_header), 1, fp) != 1 )
   {
      errorMessage("I have trouble reading the specified image file.");
   }

   fclose(fp);

   // looks like ANALYZE 7.5, cannot determine orientation
   if( hdr.magic[0]!='n' ||  (hdr.magic[1]!='+' && hdr.magic[1]!='i') ||  hdr.magic[2]!='1')
   {
      printf("\nWarning: Could not determine %s image orientation ...\n",filename);
      return;
   }

   // if dim[0] is outside range 1..7, then the header information
   // needs to be byte swapped appropriately
   if(hdr.dim[0]<1 || hdr.dim[0]>7) 
   {
      swapflg=1;
   }

   // Here I am only byte swapping the header fields relevant for determining the image orientation
   if(swapflg)
   {
      swapByteOrder( (char *)&(hdr.qform_code), sizeof(short));
      swapByteOrder( (char *)&(hdr.sform_code), sizeof(short));
      swapByteOrder( (char *)&(hdr.quatern_b), sizeof(float));
      swapByteOrder( (char *)&(hdr.quatern_c), sizeof(float));
      swapByteOrder( (char *)&(hdr.quatern_d), sizeof(float));
      swapByteOrder( (char *)&(hdr.qoffset_x), sizeof(float));
      swapByteOrder( (char *)&(hdr.qoffset_y), sizeof(float));
      swapByteOrder( (char *)&(hdr.qoffset_z), sizeof(float));
      for(int i=0; i<4; i++)
      {
         swapByteOrder( (char *)&(hdr.srow_x[i]), sizeof(float));
         swapByteOrder( (char *)&(hdr.srow_y[i]), sizeof(float));
         swapByteOrder( (char *)&(hdr.srow_z[i]), sizeof(float));
      }
   }

   if(hdr.qform_code == 0 && hdr.sform_code == 0) 
   {
      printf("\nWarning: Could not determine %s image orientation ...\n",filename);
      printf("\nWarning: The header of this \"NIFTI\" file does not contain orientation information.\n");
      return;
   }

   if(hdr.qform_code > 0 )
   {
      mat44 R;

      R = nifti_quatern_to_mat44(hdr.quatern_b, hdr.quatern_c, hdr.quatern_d, hdr.qoffset_x, hdr.qoffset_y, 
      hdr.qoffset_z, hdr.pixdim[1], hdr.pixdim[2], hdr.pixdim[3], hdr.pixdim[0]);

      orientation[0] = directionCode(R.m[0][0],R.m[1][0],R.m[2][0]);
      orientation[1] = directionCode(R.m[0][1],R.m[1][1],R.m[2][1]);
      orientation[2] = directionCode(R.m[0][2],R.m[1][2],R.m[2][2]);
      orientation[3] = '\0';
   }
   else
   {
      orientation[0] = directionCode(hdr.srow_x[0],hdr.srow_y[0],hdr.srow_z[0]);
      orientation[1] = directionCode(hdr.srow_x[1],hdr.srow_y[1],hdr.srow_z[1]);
      orientation[2] = directionCode(hdr.srow_x[2],hdr.srow_y[2],hdr.srow_z[2]);
      orientation[3] = '\0';
   }

   return;
}

void getNiftiImageOrientation(nifti_1_header hdr, char *orientation)
{
   orientation[0]='\0';

   if(hdr.qform_code == 0 && hdr.sform_code == 0) 
   {
      printf("\nWarning: Could not determine image orientation ...\n");
      printf("\nWarning: \"NIFTI\" header does not contain orientation information.\n");
      return;
   }

   if(hdr.qform_code > 0 )
   {
      mat44 R;

      R = nifti_quatern_to_mat44(hdr.quatern_b, hdr.quatern_c, hdr.quatern_d, hdr.qoffset_x, hdr.qoffset_y, 
      hdr.qoffset_z, hdr.pixdim[1], hdr.pixdim[2], hdr.pixdim[3], hdr.pixdim[0]);

      orientation[0] = directionCode(R.m[0][0],R.m[1][0],R.m[2][0]);
      orientation[1] = directionCode(R.m[0][1],R.m[1][1],R.m[2][1]);
      orientation[2] = directionCode(R.m[0][2],R.m[1][2],R.m[2][2]);
      orientation[3] = '\0';
   }
   else
   {
      orientation[0] = directionCode(hdr.srow_x[0],hdr.srow_y[0],hdr.srow_z[0]);
      orientation[1] = directionCode(hdr.srow_x[1],hdr.srow_y[1],hdr.srow_z[1]);
      orientation[2] = directionCode(hdr.srow_x[2],hdr.srow_y[2],hdr.srow_z[2]);
      orientation[3] = '\0';
   }

   return;
}
