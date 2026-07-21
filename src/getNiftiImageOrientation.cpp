///////////////////////////////////////////////////////////////////////
// Copyright (C) 2024 Babak A. Ardekani, PhD - All Rights Reserved.
///////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include "nifti1.h"
#include "nifti1_io.h"
#include "babak_lib.h"
#include "swap.h"
#include "directionCode.h"

bool getNiftiImageOrientation(const char *filename, char *orientation)
{
   FILE *fp;
   nifti_1_header hdr;

   if (filename == nullptr || orientation == nullptr)
      return false;

   orientation[0] = '\0';

   fp = fopen(filename, "rb");

   if (fp == nullptr)
      return false;

   if (fread(&hdr, sizeof(nifti_1_header), 1, fp) != 1)
   {
      fclose(fp);
      return false;
   }

   fclose(fp);

   // If dim[0] is outside the range 1..7, then the header information
   // needs to be byte swapped appropriately.
   if (hdr.dim[0] < 1 || hdr.dim[0] > 7)
   {
      swapByteOrder((char *)&(hdr.qform_code), sizeof(short));
      swapByteOrder((char *)&(hdr.sform_code), sizeof(short));
      swapByteOrder((char *)&(hdr.quatern_b), sizeof(float));
      swapByteOrder((char *)&(hdr.quatern_c), sizeof(float));
      swapByteOrder((char *)&(hdr.quatern_d), sizeof(float));
      swapByteOrder((char *)&(hdr.qoffset_x), sizeof(float));
      swapByteOrder((char *)&(hdr.qoffset_y), sizeof(float));
      swapByteOrder((char *)&(hdr.qoffset_z), sizeof(float));

      for (int i = 0; i < 4; i++)
      {
         swapByteOrder((char *)&(hdr.srow_x[i]), sizeof(float));
         swapByteOrder((char *)&(hdr.srow_y[i]), sizeof(float));
         swapByteOrder((char *)&(hdr.srow_z[i]), sizeof(float));
      }

      swap_float_array(hdr.pixdim, 8);
   }

   if (hdr.qform_code == 0 && hdr.sform_code == 0)
   {
      // It is not possible to determine image orientation.
      return false;
   }

   // Prefer qform over sform when both are defined.
   if (hdr.qform_code > 0)
   {
      mat44 R;

      R = nifti_quatern_to_mat44(
         hdr.quatern_b,
         hdr.quatern_c,
         hdr.quatern_d,
         hdr.qoffset_x,
         hdr.qoffset_y,
         hdr.qoffset_z,
         hdr.pixdim[1],
         hdr.pixdim[2],
         hdr.pixdim[3],
         hdr.pixdim[0]
      );

      orientation[0] = directionCode(R.m[0][0], R.m[1][0], R.m[2][0]);
      orientation[1] = directionCode(R.m[0][1], R.m[1][1], R.m[2][1]);
      orientation[2] = directionCode(R.m[0][2], R.m[1][2], R.m[2][2]);
   }
   else
   {
      orientation[0] = directionCode(
         hdr.srow_x[0],
         hdr.srow_y[0],
         hdr.srow_z[0]
      );

      orientation[1] = directionCode(
         hdr.srow_x[1],
         hdr.srow_y[1],
         hdr.srow_z[1]
      );

      orientation[2] = directionCode(
         hdr.srow_x[2],
         hdr.srow_y[2],
         hdr.srow_z[2]
      );
   }

   orientation[3] = '\0';

   return true;
}

bool getNiftiImageOrientation(nifti_1_header hdr,
                              char *orientation)
{
   if (orientation == nullptr)
      return false;

   orientation[0] = '\0';

   if (hdr.qform_code == 0 && hdr.sform_code == 0)
   {
      // The header does not contain orientation information.
      return false;
   }

   // Prefer qform over sform when both are defined.
   if (hdr.qform_code > 0)
   {
      mat44 R;

      R = nifti_quatern_to_mat44(hdr.quatern_b,
                                 hdr.quatern_c,
                                 hdr.quatern_d,
                                 hdr.qoffset_x,
                                 hdr.qoffset_y,
                                 hdr.qoffset_z,
                                 hdr.pixdim[1],
                                 hdr.pixdim[2],
                                 hdr.pixdim[3],
                                 hdr.pixdim[0]);

      orientation[0] = directionCode(R.m[0][0],
                                     R.m[1][0],
                                     R.m[2][0]);

      orientation[1] = directionCode(R.m[0][1],
                                     R.m[1][1],
                                     R.m[2][1]);

      orientation[2] = directionCode(R.m[0][2],
                                     R.m[1][2],
                                     R.m[2][2]);

   }
   else
   {
      orientation[0] = directionCode(hdr.srow_x[0],
                                     hdr.srow_y[0],
                                     hdr.srow_z[0]);

      orientation[1] = directionCode(hdr.srow_x[1],
                                     hdr.srow_y[1],
                                     hdr.srow_z[1]);

      orientation[2] = directionCode(hdr.srow_x[2],
                                     hdr.srow_y[2],
                                     hdr.srow_z[2]);

   }
   orientation[3] = '\0';

   return true;
}
