#include <stdio.h>
#include <stdlib.h>
#include <climits>
#include <math.h>
#include <stdint.h>

#include "babak_lib.h"
#include "nifti1.h"
#include "swap.h"

char *read_nifti_image(const char *filename, nifti_1_header *hdr)
{
   FILE *fp;
   int swapflg = 0;
   size_t datasize = 1;
   char *im = NULL;
   long voxeloffset;
   size_t nv = 1;
   float value;

   // Validate input arguments.
   if (filename == nullptr || hdr == nullptr)
      return nullptr;

   // Ensure that the file has '.nii' extension.
   if (!check_nifti_file_extension(filename))
      return nullptr;

   // Ensure that the file has the correct NIFTI-1 'magic' code.
   if (!check_nifti1_magic(filename))
      return nullptr;

   fp = fopen(filename, "rb");

   if (fp == nullptr)
      return nullptr;

   if (fread(hdr, sizeof(nifti_1_header), 1, fp) != 1)
   {
      fclose(fp);
      return nullptr;
   }

   fclose(fp);

   // If dim[0] is outside the range [1,7], then the header information
   // needs to be byte swapped appropriately.
   if (hdr->dim[0] < 1 || hdr->dim[0] > 7)
   {
      swapflg = 1;
      swapniftiheader(hdr);
   }

   if (!isfinite(hdr->vox_offset) ||
       hdr->vox_offset < 0.0f ||
       hdr->vox_offset > 2147483647.0f)
   {
      return nullptr;
   }

   voxeloffset = (long)hdr->vox_offset;

   float slope = hdr->scl_slope;
   float inter = hdr->scl_inter;

   if (slope == 0.0f)
   {
      slope = 1.0f;
      inter = 0.0f;
   }

   if (hdr->dim[0] < 1 || hdr->dim[0] > 7)
      return nullptr;

   if (hdr->bitpix <= 0 || hdr->bitpix % 8 != 0)
      return nullptr;

   for (int i = 1; i <= hdr->dim[0]; i++)
   {
      if (hdr->dim[i] <= 0)
         return nullptr;

      if (nv > SIZE_MAX / (size_t)hdr->dim[i])
         return nullptr;

      nv *= (size_t)hdr->dim[i];
   }

   size_t bytesPerVoxel = (size_t)hdr->bitpix / 8;

   if (nv > SIZE_MAX / bytesPerVoxel)
      return nullptr;

   datasize = nv * bytesPerVoxel;

   fp = fopen(filename, "rb");

   if (fp == nullptr)
      return nullptr;

   if (fseek(fp, voxeloffset, SEEK_SET) != 0)
   {
      fclose(fp);
      return nullptr;
   }

   im = (char *)calloc(datasize, 1);

   if (im == nullptr)
   {
      fclose(fp);
      return nullptr;
   }

   if ( fread(im, 1, datasize, fp) != datasize )
   {
      free(im);
      fclose(fp);
      return nullptr;
   }

   fclose(fp);

   // If necessary swap bytes of the image data. 
   if (swapflg)
   {
      if (hdr->datatype == DT_SIGNED_SHORT ||
          hdr->datatype == DT_UINT16)
      {
         swapN(im, datasize);
      }

      if (hdr->datatype == DT_FLOAT)
      {
         swap_float_array(
            (float *)im,
            datasize / sizeof(float)
         );
      }

      if (hdr->datatype == DT_DOUBLE)
      {
         swap_double_array(
            (double *)im,
            datasize / sizeof(double)
         );
      }

      if (hdr->datatype == DT_SIGNED_INT)
      {
         swap_int_array(
            (int *)im,
            datasize / sizeof(int)
         );
      }
   }

   if (hdr->datatype == DT_SIGNED_SHORT)
   {
      short *tmp;

      tmp = (short *)im;

      for (size_t i = 0; i < nv; i++)
      {
         value = roundf( tmp[i] * slope + inter );

         if (value < SHRT_MIN)
            value = SHRT_MIN;

         if (value > SHRT_MAX)
            value = SHRT_MAX;

         tmp[i] = (short)value;
      }
   }
   else if (hdr->datatype == DT_UINT16)
   {
      unsigned short *tmp;

      tmp = (unsigned short *)im;

      for (size_t i = 0; i < nv; i++)
      {
         value = roundf( tmp[i] * slope + inter );

         if (value < 0.0f)
            value = 0.0f;

         if (value > USHRT_MAX)
            value = USHRT_MAX;

         tmp[i] = (unsigned short)value; 
      }
   }

   return(im);
}
