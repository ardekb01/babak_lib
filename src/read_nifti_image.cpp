#include <limits>
#include <cmath>
#include <cfloat>

#include "babak_lib.h"
#include "nifti1.h"
#include "swap.h"
#include "minmax.h"

bool scaleFloatToUnsignedShort(
   const float *input,
   unsigned short *output,
   size_t n)
{
   float scale;

   if(input == nullptr ||
      output == nullptr ||
      n == 0)
   {
      return false;
   }

   float minValue = input[0];
   float maxValue = input[0];

   minmax(input, n, minValue, maxValue);

   // All values are identical.
   if(fabsf(maxValue - minValue) <= FLT_EPSILON)
   {
      for(size_t i = 0; i < n; i++)
      {
         output[i] = 0;
      }

      return true;
   }

   const float shortMax = (float)std::numeric_limits<short>::max();

   if( maxValue > shortMax )
   {
      scale  = maxValue/shortMax;
      for(size_t i = 0; i < n; i++)
      {
         float s = input[i] / scale;

         if(s > shortMax) s = shortMax;

         output[i] = (unsigned short)lroundf(s);
      }
   }

   return true;
}

bool scaleFloatToShort(
   const float *input,
   short *output,
   size_t n)
{
   float scale;

   if(input == nullptr ||
      output == nullptr ||
      n == 0)
   {
      return false;
   }

   float minValue = input[0];
   float maxValue = input[0];

   minmax(input, n, minValue, maxValue);

   // All values are identical.
   if(fabsf(maxValue - minValue) <= FLT_EPSILON)
   {
      for(size_t i = 0; i < n; i++)
      {
         output[i] = 0;
      }

      return true;
   }

   const float shortMin = (float)std::numeric_limits<short>::min();
   const float shortMax = (float)std::numeric_limits<short>::max();

   if( maxValue > shortMax )
   {
      scale  = maxValue/shortMax;
      for(size_t i = 0; i < n; i++)
      {
         float s = input[i] / scale;

         if(s > shortMax) s = shortMax;

         output[i] = (short)lroundf(s);
      }

      minValue /= scale;
   }

   if( minValue < shortMin )
   for(size_t i = 0; i < n; i++)
   {
      scale  = minValue/shortMin;

      float s = input[i] / scale;

      if(s < shortMin) s = shortMin;

      output[i] = (short)lroundf(s);
   }

   return true;
}

char *read_nifti_image(const char *filename, nifti_1_header *hdr)
{
   FILE *fp;
   int swapflg = 0;
   size_t datasize = 1;
   char *im = NULL;
   long voxeloffset;
   size_t nv = 1;
   float *float_im;

   const float shortMin = (float)std::numeric_limits<short>::min();
   const float shortMax = (float)std::numeric_limits<short>::max();

   const float ushortMax = (float)std::numeric_limits<unsigned short>::max();

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

   // For a NIFTI-1 header, sizeof_hdr should be 348, so an unexpected value 
   // is a reasonable indication that byte swapping is required.
   // If sizeof_hdr is not 348, assume that the header has the
   // opposite byte order and swap the header.
   if(hdr->sizeof_hdr != 348)
   {
      swapflg = 1;
      swapniftiheader(hdr);

      // After swapping, verify that the header is now valid. Otherwise, 
      // a corrupt or non-NIFTI file could be accepted:
      if(hdr->sizeof_hdr != 348)
      {
         return nullptr;
      }
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

   float_im = (float *)calloc(nv, sizeof(float) );

   if (float_im == nullptr)
   {
      free(im);
      fclose(fp);
      return nullptr;
   }

   if ( fread(im, 1, datasize, fp) != datasize )
   {
      free(im);
      free(float_im);
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

   float min, max;

   if (hdr->datatype == DT_SIGNED_SHORT)
   {
      short *tmp;

      tmp = (short *)im;

      for (size_t i = 0; i < nv; i++)
      {
         float_im[i] = roundf( tmp[i] * slope + inter );
      }

      minmax(float_im, nv, min, max);
      if (min < shortMin || max > shortMax )
      {
         scaleFloatToShort(float_im, tmp, nv);
      }
      else
      {
         for (size_t i = 0; i < nv; i++)
            tmp[i] = (short)float_im[i];
      }
   }
   else if (hdr->datatype == DT_UINT16)
   {
      unsigned short *tmp;

      tmp = (unsigned short *)im;

      for (size_t i = 0; i < nv; i++)
      {
         float_im[i] = roundf( tmp[i] * slope + inter );
      }

      minmax(float_im, nv, min, max);
      if (max > ushortMax )
      {
         scaleFloatToUnsignedShort(float_im, tmp, nv);
      }
      else
      {
         for (size_t i = 0; i < nv; i++)
            tmp[i] = (unsigned short)float_im[i];
      }
   }

   free(float_im);
   return(im);
}

// Returns false on failure.
bool read_nifti_hdr(const char *filename, nifti_1_header *hdr)
{
   // Validate inputs.
   if(filename == nullptr ||
      hdr == nullptr)
   {
      return false;
   }

   FILE *fp = fopen(filename, "rb");

   if(fp == nullptr)
   {
      return false;
   }

   if(fread(hdr, sizeof(nifti_1_header), 1, fp) != 1)
   {
      fclose(fp);
      return false;
   }

   fclose(fp);

   // For a NIFTI-1 header, sizeof_hdr should be 348, so an unexpected value 
   // is a reasonable indication that byte swapping is required.
   // If sizeof_hdr is not 348, assume that the header has the
   // opposite byte order and swap the header.
   if(hdr->sizeof_hdr != 348)
   {
      swapniftiheader(hdr);

      // After swapping, verify that the header is now valid. Otherwise, 
      // a corrupt or non-NIFTI file could be accepted:
      if(hdr->sizeof_hdr != 348)
      {
         return false;
      }
   }

   return true;
}
