#include "nifti1.h"
#include <cstdio>

// Checks the NIFTI-1 magic string "n+1".
bool check_nifti1_magic(const char *imagefilename)
{
   FILE *fp;
   nifti_1_header hdr;

   if (imagefilename == nullptr)
   {
      return false;
   }

   fp = fopen(imagefilename, "rb");

   if (fp == nullptr)
   {
      return false;
   }

   if (fread(&hdr, sizeof(hdr), 1, fp) != 1)
   {
      fclose(fp);
      return false;
   }

   fclose(fp);

   if (hdr.magic[0] == 'n' &&
       hdr.magic[1] == '+' &&
       hdr.magic[2] == '1')
   {
      return true;
   }

   return false;
}
