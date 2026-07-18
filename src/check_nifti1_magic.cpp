#include "nifti1.h"

#include <stdio.h>

// Checks the last 4 bytes of the header
// for NIFTI-1 magic string "n+1" 
int check_nifti1_magic(const char *imagefilename)
{
   FILE *fp;
   nifti_1_header hdr;

   fp = fopen(imagefilename,"r");
   
   if(fp == NULL)
   {
      return 0;
   }

   if(fread(&hdr, sizeof(hdr), 1, fp) != 1)
   {
      fclose(fp);
      return 0;
   }

   fclose(fp);

   if( hdr.magic[0]!='n' || hdr.magic[1]!='+' || hdr.magic[2]!='1' )
   {
      return 0;
   }

   return 1;
}
