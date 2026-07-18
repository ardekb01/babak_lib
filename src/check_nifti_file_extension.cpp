#include <string.h>

// Returns 1 if filename has a .nii extension, 0 otherwise.
// Correctly recognizes:
//    image.nii
// Correctly rejects:
//    image.NII
//    image.nii.gz
//    image.hdr
int check_nifti_file_extension(const char *filename)
{
   size_t filenameLength;
   const char *extension;

   if(filename == nullptr)
      return 0;

   filenameLength = strlen(filename);

   if(filenameLength < 4)
      return 0;

   extension = filename + filenameLength - 4;

   if(strcmp(extension, ".nii") == 0)
   {
      return 1;
   }

   return 0;
}
