#include <string.h>
#include <stddef.h>

// Returns 1 if filename has a .nii extension, 0 otherwise.
// Correctly recognizes:
//    image.nii
// Correctly rejects:
//    image.NII
//    image.nii.gz
//    image.hdr
bool check_nifti_file_extension(const char *filename)
{
   size_t filenameLength;
   const char *extension;

   if(filename == nullptr)
      return false;

   filenameLength = strlen(filename);

   if(filenameLength < 4)
      return false;

   extension = filename + filenameLength - 4;

   if(strcmp(extension, ".nii") == 0)
   {
      return true;
   }

   return false;
}
