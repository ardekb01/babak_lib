#include "check_nifti_file_extension.h"

#include <cstdio>
#include <cstring>

// Extracts the basename from the full path string.
//
// Examples:
//    "/home/babak/images/test.nii" -> "test"
//    "test.nii"                    -> "test"
//    "/home/babak/my.test.nii"     -> "my.test"
//
// The ".nii" extension is not included in the output basename.
//
// Returns false on failure and true on success.
bool get_nifti_basename(char *basename,
                        size_t basenameSize,
                        const char *path)
{
   size_t i;
   size_t len;

   // Validate input variables.
   if (basename == nullptr)
   {
      fprintf(stderr,
              "Error: get_nifti_basename(): NULL pointer passed as basename.\n");
      return false;
   }

   if (path == nullptr)
   {
      fprintf(stderr,
              "Error: get_nifti_basename(): NULL pointer passed as path.\n");
      return false;
   }

   if (basenameSize == 0)
   {
      fprintf(stderr,
              "Error: get_nifti_basename(): 0 passed as basenameSize.\n");
      return false;
   }

   // Ensure that the specified image has a .nii extension.
   if (!check_nifti_file_extension(path))
   {
      fprintf(stderr,
              "%s does not have `.nii' extension.\n",
              path);
      return false;
   }

   // Find the position of the basename, immediately after the last '/'.
   i = strlen(path);

   while (i > 0 && path[i - 1] != '/')
   {
      i--;
   }

   // Determine the length of the final path component, including its extension.
   len = strlen(path + i);

   if (len + 1 > basenameSize)
   {
      fprintf(stderr,
              "Error: get_nifti_basename(): output \"basename\" buffer is too small.\n");
      return false;
   }

   // Copy the basename, including the ".nii" extension.
   strcpy(basename, path + i);

   // Remove the ".nii" extension.
   if (len < 4)
   {
      fprintf(stderr,
              "Error: get_nifti_basename(): unexpected basename length.\n");
      return false;
   }

   basename[len - 4] = '\0';

   // Reject a basename consisting only of ".nii".
   if (basename[0] == '\0')
   {
      fprintf(stderr,
              "Error: get_nifti_basename(): basename has zero length.\n");
      return false;
   }

   return true;
}
