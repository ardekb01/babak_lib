#include "check_nifti_file_extension.h"

#include <cstdio>
#include <cstring>

// Extracts the filename from the full path string.
//
// Examples:
//    "/home/babak/images/test.nii" -> "test"
//    "test.nii"                    -> "test"
//    "/home/babak/my.test.nii"     -> "my.test"
//
// The ".nii" extension is not included in the output filename.
//
// Returns false on failure and true on success.
bool get_nifti_filename(char *filename,
                        size_t filenameSize,
                        const char *path)
{
   size_t i;
   size_t len;
   size_t pos;

   // Validate input variables.
   if (filename == nullptr)
   {
      fprintf(stderr,
              "Error: get_nifti_filename(): NULL pointer passed as filename.\n");
      return false;
   }

   if (path == nullptr)
   {
      fprintf(stderr,
              "Error: get_nifti_filename(): NULL pointer passed as path.\n");
      return false;
   }

   if (filenameSize == 0)
   {
      fprintf(stderr,
              "Error: get_nifti_filename(): 0 passed as filenameSize.\n");
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

   pos = i;

   // Determine the length of the filename, including its extension.
   len = strlen(path + pos);

   if (len + 1 > filenameSize)
   {
      fprintf(stderr,
              "Error: get_nifti_filename(): output \"filename\" buffer is too small.\n");
      return false;
   }

   // Copy the basename, including the ".nii" extension.
   strcpy(filename, path + pos);

   // Remove the ".nii" extension.
   if (len < 4)
   {
      fprintf(stderr,
              "Error: get_nifti_filename(): unexpected filename length.\n");
      return false;
   }

   filename[len - 4] = '\0';

   // Reject a filename consisting only of ".nii".
   if (filename[0] == '\0')
   {
      fprintf(stderr,
              "Error: get_nifti_filename(): filename has zero length.\n");
      return false;
   }

   return true;
}
