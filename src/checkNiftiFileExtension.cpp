#include <string.h>

// Returns 1 if filename has a .hdr or .nii extension, 0 otherwise.
int checkNiftiFileExtension(const char *filename)
{
   int L; // number of characters in the filename

   L = strlen(filename);

   if(L < 4)
   {
      return(0);
   }

   if(filename[L - 4] == '.' &&
      filename[L - 3] == 'h' &&
      filename[L - 2] == 'd' &&
      filename[L - 1] == 'r')
   {
      return(1);
   }

   if(filename[L - 4] == '.' &&
      filename[L - 3] == 'n' &&
      filename[L - 2] == 'i' &&
      filename[L - 1] == 'i')
   {
      return(1);
   }

   return(0);
}
