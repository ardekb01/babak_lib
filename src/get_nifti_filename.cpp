#include "check_nifti_file_extension.h"
#include "check_nifti1_magic.h"

#include <cstdio>
#include <cstring>

// extracts the "filename" from the full "path" string
// Example: If path="/home/babak/images/test.nii", then filename="test".
// Note that the extension is not included in the output filename
// Returns 0 on failture and 1 on success
bool get_nifti_filename(char *filename,
                        size_t filenameSize,
                        const char *path)
{
   int i;
   size_t pathsize;  // length of the path string
   size_t len;
   int pos;  // position of the filename

   if (filename == nullptr || path == nullptr || filenameSize == 0)
   {
      return false;
   }

   // ensure that the specified image has .nii extension
   if (check_nifti_file_extension(path) == false)
   {
      printf("%s does not have `.nii' extension\n", path);
      return false;
   }

   if (check_nifti1_magic(path) == false)
   {
      return false;
   }

   pathsize = (int)strlen(path);

   if (pathsize <= 0)
   {
      printf("Error: unexpected string length for the NIFTI image path, aborting ...\n");
      return false;
   }

   // finds the position of the first '/' character from right if any
   // returns 0 of no '/' charater is there
   // Examples: path=/sss/yyy would give pos=5
   // path=sss would give pos=0
   // path=/x/ would give pos=3
   i = pathsize - 1;

   while (i >= 0 && path[i] != '/')
   {
      i--;
   }

   pos = i + 1;

   if (strlen(path + pos) + 1 > filenameSize)
   {
      fprintf(stderr,"Error: output filename buffer is too small, aborting ...\n");
      return false;
   }

   // copy everything to the right of the first '/' character from right
   // into filename
   // Examples: path=/sss/yyy would give filename=yyy  pathsize=3
   // path=sss would give filename=sss  pathsize=3
   // path=/x/ would give filename=""  pathsize=0
   //strcpy(filename, path + pos);
   memcpy(filename, path + pos, strlen(path + pos) + 1);

   len = (int)strlen(filename);

   if (len <= 0)
   {
      printf("Error: unexpected string length for the NIFTI image filename, aborting ...\n");
      return false;
   }

   // finds the position of the first '.' character from right if any
   // returns 0 of no '.' charater is there
   // Examples: path=/sss.yyy would give pos=5
   // path=sss would give pos=0
   // path=/x. would give pos=3
   i = len - 1;

   while (i >= 0 && filename[i] != '.')
   {
      i--;
   }

   pos = i + 1;

   if (pos > 0)
   {
      filename[pos - 1] = '\0';
   }

   len = (int)strlen(filename);

   if (len <= 0)
   {
      printf("Error: unexpected string length for the NIFTI image prefix, aborting ...\n");
      return false;
   }

   return true;
}
