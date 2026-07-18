#include <stddef.h>
#include <string.h>

/*
   Extract the directory portion of a pathname.

   Returns:
      1 on success
      0 on failure

   Examples:
      "/home/user/file.txt" -> "/home/user"
      "file.txt"            -> "."
      "/"                   -> "/"

   The function does not modify pathname.
*/
int get_directory_name(
   const char *pathname,
   char *dirname,
   size_t dirnameSize
)
{
   size_t n;
   size_t i;
   size_t directoryLength;

   // Validate input arguments.
   if(pathname == nullptr || dirname == nullptr || dirnameSize == 0)
      return 0;

   n = strlen(pathname);

   // Handle the case pathname = "".
   if(n == 0)
      return 0;

   // Handle the case pathname = "/".
   if(n == 1 && pathname[0] == '/')
   {
      if(dirnameSize < 2)
         return 0;

      dirname[0] = '/';
      dirname[1] = '\0';

      return 1;
   }

   for(i = n; i > 0; i--)
   {
      if(pathname[i - 1] == '/')
         break;
   }

   // If no directory separator is present, for example, when
   // pathname = "file.nii", use the current directory.
   if(i == 0)
   {
      if(dirnameSize < 2)
         return 0;

      dirname[0] = '.';
      dirname[1] = '\0';

      return 1;
   }

   directoryLength = i - 1;

   if(directoryLength >= dirnameSize)
      return 0;

   memcpy(dirname, pathname, directoryLength);

   dirname[directoryLength] = '\0';

   return 1;
}
