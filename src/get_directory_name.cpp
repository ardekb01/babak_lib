#include <stddef.h>
#include <string.h>

/*
   Extract the directory portion of a pathname.

   Returns:
      true on success
      false on failure

   Examples:
      "/home/user/file.txt" -> "/home/user"
      "file.txt"            -> "."
      "/"                   -> "/"

   The function does not modify pathname.
*/
bool get_directory_name(
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
      return false;

   n = strlen(pathname);

   // Handle the case pathname = "".
   if(n == 0)
      return false;

   // Handle the case pathname = "/".
   if(n == 1 && pathname[0] == '/')
   {
      if(dirnameSize < 2)
         return false;

      dirname[0] = '/';
      dirname[1] = '\0';

      return true;
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
         return false;

      dirname[0] = '.';
      dirname[1] = '\0';

      return true;
   }

   directoryLength = i - 1;

   if(directoryLength >= dirnameSize)
      return false;

   memcpy(dirname, pathname, directoryLength);

   dirname[directoryLength] = '\0';

   return true;
}
