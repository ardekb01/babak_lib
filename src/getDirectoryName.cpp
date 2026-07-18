#include <stddef.h>
#include <string.h>

int getDirectoryName(
   const char *pathname,
   char *dirname,
   size_t dirnameSize
)
{
   size_t n;
   size_t i;
   size_t directoryLength;

   if(pathname == NULL || dirname == NULL || dirnameSize == 0)
      return 0;

   n = strlen(pathname);

   for(i = n; i > 0; i--)
   {
      if(pathname[i - 1] == '/')
         break;
   }

   if(i == 0)
   {
      if(dirnameSize < 2)
         return 0;

      dirname[0] = '.';
      dirname[1] = '\0';

      return 1;
   }

   directoryLength = i - 1;

   if(directoryLength + 1 > dirnameSize)
      return 0;

   memcpy(dirname, pathname, directoryLength);

   dirname[directoryLength] = '\0';

   return 1;
}
