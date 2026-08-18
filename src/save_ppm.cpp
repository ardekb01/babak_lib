#include <cstring>
#include <cstdio>
#include "babak_lib.h"

bool save_ppm(const char *filename, const int nx, const int ny,
              const unsigned char *R, 
              const unsigned char *G, 
              const unsigned char *B)
{
   if(filename == nullptr || filename[0] == '\0')
      return false;

   if(nx <= 0 || ny <= 0)
      return false;

   if(R == nullptr || G == nullptr || B == nullptr)
      return false;

  const size_t L = strlen(filename);

   if(L < 4 ||
      filename[L - 4] != '.' ||
      filename[L - 3] != 'p' ||
      filename[L - 2] != 'p' ||
      filename[L - 1] != 'm')
   {
      return false;
   }

   const size_t np = static_cast<size_t>(nx) *
                     static_cast<size_t>(ny);

   FILE *fp = fopen(filename, "wb");

   if(fp == nullptr)
      return false;

   fprintf(fp, "P6\n");
   fprintf(fp, "# Created by the Automatic Registration Toolbox (ART) package.\n");
   fprintf(fp, "%d %d\n", nx, ny);
   fprintf(fp, "255\n");

   unsigned char rgb[3];

   for(size_t i = 0; i < np; ++i)
   {
      rgb[0] = R[i];
      rgb[1] = G[i];
      rgb[2] = B[i];

      if(fwrite(rgb, 1, 3, fp) != 3)
      {
         fclose(fp);
         return false;
      }
   }

   if(fclose(fp) != 0)
      return false;

   if(opt_png)
   {
      char *pngfilename;
      char cmnd[DEFAULT_STRING_LENGTH];
      int system_return_value;

      // Allocate L + 1 bytes to include the terminating null character.
      // Previously, only L bytes were allocated, causing stpcpy() to
      // write one byte beyond the allocated buffer and corrupt the heap.
      pngfilename = static_cast<char *>(calloc(L + 1, sizeof(char)));

      if(pngfilename == nullptr)
         return false;

      stpcpy(pngfilename, filename);
      pngfilename[L - 3] = 'p';
      pngfilename[L - 2] = 'n';
      pngfilename[L - 1] = 'g';

      snprintf(cmnd, sizeof(cmnd), "pnmtopng %s > %s",
               filename, pngfilename);

      system_return_value = system(cmnd);

      if(system_return_value == -1 || system_return_value == 127)
      {
         printf("Warning: %s failed.\n", cmnd);
         printf("Install \'pnmtopng\'.\n");
      }

      free(pngfilename);
   }

   if(opt_ppm == NO)
      remove(filename);

   return true;
}
