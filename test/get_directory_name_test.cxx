#include "babak_lib.h"

#include <stdio.h>
#include <string.h>

int main()
{
   char pathname[256];
   char dirname[256];
   int exit_st=0;
   
   snprintf(pathname, sizeof(pathname), "%s", "/usr/home/file.nii");
   get_directory_name(pathname, dirname, sizeof(dirname));
   printf("pathname = \"%s\"  dirname = \"%s\"", pathname,dirname);
   if (strcmp(dirname, "/usr/home") == 0) {
      printf("   (PASSED)\n");
   } else {
      printf("   (FAILED)\n");
      exit_st=1;
   }

   snprintf(pathname, sizeof(pathname), "%s", "./home/file.nii");
   get_directory_name(pathname, dirname, sizeof(dirname));
   printf("pathname = \"%s\"  dirname = \"%s\"", pathname,dirname);
   if (strcmp(dirname, "./home") == 0) {
      printf("   (PASSED)\n");
   } else {
      printf("   (FAILED)\n");
      exit_st=1;
   }

   snprintf(pathname, sizeof(pathname), "%s", "file.nii");
   get_directory_name(pathname, dirname, sizeof(dirname));
   printf("pathname = \"%s\"  dirname = \"%s\"", pathname,dirname);
   if (strcmp(dirname, ".") == 0) {
      printf("   (PASSED)\n");
   } else {
      printf("   (FAILED)\n");
      exit_st=1;
   }

   snprintf(pathname, sizeof(pathname), "%s", "file.nii");
   snprintf(pathname, sizeof(pathname), "%s", "/");
   get_directory_name(pathname, dirname, sizeof(dirname));
   printf("pathname = \"%s\"  dirname = \"%s\"", pathname,dirname);
   if (strcmp(dirname, "/") == 0) {
      printf("   (PASSED)\n");
   } else {
      printf("   (FAILED)\n");
      exit_st=1;
   }

   snprintf(pathname, sizeof(pathname), "%s", "/usr/home/");
   get_directory_name(pathname, dirname, sizeof(dirname));
   printf("pathname = \"%s\"  dirname = \"%s\"", pathname,dirname);
   if (strcmp(dirname, "/usr/home") == 0) {
      printf("   (PASSED)\n");
   } else {
      printf("   (FAILED)\n");
      exit_st=1;
   }

   exit(exit_st);
}
