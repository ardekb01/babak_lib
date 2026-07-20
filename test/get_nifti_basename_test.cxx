#include "babak_lib.h"

#include <stdio.h>
#include <string.h>

int main()
{
   char path[256]="";
   char filename[256]="";
   int exit_st=0;
   
   snprintf(path, sizeof(path), "%s", "/usr/home/image.nii");
   get_nifti_basename(filename, sizeof(filename), path);
   printf("path = \"%s\"  filename = \"%s\"", path,filename);
   if (strcmp(filename, "image") == 0) {
      printf("   (PASSED)\n");
   } else {
      printf("   (FAILED)\n");
      exit_st=1;
   }
   
   snprintf(path, sizeof(path), "%s", "/home/file.nii");
   get_nifti_basename(filename, sizeof(filename), path);
   printf("path = \"%s\"  filename = \"%s\"", path,filename);
   if (strcmp(filename, "file") == 0) {
      printf("   (PASSED)\n");
   } else {
      printf("   (FAILED)\n");
      exit_st=1;
   }
   
   snprintf(path, sizeof(path), "%s", "file.nii");
   get_nifti_basename(filename, sizeof(filename), path);
   printf("path = \"%s\"  filename = \"%s\"", path,filename);
   if (strcmp(filename, "file") == 0) {
      printf("   (PASSED)\n");
   } else {
      printf("   (FAILED)\n");
      exit_st=1;
   }

   snprintf(path, sizeof(path), "%s", "//file.nii");
   get_nifti_basename(filename, sizeof(filename), path);
   printf("path = \"%s\"  filename = \"%s\"", path,filename);
   if (strcmp(filename, "file") == 0) {
      printf("   (PASSED)\n");
   } else {
      printf("   (FAILED)\n");
      exit_st=1;
   }

   snprintf(path, sizeof(path), "%s", ".nii");
   get_nifti_basename(filename, sizeof(filename), path);
   printf("path = \"%s\"  filename = \"%s\"", path,filename);
   if (strcmp(filename, "") == 0) {
      printf("   (PASSED)\n");
   } else {
      printf("   (FAILED)\n");
      exit_st=1;
   }

   snprintf(path, sizeof(path), "%s", "/usr/home/my.image.nii");
   get_nifti_basename(filename, sizeof(filename), path);
   printf("path = \"%s\"  filename = \"%s\"", path,filename);
   if (strcmp(filename, "my.image") == 0) {
      printf("   (PASSED)\n");
   } else {
      printf("   (FAILED)\n");
      exit_st=1;
   }
   
   exit(exit_st);
}
