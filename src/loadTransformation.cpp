///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// Copyright (C) 2026 Babak A. Ardekani, PhD - All Rights Reserved.  //
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

#include "loadTransformation.h"
#include <ctype.h>
#include <stdio.h>

// Loads a 4x4 transformation matrix from a text file.
// The file is expected to contain four lines of four floating-point numbers.
// Lines beginning with '#' or blank lines are ignored.
// If no filename is provided, returns the identity transformation.
//
// Returns:
//    LOADTRANSFORM_OK            : success
//    LOADTRANSFORM_NULL_POINTER  : T is NULL
//    LOADTRANSFORM_FILE_ERROR    : file cannot be opened
//    LOADTRANSFORM_PARSE_ERROR   : invalid or incomplete matrix

int loadTransformation(const char *filename, float *T)
{
   FILE *fp;
   char line[256];
   int i;

   // Verify that the output matrix pointer is valid.
   if(T == NULL)
   {
      return LOADTRANSFORM_NULL_POINTER;
   }

   // If no filename is provided, return the identity transformation.
   if(filename == NULL || filename[0] == '\0')
   {
      for(int i=0; i<16; i++)
         T[i] = 0.0f;

      T[0] = T[5] = T[10] = T[15] = 1.0f;

      return LOADTRANSFORM_OK;
   }

   // Open the transformation file.
   fp = fopen(filename, "r");

   if(fp == NULL)
   {
      return LOADTRANSFORM_FILE_ERROR;
   }

   // Read up to four rows of the transformation matrix.
   // Skip comments and blank lines.
   i = 0;
   while (i < 4 && fgets(line, sizeof(line), fp) != NULL)
   {
      char *p = line;
      while (isspace((unsigned char)*p))
         ++p;

      if (*p == '\0' || *p == '#')
         continue;

      float *row = &T[4*i];
      if (sscanf(line, "%f %f %f %f",
         &row[0], &row[1], &row[2], &row[3]) != 4)
      {
         fclose(fp);
         return LOADTRANSFORM_PARSE_ERROR;
      }
      i++;
   }

   fclose(fp);

   // Verify that all four rows were successfully read.
   if(i != 4)
   {
      return LOADTRANSFORM_PARSE_ERROR;
   }

   return LOADTRANSFORM_OK;
}
