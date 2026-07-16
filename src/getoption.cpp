///////////////////////////////////////////////////////////////////////
// Copyright (C) 2024 Babak A. Ardekani, PhD - All Rights Reserved.
///////////////////////////////////////////////////////////////////////

#include "getoption.h"
#include <stdio.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

int optInd = 1;
const char *optArg = NULL;

////////////////////////////////////////////////////////////////////////
// getoption()
//
// Simple command-line option parser.
//
// Returns:
//      option value   - recognized option
//      '?'            - error
//      -1             - no more options
////////////////////////////////////////////////////////////////////////

int getoption(int argc, char *const argv[], const struct CmdOption *options)
{
   size_t i;

   optArg = NULL;

   if (argv == NULL || options == NULL)
      return '?';

   if (optInd >= argc || argv[optInd] == NULL)
      return -1;

   if (strcmp(argv[optInd], "--") == 0)
   {
      optInd++;
      return -1;
   }

   if (argv[optInd][0] != '-' || argv[optInd][1] == '\0')
      return -1;

   for (i = 0; options[i].val != 0; i++)
   {
      if (options[i].name != NULL && strcmp(options[i].name, argv[optInd]) == 0)
      {
         optInd++;

         if (options[i].has_arg)
         {
            if (optInd >= argc)
            {
               fprintf(stderr,
                       "Option %s requires an argument.\n",
                       options[i].name);
               return '?';
            }

            optArg = argv[optInd++];
         }

         return options[i].val;
      }
   }

   fprintf(stderr,
           "Option %s not recognized.\n",
           argv[optInd]);

   optInd++;

   return '?';
}
