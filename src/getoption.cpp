///////////////////////////////////////////////////////////////////////
// Copyright (C) 2024 Babak A. Ardekani, PhD - All Rights Reserved.
///////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

int optind = 1;
char *optarg = NULL;

////////////////////////////////////////////////////////////////////////
// Option structure
////////////////////////////////////////////////////////////////////////

struct CmdOption
{
    const char *name;
    int has_arg;
    int val;
};

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

int getoption(int argc, char **argv, struct CmdOption *options)
{
   int i;

   // Guard against NULL options
   if (options == NULL)
      return '?';

   for (int j = optind; j < argc; j++)
   {
      /* Skip non-option arguments */
      if (argv[j][0] != '-')
         continue;

      optind = j + 1;

      /* Search for the option in the option table */
      for (i = 0; options[i].val != 0; i++)
      {
         if (strcmp(options[i].name, argv[j]) != 0)
            continue;

         /* Option requires an argument */
         if (options[i].has_arg)
         {
            if (optind >= argc)
            {
               printf("\nOption %s requires an argument.\n\n",
                      options[i].name);
               return '?';
            }

            optarg = argv[optind++];
         }

         return options[i].val;
      }

      printf("\nOption %s not recognized.\n\n", argv[j]);
      return '?';
   }

   /* No more options */
   return -1;
}
