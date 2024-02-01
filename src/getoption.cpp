///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// Copyright (C) 2024 Babak A. Ardekani, PhD - All Rights Reserved.  //
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

#define _getoption

#include <stdio.h>
#include <stdlib.h> 
#include <string.h>

//////////////////////////////////////////////////////////////////////////////////////////////
int optind=0;
char *optarg;

struct option
{
	char *name;
	int has_arg;
	int val;
};
//////////////////////////////////////////////////////////////////////////////////////////////

int getoption(int argc, char **argv, struct option *options)
{
	int i;

	for(int j=optind; j<argc; j++)
	if( argv[j][0] == '-' )					// here is an option
	{
		optind = j+1;

		i=0;
		while( options[i].val != 0 ) 		// find the option in options structure
		{
			if( strcmp(options[i].name , argv[j])==0 ) 
			{	
				if( options[i].has_arg && optind >= argc ) 
				{ 
					// option requires an argument but there are no more arguments left
					printf("\nOption %s requires an argument.\n\n",options[i].name);
					return('?');
				}

				if( options[i].has_arg && optind < argc ) { optarg = argv[optind++]; }

				return(options[i].val);
			}
			i++;
		}
		
		printf("\nOption %s not recognized.\n\n",argv[j]);
		return('?');

		break;
	}

	// no more options to be processed
	return(-1);
}
