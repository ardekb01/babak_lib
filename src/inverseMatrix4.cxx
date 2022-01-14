#include <stdlib.h>
//#include <malloc.h>  
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include "babak_lib.h"

void print_matrix(char * title, float *T)
{
    printf("\n%s:",title);
    printf("\n%f\t%f\t%f\t%f",T[0],T[1],T[2],T[3]);
    printf("\n%f\t%f\t%f\t%f",T[4],T[5],T[6],T[7]);
    printf("\n%f\t%f\t%f\t%f",T[8],T[9],T[10],T[11]);
    printf("\n%f\t%f\t%f\t%f",T[12],T[13],T[14],T[15]);
    printf("\n");
}

int main(int argc, char **argv)
{
	float Tin[16];
	float *Tout;
	FILE *fp;

  if(argc==1) exit(0);

	loadTransformation(argv[1],Tin);

	Tout = inv4(Tin);

	fp = fopen(argv[2],"w");
    fprintf(fp,"%f\t%f\t%f\t%f\n",Tout[0],Tout[1],Tout[2],Tout[3]);
    fprintf(fp,"%f\t%f\t%f\t%f\n",Tout[4],Tout[5],Tout[6],Tout[7]);
    fprintf(fp,"%f\t%f\t%f\t%f\n",Tout[8],Tout[9],Tout[10],Tout[11]);
    fprintf(fp,"%f\t%f\t%f\t%f\n",Tout[12],Tout[13],Tout[14],Tout[15]);
	fclose(fp);

//	float I[16];
//	multi(Tin, 4,4, Tout, 4,4, I);
//	print_matrix("Output",I);
}
