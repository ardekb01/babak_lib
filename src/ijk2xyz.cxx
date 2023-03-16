#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>  //      required by time(), stat()
#include <time.h>       //      required by time()
#include <sys/stat.h>   //      required by stat() 
#include <unistd.h>
#include <spm_analyze.h>
#include <babak_lib.h>
#include <sph.h>
#include <landmarks.h>
#include <minmax.h>
#include <ctype.h>

#define YES 1
#define NO 0

int opt;

static struct option options[] =
{
   {"-T", 1, 'T'},
   {"-i", 1, 'i'},
   {"-o", 1, 'o'},
   {"-v", 0, 'v'},
   {"-h", 0, 'h'},
   {"-Sdx",1,'1'},
   {"-Sdy",1,'2'},
   {"-Sdz",1,'3'},
   {"-Snx",1,'4'},
   {"-Sny",1,'5'},
   {"-Snz",1,'6'},
   {"-help", 0, 'h'},
   {0, 0, 0}
};

void print_help_and_exit()
{
   printf("\nUsage: art2fsl [-v] -T <transformation matrix> -o <FSL matrix> "
   "-Snx <int> -Sny <int> -Snz <int> -Sdx <float> -Sdy <float> -Sdz <float>\n"
   "\nRequired:\n"
   "\t-i <ART matrix file>: Input FSL matrix file usually 'something'.mrx\n"
   "\t-o <FSL matrix file>: Output ART matrix file usually 'something'.mat\n"
   "\t-Snx <int> -Sny <int> -Snz <int>: 'Subject' (aka 'moving' or FSL's 'input') image matrix dimensions\n"
   "\t-Sdx <float> -Sdy <float> -Sdz <float>: 'Subject' (aka 'moving' or FSL's 'input') image voxel dimensions (mm)\n"
   "\nOptions:\n"
   "\t-v Enables verbose mode\n" 
   );

   exit(0);
}

void readinputlm( char *filename , int &n, float * &inputLM)
{
   FILE *fp;
   float dum;

   n = 0;

   fp=fopen(filename,"r");
   while( fscanf(fp,"%f", &dum) != EOF ) n++;
   fclose(fp);

   inputLM = (float *)calloc(n, sizeof(float) );
   fp=fopen(filename,"r");
   for(int i=0; i<n; i++) fscanf(fp,"%f", inputLM+i);
   fclose(fp);

   n /= 3;
}

int main(int argc, char **argv)
{
   char inputmatrixfile[DEFAULT_STRING_LENGTH]="";
   char inputlmfile[DEFAULT_STRING_LENGTH]="";
   char outputmatrixfile[DEFAULT_STRING_LENGTH]="";
   float T[16];
   DIM sub_dim;
   FILE *fp;

   //initialization to avoid complaining from the compiler
   sub_dim.dx=sub_dim.dy=sub_dim.dz=0.0;
   sub_dim.nx=sub_dim.ny=sub_dim.nz=0;

   if(argc==1) print_help_and_exit();

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) 
      {
         case '1':
            sub_dim.dx = atof(optarg);
            break;
         case '2':
            sub_dim.dy = atof(optarg);
            break;
         case '3':
            sub_dim.dz = atof(optarg);
            break;
         case '4':
            sub_dim.nx = atoi(optarg);
            break;
         case '5':
            sub_dim.ny = atoi(optarg);
            break;
         case '6':
            sub_dim.nz = atoi(optarg);
            break;
         case 'T':
            sprintf(inputmatrixfile,"%s",optarg);
            break;
         case 'i':
            sprintf(inputlmfile,"%s",optarg);
            break;
         case 'o':
            sprintf(outputmatrixfile,"%s",optarg);
            break;
         case 'v':
            opt_v=YES;
            break;
         case 'h':
            print_help_and_exit();
            break;
         case '?':
            print_help_and_exit();
      }
   }

   if(opt_v)
   {
      printf("Input (ART) transformation matrix file: %s\n",inputmatrixfile);
      printf("Input landmarks file: %s\n",inputlmfile);
      printf("Image matrix dimensions: %d x %d x %d\n",sub_dim.nx,sub_dim.ny,sub_dim.nz);
      printf("Image voxel dimensions: %6.4f x %6.4f x %6.4f\n",sub_dim.dx,sub_dim.dy,sub_dim.dz);
   }

   loadTransformation(inputmatrixfile, T);

   if(opt_v)
   {
      printMatrix(T, 4, 4, "Input ART transformation:", NULL);
   }

   int n; // number of landmarks
   float *inputLM;
   float *outputLM;
   readinputlm( inputlmfile, n, inputLM);
   outputLM = (float *)calloc(3*n, sizeof(float) );

   printf("Number of landmarks =%d\n",n);
   for(int i=0; i<n; i++)
   {
      printf("%5.1f %5.1f %5.1f\n", inputLM[i*3 + 0], inputLM[i*3 + 1], inputLM[i*3 + 2]);
   }

   printf("Input landmarks in xyz coordinates:\n");
   for(int i=0; i<n; i++)
   {
      inputLM[i*3 + 0] = inputLM[i*3 + 0]*sub_dim.dx - sub_dim.dx*(sub_dim.nx-1.0)/2.0;
      inputLM[i*3 + 1] = inputLM[i*3 + 1]*sub_dim.dy - sub_dim.dy*(sub_dim.ny-1.0)/2.0;
      inputLM[i*3 + 2] = inputLM[i*3 + 2]*sub_dim.dz - sub_dim.dz*(sub_dim.nz-1.0)/2.0;
      printf("%7.3f %7.3f %7.3f\n", inputLM[i*3 + 0], inputLM[i*3 + 1], inputLM[i*3 + 2]);
   }

   printf("Output landmarks in xyz coordinates:\n");
   for(int i=0; i<n; i++)
   {
      outputLM[i*3+0] = T[0]*inputLM[i*3+0] + T[1]*inputLM[i*3+1] + T[2]*inputLM[i*3+2] + T[3];
      outputLM[i*3+1] = T[4]*inputLM[i*3+0] + T[5]*inputLM[i*3+1] + T[6]*inputLM[i*3+2] + T[7];
      outputLM[i*3+2] = T[8]*inputLM[i*3+0] + T[9]*inputLM[i*3+1] + T[10]*inputLM[i*3+2] + T[11];
      printf("%7.3f %7.3f %7.3f\n", outputLM[i*3 + 0], outputLM[i*3 + 1], outputLM[i*3 + 2]);
   }

   free(inputLM);
   free(outputLM);

   return 0;
}

