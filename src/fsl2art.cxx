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
   {"-i", 1, 'i'},
   {"-o", 1, 'o'},
   {"-v", 0, 'v'},
   {"-h", 0, 'h'},
   {"-Tdx",1,'X'},
   {"-Tdy",1,'Y'},
   {"-Tdz",1,'Z'},
   {"-Tnx",1,'x'},
   {"-Tny",1,'y'},
   {"-Tnz",1,'z'},
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
   printf("\nUsage: fsl2art [-v] -i <FSL matrix> -o <ART matrix> "
   "-Tnx <int> -Tny <int> -Tnz <int> -Tdx <float> -Tdy <float> -Tdz <float> "
   "-Snx <int> -Sny <int> -Snz <int> -Sdx <float> -Sdy <float> -Sdz <float>\n"
   "\nRequired:\n"
   "\t-i <FSL matrix file>: Input FSL matrix file usually 'something'.mat\n"
   "\t-o <ART matrix file>: Output ART matrix file usually 'something'.mrx\n"
   "\t-Tnx <int> -Tny <int> -Tnz <int>: 'Target' (aka 'reference') image matrix dimensions\n"
   "\t-Tdx <float> -Tdy <float> -Tdz <float>: 'Target' (aka 'reference') image voxel dimensions (mm)\n"
   "\t-Snx <int> -Sny <int> -Snz <int>: 'Subject' (aka 'moving' or FSL's 'input') image matrix dimensions\n"
   "\t-Sdx <float> -Sdy <float> -Sdz <float>: 'Subject' (aka 'moving' or FSL's 'input') image voxel dimensions (mm)\n"
   "\nOptions:\n"
   "\t-v Enables verbose mode\n" 
   );

   exit(0);
}

int main(int argc, char **argv)
{
   char inputmatrixfile[DEFAULT_STRING_LENGTH]="";
   char outputmatrixfile[DEFAULT_STRING_LENGTH]="";
   float Mart[16];
   float Mfsl[16];
   DIM sub_dim, trg_dim;
   FILE *fp;

   // initialization to avoid complaining from the compiler
   sub_dim.dx=sub_dim.dy=sub_dim.dz=0.0;
   trg_dim.dx=trg_dim.dy=trg_dim.dz=0.0;
   sub_dim.nx=sub_dim.ny=sub_dim.nz=0;
   trg_dim.nx=trg_dim.ny=trg_dim.nz=0;

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
         case 'X':
            trg_dim.dx = atof(optarg);
            break;
         case 'Y':
            trg_dim.dy = atof(optarg);
            break;
         case 'Z':
            trg_dim.dz = atof(optarg);
            break;
         case 'x':
            trg_dim.nx = atoi(optarg);
            break;
         case 'y':
            trg_dim.ny = atoi(optarg);
            break;
         case 'z':
            trg_dim.nz = atoi(optarg);
            break;
         case 'i':
            sprintf(inputmatrixfile,"%s",optarg);
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
      printf("Input (FSL) transformation matrix file: %s\n",inputmatrixfile);
      printf("Output (ART) transformation matrix file: %s\n",outputmatrixfile);
      printf("Target (reference) image matrix dimensions: %d x %d x %d\n",trg_dim.nx,trg_dim.ny,trg_dim.nz);
      printf("Target (reference) image voxel dimensions: %6.4f x %6.4f x %6.4f\n",trg_dim.dx,trg_dim.dy,trg_dim.dz);
      printf("Subject (moving) image matrix dimensions: %d x %d x %d\n",sub_dim.nx,sub_dim.ny,sub_dim.nz);
      printf("Subject (moving) image voxel dimensions: %6.4f x %6.4f x %6.4f\n",sub_dim.dx,sub_dim.dy,sub_dim.dz);
   }

   loadTransformation(inputmatrixfile, Mfsl);

   if(opt_v)
   {
      printMatrix(Mfsl, 4, 4, "Input FSL transformation:", NULL);
   }

   fsl_to_art(Mfsl, Mart, sub_dim, trg_dim);

   if(opt_v)
   {
      printMatrix(Mart, 4, 4, "Output ART transformation:", NULL);
   }

   fp = fopen(outputmatrixfile,"w");
   if(fp==NULL) file_open_error(outputmatrixfile);
   fprintf(fp,"%f %f %f %f\n",Mart[0],Mart[1],Mart[2],Mart[3]);
   fprintf(fp,"%f %f %f %f\n",Mart[4],Mart[5],Mart[6],Mart[7]);
   fprintf(fp,"%f %f %f %f\n",Mart[8],Mart[9],Mart[10],Mart[11]);
   fprintf(fp,"%f %f %f %f\n",Mart[12],Mart[13],Mart[14],Mart[15]);
   fclose(fp);

   return 0;
}
