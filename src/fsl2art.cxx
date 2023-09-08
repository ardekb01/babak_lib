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
   {"-t", 1, 't'},
   {"-trg", 1, 't'},
   {"-s", 1, 's'},
   {"-sub", 1, 's'},
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
   char trgImFile[DEFAULT_STRING_LENGTH]="";
   char subImFile[DEFAULT_STRING_LENGTH]="";
   char inputmatrixfile[DEFAULT_STRING_LENGTH]="";
   char outputmatrixfile[DEFAULT_STRING_LENGTH]="";
   float Mart[16];
   float Mfsl[16];
   DIM sdim, tdim;
   FILE *fp;

   // initialization to avoid complaining from the compiler
   sdim.dx=sdim.dy=sdim.dz=0.0;
   tdim.dx=tdim.dy=tdim.dz=0.0;
   sdim.nx=sdim.ny=sdim.nz=0;
   tdim.nx=tdim.ny=tdim.nz=0;

   if(argc==1) print_help_and_exit();

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) 
      {
         case '1':
            sdim.dx = atof(optarg);
            break;
         case '2':
            sdim.dy = atof(optarg);
            break;
         case '3':
            sdim.dz = atof(optarg);
            break;
         case '4':
            sdim.nx = atoi(optarg);
            break;
         case '5':
            sdim.ny = atoi(optarg);
            break;
         case '6':
            sdim.nz = atoi(optarg);
            break;
         case 'X':
            tdim.dx = atof(optarg);
            break;
         case 'Y':
            tdim.dy = atof(optarg);
            break;
         case 'Z':
            tdim.dz = atof(optarg);
            break;
         case 'x':
            tdim.nx = atoi(optarg);
            break;
         case 'y':
            tdim.ny = atoi(optarg);
            break;
         case 'z':
            tdim.nz = atoi(optarg);
            break;
         case 't':
            sprintf(trgImFile,"%s",optarg);
            break;
         case 's':
            sprintf(subImFile,"%s",optarg);
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
      }
   }

   if(trgImFile[0] != '\0')
   {
      nifti_1_header trghdr;
      trghdr = read_NIFTI_hdr(trgImFile);

      tdim.nx=trghdr.dim[1]; tdim.ny=trghdr.dim[2]; tdim.nz=trghdr.dim[3];
      tdim.dx=trghdr.pixdim[1]; tdim.dy=trghdr.pixdim[2]; tdim.dz=trghdr.pixdim[3];

      // to deal with the sometimes -tive voxel dimensions in SPM/FSL data sets
      if(tdim.dx<0.0) tdim.dx *= -1.0; 
      if(tdim.dy<0.0) tdim.dy *= -1.0; 
      if(tdim.dz<0.0) tdim.dz *= -1.0;
   }

   if(subImFile[0] != '\0')
   {
      nifti_1_header subhdr;
      subhdr = read_NIFTI_hdr(subImFile);

      sdim.nx=subhdr.dim[1]; sdim.ny=subhdr.dim[2]; sdim.nz=subhdr.dim[3];
      sdim.dx=subhdr.pixdim[1]; sdim.dy=subhdr.pixdim[2]; sdim.dz=subhdr.pixdim[3];

      // to deal with the sometimes -tive voxel dimensions in SPM/FSL data sets
      if(sdim.dx<0.0) sdim.dx *= -1.0; 
      if(sdim.dy<0.0) sdim.dy *= -1.0; 
      if(sdim.dz<0.0) sdim.dz *= -1.0;
   }

   if(opt_v)
   {
      printf("Input (FSL) transformation matrix: %s\n",inputmatrixfile);
      printf("Output (ART) transformation matrix: %s\n",outputmatrixfile);
      printf("Target image (reference) matrix: %d x %d x %d\n",tdim.nx,tdim.ny,tdim.nz);
      printf("Target image (reference) voxel: %6.4f x %6.4f x %6.4f\n",tdim.dx,tdim.dy,tdim.dz);
      printf("Subject image (moving) matrix: %d x %d x %d\n",sdim.nx,sdim.ny,sdim.nz);
      printf("Subject image (moving) voxel: %6.4f x %6.4f x %6.4f\n",sdim.dx,sdim.dy,sdim.dz);
   }

   loadTransformation(inputmatrixfile, Mfsl);

   if(opt_v)
   {
      printMatrix(Mfsl, 4, 4, "Input FSL transformation:", NULL);
   }

   fsl_to_art(Mfsl, Mart, sdim, tdim);

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
