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
#include <nifti1_io.h>

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
   {"-FSL", 1, 'i'},
   {"-o", 1, 'o'},
   {"-ART", 1, 'o'},
   {"-v", 0, 'v'},
   {"-h", 0, 'h'},
   {"-help", 0, 'h'},
   {0, 0, 0}
};

void print_help_and_exit()
{
   printf("\nUsage: fsl2art [-v] -FSL <matrix.mat> -ART <matrix.mrx> -t <image.nii> -s <image.nii>"
   "\nRequired:\n"
   "\t-FSL <matrix.mat>: Input FSL matrix file\n"
   "\t-ART <matrix.mrx>: Output ART matrix file\n"
   "\t-t <image.nii>: Target (reference) image for registration\n"
   "\t-s <image.nii>: Subject (moving) image for registration\n"
   "\nOptions:\n"
   "\t-v Enables verbose mode\n" 
   );

   exit(0);
}

void orientation_code(nifti_1_header hdr, char *code)
{
   if(hdr.qform_code > 0 )
   {
      mat44 R;

      R = nifti_quatern_to_mat44(hdr.quatern_b, hdr.quatern_c, hdr.quatern_d, 
      hdr.qoffset_x, hdr.qoffset_y, hdr.qoffset_z, hdr.pixdim[1], hdr.pixdim[2], 
      hdr.pixdim[3], hdr.pixdim[0]);

      code[0] = directionCode(R.m[0][0],R.m[1][0],R.m[2][0]);
      code[1] = directionCode(R.m[0][1],R.m[1][1],R.m[2][1]);
      code[2] = directionCode(R.m[0][2],R.m[1][2],R.m[2][2]);
      code[3] = '\0';
   }
   else if(hdr.sform_code > 0 )
   {
      code[0] = directionCode(hdr.srow_x[0],hdr.srow_y[0],hdr.srow_z[0]);
      code[1] = directionCode(hdr.srow_x[1],hdr.srow_y[1],hdr.srow_z[1]);
      code[2] = directionCode(hdr.srow_x[2],hdr.srow_y[2],hdr.srow_z[2]);
      code[3] = '\0';
   }

   return;
}

int main(int argc, char **argv)
{
   char trgImFile[DEFAULT_STRING_LENGTH]="";
   char subImFile[DEFAULT_STRING_LENGTH]="";
   char FSLmatrixfile[DEFAULT_STRING_LENGTH]="";
   char ARTmatrixfile[DEFAULT_STRING_LENGTH]="";
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
         case 't':
            sprintf(trgImFile,"%s",optarg);
            break;
         case 's':
            sprintf(subImFile,"%s",optarg);
            break;
         case 'i':
            sprintf(FSLmatrixfile,"%s",optarg);
            break;
         case 'o':
            sprintf(ARTmatrixfile,"%s",optarg);
            break;
         case 'v':
            opt_v=YES;
            break;
         case 'h':
            print_help_and_exit();
            break;
      }
   }

   ////////////////////////////////////////////////////////////////////////////////////
   if(trgImFile[0] == '\0')
   {
     printf("Please specify the target images using the \"-t <image.nii>\" flag\n");
     exit(0);
   }

   nifti_1_header trghdr;
   char trg_orientation_code[4];
   trghdr = read_NIFTI_hdr(trgImFile);

   tdim.nx=trghdr.dim[1]; tdim.ny=trghdr.dim[2]; tdim.nz=trghdr.dim[3];
   tdim.dx=trghdr.pixdim[1]; tdim.dy=trghdr.pixdim[2]; tdim.dz=trghdr.pixdim[3];

   // to deal with the sometimes -tive voxel dimensions in SPM/FSL data sets
   if(tdim.dx<0.0) tdim.dx *= -1.0; 
   if(tdim.dy<0.0) tdim.dy *= -1.0; 
   if(tdim.dz<0.0) tdim.dz *= -1.0;

   if(trghdr.qform_code == 0 && trghdr.sform_code == 0)
   {
     printf("I cannot determine the orientation of the target image, because\n"
     "the NIFTI qform_code sform_code are both zero, aborting ...\n");
     exit(0);
   }

   orientation_code(trghdr,trg_orientation_code);

   if(opt_v)
   {
     printf("Target image: %s\n",trgImFile);
     printf("\tMatrix = %d x %d x %d\n",tdim.nx,tdim.ny,tdim.nz);
     printf("\tVoxel = %6.4f x %6.4f x %6.4f\n",tdim.dx,tdim.dy,tdim.dz);
     printf("\tOrientation = %s\n", trg_orientation_code);
   }

   ////////////////////////////////////////////////////////////////////////////////////

   if(subImFile[0] == '\0')
   {
     printf("Please specify the subject images using the \"-s <image.nii>\" flag\n");
     exit(0);
   }

   nifti_1_header subhdr;
   char sub_orientation_code[4];
   subhdr = read_NIFTI_hdr(subImFile);

   sdim.nx=subhdr.dim[1]; sdim.ny=subhdr.dim[2]; sdim.nz=subhdr.dim[3];
   sdim.dx=subhdr.pixdim[1]; sdim.dy=subhdr.pixdim[2]; sdim.dz=subhdr.pixdim[3];

   // to deal with the sometimes -tive voxel dimensions in SPM/FSL data sets
   if(sdim.dx<0.0) sdim.dx *= -1.0; 
   if(sdim.dy<0.0) sdim.dy *= -1.0; 
   if(sdim.dz<0.0) sdim.dz *= -1.0;

   if(subhdr.qform_code == 0 && subhdr.sform_code == 0)
   {
     printf("I cannot determine the orientation of the subject image, because\n"
     "the NIFTI qform_code sform_code are both zero, aborting ...\n");
     exit(0);
   }

   orientation_code(subhdr,sub_orientation_code);

   if(opt_v)
   {
     printf("Subject image: %s\n",subImFile);
     printf("\tMatrix = %d x %d x %d\n",sdim.nx,sdim.ny,sdim.nz);
     printf("\tVoxel = %6.4f x %6.4f x %6.4f\n",sdim.dx,sdim.dy,sdim.dz);
     printf("\tOrientation = %s\n", sub_orientation_code);
   }
   ////////////////////////////////////////////////////////////////////////////////////

   if(FSLmatrixfile[0] == '\0')
   {
     printf("Please the FSL matrix using the \"-FSL <matrix.mat>\" flag\n");
     exit(0);
   }

   loadTransformation(FSLmatrixfile, Mfsl);

   if(opt_v)
   {
      printMatrix(Mfsl, 4, 4, "FSL transformation matrix:", NULL);
   }

   ////////////////////////////////////////////////////////////////////////////////////

   if(ARTmatrixfile[0] == '\0')
   {
     printf("Please the ART matrix using the \"-ART <matrix.mrx>\" flag\n");
     exit(0);
   }

   if(opt_v)
   {
      printf("ART matrix file: %s\n",ARTmatrixfile);
   }
 
   int trgSysFlg;
   int subSysFlg;

   trgSysFlg = hand_system(trg_orientation_code);
   subSysFlg = hand_system(sub_orientation_code);

   fsl_to_art(Mfsl, Mart, sdim, tdim, subSysFlg, trgSysFlg);

   if(opt_v)
   {
      printMatrix(Mart, 4, 4, "ART transformation matrix:", NULL);
   }

   fp = fopen(ARTmatrixfile,"w");
   if(fp==NULL) file_open_error(ARTmatrixfile);
   fprintf(fp,"%f %f %f %f\n",Mart[0],Mart[1],Mart[2],Mart[3]);
   fprintf(fp,"%f %f %f %f\n",Mart[4],Mart[5],Mart[6],Mart[7]);
   fprintf(fp,"%f %f %f %f\n",Mart[8],Mart[9],Mart[10],Mart[11]);
   fprintf(fp,"%f %f %f %f\n",Mart[12],Mart[13],Mart[14],Mart[15]);
   fclose(fp);

   return 0;
}
