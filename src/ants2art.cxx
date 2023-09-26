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
   {"-ANTS", 1, 'i'},
   {"-o", 1, 'o'},
   {"-ART", 1, 'o'},
   {"-v", 0, 'v'},
   {"-h", 0, 'h'},
   {"-help", 0, 'h'},
   {0, 0, 0}
};

void print_help_and_exit()
{
   printf("\nUsage: ants2art [-v] -ANTS <matrix.txt> -ART <matrix.mrx> -t <image.nii> -s <image.nii>"
   "\nRequired:\n"
   "\t-ANTS <matrix.txt>: Input ANTS matrix file (after running ConvertTransformFile)\n"
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
   char ANTSmatrixfile[DEFAULT_STRING_LENGTH]="";
   char ARTmatrixfile[DEFAULT_STRING_LENGTH]="";
   float *Mart;
   float Mants[16];
   DIM sdim, tdim;

   // initialization to avoid complaining from the compiler
   sdim.dx=sdim.dy=sdim.dz=0.0;
   tdim.dx=tdim.dy=tdim.dz=0.0;
   sdim.nx=sdim.ny=sdim.nz=0;
   tdim.nx=tdim.ny=tdim.nz=0;

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
            sprintf(ANTSmatrixfile,"%s",optarg);
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
     printf("Please specify the target image using the \"-t <image.nii>\" flag\n");
     exit(0);
   }

   nifti_1_header trghdr;
   char trg_orientation_code[4];
   float txc, tyc, tzc;

   trghdr = read_NIFTI_hdr(trgImFile);

   tdim.nx=trghdr.dim[1]; tdim.ny=trghdr.dim[2]; tdim.nz=trghdr.dim[3];
   tdim.dx=trghdr.pixdim[1]; tdim.dy=trghdr.pixdim[2]; tdim.dz=trghdr.pixdim[3];

   // to deal with the sometimes -tive voxel dimensions in SPM/FSL data sets
   if(tdim.dx<0.0) tdim.dx *= -1.0; 
   if(tdim.dy<0.0) tdim.dy *= -1.0; 
   if(tdim.dz<0.0) tdim.dz *= -1.0;

   txc = (tdim.nx-1.0)*tdim.dx/2.0;
   tyc = (tdim.ny-1.0)*tdim.dy/2.0;
   tzc = (tdim.nz-1.0)*tdim.dz/2.0;

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
     printf("Please specify the subject image using the \"-s <image.nii>\" flag\n");
     exit(0);
   }

   nifti_1_header subhdr;
   char sub_orientation_code[4];
   float sxc, syc, szc;

   subhdr = read_NIFTI_hdr(subImFile);

   sdim.nx=subhdr.dim[1]; sdim.ny=subhdr.dim[2]; sdim.nz=subhdr.dim[3];
   sdim.dx=subhdr.pixdim[1]; sdim.dy=subhdr.pixdim[2]; sdim.dz=subhdr.pixdim[3];

   // to deal with the sometimes -tive voxel dimensions in SPM/FSL data sets
   if(sdim.dx<0.0) sdim.dx *= -1.0; 
   if(sdim.dy<0.0) sdim.dy *= -1.0; 
   if(sdim.dz<0.0) sdim.dz *= -1.0;

   sxc = (sdim.nx-1.0)*sdim.dx/2.0;
   syc = (sdim.ny-1.0)*sdim.dy/2.0;
   szc = (sdim.nz-1.0)*sdim.dz/2.0;

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

   if(ANTSmatrixfile[0] == '\0')
   {
     printf("Please the ANTS matrix using the \"-ANTS <matrix.txt>\" flag\n");
     exit(0);
   }

   if(opt_v)
   {
     printf("ANTS matrix file: %s\n",ANTSmatrixfile);
   }

   FILE *fp;
   char line[1000];
   char dum[1024];
   float R[9];
   float Rc[3];
   float t[3];
   float c[3];

   fp=fopen(ANTSmatrixfile,"r");
   fgets(line,1000,fp);
   fgets(line,1000,fp);
   fgets(line,1000,fp);
   fgets(line,1000,fp);
   sscanf(line,"%s %f %f %f %f %f %f %f %f %f %f %f %f",dum,&R[0],&R[1],&R[2],&R[3],&R[4],&R[5],
		   &R[6],&R[7],&R[8],&t[0],&t[1],&t[2]);
   fgets(line,1000,fp);
   sscanf(line,"%s %f %f %f",dum,&c[0],&c[1],&c[2]);

   if(opt_v)
   {
      printf("\tRotation:\n");
      printf("\t\t%f %f %f\n",R[0],R[1],R[2]);
      printf("\t\t%f %f %f\n",R[3],R[4],R[5]);
      printf("\t\t%f %f %f\n",R[6],R[7],R[8]);

      printf("\tTranslation:\n");
      printf("\t\t%f %f %f\n",t[0],t[1],t[2]);

      printf("\tCenter of rotation:\n");
      printf("\t\t%f %f %f\n",c[0],c[1],c[2]);
   }
   fclose(fp);
  
   // This gives Mants in LPS->LPS coordinates
   multi(R,3,3,c,3,1,Rc);
   Mants[3] =  c[0]+t[0]-Rc[0];
   Mants[7] =  c[1]+t[1]-Rc[1];
   Mants[11] = c[2]+t[2]-Rc[2];

   Mants[0]=R[0]; Mants[1]=R[1]; Mants[2]=R[2]; 
   Mants[4]=R[3]; Mants[5]=R[4]; Mants[6]=R[5];
   Mants[8]=R[6]; Mants[9]=R[7]; Mants[10]=R[8]; 
  
   Mants[12]=0.0; Mants[13]=0.0; Mants[14]=0.0; Mants[15]=1.0;

   // These negations convert Mants to RAS->RAS coordinates
   Mants[2] *= (-1);
   Mants[3] *= (-1);
   Mants[6] *= (-1);
   Mants[7] *= (-1);
   Mants[8] *= (-1);
   Mants[9] *= (-1);
   
   if(opt_v)
   {
      printf("\tRAS-to-RAS Matrix\n");
      printf("\t\t%9.6f  %9.6f  %9.6f  %9.6f\n",Mants[0],Mants[1],Mants[2],Mants[3]);
      printf("\t\t%9.6f  %9.6f  %9.6f  %9.6f\n",Mants[4],Mants[5],Mants[6],Mants[7]);
      printf("\t\t%9.6f  %9.6f  %9.6f  %9.6f\n",Mants[8],Mants[9],Mants[10],Mants[11]);
      printf("\t\t%9.6f  %9.6f  %9.6f  %9.6f\n",Mants[12],Mants[13],Mants[14],Mants[15]);
   }

   // The following converts it to ART matrix format
   Mants[3] -= (sxc + subhdr.qoffset_x);
   Mants[7] -= (syc + subhdr.qoffset_y);
   Mants[11] -= (szc + subhdr.qoffset_z);
   Mart=inv4(Mants);
   Mart[3] -= (txc + trghdr.qoffset_x);
   Mart[7] -= (tyc + trghdr.qoffset_y);
   Mart[11] -= (tzc + trghdr.qoffset_z);

   if(opt_v)
   {
      printf("ART transformation matrix:\n");
      printf("%9.6f  %9.6f  %9.6f  %9.6f\n",Mart[0],Mart[1],Mart[2],Mart[3]);
      printf("%9.6f  %9.6f  %9.6f  %9.6f\n",Mart[4],Mart[5],Mart[6],Mart[7]);
      printf("%9.6f  %9.6f  %9.6f  %9.6f\n",Mart[8],Mart[9],Mart[10],Mart[11]);
      printf("%9.6f  %9.6f  %9.6f  %9.6f\n",Mart[12],Mart[13],Mart[14],Mart[15]);
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

   fp = fopen(ARTmatrixfile,"w");
   if(fp==NULL) file_open_error(ARTmatrixfile);
   fprintf(fp,"%f %f %f %f\n",Mart[0],Mart[1],Mart[2],Mart[3]);
   fprintf(fp,"%f %f %f %f\n",Mart[4],Mart[5],Mart[6],Mart[7]);
   fprintf(fp,"%f %f %f %f\n",Mart[8],Mart[9],Mart[10],Mart[11]);
   fprintf(fp,"%f %f %f %f\n",Mart[12],Mart[13],Mart[14],Mart[15]);
   fclose(fp);

   return 0;
}
