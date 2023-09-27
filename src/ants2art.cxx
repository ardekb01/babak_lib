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
   float Qart[16];  // This is the matrix that converts target ijk to xyz in ART;
   float Qant[16];  // This is the matrix that converts target ijk to xyz in ANTS;
   float *iQart; // inverse of Qart
   float *iQant; // inverse of Qant
   float Q[16];  // Qart * iQant

   if(trgImFile[0] == '\0')
   {
     printf("Please specify the target image using the \"-t <image.nii>\" flag\n");
     exit(0);
   }

   nifti_1_header thdr;
   char trg_orientation_code[4];
   float txc, tyc, tzc;

   thdr = read_NIFTI_hdr(trgImFile);

   tdim.nx=thdr.dim[1]; tdim.ny=thdr.dim[2]; tdim.nz=thdr.dim[3];
   tdim.dx=thdr.pixdim[1]; tdim.dy=thdr.pixdim[2]; tdim.dz=thdr.pixdim[3];

   // to deal with the sometimes -tive voxel dimensions in SPM/FSL data sets
   if(tdim.dx<0.0) tdim.dx *= -1.0; 
   if(tdim.dy<0.0) tdim.dy *= -1.0; 
   if(tdim.dz<0.0) tdim.dz *= -1.0;

   txc = (tdim.nx-1.0)*tdim.dx/2.0;
   tyc = (tdim.ny-1.0)*tdim.dy/2.0;
   tzc = (tdim.nz-1.0)*tdim.dz/2.0;

   if(thdr.qform_code == 0 && thdr.sform_code == 0)
   {
     printf("I cannot determine the orientation of the target image, because\n"
     "the NIFTI qform_code sform_code are both zero, aborting ...\n");
     exit(0);
   }

   orientation_code(thdr,trg_orientation_code);

   if(opt_v)
   {
     printf("Target image: %s\n",trgImFile);
     printf("\tMatrix = %d x %d x %d\n",tdim.nx,tdim.ny,tdim.nz);
     printf("\tVoxel = %6.4f x %6.4f x %6.4f\n",tdim.dx,tdim.dy,tdim.dz);
     printf("\tOrientation = %s\n", trg_orientation_code);
   }

   Qart[0]=tdim.dx; Qart[1]=0.0;     Qart[2]=0.0;      Qart[3]=-txc;
   Qart[4]=0.0;     Qart[5]=tdim.dy; Qart[6]=0.0;      Qart[7]=-tyc;
   Qart[8]=0.0;     Qart[9]=0.0;     Qart[10]=tdim.dz; Qart[11]=-tzc;
   Qart[12]=0.0;    Qart[13]=0.0;    Qart[14]=0.0;     Qart[15]=1.0;

   Qant[0]=thdr.srow_x[0]; Qant[1]=thdr.srow_x[1]; Qant[2]=thdr.srow_x[2]; Qant[3]=thdr.srow_x[3];
   Qant[4]=thdr.srow_y[0]; Qant[5]=thdr.srow_y[1]; Qant[6]=thdr.srow_y[2]; Qant[7]=thdr.srow_y[3];
   Qant[8]=thdr.srow_z[0]; Qant[9]=thdr.srow_z[1]; Qant[10]=thdr.srow_z[2];Qant[11]=thdr.srow_z[3];
   Qant[12]=0.0;    Qant[13]=0.0;    Qant[14]=0.0;     Qant[15]=1.0;

   iQart=inv4(Qart);
   iQant=inv4(Qant);
   
   multi(Qart,4,4,iQant,4,4,Q);

   ////////////////////////////////////////////////////////////////////////////////////

   float Sart[16];  // This is the matrix that converts subject ijk to xyz in ART;
   float Sant[16];  // This is the matrix that converts subject ijk to xyz in ANTS;
   float *iSart; // inverse of Sart
   float *iSant; // inverse of Sant
   float S[16];  // Sant * iSart

   if(subImFile[0] == '\0')
   {
     printf("Please specify the subject image using the \"-s <image.nii>\" flag\n");
     exit(0);
   }

   nifti_1_header shdr;
   char sub_orientation_code[4];
   float sxc, syc, szc;

   shdr = read_NIFTI_hdr(subImFile);

   sdim.nx=shdr.dim[1]; sdim.ny=shdr.dim[2]; sdim.nz=shdr.dim[3];
   sdim.dx=shdr.pixdim[1]; sdim.dy=shdr.pixdim[2]; sdim.dz=shdr.pixdim[3];

   // to deal with the sometimes -tive voxel dimensions in SPM/FSL data sets
   if(sdim.dx<0.0) sdim.dx *= -1.0; 
   if(sdim.dy<0.0) sdim.dy *= -1.0; 
   if(sdim.dz<0.0) sdim.dz *= -1.0;

   sxc = (sdim.nx-1.0)*sdim.dx/2.0;
   syc = (sdim.ny-1.0)*sdim.dy/2.0;
   szc = (sdim.nz-1.0)*sdim.dz/2.0;

   if(shdr.qform_code == 0 && shdr.sform_code == 0)
   {
     printf("I cannot determine the orientation of the subject image, because\n"
     "the NIFTI qform_code sform_code are both zero, aborting ...\n");
     exit(0);
   }

   orientation_code(shdr,sub_orientation_code);

   if(opt_v)
   {
     printf("Subject image: %s\n",subImFile);
     printf("\tMatrix = %d x %d x %d\n",sdim.nx,sdim.ny,sdim.nz);
     printf("\tVoxel = %6.4f x %6.4f x %6.4f\n",sdim.dx,sdim.dy,sdim.dz);
     printf("\tOrientation = %s\n", sub_orientation_code);
   }

   Sart[0]=sdim.dx; Sart[1]=0.0;     Sart[2]=0.0;      Sart[3]=-sxc;
   Sart[4]=0.0;     Sart[5]=sdim.dy; Sart[6]=0.0;      Sart[7]=-syc;
   Sart[8]=0.0;     Sart[9]=0.0;     Sart[10]=sdim.dz; Sart[11]=-szc;
   Sart[12]=0.0;    Sart[13]=0.0;    Sart[14]=0.0;     Sart[15]=1.0;

   Sant[0]=shdr.srow_x[0]; Sant[1]=shdr.srow_x[1]; Sant[2]=shdr.srow_x[2]; Sant[3]=shdr.srow_x[3];
   Sant[4]=shdr.srow_y[0]; Sant[5]=shdr.srow_y[1]; Sant[6]=shdr.srow_y[2]; Sant[7]=shdr.srow_y[3];
   Sant[8]=shdr.srow_z[0]; Sant[9]=shdr.srow_z[1]; Sant[10]=shdr.srow_z[2];Sant[11]=shdr.srow_z[3];
   Sant[12]=0.0;    Sant[13]=0.0;    Sant[14]=0.0;     Sant[15]=1.0;

   iSart=inv4(Sart);
   iSant=inv4(Sant);

   multi(Sant,4,4,iSart,4,4,S);

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
   Mart=inv4(Mants);
   multi(Q,4,4,Mart,4,4,Mart);
   multi(Mart,4,4,S,4,4,Mart);

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
