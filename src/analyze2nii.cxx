
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include "volume.h"
#include <ctype.h>

#include <nifti1_io.h>
#include "babak_lib.h"

#define YES 1
#define NO 0

//////////////////////////////////////////////////////////////////////////////////////////////////

int opt;

static struct option options[] =
{
	{"-v",0,'v'},
	{"--verbose",0,'v'},

	{"-i",1,'i'},
	{"--image",1,'i'},

	{"--orient",1,'O'},
	{"-O",1,'O'},

	{"-h",0,'h'},
	{"--help",0,'h'},

	{0,0,0}
};

int opt_D=NO;
int opt_i=NO;
int opt_orient=NO;
//////////////////////////////////////////////////////////////////////////////////////////////////

void print_help_and_exit()
{
   printf("\nUsage: analyze2nii [-h/--help -v/--verbose ] -O/--orient <code> -i/--image <input volume>.hdr \n\n"
   "Required arguments:\n"
   "-O or --orient <code>: Three-letter orientation code. Examples:\n"
   "\t\tPIL for Posterior-Inferior-Left\n"
   "\t\tRAS for Right-Anterior-Superior\n"
   "-i or --image <input volume>: Input image volume in ANALYZE format\n\n"
   "Optional arguments:\n"
   "-h or --help: Prints help information.\n"
   "-v or --verbose : Enables verbose mode\n"
   );

   exit(0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
void set_srow(char code, float *x, float *y, float *z, float d)
{
   if(code=='P' || code=='p')
   {
      *x=0.0; *y=-d; *z=0.0;
   }
   else if(code=='A' || code=='a' )
   {
      *x=0.0; *y=d; *z=0.0;
   }
   else if(code=='L' ||  code=='l')
   {
      *x=-d; *y=0.0; *z=0.0;
   }
   else if(code=='R' ||  code=='r')
   {
      *x=d; *y=0.0; *z=0.0;
   }
   else if(code=='I' ||  code=='i')
   {
      *x=0.0; *y=0.0; *z=-d;
   }
   else if(code=='S' ||  code=='s')
   {
      *x=0.0; *y=0.0; *z=d;
   }
}

void set_shift(char code, float *x, float *y, float *z, float d)
{
   if(code=='P' || code=='p')
   {
      *y=-d;
   }
   else if(code=='A' || code=='a' )
   {
      *y=d;
   }
   else if(code=='L' ||  code=='l')
   {
      *x=-d;
   }
   else if(code=='R' ||  code=='r')
   {
      *x=d;
   }
   else if(code=='I' ||  code=='i')
   {
      *z=-d;
   }
   else if(code=='S' ||  code=='s')
   {
      *z=d;
   }
}

//////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   FILE *fp;
   nifti_1_header hdr;
   nifti1_extender ext;
   char *imgname;
   char *im=NULL;
   unsigned int nbytes;

   char analyzefile[512];
   char niftifile[512];
   char orientation[4];
   int swapflg=0;

   while ((opt = getoption(argc, argv, options)) != -1 )
   {
      switch (opt) 
      {
         case 'i':
            snprintf(analyzefile,sizeof(analyzefile),"%s",optarg);
            opt_i=YES;
            break;
         case 'D':
            opt_D=YES;
            break;
         case 'v':
            opt_v=YES;
            break;
         case 'h':
            print_help_and_exit();
            break;
         case 'O':
            snprintf(orientation,sizeof(orientation),"%s",optarg);
            opt_orient=YES;
            break;
         case '?':
            print_help_and_exit();
      }
   }

   if(!opt_i)
   {
      errorMessage("Error: No input image.  Please specify an image using -i or --image argument.");
   }

   if(opt_v)
   {
      printf("\nInput image: %s\n",analyzefile);
   }

   if(!opt_orient)
   {
      errorMessage("Error: No orientation code.  Please specify a code using -O or --orient argument.");
   }

   if(opt_v)
   {
      printf("\nInput image orientation: %s\n",orientation);
   }

   if ( isOrientationCodeValid(orientation) == 0)
   {
      errorMessage("Error: Invalid orientation code. The specified code is not one of the 48 legal ones.");
   }

   // ensure that the specified image has extension: hdr
   if( !extension_is_hdr(analyzefile) )
   {
      errorMessage("Error: The image filename must have an `hdr' extension.");
   }

   fp = fopen(analyzefile,"r");

   if(fp==NULL)
   {
      errorMessage("Error: I have trouble opening the specified image file.");
   }

   if( fread(&hdr, sizeof(nifti_1_header), 1, fp) != 1 )
   {
      errorMessage("Error: I have trouble reading the specified image file.");
   }

   fclose(fp);

   // if dim[0] is outside range 1..7, then the header information
   // needs to be byte swapped appropriately
   if( (hdr.dim[0]<1 || hdr.dim[0]>7) && hdr.sizeof_hdr != 348) 
   {
      swapflg=1;
      swapniftiheader(&hdr);
   }

   hdr.dim_info=0;
   hdr.intent_p1=0.0;
   hdr.intent_p2=0.0;
   hdr.intent_p3=0.0;
   hdr.intent_code=NIFTI_INTENT_NONE; //no intention is indicated in the header
   hdr.slice_start=0;
   hdr.pixdim[0]=1.0;
   hdr.vox_offset=352;
   hdr.scl_slope=1.0;
   hdr.scl_inter=0.0;
   hdr.slice_end=0;
   hdr.slice_code=0;
   hdr.xyzt_units = NIFTI_UNITS_UNKNOWN;
   hdr.cal_max=0.0;
   hdr.cal_min=0.0;
   hdr.slice_duration=0.0;
   hdr.toffset=0.0;
   hdr.descrip[0]=0;
   hdr.aux_file[0]=0;
   hdr.qform_code=1;
   hdr.sform_code=1;

   set_srow(orientation[0], &(hdr.srow_x[0]), &(hdr.srow_y[0]), &(hdr.srow_z[0]), hdr.pixdim[1]);
   set_srow(orientation[1], &(hdr.srow_x[1]), &(hdr.srow_y[1]), &(hdr.srow_z[1]), hdr.pixdim[2]);
   set_srow(orientation[2], &(hdr.srow_x[2]), &(hdr.srow_y[2]), &(hdr.srow_z[2]), hdr.pixdim[3]);

   set_shift(orientation[0], &(hdr.srow_x[3]), &(hdr.srow_y[3]), &(hdr.srow_z[3]), -hdr.pixdim[1]*(hdr.dim[1]-1.0)/2.0);
   set_shift(orientation[1], &(hdr.srow_x[3]), &(hdr.srow_y[3]), &(hdr.srow_z[3]), -hdr.pixdim[2]*(hdr.dim[2]-1.0)/2.0);
   set_shift(orientation[2], &(hdr.srow_x[3]), &(hdr.srow_y[3]), &(hdr.srow_z[3]), -hdr.pixdim[3]*(hdr.dim[3]-1.0)/2.0);

   hdr.intent_name[0]=0;
   snprintf(hdr.magic,sizeof(hdr.magic),"n+1"); 

   // sets the Quaternions
   {
      float dx, dy, dz;
      mat44 R;

      R.m[0][0]=hdr.srow_x[0]; R.m[0][1]=hdr.srow_x[1]; R.m[0][2]=hdr.srow_x[2]; R.m[0][3]=hdr.srow_x[3];
      R.m[1][0]=hdr.srow_y[0]; R.m[1][1]=hdr.srow_y[1]; R.m[1][2]=hdr.srow_y[2]; R.m[1][3]=hdr.srow_y[3];
      R.m[2][0]=hdr.srow_z[0]; R.m[2][1]=hdr.srow_z[1]; R.m[2][2]=hdr.srow_z[2]; R.m[2][3]=hdr.srow_z[3];
   
      nifti_mat44_to_quatern( R , &(hdr.quatern_b), &(hdr.quatern_c), &(hdr.quatern_d),
      &(hdr.qoffset_x), &(hdr.qoffset_y), &(hdr.qoffset_z), &dx, &dy, &dz, &(hdr.pixdim[0])) ;
   }

   // just for testing
/*
   {
      mat44 S;

      S = nifti_quatern_to_mat44( hdr.quatern_b, hdr.quatern_c, hdr.quatern_d,
      hdr.qoffset_x, hdr.qoffset_y, hdr.qoffset_z, hdr.pixdim[1], hdr.pixdim[2], hdr.pixdim[3], hdr.pixdim[0]);

      printf("%f %f %f %f\n",hdr.srow_x[0], hdr.srow_x[1],hdr.srow_x[2], hdr.srow_x[3]);
      printf("%f %f %f %f\n",hdr.srow_y[0], hdr.srow_y[1],hdr.srow_y[2], hdr.srow_y[3]);
      printf("%f %f %f %f\n\n",hdr.srow_z[0], hdr.srow_z[1],hdr.srow_z[2], hdr.srow_z[3]);

      printf("%f %f %f %f\n", S.m[0][0], S.m[0][1], S.m[0][2], S.m[0][3]);
      printf("%f %f %f %f\n", S.m[1][0], S.m[1][1], S.m[1][2], S.m[1][3]);
      printf("%f %f %f %f\n\n", S.m[2][0], S.m[2][1], S.m[2][2], S.m[2][3]);
   }
*/

   ext.extension[0]=0;
   ext.extension[1]=0;
   ext.extension[2]=0;
   ext.extension[3]=0;

   {
      int L;

      L = strlen(analyzefile);

      imgname = (char *)calloc(L+1, 1);

      strcpy(imgname, analyzefile);

      imgname[L-3]='i';
      imgname[L-2]='m';
      imgname[L-1]='g';

      //printf("%s\n",imgname);
   }

   nbytes = hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.bitpix/8;

   im = (char *)calloc(nbytes, 1);
   
   if(im==NULL)
   {
      errorMessage("Error: Memory allocation problem.");
   }

   fp = fopen(imgname,"r");

   if(fp==NULL)
   {
      errorMessage("Error: I have trouble opening the specified image file.");
   }

   if( fread(im, 1, nbytes, fp) != nbytes )
   {
      errorMessage("Error: I have trouble reading the specified image file.");
   }

   fclose(fp);

   if(swapflg)
   {
      if( hdr.datatype == DT_SIGNED_SHORT || hdr.datatype == DT_UINT16) 
      {
         nifti_swap_2bytes( nbytes/2 , im );
      }
      else if( hdr.datatype == DT_SIGNED_INT || hdr.datatype == DT_FLOAT ) 
      {
         nifti_swap_4bytes( nbytes/4 , im );
      }
      else if( hdr.datatype == DT_DOUBLE ) 
      {
         nifti_swap_8bytes( nbytes/8 , im );
      }
   }

   getfilename(niftifile, analyzefile);
   {
      int L;

      L = strlen(niftifile);

      if( (niftifile[L-3]=='i' && niftifile[L-2]=='m' && niftifile[L-1]=='g') 
      ||  (niftifile[L-3]=='h' && niftifile[L-2]=='d' && niftifile[L-1]=='r') )
      {
         niftifile[L-3]='n'; 
         niftifile[L-2]='i';  
         niftifile[L-1]='i';
      }
      else
      {
         snprintf(niftifile,sizeof(niftifile),"%s.nii",niftifile);
      }
   }

   fp = fopen(niftifile,"w");

   if(fp==NULL)
   {
      errorMessage("Error: I have trouble opening the output NIFIT file.");
   }

   if( fwrite(&hdr, 348, 1, fp) != 1 )
   {
      errorMessage("Error: I have trouble writing to the output NIFIT file.");
   }

   if( fwrite(&ext, 4, 1, fp) != 1 )
   {
      errorMessage("Error: I have trouble writing to the output NIFIT file.");
   }

   if( fwrite(im, 1, nbytes, fp) != nbytes )
   {
      errorMessage("Error: I have trouble writing to the output NIFIT file.");
   }

   fclose(fp);
}
