#include <stdio.h>
#include <string.h>
#include <nifti1_io.h>
#include <babak_lib.h>

//////////////////////////////////////////////////////////////////
nifti_1_header read_NIFTI_hdr(const char *filename)
{
   FILE *fp;
   nifti_1_header hdr;

   fp = fopen(filename,"r");

   if(fp==NULL)
   {
      printf("read_NIFTI_hdr(): I could not open file: %s\n", filename);
      return(hdr);
   }

   if( fread(&hdr, sizeof(hdr), 1, fp)!=1 )
   {
      printf("read_NIFTI_hdr(): I could not read the NIFTI header from file: %s\n", filename);
      fclose(fp);
      return(hdr);
   }

   if(hdr.dim[0]<1 || hdr.dim[0]>7)
   {
      swapniftiheader(&hdr);
   }

   return(hdr);
}

// returns 0 on failure
int read_NIFTI_hdr(const char *filename, nifti_1_header *hdr)
{
   FILE *fp;

   fp = fopen(filename,"r");

   if(fp==NULL)
   {
      printf("read_NIFTI_hdr(): I could not open file: %s\n", filename);
      return(0);
   }

   if( fread(hdr, sizeof(nifti_1_header), 1, fp)!=1 )
   {
      printf("read_NIFTI_hdr(): I could not read the NIFTI header from file: %s\n", filename);
      fclose(fp);
      return(0);
   }

   if(hdr->dim[0]<1 || hdr->dim[0]>7)
   {
      swapniftiheader(hdr);
   }

   return(1);
}
//////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
int not_magical_nifti(const char *imagefilename)
{
   FILE *fp;
   nifti_1_header hdr;

   fp = fopen(imagefilename,"r");
   
   if(fp == NULL)
   {
      printf("not_magical_nifti(): could not open file %s\n", imagefilename);
      return(1);
   }

   if( fread(&hdr, sizeof(hdr), 1, fp)!=1 )
   {
      printf("not_magical_nifti(): could not read the NIFTI header from %s\n", imagefilename);
      fclose(fp);
      return(1);
   }

   fclose(fp);

   if( hdr.magic[0]!='n' || (hdr.magic[1]!='+' && hdr.magic[1]!='i') || hdr.magic[2]!='1' )
   {
      printf("hdr.magic = %s\n", hdr.magic);
      printf("not_magical_nifti(): NIFTI file magic field is neither \"n+1\" nor \"ni1\"\n"); 
      return(1);
   }

   return(0);
}

int not_magical_nifti(const char *imagefilename, int verbose)
{
   FILE *fp;
   nifti_1_header hdr;

   fp = fopen(imagefilename,"r");
   
   if(fp == NULL)
   {
      if(verbose) printf("not_magical_nifti(): could not open file %s\n", imagefilename);
      return(1);
   }

   if( fread(&hdr, sizeof(hdr), 1, fp)!=1 )
   {
      if(verbose) printf("not_magical_nifti(): could not read the NIFTI header from %s\n", imagefilename);
      fclose(fp);
      return(1);
   }

   fclose(fp);

   if( hdr.magic[0]!='n' || (hdr.magic[1]!='+' && hdr.magic[1]!='i') || hdr.magic[2]!='1' )
   {
      if(verbose) printf("hdr.magic = %s\n", hdr.magic);
      if(verbose) printf("not_magical_nifti(): NIFTI file magic field is neither \"n+1\" nor \"ni1\"\n"); 
      return(1);
   }

   return(0);
}

////////////////////////////////////////////////////////////////////////////////////

// returns 1 if filename has a .hdr or .nii extnension, 0 otherwise
/*
This function was replaced by checkNiftiFileExtension() function in checkNiftiFileExtension.cxx
int niftiextension(const char *filename)
{
   int L; // number of characters in the filename

   L=strlen(filename);

   if(L<4)
   {
      return(0);
   }

   if( filename[L-4]=='.' && filename[L-3]=='h' && filename[L-2]=='d' && filename[L-1]=='r')
   {
      return(1);
   }

   if( filename[L-4]=='.' && filename[L-3]=='n' && filename[L-2]=='i' && filename[L-1]=='i')
   {
      return(1);
   }

   return(0);
}
*/

void swapniftiheader(nifti_1_header *hdr)
{
   swapByteOrder( (char *)&(hdr->sizeof_hdr), sizeof(int));

   swapByteOrder( (char *)&(hdr->extents), sizeof(int));
   swapByteOrder( (char *)&(hdr->session_error), sizeof(short));

   for(int i=0; i<8; i++)
   {
      swapByteOrder( (char *)&(hdr->dim[i]), sizeof(short));
   }

   swapByteOrder( (char *)&(hdr->intent_p1), sizeof(float));
   swapByteOrder( (char *)&(hdr->intent_p2), sizeof(float));
   swapByteOrder( (char *)&(hdr->intent_p3), sizeof(float));

   swapByteOrder( (char *)&(hdr->intent_code), sizeof(short));
   swapByteOrder( (char *)&(hdr->datatype), sizeof(short));
   swapByteOrder( (char *)&(hdr->bitpix), sizeof(short));
   swapByteOrder( (char *)&(hdr->slice_start), sizeof(short));

   for(int i=0; i<8; i++)
   {
      swapByteOrder( (char *)&(hdr->pixdim[i]), sizeof(float));
   }

   swapByteOrder( (char *)&(hdr->vox_offset), sizeof(float));
   swapByteOrder( (char *)&(hdr->scl_slope), sizeof(float));
   swapByteOrder( (char *)&(hdr->scl_inter), sizeof(float));

   swapByteOrder( (char *)&(hdr->slice_end), sizeof(short));

   swapByteOrder( (char *)&(hdr->cal_max), sizeof(float));
   swapByteOrder( (char *)&(hdr->cal_min), sizeof(float));
   swapByteOrder( (char *)&(hdr->slice_duration), sizeof(float));
   swapByteOrder( (char *)&(hdr->toffset), sizeof(float));

   swapByteOrder( (char *)&(hdr->glmax), sizeof(int));
   swapByteOrder( (char *)&(hdr->glmin), sizeof(int));

   swapByteOrder( (char *)&(hdr->qform_code), sizeof(short));
   swapByteOrder( (char *)&(hdr->sform_code), sizeof(short));

   swapByteOrder( (char *)&(hdr->quatern_b), sizeof(float));
   swapByteOrder( (char *)&(hdr->quatern_c), sizeof(float));
   swapByteOrder( (char *)&(hdr->quatern_d), sizeof(float));
   swapByteOrder( (char *)&(hdr->qoffset_x), sizeof(float));
   swapByteOrder( (char *)&(hdr->qoffset_y), sizeof(float));
   swapByteOrder( (char *)&(hdr->qoffset_z), sizeof(float));

   for(int i=0; i<4; i++)
   {
      swapByteOrder( (char *)&(hdr->srow_x[i]), sizeof(float));
      swapByteOrder( (char *)&(hdr->srow_y[i]), sizeof(float));
      swapByteOrder( (char *)&(hdr->srow_z[i]), sizeof(float));
   }
}

int same_nifti_image_size(int N, char **imagefile, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz)
{
   nifti_1_header hdr;

   hdr=read_NIFTI_hdr(imagefile[0]);

   *nx = hdr.dim[1];
   *ny = hdr.dim[2];
   *nz = hdr.dim[3];
   *dx = hdr.pixdim[1];
   *dy = hdr.pixdim[2];
   *dz = hdr.pixdim[3];

   for(int i=1; i<N; i++)
   {
      hdr=read_NIFTI_hdr(imagefile[i]);

      if( (*nx != hdr.dim[1]) ||  (*ny != hdr.dim[2]) ||  (*nz != hdr.dim[3]) )
      {
         printf("\n\nImage %d: %s",i+1,imagefile[i]);
         printf("\n\tMatrix size = %d x %d x %d", hdr.dim[1], hdr.dim[2], hdr.dim[3]);
         return(0);
      }
   }

   return(1);
}


void read_nifti_image(const char *filename, unsigned char **im, nifti_1_header *hdr)
{
   FILE *fp;
   int nv;

   fp = fopen(filename,"r");
   if(fp==NULL) file_open_error(filename);
   if( fread(hdr, sizeof(nifti_1_header), 1, fp) != 1 )
   {
      printf("Failed to read variable: \"hdr\"\n");
   }
   fclose(fp);

   if(hdr->dim[0]<1 || hdr->dim[0]>7)
   {
      swapniftiheader(hdr);
   }

   nv=1;
   for( int i=1; i<=hdr->dim[0]; i++)
   {
      nv *= hdr->dim[i];
   }

   *im = (unsigned char *)calloc(nv,sizeof(unsigned char));

   fp = fopen(filename,"r");
   if(fp==NULL) file_open_error(filename);
   fseek(fp, (long)hdr->vox_offset, SEEK_SET);
   if( fread(*im, sizeof(unsigned char), nv, fp) != (size_t)nv )
   {
      printf("Failed to read variable: \"im\"\n");
   }
   fclose(fp);
}

void read_nifti_image(const char *filename, short **im, nifti_1_header *hdr)
{
   FILE *fp;
   int nv;

   fp = fopen(filename,"r");
   if(fp==NULL) file_open_error(filename);
   if( fread(hdr, sizeof(nifti_1_header), 1, fp) != 1 )
   {
      printf("Failed to read variable: \"hdr\"\n");
   }
   fclose(fp);

   if(hdr->dim[0]<1 || hdr->dim[0]>7)
   {
      swapniftiheader(hdr);
   }

   nv=1;
   for( int i=1; i<=hdr->dim[0]; i++)
   {
      nv *= hdr->dim[i];
   }

   *im = (short *)calloc(nv,sizeof(short));

   fp = fopen(filename,"r");
   if(fp==NULL) file_open_error(filename);
   fseek(fp, (long)hdr->vox_offset, SEEK_SET);
   if( fread(*im, sizeof(short), nv, fp) != (size_t)nv )
   {
      printf("Failed to read variable: \"im\"\n");
   }
   fclose(fp);
}

void save_nifti_image(const char *filename, char *im, nifti_1_header *hdr)
{
   FILE *fp;
   nifti1_extender extender;
   int nv;

   nv = hdr->dim[1] * hdr->dim[2] * hdr->dim[3];

   hdr->dim[0]=3;
   hdr->sizeof_hdr=348;
   strcpy(hdr->magic,"n+1");
   hdr->vox_offset=352.0;
   extender.extension[0]=0;

   if(hdr->datatype==DT_SIGNED_SHORT)
   {
      hdr->bitpix=16;
   }

   if(hdr->datatype==DT_FLOAT)
   {
      hdr->bitpix=32;
   }

   if(hdr->datatype==DT_UNSIGNED_CHAR)
   {
      hdr->bitpix=8;
   }

   fp=fopen(filename,"w");
   if(fp==NULL) file_open_error(filename);
   fwrite(hdr, sizeof(nifti_1_header),1,fp);
   fwrite(&extender, sizeof(nifti1_extender),1,fp);
   fwrite(im, 1, nv*hdr->bitpix/8, fp);
   fclose(fp);
}

void save_nifti_image(const char *filename, unsigned char *im, nifti_1_header *hdr)
{
   FILE *fp;
   nifti1_extender extender;
   int nv;

   nv = hdr->dim[1] * hdr->dim[2] * hdr->dim[3];

   hdr->dim[0]=3;
   hdr->sizeof_hdr=348;
   strcpy(hdr->magic,"n+1");
   hdr->datatype=DT_UNSIGNED_CHAR;
   hdr->bitpix=8;
   hdr->vox_offset=352.0;
   extender.extension[0]=0;

   fp=fopen(filename,"w");
   if(fp==NULL) file_open_error(filename);
   fwrite(hdr, sizeof(nifti_1_header),1,fp);
   fwrite(&extender, sizeof(nifti1_extender),1,fp);
   fwrite(im, sizeof(unsigned char),nv,fp);
   fclose(fp);
}

void save_nifti_image(const char *filename, short *im, nifti_1_header *hdr)
{
   FILE *fp;
   nifti1_extender extender;
   int nv=1;

   for( int i=1; i<=hdr->dim[0]; i++)
   {
      nv *= hdr->dim[i];
   }
   
   hdr->sizeof_hdr=348;
   strcpy(hdr->magic,"n+1");
   hdr->datatype=DT_SIGNED_SHORT;
   hdr->bitpix=16;
   hdr->vox_offset=352.0;
   extender.extension[0]=0;
   extender.extension[1]=0;
   extender.extension[2]=0;
   extender.extension[3]=0;

   fp=fopen(filename,"w");
   if(fp==NULL) file_open_error(filename);
   fwrite(hdr, sizeof(nifti_1_header),1,fp);
   fwrite(&extender, sizeof(nifti1_extender),1,fp);
   fwrite(im, sizeof(short),nv,fp);
   fclose(fp);
}

void save_nifti_image(const char *filename, float *im, nifti_1_header *hdr)
{
   FILE *fp;
   nifti1_extender extender;
   int nv;

   nv = hdr->dim[1] * hdr->dim[2] * hdr->dim[3];
   
   hdr->dim[0]=3;
   hdr->sizeof_hdr=348;
   strcpy(hdr->magic,"n+1");
   hdr->datatype=DT_FLOAT;
   hdr->bitpix=32;
   hdr->vox_offset=352.0;
   extender.extension[0]=0;

   fp=fopen(filename,"w");
   if(fp==NULL) file_open_error(filename);
   fwrite(hdr, sizeof(nifti_1_header),1,fp);
   fwrite(&extender, sizeof(nifti1_extender),1,fp);
   fwrite(im, sizeof(float),nv,fp);
   fclose(fp);
}


nifti_1_header read_NIFTI_hdr(const char *filename, nifti1_extender *extender, char **extension)
{
   FILE *fp;
   nifti_1_header hdr;
   int extension_size;

   fp = fopen(filename,"r");
   if(fp==NULL) file_open_error(filename);
   if( fread(&hdr, sizeof(nifti_1_header), 1, fp) != 1 )
   {
      printf("Failed to read variable: \"hdr\"\n");
   }
   if( fread(extender, sizeof(nifti1_extender), 1, fp) != 1 )
   {
      printf("Failed to read variable: \"extender\"\n");
   }

   if(hdr.dim[0]<1 || hdr.dim[0]>7)
   {
      swapniftiheader(&hdr);
   }

   extension_size = (int)(hdr.vox_offset-352);
   if(extension_size>0)
   {
      *extension = (char *)calloc(extension_size,1);
      if( fread(*extension, extension_size, 1, fp) != 1 )
      {
         printf("Failed to read variable: \"extension\"\n");
      }
   }
   else
   {
      *extension=NULL;
   }

   fclose(fp);

   return(hdr);
}

void print_NIFTI_hdr(const char *filename)
{
   FILE *fp;
   nifti_1_header hdr;

   fp = fopen(filename,"r");
   if(fp==NULL) file_open_error(filename);
   if( fread(&hdr, sizeof(nifti_1_header), 1, fp) != 1 )
   {
       printf("Failed to read variable: \"hdr\"\n");
   }
   fclose(fp);

   if(hdr.dim[0]<1 || hdr.dim[0]>7)
   {
      swapniftiheader(&hdr);
   }

   printf("sizeof_hdr = %d\n",hdr.sizeof_hdr);

   printf("dim_info = %d\n",hdr.dim_info);
   printf("\tFrequency dimension = %d\n", DIM_INFO_TO_FREQ_DIM(hdr.dim_info) );
   printf("\tPhase dimension = %d\n", DIM_INFO_TO_PHASE_DIM(hdr.dim_info) );
   printf("\tSlice dimension = %d\n", DIM_INFO_TO_SLICE_DIM(hdr.dim_info) );

   printf("dim:\n");
   for(int i=0; i<8; i++)
   {
      printf("\tdim[%d] = %d\n",i,hdr.dim[i]);
   }

   printf("intent_p1 = %f\n",hdr.intent_p1);
   printf("intent_p2 = %f\n",hdr.intent_p2);
   printf("intent_p3 = %f\n",hdr.intent_p3);

   printf("intent_code = %d\n",hdr.intent_code);

   printf("datatype = %d\n",hdr.datatype);

   printf("bitpix = %d\n",hdr.bitpix);
   printf("slice_start = %d\n",hdr.slice_start);

   printf("pixdim:\n");
   for(int i=0; i<8; i++)
   {
      printf("\tpixdim[%d] = %f\n",i,hdr.pixdim[i]);
   }

   printf("vox_offset = %f\n",hdr.vox_offset);
   printf("scl_slope = %f\n",hdr.scl_slope);
   printf("scl_inter = %f\n",hdr.scl_inter);
   printf("slice_end = %d\n",hdr.slice_end);
   printf("slice_code = %d\n",hdr.slice_code);
   printf("xyzt_units = %d\n",hdr.xyzt_units);
   printf("cal_max = %f\n",hdr.cal_max);
   printf("cal_min = %f\n",hdr.cal_min);
   printf("slice_duration = %f\n",hdr.slice_duration);
   printf("toffset = %f\n",hdr.toffset);
   printf("descrip = %s\n",hdr.descrip);
   printf("aux_file = %s\n",hdr.aux_file);
   printf("qform_code = %d\n",hdr.qform_code);
   printf("sform_code = %d\n",hdr.sform_code);
   printf("quatern_b = %f\n",hdr.quatern_b);
   printf("quatern_c = %f\n",hdr.quatern_c);
   printf("quatern_d = %f\n",hdr.quatern_d);
   printf("qoffset_x = %f\n",hdr.qoffset_x);
   printf("qoffset_y = %f\n",hdr.qoffset_y);
   printf("qoffset_z = %f\n",hdr.qoffset_z);
   printf("srow_x = %f %f %f %f\n",hdr.srow_x[0],hdr.srow_x[1],hdr.srow_x[2],hdr.srow_x[3]);
   printf("srow_y = %f %f %f %f\n",hdr.srow_y[0],hdr.srow_y[1],hdr.srow_y[2],hdr.srow_y[3]);
   printf("srow_z = %f %f %f %f\n",hdr.srow_z[0],hdr.srow_z[1],hdr.srow_z[2],hdr.srow_z[3]);
   printf("intent_name = %s\n",hdr.intent_name);
   printf("magic = %s\n",hdr.magic);
}

void print_NIFTI_hdr(nifti_1_header hdr)
{
   if(hdr.dim[0]<1 || hdr.dim[0]>7)
   {
      swapniftiheader(&hdr);
   }

   printf("sizeof_hdr = %d\n",hdr.sizeof_hdr);

   printf("dim_info = %d\n",hdr.dim_info);
   printf("\tFrequency dimension = %d\n", DIM_INFO_TO_FREQ_DIM(hdr.dim_info) );
   printf("\tPhase dimension = %d\n", DIM_INFO_TO_PHASE_DIM(hdr.dim_info) );
   printf("\tSlice dimension = %d\n", DIM_INFO_TO_SLICE_DIM(hdr.dim_info) );

   printf("dim:\n");
   for(int i=0; i<8; i++)
   {
      printf("\tdim[%d] = %d\n",i,hdr.dim[i]);
   }

   printf("intent_p1 = %f\n",hdr.intent_p1);
   printf("intent_p2 = %f\n",hdr.intent_p2);
   printf("intent_p3 = %f\n",hdr.intent_p3);

   printf("intent_code = %d\n",hdr.intent_code);

   printf("datatype = %d\n",hdr.datatype);

   printf("bitpix = %d\n",hdr.bitpix);
   printf("slice_start = %d\n",hdr.slice_start);

   printf("pixdim:\n");
   for(int i=0; i<8; i++)
   {
      printf("\tpixdim[%d] = %f\n",i,hdr.pixdim[i]);
   }

   printf("vox_offset = %f\n",hdr.vox_offset);
   printf("scl_slope = %f\n",hdr.scl_slope);
   printf("scl_inter = %f\n",hdr.scl_inter);
   printf("slice_end = %d\n",hdr.slice_end);
   printf("slice_code = %d\n",hdr.slice_code);
   printf("xyzt_units = %d\n",hdr.xyzt_units);
   printf("cal_max = %f\n",hdr.cal_max);
   printf("cal_min = %f\n",hdr.cal_min);
   printf("slice_duration = %f\n",hdr.slice_duration);
   printf("toffset = %f\n",hdr.toffset);
   printf("glmax = %d\n",hdr.glmax);
   printf("glmin = %d\n",hdr.glmin);
   printf("descrip = %s\n",hdr.descrip);
   printf("aux_file = %s\n",hdr.aux_file);
   printf("qform_code = %d\n",hdr.qform_code);
   printf("sform_code = %d\n",hdr.sform_code);
   printf("quatern_b = %f\n",hdr.quatern_b);
   printf("quatern_c = %f\n",hdr.quatern_c);
   printf("quatern_d = %f\n",hdr.quatern_d);
   printf("qoffset_x = %f\n",hdr.qoffset_x);
   printf("qoffset_y = %f\n",hdr.qoffset_y);
   printf("qoffset_z = %f\n",hdr.qoffset_z);
   printf("srow_x = %f %f %f %f\n",hdr.srow_x[0],hdr.srow_x[1],hdr.srow_x[2],hdr.srow_x[3]);
   printf("srow_y = %f %f %f %f\n",hdr.srow_y[0],hdr.srow_y[1],hdr.srow_y[2],hdr.srow_y[3]);
   printf("srow_z = %f %f %f %f\n",hdr.srow_z[0],hdr.srow_z[1],hdr.srow_z[2],hdr.srow_z[3]);
   printf("intent_name = %s\n",hdr.intent_name);
   printf("magic = %s\n",hdr.magic);
}

#if 0
// Replaced by the functions in directionCode.cxx
// Input: (x,y,z) a vector defined in RAS system
// Output: One of six charaters {R,L,A,P,S,I}
char directionCode(float x, float y, float z) 
{
   char code;
   float dum;

   dum=x;
   code ='R';

   if( -x > dum )
   {
      dum = -x;
      code ='L';
   }

   if( y > dum )
   {
      dum = y;
      code='A';
   }

   if( -y > dum )
   {
      dum = -y;
      code='P';
   }

   if( z > dum )
   {
      dum = z;
      code='S';
   }

   if( -z > dum )
   {
      dum = -z;
      code='I';
   }

   return(code);
}
#endif

/*
Superceeded by getNiftiImageOrientation() in getNiftiImageOrientation.cxx
void readOrientationFromFile(const char *filename, char *orientation)
{
   FILE *fp;
   nifti_1_header hdr;
   int swapflg=0;

   // ensure that the specified image has either a .hdr or a .nii extension
   if( !checkNiftiFileExtension(filename) )
   {
      errorMessage("Error: The image filename must have a `.hdr' or `.nii' extension.");
   }

   orientation[0]='\0';

   fp = fopen(filename,"r");

   if(fp==NULL)
   {
      errorMessage("Error: I have trouble opening the specified image file.");
   }

   if( fread(&hdr, sizeof(nifti_1_header), 1, fp) != 1 )
   {
      errorMessage("Error: I have trouble reading the specified image file.");
   }

   fclose(fp);

   // looks like ANALYZE 7.5, cannot determine orientation
   if( hdr.magic[0]!='n' ||  (hdr.magic[1]!='+' && hdr.magic[1]!='i') ||  hdr.magic[2]!='1')
   {
      return;
   }

   // if dim[0] is outside range 1..7, then the header information
   // needs to be byte swapped appropriately
   if(hdr.dim[0]<1 || hdr.dim[0]>7) 
   {
      swapflg=1;
   }

   // Here I am only byte swapping the header fields relevant for determining the image orientation
   if(swapflg)
   {
      swapByteOrder( (char *)&(hdr.qform_code), sizeof(short));
      swapByteOrder( (char *)&(hdr.sform_code), sizeof(short));

      swapByteOrder( (char *)&(hdr.quatern_b), sizeof(float));
      swapByteOrder( (char *)&(hdr.quatern_c), sizeof(float));
      swapByteOrder( (char *)&(hdr.quatern_d), sizeof(float));
      swapByteOrder( (char *)&(hdr.qoffset_x), sizeof(float));
      swapByteOrder( (char *)&(hdr.qoffset_y), sizeof(float));
      swapByteOrder( (char *)&(hdr.qoffset_z), sizeof(float));

      for(int i=0; i<4; i++)
      {
         swapByteOrder( (char *)&(hdr.srow_x[i]), sizeof(float));
         swapByteOrder( (char *)&(hdr.srow_y[i]), sizeof(float));
         swapByteOrder( (char *)&(hdr.srow_z[i]), sizeof(float));
      }
   }

   if(hdr.qform_code == 0 && hdr.sform_code == 0) 
   {
      printf("\nThe header of this so called \"NIFTI\" file does not contain orientation information.\n");
      return;
   }

   if(hdr.qform_code > 0 )
   {
      mat44 R;

      R = nifti_quatern_to_mat44(hdr.quatern_b, hdr.quatern_c, hdr.quatern_d, hdr.qoffset_x, hdr.qoffset_y, 
      hdr.qoffset_z, hdr.pixdim[1], hdr.pixdim[2], hdr.pixdim[3], hdr.pixdim[0]);

      orientation[0] = directionCode(R.m[0][0],R.m[1][0],R.m[2][0]);
      orientation[1] = directionCode(R.m[0][1],R.m[1][1],R.m[2][1]);
      orientation[2] = directionCode(R.m[0][2],R.m[1][2],R.m[2][2]);
      orientation[3] = '\0';
   }
   else
   {
      orientation[0] = directionCode(hdr.srow_x[0],hdr.srow_y[0],hdr.srow_z[0]);
      orientation[1] = directionCode(hdr.srow_x[1],hdr.srow_y[1],hdr.srow_z[1]);
      orientation[2] = directionCode(hdr.srow_x[2],hdr.srow_y[2],hdr.srow_z[2]);
      orientation[3] = '\0';
   }

   return;
}
*/

// returns the orientations vectors xvec, yvec, and zvec in NIFTI's RAS system
void readOrientationVectorsFromFile(const char *filename, float *xvec, float *yvec, float *zvec)
{
   FILE *fp;
   nifti_1_header hdr;
   int swapflg=0;

   // ensure that the specified image has either a .hdr or a .nii extension
   if( !checkNiftiFileExtension(filename) )
   {
      errorMessage("Error: The image filename must have a `.hdr' or `.nii' extension.");
   }

   fp = fopen(filename,"r");

   if(fp==NULL)
   {
      errorMessage("Error: I have trouble opening the specified image file.");
   }

   if( fread(&hdr, sizeof(nifti_1_header), 1, fp) != 1 )
   {
      errorMessage("Error: I have trouble reading the specified image file.");
   }

   fclose(fp);

   // looks like ANALYZE 7.5, cannot determine orientation
   if( hdr.magic[0]!='n' ||  (hdr.magic[1]!='+' && hdr.magic[1]!='i') ||  hdr.magic[2]!='1')
   {
      return;
   }

   // if dim[0] is outside range 1..7, then the header information
   // needs to be byte swapped appropriately
   if(hdr.dim[0]<1 || hdr.dim[0]>7) 
   {
      swapflg=1;
   }

   // Here I am only byte swapping the header fields relevant for determining the image orientation
   if(swapflg)
   {
      swapByteOrder( (char *)&(hdr.qform_code), sizeof(short));
      swapByteOrder( (char *)&(hdr.sform_code), sizeof(short));

      swapByteOrder( (char *)&(hdr.quatern_b), sizeof(float));
      swapByteOrder( (char *)&(hdr.quatern_c), sizeof(float));
      swapByteOrder( (char *)&(hdr.quatern_d), sizeof(float));
      swapByteOrder( (char *)&(hdr.qoffset_x), sizeof(float));
      swapByteOrder( (char *)&(hdr.qoffset_y), sizeof(float));
      swapByteOrder( (char *)&(hdr.qoffset_z), sizeof(float));

      for(int i=0; i<4; i++)
      {
         swapByteOrder( (char *)&(hdr.srow_x[i]), sizeof(float));
         swapByteOrder( (char *)&(hdr.srow_y[i]), sizeof(float));
         swapByteOrder( (char *)&(hdr.srow_z[i]), sizeof(float));
      }
   }

   if(hdr.qform_code == 0 && hdr.sform_code == 0) 
   {
      printf("\nThe header of this so called \"NIFTI\" file does not contain orientation information.\n");
      return;
   }

   if(hdr.qform_code > 0 )
   {
      mat44 R;

      R = nifti_quatern_to_mat44(hdr.quatern_b, hdr.quatern_c, hdr.quatern_d, hdr.qoffset_x, hdr.qoffset_y, 
      hdr.qoffset_z, hdr.pixdim[1], hdr.pixdim[2], hdr.pixdim[3], hdr.pixdim[0]);

      xvec[0] = R.m[0][0];
      xvec[1] = R.m[1][0];
      xvec[2] = R.m[2][0];

      yvec[0] = R.m[0][1];
      yvec[1] = R.m[1][1];
      yvec[2] = R.m[2][1];

      zvec[0] = R.m[0][2];
      zvec[1] = R.m[1][2];
      zvec[2] = R.m[2][2];
   }
   else
   {
      xvec[0] = hdr.srow_x[0];
      xvec[1] = hdr.srow_y[0];
      xvec[2] = hdr.srow_z[0];

      yvec[0] = hdr.srow_x[1];
      yvec[1] = hdr.srow_y[1];
      yvec[2] = hdr.srow_z[1];

      zvec[0] = hdr.srow_x[2];
      zvec[1] = hdr.srow_y[2];
      zvec[2] = hdr.srow_z[2];
   }

   normalizeVector(xvec, 3);
   normalizeVector(yvec, 3);
   normalizeVector(zvec, 3);

   return;
}

// extracts the "filename" from the full "path" string
// Example: If path="/home/babak/images/test.nii", then filename="test".
// Note that the extension is not included in the output filename
// Returns 0 on failture and 1 on success
int niftiFilename(char *filename, const char *path)
{
   int i;
   int len;	// length of the path string
   int pos;	// position of the filename

   if( !checkNiftiFileExtension(path) )
   {
      printf("%s does not have a `.hdr' or `.nii' extension\n",path);
      return(0);
   }

   if( not_magical_nifti(path) )
   {
      return(0);
   }

   len=(int)strlen(path);

   if(len<=0)
   {
      printf("Error: unexpected string length for the NIFTI image path, aborting ...\n");
      return(0);
   }

   // finds the position of the first '/' character from right if any
   // returns 0 of no '/' charater is there
   // Examples: path=/sss/yyy would give pos=5
   // path=sss would give pos=0
   // path=/x/ would give pos=3
   i=len-1;
   while( i>=0 && path[i] != '/' )
   {
      i--;
   }
   pos=i+1;

   // copy everything to the right of the first '/' character from right
   // into filename
   // Examples: path=/sss/yyy would give filename=yyy  len=3
   // path=sss would give filename=sss  len=3
   // path=/x/ would give filename=""  len=0
   strcpy(filename,path+pos);
   len=(int)strlen(filename);

   if(len<=0)
   {
      printf("Error: unexpected string length for the NIFTI image filename, aborting ...\n");
      return(0);
   }

   if( len>=2 && filename[len-2]=='g' && filename[len-1]=='z' )
   {
      printf("Sorry but this program currently does not handle gzipped images, aborting ...\n");
      return(0);
   }

   // finds the position of the first '.' character from right if any
   // returns 0 of no '.' charater is there
   // Examples: path=/sss.yyy would give pos=5
   // path=sss would give pos=0
   // path=/x. would give pos=3
   i=len-1;
   while( i>=0 && filename[i] != '.' )
   {
      i--;
   }
   pos=i+1;

   if(pos>0) filename[pos-1]='\0';

   len=(int)strlen(filename);

   if(len<=0)
   {
      printf("Error: unexpected string length for the NIFTI image prefix, aborting ...\n");
      return(0);
   }

   return(1);
}

short *readNiftiImage(const char *filename, DIM *dim, int flg)
{
   FILE *fp;
   nifti_1_header hdr;
   int swapflg=0;
   int nv;
   char *imgname;
   short *im=NULL;

   // ensure that the specified image has either a .hdr or a .nii extension
   if( !checkNiftiFileExtension(filename) )
   {
      errorMessage("Error: The image filename must have a `.hdr' or `.nii' extension.");
   }

   fp = fopen(filename,"r");

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
   if(hdr.dim[0]<1 || hdr.dim[0]>7) 
   {
      swapflg=1;
   }

   if(swapflg)
   {
      swapByteOrder( (char *)&(hdr.sizeof_hdr), sizeof(int));

      swapByteOrder( (char *)&(hdr.dim[0]), sizeof(short));
      for(int i=1; i<=hdr.dim[0]; i++)
      {
         swapByteOrder( (char *)&(hdr.dim[i]), sizeof(short));
      }

      swapByteOrder( (char *)&(hdr.datatype), sizeof(short));

      for(int i=0; i<=hdr.dim[0]; i++)
      {
         swapByteOrder( (char *)&(hdr.pixdim[i]), sizeof(float));
      }

      swapByteOrder( (char *)&(hdr.vox_offset), sizeof(float));
   }

   dim->nx = hdr.dim[1];
   dim->ny = hdr.dim[2];
   dim->nz = hdr.dim[3];
   dim->dx = hdr.pixdim[1];
   dim->dy = hdr.pixdim[2];
   dim->dz = hdr.pixdim[3];
   dim->np = dim->nx*dim->ny;
   dim->nv = dim->np*dim->nz;

   nv = hdr.dim[1]*hdr.dim[2]*hdr.dim[3];

   im = (short *)calloc(nv, sizeof(short));
   
   if(im==NULL)
   {
      errorMessage("Error: Memory allocation problem in readNiftiImage() function.");
   }

   {
      int L;

      L = strlen(filename);

      imgname = (char *)calloc(L+1, 1);

      strcpy(imgname, filename);

      if(imgname[L-3]=='h' && imgname[L-2]=='d' && imgname[L-1]=='r')
      {
         imgname[L-3]='i';
         imgname[L-2]='m';
         imgname[L-1]='g';
      }
   }

   if(flg)
   {
      printf("Reading %s ...\n", imgname);
      printf("\tmatrix size = %d x %d x %d\n", dim->nx, dim->ny, dim->nz);
      printf("\tvoxel size = %5.3f x %5.3f x %5.3f\n", dim->dx, dim->dy, dim->dz);
      printf("\tdata type = %d\n", hdr.datatype);
   }

   // if it is ANALYZE format, we want vox_offset to be zero
   if( hdr.magic[0]!='n' ||  (hdr.magic[1]!='+' && hdr.magic[1]!='i') ||  hdr.magic[2]!='1')
   {
      hdr.vox_offset = 0;
   }

   fp = fopen(imgname,"r");

   if(fp==NULL)
   {
      errorMessage("Error: I have trouble opening the specified image file.");
   }

   if( fseek(fp, (long)hdr.vox_offset, SEEK_SET) != 0 )
   {
      errorMessage("Error: I have trouble reading the specified image file.");
   }

   if( hdr.datatype == DT_SIGNED_SHORT || hdr.datatype == DT_UINT16) 
   {
      if( (int)fread(im, sizeof(short), nv, fp) != nv )
      {
         errorMessage("Error: I have trouble reading the specified image file.");
      }

      if(swapflg)
      {
         swapN( (char *)im, nv*2);
      }

   }
   else if( hdr.datatype == DT_UNSIGNED_CHAR || hdr.datatype == DT_INT8 ) 
   {
      unsigned char *dum;

      dum = (unsigned char *)calloc( nv, sizeof(unsigned char));

      if(dum==NULL)
      {
         errorMessage("Error: Memory allocation problem in readNiftiImage() function.");
      }

      if ( (int)fread(dum, sizeof(unsigned char), nv, fp) != nv )
      {
         errorMessage("Error: I have trouble reading the specified image file.");
      }

      for(int i=0; i<nv; i++) im[i]=(short)dum[i];

      free(dum);
   }
   else
   {
      errorMessage("Error: Sorry, but this program can only handle image types `short' and `unsigned char'.");
   }

   fclose(fp);

  if( hdr.datatype == DT_SIGNED_SHORT || hdr.datatype == DT_UINT16) 
  {
    for(int i=0; i<nv; i++)
      im[i] = (short)( im[i]*hdr.scl_slope + hdr.scl_inter + 0.5);
  }

  return(im);
}

#if 0
char *read_nifti_image(const char *filename, nifti_1_header *hdr)
{
   FILE *fp;
   int swapflg=0;
   int nv;
   char *imgname;
   char *im=NULL;

   // ensure that the specified image has either a .hdr or a .nii extension
   if( !checkNiftiFileExtension(filename) )
   {
      errorMessage("Error: The image filename must have a `.hdr' or `.nii' extension.");
   }

   fp = fopen(filename,"r");

   if(fp==NULL)
   {
      errorMessage("Error: I have trouble opening the specified image file.");
   }

   if( fread(hdr, sizeof(nifti_1_header), 1, fp) != 1 )
   {
      errorMessage("Error: I have trouble reading the specified image file.");
   }

   fclose(fp);

   // if dim[0] is outside range 1..7, then the header information
   // needs to be byte swapped appropriately
   if(hdr->dim[0]<1 || hdr->dim[0]>7) 
   {
      swapflg=1;
   }

   if(swapflg)
   {
      swapByteOrder( (char *)&(hdr->sizeof_hdr), sizeof(int));

      swapByteOrder( (char *)&(hdr->dim[0]), sizeof(short));
      for(int i=1; i<=hdr->dim[0]; i++)
      {
         swapByteOrder( (char *)&(hdr->dim[i]), sizeof(short));
      }

      swapByteOrder( (char *)&(hdr->datatype), sizeof(short));

      for(int i=0; i<=hdr->dim[0]; i++)
      {
         swapByteOrder( (char *)&(hdr->pixdim[i]), sizeof(float));
      }

      swapByteOrder( (char *)&(hdr->vox_offset), sizeof(float));
   }

   nv = hdr->dim[1]*hdr->dim[2]*hdr->dim[3];

   im = (char *)calloc(nv*hdr->bitpix/8, sizeof(char));
   
   if(im==NULL)
   {
      errorMessage("Error: Memory allocation problem in readNiftiImage() function.");
   }

   {
      int L;

      L = strlen(filename);

      imgname = (char *)calloc(L+1, 1);

      strcpy(imgname, filename);

      if(imgname[L-3]=='h' && imgname[L-2]=='d' && imgname[L-1]=='r')
      {
         imgname[L-3]='i';
         imgname[L-2]='m';
         imgname[L-1]='g';
      }
   }

   // if it is ANALYZE format, we want vox_offset to be zero
   if( hdr->magic[0]!='n' ||  (hdr->magic[1]!='+' && hdr->magic[1]!='i') ||  hdr->magic[2]!='1')
   {
      hdr->vox_offset = 0;
   }

   fp = fopen(imgname,"r");

   if(fp==NULL)
   {
      errorMessage("Error: I have trouble opening the specified image file.");
   }

   if( fseek(fp, (long)hdr->vox_offset, SEEK_SET) != 0 )
   {
      errorMessage("Error: I have trouble reading the specified image file.");
   }

   if( hdr->datatype == DT_SIGNED_SHORT || hdr->datatype == DT_UINT16) 
   {
      if( fread(im, sizeof(short), nv, fp) != nv )
      {
         errorMessage("Error: I have trouble reading the specified image file.");
      }

      if(swapflg)
      {
         swapN( im, nv*2);
      }
   }
   else if( hdr->datatype == DT_UNSIGNED_CHAR || hdr->datatype == DT_INT8 ) 
   {
      if (fread(im, sizeof(unsigned char), nv, fp) != nv )
      {
         errorMessage("Error: I have trouble reading the specified image file.");
      }
   }
   else
   {
      errorMessage("Error: Sorry, but this program can only handle image types `short' and `unsigned char'.");
   }

   fclose(fp);

   return(im);
}
#endif

char *read_nifti_image(const char *filename, nifti_1_header *hdr)
{
   FILE *fp;
   int swapflg=0;
   int datasize;
   char *imgname;
   char *im=NULL;
   long voxeloffset;

/*
   // ensure that the specified image has either a .hdr or a .nii extension
   if( !checkNiftiFileExtension(filename) )
   {
      printf("\nread_nifti_image(): %s does not have a `.hdr' or `.nii' extension.\n\n",filename);
      return(NULL);
   }
*/

   fp = fopen(filename,"r");

   if(fp==NULL)
   {
      printf("\nread_nifti_image(): Opening %s for reading failed.\n\n",filename);
      return(NULL);
   }

   if( fread(hdr, sizeof(nifti_1_header), 1, fp) != 1 )
   {
      printf("\nread_nifti_image(): Reading %s failed.\n\n",filename);
      return(NULL);
   }

   fclose(fp);

   if( hdr->magic[0]!='n' ||  (hdr->magic[1]!='+' && hdr->magic[1]!='i') ||  hdr->magic[2]!='1')
   {
      printf("\nread_nifti_image(): %s does not have the NIFTI magic.\n\n",filename);
      return(NULL);
   }

   // if dim[0] is outside range 1..7, then the header information
   // needs to be byte swapped appropriately
   if(hdr->dim[0]<1 || hdr->dim[0]>7) 
   {
      swapflg=1;
      swapniftiheader(hdr);
   }

   {
      int L;

      L = strlen(filename);

      imgname = (char *)calloc(L+1, 1);

      strcpy(imgname, filename);

      if(imgname[L-3]=='h' && imgname[L-2]=='d' && imgname[L-1]=='r')
      {
         imgname[L-3]='i';
         imgname[L-2]='m';
         imgname[L-1]='g';

         voxeloffset = 0;
      }
      else
      {
         voxeloffset = (long)hdr->vox_offset;
      }
   }

   fp = fopen(imgname,"r");

   if(fp==NULL)
   {
      printf("\nread_nifti_image(): Opening %s for reading failed.\n\n",imgname);
      free(imgname);
      return(NULL);
   }

   if( fseek(fp, voxeloffset, SEEK_SET) != 0 )
   {
      printf("\nread_nifti_image(): Reading %s failed.\n\n",imgname);
      fclose(fp);
      free(imgname);
      return(NULL);
   }

   datasize = 1;
   for(int i=1; i<=hdr->dim[0]; i++)
   {
      datasize *= hdr->dim[i];
   }
   datasize *= (hdr->bitpix/8);

   im = (char *)calloc(datasize, sizeof(char));
   if(im==NULL)
   {
      printf("\nread_nifti_image(): Memory allocation problem.\n\n");
      fclose(fp);
      free(imgname);
      return(NULL);
   }

   if( (int)fread(im, 1, datasize, fp) != datasize )
   {
      printf("\nread_nifti_image(): Reading %s failed.\n\n",imgname);
      free(im);
      fclose(fp);
      free(imgname);
      return(NULL);
   }

   fclose(fp);

   if(swapflg)
   {
      if( hdr->datatype == DT_SIGNED_SHORT || hdr->datatype == DT_UINT16) 
      { 
         swapN(im, datasize);
      }

      if( hdr->datatype == DT_FLOAT ) 
      { 
         swap_float_array( (float *)im, datasize/sizeof(float));
      }

      if( hdr->datatype == DT_DOUBLE) 
      { 
         swap_double_array( (float8 *)im, datasize/sizeof(float8));
      }

      if( hdr->datatype == DT_SIGNED_INT ) 
      { 
         swap_int_array( (int *)im, datasize/sizeof(int));
      }
   }

  int nv=1;
  for(int i=1; i<=hdr->dim[0]; i++)
  {
    nv *= hdr->dim[i];
  }

  if( hdr->datatype == DT_SIGNED_SHORT || hdr->datatype == DT_UINT16) 
  {
    short *tmp;

    tmp = (short *)im;
    for(int i=0; i<nv; i++)
      tmp[i] = (short)( tmp[i]*hdr->scl_slope + hdr->scl_inter + 0.5);    
  }

  free(imgname);
  return(im);
}
