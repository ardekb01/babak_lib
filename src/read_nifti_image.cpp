#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "babak_lib.h"
#include "nifti1.h"
#include "swap.h"
#include "minmax.h"

char *read_nifti_image(const char *filename, nifti_1_header *hdr)
{
   FILE *fp;
   int swapflg=0;
   int datasize;
   char *imgname;
   char *im=NULL;
   long voxeloffset;

/*
   // ensure that the specified image has either .nii extension
   if( check_nifti_file_extension(filename) == false )
   {
      printf("\nread_nifti_image(): %s does not have `.nii' extension.\n\n",filename);
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


  float *floatim;
  float max=0.0;
  floatim = (float *)calloc(nv, sizeof(float));

  if( hdr->datatype == DT_SIGNED_SHORT) 
  {
    short *tmp;
    tmp = (short *)im;

    for(int i=0; i<nv; i++)
      floatim[i] = tmp[i]*hdr->scl_slope + hdr->scl_inter + .5;    

    arraymax(floatim, nv, max);
    if(max>32767) for(int i=0; i<nv; i++) floatim[i] *= (32767/max);    

    for(int i=0; i<nv; i++)
      tmp[i] = (short)(floatim[i]);    
  }

  if( hdr->datatype == DT_UINT16) 
  {
    unsigned short *tmp;
    tmp = (unsigned short *)im;

    for(int i=0; i<nv; i++)
      floatim[i] = tmp[i]*hdr->scl_slope + hdr->scl_inter + .5;    

    arraymax(floatim, nv, max);
    if(max>32767) for(int i=0; i<nv; i++) floatim[i] *= (32767/max);    

    for(int i=0; i<nv; i++)
      tmp[i] = (unsigned short)(floatim[i]);    
  }

  free(imgname);
  free(floatim);
  return(im);
}
