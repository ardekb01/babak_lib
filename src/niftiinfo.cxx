#include <ftw.h> // ftw()
#include <sys/types.h> // stat()
#include <sys/stat.h> // stat()
#include <unistd.h>	// stat()
#include <stdio.h> // printf()
#include <string.h> // strcmp()
#include <stdlib.h> // atoi()
#include <math.h> // sqrt()
#include <nifti1_io.h>

#include "babak_lib.h"

#define IMPLICIT_LITTLE_ENDIAN 0
#define EXPLICIT_LITTLE_ENDIAN 1
#define EXPLICIT_BIG_ENDIAN 2

int main(int argc, char **argv)
{
   FILE *fp;
   nifti_1_header hdr;
   nifti1_extender hdr_extender;
   nifti1_extension hdr_extension;
   dicominfo extension_data;
   char orientation_code[4];

   int swapflg=0;

   if(argc==1)
   {
      exit(0);
   }

   fp = fopen(argv[1],"r");
   fread(&hdr, sizeof(nifti_1_header), 1, fp);
   fread(&hdr_extender, sizeof(nifti1_extender), 1, fp);
   if(hdr_extender.extension[0] == 1)
   {
      fread(&hdr_extension.esize, 4, 1, fp);
      fread(&hdr_extension.ecode, 4, 1, fp);

      if(hdr_extension.ecode == 1022)
      {
         fread(&extension_data, sizeof(dicominfo), 1, fp);
      }
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
      swapByteOrder( (char *)&(hdr.bitpix), sizeof(short));

      for(int i=0; i<=hdr.dim[0]; i++)
      {
         swapByteOrder( (char *)&(hdr.pixdim[i]), sizeof(float));
      }

      swapByteOrder( (char *)&(hdr.vox_offset), sizeof(float));
      swapByteOrder( (char *)&(hdr.scl_slope), sizeof(float));
      swapByteOrder( (char *)&(hdr.scl_inter), sizeof(float));
      swapByteOrder( (char *)&(hdr.cal_max), sizeof(float));
      swapByteOrder( (char *)&(hdr.cal_min), sizeof(float));
      swapByteOrder( (char *)&(hdr.toffset), sizeof(float));
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

   printf("sizeof_hdr = %d\n", hdr.sizeof_hdr);
   printf("freq_dim = %d\n", DIM_INFO_TO_FREQ_DIM(hdr.dim_info)); 
   printf("phase_dim = %d\n", DIM_INFO_TO_PHASE_DIM(hdr.dim_info)); 
   printf("slice_dim = %d\n", DIM_INFO_TO_SLICE_DIM(hdr.dim_info)); 

   for(int i=0; i<=hdr.dim[0]; i++)
   {
      printf("dim[%d] = %d\n",i , hdr.dim[i]);
   }

   printf("datatype = %d\n",hdr.datatype);
   printf("bitpix= %d\n",hdr.bitpix);

   for(int i=0; i<=hdr.dim[0]; i++)
   {
      printf("pixdim[%d] = %f\n",i , hdr.pixdim[i]);
   }

   printf("vox_offset = %d\n",(int)hdr.vox_offset);
   printf("scl_slope = %f\n",hdr.scl_slope);
   printf("scl_inter = %f\n",hdr.scl_inter);

   printf("xyzt_units = %d\n",hdr.xyzt_units);
   printf("cal_max = %f\n",hdr.cal_max);
   printf("cal_min = %f\n",hdr.cal_min);
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

   printf("srow_x = (%f, %f, %f, %f)\n",hdr.srow_x[0], hdr.srow_x[1], hdr.srow_x[2], hdr.srow_x[3]);
   printf("srow_y = (%f, %f, %f, %f)\n",hdr.srow_y[0], hdr.srow_y[1], hdr.srow_y[2], hdr.srow_y[3]);
   printf("srow_z = (%f, %f, %f, %f)\n",hdr.srow_z[0], hdr.srow_z[1], hdr.srow_z[2], hdr.srow_z[3]);
   printf("indent_name = %s\n",hdr.intent_name);
   printf("magic = %s\n",hdr.magic);

   if(hdr_extender.extension[0] == 1 && hdr_extension.ecode == 1022)
   {
         printf("Patient name:  %s\n", extension_data.patientName);
         printf("DOB: %s\n", extension_data.DOB);
         printf("Patient ID: %s\n", extension_data.patientID);
         printf("Study date: %s\n", extension_data.studyDate);
         printf("TE = %s\n", extension_data.TE);
         printf("TR = %s\n", extension_data.TR);
         if(extension_data.TI[0] != '\0')
            printf("TI = %s\n", extension_data.TI);
         printf("ETL = %s\n", extension_data.ETL);
         if(extension_data.NEX[0] != '\0')
            printf("NEX = %s\n", extension_data.NEX);
         if(extension_data.flipAngle[0] != '\0')
            printf("Flip angle = %s\n", extension_data.flipAngle);
         if(extension_data.bandwidth[0] != '\0')
            printf("Pixel bandwidth = %s Hz/Pixel\n", extension_data.bandwidth);
         if(extension_data.freq[0] != '\0')
            printf("Imaging frequency = %s MHz\n", extension_data.freq);
         if(extension_data.acquisitionMatrix[0] != 0)
            printf("Frequency rows = %d\n", extension_data.acquisitionMatrix[0]);
         if(extension_data.acquisitionMatrix[1] != 0)
            printf("Frequency columns = %d\n", extension_data.acquisitionMatrix[1]);
         if(extension_data.acquisitionMatrix[2] != 0)
            printf("Phase rows = %d\n", extension_data.acquisitionMatrix[2]);
         if(extension_data.acquisitionMatrix[3] != 0)
            printf("Phase columns = %d\n", extension_data.acquisitionMatrix[3]);
         if(extension_data.phaseFOV[0] != '\0')
            printf("Phase FOV = %s\n", extension_data.phaseFOV);
   }

   if(hdr.qform_code == 0 && hdr.sform_code == 0)
   {
      printf("No particular spatial orientation is specified in the header (qform_code=0 and sform_code=0).\n");
   }

   if(hdr.qform_code > 0 )
   {
      mat44 R;

      R = nifti_quatern_to_mat44(hdr.quatern_b, hdr.quatern_c, hdr.quatern_d, hdr.qoffset_x, hdr.qoffset_y, 
      hdr.qoffset_z, hdr.pixdim[1], hdr.pixdim[2], hdr.pixdim[3], hdr.pixdim[0]);

      orientation_code[0] = directionCode(R.m[0][0],R.m[1][0],R.m[2][0]);
      orientation_code[1] = directionCode(R.m[0][1],R.m[1][1],R.m[2][1]);
      orientation_code[2] = directionCode(R.m[0][2],R.m[1][2],R.m[2][2]);
      orientation_code[3] = '\0';

      printf("Image orientation: %s\n", orientation_code);
   }
   else if(hdr.sform_code > 0 )
   {
      orientation_code[0] = directionCode(hdr.srow_x[0],hdr.srow_y[0],hdr.srow_z[0]);
      orientation_code[1] = directionCode(hdr.srow_x[1],hdr.srow_y[1],hdr.srow_z[1]);
      orientation_code[2] = directionCode(hdr.srow_x[2],hdr.srow_y[2],hdr.srow_z[2]);
      orientation_code[3] = '\0';

      printf("Image orientation: %s\n", orientation_code);
   }
}
