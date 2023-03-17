#include <niftiimage.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// destructor
NIFTIIMAGE::~NIFTIIMAGE()
{
   if(extension.edata != NULL)
      free(extension.edata);

   if(data != NULL)
      free(data);

   //printf("destructor called\n");
}

NIFTIIMAGE::NIFTIIMAGE(int nx, int ny, int nz, float dx, float dy, float dz, short datatype)   // constructor
{
   hdr.sizeof_hdr = 348;
   hdr.data_type[0] = '\0';
   hdr.db_name[0] = '\0';
   hdr.extents = 0;
   hdr.session_error = 0;
   hdr.regular = 0;

   // 57 or 111001 sets freq_dim=1, phase_dim=2, slice_dim=3
   hdr.dim_info = 0;

   hdr.dim[0] = 4;
   hdr.dim[1] = nx;
   hdr.dim[2] = ny;
   hdr.dim[3] = nz;
   hdr.dim[4] = 1;

   hdr.intent_p1 = 0.0;
   hdr.intent_p2 = 0.0;
   hdr.intent_p3 = 0.0;

   hdr.intent_code = NIFTI_INTENT_NONE;
   hdr.datatype = datatype;

   if(datatype == NIFTI_TYPE_UINT8)
      hdr.bitpix = 8;
   else if(datatype == NIFTI_TYPE_INT16)
      hdr.bitpix = 16;
   else if(datatype == NIFTI_TYPE_INT32)
      hdr.bitpix = 32;
   else if(datatype == NIFTI_TYPE_FLOAT32)
      hdr.bitpix = 32;

   hdr.slice_start = 0;

   hdr.pixdim[0] = 1.0;
   hdr.pixdim[1] = dx;
   hdr.pixdim[2] = dy;
   hdr.pixdim[3] = dz;

   hdr.vox_offset = 352.0;
   hdr.scl_slope = 1.0;
   hdr.scl_inter = 0.0;
   hdr.slice_end = 0;
   hdr.slice_code = 0;
   hdr.xyzt_units = 0;
   hdr.cal_max = 0.0;
   hdr.cal_min = 0.0;
   hdr.slice_duration= 0.0;
   hdr.toffset= 0.0;
   hdr.glmax = 0;
   hdr.glmin = 0;
   strcpy(hdr.descrip, "Created by ART software.");
   strcpy(hdr.aux_file, "");
   hdr.qform_code = 0;
   hdr.sform_code = 0;

   hdr.quatern_b = 0.0;
   hdr.quatern_c = 0.0;
   hdr.quatern_d = 0.0;
   hdr.qoffset_x = 0.0;
   hdr.qoffset_y = 0.0;
   hdr.qoffset_z = 0.0;

   for(int i=0; i<4; i++) 
   {
      hdr.srow_x[i] = hdr.srow_y[i] = hdr.srow_z[i] = 0.0;
   }

   strcpy(hdr.intent_name, "");

   strcpy(hdr.magic, "n+1");

   ext.extension[0]=0;
   ext.extension[1]=0;
   ext.extension[2]=0;
   ext.extension[3]=0;

   extension.esize = 0;
   extension.ecode = 0;
   extension.edata = NULL; // important to initialize this to NULL

   data = (char *)calloc(nx*ny*nz*hdr.bitpix/8, 1);
}

NIFTIIMAGE::NIFTIIMAGE(int nx, int ny, int nz, int nt, float dx, float dy, float dz, short datatype)   // constructor
{
   hdr.sizeof_hdr = 348;
   hdr.data_type[0] = '\0';
   hdr.db_name[0] = '\0';
   hdr.extents = 0;
   hdr.session_error = 0;
   hdr.regular = 0;

   // 57 or 111001 sets freq_dim=1, phase_dim=2, slice_dim=3
   hdr.dim_info = 0;

   hdr.dim[0] = 4;
   hdr.dim[1] = nx;
   hdr.dim[2] = ny;
   hdr.dim[3] = nz;
   hdr.dim[4] = nt;

   hdr.intent_p1 = 0.0;
   hdr.intent_p2 = 0.0;
   hdr.intent_p3 = 0.0;

   hdr.intent_code = NIFTI_INTENT_NONE;
   hdr.datatype = datatype;

   if(datatype == NIFTI_TYPE_UINT8)
      hdr.bitpix = 8;
   else if(datatype == NIFTI_TYPE_INT16)
      hdr.bitpix = 16;
   else if(datatype == NIFTI_TYPE_INT32)
      hdr.bitpix = 32;
   else if(datatype == NIFTI_TYPE_FLOAT32)
      hdr.bitpix = 32;

   hdr.slice_start = 0;

   hdr.pixdim[0] = 1.0;
   hdr.pixdim[1] = dx;
   hdr.pixdim[2] = dy;
   hdr.pixdim[3] = dz;

   hdr.vox_offset = 352.0;
   hdr.scl_slope = 1.0;
   hdr.scl_inter = 0.0;
   hdr.slice_end = 0;
   hdr.slice_code = 0;
   hdr.xyzt_units = 0;
   hdr.cal_max = 0.0;
   hdr.cal_min = 0.0;
   hdr.slice_duration= 0.0;
   hdr.toffset= 0.0;
   hdr.glmax = 0;
   hdr.glmin = 0;
   strcpy(hdr.descrip, "Created by ART software.");
   strcpy(hdr.aux_file, "");
   hdr.qform_code = 0;
   hdr.sform_code = 0;

   hdr.quatern_b = 0.0;
   hdr.quatern_c = 0.0;
   hdr.quatern_d = 0.0;
   hdr.qoffset_x = 0.0;
   hdr.qoffset_y = 0.0;
   hdr.qoffset_z = 0.0;

   for(int i=0; i<4; i++) 
   {
      hdr.srow_x[i] = hdr.srow_y[i] = hdr.srow_z[i] = 0.0;
   }

   strcpy(hdr.intent_name, "");

   strcpy(hdr.magic, "n+1");

   ext.extension[0]=0;
   ext.extension[1]=0;
   ext.extension[2]=0;
   ext.extension[3]=0;

   extension.esize = 0;
   extension.ecode = 0;
   extension.edata = NULL; // important to initialize this to NULL

   data = (char *)calloc(nx*ny*nz*nt*hdr.bitpix/8, 1);
}

// constructor
NIFTIIMAGE::NIFTIIMAGE()
{
   hdr.sizeof_hdr = 348;
   hdr.data_type[0] = '\0';
   hdr.db_name[0] = '\0';
   hdr.extents = 0;
   hdr.session_error = 0;
   hdr.regular = 0;

   // 57 or 111001 sets freq_dim=1, phase_dim=2, slice_dim=3
   hdr.dim_info = 0;

   hdr.dim[0] = 1;
   hdr.dim[1] = 0;

   hdr.intent_p1 = 0.0;
   hdr.intent_p2 = 0.0;
   hdr.intent_p3 = 0.0;

   hdr.intent_code = NIFTI_INTENT_NONE;
   hdr.datatype = DT_NONE;
   hdr.bitpix = 0;
   hdr.slice_start = 0;

   hdr.pixdim[0] = 1.0;
   hdr.pixdim[1] = 0.0;

   hdr.vox_offset = 352.0;
   hdr.scl_slope = 1.0;
   hdr.scl_inter = 0.0;
   hdr.slice_end = 0;
   hdr.slice_code = 0;
   hdr.xyzt_units = 0;
   hdr.cal_max = 0.0;
   hdr.cal_min = 0.0;
   hdr.slice_duration= 0.0;
   hdr.toffset= 0.0;
   hdr.glmax = 0;
   hdr.glmin = 0;
   strcpy(hdr.descrip, "Created by ART software.");
   strcpy(hdr.aux_file, "");
   hdr.qform_code = 0;
   hdr.sform_code = 0;

   hdr.quatern_b = 0.0;
   hdr.quatern_c = 0.0;
   hdr.quatern_d = 0.0;
   hdr.qoffset_x = 0.0;
   hdr.qoffset_y = 0.0;
   hdr.qoffset_z = 0.0;

   for(int i=0; i<4; i++) 
   {
      hdr.srow_x[i] = hdr.srow_y[i] = hdr.srow_z[i] = 0.0;
   }

   strcpy(hdr.intent_name, "");

   strcpy(hdr.magic, "n+1");

   ext.extension[0]=0;
   ext.extension[1]=0;
   ext.extension[2]=0;
   ext.extension[3]=0;

   extension.esize = 0;
   extension.ecode = 0;
   extension.edata = NULL; // important to initialize this to NULL

   data = NULL;
}

void NIFTIIMAGE::printheader()
{
   printf("nifti_1_header:\n");
   printf("sizeof_hdr = %d (Must be 348)\n", hdr.sizeof_hdr);
   printf("data_type[10] = %s (UNUSED)\n", hdr.data_type);
   printf("db_name[18] = %s (UNUSED)\n", hdr.db_name);
   printf("extents = %d (UNUSED)\n", hdr.extents);
   printf("session_error = %d (UNUSED)\n", hdr.session_error);
   printf("regular = %d (UNUSED)\n", hdr.regular);

   printf("dim_info (MRI slice ordering):\n");
   printf("\tfreq_dim = %d\n", DIM_INFO_TO_FREQ_DIM( hdr.dim_info ) );   
   printf("\tphase_dim = %d\n", DIM_INFO_TO_PHASE_DIM( hdr.dim_info ) );
   printf("\tslice_dim = %d\n", DIM_INFO_TO_SLICE_DIM( hdr.dim_info ) );

   printf("dim[8] (Data array dimensions):\n");
   printf("\tdim[0] = %d (number of dimensions)\n", hdr.dim[0]);   
   for(int i=1; i<=hdr.dim[0]; i++)
   {
      printf("\tdim[%d] = %d\n", i, hdr.dim[i]);   
   }

   printf("intent_p1 = %f (1st intent parameter)\n", hdr.intent_p1 );
   printf("intent_p2 = %f (2nd intent parameter)\n", hdr.intent_p2 );
   printf("intent_p3 = %f (3rd intent parameter)\n", hdr.intent_p3 );

   printf("intent_code = %d (NIFTI_INTENT_* code)\n", hdr.intent_code);
   printf("datatype = %d (Defines data type!)\n", hdr.datatype);
   printf("bitpix = %d (Number bits/voxel)\n", hdr.bitpix);
   printf("slice_start = %d (First slice index)\n", hdr.slice_start);

   printf("pixdim[8] (Grid spacings):\n");
   for(int i=1; i<=hdr.dim[0]; i++)
   {
      printf("\tpixdim[%d] = %f\n", i, hdr.pixdim[i]);   
   }

   printf("vox_offset = %f (Offset into .nii file)\n", hdr.vox_offset);
   printf("scl_slope = %f (Data scaling: slope)\n", hdr.scl_slope);
   printf("scl_inter = %f (Data scaling: offset)\n", hdr.scl_inter);
   printf("slice_end = %d (Last slice index)\n", hdr.slice_end);
   printf("slice_code = %d (Slice timing order)\n", hdr.slice_code);
   printf("xyzt_units = %d (Units of pixdim[1..4])\n", hdr.xyzt_units);
   printf("cal_max = %f (Max display intensity)\n", hdr.cal_max);
   printf("cal_min = %f (Min display intensity)\n", hdr.cal_min);
   printf("slice_duration = %f (Time for 1 slice)\n", hdr.slice_duration);
   printf("toffset = %f (Time axis shift)\n", hdr.toffset);
   printf("glmax = %d (UNUSED)\n", hdr.glmax);
   printf("glmin = %d (UNUSED)\n", hdr.glmin);

   printf("descrip[80] = %s\n", hdr.descrip);
   printf("aux_file[24] = %s (auxiliary filename)\n", hdr.aux_file);

   printf("qform_code = %d (NIFTI_XFORM_* code)\n", hdr.qform_code);
   printf("sform_code = %d (NIFTI_XFORM_* code)\n", hdr.sform_code);

   printf("quatern_b = %f (Quaternion b param)\n", hdr.quatern_b);
   printf("quatern_c = %f (Quaternion c param)\n", hdr.quatern_c);
   printf("quatern_d = %f (Quaternion d param)\n", hdr.quatern_d);

   printf("qoffset_x = %f (Quaternion x shift)\n", hdr.qoffset_x);
   printf("qoffset_y = %f (Quaternion y shift)\n", hdr.qoffset_y);
   printf("qoffset_z = %f (Quaternion z shift)\n", hdr.qoffset_z);

   printf("srow_x[4] (1st row affine transform):\n");
   printf("\t%f\t%f\t%f\t%f\n",hdr.srow_x[0], hdr.srow_x[1], hdr.srow_x[2], hdr.srow_x[3]);

   printf("srow_y[4] (2nd row affine transform):\n");
   printf("\t%f\t%f\t%f\t%f\n",hdr.srow_y[0], hdr.srow_y[1], hdr.srow_y[2], hdr.srow_y[3]);

   printf("srow_z[4] (3rd row affine transform):\n");
   printf("\t%f\t%f\t%f\t%f\n",hdr.srow_z[0], hdr.srow_z[1], hdr.srow_z[2], hdr.srow_z[3]);

   printf("intent_name[16] = %s ('name' or meaning of data)\n", hdr.intent_name);

   printf("magic[4] = %s (Must be \"nil\" or \"n+1\")\n", hdr.magic);

   printf("nifti1_extender:\n");
   printf("\t%d\t%d\t%d\t%d\n", ext.extension[0], ext.extension[1], ext.extension[2], ext.extension[3]);

   if(ext.extension[0] == 1)
   {
      printf("esize = %d (size of extension, in bytes (must be multiple of 16))\n", extension.esize);
      printf("ecode = %d (extension code, one of the NIFTI_ECODE_ values)\n", extension.ecode);
   }

   return;
}

void NIFTIIMAGE::setheader(nifti_1_header newhdr)
{
   hdr = newhdr;

   return;
}

nifti_1_header NIFTIIMAGE::getheader()
{
   return(hdr);
}

char *NIFTIIMAGE::getdata()
{
   return(data);
}

void NIFTIIMAGE::read(const char *filename)
{
   FILE *fp;
   int datasize;

   fp = fopen(filename,"r");

   if(fp==NULL)
   {
      printf("\n\nWARNING: could not open %s for reading!\n\n", filename);
      return;
   }

   if( fread(&hdr, sizeof(nifti_1_header), 1, fp) != 1)
   {
      printf("Failed to read \"hdr\" variable.\n");
   }

   if( fread(&ext, sizeof(nifti1_extender), 1, fp) != 1)
   {
      printf("Failed to read \"ext\" variable.\n");
   }

   if(ext.extension[0] == 1)
   {
      if( fread(&(extension.esize), sizeof(int), 1, fp) != 1)
      {
         printf("Failed to read \"extension.esize\" variable.\n");
      }

      if( fread(&(extension.ecode), sizeof(int), 1, fp) != 1)
      {
         printf("Failed to read \"extension.ecode\" variable.\n");
      }

      extension.edata = (char *)calloc(extension.esize-2*sizeof(int), 1);

      if( fread(extension.edata, 1, extension.esize-2*sizeof(int), fp) != extension.esize-2*sizeof(int) )
      {
         printf("Failed to read \"extension.edata\" variable.\n");
      }
   }

   datasize = 1;
   for(int i=1; i<=hdr.dim[0]; i++)
   {
      datasize *= hdr.dim[i];
   }
   datasize *= hdr.bitpix/8;

   data = (char *)calloc(datasize, 1);
   if( fread(data, 1, datasize, fp) != (size_t)datasize )
   {
      printf("Failed to read \"data\" variable.\n");
   }

   fclose(fp);

   return;
}

void NIFTIIMAGE::readheader(const char *filename)
{
   FILE *fp;
   //int datasize;

   fp = fopen(filename,"r");

   if(fp==NULL)
   {
      printf("\n\nWARNING: could not open %s for reading!\n\n", filename);
      return;
   }

   if( fread(&hdr, sizeof(nifti_1_header), 1, fp) != 1)
   {
      printf("Failed to read \"hdr\" variable.\n");
   }

   if( fread(&ext, sizeof(nifti1_extender), 1, fp) != 1)
   {
      printf("Failed to read \"ext\" variable.\n");
   }

   if(ext.extension[0] == 1)
   {
      if( fread(&(extension.esize), sizeof(int), 1, fp) != 1)
      {
         printf("Failed to read \"extension.esize\" variable.\n");
      }

      if( fread(&(extension.ecode), sizeof(int), 1, fp) != 1)
      {
         printf("Failed to read \"extension.ecode\" variable.\n");
      }

      extension.edata = (char *)calloc(extension.esize-2*sizeof(int), 1);
      if( fread(extension.edata, 1, extension.esize-2*sizeof(int), fp) != extension.esize-2*sizeof(int) )
      {
         printf("Failed to read \"extension.edata\" variable.\n");
      }
   }

   fclose(fp);

   return;
}

int NIFTIIMAGE::nv()
{
   int number_of_voxels;

   number_of_voxels = 1;

   for(int i=1; i<=hdr.dim[0]; i++)
   {
      number_of_voxels *= hdr.dim[i];
   }

   return(number_of_voxels);
}

int NIFTIIMAGE::datasize()
{
   return( nv() * hdr.bitpix / 8);
}

int NIFTIIMAGE::nx()
{
   return(hdr.dim[1]);
}

int NIFTIIMAGE::ny()
{
   return(hdr.dim[2]);
}

int NIFTIIMAGE::nz()
{
   return(hdr.dim[3]);
}

float NIFTIIMAGE::dx()
{
   return(hdr.pixdim[1]);
}

float NIFTIIMAGE::dy()
{
   return(hdr.pixdim[2]);
}

float NIFTIIMAGE::dz()
{
   return(hdr.pixdim[3]);
}

short NIFTIIMAGE::datatype()
{
   return(hdr.datatype);
}

void NIFTIIMAGE::write(const char *filename)
{
   FILE *fp;
   int datasize;

   fp = fopen(filename,"w");

   if(fp==NULL)
   {
      printf("\n\nWARNING: could not open %s for writing!\n\n", filename);
      return;
   }

   fwrite(&hdr, sizeof(nifti_1_header), 1, fp);
   fwrite(&ext, sizeof(nifti1_extender), 1, fp);

   if(ext.extension[0] == 1)
   {
      fwrite(&(extension.esize), sizeof(int), 1, fp);
      fwrite(&(extension.ecode), sizeof(int), 1, fp);
      fwrite(extension.edata, 1, extension.esize-2*sizeof(int), fp);
   }

   datasize = 1;
   for(int i=1; i<=hdr.dim[0]; i++)
   {
      datasize *= hdr.dim[i];
   }
   datasize *= hdr.bitpix/8;

   fwrite(data, 1, datasize, fp);

   fclose(fp);

   return;
}
