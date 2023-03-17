// Acronyms
// VR = Value Representation
// IS = Integer String
// GN = Group Number 
// EN = Element Number 
// VL = Value Length

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \
*  program: dicomIO.c                                        *
*  Copyright 2005 by Babak A. Ardekani                       *
*  ALL RIGHTS RESERVED.  No part of this program may be      *
*  used, transferred, or modified by any means without       *
*  prior written permission.  Making copies of any part      *
*  of this program for any purpose is a violation of         *
*  copyright laws.                                           *
\ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define _dicomIO

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <babak_lib.h>

#define IMPLICIT_LITTLE_ENDIAN 0
#define EXPLICIT_LITTLE_ENDIAN 1
#define EXPLICIT_BIG_ENDIAN 2

int readNEX(const char *file, int *nex);
int readSliceThickness(const char *file, float *sliceThickness);
int readTR(const char *file, int *TR);
int readMetaFileInfo(const char *file, long *offset, int *syntax);
int readTag(const char *file, unsigned short iGN, unsigned short iEN, long byteOffset, int transferSyntax,
unsigned long *oVL, char *oV, long *valueOffset);

///////////////////////////////////////////////////////////////////////////////////////////////
// This function returns a value of 1 if the input "file" is recognized as a DICOM image, and
// 0 otherwise.
// Every DICOM file must contain a 128 byte File Preamble following by a 4 byte prefix 
// containing the character string "DICM" encoded as uppercase characters.
// Ref: PS 3.10-1996 Page 12.
///////////////////////////////////////////////////////////////////////////////////////////////
int dicomFormat(const char *file)
{
   FILE *fp;
   char dicomPrefix[4];

   // open the file for reading
   fp = fopen(file,"r");

   if(fp==NULL) 
   {
      printf("\ndicomFormat() error: cannot open %s\n\n",file);
      return(0);
   }

   // Skip over the 128 byte File Preamble 
   if( fseek(fp,128,SEEK_SET) == -1) 
   {
      printf("\ndicomFormat() error: 128-byte File Preamble fseek() failed for %s\n\n",file);
      fclose(fp);
      return(0);
   }

   // The four byte DICOM Prefix must contain the character string "DICM"
   if( fread(dicomPrefix,1,4,fp) != 4) 
   {
      printf("\ndicomFormat() error: DICOM Prefix fread() failed for %s\n\n",file);
      fclose(fp);
      return(0);
   }

   fclose(fp);

   if(dicomPrefix[0]!='D' || dicomPrefix[1]!='I' || dicomPrefix[2]!='C' || dicomPrefix[3]!='M')
   {
      printf("\ndicomFormat() error: %s does not contain the required \"DICM\" prefix.\n\n", file);
      return(0);
   }

   return(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////

// Reads the DICOM File Meta Information.
// By definition, File Meta Information is encoded uisng the Explicit VR Little Endian
// Transfer Syntax.  So if the computer is Big Endian, this function swaps bytes at appropriate places.
// v=1 puts it in verbose mode
// Returns 0 on failure, 1 on success.
int readFileMetaInfo(const char *filename, DICOM_file_meta_info *file_meta_info, char v)
{
   FILE *fp;
   char dicomPrefix[4];
   char VR[2];
   uint2 GN=0, EN=0;
   uint2 VL2=0;
   uint4 VL4=0;
   uint4 group_length=0;
   size_t VL=0;

   file_meta_info->media_storage_SOP_class[0]='\0';
   file_meta_info->transfer_syntax[0]='\0';
   file_meta_info->dataset_offset=0;

   // open the filename for reading
   fp = fopen(filename,"r");

   if(fp==NULL) 
   {
      printf("\nreadFileMetaInfo() error: cannot open %s\n\n",filename);
      return(0);
   }

   // Skip over the 128 byte File Preamble 
   if( fseek(fp,128,SEEK_SET) == -1) 
   {
      printf("\nreadFileMetaInfo() error: 128-byte File Preamble fseek() failed for %s\n\n",filename);
      fclose(fp);
      return(0);
   }

   // The four byte DICOM Prefix must contain the character string "DICM"
   if( fread(dicomPrefix,1,4,fp) != 4) 
   {
      printf("\nreadFileMetaInfo() error: DICOM Prefix fread() failed for %s\n\n",filename);
      fclose(fp);
      return(0);
   }

   if(dicomPrefix[0]!='D' || dicomPrefix[1]!='I' || dicomPrefix[2]!='C' || dicomPrefix[3]!='M')
   {
      printf("\nreadFileMetaInfo() error: %s does not contain the required \"DICM\" prefix.\n\n", filename);
      fclose(fp);
      return(0);
   }

   if(v)
   {
      printf("Prefix = %c%c%c%c\n", dicomPrefix[0],dicomPrefix[1],dicomPrefix[2],dicomPrefix[3]);
   }

   for( ; ;)
   {
      // read the GN
      if(fread(&GN, 2, 1, fp)!=1)
      {
         printf("readFileMetaInfo() error: GN fread() failed.\n");
         fclose(fp);
         return(0);
      }

      if( bigEndian() ) swapByteOrder( (char *)&GN, 2);

      // we have reached the next group if GN is not 0x0002
      if( GN != 0x0002)
      {
         break;
      }

      // read the EN
      if(fread(&EN, 2, 1, fp)!=1)
      {
         printf("readFileMetaInfo() error: EN fread() failed.\n");
         fclose(fp);
         return(0);
      }

      if( bigEndian() ) swapByteOrder( (char *)&EN, 2);

      // read the VR 
      if(fread(VR, 1, 2, fp)!=2)
      {
         printf("readFileMetaInfo() error: VR fread() failed.\n");
         fclose(fp);
         return(0);
      }

      // see P. 36, PS 3.5-2008
      if( (VR[0]=='O' && VR[1]=='B') 
      ||  (VR[0]=='O' && VR[1]=='W') 
      ||  (VR[0]=='O' && VR[1]=='F') 
      ||  (VR[0]=='S' && VR[1]=='Q') 
      ||  (VR[0]=='U' && VR[1]=='T') 
      ||  (VR[0]=='U' && VR[1]=='N') 
      )
      {
         if( fseek(fp,2,SEEK_CUR) == -1) // jump over the reserved 2 bytes
         {
            printf("readFileMetaInfo() error: fseek() failed\n");
            fclose(fp);
            return(0);
         }

         if( fread(&VL4,4,1,fp) != 1)
         {
            printf("readFileMetaInfo() error: VL4 fread() failed.\n");
            fclose(fp);
            return(0);
         }

         if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

         VL = VL4;
      }
      else
      {
         if( fread(&VL2,2,1,fp) != 1)
         {
            printf("readFileMetaInfo() error: VL2 fread() failed.\n");
            fclose(fp);
            return(0);
         }

         if( bigEndian() ) swapByteOrder( (char *)&VL2, 2);

         VL = VL2;
      }

      if(GN==0x0002 && EN==0x0000)
      {
         if( fread(&group_length,4,1,fp) != 1)
         {
            printf("readFileMetaInfo() error: group_length fread() failed.\n");
            fclose(fp);
            return(0);
         }

         if( bigEndian() ) swapByteOrder( (char *)&group_length, 4);
        
         if(v)
         {
            printf("\nTag (%04X,%04X)\n", GN,EN);
            printf("VR = %c%c\n", VR[0], VR[1]);
            printf("Value Length = %lu\n", VL);
            printf("Group Length = %d\n", group_length);
         }
      }
      else if(GN==0x0002 && EN==0x0002)
      {
         if( fread(file_meta_info->media_storage_SOP_class, 1, VL, fp) != VL )
         {
            printf("readFileMetaInfo() error: media_storage_SOP_class fread() failed.\n");
            fclose(fp);
            return(0);
         }
         file_meta_info->media_storage_SOP_class[VL]='\0';

         if(v)
         {
            printf("\nTag (%04X,%04X)\n", GN,EN);
            printf("VR = %c%c\n", VR[0], VR[1]);
            printf("Value Length = %lu\n", VL);
            printf("Media Storage SOP Class UID = %s\n", file_meta_info->media_storage_SOP_class);
         }
      }
      else if(GN==0x0002 && EN==0x0010)
      {
         if( fread(file_meta_info->transfer_syntax, 1, VL, fp) != VL )
         {
            printf("readFileMetaInfo() error: transfer_syntax fread() failed.\n");
            fclose(fp);
            return(0);
         }
         file_meta_info->transfer_syntax[VL]='\0';

         if(v)
         {
            printf("\nTag (%04X,%04X)\n", GN,EN);
            printf("VR = %c%c\n", VR[0], VR[1]);
            printf("Value Length = %lu\n", VL);
            printf("Transfer Syntax UID = %s\n", file_meta_info->transfer_syntax);
         }
      }
      else
      {
         if( fseek(fp, VL, SEEK_CUR) == -1)
         {
            printf("readFileMetaInfo() error: fseek() failed\n");
            fclose(fp);
            return(0);
         }
      }
   }
   fclose(fp);

   if(group_length == 0)
   {
      printf("readFileMetaInfo() error: could not determine (0002) group length.\n");
      return(0);
   }

   if(file_meta_info->media_storage_SOP_class[0]=='\0')
   {
      printf("readFileMetaInfo() error: could not determine media storage SOP class.\n");
      return(0);
   }

   if(file_meta_info->transfer_syntax[0]=='\0')
   {
      printf("readFileMetaInfo() error: could not determine transfer syntax.\n");
      return(0);
   }

   // 144 = 128 (Preamble) + 4 (DICM) + 2 (Group #) + 2 (Element #) + 2 (VR) + 2 (VL) + 4 (Group Length).
   file_meta_info->dataset_offset = group_length + 144;

   return(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////

int readDataSet(const char *filename, DICOM_hdr *hdr, char v)
{
   FILE *fp;
   DICOM_file_meta_info file_meta_info;
   char VR[2];
   uint2 GN=0, EN=0;
   uint2 VL2=0;
   uint4 VL4=0;
   size_t VL=0;

   char DS[17]; // any data with VR of DS
   char IS[13]; // any data with VR of IS
   uint2 US; // any data with VR of US

   if( readFileMetaInfo(filename, &file_meta_info, 0) == 0)
   {
      printf("readDataSet() error: readFileMetaInfo() failed.\n");
      return(0);
   }

   fp = fopen(filename,"r");

   if(fp==NULL) 
   {
      printf("\nreadDataSet() error: cannot open %s\n\n",filename);
      return(0);
   }

   if( fseek(fp, file_meta_info.dataset_offset, SEEK_CUR) == -1)
   {
      printf("readDataSet() error: fseek() failed\n");
      fclose(fp);
      return(0);
   }

   // initializations
   hdr->patientID[0]='\0';
   hdr->slice_thickness=0.0;
   hdr->TR=0;
   hdr->TE=0;
   hdr->TI=0;
   hdr->NEX=0;
   hdr->dz=0.0;
   hdr->flip_angle=0.0;
   hdr->seriesID[0]='\0';
   hdr->series_number=0;
   hdr->TLHC[0] = hdr->TLHC[1] = hdr->TLHC[2] = 0.0;
   hdr->rowvec[0] = hdr->rowvec[1] = hdr->rowvec[2] = 0.0;
   hdr->colvec[0] = hdr->colvec[1] = hdr->colvec[2] = 0.0;
   hdr->frame_of_referenceID[0]='\0';
   hdr->nz=0;
   hdr->slice_location=0.0;
   hdr->ny=0;
   hdr->nx=0;
   hdr->dy=0.0;
   hdr->dx=0.0;

int count=0;
   for( ; ; )
   {
count++;
      printf("\ncount=%d\n",count);

      // read the GN
      if(fread(&GN, 2, 1, fp)!=1)
      {
         printf("readDataSet() error: GN fread() failed.\n");
         fclose(fp);
         return(0);
      }
      // read the EN
      if(fread(&EN, 2, 1, fp)!=1)
      {
         printf("readDataSet() error: EN fread() failed.\n");
         fclose(fp);
         return(0);
      }

      printf("Tag (%04X,%04X)\n", GN,EN);

      if( GN==0xFFFE && (EN==0xE000 || EN==0xE00D || EN==0xE0DD) ) 
      {
         // jump over the 4-byte Item Length
         if ( fseek(fp, 4, SEEK_CUR) == -1 ) 
         {
            printf("readDataSet() error: fseek() failed\n");
            fclose(fp);
            return(0);
         }
         continue;
      }

      // read the VR 
      if(fread(VR, 1, 2, fp)!=2)
      {
         printf("readDataSet() error: VR fread() failed.\n");
         fclose(fp);
         return(0);
      }

      printf("VR = %c%c\n", VR[0], VR[1]);

      // see P. 36, PS 3.5-2008
      if( (VR[0]=='O' && VR[1]=='B') 
      ||  (VR[0]=='O' && VR[1]=='W') 
      ||  (VR[0]=='O' && VR[1]=='F') 
      ||  (VR[0]=='S' && VR[1]=='Q') 
      ||  (VR[0]=='U' && VR[1]=='T') 
      ||  (VR[0]=='U' && VR[1]=='N') 
      )
      {
         if( fseek(fp,2,SEEK_CUR) == -1) // jump over the reserved 2 bytes
         {
            printf("readDataSet() error: fseek() failed\n");
            fclose(fp);
            return(0);
         }

         if( fread(&VL4,4,1,fp) != 1)
         {
            printf("readDataSet() error: VL4 fread() failed.\n");
            fclose(fp);
            return(0);
         }

         VL = VL4;
      }
      else
      {
         if( fread(&VL2,2,1,fp) != 1)
         {
            printf("readDataSet() error: VL2 fread() failed.\n");
            fclose(fp);
            return(0);
         }

         VL = VL2;
      }

      if(VR[0]=='S' && VR[1]=='Q')
      {
         continue;
      }

      printf("Value Length = %lu\n", VL);

      if( GN==0x0008 && EN==0x0008) // Patient ID, VL=LO
      {
         char *dum;

         dum = (char *)calloc(VL+1, 1);

         if( fread(dum, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: patient ID fread() failed.\n");
            fclose(fp);
            return(0);
         }
         dum[VL]='\0';
         if(v) printf("Image Type = %s\n", dum);
         free(dum);

         continue;
      }

      if( GN==0x0010 && EN==0x0020) // Patient ID, VL=LO
      {
         if( fread(hdr->patientID, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: patient ID fread() failed.\n");
            fclose(fp);
            return(0);
         }
         hdr->patientID[VL]='\0';
         if(v) printf("Patient ID = %s\n", hdr->patientID);
         continue;
      }

      if( GN==0x0018 && EN==0x0050) // Slice Thickness, VL=DS
      {
         if( fread(DS, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: Slice Thickness (0018,0050) fread() failed.\n");
            fclose(fp);
            return(0);
         }
         DS[VL]='\0';
         hdr->slice_thickness = atof(DS);
         if(v) printf("Slice Thickness = %f\n", hdr->slice_thickness);
         continue;
      }

      if( GN==0x0018 && EN==0x0080) // Repetition Time, VL=DS
      {
         if( fread(DS, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: Repetition Time (0018,0080) fread() failed.\n");
            fclose(fp);
            return(0);
         }
         DS[VL]='\0';
         hdr->TR = atoi(DS);
         if(v) printf("TR = %d\n", hdr->TR);
         continue;
      }

      if( GN==0x0018 && EN==0x0081) // Echo Time, VL=DS
      {
         if( fread(DS, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: Echo Time (0018,0081) fread() failed.\n");
            fclose(fp);
            return(0);
         }
         DS[VL]='\0';
         hdr->TE = atoi(DS);
         if(v) printf("TE = %d\n", hdr->TE);
         continue;
      }

      if( GN==0x0018 && EN==0x0082) // Inversion Time, VL=DS
      {
         if( fread(DS, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: Inversion Time (0018,0082) fread() failed.\n");
            fclose(fp);
            return(0);
         }
         DS[VL]='\0';
         hdr->TI = atoi(DS);
         if(v) printf("TI = %d\n", hdr->TI);
         continue;
      }

      if( GN==0x0018 && EN==0x0083) // Number of Averages , VL=DS
      {
         if( fread(DS, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: Number of Averages (0018,0083) fread() failed.\n");
            fclose(fp);
            return(0);
         }
         DS[VL]='\0';
         hdr->NEX = atoi(DS);
         if(v) printf("Number of Averages = %d\n", hdr->NEX);
         continue;
      }

      if( GN==0x0018 && EN==0x0088) // Spacing Between Slices, VL=DS
      {
         if( fread(DS, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: Spacing Between Slices (0018,0088) fread() failed.\n");
            fclose(fp);
            return(0);
         }
         DS[VL]='\0';
         hdr->dz = atof(DS);
         if(v) printf("Spacing Between Slices = %f\n", hdr->dz);
         continue;
      }

      if( GN==0x0018 && EN==0x1314) // Flip Angle, VL=DS
      {
         if( fread(DS, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: Flip Angle (0018,1314) fread() failed.\n");
            fclose(fp);
            return(0);
         }
         DS[VL]='\0';
         hdr->flip_angle = atof(DS);
         if(v) printf("Flip Angle = %f\n", hdr->flip_angle);
         continue;
      }

      if( GN==0x0020 && EN==0x000E) // Series ID, VL=UI
      {
         if( fread(hdr->seriesID, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: Series Instance UID (0020,000E) fread() failed.\n");
            fclose(fp);
            return(0);
         }
         hdr->seriesID[VL]='\0';
         if(v) printf("Series Instance UID = %s\n", hdr->seriesID);
         continue;
      }

      if( GN==0x0020 && EN==0x0011) // Series Number, VL=IS
      {
         if( fread(IS, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: Series Number (0020,0011) fread() failed.\n");
            fclose(fp);
            return(0);
         }
         IS[VL]='\0';
         hdr->series_number = atoi(IS);
         if(v) printf("Series Number = %d\n", hdr->series_number);
         continue;
      }

      // The x, y, and z coordinates of the upper left hand corner (center of the first voxel
      // transmitted) of the image, in mm.  See PS 3.3-2008 C7.6.2.1.1. for further explanation.
      // (x,y,z)=(L,P,S) See PS 3.3-2008 Page 318
      if( GN==0x0020 && EN==0x0032) // Image Position, VL=DS
      {
         if( fread(DS, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: Image Position (0020,0032) fread() failed.\n");
            fclose(fp);
            return(0);
         }
         DS[VL]='\0';
         sscanf(DS,"%f\\%f\\%f",hdr->TLHC,hdr->TLHC+1,hdr->TLHC+2);
         if(v) printf("Image Position = (%f, %f, %f)\n", hdr->TLHC[0], hdr->TLHC[1], hdr->TLHC[2]);
         continue;
      }

      // (x,y,z)=(L,P,S) See PS 3.3-2008 Page 318
      if( GN==0x0020 && EN==0x0037) // Image Orientation , VL=DS
      {
         if( fread(DS, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: Image Orientation (0020,0037) fread() failed.\n");
            fclose(fp);
            return(0);
         }
         DS[VL]='\0';

         sscanf(DS,"%f\\%f\\%f\\%f\\%f\\%f",hdr->rowvec,hdr->rowvec+1,hdr->rowvec+2,
         hdr->colvec,hdr->colvec+1,hdr->colvec+2);

         if(v) printf("Image Row Orientaton = (%f, %f, %f)\n", hdr->rowvec[0],hdr->rowvec[1],hdr->rowvec[2]);
         if(v) printf("Image Column Orientaton = (%f, %f, %f)\n", hdr->colvec[0],hdr->colvec[1],hdr->colvec[2]);
         continue;
      }

      if( GN==0x0020 && EN==0x0052) // Frame of Reference ID, VL=UI
      {
         if( fread(hdr->frame_of_referenceID, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: Frame of Reference ID fread() failed.\n");
            fclose(fp);
            return(0);
         }
         hdr->frame_of_referenceID[VL]='\0';
         if(v) printf("Frame of Reference ID = %s\n", hdr->frame_of_referenceID);
         continue;
      }

      if( GN==0x0020 && EN==0x1002) // Images in Acquistion, VL=IS
      {
         if( fread(IS, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: Images in Acquisition (0020,1002) fread() failed.\n");
            fclose(fp);
            return(0);
         }
         IS[VL]='\0';
         hdr->nz = atoi(IS);
         if(v) printf("Number of Images in Acquisition (nz) = %d\n", hdr->nz);
         continue;
      }

      if( GN==0x0020 && EN==0x1041) // Slice Location
      {
         if( fread(DS, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: Slice Location (0020,1041) fread() failed.\n");
            fclose(fp);
            return(0);
         }
         DS[VL]='\0';
         hdr->slice_location = atof(DS);
         if(v) printf("Slice Location = %f\n", hdr->slice_location);
         continue;
      }

      if( GN==0x0028 && EN==0x0010) // Rows 
      {
         if( fread(&US, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: Rows (0028,0010) fread() failed.\n");
            fclose(fp);
            return(0);
         }
         hdr->ny = US;
         if(v) printf("Number of rows (ny) = %d\n", hdr->ny);
         continue;
      }

      if( GN==0x0028 && EN==0x0011) // Columns
      {
         if( fread(&US, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: Columns (0028,0011) fread() failed.\n");
            fclose(fp);
            return(0);
         }
         hdr->nx = US;
         if(v) printf("Number of columns (nx) = %d\n", hdr->nx);
         continue;
      }

      if( GN==0x0028 && EN==0x0030) // Pixel Spacing, VL=DS
      {
         if( fread(DS, VL, 1, fp) != 1)
         {
            printf("readDataSet() error: Pixel Spacing (0028,0030) fread() failed.\n");
            fclose(fp);
            return(0);
         }
         DS[VL]='\0';
         sscanf(DS,"%f\\%f",&(hdr->dy),&(hdr->dx));
         if(v) printf("Row spacing (dy) = %f\n", hdr->dy);
         if(v) printf("Column spacing (dx) = %f\n", hdr->dx);
         continue;
      }

      if(GN==0x7FE0 && EN==0x0010)
      {
         break;
      }
/*****************************
// somethings to consider later
GN==0x0018 && EN==0x0093) 
GN==0x0018 && EN==0x0094) 
Tag (0018,1310)  Acquisition Matrix
Tag (0018,1312)  Phase encoding axis
Tag (0020,0012)  Acquisition number
*****************************/

      if( fseek(fp,VL,SEEK_CUR) == -1) // jump over the reserved 2 bytes
      {
         printf("readDataSet() error: fseek() failed\n");
         fclose(fp);
         return(0);
      }
   }

   fclose(fp);

   return(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////

int readNEX(const char *file, int *nex)
{
	int errorFlag;
	unsigned long VL;
    int transferSyntax = IMPLICIT_LITTLE_ENDIAN;  // DICOM default
    long byteOffset = 0;
	long valueOffset = 0;
	char DS[64];

	errorFlag=readMetaFileInfo(file, &byteOffset, &transferSyntax);
	if(errorFlag) return(0);

	errorFlag=readTag(file, 0x0018, 0x0050, byteOffset, transferSyntax, &VL, DS, &valueOffset);
	if(errorFlag) return(0);
	DS[VL]='\0';
	*nex= atoi(DS);

	return(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////
int readSliceThickness(const char *file, float *sliceThickness)
{
	int errorFlag;
	unsigned long VL;
    int transferSyntax = IMPLICIT_LITTLE_ENDIAN;  // DICOM default
    long byteOffset = 0;
	long valueOffset = 0;
	char DS[64];

	errorFlag=readMetaFileInfo(file, &byteOffset, &transferSyntax);
	if(errorFlag) return(0);

	errorFlag=readTag(file, 0x0018, 0x0050, byteOffset, transferSyntax, &VL, DS, &valueOffset);
	if(errorFlag) return(0);
	DS[VL]='\0';
	*sliceThickness = atof(DS);

	return(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////

extern int readSeriesNumber(const char *file, int *seriesNumber, int np)
{
	int filesize, hdrsize;
	int seriesNumberflag=0;

	unsigned short GN; 
	unsigned short EN;
	unsigned short VL2;	// Value Length see: PS 3.5-1998 Page 27
	unsigned long VL4; 	// Value Length see: PS 3.5-1998 Page 27

	FILE *fp;

	int errorFlag;
    int transferSyntax = IMPLICIT_LITTLE_ENDIAN;  // DICOM default
    long byteOffset = 0;

	// The VR of Series Number is IS () (Page 18, PS 3.6-1998).
	char IS[13]; // IS has 12 bytes maximum (Page 16, PS 3.5-1998).
	char VR[3]; // DICOM Value Representation   see: PS 3.5-1998 Page 15 and 27

	filesize = getFileSize(file);
	if(filesize==0) return(0);

	hdrsize = filesize-2*np;

	errorFlag=readMetaFileInfo(file, &byteOffset, &transferSyntax);
	if(errorFlag) return(0);

   fp=fopen(file,"r");
   if(fp==NULL) file_open_error(file);

	for(int b=byteOffset; b<hdrsize; b++)
	{
		fseek(fp,b,SEEK_SET);

		if( fread(&GN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&GN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);

		if( fread(&EN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&EN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);

		if(GN==0x0020 && EN==0x0011) 
		{ 
			if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
			{
				if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
				if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

				if(VL4>13) {errorFlag=1; break;}

		    	if( fread(IS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
				IS[VL4]='\0';
			}
			else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
			{
				if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
				VR[2]='\0';

				if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
				if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

				if(VL2>13) {errorFlag=1; break;}

		    	if( fread(IS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
				IS[VL2]='\0';
			}

			*seriesNumber=atoi(IS);
			seriesNumberflag=1;
			break;
		}
	};	

	fclose(fp);

	if(errorFlag || !seriesNumberflag) return(0);

	return(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////

int readImageNumber(const char *file, int *imageNumber, int np)
{
	int filesize, hdrsize;
	int imageNumberflag=0;

	unsigned short GN; 
	unsigned short EN;
	unsigned short VL2;	// Value Length see: PS 3.5-1998 Page 27
	unsigned long VL4; 	// Value Length see: PS 3.5-1998 Page 27

	FILE *fp;

	int errorFlag;
    int transferSyntax = IMPLICIT_LITTLE_ENDIAN;  // DICOM default
    long byteOffset = 0;

	// The VR of Series Number is IS () (Page 18, PS 3.6-1998).
	char IS[13]; // IS has 12 bytes maximum (Page 16, PS 3.5-1998).
	char VR[3]; // DICOM Value Representation   see: PS 3.5-1998 Page 15 and 27

	filesize = getFileSize(file);
	if(filesize==0) return(0);

	hdrsize = filesize-2*np;

	errorFlag=readMetaFileInfo(file, &byteOffset, &transferSyntax);
	if(errorFlag) return(0);

   fp=fopen(file,"r");
   if(fp==NULL) file_open_error(file);

	for(int b=byteOffset; b<hdrsize; b++)
	{
		fseek(fp,b,SEEK_SET);

		if( fread(&GN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&GN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);

		if( fread(&EN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&EN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);

		if(GN==0x0020 && EN==0x0013) 
		{ 
			if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
			{
				if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
				if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

				if(VL4>13) {errorFlag=1; break;}

		    	if( fread(IS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
				IS[VL4]='\0';
			}
			else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
			{
				if( fread(VR, 1, 2, fp) != 2) { errorFlag=1; break;}
				VR[2]='\0';

				if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
				if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

				if(VL2>13) {errorFlag=1; break;}

		    	if( fread(IS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
				IS[VL2]='\0';
			}

			*imageNumber=atoi(IS);
			imageNumberflag=1;
			break;
		}
	};	

	fclose(fp);

	if(errorFlag || !imageNumberflag) return(0);

	return(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
void readDicomInfo(const char *file, int np, dicominfo *info)
{ 
   int filesize, hdrsize;
   char patientName_found=0;
   char patientID_found=0;
   char DOB_found=0;
   char studyDate_found=0;
   char TE_found=0;
   char TR_found=0;
   char ETL_found=0;
   char flipAngle_found=0;

   unsigned short GN; 
   unsigned short EN;
   unsigned short VL2;	// Value Length see: PS 3.5-1998 Page 27
   unsigned long VL4; 	// Value Length see: PS 3.5-1998 Page 27

   FILE *fp;

   int errorFlag;
   int transferSyntax = IMPLICIT_LITTLE_ENDIAN;  // DICOM default

   long byteOffset = 0;

   // The VR of "Series Number" and "Image Number" is IS (Page 18, PS 3.6-1998).
   char VR[3]; // DICOM Value Representation   see: PS 3.5-1998 Page 15 and 27

   info->patientName[0]='\0';
   info->DOB[0]='\0';
   info->patientID[0]='\0';
   info->studyDate[0]='\0';
   info->TE[0]='\0';
   info->TI[0]='\0';
   info->TR[0]='\0';
   info->ETL[0]='\0';
   info->NEX[0]='\0';
   info->flipAngle[0]='\0';
   info->bandwidth[0]='\0';
   info->freq[0]='\0';
   for(int i=0; i<4; i++)
   {
      info->acquisitionMatrix[i]=0;
   }
   info->phaseFOV[0]='\0';

   filesize = getFileSize(file);
   if(filesize==0) 
   {
      return;
   }

   hdrsize = filesize-2*np;

   errorFlag=readMetaFileInfo(file, &byteOffset, &transferSyntax);
   if(errorFlag) 
   {
      return;
   }

   fp=fopen(file,"r");
   if(fp==NULL) file_open_error(file);

   for(int b=byteOffset; b<hdrsize; b++)
   {
      fseek(fp,b,SEEK_SET);

      if( fread(&GN, 2, 1, fp) != 1)  { errorFlag=1; break; }
      if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&GN, 2);
      if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);

      if( fread(&EN, 2, 1, fp) != 1)  { errorFlag=1; break; }
      if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&EN, 2);
      if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);

      if(GN==0x0008 && EN==0x0020) 
      { 
         studyDate_found=1;

         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; continue; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);
      
            if(VL4>10) {errorFlag=1; continue;}
      
            if( fread(info->studyDate, 1, VL4, fp) != VL4) {errorFlag=1; continue;}
            info->studyDate[VL4]='\0';

         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; continue;}
            VR[2]='\0';
      
            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; continue;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);
      
            if(VL2>10) {errorFlag=1; continue;}
      
            if( fread(info->studyDate, 1, VL2, fp) != VL2) {errorFlag=1; continue;}
            info->studyDate[VL2]='\0';
         }
      }

      if(GN==0x0010 && EN==0x0010) 
      { 
         patientName_found=1;

         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; continue; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);
      
            if(VL4>320) {errorFlag=1; continue;}
      
            if( fread(info->patientName, 1, VL4, fp) != VL4) {errorFlag=1; continue;}
            info->patientName[VL4]='\0';

         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; continue;}
            VR[2]='\0';
      
            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; continue;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);
      
            if(VL2>320) {errorFlag=1; continue;}
      
            if( fread(info->patientName, 1, VL2, fp) != VL2) {errorFlag=1; continue;}
            info->patientName[VL2]='\0';
         }

      }

      if(GN==0x0010 && EN==0x0020) 
      { 
         patientID_found=1;

         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; continue; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);
      
            if(VL4>64) {errorFlag=1; continue;}
      
            if( fread(info->patientID, 1, VL4, fp) != VL4) {errorFlag=1; continue;}
            info->patientID[VL4]='\0';

         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; continue;}
            VR[2]='\0';
      
            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; continue;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);
      
            if(VL2>64) {errorFlag=1; continue;}
      
            if( fread(info->patientID, 1, VL2, fp) != VL2) {errorFlag=1; continue;}
            info->patientID[VL2]='\0';
         }
      }

      if(GN==0x0010 && EN==0x0030) 
      { 
         DOB_found=1;

         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; continue; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);
      
            if(VL4>10) {errorFlag=1; continue;}
      
            if( fread(info->DOB, 1, VL4, fp) != VL4) {errorFlag=1; continue;}
            info->DOB[VL4]='\0';

         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; continue;}
            VR[2]='\0';
      
            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; continue;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);
      
            if(VL2>10) {errorFlag=1; continue;}
      
            if( fread(info->DOB, 1, VL2, fp) != VL2) {errorFlag=1; continue;}
            info->DOB[VL2]='\0';
         }
      }

      if(GN==0x0018 && EN==0x0080) 
      { 
         TR_found=1;

         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; continue; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

            if(VL4>16) {errorFlag=1; continue;}

            if( fread(info->TR, 1, VL4, fp) != VL4) {errorFlag=1; continue;}
            info->TR[VL4]='\0';
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; continue;}
            VR[2]='\0';

            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; continue;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

            if(VL2>16) {errorFlag=1; continue;}

            if( fread(info->TR, 1, VL2, fp) != VL2) {errorFlag=1; continue;}
            info->TR[VL2]='\0';
         }
      }

      if(GN==0x0018 && EN==0x0081) 
      { 
         TE_found=1;

         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; continue; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

            if(VL4>16) {errorFlag=1; continue;}

            if( fread(info->TE, 1, VL4, fp) != VL4) {errorFlag=1; continue;}
            info->TE[VL4]='\0';
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; continue;}
            VR[2]='\0';

            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; continue;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

            if(VL2>16) {errorFlag=1; continue;}

            if( fread(info->TE, 1, VL2, fp) != VL2) {errorFlag=1; continue;}
            info->TE[VL2]='\0';
         }
      }

      if(GN==0x0018 && EN==0x0082) 
      { 
         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; continue; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

            if(VL4>16) {errorFlag=1; continue;}

            if( fread(info->TI, 1, VL4, fp) != VL4) {errorFlag=1; continue;}
            info->TI[VL4]='\0';
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; continue;}
            VR[2]='\0';

            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; continue;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

            if(VL2>16) {errorFlag=1; continue;}

            if( fread(info->TI, 1, VL2, fp) != VL2) {errorFlag=1; continue;}
            info->TI[VL2]='\0';
         }
      }

      if(GN==0x0018 && EN==0x0083) 
      { 
         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; continue; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

            if(VL4>16) {errorFlag=1; continue;}

            if( fread(info->NEX, 1, VL4, fp) != VL4) {errorFlag=1; continue;}
            info->NEX[VL4]='\0';
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; continue;}
            VR[2]='\0';

            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; continue;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

            if(VL2>16) {errorFlag=1; continue;}

            if( fread(info->NEX, 1, VL2, fp) != VL2) {errorFlag=1; continue;}
            info->NEX[VL2]='\0';
         }
      }

      if(GN==0x0018 && EN==0x0084) 
      { 
         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; continue; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

            if(VL4>16) {errorFlag=1; continue;}

            if( fread(info->freq, 1, VL4, fp) != VL4) {errorFlag=1; continue;}
            info->freq[VL4]='\0';
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; continue;}
            VR[2]='\0';

            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; continue;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

            if(VL2>16) {errorFlag=1; continue;}

            if( fread(info->freq, 1, VL2, fp) != VL2) {errorFlag=1; continue;}
            info->freq[VL2]='\0';
         }
      }

      if(GN==0x0018 && EN==0x0091) 
      { 
         ETL_found=1;

         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; continue; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

            if(VL4>12) {errorFlag=1; continue;}

            if( fread(info->ETL, 1, VL4, fp) != VL4) {errorFlag=1; continue;}
            info->ETL[VL4]='\0';
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; continue;}
            VR[2]='\0';

            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; continue;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

            if(VL2>12) {errorFlag=1; continue;}

            if( fread(info->ETL, 1, VL2, fp) != VL2) {errorFlag=1; continue;}
            info->ETL[VL2]='\0';
         }
      }

      if(GN==0x0018 && EN==0x0094) 
      { 
         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; continue; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

            if(VL4>16) {errorFlag=1; continue;}

            if( fread(info->phaseFOV, 1, VL4, fp) != VL4) {errorFlag=1; continue;}
            info->phaseFOV[VL4]='\0';
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; continue;}
            VR[2]='\0';

            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; continue;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

            if(VL2>16) {errorFlag=1; continue;}

            if( fread(info->phaseFOV, 1, VL2, fp) != VL2) {errorFlag=1; continue;}
            info->phaseFOV[VL2]='\0';
         }
      }

      if(GN==0x0018 && EN==0x0095) 
      { 
         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; continue; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

            if(VL4>16) {errorFlag=1; continue;}

            if( fread(info->bandwidth, 1, VL4, fp) != VL4) {errorFlag=1; continue;}
            info->bandwidth[VL4]='\0';
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; continue;}
            VR[2]='\0';

            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; continue;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

            if(VL2>16) {errorFlag=1; continue;}

            if( fread(info->bandwidth, 1, VL2, fp) != VL2) {errorFlag=1; continue;}
            info->bandwidth[VL2]='\0';
         }
      }

      if(GN==0x0018 && EN==0x1310) 
      { 
         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; continue; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

            if(VL4 != 8) {errorFlag=1; continue;}

            if( fread(info->acquisitionMatrix, 1, VL4, fp) != VL4) {errorFlag=1; continue;}
            if( bigEndian() ) swapN( (char *)info->acquisitionMatrix, 8);
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; continue;}
            VR[2]='\0';

            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; continue;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

            if(VL2 != 8) {errorFlag=1; continue;}

            if( fread(info->acquisitionMatrix, 1, VL2, fp) != VL2) {errorFlag=1; continue;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapN( (char *)info->acquisitionMatrix, 8);
            if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapN( (char *)info->acquisitionMatrix, 8);
         }
      }

      if(GN==0x0018 && EN==0x1314) 
      { 
         flipAngle_found=1;

         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; continue; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

            if(VL4>16) {errorFlag=1; continue;}

            if( fread(info->flipAngle, 1, VL4, fp) != VL4) {errorFlag=1; continue;}
            info->flipAngle[VL4]='\0';
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; continue;}
            VR[2]='\0';

            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; continue;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

            if(VL2>16) {errorFlag=1; continue;}

            if( fread(info->flipAngle, 1, VL2, fp) != VL2) {errorFlag=1; continue;}
            info->flipAngle[VL2]='\0';
         }
      }

      // no need to pass this element (0020,000D) since all need information should be before this
      if(GN==0x0020 && EN==0x000D) 
      {
         break;
      }

      // found everthing I was looking for
      if(patientName_found && DOB_found && patientID_found && studyDate_found && TE_found && TR_found && ETL_found && flipAngle_found) 
      {
         break;
      }
   };	

   fclose(fp);
} 

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

int readImageParams(const char *file, float *TLHC, float *rowvec, float *colvec, 
float *dx, float *dy, float *dz, int *TE, char *patientID, int *imageNumber, int *seriesNumber, int np)
{
   int filesize, hdrsize;
   int imageNumberflag=0;
   int seriesNumberflag=0;
   int patientIDflag=0;
   int TEflag=0;
   int dx_dy_flag=0;
   int rowcolvecflag=0;
   int TLHCflag=0;
   int dz_flag = 0;

   unsigned short GN; 
   unsigned short EN;
   unsigned short VL2;	// Value Length see: PS 3.5-1998 Page 27
   unsigned long VL4; 	// Value Length see: PS 3.5-1998 Page 27

   FILE *fp;

   int errorFlag;
   int transferSyntax = IMPLICIT_LITTLE_ENDIAN;  // DICOM default

   long byteOffset = 0;

   // The VR of "Series Number" and "Image Number" is IS (Page 18, PS 3.6-1998).
   char IS[13]; // IS has 12 bytes maximum (Page 16, PS 3.5-1998).
   char VR[3]; // DICOM Value Representation   see: PS 3.5-1998 Page 15 and 27

   // The VR of TE is DS  (Page 10, PS 3.6-1998).
   char DS[256];

   char dum[64];
   int L1,L2;

   filesize = getFileSize(file);
   if(filesize==0) return(0);

   hdrsize = filesize-2*np;

   errorFlag=readMetaFileInfo(file, &byteOffset, &transferSyntax);
   if(errorFlag) return(0);

   fp=fopen(file,"r");
   if(fp==NULL) file_open_error(file);

   for(int b=byteOffset; b<hdrsize; b++)
   {
      fseek(fp,b,SEEK_SET);

      if( fread(&GN, 2, 1, fp) != 1)  { errorFlag=1; break; }
      if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&GN, 2);
      if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);

      if( fread(&EN, 2, 1, fp) != 1)  { errorFlag=1; break; }
      if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&EN, 2);
      if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);

      if(GN==0x0010 && EN==0x0020) 
      { 
         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);
      
            if(VL4>65) {errorFlag=1; break;}
      
            if( fread(patientID, 1, VL4, fp) != VL4) {errorFlag=1; break;}
            patientID[VL4]='\0';
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
            VR[2]='\0';
      
            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);
      
            if(VL2>65) {errorFlag=1; break;}
      
            if( fread(patientID, 1, VL2, fp) != VL2) {errorFlag=1; break;}
            patientID[VL2]='\0';
         }

         patientIDflag=1;
      }

      // Some 3D sequences may not have element 88 of group 18.
      if(GN==0x0018 && EN==0x0050) 
      { 
         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

            if(VL4>64) {errorFlag=1; break;}

            if( fread(DS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
            DS[VL4]='\0';
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
            VR[2]='\0';

            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

            if(VL2>64) {errorFlag=1; break;}

            if( fread(DS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
            DS[VL2]='\0';
         }

         *dz = atof(DS);

         dz_flag=1;
      }

      if(GN==0x0018 && EN==0x0081) 
      { 
         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

            if(VL4>16) {errorFlag=1; break;}

            if( fread(DS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
            DS[VL4]='\0';
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
            VR[2]='\0';

            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

            if(VL2>16) {errorFlag=1; break;}

            if( fread(DS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
            DS[VL2]='\0';
         }

         *TE=atoi(DS);
         TEflag=1;
      }

      if(GN==0x0018 && EN==0x0088) 
      {       
         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

            if(VL4>64) {errorFlag=1; break;}

            if( fread(DS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
            DS[VL4]='\0';
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
            VR[2]='\0';

            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

            if(VL2>64) {errorFlag=1; break;}

            if( fread(DS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
            DS[VL2]='\0';
         }

         *dz = atof(DS);

         dz_flag=1;
      }

      if(GN==0x0020 && EN==0x0011) 
      { 
         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

            if(VL4>13) {errorFlag=1; break;}

            if( fread(IS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
            IS[VL4]='\0';
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
            VR[2]='\0';

            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

            if(VL2>13) {errorFlag=1; break;}

            if( fread(IS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
            IS[VL2]='\0';
         }

         *seriesNumber=atoi(IS);
         seriesNumberflag=1;
      }

      if(GN==0x0020 && EN==0x0013) 
      { 
         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);
      
            if(VL4>13) {errorFlag=1; break;}
      
            if( fread(IS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
            IS[VL4]='\0';
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) { errorFlag=1; break;}
            VR[2]='\0';
      
            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);
      
            if(VL2>13) {errorFlag=1; break;}
      
            if( fread(IS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
            IS[VL2]='\0';
         }

         *imageNumber=atoi(IS);
         imageNumberflag=1;
      }

      if(GN==0x0020 && EN==0x0032) 
      { 
         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);
            
            if(VL4>256) {errorFlag=1; break;}

            if( fread(DS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
            DS[VL4] = '\0';
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
            VR[2]='\0';

            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

            if(VL2>256) {errorFlag=1; break;}
            if( fread(DS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
            DS[VL2] = '\0';
         }

         sscanf(DS,"%f\\%f\\%f",TLHC,TLHC+1,TLHC+2);

         TLHCflag=1;
      }

      if(GN==0x0020 && EN==0x0037) 
      { 
         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

            if(VL4>256) {errorFlag=1; break;}

            if( fread(DS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
            VR[2]='\0';

            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

            if(VL2>256) {errorFlag=1; break;}

            if( fread(DS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
         }

         sscanf(DS,"%f\\%f\\%f\\%f\\%f\\%f",rowvec,rowvec+1,rowvec+2,colvec,colvec+1,colvec+2);

         rowcolvecflag=1;
      }

      if(GN==0x0028 && EN==0x0030) 
      { 
         if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
         {
            if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
            if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

            if(VL4>64) {errorFlag=1; break;}

            if( fread(DS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
         }
         else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
         {
            if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
            VR[2]='\0';

            if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
            if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

            if(VL2>64) {errorFlag=1; break;}

            if( fread(DS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
         }

         L1=strcspn(DS,"\\");
         strncpy(dum,DS,L1); dum[L1]='\0';
         *dx=atof(dum);

         strcpy(dum,DS+L1+1);
         L2=strlen(DS+L1+1);
         dum[L2]='\0';
         *dy=atof(dum);

         dx_dy_flag=1;
      }

      // found everthing I was looking for
      if(imageNumberflag && seriesNumberflag && patientIDflag && TEflag && dx_dy_flag 
      && rowcolvecflag && TLHCflag && dz_flag) 
      {
         break;
      }
   };	

   fclose(fp);

   if(errorFlag || !imageNumberflag || !seriesNumberflag || !patientIDflag || !TEflag 
   || !dx_dy_flag || !rowcolvecflag || !TLHCflag ) return(0);

   return(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////

int readTE(const char *file, int *TE, int np)
{
	int filesize, hdrsize;
	int TEflag=0;

	unsigned short GN; 
	unsigned short EN;
	unsigned short VL2;	// Value Length see: PS 3.5-1998 Page 27
	unsigned long VL4; 	// Value Length see: PS 3.5-1998 Page 27

	FILE *fp;

	int errorFlag;
    int transferSyntax = IMPLICIT_LITTLE_ENDIAN;  // DICOM default
    long byteOffset = 0;

	// The VR of Series Number is IS () (Page 18, PS 3.6-1998).
	char DS[16];
	char VR[3]; // DICOM Value Representation   see: PS 3.5-1998 Page 15 and 27

	filesize = getFileSize(file);
	if(filesize==0) return(0);

	hdrsize = filesize-2*np;

	errorFlag=readMetaFileInfo(file, &byteOffset, &transferSyntax);
	if(errorFlag) return(0);

   fp=fopen(file,"r");
   if(fp==NULL) file_open_error(file);

	for(int b=byteOffset; b<hdrsize; b++)
	{
		fseek(fp,b,SEEK_SET);

		if( fread(&GN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&GN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);

		if( fread(&EN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&EN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);

		if(GN==0x0018 && EN==0x0081) 
		{ 
			if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
			{
				if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
				if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

				if(VL4>16) {errorFlag=1; break;}

		    	if( fread(DS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
				DS[VL4]='\0';
			}
			else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
			{
				if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
				VR[2]='\0';

				if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			    if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

				if(VL2>16) {errorFlag=1; break;}

		    	if( fread(DS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
				DS[VL2]='\0';
			}

			*TE=atoi(DS);
			TEflag=1;
			break;
		}
	};	

	fclose(fp);

	if(errorFlag || !TEflag) return(0);

	return(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////

int readPhaseEncodingDirection(const char *file, char *PED, int np)
{
	int filesize, hdrsize;
	int PEDflag=0;

	unsigned short GN; 
	unsigned short EN;
	unsigned short VL2;	// Value Length see: PS 3.5-1998 Page 27
	unsigned long VL4; 	// Value Length see: PS 3.5-1998 Page 27

	FILE *fp;

	int errorFlag;
    int transferSyntax = IMPLICIT_LITTLE_ENDIAN;  // DICOM default
    long byteOffset = 0;

	// The VR of Phase Encoding Direction is CS (Page 14, PS 3.6-1998).
	char CS[16];
	char VR[3]; // DICOM Value Representation   see: PS 3.5-1998 Page 15 and 27

	filesize = getFileSize(file);
	if(filesize==0) return(0);

	hdrsize = filesize-2*np;

	errorFlag=readMetaFileInfo(file, &byteOffset, &transferSyntax);
	if(errorFlag) return(0);

   fp=fopen(file,"r");
   if(fp==NULL) file_open_error(file);

	for(int b=byteOffset; b<hdrsize; b++)
	{
		fseek(fp,b,SEEK_SET);

		if( fread(&GN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&GN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);

		if( fread(&EN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&EN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);

		if(GN==0x0018 && EN==0x1312) 
		{ 
			if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
			{
				if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
				if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

				if(VL4>16) {errorFlag=1; break;}

		    	if( fread(CS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
				CS[VL4]='\0';
			}
			else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
			{
				if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
				VR[2]='\0';

				if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			    if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

				if(VL2>16) {errorFlag=1; break;}

		    	if( fread(CS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
				CS[VL2]='\0';
			}

			strcpy(PED, CS);
			PEDflag=1;
			break;
		}
	};	

	fclose(fp);

	if(errorFlag || !PEDflag) return(0);

	return(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////

int readTR(const char *file, int *TR)
{
	int errorFlag;
    int transferSyntax = IMPLICIT_LITTLE_ENDIAN;  // DICOM default
    long byteOffset = 0;
	long valueOffset = 0;
	unsigned long VL;
	char DS[16];

	errorFlag=readMetaFileInfo(file, &byteOffset, &transferSyntax);
	if(errorFlag) return(0);

	errorFlag=readTag(file, 0x0018, 0x0080, byteOffset, transferSyntax, &VL, DS, &valueOffset);
	if(errorFlag) return(0);
	DS[VL]='\0';
	*TR=atoi(DS);

	return(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////

int readImageSliceThickness(const char *file, float *dz, int np)
{
	int filesize, hdrsize;

	int dz_flag = 0;

	unsigned short GN; 
	unsigned short EN;
	unsigned short VL2;	// Value Length see: PS 3.5-1998 Page 27
	unsigned long VL4; 	// Value Length see: PS 3.5-1998 Page 27

	char DS[64];

	FILE *fp;

	int errorFlag;
    int transferSyntax = IMPLICIT_LITTLE_ENDIAN;  // DICOM default
    long byteOffset = 0;

	// The VR of Series Number is IS () (Page 18, PS 3.6-1998).
	char VR[3]; // DICOM Value Representation   see: PS 3.5-1998 Page 15 and 27

	filesize = getFileSize(file);
	if(filesize==0) return(0);

	hdrsize = filesize-2*np;

	errorFlag=readMetaFileInfo(file, &byteOffset, &transferSyntax);
	if(errorFlag) return(0);

   fp=fopen(file,"r");
   if(fp==NULL) file_open_error(file);

	// Some 3D sequences may not have element 88 of group 18.
	for(int b=byteOffset; b<hdrsize; b++)
	{
		fseek(fp,b,SEEK_SET);

		if( fread(&GN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&GN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);

		if( fread(&EN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&EN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);

		if(GN==0x0018 && EN==0x0050) 
		{ 
			if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
			{
				if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
				if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

				if(VL4>64) {errorFlag=1; break;}

		    	if( fread(DS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
				DS[VL4]='\0';
			}
			else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
			{
				if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
				VR[2]='\0';

				if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			    if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

				if(VL2>64) {errorFlag=1; break;}

		    	if( fread(DS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
				DS[VL2]='\0';
			}

			*dz = atof(DS);

			dz_flag=1;

			break;
		}
	};	

	fclose(fp);

	if(errorFlag) return(0);

   fp=fopen(file,"r");
   if(fp==NULL) file_open_error(file);

	for(int b=byteOffset; b<hdrsize; b++)
	{
		fseek(fp,b,SEEK_SET);

		if( fread(&GN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&GN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);

		if( fread(&EN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&EN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);

		if(GN==0x0018 && EN==0x0088) 
		{ 
			if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
			{
				if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
				if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

				if(VL4>64) {errorFlag=1; break;}

	    		if( fread(DS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
				DS[VL4]='\0';
			}
			else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
			{
				if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
				VR[2]='\0';

				if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			    if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

				if(VL2>64) {errorFlag=1; break;}

	    		if( fread(DS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
				DS[VL2]='\0';
			}

			*dz = atof(DS);

			dz_flag=1;

			break;
		}
	};	

	fclose(fp);

	if(errorFlag || !dz_flag) return(0);

	return(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////

int readImageVoxelSize(const char *file, float *dx, float *dy, float *dz, int np)
{
	int filesize, hdrsize;

	int dx_dy_flag=0;
	int dz_flag = 0;

	unsigned short GN; 
	unsigned short EN;
	unsigned short VL2;	// Value Length see: PS 3.5-1998 Page 27
	unsigned long VL4; 	// Value Length see: PS 3.5-1998 Page 27

	char DS[64];
	char dum[64];
	int L1,L2;

	FILE *fp;

	int errorFlag;
    int transferSyntax = IMPLICIT_LITTLE_ENDIAN;  // DICOM default
    long byteOffset = 0;

	// The VR of Series Number is IS () (Page 18, PS 3.6-1998).
	char VR[3]; // DICOM Value Representation   see: PS 3.5-1998 Page 15 and 27

	filesize = getFileSize(file);
	if(filesize==0) return(0);

	hdrsize = filesize-2*np;

	errorFlag=readMetaFileInfo(file, &byteOffset, &transferSyntax);
	if(errorFlag) return(0);

   fp=fopen(file,"r");
   if(fp==NULL) file_open_error(file);

	for(int b=byteOffset; b<hdrsize; b++)
	{
		fseek(fp,b,SEEK_SET);

		if( fread(&GN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&GN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);

		if( fread(&EN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&EN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);

		if(GN==0x0028 && EN==0x0030) 
		{ 
			if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
			{
				if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
				if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

				if(VL4>64) {errorFlag=1; break;}

		    	if( fread(DS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
			}
			else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
			{
				if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
				VR[2]='\0';

				if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			    if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

				if(VL2>64) {errorFlag=1; break;}

		    	if( fread(DS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
			}

			L1=strcspn(DS,"\\");
			strncpy(dum,DS,L1); dum[L1]='\0';
			*dx=atof(dum);

			strcpy(dum,DS+L1+1);
			L2=strlen(DS+L1+1);
			dum[L2]='\0';
			*dy=atof(dum);

			dx_dy_flag=1;
			break;
		}
	}

	fclose(fp);

	if(errorFlag || !dx_dy_flag) return(0);

   fp=fopen(file,"r");
   if(fp==NULL) file_open_error(file);

	// Some 3D sequences may not have element 88 of group 18.
	for(int b=byteOffset; b<hdrsize; b++)
	{
		fseek(fp,b,SEEK_SET);

		if( fread(&GN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&GN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);

		if( fread(&EN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&EN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);

		if(GN==0x0018 && EN==0x0050) 
		{ 
			if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
			{
				if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
				if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

				if(VL4>64) {errorFlag=1; break;}

		    	if( fread(DS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
				DS[VL4]='\0';
			}
			else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
			{
				if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
				VR[2]='\0';

				if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			    if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

				if(VL2>64) {errorFlag=1; break;}

		    	if( fread(DS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
				DS[VL2]='\0';
			}

			*dz = atof(DS);

			dz_flag=1;

			break;
		}
	};	

	fclose(fp);

	if(errorFlag) return(0);

   fp=fopen(file,"r");
   if(fp==NULL) file_open_error(file);

	for(int b=byteOffset; b<hdrsize; b++)
	{
		fseek(fp,b,SEEK_SET);

		if( fread(&GN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&GN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);

		if( fread(&EN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&EN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);

		if(GN==0x0018 && EN==0x0088) 
		{ 
			if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
			{
				if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
				if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

				if(VL4>64) {errorFlag=1; break;}

	    		if( fread(DS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
				DS[VL4]='\0';
			}
			else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
			{
				if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
				VR[2]='\0';

				if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
			    if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

				if(VL2>64) {errorFlag=1; break;}

	    		if( fread(DS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
				DS[VL2]='\0';
			}

			*dz = atof(DS);

			dz_flag=1;

			break;
		}
	};	

	fclose(fp);

	if(errorFlag || !dz_flag) return(0);

	return(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////

int readImageMatrixSize(const char *file, unsigned short *nx, unsigned short *ny)
{
   int filesize;
   int nxflag=0;
   int nyflag=0;

   unsigned short GN; 
   unsigned short EN;
   unsigned short VL2;	// Value Length see: PS 3.5-1998 Page 27
   unsigned long VL4; 	// Value Length see: PS 3.5-1998 Page 27

   FILE *fp;

   int errorFlag;
   int transferSyntax = IMPLICIT_LITTLE_ENDIAN;  // DICOM default
   long byteOffset = 0;

   // The VR of Series Number is IS () (Page 18, PS 3.6-1998).
   char VR[3]; // DICOM Value Representation   see: PS 3.5-1998 Page 15 and 27

   filesize = getFileSize(file);
   if(filesize==0) return(0);

   errorFlag=readMetaFileInfo(file, &byteOffset, &transferSyntax);
   if(errorFlag) return(0);

   fp=fopen(file,"r");
   if(fp==NULL) file_open_error(file);

   while( byteOffset<filesize-5 ) 
   {
      fseek(fp,byteOffset,SEEK_SET);
      byteOffset++;

      if( fread(&GN, 2, 1, fp) != 1)  { errorFlag=1; break; }
      if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);
      if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);

      if( fread(&EN, 2, 1, fp) != 1)  { errorFlag=1; break; }
      if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);
      if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);

		if(GN==0x0028 && EN==0x0010) 
		{ 
			if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
			{
				if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
				if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

				if(VL4>2) {errorFlag=1; break;}

		    	if( fread(ny, 1, VL4, fp) != VL4) {errorFlag=1; break;}
				if( bigEndian() ) swapByteOrder( (char *)ny, 2);
			}
			else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
			{
				if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
				VR[2]='\0';

				if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
				if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

				if(VL2>2) {errorFlag=1; break;}

		    	if( fread(ny, 1, VL2, fp) != VL2) {errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)ny, 2);
				if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)ny, 2);
			}
			nyflag=1;
		}

		if(GN==0x0028 && EN==0x0011) 
		{ 
			if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
			{
				if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
				if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

				if(VL4>2) {errorFlag=1; break;}

		    	if( fread(nx, 1, VL4, fp) != VL4) {errorFlag=1; break;}
				if( bigEndian() ) swapByteOrder( (char *)nx, 2);
			}
			else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
			{
				if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
				VR[2]='\0';

				if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
				if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

				if(VL2>2) {errorFlag=1; break;}

		    	if( fread(nx, 1, VL2, fp) != VL2) {errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)nx, 2);
				if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)nx, 2);
			}

			nxflag=1;
		}

		if(nxflag && nyflag)
        {
			break;
        }
	};	

	fclose(fp);

	if(errorFlag || !nxflag || !nyflag) return(0);

	return(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////

int readPatientID(const char *file, char *patientID, int np)
{
	int filesize, hdrsize;
	int patientIDflag=0;

	unsigned short GN; 
	unsigned short EN;
	unsigned short VL2;	// Value Length see: PS 3.5-1998 Page 27
	unsigned long VL4; 	// Value Length see: PS 3.5-1998 Page 27

	FILE *fp;

	int errorFlag;
    int transferSyntax = IMPLICIT_LITTLE_ENDIAN;  // DICOM default
    long byteOffset = 0;

	// The VR of Series Number is IS () (Page 18, PS 3.6-1998).
	char VR[3]; // DICOM Value Representation   see: PS 3.5-1998 Page 15 and 27

	filesize = getFileSize(file);
	if(filesize==0) return(0);

	hdrsize = filesize-2*np;

	errorFlag=readMetaFileInfo(file, &byteOffset, &transferSyntax);
	if(errorFlag) return(0);

   fp=fopen(file,"r");
   if(fp==NULL) file_open_error(file);

	for(int b=byteOffset; b<hdrsize; b++)
	{
		fseek(fp,b,SEEK_SET);

		if( fread(&GN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&GN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);

		if( fread(&EN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&EN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);

		if(GN==0x0010 && EN==0x0020) 
		{ 
			if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
			{
				if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
				if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

				if(VL4>65) {errorFlag=1; break;}

		    	if( fread(patientID, 1, VL4, fp) != VL4) {errorFlag=1; break;}
				patientID[VL4]='\0';
			}
			else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
			{
				if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
				VR[2]='\0';

				if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
				if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

				if(VL2>65) {errorFlag=1; break;}

		    	if( fread(patientID, 1, VL2, fp) != VL2) {errorFlag=1; break;}
				patientID[VL2]='\0';
			}
			patientIDflag=1;
			break;
		}
	};	

	fclose(fp);

	if(errorFlag || !patientIDflag ) return(0);

	return(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////

int readRowCol(const char *file, float *rowvec, float *colvec, int np)
{
	int filesize, hdrsize;
	int rowcolvecflag=0;

	unsigned short GN; 
	unsigned short EN;
	unsigned short VL2;	// Value Length see: PS 3.5-1998 Page 27
	unsigned long VL4; 	// Value Length see: PS 3.5-1998 Page 27

	FILE *fp;

	int errorFlag;
    int transferSyntax = IMPLICIT_LITTLE_ENDIAN;  // DICOM default
    long byteOffset = 0;

	// The VR of Series Number is IS () (Page 18, PS 3.6-1998).
	char VR[3]; // DICOM Value Representation   see: PS 3.5-1998 Page 15 and 27
	char DS[256];

	filesize = getFileSize(file);
	if(filesize==0) return(0);

	hdrsize = filesize-2*np;

	errorFlag=readMetaFileInfo(file, &byteOffset, &transferSyntax);
	if(errorFlag) return(0);

   fp=fopen(file,"r");
   if(fp==NULL) file_open_error(file);

	for(int b=byteOffset; b<hdrsize; b++)
	{
		fseek(fp,b,SEEK_SET);

		if( fread(&GN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&GN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);

		if( fread(&EN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&EN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);

		if(GN==0x0020 && EN==0x0037) 
		{ 
			if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
			{
				if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
				if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

				if(VL4>256) {errorFlag=1; break;}

		    	if( fread(DS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
			}
			else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
			{
				if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
				VR[2]='\0';

				if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
				if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

				if(VL2>256) {errorFlag=1; break;}

		    	if( fread(DS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
			}

			sscanf(DS,"%f\\%f\\%f\\%f\\%f\\%f",rowvec,rowvec+1,rowvec+2,colvec,colvec+1,colvec+2);

			rowcolvecflag=1;
			break;
		}
	};	

	fclose(fp);

	if(errorFlag || !rowcolvecflag ) return(0);

	return(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////

int readTLHC(const char *file, float *TLHC, int np)
{
	int filesize, hdrsize;
	int TLHCflag=0;

	unsigned short GN; 
	unsigned short EN;
	unsigned short VL2;	// Value Length see: PS 3.5-1998 Page 27
	unsigned long VL4; 	// Value Length see: PS 3.5-1998 Page 27

	FILE *fp;

	int errorFlag;
    int transferSyntax = IMPLICIT_LITTLE_ENDIAN;  // DICOM default
    long byteOffset = 0;

	// The VR of Series Number is IS () (Page 18, PS 3.6-1998).
	char VR[3]; // DICOM Value Representation   see: PS 3.5-1998 Page 15 and 27
	char DS[256];

	filesize = getFileSize(file);
	if(filesize==0) return(0);

	hdrsize = filesize-2*np;

	errorFlag=readMetaFileInfo(file, &byteOffset, &transferSyntax);
	if(errorFlag) return(0);

   fp=fopen(file,"r");
   if(fp==NULL) file_open_error(file);

	for(int b=byteOffset; b<hdrsize; b++)
	{
		fseek(fp,b,SEEK_SET);

		if( fread(&GN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&GN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);

		if( fread(&EN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&EN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);

		if(GN==0x0020 && EN==0x0032) 
		{ 
			if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
			{
				if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
				if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

				if(VL4>256) {errorFlag=1; break;}

		    	if( fread(DS, 1, VL4, fp) != VL4) {errorFlag=1; break;}
                DS[VL4] = '\0';
			}
			else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
			{
				if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
				VR[2]='\0';

				if( fread(&VL2, 2, 1, fp) != 1) { errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL2, 2);
				if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL2, 2);

				if(VL2>256) {errorFlag=1; break;}

		    	if( fread(DS, 1, VL2, fp) != VL2) {errorFlag=1; break;}
                DS[VL2] = '\0';
			}

			sscanf(DS,"%f\\%f\\%f",TLHC,TLHC+1,TLHC+2);

			TLHCflag=1;
			break;
		}
	};	

	fclose(fp);

	if(errorFlag || !TLHCflag) return(0);

	return(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////
    
int readPixelData(const char *file, char *data, int opt_j, int np)
{
	int filesize, hdrsize;
	int dataflag=0;

	unsigned short GN; 
	unsigned short EN;
	//unsigned short VL2;	// Value Length see: PS 3.5-1998 Page 27
	unsigned long VL4; 	// Value Length see: PS 3.5-1998 Page 27

	FILE *fp;

	int errorFlag;
    int transferSyntax = IMPLICIT_LITTLE_ENDIAN;  // DICOM default
    long byteOffset = 0;

	// The VR of Series Number is IS () (Page 18, PS 3.6-1998).
	char VR[3]; // DICOM Value Representation   see: PS 3.5-1998 Page 15 and 27

	filesize = getFileSize(file);
	if(filesize==0) return(0);

	hdrsize = filesize-2*np;

	if(opt_j)
	{
		byteOffset = filesize - np*2;

		fp=fopen(file,"r");
        if(fp==NULL) file_open_error(file);

		fseek(fp,byteOffset,SEEK_SET);

		if( (int)fread(data, 1, np*2, fp) != np*2) {fclose(fp); return(0);}
		fclose(fp);

		return(1);
	}

	errorFlag=readMetaFileInfo(file, &byteOffset, &transferSyntax);
	if(errorFlag) return(0);

   fp=fopen(file,"r");
   if(fp==NULL) file_open_error(file);

	for(int b=byteOffset; b<hdrsize; b++)
	{
		fseek(fp,b,SEEK_SET);

		if( fread(&GN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&GN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);

		if( fread(&EN, 2, 1, fp) != 1)  { errorFlag=1; break; }
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&EN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);

		if(GN==0x7FE0 && EN==0x0010) 
		{ 
			if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
			{
				if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break; }
				if( bigEndian() ) swapByteOrder( (char *)&VL4, 4);

                // important: in case the wrong (7FE0,0010) is found
                if( (int)VL4 != 2*np) continue;

		    	if( fread(data, 1, VL4, fp) != VL4) {errorFlag=1; break;}
				if( bigEndian() ) swapN( data, VL4);
			}
			else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
			{
				if( fread(VR, 1, 2, fp) != 2) {errorFlag=1; break;}
				VR[2]='\0';

				// skip the two reserved bytes, see Table 7.1-1, PS 3.5-1998, Page 27
				if( fseek(fp,2,SEEK_CUR) == -1) { errorFlag=1; break;}

				if( fread(&VL4, 4, 1, fp) != 1) { errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)&VL4, 4);
				if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)&VL4, 4);

                // important: in case the wrong (7FE0,0010) is found
                if( (int)VL4 != 2*np) continue;

		    	if( fread(data, 1, VL4, fp) != VL4) {errorFlag=1; break;}
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapN( data, VL4);
				if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapN( data, VL4);
			}

			dataflag=1;
			break;
		}
	};	

	fclose(fp);

	if(errorFlag || !dataflag) return(0);

	return(1);
}
///////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////
// Reads DICOM File Meta Information
// The File Meta Information consists of a 128 byte File Preamble, followed by a 4 byte
// DICOM Prefix, followed by group 0002 containing the File Meta Elements
// returns 0 if no error occurs, 1 otherwise
///////////////////////////////////////////////////////////////////////////////////////////////
int readMetaFileInfo(const char *file, long *offset, int *syntax)
{
   FILE *fp;
   int errorFlag;
   char dicomPrefix[4];
   char UI[65]; // unique identifier
   char VR[3]; // DICOM Value Representation	see: PS 3.5-1998 Page 15 and 27
   unsigned short *GN = new unsigned short; // group number
   unsigned short *EN = new unsigned short; // element number
   unsigned short *VL2 = new unsigned short; // Value Length see: PS 3.5-1998 Page 27
   unsigned long *VL4 = new unsigned long; // Value Length see: PS 3.5-1998 Page 27
   unsigned long *GL = new unsigned long; // Group Length see Tag:(0002,0000)

   // open the file for reading
   fp = fopen(file,"r");

   if(fp==NULL) 
   {
      printf("\nError: cannot open %s\n\n",file);
      delete GN; 
      delete EN; 
      delete VL2; 
      delete VL4; 
      delete GL;
      return(1);
   }

   // Skip over the 128 byte File Preamble 
   if( fseek(fp,128,SEEK_SET) == -1) 
   {
      printf("\n128 byte File Preamble fseek() failed for %s\n\n",file);
      fclose(fp);
      delete GN; 
      delete EN; 
      delete VL2; 
      delete VL4; 
      delete GL;
      return(1);
   }

   // The four byte DICOM Prefix must contain the character string "DICM"
   if( fread(dicomPrefix,1,4,fp) != 4) 
   {
      printf("\nDICOM Prefix fread() failed for %s\n\n",file);
      fclose(fp);
      delete GN; 
      delete EN; 
      delete VL2; 
      delete VL4; 
      delete GL;
      return(1);
   }

   // we take this condition to mean that the File Meta Information header does not exist
   if(dicomPrefix[0]!='D' || dicomPrefix[1]!='I' || dicomPrefix[2]!='C' || dicomPrefix[3]!='M')
   {
      *offset=0;
      *syntax=IMPLICIT_LITTLE_ENDIAN;		// This is DICOM default
      fclose(fp);
      delete GN; 
      delete EN; 
      delete VL2; 
      delete VL4; 
      delete GL;
      return(1);
   }

   //////////////////////////////////////////////////////////////////////
   // read Tag: (0002,0000)	Group Length
   errorFlag=0;

   if( fread(GN, 2, 1, fp) != 1)  errorFlag=1;
   if( bigEndian() ) swapByteOrder( (char *)GN, 2);

   if( fread(EN, 2, 1, fp) != 1)  errorFlag=1;
   if( bigEndian() ) swapByteOrder( (char *)EN, 2);

   if( fread(VR, 1, 2, fp) != 2) errorFlag=1;
   VR[2]='\0';

   if( fread(VL2, 2, 1, fp) != 1) errorFlag=1;
   if( bigEndian() ) swapByteOrder( (char *)VL2, 2);

   if( fread(GL, 4, 1, fp) != 1) errorFlag=1;
   if( bigEndian() ) swapByteOrder( (char *)GL, 4);

   if(*GN!=0x0002 || *EN!=0x0000 || errorFlag==1)
   {
      printf("\nError reading Tag:(0002,0000) from %s\n\n",file);
      fclose(fp);
      delete GN; 
      delete EN; 
      delete VL2; 
      delete VL4; 
      delete GL;
      return(1);
   }

   //printf("GN=%04x EN=%04x VR=%s VL=%d GL=%d\n",*GN,*EN,VR,*VL2,*GL);

   *offset = *GL + ftell(fp);
   //////////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////////
   // read Tag: (0002,0001)	File Meta Information Version
   errorFlag=0;

   if( fread(GN, 2, 1, fp) != 1)  errorFlag=1;
   if( bigEndian() ) swapByteOrder( (char *)GN, 2);

   if( fread(EN, 2, 1, fp) != 1)  errorFlag=1;
   if( bigEndian() ) swapByteOrder( (char *)EN, 2);

   if( fread(VR, 1, 2, fp) != 2) errorFlag=1;
   VR[2]='\0';

   // skip the two reserved bytes, see Table 7.1-1, PS 3.5-1998, Page 27
   if( fseek(fp,2,SEEK_CUR) == -1) errorFlag=1;

   if( fread(VL4, 4, 1, fp) != 1) errorFlag=1;
   if( bigEndian() ) swapByteOrder( (char *)VL4, 4);

   // skip over the value of this tag
   if(*VL4 != 0xFFFFFFFF )
   if( fseek(fp,*VL4,SEEK_CUR) == -1) errorFlag=1;

   if(*GN!=0x0002 || *EN!=0x0001 || errorFlag==1)
   {
      printf("\nError reading Tag:(0002,0001) from %s\n\n",file);
      fclose(fp);
      delete GN; 
      delete EN; 
      delete VL2; 
      delete VL4; 
      delete GL;
      return(1);
   }
   //printf("GN=%04x EN=%04x VR=%s VL=%08x\n",*GN,*EN,VR,*VL4);
   //////////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////////
   // read Tag: (0002,0002)  Media Storage SOP Class UID	
   errorFlag=0;

   if( fread(GN, 2, 1, fp) != 1)  errorFlag=1;
   if( bigEndian() ) swapByteOrder( (char *)GN, 2);

   if( fread(EN, 2, 1, fp) != 1)  errorFlag=1;
   if( bigEndian() ) swapByteOrder( (char *)EN, 2);

   if( fread(VR, 1, 2, fp) != 2) errorFlag=1;
   VR[2]='\0';

   if( fread(VL2, 2, 1, fp) != 1) errorFlag=1;
   if( bigEndian() ) swapByteOrder( (char *)VL2, 2);

   // skip over the value of this tag
   if(*VL2>0)
   if( fseek(fp,*VL2,SEEK_CUR) == -1) errorFlag=1;
   
   if(*GN!=0x0002 || *EN!=0x0002 || errorFlag==1)
   {
      printf("\nError reading Tag:(0002,0002) from %s\n\n",file);
      fclose(fp);
      delete GN; 
      delete EN; 
      delete VL2; 
      delete VL4; 
      delete GL;
      return(1);
   }
   //printf("GN=%04x EN=%04x VR=%s VL=%d\n",*GN,*EN,VR,*VL2);
   //////////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////////
   // read Tag: (0002,0003)  Media Storage SOP Instance UID	
   errorFlag=0;

   if( fread(GN, 2, 1, fp) != 1)  errorFlag=1;
   if( bigEndian() ) swapByteOrder( (char *)GN, 2);

   if( fread(EN, 2, 1, fp) != 1)  errorFlag=1;
   if( bigEndian() ) swapByteOrder( (char *)EN, 2);

   if( fread(VR, 1, 2, fp) != 2) errorFlag=1;
   VR[2]='\0';

   if( fread(VL2, 2, 1, fp) != 1) errorFlag=1;
   if( bigEndian() ) swapByteOrder( (char *)VL2, 2);

   // skip over the value of this tag
   if(*VL2>0)
   if( fseek(fp,*VL2,SEEK_CUR) == -1) errorFlag=1;

   if(*GN!=0x0002 || *EN!=0x0003 || errorFlag==1)
   {
      printf("\nError reading Tag:(0002,0003) from %s\n\n",file);
      fclose(fp);
      delete GN; 
      delete EN; 
      delete VL2; 
      delete VL4; 
      delete GL;
      return(1);
   }
   //printf("GN=%04x EN=%04x VR=%s VL=%d\n",*GN,*EN,VR,*VL2);
   //////////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////////
   // read Tag: (0002,0010)  Transfer Syntax UID
   errorFlag=0;

   if( fread(GN, 2, 1, fp) != 1)  errorFlag=1;
   if( bigEndian() ) swapByteOrder( (char *)GN, 2);

   if( fread(EN, 2, 1, fp) != 1)  errorFlag=1;
   if( bigEndian() ) swapByteOrder( (char *)EN, 2);

   if( fread(VR, 1, 2, fp) != 2) errorFlag=1;
   VR[2]='\0';

   if( fread(VL2, 2, 1, fp) != 1) errorFlag=1;
   if( bigEndian() ) swapByteOrder( (char *)VL2, 2);

   // skip over the value of this tag
   if( fread(UI,1,*VL2,fp) != *VL2) errorFlag=1;
   UI[*VL2]='\0';

   if(*GN!=0x0002 || *EN!=0x0010 || errorFlag==1)
   {
      printf("\nError reading Tag:(0002,0010) from %s\n\n",file);
      fclose(fp);
      delete GN; 
      delete EN; 
      delete VL2; 
      delete VL4; 
      delete GL;
      return(1);
   }
   // printf("GN=%04x EN=%04x VR=%s VL=%d UI=%s\n",*GN,*EN,VR,*VL2,UI);
   //////////////////////////////////////////////////////////////////////

   fclose(fp);

   if( strcmp(UI, "1.2.840.10008.1.2.1") == 0)
      *syntax=EXPLICIT_LITTLE_ENDIAN;
   else if( strcmp(UI, "1.2.840.10008.1.2.2") == 0)
      *syntax=EXPLICIT_BIG_ENDIAN;
   else
   {
      printf("\nError: Transfer Syntax not recognized in %s\n\n",file);
      delete GN; 
      delete EN; 
      delete VL2; 
      delete VL4; 
      delete GL;
      return(1);
   }

   delete GN; 
   delete EN; 
   delete VL2; 
   delete VL4; 
   delete GL;

   return(0);
}
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
// oVL is output value lenght
// iGN is input group number
// iEN is input element number
// oV is output value
///////////////////////////////////////////////////////////////////////////////////////////////
int readTag(const char *file, unsigned short iGN, unsigned short iEN, long byteOffset, int transferSyntax,
unsigned long *oVL, char *oV, long *valueOffset)
{
	int errorFlag=0;
	FILE *fp;
	char VR[3]; // DICOM Value Representation	see: PS 3.5-1998 Page 15 and 27
	unsigned short *GN = new unsigned short; // group number
	unsigned short *EN = new unsigned short; // element number
	unsigned short *VL2 = new unsigned short; // Value Length see: PS 3.5-1998 Page 27
	unsigned long *VL4 = new unsigned long; // Value Length see: PS 3.5-1998 Page 27

	fp=fopen(file,"r");
	if(fp==NULL)
	{
		printf("\nError: cannot open %s\n\n",file);
		delete GN; delete EN; delete VL2; delete VL4;
		return(1);
	}

	if( fseek(fp,byteOffset,SEEK_SET)==-1)
	{
		printf("\n%ld byte fseek() failed for %s\n\n",byteOffset,file);
		fclose(fp);
		delete GN; delete EN; delete VL2; delete VL4;
		return(1);
	}

	do 
	{
		if( fread(GN, 2, 1, fp) != 1)  errorFlag=1;
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)GN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&GN, 2);

		if( fread(EN, 2, 1, fp) != 1)  errorFlag=1;
		if( bigEndian() && transferSyntax!=EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)EN, 2);
        if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN ) swapByteOrder( (char *)&EN, 2);

		//printf("GN=%04x EN=%04x\n",*GN,*EN);

		if(transferSyntax==IMPLICIT_LITTLE_ENDIAN)
		{
			if( fread(VL4, 4, 1, fp) != 1) errorFlag=1;
			if( bigEndian() ) swapByteOrder( (char *)VL4, 4);

			*oVL = *VL4;

			// printf("\tVL=%08x\n",*VL4);

			if(*GN==iGN && *EN==iEN) break;

			// skip over the value of this tag
			if(*VL4 != 0xFFFFFFFF )
			if( fseek(fp,*VL4,SEEK_CUR) == -1) errorFlag=1;
		}
		else if( transferSyntax==EXPLICIT_LITTLE_ENDIAN || transferSyntax==EXPLICIT_BIG_ENDIAN)
		{
			if(*GN==0xFFFE && (*EN==0xE000 || 0xE0DD) )
			{
				if( fread(VL4, 4, 1, fp) != 1) errorFlag=1;

				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)VL4, 4);
				if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)VL4, 4);

				if( fseek(fp,*VL4,SEEK_CUR) == -1) errorFlag=1;
				continue;
			}

			if( fread(VR, 1, 2, fp) != 2) errorFlag=1;
			VR[2]='\0';
			//printf("\tVR=%s \n",VR);

			if(strcmp(VR,"OB")==0 || strcmp(VR,"OW")==0 || strcmp(VR,"SQ")==0 || strcmp(VR,"UN")==0)
			{
				// skip the two reserved bytes, see Table 7.1-1, PS 3.5-1998, Page 27
				if( fseek(fp,2,SEEK_CUR) == -1) errorFlag=1;

				if( fread(VL4, 4, 1, fp) != 1) errorFlag=1;
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)VL4, 4);
				if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)VL4, 4);

				*oVL = *VL4;

				//printf("\tVL4=%08x\n",*VL4);

				if(*GN==iGN && *EN==iEN) break;

				// skip over the value of this tag
				if(*VL4 != 0xFFFFFFFF )
				if( fseek(fp,*VL4,SEEK_CUR) == -1) errorFlag=1;
			}
			else
			{
				if( fread(VL2, 2, 1, fp) != 1) errorFlag=1;
				if( bigEndian() && transferSyntax==EXPLICIT_LITTLE_ENDIAN) swapByteOrder( (char *)VL2, 2);
				if( !bigEndian() && transferSyntax==EXPLICIT_BIG_ENDIAN) swapByteOrder( (char *)VL2, 2);

				*oVL = *VL2;

				//printf("\tVL2=%d\n",*VL2);

				if(*GN==iGN && *EN==iEN) break;

				// skip over the value of this tag
				if(*VL2>0)
				if( fseek(fp,*VL2,SEEK_CUR) == -1) errorFlag=1;
			}
		}
	}
	while( !errorFlag && (iGN != *GN || iEN != *EN) );

	*valueOffset = ftell(fp);

	// There was a problem here. In some example images *oVL was 100 and the
	// size of oV was 64.  So we changed the size of oV to 256 in 
	// the calling function.
	if( fread(oV, 1, *oVL, fp) != *oVL) errorFlag=1;

	fclose(fp);


	delete GN; delete EN; delete VL2; delete VL4;

	if(errorFlag)
		return(1);
	else
		return(0);
}

///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
// VL: value length	(see PS 3.5-1988, P. 27, Table 7-1-3)
// V: element value
int read_element(char *filename, short S_GN, short S_EN, char *V, int *VL)
{
	FILE *fp;
	short GN;	// group number
	short EN;	// element number
	
	int i;

	fp=fopen(filename,"r");
	if(fp==NULL) return(0);

	*VL=0;
	i=0;
	do {
		i++;

		if( fread (&GN, 2, 1, fp) == 1 ) {};
		if( bigEndian() ) swapByteOrder( (char *)(&GN), 2);

		if( fread (&EN, 1, 2, fp) == 2 ) {};
		if( bigEndian() ) swapByteOrder( (char *)(&EN), 2);

		if ( fread (VL, 1, 4, fp) == 4 ) {};
		if( bigEndian() ) swapByteOrder( (char *)VL, 4);

		if(EN != S_EN || GN != S_GN)
			fseek(fp,*VL,1);
		else
			if ( fread (V, 1, *VL, fp) == (size_t)(*VL) ) {};
      
		// printf("Group=%04x   Element Number=%04x   Element Length=%4d\n",GN,EN,*VL); 

	} while( i<1000 && !(GN==0x7fe0 && EN==0x10) && (EN != S_EN || GN != S_GN) );

	fclose(fp);   

	if(EN != S_EN || GN != S_GN)
	{
		V=NULL;
		return(0);
	}

	return(1);
}

void readMatrixSize(char *filename, int *nx, int *ny)
{
	char US[2];
	int VL;

	read_element(filename,0x28,0x11,US,&VL);
	if( bigEndian() ) swapByteOrder(US,2);
	*nx=*(short *)US;

	read_element(filename,0x28,0x10,US,&VL);
	if( bigEndian() ) swapByteOrder(US,2);
	*ny=*(short *)US;

	return;
}

void readVoxelSize(char *filename, float *dx, float *dy, float *dz)
{
	char DS[256];
	char dum[32];
	int L1,L2;
	int VL;		// value length

	read_element(filename,0x28,0x30,DS,&VL);
	L1=strcspn(DS,"\\");
	strncpy(dum,DS,L1); dum[L1]='\0';
	*dx=atof(dum);

	strcpy(dum,DS+L1+1);
	L2=strlen(DS+L1+1);
	dum[L2]='\0';
	*dy=atof(dum);

	if( read_element(filename,0x18,0x88,DS,&VL) )
	{
		DS[VL]='\0';
		*dz = atof(DS);
	}
	else // Some 3D sequences may not have element 88 of group 18.
	{
		read_element(filename,0x18,0x50,DS,&VL);
		DS[VL]='\0';
		*dz = atof(DS);
	}
}
