///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// Copyright (C) 2024 Babak A. Ardekani, PhD - All Rights Reserved.  //
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

// This is The INTERFACE for the class NIFTIIMAGE

#ifndef NIFTIIMAGE_H 
#define NIFTIIMAGE_H

#include <nifti1.h>

class NIFTIIMAGE
{
public:
   NIFTIIMAGE();	// constructor
   NIFTIIMAGE(int nx, int ny, int nz, float dx, float dy, float dz, short datatype);   // constructor
   NIFTIIMAGE(int nx, int ny, int nz, int nt, float dx, float dy, float dz, short datatype);   // constructor

   ~NIFTIIMAGE();	// destructor

   void printheader();
   void read(const char *filename);
   void readheader(const char *filename);
   void write(const char *filename);
   void setheader(nifti_1_header newhdr);
   nifti_1_header getheader();
   int datasize();
   int nv();
   int nx();
   int ny();
   int nz();
   float dx();
   float dy();
   float dz();
   short datatype();
   char *getdata();

private:
   nifti_1_header hdr;
   nifti1_extender ext;
   nifti1_extension extension;
   char *data;
};

#endif // NIFTIIMAGE_H
