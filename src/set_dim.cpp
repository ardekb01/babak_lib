///////////////////////////////////////////////////////////////////////
// Copyright (C) 2024 Babak A. Ardekani, PhD - All Rights Reserved.
///////////////////////////////////////////////////////////////////////

#include <babak_lib.h>

void set_dim(SHORTIM &im, const SHORTIM &sourceim)
{
   im.nx = sourceim.nx;
   im.ny = sourceim.ny;
   im.nz = sourceim.nz;
   im.nt = sourceim.nt;

   im.dx = sourceim.dx;
   im.dy = sourceim.dy;
   im.dz = sourceim.dz;

   im.np = sourceim.np;
   im.nv = sourceim.nv;
}

void set_dim(DIM &dim, const SHORTIM &im)
{
   dim.nx = im.nx;
   dim.ny = im.ny;
   dim.nz = im.nz;
   dim.nt = im.nt;

   dim.dx = im.dx;
   dim.dy = im.dy;
   dim.dz = im.dz;

   dim.np = im.np;
   dim.nv = im.nv;
}

void set_dim(SHORTIM &im, const DIM &dim)
{
   im.nx = dim.nx;
   im.ny = dim.ny;
   im.nz = dim.nz;
   im.nt = dim.nt;

   im.dx = dim.dx;
   im.dy = dim.dy;
   im.dz = dim.dz;

   im.np = dim.np;
   im.nv = dim.nv;
}

void set_dim(SHORTIM &im, const nifti_1_header &hdr)
{
   im.nx = hdr.dim[1];
   im.ny = hdr.dim[2];
   im.nz = hdr.dim[3];
   im.nt = hdr.dim[4];

   im.dx = hdr.pixdim[1];
   im.dy = hdr.pixdim[2];
   im.dz = hdr.pixdim[3];

   im.np = im.nx * im.ny;
   im.nv = im.np * im.nz;
}

void set_dim(nifti_1_header &hdr, const DIM &dim)
{
   hdr.dim[0] = 3;
   hdr.dim[1] = dim.nx;
   hdr.dim[2] = dim.ny;
   hdr.dim[3] = dim.nz;

   hdr.pixdim[1] = dim.dx;
   hdr.pixdim[2] = dim.dy;
   hdr.pixdim[3] = dim.dz;
}

void set_dim(DIM &dim,
             int nx,
             int ny,
             int nz,
             float dx,
             float dy,
             float dz)
{
   dim.nx = nx;
   dim.ny = ny;
   dim.nz = nz;

   dim.dx = dx;
   dim.dy = dy;
   dim.dz = dz;

   dim.np = nx * ny;
   dim.nv = nx * ny * nz;
}

void set_dim(DIM &dim, const nifti_1_header *hdr)
{
   if (hdr == nullptr)
      return;

   dim.nx = hdr->dim[1];
   dim.ny = hdr->dim[2];
   dim.nz = hdr->dim[3];
   dim.nt = hdr->dim[4];

   dim.dx = hdr->pixdim[1];
   dim.dy = hdr->pixdim[2];
   dim.dz = hdr->pixdim[3];

   dim.np = dim.nx * dim.ny;
   dim.nv = dim.np * dim.nz;
}

void set_dim(DIM &dim, const nifti_1_header hdr)
{
   dim.nx = hdr.dim[1]; 
   dim.ny = hdr.dim[2]; 
   dim.nz = hdr.dim[3];
   dim.nt = hdr.dim[4];
   dim.dx = hdr.pixdim[1]; 
   dim.dy = hdr.pixdim[2]; 
   dim.dz = hdr.pixdim[3];
   dim.np = hdr.dim[1] * hdr.dim[2];
   dim.nv = hdr.dim[1] * hdr.dim[2] * hdr.dim[3];
}
