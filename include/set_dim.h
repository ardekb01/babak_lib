#ifndef SET_DIM_H
#define SET_DIM_H

void set_dim(SHORTIM &im, const SHORTIM &sourceim);
void set_dim(DIM &dim, const SHORTIM &im);
void set_dim(SHORTIM &im, const DIM &dim);
void set_dim(SHORTIM &im, const nifti_1_header &hdr);
void set_dim(nifti_1_header &hdr, const DIM &dim);
void set_dim(DIM &dim,
             int nx,
             int ny,
             int nz,
             float dx,
             float dy,
             float dz);
void set_dim(DIM &dim, const nifti_1_header *hdr);
void set_dim(DIM &dim, const nifti_1_header hdr);

#endif
