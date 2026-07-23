#ifndef REORIENTVOLUME_H
#define REORIENTVOLUME_H

// The upper-left 3x3 submatrix of T is a signed permutation matrix.
// Callers must release the returned image with free()
short *reorientVolume(short *input_image,
                       int nx1,
                       int ny1,
                       int nz1,
                       float dx1,
                       float dy1,
                       float dz1,
                       float *T,
                       int &nx2,
                       int &ny2,
                       int &nz2,
                       float &dx2,
                       float &dy2,
                       float &dz2);

#endif
