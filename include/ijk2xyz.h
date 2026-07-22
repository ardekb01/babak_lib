#ifndef IJK2XYZ_H 
#define IJK2XYZ_H 

// Convert P from (i, j, k) to (x, y, z) coordinates.
// P is a 3 x n matrix. Each column is a point in the
// (i, j, k) coordinate system. Each point is converted
// in-place to the (x, y, z) coordinate system.
bool convert_to_xyz(float *P, int n, SHORTIM im);

// Convert P from (x, y, z) to (i, j, k) coordinates.
// P is a 3 x n matrix. Each column is a point in the
// (x, y, z) coordinate system. Each point is converted
// in-place to the (i, j, k) coordinate system.
// Core operation: P -> point -> matrix multiplication -> result -> P
bool convert_to_ijk(float *P, int n, SHORTIM im);

bool xyz2ijk(float *T,
             int nx,
             int ny,
             int nz,
             float dx,
             float dy,
             float dz);

bool ijk2xyz(float *T,
             int nx,
             int ny,
             int nz,
             float dx,
             float dy,
             float dz);

bool xyz2ijk(float *T, const DIM &dim);
bool ijk2xyz(float *T, const DIM &dim);

#endif
