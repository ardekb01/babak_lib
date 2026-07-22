#ifndef IJK2XYZ_H 
#define IJK2XYZ_H 

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
