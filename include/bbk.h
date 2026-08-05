#ifndef BBK_H 
#define BBK_H

float symm_objective_func(const short *image, const DIM &dim, float A, float B, float C);

bool intensity_weighted_centroid(
   const short *image,
   int nx,
   int ny,
   int nz,
   float dx,
   float dy,
   float dz,
   float &x,
   float &y,
   float &z);

#endif
