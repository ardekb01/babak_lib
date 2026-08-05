#ifndef INTENSITY_WEIGHTED_CENTROID_H 
#define INTENSITY_WEIGHTED_CENTROID_H

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
