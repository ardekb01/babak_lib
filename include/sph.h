#ifndef _sph_h

class SPH 
{
   public:
   int *i;   // i voxel coordinates
   int *j;   // j voxel coordinates
   int *k;   // k voxel coordinates
   int n;    // size in voxels
   int r;    // radius in voxels
   float *v; // values at (i,j,k) - can be set by set function
   SPH(int radius);
   void reset();
   void set(SHORTIM im, int ic, int jc, int kc);
   void set(float *im, int nx, int ny, int nz, int ic, int jc, int kc);

   void get(float *array);

   // Subtracts the mean value of array v from itself
   // v[i] = v[i] - mean(v[.])
   // The resulting v[.] will have zero mean.
   void zeromean();

   // returns the mean value of array v[.]
   float mean();

   // returns the L2 norm of array v[.]
   float norm();

   void normalize();
   float dot(float *u);
   ~SPH();
};

#define _sph_h

#endif
