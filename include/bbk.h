#ifndef BBK_H 
#define BBK_H

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

bool inversePILtransform(const char *orientCode, float *T);
bool PILtransform(const char *orientCode, float *T);

// orientation must point to a buffer of at least 4 bytes.
bool getNiftiImageOrientation(const char *filename,
                              char *orientation);

// orientation must point to a buffer of at least 4 bytes.
bool getNiftiImageOrientation(nifti_1_header hdr,
                              char *orientation);

bool get_directory_name(const char *pathname, char *dirname, size_t dirnameSize);

bool check_nifti_file_extension(const char *filename);

bool check_nifti1_magic(const char *imagefilename);

bool get_nifti_basename(char *filename,
                        size_t filenameSize,
                        const char *path);


bool valid_orientation_code(const char *orientCode);

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
