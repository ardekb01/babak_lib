#ifndef BBK_H 
#define BBK_H

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
