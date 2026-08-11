#ifndef BBK_H 
#define BBK_H

#include "nifti1.h"

// Definitions associated with loadTransformation()
#define LOADTRANSFORM_OK 0
#define LOADTRANSFORM_NULL_POINTER 1
#define LOADTRANSFORM_FILE_ERROR 2
#define LOADTRANSFORM_PARSE_ERROR 3

// Global variables associated with getoption()
extern int optInd;
extern const char *optArg;

struct model_file_hdr
{
   int nxHR;
   int nzHR;
   float dxHR;
   int nxLR;
   int nvol; // number of image volumes in the training set 
   int RPtemplateradius;
   int RPtemplateheight;
   int RPtemplatesize;
   int ACtemplateradius;
   int ACtemplateheight;
   int ACtemplatesize;
   int PCtemplateradius;
   int PCtemplateheight;
   int PCtemplatesize;
   int nangles; // number of angles, each template is rotated by this many angles and saved
};

typedef struct model_file_hdr model_file_hdr;

struct model_file_tail
{
   float RPPCmean[2]; // RPPC is a vector on the MSP that points from the RP point to the PC.   RP------->PC
   float parcomMean;
   float percomMean; 
   float RPmean[2];
};

typedef struct model_file_tail model_file_tail;

typedef struct shortim
{
   int nx;
   int ny;
   int nz;
   int nt;
   int np;
   int nv;
   float dx;
   float dy;
   float dz;
   short *v; // image values
} SHORTIM;

typedef struct dim {
   int nx; // number of columns 
   int ny; // number of rows
   int nz; // number of slices
   int nt; // number of frames (epochs)
   int np; // nx*ny
   int nv; // nx*ny*nz
   float dx; // x voxel dimension (mm)
   float dy; // y voxel dimension (mm)
   float dz; // z voxel dimension (mm)
   float dt; // t between frames (sec)
} DIM;

// Structure associated with getoption()
struct CmdOption
{
    const char *name;
    int has_arg;
    int val;
};

char *circularMask(short *mask, 
                   const DIM &dim, 
                   const float *center, 
                   const double radius);

bool initialMSPestimation(const short *image, 
                          const DIM &dim, 
                          float &A, 
                          float &B, 
                          float &C);

int getoption(int argc, char *const argv[], const struct CmdOption *options);

int loadTransformation(const char *filename, float *T);

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

bool vectorNorm(const float *x, int n, double &norm);
bool normalizeVector(float *x, int n);

bool MSPtransformation(const char *filename, const char *orient, float *Tmsp, DIM &dim);
bool MSPtransformation(const char *filename, const char *orient, const char *lmfile, float *Tmsp, DIM &dim);

char directionCode(float x, float y, float z);

bool bigEndian();
bool swapByteOrder(char *in, size_t N);
bool swapN(char *in, size_t N);
bool swap_float_array(float *x, size_t n);
bool swap_double_array(double *x, size_t n);
bool swap_int_array(int *x, size_t n); 
void swap_model_file_hdr(model_file_hdr *hdr);
void swap_model_file_tail(model_file_tail *tail);
void swapniftiheader(nifti_1_header *hdr);

bool rotate(float *R, float alpha, float x, float y, float z);

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

char *read_nifti_image(const char *filename, nifti_1_header *hdr);

bool read_nifti_hdr(const char *filename, nifti_1_header *hdr);

void setLowHigh(const short *image, int nv, int &low, int &high, float percent);

float *gaussian_kernel(const float sd, int *n);

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
