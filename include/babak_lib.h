#ifndef _babak_lib_h
#define _babak_lib_h

#include <nifti1.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef _PILTRANSFORM
extern double searchradius[3]; // in units of mm
extern int opt_CENTER_AC;
#endif

typedef unsigned char uchar;
typedef unsigned short uint2;
typedef unsigned short uint2;
typedef unsigned int uint4;
typedef short int2;
typedef int int4;
typedef float float4;
typedef double float8;

#define YES 1
#define NO 0
#define ESMALL 1e-10
#define DEFAULT_STRING_LENGTH 512
#define MAXIM 256 // maximum number of images allowed 

#ifndef _TPS
extern float4 tpsU(float4 x, float4 y);
extern float4 *tpsK(float4 *x, float4 *y, int n);
#endif

#ifndef _update_qsform
void update_qsform(nifti_1_header &hdr, float *matrix);
void update_qsform( const char *imagefilename , float *matrix);
void update_qsform( const char *imagefile1, const char *imagefile2);
extern int opt_sform;
extern int opt_qform;
#endif

#ifndef _getARTHOME
extern char *ARTHOME;
extern void getARTHOME();
#endif

#ifndef _artlib
// Turns verbose mode on (opt_v=YES) or off (opt_v=NO)
// Initialized to NO in artlib.c
extern char opt_v;
#endif

#ifndef _SHORTIM
#define _SHORTIM
typedef struct shortim
{
   int4 nx;
   int4 ny;
   int4 nz;
   int4 nt;
   int4 np;
   int4 nv;
   float4 dx;
   float4 dy;
   float4 dz;
   int2 *v; // image values
} SHORTIM;
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef _DIM
#define _DIM
typedef struct dim {
   int4 nx; // number of columns 
   int4 ny; // number of rows
   int4 nz; // number of slices
   int4 nt; // number of frames (epochs)
   int4 np; // nx*ny
   int4 nv; // nx*ny*nz
   float4 dx; // x voxel dimension (mm)
   float4 dy; // y voxel dimension (mm)
   float4 dz; // z voxel dimension (mm)
   float4 dt; // t between frames (sec)
} DIM;
#endif
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
struct DICOM_file_meta_info
{
   // maximum bytes for VR=UI is 64,  1 byte is added for the string terminator '\0'
   char media_storage_SOP_class[65];  // Tag (0002,0002)
   char transfer_syntax[65]; // Tag (0002,0010)
   int4 dataset_offset;
};
typedef struct DICOM_file_meta_info DICOM_file_meta_info;

struct DICOM_hdr
{
   char patientID[65];    // Tag (0010,0020)  VR=LO (Long String 64 chars max.)
   float4 slice_thickness; // Tag (0018,0050)
   int4 TR;             // Tag (018,0080)
   int4 TE;             // Tag (018,0081)
   int4 TI;             // Tag (018,0082)
   int4 NEX;            // Tag (018,0083)
   float4 dz;           // Tag (018,0088)
   float4 flip_angle;   // Tag (018,1314)
   char seriesID[65];  // Tag (020,000E)
   int4 series_number;  // Tag (020,0011)
   float4 TLHC[3];      // Tag (020,0032)
   float4 rowvec[3];    // Tag (020,0037)
   float4 colvec[3];    // Tag (020,0037)
   char frame_of_referenceID[65]; // Tag (020,0052)
   int4 nz; // Tag (020,1002)
   float4 slice_location; // Tag (020,1041)
   int4 ny; // Tag (028,0010)
   int4 nx; // Tag (028,0011)
   float4 dy; // Tag (028,0030)
   float4 dx; // Tag (028,0030)
};
typedef struct DICOM_hdr DICOM_hdr;

//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
struct model_file_hdr
{
   int4 nxHR;
   int4 nzHR;
   float4 dxHR;
   int4 nxLR;
   int4 nvol; // number of image volumes in the training set 
   int4 RPtemplateradius;
   int4 RPtemplateheight;
   int4 RPtemplatesize;
   int4 ACtemplateradius;
   int4 ACtemplateheight;
   int4 ACtemplatesize;
   int4 PCtemplateradius;
   int4 PCtemplateheight;
   int4 PCtemplatesize;
   int4 nangles; // number of angles, each template is rotated by this many angles and saved
};

typedef struct model_file_hdr model_file_hdr;

struct model_file_tail
{
   float4 RPPCmean[2]; // RPPC is a vector on the MSP that points from the RP point to the PC.   RP------->PC
   float4 parcomMean;
   float4 percomMean; 
   float4 RPmean[2];
};

typedef struct model_file_tail model_file_tail;
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
struct im_params {

   int4 MRimageId;
   int4 PETimageId;

	/* maximum values of the image 1 */	
	uchar max1; 	

	int4 	*Size;   /* size of non-zero CC's in CCI */
	float4 *S,*SS;  

	int4 	FROM;
	int4 	TO;

	/* number of image 1 columns, rows and slices  */
	int4	nx1,ny1,nz1; 

	/* number of image 1 voxels, number of image 1 pixels per slice */
	int4	nv1,np1;      	

	float4  dx1,dy1,dz1; 	/* image 1 voxel size (mm) */

	/* coordinates of image 1 volume center wrt pixel (0,0,0) in (mm) */
	float4  xc1,yc1,zc1; 	

	// number of image 2 columns, rows and slices  
	int4	nx2,ny2,nz2; 
	float4  dx2,dy2,dz2; 	// image 2 voxel size (mm) 

	/* number of image 2 voxels, number of image 2 pixels per slice */
	int4	nv2,np2;      

	/* coordinates of image 2 volume center wrt pixel (0,0,0) in (mm) */
	float4  xc2,yc2,zc2; 	

	int4 NCC;

	int4 NP;  // number of non-zero pixels in CCI

	int4 sf;  // speed factor 

	unsigned *CCI;  	/* labeled connected compontents' image */

	uchar *data1;  /* image 1 */
	uchar *data2;  /* image 2 */
};
//////////////////////////////////////////////////////////////////////////////////////////////////

struct dicominfo
{
   char patientName[321];
   char DOB[11];
   char patientID[65];
   char studyDate[11];
   char TE[17];
   char TR[17];
   char TI[17];
   char ETL[13];
   char NEX[17];
   char flipAngle[17];
   char bandwidth[17];
   char freq[17];
   uint2 acquisitionMatrix[4];
   char phaseFOV[17];
   char dum[3];
};
typedef struct dicominfo dicominfo;
//////////////////////////////////////////////////////////////////////////////////////////////////

void leastSquaresAffineTrans(float *P, float *Q, int n, float *A);

void hist1D_plot(const char *name, int n, int *bin, float *data1, float *data2);
void hist1D_plot(const char *name, int n, int *bin, float *data1, float *data2, int max_x);
void hist1D_plot(const char *name, int n, int *bin, float *data1, float *data2, int max_x, int thr);

// The following functions are defined in errorMessage.cxx
void memory_allocation_error(const char *variablename);
void file_open_error(const char *filename);
void errorMessage(const char *message);

void set_dim(nifti_1_header &hdr, DIM dim);
void set_dim(DIM &dim, int4 nx, int4 ny, int4 nz, float4 dx, float4 dy, float4 dz);
void set_dim(DIM &dim, nifti_1_header hdr);
void set_dim(DIM &dim, nifti_1_header *hdr);
void set_dim(SHORTIM &im, nifti_1_header hdr);
void set_dim(SHORTIM &im, DIM dim);
void set_dim(DIM &dim, SHORTIM im);
void set_dim(SHORTIM &im, SHORTIM sourceim);
void findMSP(const char *filename, char *orient, const char *lmfile, float4 *Tmsp, int4 verbose, DIM &dim);
// Input: (x,y,z) a vector defined in RAS system
// Output: One of six charaters {R,L,A,P,S,I}
char directionCode(float4 x, float4 y, float4 z);
void getNiftiImageOrientation(const char *filename, char *orientation);
void getNiftiImageOrientation(nifti_1_header hdr, char *orientation);
int4 checkNiftiFileExtension(const char *filename);
int4 isOrientationCodeValid(const char *orientCode);
void swap_model_file_hdr(model_file_hdr *hdr);
void swap_model_file_tail(model_file_tail *tail);
void new_PIL_transform(const char *subfile, const char *lmfile, char *orient, float4 *T, int SAVE_PIL_FLAG);
void standard_PIL_transformation(const char *imfile, const char *lmfile, char *orient, int4 verbose, float4 *TPIL);
void convert_to_xyz(float *P, int n, SHORTIM im);
void Procrustes(float *Q, float *Qavg, int n, float *P, float *Pavg, float *TLM);
void Procrustes(float *Q, int n, float *P, float *TLM);
void PILtransform(const char *orientCode, float4 *orientMat);
void inversePILtransform(const char *orientCode, float4 *orientMat);
int2 *reorientVolume(int2 *v1, int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1, float4 *orientMat,
int4 &nx2, int4 &ny2, int4 &nz2, float4 &dx2, float4 &dy2, float4 &dz2);
int2 *reorientVolume(int2 *input_image, nifti_1_header oldhdr, const char *neworient, nifti_1_header &newhdr, 
float4 *T_oldorient_to_neworient);
void rotate(float4 *T, float4 alpha, float4 x, float4 y, float4 z);
float4 *rotate(float4 alpha, float4 x, float4 y, float4 z);
void rotate(float4 *T, float4 CosAlpha, float4 SinAlpha, float4 x, float4 y, float4 z);
void compute_cm(int2 *image, int4 nx, int4 ny, int4 nz, float4 dx, float4 dy, float4 dz, float4 *x, float4 *y, float4 *z);
void compute_cm(short *image, DIM dim, float *cm);
void compute_cov(short *image, DIM dim, float *cm, double *I);
void standardize(float4 *x, int4 n);
void standardize(float *x, float *mask, int n);
void irodrigues_formula(float4 *R, float4 *w, float4 &theta);
void rodrigues_formula(float4 *R, float4 *w, float4 theta);
void rodrigues_formula4x4(float4 *R, float4 *w, float4 theta);
void se3_to_SE3(float4 *M, float4 *w, float4 *v, float4 theta);
void SE3_to_se3(float4 *M, float4 *w, float4 *v, float4 &theta);

////////////////////////////////////////////////////////////////////////////////
//defined in setLowHigh.cxx
void setLowHigh(int2 *image, int4 nv, int4 *low, int4 *high);	
void setLowHigh(int2 *image, int4 nv, int4 *low, int4 *high, float4 percent);
void setMX(int2 *image, int2 *msk, int4 nv, int4 &high, float4 percent);
////////////////////////////////////////////////////////////////////////////////

#ifndef _singular_value_decomposition
extern int4 Singular_Value_Decomposition(float8* A, int4 nrows, int4 ncols, float8* U, 
                      float8* singular_values, float8* V, float8* dummy_array);
extern void Singular_Value_Decomposition_Inverse(float8* U, float8* D, float8* V,  
                        float8 tolerance, int4 nrows, int4 ncols, float8 *Astar);
#endif

#ifndef _DKI
extern void formA7(float4 *A, int4 n, float4 *u1, float4 *u2, float4 *u3, float4 *b, float4 *w);

extern void restricted_tensor_model(float4 *x, float4 *y);
extern void second_and_fourth_order_moments(float4 *x, float4 *y, float4 *E, float4 *F, int4 v);
extern float4 mardia_kurtosis(float4 *x1, float4 *x2, int4 v);
extern float4 mardia_kurtosis(float4 *x);

extern float4 estimate_K(float4 *A, float4 *x, float4 *y, int4 n);
extern void update_x(float4 *A, float4 *x, float4 *y, float4 K, int4 n);
extern void update_v(float4 *u0, float4 *u1, float4 *u2, float4 *b, float4 *y, int4 n, float4 *v, float4 scale);
extern void update_v(float4 *u0, float4 *u1, float4 *u2, float4 *b, float4 *y, int4 n, float4 *v1, float4 *v2, float4 *v3, float4 scale);
extern float4 quadratic_form_1(float4 *x, float4 *v);
extern float4 quadratic_form_1(float4 x0, float4 x1, float4 x2, float4 *v);
extern void vector_to_symmetric_matrix_form(float4 *u, float4 *D);
extern void symmetric_matrix_to_vector_form(float4 *D, float4 *u);
extern void cholesky_composition_1(float4 *v, float4 *u);
extern void cholesky_composition_2(float4 *L, float4 *D);
extern int4 cholesky_decomposition_1(float4 *u, float4 *v);
extern int4 cholesky_decomposition_2(float4 *D, float4 *L);
extern void positive_definite_tensor_model(float4 *v1, float4 *v2, float4 *y);
#endif

#ifndef _artlib
extern void find_pil_transformation(char *imfile, DIM dim, float4 *pilT, float4 *AC, float4 *PC, float4 *VSPS);
extern void find_pil_transformation(char *imfile, DIM dim, float4 *pilT);
extern void update_qsform(nifti_1_header &hdr, const char *orientationcode);
extern char opt_ppm;
extern char opt_png;
extern char opt_txt;
extern char opt_AC; // if YES AC will be detected automatically
extern char opt_PC; // if YES PC will be detected automatically
extern char opt_RP; // if YES RP will be detected automatically
extern char opt_MSP; // if YES MSP will be detected automatically
extern void sub2trg_rigid_body_transformation(float4 *sub2trg, const char *subfile, const char *trgfile);
extern void getDirectoryName(const char *pathname, char *dirname);
extern void forwardTCSAP(float4 *xvec, float4 *yvec, float4 *zvec, float4 *TLHC, float4 *angle, float4 *translation, DIM dim);
extern void backwardTCSAP(float4 *xvec, float4 *yvec, float4 *angle);
extern void orig_ijk_to_pil_xyz(float4 *Tmsp, DIM orig_dim, float4 *AC, float4 *PC);
extern void initialAC(float4 Ax, float4 Ay, float4 Bx, float4 By, float4 *Cx, float4 *Cy, float4 parcomMean, float4 percomMean);
extern int2 *thresholdImageOtsu(int2 *im, int4 nv, int4 *nbv);
extern void defineTemplate(int4 r, int4 h, int2 *x, int2 *y, int2 *z);
extern char *defineACregion(DIM dim, float4 *RP, float4 *PC, float4 parcomMean, float4 percomMean, float8 ACsr);
extern char *definePCregion(DIM HR, float4 *RP, float4 *RPPCmean, float8 PCsr);
extern char *expandMask(int2 *mask_HR, DIM HR, float4 *RPmean, float8 RPsr);
extern void ACPCtransform(float4 *Tacpc, float4 *Tmsp, float4 *AC, float4 *PC, char flg);
extern void compute_MSP_parameters_from_Tmsp(float4 *Tmsp, float4 *n, float4 *d);
extern int4 detect_AC_PC_MSP( const char *imagefilename, char *orientation, char *modelfile,
float4 *AC, float4 *PC, float4 *RP, float4 *Tmsp, int4 opt_v, int4 opt_T2);
extern float4 optimizeNormalVector(int2 *image,DIM dim,float4 *A, float4 *B, float4 *C);
extern float4 reflection_cross_correlation2(int2 *image, DIM dim, float4 A, float4 B, float4 C);
extern float4 reflection_cross_correlation(int2 *image, DIM dim, float4 a, float4 b, float4 c, float4 d);
extern void findInitialNormalVector(int2 *image, DIM dim, float4 *A, float4 *B,float4 *C);
extern float4 msp(int2 *im_in, int4 nx, int4 ny, int4 nz, float4 dx, float4 dy, float4 dz, float4 *A, float4 *B, float4 *C);
extern void computeTmsp(char *orientation, int2 *volOrig, DIM dim, float4 *Tmsp);
extern int4 save_as_ppm(const char *filename, int4 nx, int4 ny, uchar *R, uchar *G, uchar *B);
extern void combine_warps_and_trans(int4 nx, int4 ny, int4 nz, float4 dx, float4 dy, float4 dz, float4 *Xwarp, float4 *Ywarp, float4 *Zwarp, float4 *T);
#endif

#ifndef _binomial
extern float8 BinomialCoeff( int4 n, int4 k );
extern void k_subset ( int4 k, float8 rank, int4 a[] );
#endif

#ifndef _maskOps
extern char *find_foreground_mask(int2 *im, int4 nv, int4 nb, int4 nclass, int4 niter, int2 *thresh);
extern int2 *mask_image(int2 *im, char *mask, int4 nv, int4 *nbv_out);
extern void cc_filter(char *mask, int4 nx, int4 ny, int4 nz);
#endif

#ifndef _EMFIT
extern void EMFIT1d(float8 *hist, float8 *fit, int2 *label, int4 nb, float8 *mean, float8 *var, float8 *p, int4 nclass, int4 niter);
#endif

#ifndef _max_cc

extern void Connected_Components( char *im, int4 nx, int4 ny, int4 nz, int4 *LABEL, int4 *N, int4 **Clabel, int4 **Size);
extern void max_Connected_Component(char *im, int4 nx, int4 ny, int4 nz,int4 *ncc, int4 *maxsize);
extern void thr_Connected_Component(char *im,int4 thr, int4 nx, int4 ny, int4 nz,int4 *ncc, int4 *ntcc);
extern void heightthr_Connected_Component(float4 *im, float4 thr, int4 nx, int4 ny, int4 nz, int4 *ncc, int4 *ntcc);
extern void Connected_Component_location(char *im, int4 nx, int4 ny, int4 nz, int4 *ncc, int4 **S, int4 **L);
extern void Connected_Component_location(char *im, int4 nx, int4 ny, int4 nz, int4 *ncc, int4 **S, int4 **L, int4 **CCsize);

#endif

#ifndef _statistics
extern float4 imageMean(int2 *im, int2 *msk, int4 nv);
extern float4 imageMean(int2 *im, int4 nv);
extern float8 one_sample_t(float8 *x, int4 n);

extern float4 median(float4 *x, char *mask, int4 n);
extern void scaleAbsToOne(float8 *y, int4 n, int4 p);
extern void scaleAbsToOne(float4 *y, int4 n, int4 p);
extern float8 scaleAbsToOne(float8 *y, int4 n);
extern void decomposeVector(float8 *x, float8 *xpar, float8 *xper, float8 *u, int4 n);

extern void removeComponent(float8 *x, float8 *y, int4 n);
extern void removeComponent(float4 *x, float4 *y, int4 n);

extern float8 removeVectorMean(float4 *y, int4 n);
extern float8 removeVectorMean(float8 *y, int4 n);
extern float8 removeVectorMean(float4 *x, float4 *y, int4 n);
extern float8 removeVectorMean(float8 *x, float8 *y, int4 n);
extern float8 removeVectorMean(int2 *x, float8 *y, int4 n);
extern void removeVectorMean(float8 *y, int4 n, int4 p);
extern void removeVectorMean(float4 *y, int4 n, int4 p);

extern void partialCorrelation(float8 *Y, float8 *X1, float8 *X2, int4 n, float8 *pr1, float8 *pr2);
extern void partialCorrelation(float4 *Y, float4 *X1, float4 *X2, int4 n, float8 *pr1, float8 *pr2);
extern void partialCorrelation(int2 *Y, float4 *X1, float4 *X2, int4 n, float8 *pr1, float8 *pr2);

#endif

#ifndef _ginverse
extern int4 ginverse(float4 *X, int4 N, int4 p, float4 *G);
extern int4 ginverse(float8 *X, int4 N, int4 p, float4 *G);
#endif

#ifndef _convolution
extern float4 conv_pnt_sk(uchar *x,int4 sx,float4 *h,int4 sh,int4 i0);
extern float4 conv_pnt_sk(int2 *x,int4 sx,float4 *h,int4 sh,int4 i0);
extern float4 conv_pnt_sk(float4 *x,int4 sx,float4 *h,int4 sh,int4 i0);
extern float4 *conv_sk(int2 *x,int4 sx,float4 *h,int4 sh);
extern float4 *conv_sk(float4 *x,int4 sx,float4 *h,int4 sh);
extern void conv_sk(float4 *x, float4 *y, int4 sx,float4 *h,int4 sh);
extern void conv_sk_inplace(float4 *x, int4 sx, float4 *h, int4 sh);
#endif

#ifndef _gaussian_kernel
extern float4 *gaussian_kernel(float4 sd, int4 *n);
#endif

#ifndef _hpsort
extern void hpsort(unsigned long n, float4 *ra);
extern void hpsort(unsigned long n, float4 *ra, int4 *indx);
#endif

#ifndef _random
extern void initializeRandomNumberGenerator();
extern void sampleWithReplacementIndex(int4 *I, int4 N);
#endif

#ifndef _dicomIO

#define IMPLICIT_LITTLE_ENDIAN 0
#define EXPLICIT_LITTLE_ENDIAN 1
#define EXPLICIT_BIG_ENDIAN 2

extern int4 readFileMetaInfo(const char *filename, DICOM_file_meta_info *file_meta_info, char v);
extern int4 readDataSet(const char *filename, DICOM_hdr *hdr, char v);

extern void readDicomInfo(const char *file, int4 np, dicominfo *info);
extern int4 readImageParams(const char *file, float4 *TLHC, float4 *rowvec, float4 *colvec, 
float4 *dx, float4 *dy, float4 *dz, int4 *TE, char *patientID, int4 *imageNumber, int4 *seriesNumber, int4 np);
extern int4 readPhaseEncodingDirection(const char *file, char *PED, int4 np);
extern int4 readPixelData(const char *file, char *data, int4 opt_j, int4 np);
extern int4 readTLHC(const char *file, float4 *TLHC, int4 np);
extern int4 readRowCol(const char *file, float4 *rowvec, float4 *colvec, int4 np);
extern int4 readNEX(const char *file, int4 *nex);
extern int4 readSliceThickness(const char *file, float4 *sliceThickness);
extern int4 readSeriesNumber(const char *file, int4 *seriesNumber, int4 np);
extern int4 readImageNumber(const char *file, int4 *imageNumber, int4 np);
extern int4 readTE(const char *file, int4 *TE, int4 np);
extern int4 readTR(const char *file, int4 *TR);
extern int4 readImageVoxelSize(const char *file, float4 *dx, float4 *dy, float4 *dz, int4 np);
extern int4 readImageMatrixSize(const char *file, uint2 *nx, uint2 *ny);
extern int4 readPatientID(const char *file, char *patientID, int4 np);
extern int4 dicomFormat(const char *file);
extern int4 readMetaFileInfo(const char *file, long *offset, int4 *syntax);
extern int4 readTag(const char *file, uint2 iGN, uint2 iEN, long byteOffset, int4 transferSyntax,
unsigned long *oVL, char *oV, long *valueOffset);
extern int4 read_element(char *filename, int2 S_GN, int2 S_EN, char *V, int4 *VL);
extern void readMatrixSize(char *filename, int4 *nx, int4 *ny);
extern void readVoxelSize(char *filename, float4 *dx, float4 *dy, float4 *dz);
extern int4 readImageSliceThickness(const char *file, float4 *dz, int4 np);

#endif

#ifndef _subsets
extern int4 *subsets(int4 N, int4 M);
extern int4 binomialCoeff(int4 N,int4 M);
#endif

#ifndef _cubicspline
extern float4 cubicSplineSynthesis(float4 *c, int4 nx, int4 ny, int4 nz, float4 x, float4 y, float4 z, float4 *beta, float4 del);

extern int2 *resliceImageCubicSpline(int2 *im1, int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1,
int4 nx2, int4 ny2, int4 nz2, float4 dx2, float4 dy2, float4 dz2, float4 *T);

extern float4 *resliceImageCubicSpline(float4 *im1, int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1,
int4 nx2, int4 ny2, int4 nz2, float4 dx2, float4 dy2, float4 dz2, float4 *T);

extern float4 *computeBeta(float4 *del);

extern void cubicSplineAnalysis(uchar *s, float4 *c, int4 nx, int4 ny, int4 nz);
extern void cubicSplineAnalysis(int2 *s, float4 *c, int4 nx, int4 ny, int4 nz);
extern void cubicSplineAnalysis(float4 *s, float4 *c, int4 nx, int4 ny, int4 nz);

extern void cubicSplineAnalysis(float4 *s, float4 *c, int4 N);

#endif

#ifndef _registration
extern void label_3d_cc(int2 *KMI,uint2 label,int4 i,int4 j, int4 k,int4 *size,int2 CC, struct im_params *IP);
extern void set_transformation(float4 x, float4 y, float4 z, float4 ax, float4 ay, float4 az, const char *code, float4 *T);
extern int4 label_CCI(int2 *KMI, int4 size_thresh,struct im_params * IP, int4 nvoxels);
extern float (*interpolator)(float4 x, float4 y, float4 z, int2 *array, int4 nx, int4 ny, int4 nz, int4 np);
extern float4 P[12];
extern struct im_params IP;

extern uchar linearInterpolatorUC(float4 x, float4 y, float4 z, uchar *array, int4 nx, int4 ny, int4 nz, int4 np);
extern void scale_short_minmax(int2 *imagein, uchar **imageout, int4 np, int4 min,int4 max);

extern void testCostFunc1(int2 *trg, int4 Tnx, int4 Tny, int4 Tnz, float4 Tdx, float4 Tdy, float4 Tdz,
int2 *obj, int4 Onx, int4 Ony, int4 Onz, float4 Odx, float4 Ody, float4 Odz);

extern float4 *transformation(float4 x, float4 y, float4 z, float4 ax, float4 ay,
float4 az, float4 sx, float4 sy, float4 sz, int4 rX, int4 rY, int4 rZ, char *code);

extern void cca(int2 *im, int4 nx, int4 ny, int4 nz);

extern int4 loadTransformation( char *filename, float4 *T);
extern int4 findThresholdLevel(int2 *image_in, int4 nv);
extern float4 *findTransMatrix(int2 *trg, int4 Tnx, int4 Tny, int4 Tnz, float4 Tdx, float4 Tdy, float4 Tdz,
int2 *obj, int4 Onx, int4 Ony, int4 Onz, float4 Odx, float4 Ody, float4 Odz);

extern int2 *KMcluster(int2 *ccImage, int2 *im_in, int4 nclass, int4 maxiter, int4 thresh, int4 Low, int4 High, int4 nv);
#endif

#ifndef _legendre
//extern float4 *FourierLegendreSynthesis(float8 *c, int4 nx, int4 ny, int4 nz, int4 mx, int4 my, int4 mz);
//extern float8 *FourierLegendreAnalysis(float4 *f, int4 nx, int4 ny, int4 nz, int4 mx, int4 my, int4 mz, int4 N);
extern float8 *LegendreAnalysis(float4 *image, int4 nx, int4 ny, int4 nz, int4 mx, int4 my, int4 mz);
extern	float8 *LegendreAnalysis(float4 *image, int4 nx, int4 ny, int4 mx, int4 my);
extern void LegendreSynthesis(float8 *c, int4 mx, int4 my, float4 *image, int4 nx, int4 ny);
extern float4 *LegendreSynthesis(float8 *c, int4 nx, int4 ny, int4 nz, int4 mx, int4 my, int4 mz);
void LegendrePoly(float8 *p0, float8 *q0, float8 *p1, float8 *q1, float8 *x, int4 N, int4 n);
void LegendrePoly(float4 *p0, float4 *q0, float4 *p1, float4 *q1, float4 *x, int4 N, int4 n);
void integral_1d(float8 *a, int2 *b, int4 n, float8 *d);
void integral_1d(float8 *a, float4 *b, int4 n, float8 *d);
#endif

#ifndef _matrixCom

extern float8 Frobenius_s3(float8 *qqT, float4 *D);
extern float8 Frobenius_s3(float8 *qqT, float8 *D);
extern void p3update(float8 *D, float8 *W, float8 epsilon);
extern void p3update(float8 *D, float8 *invsqrtD, float8 *sqrtD, float8 *W, float8 epsilon);
extern void s3multi(float8 *A, float8 *B, float8 *AB);
extern void p3invsqrt(float8 *A, float8 *invsqrtA, float8 *sqrtA);
extern void p3invsqrt(float8 *A, float8 *invsqrtA);
extern int4 p3invsqrt(float4 *A, float4 *invsqrtA);
extern float8 p3RiemannianDistance(float8 *D, float8 *invsqrtF);
extern float8 p3RiemannianDistance(float8 *L);
extern float4 p3RiemannianDistance(float4 *D, float4 *invsqrtF);
extern void s3ABA(float8 *A, float8 *B, float8 *ABA);
extern void s3ULUT(float8 *L, float8 *UT, float8 *ULUT);
extern void s3eigenvec(float8 *A, float8 *evalue, float8 *evector);
extern void getcomplement1(float8 *a, float8 *b, float8 *c);
extern void getcomplement2 (float8 *U, float8 *V, float8 *W);
extern void s3inv(float8 *A, float8 *invA);
extern int4 s3inv(float4 *A, float4 *invA);
extern void s3multi(float8 *A, float8 *B, float8 *AB);
extern void differentiate_distance(float8 *D, float8 *F, float8 *L, float8 *dLdD);

extern void subtractRowAvg(float4 *X, int4 N, int4 P, float4 *X0);

extern void crossProduct(float8 *a, float8 *b, float8 *c);
extern void crossProduct(float4 *a, float4 *b, float4 *c);

extern void copyVector(float4 *v1, float4 *v2, int4 n);
extern void subtractVector(float4 *v1, float4 *v2, int4 n);
void normalizeVector(float4 *x, int4 n, float8 *norm);

extern int4 centerMatrixRow(float4 *X, int4 N, int4 P, float4 *avg);
extern int4 centerMatrixRow(float4 *X, int4 N, int4 P);

// Returns 1 if an error condition occurs, 0 otherwise

// returns an N-vector in "avg" containing the average of the P columns of X
int4 avgCol(float4 *X, int4 N, int4 P, float4 *avg);
int subtractAvgCol(float *X, int N, int P, float *avg);

extern int4 avgRow(float4 *X, int4 N, int4 P, float4 *avg, char *rowmask);
extern int4 avgRow(float8 *X, int4 N, int4 P, float8 *avg, char *rowmask);
extern int4 avgRow(float4 *X, int4 N, int4 P, float4 *avg);
extern int4 avgRow(float8 *X, int4 N, int4 P, float8 *avg);
extern int4 varRow(float4 *X, int4 N, int4 P, float4 *avg, float4 *var);
extern int4 varRow(float8 *X, int4 N, int4 P, float8 *avg, float8 *var);
extern int4 varRow(float8 *X, int4 N, int4 P, float8 *avg, float8 *var, char *rowmask);
extern int4 varRow(float4 *X, int4 N, int4 P, float4 *avg, float4 *var, char *rowmask);
extern int4 ssdRow(float4 *X, int4 N, int4 P, float4 *avg, float4 *ssd, char *rowmask);
extern int4 ssdRow(float4 *X, int4 N, int4 P, float4 *avg, float4 *ssd);

// compute the Euclidian distance between two vectors r0 and r1
extern float8 euclideandistance(float4 *r0, float4 *r1, int4 n);
extern float8 xtAx(float4 *A, float8 *x, int4 p);
extern float8 vectorNorm(float4 *x, int4 n);
extern void normalizeVector(float4 *x, int4 n);
extern void transpose_matrix(float4 *A, int4 N,  int4 M);
extern void transpose_matrix(float4 *A, int4 N,  int4 M, float4 *AT);
/* computes the transpose of the NxM matrix A */
extern float *trans(float *A, int N,  int M);
extern double *trans(double *A, int N,  int M);

extern float4 normalize(float4 *s, int4 n);
extern float8 normalize(float8 *s, int4 n);

extern int4 ComputeRank(float8 *M);
extern void s3eigenval(float8 *A, float8 *L);
extern float8 s3tr(float8 *A, float8 *B);
extern float4 s3tr(float4 *A, float4 *B);
extern void s3vec_to_mat(float4 *M, float4 *V);
extern void s3vec_to_mat(float8 *M, float8 *V);
extern void s3mat_to_vec(float4 *M, float4 *V);
extern void s3mat_to_vec(float4 *M, float8 *V);
extern void s3mat_to_vec(float8 *M, float8 *V);
extern void s3adjoint(float8 *A, float8 *ADJ);
extern void s3adjoint(float4 *A, float4 *ADJ);
extern float8 s3det(float8 *A);
extern float4 s3det(float4 *A);
void ds3det(float8 *A, float8 *B);
void ds3det(float4 *A, float4 *B);
extern float8 det3(float8 *A);
extern float4 det3(float4 *A);
extern float4 det4(float4 *A);
extern float8 det4(float8 *A);
extern float8 *inv3(float8 *A);
extern float4 *inv3(float4 *A);
extern void inv3(float4 *A, float4 *invA);
extern float4 *inv2(float4 *A);
extern float8 *inv2(float8 *A);
extern float4 *inv4(float4 *A);
extern float8 *inv4(float8 *A);
extern void multi(float4 *A,int4 iA,int4 jA,float4 *B,int4 iB,int4 jB,float4 *C);
extern void multi(float8 *A,int4 iA,int4 jA,float8 *B,int4 iB,int4 jB,float8 *C);
extern void multi(float4 *A,int4 iA,int4 jA, float8 *B,int4 iB,int4 jB,float8 *C);
extern void multi(float8 *A,int4 iA,int4 jA, float4 *B,int4 iB,int4 jB,float4 *C);
#endif

//////////////////////////////////////////////////////////////////////////////
// The following are define in reslice.c
#define LIN 1
#define NEARN 2
#define SINC 3
#define CUBICSPLINE 4	

void resliceImage(SHORTIM im1, SHORTIM &im2, float4 *T, int4 interpolation_method);

int2 *resliceImage(int2 *im1, DIM dim1, DIM dim2, float4 *T, int4 interpolation_method);
float4 *resliceImage(float4 *im1, DIM dim1, DIM dim2, float4 *T, int4 interpolation_method);

float4 *resliceImage(float4 *im1, int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1,
int4 nx2, int4 ny2, int4 nz2, float4 dx2, float4 dy2, float4 dz2, float4 *T);

float4 *resliceImage(float4 *im1, int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1,
int4 nx2, int4 ny2, int4 nz2, float4 dx2, float4 dy2, float4 dz2, float4 *T, float4 *xjit, float4 *yjit);

float4 partial_var(float4 x, float4 y, float4 z, uchar *array, int4 nx, int4 ny, int4 nz, int4 np, float4 mu);

// you must initialize drand48 before using this function
uchar PNN(float4 x, float4 y, float4 z, uchar *array, int4 nx, int4 ny, int4 nz);

char *resliceImage(char *obj, int4 Onx, int4 Ony, float4 Odx, float4 Ody, int4 Tnx, int4 Tny, float4 Tdx, float4 Tdy, float4 *T);
int2 *resliceImage(int2 *im1, int4 nx1, int4 ny1, float4 dx1, float4 dy1, int4 nx2, int4 ny2, float4 dx2, float4 dy2, float4 *T);
int2 *resliceImage(float4 *im1, int4 nx1, int4 ny1, float4 dx1, float4 dy1, int4 nx2, int4 ny2, float4 dx2, float4 dy2, float4 *T);

int2 *resliceImage(int2 *im1, int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1,
int4 nx2, int4 ny2, int4 nz2, float4 dx2, float4 dy2, float4 dz2, float4 *T, int4 interpolation_method);

uchar *resliceImage(uchar *im1, int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1,
int4 nx2, int4 ny2, int4 nz2, float4 dx2, float4 dy2, float4 dz2, float4 *T);

float4 *resliceImage(float4 *im1, int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1,
int4 nx2, int4 ny2, int4 nz2, float4 dx2, float4 dy2, float4 dz2, float4 *T, float4 *w);

int2 *resliceImage(int2 *im1, int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1,
int4 nx2, int4 ny2, int4 nz2, float4 dx2, float4 dy2, float4 dz2, float4 *Xwarp, float4 *Ywarp, float4 *Zwarp);

float4 linearInterpolator(float4 x, float4 y, float4 z, float4 *array, int4 nx, int4 ny, int4 nz, int4 np, float4 *w);
uchar  linearInterpolator(float4 x, float4 y, float4 z, uchar  *array, int4 nx, int4 ny, int4 nz, int4 np, float4 *w);

int2 *computeReslicedImage(int2 *im1, int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1,
int4 nx2, int4 ny2, int4 nz2, float4 dx2, float4 dy2, float4 dz2, float4 *Xwarp, float4 *Ywarp, float4 *Zwarp);

float4 *computeReslicedImage(float4 *im1, int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1,
int4 nx2, int4 ny2, int4 nz2, float4 dx2, float4 dy2, float4 dz2, float4 *Xwarp, float4 *Ywarp, float4 *Zwarp);

int2 *computeReslicedImage(float4 *im1, int4 nx1, int4 ny1, float4 dx1, float4 dy1,
int4 nx2, int4 ny2, float4 dx2, float4 dy2, float4 *Xwarp, float4 *Ywarp);

int2 *computeReslicedImage(int2 *im1, int4 nx1, int4 ny1, float4 dx1, float4 dy1,
int4 nx2, int4 ny2, float4 dx2, float4 dy2, float4 *Xwarp, float4 *Ywarp);
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// defined in resize.c
int2 *resizeXYZ(int2 *image1,  DIM dim1, DIM dim2);

int2 *resizeXYZ(int2 *image1,
int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1,
int4 nx2, int4 ny2, int4 nz2, float4 dx2, float4 dy2, float4 dz2);

uchar *resizeXYZ(uchar *image1,
int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1,
int4 nx2, int4 ny2, int4 nz2, float4 dx2, float4 dy2, float4 dz2);

int2 *resizeXYZ(char *image1,
int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1,
int4 nx2, int4 ny2, int4 nz2, float4 dx2, float4 dy2, float4 dz2);

float4 *resizeXYZ(float4 *image1,
int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1,
int4 nx2, int4 ny2, int4 nz2, float4 dx2, float4 dy2, float4 dz2);

float4 *resizeXY(float4 *image1, int4 nx1, int4 ny1, float4 dx1, float4 dy1, int4 nx2, int4 ny2, float4 dx2, float4 dy2);
int2 *resizeXY(int2 *image1, int4 nx1, int4 ny1, float4 dx1, float4 dy1, int4 nx2, int4 ny2, float4 dx2, float4 dy2);
//////////////////////////////////////////////////////////////////////////////

#ifndef _getoption

extern int4 optind;
extern char *optarg;

struct option
{
	const char *name;
	int4 has_arg;
	int4 val;
};

extern int4 getoption(int4 argc, char **argv, struct option *options);
#endif

////////////////////////////////////////////////////////////////////////////////////////
// The following functions are defined in analyzeio.c
int4 extension_is_hdr(const char *filename);
void read_analyze_image(const char *filename, int2 *im);
char *read_analyze_image(const char *filename, DIM *dim, int4 *type, int4 v);
float4 read_dx(const char *hdrfile);
float4 read_dy(const char *hdrfile);
float4 read_dz(const char *hdrfile);
int4 read_nt(const char *hdrfile);
int4 read_nx(const char *hdrfile);
int4 read_ny(const char *hdrfile);
int4 read_nz(const char *hdrfile);
int4 read_datatype(char *hdrfile);
char *read_analyze_image(const char *filename, int4 *nx, int4 *ny, int4 *nz, float4 *dx, float4 *dy, float4 *dz);
char *read_analyze_image(const char *filename, int4 *nx, int4 *ny, int4 *nz, int4 *nt, float4 *dx, float4 *dy, float4 *dz, int4 *type, int4 v);
char *read_analyze_image(const char *filename, int4 *nx, int4 *ny, int4 *nz, float4 *dx, float4 *dy, float4 *dz, int4 *type, int4 v);
char *read_analyze_image(const char *filename, int4 *nx, int4 *ny, int4 *nz, float4 *dx, float4 *dy, float4 *dz, int4 *type);
char *read_image(char *file,int4 n);
void get_analyze_file_names(const char *filename, char *basename_hdr, char *basename_img);
void read_analyze_hdr(struct dsr *hdr, char *filename);
void setDimensions(struct dsr hdr, int4 *nx, int4 *ny, int4 *nz, float8 *dx, float8 *dy, float8 *dz, int2 *dataType);
void setDimensions(struct dsr hdr, int4 *nx, int4 *ny, int4 *nz, float4 *dx, float4 *dy, float4 *dz, int2 *dataType);
void setDimensions(struct dsr hdr, int4 *nx, int4 *ny, int4 *nz, float4 *dx, float4 *dy, float4 *dz, int2 *dataType, int4 v);
void setDimensions(struct dsr hdr, int4 *nx, int4 *ny, int4 *nz, int4 *nt, float4 *dx, float4 *dy, float4 *dz, int2 *dataType, int4 v);
void create_analyze_hdr(struct dsr *hdr, int4 nx, int4 ny, int4 nz, int4 dt, float4 dx, float4 dy, float4 dz);
void create_analyze_hdr(struct dsr *hdr, int4 nx, int4 ny, int4 nz, int4 nt, int4 datatype, float4 dx, float4 dy, float4 dz);
void write_analyze_image(const char *filename, int2 *im, int4 nx, int4 ny, int4 nz, float4 dx, float4 dy, float4 dz); 
void write_analyze_image(const char *filename, float4 *im, int4 nx, int4 ny, int4 nz, float4 dx, float4 dy, float4 dz); 
void write_analyze_image(const char *filename, uchar *im, int4 nx, int4 ny, int4 nz, float4 dx, float4 dy, float4 dz); 
void write_analyze_image(const char *filename, uchar *im, int4 nx, int4 ny, int4 nz, float4 dx, float4 dy, float4 dz, int4 v); 
void write_analyze_image(const char *filename, int2 *im, int4 nx, int4 ny, int4 nz, float4 dx, float4 dy, float4 dz,int4 v); 
void write_analyze_image(const char *filename, float4 *im, int4 nx, int4 ny, int4 nz, float4 dx, float4 dy, float4 dz,int4 v); 
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
// functions defined in swap.cxx
extern int4 bigEndian();
extern void swapByteOrder(char *in, int4 N);
extern void swapN(char *in, int4 N);
extern void swap_float_array(float4 *x, int4 n);
extern void swap_double_array(float8 *x, int4 n);
extern void swap_int_array(int4 *x, int4 n);
////////////////////////////////////////////////////////////////////////////////////////

#ifndef _medianfilter
extern void medianFilter(float4 *image1, int4 nx, int4 ny, int4 nz, int4 Wx, int4 Wy, int4 Wz);
#endif

#ifndef _fileinfo
extern int4 isregular(const char *file);
extern int4 getFileSize(const char *file);
extern int4 checkWriteAccess(char *file);
extern int4 checkReadAccess(char *file);
extern int4 checkFileExistence(char *file);
extern int4 checkFileReadOK(char *file);
extern int4 checkFileWriteOK(char *file);
extern int4 check_F_R_permission(char *file);
#endif

#ifndef _histogram
extern int4 otsu(float8 *histogram, int4 numberOfBins);
extern float8 *findHistogram(int2 *im, int4 nv, int4 *nb, int2 &min, int2 &max);
extern float8 *findHistogram(int2 *im, int4 nv, int4 nb, int4 low, int4 high, int4 *bw_return);
extern float8 *findHistogram(int2 *im1, int2 *im2, int4 nv, int4 nb1, int4 nb2, int4 *bw1_r, int4 *bw2_r, int4 low1, int4 high1, int4 low2, int4 high2);
extern void trimExtremes(int2 *image, int2 *msk, int4 nv, float4 percent);
#endif

#ifndef _matrixops
extern void zeroVector(float4 *v, int4 n);
extern void zeroVector(char *v, int4 n);
extern void oneVector(float4 *v, int4 n);
extern void oneVector(char *v, int4 n);
extern void svd(float4 *At, int4 M, int4 N, float4 *Ut, float4 *V, float4 *S);
extern int4 zeroRowCol(float4 *A, int4 N, int4 n);
extern int4 setRowCol(float4 *A, int4 N, int4 n, float4 *a);
extern void projectVector(float8 *x, float8 *xpar, float8 *xper, float4 *Pz, int4 n);
extern float4 *projectionMatrix(float8 *X, int4 N, int4 p, int4 *rank);
extern float4 *projectionMatrix(float4 *X, int4 N, int4 p);
extern void mat_mat_trans(float4 *A,int4 Ar,int4 Ac,float4 *B,int4 Br, float4 *C);
extern float4 *diagATA_float(float4 *ATA, int4 n, char uplo);
extern float4 *AAT_float(float4 *A,int4 nr,int4 nc, char uplo);
extern void mat_trans_mat(float4 *A, int4 Ar, int4 Ac, float4 *B, int4 Bc, float4 *C);
extern void mat_trans_mat(float8 *A, int4 Ar, int4 Ac, float8 *B, int4 Bc, float8 *C);
#endif

///////////////////////////////////////////////////////////////
// The following functions are defined in nifti.cxx
int4 not_magical_nifti(const char *imagefilename);
int4 not_magical_nifti(const char *imagefilename, int verbose);
char *read_nifti_image(const char *filename, nifti_1_header *hdr);
int4 same_nifti_image_size(int4 N, char **imagefile, int4 *nx, int4 *ny, int4 *nz, float4 *dx, float4 *dy, float4 *dz);
void read_nifti_image(const char *filename, uchar **im, nifti_1_header *hdr);
void read_nifti_image(const char *filename, int2 **im, nifti_1_header *hdr);
void print_NIFTI_hdr(const char *filename);
void print_NIFTI_hdr(nifti_1_header hdr);
nifti_1_header read_NIFTI_hdr(const char *filename);
int4 read_NIFTI_hdr(const char *filename, nifti_1_header *hdr);
nifti_1_header read_NIFTI_hdr(const char *filename, nifti1_extender *extender, char **extension);
void save_nifti_image(const char *filename, uchar *im, nifti_1_header *hdr);
void save_nifti_image(const char *filename, int2 *im, nifti_1_header *hdr);
void save_nifti_image(const char *filename, float4 *im, nifti_1_header *hdr);
void save_nifti_image(const char *filename, char *im, nifti_1_header *hdr);

// returns the orientations vectors xvec, yvec, and zvec in NIFTI's RAS system
void readOrientationVectorsFromFile(const char *filename, float4 *xvec, float4 *yvec, float4 *zvec);

int4  niftiFilename(char *filename, const char *path);
void swapniftiheader(nifti_1_header *hdr);
int2 *readNiftiImage(const char *filename, DIM *dim, int4 flg);
///////////////////////////////////////////////////////////////


#ifndef _utils

extern int hand_system(char *code);

extern void art_to_fsl(float *Mart, float *Mfsl, DIM sub_dim, DIM trg_dim);
extern void fsl_to_art(float *Mfsl, float *Mart, DIM sub_dim, DIM trg_dim);
extern void fsl_to_art(float *Mfsl, float *Mart, DIM sub_dim, DIM trg_dim, int subflg, int trgflg);

// Set the nxn matrix A equal to the identity matrix
extern void set_to_I( float4 *A, int4 n);

extern void sobel_edge_x(int2 *in, float4 *out, int4 nx, int4 ny);
extern void sobel_edge_y(int2 *in, float4 *out, int4 nx, int4 ny);
extern void sobel_edge(int2 *in, float4 *out, int4 nx, int4 ny);
extern void sobel_edge(int2 *in, int2 *out, int4 nx, int4 ny);

extern int4 ccsize(int2 *im, int4 nv);

extern void copyarray(float4 *source, float4 *destination, int4 size);
extern void copyarray(int2 *source, char *destination, int4 size);
extern void copyarray(char *source, int2 *destination, int4 size);

extern void zeroarray(float4 *y, int4 size);
extern float4 diceindex(int2 *setA, int2 *setB, int4 n);
extern void remove_space(char *inp);
extern void orientationCodeConverter(int4 integerCode, char *stringCode);
extern int4 orientationCodeOK(char *stringCode);

extern void ijk2xyz(float4 *T, DIM dim);
extern void ijk2xyz(float4 *T, int4 nx, int4 ny, int4 nz, float4 dx, float4 dy, float4 dz);
extern void xyz2ijk(float4 *T, int4 nx, int4 ny, int4 nz, float4 dx, float4 dy, float4 dz);
extern void xyz2ijk(float4 *T, DIM dim);
extern void saveMatrix(float4 *A, int4 n, int4 m, char *filename);
extern void saveMatrix(int2 *A, int4 n, int4 m, char *filename);
extern float4 *readMatrix(int4 *n, int4 *m, char *filename);

extern float4 *readDataMatrix(char **imageList, int4 n, int4 p, int2 *mask);
extern float4 *readDataMatrix_nifti(char **imageList, int4 n, int4 p, int2 *mask);

extern int2 *readDataMatrixShort(char **imageList, int4 n, int4 p, int2 *mask);
extern int2 *readDataMatrixShort_nifti(char **imageList, int4 n, int4 p, int2 *mask);

extern int2 *readMask(const char *filename, int4 *nx, int4 *ny, int4 *nz);
extern int2 *readMask_nifti(const char *filename, int4 *nx, int4 *ny, int4 *nz);


extern void checkDimension(int4 N, char **imagefile, int4 nx, int4 ny, int4 nz);
extern void checkDimension_nifti(int4 N, char **imagefile, int4 nx, int4 ny, int4 nz);
extern void checkDimension(int4 N, char **imagefile, int4 *nx, int4 *ny, int4 *nz, float4 *dx, float4 *dy, float4 *dz);


extern void affineLSE(int2 *msk, int4 nx, int4 ny, int4 nz, float4 dx, float4 dy, float4 dz, float4 *Xwarp, float4 *Ywarp, float4 *Zwarp, float4 *T);
extern void affineLSE(char *msk, int4 nx, int4 ny, int4 nz, float4 dx, float4 dy, float4 dz, float4 *Xwarp, float4 *Ywarp, float4 *Zwarp, float4 *T);
extern float4 *affineLSE(char *msk, int4 nx, int4 ny, float4 dx, float4 dy, float4 *Xwarp, float4 *Ywarp);
extern void affineLSE(char *msk, int4 nx, int4 ny, float4 dx, float4 dy, float4 *Xwarp, float4 *Ywarp, float4 *T);

extern void affineLSE(int2 *msk, int4 nx, int4 ny, float4 dx, float4 dy, float4 *Xwarp, float4 *Ywarp, float4 *T);

extern void extractArray(float4 *im, int4 nx, int4 ny, int4 nz, int4 np, int4 nx0, int4 ny0, int4 nz0, int4 Lx, int4 Ly, int4 Lz, float4 *array);
extern void extractArray(uchar *im, int4 nx, int4 ny, int4 nz, int4 nx0, int4 ny0, int4 nz0, int4 Lx, int4 Ly, int4 Lz, float4 *array);
extern void extractArray(float4 *im, int4 nx, int4 ny, int4 nz, int4 nx0, int4 ny0, int4 nz0, int4 Lx, int4 Ly, int4 Lz, float4 *array);
extern void extractArray(int2 *im, int4 nx, int4 i,  int4 j,  int4 L, float4 *array);
extern void extractArray(int2 *im, int4 nx, int4 ny, int4 nz, int4 x0, int4 y0, int4 z0, int2 *x,int2 *y,int2 *z,int4 n,float4 *array);
extern void extractArray(int2 *im, int4 nx, int4 ny, int4 nx0,int4 ny0,int4 Lx, int4 Ly, float4 *array);
extern void extractArray(int2 *im, int4 nx, int4 ny, int4 nz, int4 nx0,int4 ny0, int4 nz0, int4 Lx, int4 Ly, int4 Lz, float4 *array);
extern int4 extractArray(int2 *im, int4 nx, int4 ny, int4 nz, int4 np, int4 nx0,int4 ny0, int4 nz0, int4 Lx, int4 Ly, int4 Lz, int2 *array);
extern void extractArray(int2 *im, int4 nx, int4 ny, int4 nz, int4 nx0,int4 ny0,int4 nz0, int4 Lx, int4 Ly, int4 Lz, int2 *array);
extern void extractArray(int2 *im, int4 nx, int4 ny, int4 nz, int4 np, int4 nx0,int4 ny0, int4 nz0, int4 Lx, int4 Ly, int4 Lz, float4 *array);

// This function extracts the "filename" from the full "path" string.
// For example, if path="/home/babak/testimages/test.img", then filname="test.img".
extern void getfilename(char *filename, const char *path);

extern void printMatrix(int4 *mat, int4 n, int4 p, const char *s, FILE *fp);
extern void printMatrix(float4 *mat, int4 n, int4 p, const char *s, FILE *fp);
extern void printMatrix(float8 *mat, int4 n, int4 p, const char *s, FILE *fp);
extern void get_temp_filename(char *filename);
extern void mask_and_save(const char *inputfile, const char *outputfile, int2 *mask, int2 *masked_image, int4 nbv, float4 FWHM);
extern void mask_and_save_nii(const char *inputfile, const char *outputfile, int2 *mask, int2 *masked_image, int4 nbv, float4 FWHM);
extern void read_transpose_save(char *inputfile, char *outputfile, int4 nr, int4 v);
extern void centerOfMass(int2 *im, int4 nx, int4 ny, int4 nz, float4 dx, float4 dy, float4 dz, float4 *CM);

#endif

////////////////////////////////////////////////////////////////////////
// defined in smooth.cxx
float4 *smoothY(float4 *image, int4 nx, int4 ny, int4 nz, float4 sd);
float4 *smoothZ(float4 *image, int4 nx, int4 ny, int4 nz, float4 sd);
////////////////////////////////////////////////////////////////////////


#endif
