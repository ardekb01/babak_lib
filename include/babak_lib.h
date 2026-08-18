///////////////////////////////////////////////////////////////////////
// Copyright (C) 2024 Babak A. Ardekani, PhD - All Rights Reserved.
///////////////////////////////////////////////////////////////////////

#ifndef BABAK_LIB_H
#define BABAK_LIB_H

#include <stdio.h>
#include <stdlib.h>

#include "nifti1.h"
#include "bbk.h"

#ifndef _PILTRANSFORM
extern int opt_CENTER_AC;
#endif

#ifndef GETARTHOME
extern const char *ARTHOME;
extern bool getARTHOME();
#endif

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned short ushort;
typedef unsigned int uint;

#define YES 1
#define NO 0
#define ESMALL 1e-10
#define DEFAULT_STRING_LENGTH 1024
#define MAXIM 256 // maximum number of images allowed 

#ifndef _TPS
extern float tpsU(float x, float y);
extern float *tpsK(float *x, float *y, int n);
#endif

#ifndef _update_qsform
void update_qsform(nifti_1_header &hdr, float *matrix);
void update_qsform( const char *imagefilename , float *matrix);
void update_qsform( const char *imagefile1, const char *imagefile2);
extern int opt_sform;
extern int opt_qform;
#endif

#ifndef _artlib
// Turns verbose mode on (opt_v=YES) or off (opt_v=NO)
// Initialized to NO in artlib.c
extern char opt_v;
extern double searchradius[3]; // in units of mm
#endif

struct DICOM_file_meta_info
{
   // maximum bytes for VR=UI is 64,  1 byte is added for the string terminator '\0'
   char media_storage_SOP_class[65];  // Tag (0002,0002)
   char transfer_syntax[65]; // Tag (0002,0010)
   int dataset_offset;
};
typedef struct DICOM_file_meta_info DICOM_file_meta_info;

struct DICOM_hdr
{
   char patientID[65];    // Tag (0010,0020)  VR=LO (Long String 64 chars max.)
   float slice_thickness; // Tag (0018,0050)
   int TR;             // Tag (018,0080)
   int TE;             // Tag (018,0081)
   int TI;             // Tag (018,0082)
   int NEX;            // Tag (018,0083)
   float dz;           // Tag (018,0088)
   float flip_angle;   // Tag (018,1314)
   char seriesID[65];  // Tag (020,000E)
   int series_number;  // Tag (020,0011)
   float TLHC[3];      // Tag (020,0032)
   float rowvec[3];    // Tag (020,0037)
   float colvec[3];    // Tag (020,0037)
   char frame_of_referenceID[65]; // Tag (020,0052)
   int nz; // Tag (020,1002)
   float slice_location; // Tag (020,1041)
   int ny; // Tag (028,0010)
   int nx; // Tag (028,0011)
   float dy; // Tag (028,0030)
   float dx; // Tag (028,0030)
};
typedef struct DICOM_hdr DICOM_hdr;

//////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////
struct im_params {

   int MRimageId;
   int PETimageId;

	/* maximum values of the image 1 */	
	uchar max1; 	

	int 	*Size;   /* size of non-zero CC's in CCI */
	float *S,*SS;  

	int 	FROM;
	int 	TO;

	/* number of image 1 columns, rows and slices  */
	int	nx1,ny1,nz1; 

	/* number of image 1 voxels, number of image 1 pixels per slice */
	int	nv1,np1;      	

	float  dx1,dy1,dz1; 	/* image 1 voxel size (mm) */

	/* coordinates of image 1 volume center wrt pixel (0,0,0) in (mm) */
	float  xc1,yc1,zc1; 	

	// number of image 2 columns, rows and slices  
	int	nx2,ny2,nz2; 
	float  dx2,dy2,dz2; 	// image 2 voxel size (mm) 

	/* number of image 2 voxels, number of image 2 pixels per slice */
	int	nv2,np2;      

	/* coordinates of image 2 volume center wrt pixel (0,0,0) in (mm) */
	float  xc2,yc2,zc2; 	

	int NCC;

	int NP;  // number of non-zero pixels in CCI

	int sf;  // speed factor 

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
   ushort acquisitionMatrix[4];
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

// Input: (x,y,z) a vector defined in RAS system
// Output: One of six charaters {R,L,A,P,S,I}
void new_PIL_transform(const char *subfile, const char *lmfile, char *orient, float *T);
void standard_PIL_transformation(const char *imfile, const char *lmfile, char *orient, float *TPIL);
void Procrustes(float *Q, float *Qavg, int n, float *P, float *Pavg, float *TLM);
void Procrustes(float *Q, int n, float *P, float *TLM);
void compute_cm(short *image, DIM dim, float *cm);
void compute_cov(short *image, DIM dim, float *cm, double *I);
void standardize(float *x, int n);
void standardize(float *x, float *mask, int n);
void irodrigues_formula(float *R, float *w, float &theta);
void rodrigues_formula(float *R, float *w, float theta);
void se3_to_SE3(float *M, float *w, float *v, float theta);
void SE3_to_se3(float *M, float *w, float *v, float &theta);

#ifndef _singular_value_decomposition
extern int Singular_Value_Decomposition(double* A, int nrows, int ncols, double* U, 
                      double* singular_values, double* V, double* dummy_array);
extern void Singular_Value_Decomposition_Inverse(double* U, double* D, double* V,  
                        double tolerance, int nrows, int ncols, double *Astar);
#endif

#ifndef _DKI
extern void formA7(float *A, int n, float *u1, float *u2, float *u3, float *b, float *w);

extern void restricted_tensor_model(float *x, float *y);
extern void second_and_fourth_order_moments(float *x, float *y, float *E, float *F, int v);
extern float mardia_kurtosis(float *x1, float *x2, int v);
extern float mardia_kurtosis(float *x);

extern float estimate_K(float *A, float *x, float *y, int n);
extern void update_x(float *A, float *x, float *y, float K, int n);
extern void update_v(float *u0, float *u1, float *u2, float *b, float *y, int n, float *v, float scale);
extern void update_v(float *u0, float *u1, float *u2, float *b, float *y, int n, float *v1, float *v2, float *v3, float scale);
extern float quadratic_form_1(float *x, float *v);
extern float quadratic_form_1(float x0, float x1, float x2, float *v);
extern void vector_to_symmetric_matrix_form(float *u, float *D);
extern void symmetric_matrix_to_vector_form(float *D, float *u);
extern void cholesky_composition_1(float *v, float *u);
extern void cholesky_composition_2(float *L, float *D);
extern int cholesky_decomposition_1(float *u, float *v);
extern int cholesky_decomposition_2(float *D, float *L);
extern void positive_definite_tensor_model(float *v1, float *v2, float *y);
#endif

#ifndef _artlib
extern void find_pil_transformation(char *imfile, DIM dim, float *pilT, float *AC, float *PC, float *VSPS);
extern void find_pil_transformation(char *imfile, DIM dim, float *pilT);
extern void update_qsform(nifti_1_header &hdr, const char *orientationcode);
extern char opt_ppm;
extern char opt_png;
extern char opt_txt;
extern char opt_AC; // if YES AC will be detected automatically
extern char opt_PC; // if YES PC will be detected automatically
extern char opt_RP; // if YES RP will be detected automatically
extern char opt_MSP; // if YES MSP will be detected automatically
extern void sub2trg_rigid_body_transformation(float *sub2trg, const char *subfile, const char *trgfile);
extern void forwardTCSAP(float *xvec, float *yvec, float *zvec, float *TLHC, float *angle, float *translation, DIM dim);
extern void backwardTCSAP(float *xvec, float *yvec, float *angle);
extern void orig_ijk_to_pil_xyz(float *Tmsp, DIM orig_dim, float *AC, float *PC);
extern void initialAC(float Ax, float Ay, float Bx, float By, float *Cx, float *Cy, float parcomMean, float percomMean);
extern short *thresholdImageOtsu(short *im, int nv, int *nbv);
extern void defineTemplate(int r, int h, short *x, short *y, short *z);
extern char *defineACregion(DIM dim, float *RP, float *PC, float parcomMean, float percomMean, double ACsr);
extern char *definePCregion(DIM HR, float *RP, float *RPPCmean, double PCsr);
extern void ACPCtransform(float *Tacpc, float *Tmsp, float *AC, float *PC, char flg);
extern void compute_MSP_parameters_from_Tmsp(float *Tmsp, float *n, float *d);

bool detect_AC_PC_MSP(const char *imagefilename,
                     char *orientation,
                     float *AC, 
                     float *PC, 
                     float *RP, 
                     float *Tmsp,
                     int opt_T2);

extern float optimizeNormalVector(short *image,DIM dim,float &A, float &B, float &C);
extern float msp(short *im_in, int nx, int ny, int nz, float dx, float dy, float dz, float &A, float &B, float &C);
extern void computeTmsp(char *orientation, short *volOrig, DIM dim, float *Tmsp);
extern void combine_warps_and_trans(int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T);
#endif

#ifndef _binomial
extern double BinomialCoeff( int n, int k );
extern void k_subset ( int k, double rank, int a[] );
#endif

#ifndef _maskOps
extern char *find_foreground_mask(short *im, int nv, int nb, int nclass, int niter, short *thresh);
extern short *mask_image(short *im, char *mask, int nv, int *nbv_out);
extern void cc_filter(char *mask, int nx, int ny, int nz);
#endif

#ifndef _EMFIT
extern void EMFIT1d(double *hist, double *fit, short *label, int nb, double *mean, double *var, double *p, int nclass, int niter);
#endif

#ifndef _max_cc

extern void Connected_Components( char *im, int nx, int ny, int nz, int *LABEL, int *N, int **Clabel, int **Size);
extern void max_Connected_Component(char *im, int nx, int ny, int nz,int *ncc, int *maxsize);
extern void thr_Connected_Component(char *im,int thr, int nx, int ny, int nz,int *ncc, int *ntcc);
extern void heightthr_Connected_Component(float *im, float thr, int nx, int ny, int nz, int *ncc, int *ntcc);
extern void Connected_Component_location(char *im, int nx, int ny, int nz, int *ncc, int **S, int **L);
extern void Connected_Component_location(char *im, int nx, int ny, int nz, int *ncc, int **S, int **L, int **CCsize);

#endif

#ifndef _statistics
extern float imageMean(short *im, short *msk, int nv);
extern float imageMean(short *im, int nv);
extern double one_sample_t(double *x, int n);

extern float median(float *x, char *mask, int n);
extern void scaleAbsToOne(double *y, int n, int p);
extern void scaleAbsToOne(float *y, int n, int p);
extern double scaleAbsToOne(double *y, int n);
extern void decomposeVector(double *x, double *xpar, double *xper, double *u, int n);

extern void removeComponent(double *x, double *y, int n);
extern void removeComponent(float *x, float *y, int n);

extern double removeVectorMean(float *y, int n);
extern double removeVectorMean(double *y, int n);
extern double removeVectorMean(float *x, float *y, int n);
extern double removeVectorMean(double *x, double *y, int n);
extern double removeVectorMean(short *x, double *y, int n);
extern void removeVectorMean(double *y, int n, int p);
extern void removeVectorMean(float *y, int n, int p);

extern void partialCorrelation(double *Y, double *X1, double *X2, int n, double *pr1, double *pr2);
extern void partialCorrelation(float *Y, float *X1, float *X2, int n, double *pr1, double *pr2);
extern void partialCorrelation(short *Y, float *X1, float *X2, int n, double *pr1, double *pr2);

#endif

#ifndef _ginverse
extern int ginverse(float *X, int N, int p, float *G);
extern int ginverse(double *X, int N, int p, float *G);
#endif

#ifndef _convolution
extern float conv_pnt_sk(uchar *x,int sx,float *h,int sh,int i0);
extern float conv_pnt_sk(short *x,int sx,float *h,int sh,int i0);
extern float conv_pnt_sk(float *x,int sx,float *h,int sh,int i0);
extern float *conv_sk(short *x,int sx,float *h,int sh);
extern float *conv_sk(float *x,int sx,float *h,int sh);
extern void conv_sk(float *x, float *y, int sx,float *h,int sh);
extern void conv_sk_inplace(float *x, int sx, float *h, int sh);
#endif

#ifndef _hpsort
extern void hpsort(unsigned long n, float *ra);
extern void hpsort(unsigned long n, float *ra, int *indx);
#endif

#ifndef _random
extern void initializeRandomNumberGenerator();
extern void sampleWithReplacementIndex(int *I, int N);
#endif

#ifndef _dicomIO

#define IMPLICIT_LITTLE_ENDIAN 0
#define EXPLICIT_LITTLE_ENDIAN 1
#define EXPLICIT_BIG_ENDIAN 2

extern int readFileMetaInfo(const char *filename, DICOM_file_meta_info *file_meta_info, char v);
extern int readDataSet(const char *filename, DICOM_hdr *hdr, char v);

extern void readDicomInfo(const char *file, int np, dicominfo *info);
extern int readImageParams(const char *file, float *TLHC, float *rowvec, float *colvec, 
float *dx, float *dy, float *dz, int *TE, char *patientID, int *imageNumber, int *seriesNumber, int np);
extern int readPhaseEncodingDirection(const char *file, char *PED, int np);
extern int readPixelData(const char *file, char *data, int opt_j, int np);
extern int readTLHC(const char *file, float *TLHC, int np);
extern int readRowCol(const char *file, float *rowvec, float *colvec, int np);
extern int readNEX(const char *file, int *nex);
extern int readSliceThickness(const char *file, float *sliceThickness);
extern int readSeriesNumber(const char *file, int *seriesNumber, int np);
extern int readImageNumber(const char *file, int *imageNumber, int np);
extern int readTE(const char *file, int *TE, int np);
extern int readTR(const char *file, int *TR);
extern int readImageVoxelSize(const char *file, float *dx, float *dy, float *dz, int np);
extern int readImageMatrixSize(const char *file, ushort *nx, ushort *ny);
extern int readPatientID(const char *file, char *patientID, int np);
extern int dicomFormat(const char *file);
extern int readMetaFileInfo(const char *file, long *offset, int *syntax);
extern int readTag(const char *file, ushort iGN, ushort iEN, long byteOffset, int transferSyntax,
unsigned long *oVL, char *oV, long *valueOffset);
extern int read_element(char *filename, short S_GN, short S_EN, char *V, int *VL);
extern void readMatrixSize(char *filename, int *nx, int *ny);
extern void readVoxelSize(char *filename, float *dx, float *dy, float *dz);
extern int readImageSliceThickness(const char *file, float *dz, int np);

#endif

#ifndef _subsets
extern int *subsets(int N, int M);
extern int binomialCoeff(int N,int M);
#endif

#ifndef _cubicspline
extern float cubicSplineSynthesis(float *c, int nx, int ny, int nz, float x, float y, float z, float *beta, float del);

extern short *resliceImageCubicSpline(short *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *T);

extern float *resliceImageCubicSpline(float *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *T);

extern float *computeBeta(float *del);

extern void cubicSplineAnalysis(unsigned char *s, float *c, int nx, int ny, int nz);
extern void cubicSplineAnalysis(short *s, float *c, int nx, int ny, int nz);
extern void cubicSplineAnalysis(float *s, float *c, int nx, int ny, int nz);

extern void cubicSplineAnalysis(float *s, float *c, int N);

#endif

#ifndef _registration
extern void label_3d_cc(short *KMI,ushort label,int i,int j, int k,int *size,short CC, struct im_params *IP);
extern void set_transformation(float x, float y, float z, float ax, float ay, float az, const char *code, float *T);
extern int label_CCI(short *KMI, int size_thresh,struct im_params * IP, int nvoxels);
extern float (*interpolator)(float x, float y, float z, short *array, int nx, int ny, int nz, int np);
extern float P[12];
extern struct im_params IP;

extern void scale_short_minmax(short *imagein, unsigned char **imageout, int np, int min, int max);

extern void testCostFunc1(short *trg, int Tnx, int Tny, int Tnz, float Tdx, float Tdy, float Tdz,
short *obj, int Onx, int Ony, int Onz, float Odx, float Ody, float Odz);

extern float *transformation(float x, float y, float z, float ax, float ay,
float az, float sx, float sy, float sz, int rX, int rY, int rZ, char *code);

extern void cca(short *im, int nx, int ny, int nz);

extern int findThresholdLevel(short *image_in, int nv);
extern float *findTransMatrix(short *trg, int Tnx, int Tny, int Tnz, float Tdx, float Tdy, float Tdz,
short *obj, int Onx, int Ony, int Onz, float Odx, float Ody, float Odz);

extern short *KMcluster(short *ccImage, short *im_in, int nclass, int maxiter, int thresh, int Low, int High, int nv);
#endif

#ifndef _legendre
//extern float *FourierLegendreSynthesis(double *c, int nx, int ny, int nz, int mx, int my, int mz);
//extern double *FourierLegendreAnalysis(float *f, int nx, int ny, int nz, int mx, int my, int mz, int N);
extern double *LegendreAnalysis(float *image, int nx, int ny, int nz, int mx, int my, int mz);
extern	double *LegendreAnalysis(float *image, int nx, int ny, int mx, int my);
extern void LegendreSynthesis(double *c, int mx, int my, float *image, int nx, int ny);
extern float *LegendreSynthesis(double *c, int nx, int ny, int nz, int mx, int my, int mz);
void LegendrePoly(double *p0, double *q0, double *p1, double *q1, double *x, int N, int n);
void LegendrePoly(float *p0, float *q0, float *p1, float *q1, float *x, int N, int n);
void integral_1d(double *a, short *b, int n, double *d);
void integral_1d(double *a, float *b, int n, double *d);
#endif

#ifndef _matrixCom

extern double Frobenius_s3(double *qqT, float *D);
extern double Frobenius_s3(double *qqT, double *D);
extern void p3update(double *D, double *W, double epsilon);
extern void p3update(double *D, double *invsqrtD, double *sqrtD, double *W, double epsilon);
extern void s3multi(double *A, double *B, double *AB);
extern void p3invsqrt(double *A, double *invsqrtA, double *sqrtA);
extern void p3invsqrt(double *A, double *invsqrtA);
extern int p3invsqrt(float *A, float *invsqrtA);
extern double p3RiemannianDistance(double *D, double *invsqrtF);
extern double p3RiemannianDistance(double *L);
extern float p3RiemannianDistance(float *D, float *invsqrtF);
extern void s3ABA(double *A, double *B, double *ABA);
extern void s3ULUT(double *L, double *UT, double *ULUT);
extern void s3eigenvec(double *A, double *evalue, double *evector);
extern void getcomplement1(double *a, double *b, double *c);
extern void getcomplement2 (double *U, double *V, double *W);
extern void s3inv(double *A, double *invA);
extern int s3inv(float *A, float *invA);
extern void s3multi(double *A, double *B, double *AB);
extern void differentiate_distance(double *D, double *F, double *L, double *dLdD);

extern void subtractRowAvg(float *X, int N, int P, float *X0);
extern void copyVector(float *v1, float *v2, int n);
extern void subtractVector(float *v1, float *v2, int n);
extern int centerMatrixRow(float *X, int N, int P, float *avg);
extern int centerMatrixRow(float *X, int N, int P);

// Returns 1 if an error condition occurs, 0 otherwise

// returns an N-vector in "avg" containing the average of the P columns of X
int avgCol(float *X, int N, int P, float *avg);
int subtractAvgCol(float *X, int N, int P, float *avg);

extern int avgRow(float *X, int N, int P, float *avg, char *rowmask);
extern int avgRow(double *X, int N, int P, double *avg, char *rowmask);
extern int avgRow(float *X, int N, int P, float *avg);
extern int avgRow(double *X, int N, int P, double *avg);
extern int varRow(float *X, int N, int P, float *avg, float *var);
extern int varRow(double *X, int N, int P, double *avg, double *var);
extern int varRow(double *X, int N, int P, double *avg, double *var, char *rowmask);
extern int varRow(float *X, int N, int P, float *avg, float *var, char *rowmask);
extern int ssdRow(float *X, int N, int P, float *avg, float *ssd, char *rowmask);
extern int ssdRow(float *X, int N, int P, float *avg, float *ssd);

// compute the Euclidian distance between two vectors r0 and r1
extern double euclideandistance(float *r0, float *r1, int n);
extern double xtAx(float *A, double *x, int p);
extern float normalize(float *s, int n);
extern double normalize(double *s, int n);

extern int ComputeRank(double *M);
extern void s3eigenval(double *A, double *L);
extern double s3tr(double *A, double *B);
extern float s3tr(float *A, float *B);
extern void s3vec_to_mat(float *M, float *V);
extern void s3vec_to_mat(double *M, double *V);
extern void s3mat_to_vec(float *M, float *V);
extern void s3mat_to_vec(float *M, double *V);
extern void s3mat_to_vec(double *M, double *V);
extern void s3adjoint(double *A, double *ADJ);
extern void s3adjoint(float *A, float *ADJ);
extern double s3det(double *A);
extern float s3det(float *A);
void ds3det(double *A, double *B);
void ds3det(float *A, float *B);
extern void multi(float *A,int iA,int jA,float *B,int iB,int jB,float *C);
extern void multi(double *A,int iA,int jA,double *B,int iB,int jB,double *C);
extern void multi(float *A,int iA,int jA, double *B,int iB,int jB,double *C);
extern void multi(double *A,int iA,int jA, float *B,int iB,int jB,float *C);
#endif

//////////////////////////////////////////////////////////////////////////////
// The following are define in reslice.c
#define LIN 1
#define NEARN 2
#define SINC 3
#define CUBICSPLINE 4	

void resliceImage(SHORTIM im1, SHORTIM &im2, float *T, int interpolation_method);

short *resliceImage(short *im1, DIM dim1, DIM dim2, float *T, int interpolation_method);
float *resliceImage(float *im1, DIM dim1, DIM dim2, float *T, int interpolation_method);

float *resliceImage(float *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *T);

float *resliceImage(float *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *T, float *xjit, float *yjit);

float partial_var(float x, float y, float z, uchar *array, int nx, int ny, int nz, int np, float mu);

// you must initialize drand48 before using this function
uchar PNN(float x, float y, float z, uchar *array, int nx, int ny, int nz);

char *resliceImage(char *obj, int Onx, int Ony, float Odx, float Ody, int Tnx, int Tny, float Tdx, float Tdy, float *T);
short *resliceImage(short *im1, int nx1, int ny1, float dx1, float dy1, int nx2, int ny2, float dx2, float dy2, float *T);
short *resliceImage(float *im1, int nx1, int ny1, float dx1, float dy1, int nx2, int ny2, float dx2, float dy2, float *T);

short *resliceImage(short *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *T, int interpolation_method);

uchar *resliceImage(uchar *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *T);

float *resliceImage(float *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *T, float *w);

short *resliceImage(short *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *Xwarp, float *Ywarp, float *Zwarp);

float linearInterpolator(float x, float y, float z, float *array, int nx, int ny, int nz, int np, float *w);
uchar  linearInterpolator(float x, float y, float z, uchar  *array, int nx, int ny, int nz, int np, float *w);

short *computeReslicedImage(short *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *Xwarp, float *Ywarp, float *Zwarp);

float *computeReslicedImage(float *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *Xwarp, float *Ywarp, float *Zwarp);

short *computeReslicedImage(float *im1, int nx1, int ny1, float dx1, float dy1,
int nx2, int ny2, float dx2, float dy2, float *Xwarp, float *Ywarp);

short *computeReslicedImage(short *im1, int nx1, int ny1, float dx1, float dy1,
int nx2, int ny2, float dx2, float dy2, float *Xwarp, float *Ywarp);
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// defined in resize.c
short *resizeXYZ(short *image1,  DIM dim1, DIM dim2);

short *resizeXYZ(short *image1,
int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2);

uchar *resizeXYZ(uchar *image1,
int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2);

short *resizeXYZ(char *image1,
int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2);

float *resizeXYZ(float *image1,
int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2);

float *resizeXY(float *image1, int nx1, int ny1, float dx1, float dy1, int nx2, int ny2, float dx2, float dy2);
short *resizeXY(short *image1, int nx1, int ny1, float dx1, float dy1, int nx2, int ny2, float dx2, float dy2);
//////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
// The following functions are defined in analyzeio.c
int extension_is_hdr(const char *filename);
void read_analyze_image(const char *filename, short *im);
char *read_analyze_image(const char *filename, DIM *dim, int *type, int v);
float read_dx(const char *hdrfile);
float read_dy(const char *hdrfile);
float read_dz(const char *hdrfile);
int read_nt(const char *hdrfile);
int read_nx(const char *hdrfile);
int read_ny(const char *hdrfile);
int read_nz(const char *hdrfile);
int read_datatype(char *hdrfile);
char *read_analyze_image(const char *filename, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz);
char *read_analyze_image(const char *filename, int *nx, int *ny, int *nz, int *nt, float *dx, float *dy, float *dz, int *type, int v);
char *read_analyze_image(const char *filename, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz, int *type, int v);
char *read_analyze_image(const char *filename, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz, int *type);
char *read_image(char *file,int n);
void get_analyze_file_names(const char *filename, char *basename_hdr, char *basename_img);
void read_analyze_hdr(struct dsr *hdr, char *filename);
void setDimensions(struct dsr hdr, int *nx, int *ny, int *nz, double *dx, double *dy, double *dz, short *dataType);
void setDimensions(struct dsr hdr, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz, short *dataType);
void setDimensions(struct dsr hdr, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz, short *dataType, int v);
void setDimensions(struct dsr hdr, int *nx, int *ny, int *nz, int *nt, float *dx, float *dy, float *dz, short *dataType, int v);
void create_analyze_hdr(struct dsr *hdr, int nx, int ny, int nz, int dt, float dx, float dy, float dz);
void create_analyze_hdr(struct dsr *hdr, int nx, int ny, int nz, int nt, int datatype, float dx, float dy, float dz);
void write_analyze_image(const char *filename, short *im, int nx, int ny, int nz, float dx, float dy, float dz); 
void write_analyze_image(const char *filename, float *im, int nx, int ny, int nz, float dx, float dy, float dz); 
void write_analyze_image(const char *filename, uchar *im, int nx, int ny, int nz, float dx, float dy, float dz); 
void write_analyze_image(const char *filename, uchar *im, int nx, int ny, int nz, float dx, float dy, float dz, int v); 
void write_analyze_image(const char *filename, short *im, int nx, int ny, int nz, float dx, float dy, float dz,int v); 
void write_analyze_image(const char *filename, float *im, int nx, int ny, int nz, float dx, float dy, float dz,int v); 
////////////////////////////////////////////////////////////////////////////////////////

#ifndef _medianfilter
extern void medianFilter(float *image1, int nx, int ny, int nz, int Wx, int Wy, int Wz);
#endif

#ifndef _fileinfo
extern int isregular(const char *file);
extern int getFileSize(const char *file);
extern int checkWriteAccess(char *file);
extern int checkReadAccess(char *file);
extern int checkFileExistence(char *file);
extern int checkFileReadOK(char *file);
extern int checkFileWriteOK(char *file);
extern int check_F_R_permission(char *file);
#endif

#ifndef _histogram
extern int otsu(double *histogram, int numberOfBins);
extern double *findHistogram(short *im, int nv, int *nb, short &min, short &max);
extern double *findHistogram(short *im, int nv, int nb, int low, int high, int *bw_return);
extern double *findHistogram(short *im1, short *im2, int nv, int nb1, int nb2, int *bw1_r, int *bw2_r, int low1, int high1, int low2, int high2);
extern void trimExtremes(short *image, short *msk, int nv, float percent);
#endif

#ifndef _matrixops
extern void zeroVector(float *v, int n);
extern void zeroVector(char *v, int n);
extern void oneVector(float *v, int n);
extern void oneVector(char *v, int n);
extern void svd(float *At, int M, int N, float *Ut, float *V, float *S);
extern int zeroRowCol(float *A, int N, int n);
extern int setRowCol(float *A, int N, int n, float *a);
extern void projectVector(double *x, double *xpar, double *xper, float *Pz, int n);
extern float *projectionMatrix(double *X, int N, int p, int *rank);
extern float *projectionMatrix(float *X, int N, int p);
extern void mat_mat_trans(float *A,int Ar,int Ac,float *B,int Br, float *C);
extern float *diagATA_float(float *ATA, int n, char uplo);
extern float *AAT_float(float *A,int nr,int nc, char uplo);
extern void mat_trans_mat(float *A, int Ar, int Ac, float *B, int Bc, float *C);
extern void mat_trans_mat(double *A, int Ar, int Ac, double *B, int Bc, double *C);
#endif

///////////////////////////////////////////////////////////////
// The following functions are defined in nifti.cxx
int same_nifti_image_size(int N, char **imagefile, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz);
void read_nifti_image(const char *filename, uchar **im, nifti_1_header *hdr);
void read_nifti_image(const char *filename, short **im, nifti_1_header *hdr);
void print_NIFTI_hdr(const char *filename);
void print_NIFTI_hdr(nifti_1_header hdr);
nifti_1_header read_NIFTI_hdr(const char *filename, nifti1_extender *extender, char **extension);
void save_nifti_image(const char *filename, uchar *im, nifti_1_header *hdr);
void save_nifti_image(const char *filename, short *im, nifti_1_header *hdr);
void save_nifti_image(const char *filename, float *im, nifti_1_header *hdr);
void save_nifti_image(const char *filename, char *im, nifti_1_header *hdr);
// returns the orientations vectors xvec, yvec, and zvec in NIFTI's RAS system
void readOrientationVectorsFromFile(const char *filename, float *xvec, float *yvec, float *zvec);
short *readNiftiImage(const char *filename, DIM *dim, int flg);
///////////////////////////////////////////////////////////////

#ifndef _utils

extern int hand_system(char *code);

extern void art_to_fsl(float *Mart, float *Mfsl, DIM sub_dim, DIM trg_dim);
extern void fsl_to_art(float *Mfsl, float *Mart, DIM sub_dim, DIM trg_dim, int subflg, int trgflg);

extern void sobel_edge_x(short *in, float *out, int nx, int ny);
extern void sobel_edge_y(short *in, float *out, int nx, int ny);
extern void sobel_edge(short *in, float *out, int nx, int ny);
extern void sobel_edge(short *in, short *out, int nx, int ny);

extern int ccsize(short *im, int nv);

extern void copyarray(float *source, float *destination, int size);
extern void copyarray(short *source, char *destination, int size);
extern void copyarray(char *source, short *destination, int size);

extern void zeroarray(float *y, int size);
extern float diceindex(short *setA, short *setB, int n);
extern void remove_space(char *inp);
extern void orientationCodeConverter(int integerCode, char *stringCode);
extern int orientationCodeOK(char *stringCode);

extern void saveMatrix(float *A, int n, int m, char *filename);
extern void saveMatrix(short *A, int n, int m, char *filename);
extern float *readMatrix(int *n, int *m, char *filename);

extern float *readDataMatrix(char **imageList, int n, int p, short *mask);
extern float *readDataMatrix_nifti(char **imageList, int n, int p, short *mask);

extern short *readDataMatrixShort(char **imageList, int n, int p, short *mask);
extern short *readDataMatrixShort_nifti(char **imageList, int n, int p, short *mask);

extern short *readMask(const char *filename, int *nx, int *ny, int *nz);
extern short *readMask_nifti(const char *filename, int *nx, int *ny, int *nz);


extern void checkDimension(int N, char **imagefile, int nx, int ny, int nz);
extern void checkDimension_nifti(int N, char **imagefile, int nx, int ny, int nz);
extern void checkDimension(int N, char **imagefile, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz);


extern void affineLSE(short *msk, int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T);
extern void affineLSE(char *msk, int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T);
extern float *affineLSE(char *msk, int nx, int ny, float dx, float dy, float *Xwarp, float *Ywarp);
extern void affineLSE(char *msk, int nx, int ny, float dx, float dy, float *Xwarp, float *Ywarp, float *T);

extern void affineLSE(short *msk, int nx, int ny, float dx, float dy, float *Xwarp, float *Ywarp, float *T);

extern void extractArray(float *im, int nx, int ny, int nz, int np, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, float *array);
extern void extractArray(uchar *im, int nx, int ny, int nz, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, float *array);
extern void extractArray(float *im, int nx, int ny, int nz, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, float *array);
extern void extractArray(short *im, int nx, int i,  int j,  int L, float *array);
extern void extractArray(short *im, int nx, int ny, int nz, int x0, int y0, int z0, short *x,short *y,short *z,int n,float *array);
extern void extractArray(short *im, int nx, int ny, int nx0,int ny0,int Lx, int Ly, float *array);
extern void extractArray(short *im, int nx, int ny, int nz, int nx0,int ny0, int nz0, int Lx, int Ly, int Lz, float *array);
extern int extractArray(short *im, int nx, int ny, int nz, int np, int nx0,int ny0, int nz0, int Lx, int Ly, int Lz, short *array);
extern void extractArray(short *im, int nx, int ny, int nz, int nx0,int ny0,int nz0, int Lx, int Ly, int Lz, short *array);
extern void extractArray(short *im, int nx, int ny, int nz, int np, int nx0,int ny0, int nz0, int Lx, int Ly, int Lz, float *array);

// This function extracts the "filename" from the full "path" string.
// For example, if path="/home/babak/testimages/test.img", then filname="test.img".
extern void getfilename(char *filename, const char *path);

extern void printMatrix(int *mat, int n, int p, const char *s, FILE *fp);
extern void printMatrix(float *mat, int n, int p, const char *s, FILE *fp);
extern void printMatrix(double *mat, int n, int p, const char *s, FILE *fp);
extern void get_temp_filename(char *filename);
extern void mask_and_save(const char *inputfile, const char *outputfile, short *mask, short *masked_image, int nbv, float FWHM);
extern void mask_and_save_nii(const char *inputfile, const char *outputfile, short *mask, short *masked_image, int nbv, float FWHM);
extern void read_transpose_save(char *inputfile, char *outputfile, int nr, int v);
extern void centerOfMass(short *im, int nx, int ny, int nz, float dx, float dy, float dz, float *CM);

#endif

////////////////////////////////////////////////////////////////////////
// defined in smooth.cxx
float *smoothY(float *image, int nx, int ny, int nz, float sd);
float *smoothZ(float *image, int nx, int ny, int nz, float sd);
////////////////////////////////////////////////////////////////////////


#endif
