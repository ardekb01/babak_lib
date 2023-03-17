#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <babak_lib.h>
#include <minmax.h>
#include <stats.h>
#include <interpolator.h>
#include <nifti1_io.h>

#define _artlib

#define MAXITER 5000
#define DEL 3

// Turns verbose mode on (opt_v=YES) or off (opt_v=NO)
char opt_v= NO;

char opt_ppm = YES;
char opt_png = YES; // flag for outputing PNG images
char opt_txt = YES;
char opt_AC = YES;
char opt_PC = YES;
char opt_RP = YES;
char opt_MSP = YES;

static float Sff,Sf;
static short *gimage;
static DIM gdim;
static float vertex[4][3];
static float value[4];
static float vertexSum[3];
static float VertexNew[3];

void getDirectoryName(const char *pathname, char *dirname);
void orig_ijk_to_pil_xyz(float *Tmsp, DIM orig_dim, float *AC, float *PC);
void ACPCtransform(float *Tacpc, float *Tmsp, float *AC, float *PC, char flg);
void brandImage(unsigned char *R, unsigned char *G, unsigned char *B, int nx, int ny, int sx, int sy, int L1, int L2, unsigned char Rvalue, unsigned char Gvalue, unsigned char Bvalue);
void saveACPCimages(const char *imagefilename, char *ACregion, char *PCregion, char *RPregion,  
float *AC, float *PC, float *RP, DIM HR, DIM Orig, short *volOrig, float *Tmsp);
void compute_MSP_parameters_from_Tmsp(float *Tmsp, float *n, float *d);
void compute_Tmsp_from_MSP_parameters(const char *orientation, float *Tmsp, float *n, float d);
void updateTmsp(const char *imagefilename, float *Tmsp, float *RP, float *AC, float *PC);
float detectPC(float *PC, char *modelfile, short *volumeMSP_HR, char *PCregion, short *xPC, short *yPC, short *zPC, int opt_T2);
float detectAC(float *AC, char *modelfile, short *volumeMSP_HR, char *ACregion, short *xAC, short *yAC, short *zAC, int opt_T2);
void initialAC(float Ax, float Ay, float Bx, float By, float *Cx, float *Cy, float parcomMean, float percomMean);
char *defineACregion(DIM dim, float *RP, float *PC, float parcomMean, float percomMean, double ACsr);
char *definePCregion(DIM HR, float *RP, float *RPPCmean, double PCsr);
void detectRP(float *RP1, float *RP2, char *modelfile, short *volumeMSP_LR, short *mask_LR, short *xRP, short *yRP, short *zRP);
void defineTemplate(int r, int h, short *x, short *y, short *z);
char *expandMask(short *mask_HR, DIM HR, float *RPmean, double RPsr);
short *thresholdImageOtsu(short *im, int nv, int *nbv);
int detect_AC_PC_MSP(const char *imagefilename, char *orientation, char *modelfile,
float *AC, float *PC, float *RP, float *Tmsp, int opt_v, int opt_T2);
float reflectVertex(int pmax, float fac);
int dsm(void);
float reflection_cross_correlation2(short *image, DIM dim, float A, float B, float C);
float optimizeNormalVector(short *image,DIM dim,float *A, float *B, float *C);
float reflection_cross_correlation(short *image, DIM dim, float a, float b, float c, float d);
void findInitialNormalVector(short *image, DIM dim, float *A, float *B,float *C);
float msp(short *im_in, int nx, int ny, int nz, float dx, float dy, float dz, float *A, float *B, float *C);
void combine_warps_and_trans(int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T);

///////////////////////////////////////////////////////////////////////////////////
void sub2trg_rigid_body_transformation(float *sub2trg, const char *subfile, const char *trgfile)
{
   nifti_1_header trg_hdr; 
   nifti_1_header sub_hdr;

   // Takes a point (x, y, z) in the original orientation of the subject volume
   // to point (i, j, k), also in the original orientation of the subject volume
   float sub_xyz2ijk[16];

   // Takes a point (i, j, k) in the original orientation of the subject volume
   // to point (x, y, z) in RAS orientation. 
   float sub_ijk2RAS[16];

   // Takes a point (i, j, k) in the original orientation of the trg volume
   // to point (x, y, z) in RAS orientation.
   float trg_ijk2RAS[16];

   // Takes a point (x, y, z) in RAS orientation 
   // to point (i, j, k) in the original orientation of the trg volume
   float *trg_RAS2ijk;  // inverse of trg_ijk2RAS;

   // Takes a point (i, j, k) in the original orientation of the trg volume
   // to point (x, y, z), also in the original orientation of the trg volume
   float trg_ijk2xyz[16];

   // (i, j, k)_sub = sub_xyz2ijk * (x, y, z)_sub
   // (x, y, z)_RAS = sub_ijk2RAS * (i, j, k)_sub
   // (i, j, k)_trg = trg_RAS2ijk * (x, y, z)_RAS
   // (x, y, z)_trg = trg_ijk2xyz * (i, j, k)_trg

   // combining the two transformations gives
   // (x, y, z)_trg = trg_ijk2xyz * trg_RAS2ijk * sub_ijk2RAS * sub_xyz2ijk * (x, y, z)_sub

   sub_hdr = read_NIFTI_hdr(subfile);
   trg_hdr = read_NIFTI_hdr(trgfile);

   xyz2ijk(sub_xyz2ijk, sub_hdr.dim[1], sub_hdr.dim[2], sub_hdr.dim[3], 
   sub_hdr.pixdim[1], sub_hdr.pixdim[2], sub_hdr.pixdim[3]);

   if(sub_hdr.qform_code > 0)
   {
      mat44 R;

      R = nifti_quatern_to_mat44(sub_hdr.quatern_b, sub_hdr.quatern_c, sub_hdr.quatern_d, 
      sub_hdr.qoffset_x, sub_hdr.qoffset_y, sub_hdr.qoffset_z, sub_hdr.pixdim[1], sub_hdr.pixdim[2], 
      sub_hdr.pixdim[3], sub_hdr.pixdim[0]);

      sub_ijk2RAS[0]= R.m[0][0];
      sub_ijk2RAS[1]= R.m[0][1];
      sub_ijk2RAS[2]= R.m[0][2];
      sub_ijk2RAS[3]= R.m[0][3];

      sub_ijk2RAS[4]= R.m[1][0];
      sub_ijk2RAS[5]= R.m[1][1];
      sub_ijk2RAS[6]= R.m[1][2];
      sub_ijk2RAS[7]= R.m[1][3];

      sub_ijk2RAS[8]= R.m[2][0];
      sub_ijk2RAS[9]= R.m[2][1];
      sub_ijk2RAS[10]=R.m[2][2];
      sub_ijk2RAS[11]= R.m[2][3];
   }
   else if(sub_hdr.sform_code > 0)
   {
      sub_ijk2RAS[0]= sub_hdr.srow_x[0]; 
      sub_ijk2RAS[1]= sub_hdr.srow_x[1];
      sub_ijk2RAS[2]= sub_hdr.srow_x[2];
      sub_ijk2RAS[3]= sub_hdr.srow_x[3];

      sub_ijk2RAS[4]= sub_hdr.srow_y[0]; 
      sub_ijk2RAS[5]= sub_hdr.srow_y[1];
      sub_ijk2RAS[6]= sub_hdr.srow_y[2];
      sub_ijk2RAS[7]= sub_hdr.srow_y[3];

      sub_ijk2RAS[8]= sub_hdr.srow_z[0]; 
      sub_ijk2RAS[9]= sub_hdr.srow_z[1];
      sub_ijk2RAS[10]= sub_hdr.srow_z[2];
      sub_ijk2RAS[11]= sub_hdr.srow_z[3];
   }
   else
   {
      printf("Warning: the subject volume header lacks qform or sform information!\n");
      return;
   }

   sub_ijk2RAS[12]= 0.0; 
   sub_ijk2RAS[13]= 0.0;
   sub_ijk2RAS[14]= 0.0; 
   sub_ijk2RAS[15]= 1.0;

   if(trg_hdr.qform_code > 0)
   {
      mat44 R;

      R = nifti_quatern_to_mat44(trg_hdr.quatern_b, trg_hdr.quatern_c, trg_hdr.quatern_d, 
      trg_hdr.qoffset_x, trg_hdr.qoffset_y, trg_hdr.qoffset_z, trg_hdr.pixdim[1], trg_hdr.pixdim[2], 
      trg_hdr.pixdim[3], trg_hdr.pixdim[0]);

      trg_ijk2RAS[0]= R.m[0][0];
      trg_ijk2RAS[1]= R.m[0][1];
      trg_ijk2RAS[2]= R.m[0][2];
      trg_ijk2RAS[3]= R.m[0][3];

      trg_ijk2RAS[4]= R.m[1][0];
      trg_ijk2RAS[5]= R.m[1][1];
      trg_ijk2RAS[6]= R.m[1][2];
      trg_ijk2RAS[7]= R.m[1][3];

      trg_ijk2RAS[8]= R.m[2][0];
      trg_ijk2RAS[9]= R.m[2][1];
      trg_ijk2RAS[10]=R.m[2][2];
      trg_ijk2RAS[11]= R.m[2][3];
   }
   else if(trg_hdr.sform_code > 0)
   {
      trg_ijk2RAS[0]= trg_hdr.srow_x[0]; 
      trg_ijk2RAS[1]= trg_hdr.srow_x[1];
      trg_ijk2RAS[2]= trg_hdr.srow_x[2];
      trg_ijk2RAS[3]= trg_hdr.srow_x[3];

      trg_ijk2RAS[4]= trg_hdr.srow_y[0]; 
      trg_ijk2RAS[5]= trg_hdr.srow_y[1];
      trg_ijk2RAS[6]= trg_hdr.srow_y[2];
      trg_ijk2RAS[7]= trg_hdr.srow_y[3];

      trg_ijk2RAS[8]= trg_hdr.srow_z[0]; 
      trg_ijk2RAS[9]= trg_hdr.srow_z[1];
      trg_ijk2RAS[10]= trg_hdr.srow_z[2];
      trg_ijk2RAS[11]= trg_hdr.srow_z[3];
   }
   else
   {
      printf("Warning: the target volume header lacks qform or sform information!\n");
      return;
   }

   trg_ijk2RAS[12]= 0.0; 
   trg_ijk2RAS[13]= 0.0;
   trg_ijk2RAS[14]= 0.0; 
   trg_ijk2RAS[15]= 1.0;

   ijk2xyz(trg_ijk2xyz, trg_hdr.dim[1], trg_hdr.dim[2], trg_hdr.dim[3], 
   trg_hdr.pixdim[1], trg_hdr.pixdim[2], trg_hdr.pixdim[3]);

   float temp_mat[16];
   multi(sub_ijk2RAS, 4, 4, sub_xyz2ijk, 4, 4, temp_mat);

   trg_RAS2ijk = inv4(trg_ijk2RAS);
   multi(trg_RAS2ijk, 4, 4, temp_mat, 4, 4, temp_mat);
   free(trg_RAS2ijk);

   multi(trg_ijk2xyz, 4, 4, temp_mat, 4, 4, sub2trg);
}
///////////////////////////////////////////////////////////////////////////////////

void getDirectoryName(const char *pathname, char *dirname)
{
	int n;
	int i;

	n = (int)strlen(pathname);

	for(i=n-1;i>=0;i--)
	if( pathname[i] == '/') break;

	if(i==-1)
	{
		dirname[0]='.';
		dirname[1]='\0';
	}
	else
	{
		strncpy(dirname, pathname, i);
		dirname[i]='\0';
	}

	return;
}

// The AC and PC points are passed to this routine in ijk coordinates of the original image.
// The routine returns their location in xyz coordinates of a PIL transformed image.
void orig_ijk_to_pil_xyz(float *Tmsp, DIM orig_dim, float *AC, float *PC)
{
   float I2X[16];

   ijk2xyz(I2X, orig_dim.nx, orig_dim.ny, orig_dim.nz, orig_dim.dx, orig_dim.dy, orig_dim.dz);

   multi(I2X, 4, 4,  AC, 4,  1, AC);
   multi(I2X, 4, 4,  PC, 4,  1, PC);

   multi(Tmsp, 4, 4,  AC, 4,  1, AC);
   multi(Tmsp, 4, 4,  PC, 4,  1, PC);
}

// if flg=1 then Tacpc makes AC the center of the FOV
// if flg=0 then Tacpc makes the midpoint between AC and PC the center of the FOV
void ACPCtransform(float *Tacpc, float *Tmsp, float *AC, float *PC, char flg)
{
   float ac[4], pc[4];
   float ACPC[2];
   float R[16], T[16];

   // determined the AC-PC vector (a vector pointing from AC to PC)
   ACPC[0]=PC[0]-AC[0];
   ACPC[1]=PC[1]-AC[1];
		
   normalizeVector(ACPC,2);

   // determine the rotation vector necessary to align the AC-PC vector to the +x axis
   rotate(R, ACPC[0], -ACPC[1], 0.0, 0.0, 1.0);
   multi(R, 4,4, Tmsp, 4, 4, Tacpc);

   multi(R, 4, 4,  PC, 4,  1, pc);
   multi(R, 4, 4,  AC, 4,  1, ac);

   if(flg==1)
   {
      T[0]=1.0;  T[1]=0.0;  T[2]=0.0;  T[3]=-ac[0]; // makes ac the center of the FOV
   }
   else
   {
      T[0]=1.0;  T[1]=0.0;  T[2]=0.0;  T[3]=-(ac[0]+pc[0])/2.0;
   }
   T[4]=0.0;  T[5]=1.0;  T[6]=0.0;  T[7]=-ac[1];  // ac[1] should be equal to pc[1]
   T[8]=0.0;  T[9]=0.0;  T[10]=1.0; T[11]=0;
   T[12]=0.0; T[13]=0.0; T[14]=0.0; T[15]=1.0;

   multi(T, 4,4, Tacpc, 4, 4, Tacpc);
}

void brandImage(unsigned char *R, unsigned char *G, unsigned char *B, int nx, int ny, int sx, int sy, int L1, int L2, unsigned char Rvalue, unsigned char Gvalue, unsigned char Bvalue)
{
   if(sy>=0 && sy<ny && sx>=0 && sx<nx)
   {
      R[nx*sy + sx]=Rvalue;
      G[nx*sy + sx]=Gvalue;
      B[nx*sy + sx]=Bvalue;
   }

   for(int i=sx-L1; i<=(sx+L1); i++)
   if( sy>=0 && sy<ny && i>=0 && i<nx )
   {
      R[sy*nx + i]=Rvalue;
      G[sy*nx + i]=Gvalue;
      B[sy*nx + i]=Bvalue;
   }

   for(int j=sy-L2; j<=(sy+L2); j++)
   if( sx>=0 && sx<nx && j>=0 && j<ny )
   {
      R[j*nx + sx]=Rvalue;
      G[j*nx + sx]=Gvalue;
      B[j*nx + sx]=Bvalue;
   }
}

void saveACPCimages(const char *imagefilename, char *ACregion, char *PCregion, char *RPregion,  
float *AC, float *PC, float *RP, DIM HR, DIM Orig, short *volOrig, float *Tmsp)
{
   char filename[1024]; // filename variable for reading/writing data files
   unsigned char *Rchannel, *Gchannel, *Bchannel;
   int npHR;
   float TPIL2LPS[16];
   float Tacpc[16], *invT; 
   float R[16],T[16];
   float ACPC[2];
   float X2I[16]; // 4x4 matrix that transforms a vector from xyz-coordinates to ijk-coordinates
   float ac[4], pc[4], rp[4];
   short min=0, max=0;
   short *im;

   npHR = HR.nx * HR.ny;

   Rchannel = (unsigned char *)calloc(npHR, 1);
   Gchannel = (unsigned char *)calloc(npHR, 1);
   Bchannel = (unsigned char *)calloc(npHR, 1);

   inversePILtransform("LPS", TPIL2LPS);

   xyz2ijk(X2I, HR.nx, HR.ny, 1, HR.dx, HR.dy, HR.dz);

   multi(X2I, 4, 4,  RP, 4,  1, rp);
   multi(X2I, 4, 4,  PC, 4,  1, pc);
   multi(X2I, 4, 4,  AC, 4,  1, ac);
		
   invT = inv4(Tmsp);
   im=resliceImage(volOrig,Orig.nx,Orig.ny,Orig.nz,Orig.dx,Orig.dy,Orig.dz,HR.nx,HR.ny,1, HR.dx,HR.dy,HR.dz,invT,LIN);
   free(invT);

   minmax(im, npHR, min, max);

   for(int i=0; i<npHR; i++)
   {
      Rchannel[i] = (unsigned char)(im[i]*255.0/max);
      Gchannel[i] = (unsigned char)(im[i]*255.0/max);
      Bchannel[i] = (unsigned char)(im[i]*255.0/max);
   }

   brandImage(Rchannel, Gchannel, Bchannel, HR.nx, HR.ny, (int)(rp[0]+0.5), (int)(rp[1]+0.5), 4, 4, 0, 0, 255);
   brandImage(Rchannel, Gchannel, Bchannel, HR.nx, HR.ny, (int)(pc[0]+0.5), (int)(pc[1]+0.5), 4, 4, 255, 0, 0);
   brandImage(Rchannel, Gchannel, Bchannel, HR.nx, HR.ny, (int)(ac[0]+0.5), (int)(ac[1]+0.5), 4, 4, 0, 255, 0);

	for(int i=1; i<HR.nx-1; i++)
	for(int j=1; j<HR.ny-1; j++)
	{
		if( ACregion[npHR*(HR.nz-1)/2 + j*HR.nx + i]==1 && 
		(ACregion[npHR*(HR.nz-1)/2 + j*HR.nx + i + 1]==0 || 
		ACregion[npHR*(HR.nz-1)/2 + j*HR.nx + i - 1]==0 || 
		ACregion[npHR*(HR.nz-1)/2 + (j+1)*HR.nx + i ]==0 || 
		ACregion[npHR*(HR.nz-1)/2 + (j-1)*HR.nx + i ]==0) )
		{
			Rchannel[j*HR.nx + i]=0;
			Gchannel[j*HR.nx + i]=255;
			Bchannel[j*HR.nx + i]=0;
		}

		if( PCregion[npHR*(HR.nz-1)/2 + j*HR.nx + i]==1 && 
		(PCregion[npHR*(HR.nz-1)/2 + j*HR.nx + i + 1]==0 || 
		PCregion[npHR*(HR.nz-1)/2 + j*HR.nx + i - 1]==0 || 
		PCregion[npHR*(HR.nz-1)/2 + (j+1)*HR.nx + i ]==0 || 
		PCregion[npHR*(HR.nz-1)/2 + (j-1)*HR.nx + i ]==0) )
		{
			Rchannel[j*HR.nx + i]=255;
			Gchannel[j*HR.nx + i]=0;
			Bchannel[j*HR.nx + i]=0;
		}

		if( RPregion[j*HR.nx + i]==1 && 
		(RPregion[j*HR.nx + i + 1]==0 || 
		RPregion[j*HR.nx + i - 1]==0 || 
		RPregion[(j+1)*HR.nx + i ]==0 || 
		RPregion[(j-1)*HR.nx + i ]==0) )
		{
			Rchannel[j*HR.nx + i]=0;
			Gchannel[j*HR.nx + i]=0;
			Bchannel[j*HR.nx + i]=255;
		}
	}
   
   char dirname[DEFAULT_STRING_LENGTH]; // name of the directory only
   // char fullpath[DEFAULT_STRING_LENGTH]; // directory + filename
   char fullpath[2*DEFAULT_STRING_LENGTH+1]; // directory + filename
   getDirectoryName(imagefilename, dirname);
   niftiFilename(filename,imagefilename);
   snprintf(fullpath,sizeof(fullpath),"%s/%s_ACPC_sagittal.ppm",dirname,filename);

   if(opt_ppm || opt_png)
   {
      save_as_ppm((const char *)fullpath, HR.nx, HR.ny, Rchannel, Gchannel, Bchannel);
   }

   free(im);

	/////////////////////////////////////////////////////////////////////////////////////////////////

	// determined the AC-PC vector (a vector pointing from AC to PC)
	ACPC[0]=PC[0]-AC[0];
	ACPC[1]=PC[1]-AC[1];
		
	normalizeVector(ACPC,2);

	// determine the rotation vector necessary to align the AC-PC vector to the +x axis
	rotate(R, ACPC[0], -ACPC[1], 0.0, 0.0, 1.0);

	multi(R, 4, 4,  PC, 4,  1, pc);
	multi(R, 4, 4,  AC, 4,  1, ac);

	T[0]=1.0;  T[1]=0.0;  T[2]=0.0;  T[3]=-(ac[0]+pc[0])/2.0;
	T[4]=0.0;  T[5]=1.0;  T[6]=0.0;  T[7]=-ac[1];
	T[8]=0.0;  T[9]=0.0;  T[10]=1.0; T[11]=0;
	T[12]=0.0; T[13]=0.0; T[14]=0.0; T[15]=1.0;

	multi(T, 4, 4,  R, 4,  4, T);

	multi(T, 4, 4,  PC, 4,  1, pc);
	multi(T, 4, 4,  AC, 4,  1, ac);

	multi(T, 4,4, Tmsp, 4, 4, Tacpc);
	multi(TPIL2LPS,4, 4, Tacpc, 4, 4, T);

	multi(TPIL2LPS, 4, 4,  pc, 4,  1, pc);
	multi(TPIL2LPS, 4, 4,  ac, 4,  1, ac);

	multi(X2I, 4, 4,  pc, 4,  1, pc);
	multi(X2I, 4, 4,  ac, 4,  1, ac);

	invT = inv4(T);
	im=resliceImage(volOrig,Orig.nx,Orig.ny,Orig.nz,Orig.dx,Orig.dy,Orig.dz,HR.nx,HR.ny,1, HR.dx,HR.dy,HR.dz,invT,LIN);
	free(invT);

	minmax(im, npHR, min, max);

	for(int i=0; i<npHR; i++)
	{
		Rchannel[i] = (char)(im[i]*255.0/max);
		Gchannel[i] = (char)(im[i]*255.0/max);
		Bchannel[i] = (char)(im[i]*255.0/max);
	}

	brandImage(Rchannel, Gchannel, Bchannel, HR.nx, HR.ny, (HR.nx-1)/2, (HR.ny-1)/2, 0, (HR.nx-1)/2, 255, 255, 255);
	brandImage(Rchannel, Gchannel, Bchannel, HR.nx, HR.ny, (int)(pc[0]+0.5), (int)(pc[1]+0.5), 4, 4, 255, 0, 0);
	brandImage(Rchannel, Gchannel, Bchannel, HR.nx, HR.ny, (int)(ac[0]+0.5), (int)(ac[1]+0.5), 4, 4, 0, 255, 0);


   niftiFilename(filename,imagefilename);
   snprintf(fullpath,sizeof(fullpath),"%s/%s_ACPC_axial.ppm",dirname,filename);

   if(opt_ppm || opt_png)
   {
      save_as_ppm((const char *)fullpath, HR.nx, HR.ny, Rchannel, Gchannel, Bchannel);
   }

   free(im);

   free(Rchannel); free(Gchannel); free(Bchannel);
}

void saveACPClocation(const char *imagefilename, float *Tmsp, DIM Orig, float *AC, float *PC, float *RP, int opt_v)
{
   FILE *fp;
   float ac[4], pc[4], rp[4];
   float *invT; 
   float acpc_distance;
   float X2I[16]; // 4x4 matrix that transforms a vector from xyz-coordinates to ijk-coordinates
   char filename[1024];

   invT = inv4(Tmsp);

   acpc_distance = sqrtf( (AC[0]-PC[0])*(AC[0]-PC[0]) + (AC[1]-PC[1])*(AC[1]-PC[1]) + (AC[2]-PC[2])*(AC[2]-PC[2]));

   multi(invT, 4, 4,  AC, 4,  1, ac);
   multi(invT, 4, 4,  PC, 4,  1, pc);
   multi(invT, 4, 4,  RP, 4,  1, rp);
   free(invT);
 
   xyz2ijk(X2I, Orig);

   multi(X2I, 4, 4,  ac, 4,  1, ac);
   multi(X2I, 4, 4,  pc, 4,  1, pc);
   multi(X2I, 4, 4,  rp, 4,  1, rp);

   if(opt_v) 
   {
      printf("\n%s\n",imagefilename);
      printf("\tVSPS detected at (i,j,k) = (%5.1f, %5.1f, %5.1f)\n",rp[0],rp[1],rp[2]);
      printf("\tAC detected at (i,j,k) = (%5.1f, %5.1f, %5.1f)\n",ac[0],ac[1],ac[2]);
      printf("\tPC detected at (i,j,k) = (%5.1f, %5.1f, %5.1f)\n",pc[0],pc[1],pc[2]);
      printf("\tAC-PC distance = %6.3f mm\n",acpc_distance);
   }

   char dirname[DEFAULT_STRING_LENGTH]; // name of the directory only
   char fullpath[DEFAULT_STRING_LENGTH]; // directory + filename
   getDirectoryName(imagefilename, dirname);
   niftiFilename(filename,imagefilename);
   snprintf(fullpath,sizeof(fullpath),"%s/%s_ACPC.txt",dirname, filename);

   if(opt_txt)
   {
      fp=fopen(fullpath,"a"); // note that we are appending to this file here

      if(fp==NULL)
      {
         errorMessage("Warning: I have trouble opening the *_ACPC.txt file to save the output.");
      }
      else
      {
         fprintf(fp,"# AC-PC distance = %6.3f mm\n\n",acpc_distance);

         fprintf(fp,"# VSPS (i,j,k) voxel location:\n");
         fprintf(fp,"%5.1f %5.1f %5.1f\n\n",rp[0],rp[1],rp[2]);

         fprintf(fp,"# AC (i,j,k) voxel location:\n");
         fprintf(fp,"%5.1f %5.1f %5.1f\n\n",ac[0],ac[1],ac[2]);

         fprintf(fp,"# PC (i,j,k) voxel location:\n");
         fprintf(fp,"%5.1f %5.1f %5.1f\n\n",pc[0],pc[1],pc[2]);

         fclose(fp);
      }
   }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
// Tmsp is a 4x4 transformation matrix that maps the points in the original input image space
// to a new PIL space.  The detected MSP in the original image
// space is mapped to the z=0 plane of the new PIL space.

// n is the unit normal to the MSP in the original image xyz coordinate system
// I think n will always point towards the left since it's the inverse of (0,0,1) in the PIL system

// d is the normal distance between origin and the MSP in the original image xyz system

void compute_MSP_parameters_from_Tmsp(float *Tmsp, float *n, float *d)
{
   float *invTmsp;

   // a point on the MSP in the original image xyz coordinate system
   float p[3]; 

   invTmsp = inv4(Tmsp);

   n[0] = invTmsp[2];
   n[1] = invTmsp[6];
   n[2] = invTmsp[10];

   p[0] = invTmsp[3];
   p[1] = invTmsp[7];
   p[2] = invTmsp[11];

   *d = n[0]*p[0] + n[1]*p[1] + n[2]*p[2];

   free(invTmsp);
}

// n & d are assumed to be in the original xyz system
void compute_Tmsp_from_MSP_parameters(const char *orientation, float *Tmsp, float *n, float d)
{
   float R[16], T[16];
   float TPIL[16]; // Transformation from original to PIL orientation
   float nPIL[4];
   float alpha;

   // Compute the transformation from original to PIL orientation.
   PILtransform(orientation, TPIL);

   nPIL[0]=n[0]; nPIL[1]=n[1]; nPIL[2]=n[2]; nPIL[3]=1.0;

   multi(TPIL, 4, 4, nPIL, 4, 1, nPIL);

   // ensures that the MSP in PIL pass through the center of the FOV
   T[0]=1.0;  T[1]=0.0;  T[2]=0.0;  T[3]=-nPIL[0]*d;
   T[4]=0.0;  T[5]=1.0;  T[6]=0.0;  T[7]=-nPIL[1]*d;
   T[8]=0.0;  T[9]=0.0;  T[10]=1.0; T[11]=-nPIL[2]*d;
   T[12]=0.0; T[13]=0.0; T[14]=0.0; T[15]=1.0;

   if(nPIL[2]>1.0) nPIL[2]=1.0; // just in case to prevent acos from getting into trouble
   alpha = (float)acos((double)nPIL[2]);

   // ensures that MSP in PIL is parallel to the x-y plane
   rotate(R, alpha, nPIL[1], -nPIL[0], 0.0);

   // Tmsp=R*T*TPIL transforms the volOrig to MSP aligned PIL
   multi(T, 4, 4,  TPIL, 4,  4, Tmsp);
   multi(R, 4, 4,  Tmsp, 4,  4, Tmsp);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Tmsp is a 4x4 transformation matrix that maps the points in the original input image space
// to a new PIL space.  The transformation is such that the detected MSP in the original image
// space is mapped to the z=0 plane of the new PIL space.

// AC and PC are locations of the AC and PC landmarks in the PIL space in xyz coordinates system

// After updating, the Tmsp also maps the AC and PC landmarks in the original input image space
// to z=0 plane of the PIL space.

void updateTmsp(const char *imagefilename, float *Tmsp, float *RP, float *AC, float *PC)
{
   float R[16];

   // Assume that the AC and PC points are projected onto the zx-plane. Let `ACPC_zx[3]' be a vector that goes
   // from the projected AC to the projected PC on the zx-plane.
   float ACPC_zx[3]; 

   ACPC_zx[0] = PC[0]-AC[0];
   ACPC_zx[1] = 0.0; // projection onto the zx-plane
   ACPC_zx[2] = PC[2]-AC[2];

   // length of the ACPC_zx vector
   float ACPC_zx_length;

   ACPC_zx_length = sqrtf( ACPC_zx[0]*ACPC_zx[0] + ACPC_zx[2]*ACPC_zx[2] );

   // cosine and sine of the angle alpha that the vector ACPC_zx makes with the positive x-axis
   float cosalpha, sinalpha;

   cosalpha = ACPC_zx[0] / ACPC_zx_length;
   sinalpha = ACPC_zx[2] / ACPC_zx_length;

   rotate(R, cosalpha, sinalpha, 0.0, 1.0, 0.0);

   // Note: when rotation R is applied to the AC and PC vectors, their y-coordinates do not
   // change, since the rotation is about the y-axis.
   // The z-coordinate will be equal since the rotation makes the AC-PC vector horizontal 
   // in the zx-plate

   multi(R, 4, 4,  AC, 4,  1, AC);
   multi(R, 4, 4,  PC, 4,  1, PC);
   multi(R, 4, 4,  RP, 4,  1, RP);
   multi(R, 4, 4,  Tmsp, 4,  4, Tmsp);

/*
   float T[16];
   T[0]=1.0;  T[1]=0.0;  T[2]=0.0;  T[3]=0.0;
   T[4]=0.0;  T[5]=1.0;  T[6]=0.0;  T[7]=0.0;
   T[8]=0.0;  T[9]=0.0;  T[10]=1.0; T[11]=-AC[2]; // could also use PC[2] since AC[2]=PC[2]
   T[12]=0.0; T[13]=0.0; T[14]=0.0; T[15]=1.0;
   multi(T, 4, 4,  Tmsp, 4,  4, Tmsp);
*/
   // equivalent to the above
   Tmsp[11] -= AC[2]; // could also use PC[2] since AC[2]=PC[2]

   // AC and PC must also be updated by applying T to AC and PC
   // Applying T is equivalent to making their z-coordinates z
   AC[2] = 0.0;
   PC[2] = 0.0;

//   RP[2] -= AC[2]; // could subtract PC[2] since AC[2]=PC[2]
   RP[2] = 0.0; // wanted to ensure that RP is also on the z=0 plane

   ////////////////////////////////////////////////////////////////////////////////
   {
      float n[3]; // unit normal to the MSP in the original image xyz coordinate system
      float d;
      char filename[DEFAULT_STRING_LENGTH];
      FILE *fp;

      compute_MSP_parameters_from_Tmsp(Tmsp, n, &d);

      char dirname[DEFAULT_STRING_LENGTH]; // name of the directory only
      //char fullpath[DEFAULT_STRING_LENGTH]; // directory + filename
      char fullpath[2*DEFAULT_STRING_LENGTH+1]; // directory + filename
      getDirectoryName(imagefilename, dirname);
      niftiFilename(filename,imagefilename);
      snprintf(fullpath,sizeof(fullpath),"%s/%s_ACPC.txt",dirname, filename);

      if(opt_txt)
      {
         fp=fopen(fullpath,"w");

         if(fp==NULL)
         {
            errorMessage("Error: I have trouble opening the *_ACPC.txt file to save the output.");
         }

         fprintf(fp,"# Input volume: %s\n\n",imagefilename);
         fprintf(fp,"# Estimated mid-sagittal plane: (%8.7fx) + (%8.7fy) + (%8.7fz) = %8.5f (mm)\n", n[0],n[1],n[2],d);
         fprintf(fp,"%8.7f %8.7f %8.7f %8.5f\n\n", n[0],n[1],n[2],d);
      
         fclose(fp);
      }
   }
}

float detectPC(float *PC, char *modelfile, short *volumeMSP_HR, char *PCregion, short *xPC, short *yPC, short *zPC, int opt_T2)
{	
   model_file_hdr mhdr;
   DIM HR, LR;
   float *pc_template; // PC template, array of dimension PCtemplatesize
   float *mean_pc_template; // PC template, array of dimension PCtemplatesize
   float *pc_test_vec;
   float cc, ccmax;
   float *ccmap;
   float I2X[16]; // 4x4 matrix that transfomrs a vector from ijk-coordinates to xyz-coordinates
   //int PCint[3]; // variables storing PC ijk coordinates
   int PCint[3]={0,0,0}; // variables storing PC ijk coordinates
   int npHR, nvHR;
   int *x, *y, *z; 
   int n;
   FILE *fp;

   fp = fopen(modelfile,"r"); // open setup file for reading
   if( fp == NULL )
   {
      printf("\nCould not open %s, aborting ...\n\n",modelfile);
      exit(1);
   }

   if( fread(&mhdr, sizeof(mhdr), 1, fp) != 1 ) 
   {
      printf("Failed to read \"mhdr\" variable.\n");
   }

   if(bigEndian()) 
   {
      swap_model_file_hdr(&mhdr);
   }

   HR.nx=HR.ny=mhdr.nxHR;
   LR.nz=HR.nz=mhdr.nzHR;
   LR.dz=HR.dz=HR.dy=HR.dx=mhdr.dxHR; 
   LR.nx=LR.ny=mhdr.nxLR;
   LR.dy=LR.dx=2.0*HR.dx;

   npHR = HR.nx * HR.ny;
   nvHR = HR.nx * HR.ny * HR.nz;

   // MEMORY ALLOCATION
   pc_template = (float *)calloc(mhdr.nvol*mhdr.PCtemplatesize*mhdr.nangles,sizeof(float));
   mean_pc_template = (float *)calloc(mhdr.PCtemplatesize*mhdr.nangles,sizeof(float));
   pc_test_vec = (float *)calloc(mhdr.PCtemplatesize, sizeof(float));

   n=0;
   for(int v=0; v<nvHR; v++) if(PCregion[v]==1) n++;

   x = (int *)calloc(n, sizeof(int));
   y = (int *)calloc(n, sizeof(int));
   z = (int *)calloc(n, sizeof(int));

   n=0;
   for(int i=0; i<HR.nx; i++)
   for(int j=0; j<HR.ny; j++)
   for(int k=0; k<HR.nz; k++)
   if(PCregion[npHR*k + j*HR.nx + i]==1)
   {
      x[n]=i; y[n]=j; z[n]=k;
      n++;
   }

   // MEMORY ALLOCATION
   ccmap = (float *)calloc(nvHR, sizeof(float));

   // zero ccmap
   for(int i=0; i<nvHR; i++) ccmap[i]=0.0;

   for(int i=0; i<mhdr.nvol; i++)
   {
      fseek(fp, sizeof(float)*mhdr.RPtemplatesize*mhdr.nangles, SEEK_CUR);
      fseek(fp, sizeof(float)*mhdr.ACtemplatesize*mhdr.nangles, SEEK_CUR);

      if(fread(pc_template+i*mhdr.PCtemplatesize*mhdr.nangles, sizeof(float),mhdr.PCtemplatesize*mhdr.nangles, fp) != (size_t)(mhdr.PCtemplatesize*mhdr.nangles) )
      {
         printf("Failed to read \"pc_template\" variable.\n");
      }

      if(bigEndian())
      {
         for(int j=0; j<mhdr.PCtemplatesize*mhdr.nangles; j++)
            swapByteOrder( (char *)(pc_template+i*mhdr.PCtemplatesize*mhdr.nangles+j), sizeof(float));
      }
   }

   {
      int ns;
			
      ns = mhdr.PCtemplatesize*mhdr.nangles;

      for(int s=0; s<ns; s++)
      {
         mean_pc_template[s]=0.0;

         for(int i=0; i<mhdr.nvol; i++)
            mean_pc_template[s] += pc_template[ i*ns + s];

         mean_pc_template[s] /= mhdr.nvol;
      }
   }

	for(int v=0; v<n; v++)
	{
		extractArray(volumeMSP_HR,HR.nx,HR.ny,HR.nz,x[v],y[v],z[v],xPC,yPC,zPC,mhdr.PCtemplatesize, pc_test_vec);

		removeVectorMean( pc_test_vec, mhdr.PCtemplatesize );
		normalizeVector( pc_test_vec, mhdr.PCtemplatesize );

		ccmax=0.0;
		for(int j=0; j<mhdr.nangles; j++)
		{
			cc=dot(pc_test_vec, mean_pc_template+j*mhdr.PCtemplatesize, mhdr.PCtemplatesize);

            if(opt_T2) cc *= -1.0;

			if(cc>ccmax) ccmax=cc;
		}
		
		ccmap[npHR*z[v] + y[v]*HR.nx + x[v]] = ccmax;
	}

	ccmax=0.0;
	for(int v=0; v<n; v++)
	{
		cc = ccmap[npHR*z[v] + y[v]*HR.nx + x[v]];
		if(cc>ccmax) 
		{
			ccmax=cc;
			PCint[0]=x[v]; PCint[1]=y[v]; PCint[2]=z[v];
		}
	}

	PC[0]=PCint[0]; PC[1]=PCint[1]; PC[2]=PCint[2]; PC[3]=1.0;
	ijk2xyz(I2X, HR.nx, HR.ny, HR.nz, HR.dx, HR.dy, HR.dz);
   	multi(I2X,4,4, PC, 4, 1, PC);

	fclose(fp);

	delete ccmap;
	delete pc_template;
	delete mean_pc_template;
	delete pc_test_vec;
	delete x;
	delete y;
	delete z;

	return(ccmax);
}

float detectAC(float *AC, char *modelfile, short *volumeMSP_HR, char *ACregion, short *xAC, short *yAC, short *zAC, int opt_T2)
{
   model_file_hdr mhdr;
   FILE *fp;
   DIM HR, LR;
   int npHR, nvHR;
   int *x, *y, *z; 
   int n;
   //int ACint[3]; // variables storing AC ijk coordinates
   int ACint[3]={0,0,0}; // variables storing AC ijk coordinates
   float *ac_template; // AC template, array of dimension ACtemplatesize
   float *mean_ac_template; // AC template, array of dimension ACtemplatesize
   float *ac_test_vec;
   float I2X[16]; // 4x4 matrix that transfomrs a vector from ijk-coordinates to xyz-coordinates
   float cc, ccmax;
   float *ccmap;

   fp = fopen(modelfile,"r"); // open setup file for reading
   if( fp == NULL )
   {
      printf("\nCould not open %s, aborting ...\n\n",modelfile);
      exit(1);
   }

   if( fread(&mhdr, sizeof(mhdr), 1, fp) != 1 ) 
   {
      printf("Failed to read \"mhdr\" variable.\n");
   }

   if(bigEndian()) 
   {
      swap_model_file_hdr(&mhdr);
   }

   HR.nx=HR.ny=mhdr.nxHR;
   LR.nz=HR.nz=mhdr.nzHR;
   LR.dz=HR.dz=HR.dy=HR.dx=mhdr.dxHR; 
   LR.nx=LR.ny=mhdr.nxLR;
   LR.dy=LR.dx=2.0*HR.dx;

   npHR = HR.nx*HR.ny;
   nvHR = HR.nx*HR.ny*HR.nz;

   n=0;
   for(int v=0; v<nvHR; v++) if(ACregion[v]==1) n++;

   x = (int *)calloc(n, sizeof(int));
   y = (int *)calloc(n, sizeof(int));
   z = (int *)calloc(n, sizeof(int));

	n=0;
	for(int i=0; i<HR.nx; i++)
	for(int j=0; j<HR.ny; j++)
	for(int k=0; k<HR.nz; k++)
	if(ACregion[k*npHR + j*HR.nx + i]==1)
	{
		x[n]=i; y[n]=j; z[n]=k;
		n++;
	}

   // MEMORY ALLOCATION
   ccmap = (float *)calloc(nvHR, sizeof(float));

   for(int i=0; i<nvHR; i++) ccmap[i]=0.0; // zero ccmap


   // MEMORY ALLOCATION
   ac_template = (float *)calloc(mhdr.nvol*mhdr.ACtemplatesize*mhdr.nangles,sizeof(float));
   mean_ac_template = (float *)calloc(mhdr.ACtemplatesize*mhdr.nangles,sizeof(float));
   ac_test_vec = (float *)calloc(mhdr.ACtemplatesize, sizeof(float));

   //////////////////////////////////////////////////////////////////////////////////////
   for(int i=0; i<mhdr.nvol; i++)
   {
      fseek(fp, sizeof(float)*mhdr.RPtemplatesize*mhdr.nangles, SEEK_CUR);

      if( fread(ac_template+i*mhdr.ACtemplatesize*mhdr.nangles, sizeof(float), mhdr.ACtemplatesize*mhdr.nangles, fp) != (size_t)(mhdr.ACtemplatesize*mhdr.nangles) ) 
      {
         printf("Failed to read \"ac_template\" variable.\n");
      }

      fseek(fp, sizeof(float)*mhdr.PCtemplatesize*mhdr.nangles, SEEK_CUR);

      if(bigEndian())
      {
         for(int j=0; j<mhdr.ACtemplatesize*mhdr.nangles; j++)
            swapByteOrder( (char *)(ac_template+i*mhdr.ACtemplatesize*mhdr.nangles+j), sizeof(float));
      }
   }

	{
		int ns;
			
		ns = mhdr.ACtemplatesize*mhdr.nangles;

		for(int s=0; s<ns; s++)
		{
			mean_ac_template[s]=0.0;

			for(int i=0; i<mhdr.nvol; i++)
				mean_ac_template[s] += ac_template[ i*ns + s];

			mean_ac_template[s] /= mhdr.nvol;
		}
	}

	for(int v=0; v<n; v++)
	{
		extractArray(volumeMSP_HR,HR.nx,HR.ny,HR.nz,x[v],y[v],z[v],xAC,yAC,zAC,mhdr.ACtemplatesize,ac_test_vec);

		removeVectorMean( ac_test_vec, mhdr.ACtemplatesize);
		normalizeVector(ac_test_vec, mhdr.ACtemplatesize);

		ccmax=0.0;
		for(int j=0; j<mhdr.nangles; j++)
		{
			cc=dot(ac_test_vec, mean_ac_template+j*mhdr.ACtemplatesize, mhdr.ACtemplatesize);

            if(opt_T2) cc *= -1.0;
			if(cc>ccmax) ccmax=cc;
		}
		
		ccmap[z[v]*npHR + y[v]*HR.nx + x[v]] = ccmax;
	}

	ccmax=0.0;
	for(int v=0; v<n; v++)
	{
		cc = ccmap[z[v]*npHR + y[v]*HR.nx + x[v]];
		if(cc>ccmax) 
		{
			ccmax=cc;
			ACint[0]=x[v]; ACint[1]=y[v]; ACint[2]=z[v];
		}
	}
	fclose(fp);

	AC[0]=ACint[0]; AC[1]=ACint[1]; AC[2]=ACint[2]; AC[3]=1.0;
	ijk2xyz(I2X, HR.nx, HR.ny, HR.nz, HR.dx, HR.dy, HR.dz);
   	multi(I2X,4,4, AC, 4, 1, AC);

	delete ac_template;
	delete mean_ac_template;
	delete ac_test_vec;
	delete ccmap;
	delete x;
	delete y;
	delete z;

	return(ccmax);
}

void initialAC(float Ax, float Ay, float Bx, float By, float *Cx, float *Cy, float parcomMean, float percomMean)
{
	double AB;

	AB = sqrt( (double) ((Ax-Bx)*(Ax-Bx) + (Ay-By)*(Ay-By)) );

	*Cx = (float)( Ax + parcomMean*(Bx-Ax)/AB + percomMean*(By-Ay)/AB); 
	*Cy = (float)( Ay + parcomMean*(By-Ay)/AB - percomMean*(Bx-Ax)/AB);
}

char *defineACregion(DIM dim, float *RP, float *PC, float parcomMean, float percomMean, double ACsr) 
{
   double r;
   char *ACregion;
   int np;
   float x,y;
   float ACinitial_guess[4];

   np = dim.nx*dim.ny;

   initialAC(RP[0], RP[1], PC[0], PC[1], &ACinitial_guess[0], &ACinitial_guess[1], parcomMean, percomMean);

   ACregion = (char *)calloc(dim.nx*dim.ny*dim.nz, sizeof(char));

   for(int i=0; i<dim.nx; i++)
   for(int j=0; j<dim.ny; j++)
   for(int k = (dim.nz-1)/2 - DEL; k <= (dim.nz-1)/2 + DEL; k++)
   {
      x = ACinitial_guess[0] + dim.dx*(dim.nx-1)/2.0;
      y = ACinitial_guess[1] + dim.dy*(dim.ny-1)/2.0;

      r = (i*dim.dx-x)*(i*dim.dx-x) + (j*dim.dy-y)*(j*dim.dy-y);

      if( r <= ACsr*ACsr ) ACregion[np*k + dim.nx*j + i]=1;
      else ACregion[np*k + dim.nx*j + i]=0;
   }

   return(ACregion);
}

char *definePCregion(DIM HR, float *RP, float *RPPCmean, double PCsr) 
{
   char *PCregion;
   float rp[4];
   float x, y; 
   float X2I[16]; // 4x4 matrix that transforms a vector from xyz-coordinates to ijk-coordinates
   double r;
   int np;

   xyz2ijk(X2I, HR);

   multi(X2I,4,4,RP,4,1,rp);

   np = HR.nx*HR.ny;

   PCregion = (char *)calloc(HR.nx*HR.ny*HR.nz, sizeof(char));

   for(int i=0; i<HR.nx; i++)
   for(int j=0; j<HR.ny; j++)
   for(int k = (HR.nz-1)/2 - DEL; k <= (HR.nz-1)/2 + DEL; k++)
   {
      x = rp[0]*HR.dx + RPPCmean[0];
      y = rp[1]*HR.dy + RPPCmean[1];

      r = (i*HR.dx-x)*(i*HR.dx-x) + (j*HR.dy-y)*(j*HR.dy-y);

      if( r <= PCsr*PCsr ) PCregion[np*k + HR.nx*j + i]=1;
      else PCregion[np*k + HR.nx*j + i]=0;
   }

   return(PCregion);
}

void detectRP(float *RP1, float *RP2, char *modelfile, short *volumeMSP_LR, short *mask_LR, short *xRP, 
short *yRP, short *zRP, int opt_T2)
{
   FILE *fp;
   model_file_hdr mhdr;
   DIM LR, HR;
   long templatesFilePosition;
   float *ccmap_LR;
   //int RPint[3]; // variables storing RP ijk coordinates
   int RPint[3]={0,0,0}; // variables storing RP ijk coordinates
   int npLR;
   float *rp_template;
   float *mean_rp_template;
   float *rp_test_vec;
   float cc, ccmax;
   float I2X[16];

   fp = fopen(modelfile,"r"); // open setup file for reading
   if( fp == NULL )
   {
      printf("\nCould not open %s, aborting ...\n\n",modelfile);
      exit(1);
   }

   if( fread(&mhdr, sizeof(mhdr), 1, fp) != 1 ) 
   {
      printf("Failed to read \"mhdr\" variable.\n");
   }

   if(bigEndian()) 
   {
      swap_model_file_hdr(&mhdr);
   }

   HR.nx=HR.ny=mhdr.nxHR;
   LR.nz=HR.nz=mhdr.nzHR;
   LR.dz=HR.dz=HR.dy=HR.dx=mhdr.dxHR; 
   LR.nx=LR.ny=mhdr.nxLR;
   LR.dy=LR.dx=2.0*HR.dx;

   npLR = LR.nx*LR.ny;

   // MEMORY ALLOCATION
   rp_template = (float *)calloc(mhdr.RPtemplatesize*mhdr.nangles*mhdr.nvol,sizeof(float));
   mean_rp_template = (float *)calloc(mhdr.RPtemplatesize*mhdr.nangles,sizeof(float));
   rp_test_vec = (float *)calloc(mhdr.RPtemplatesize, sizeof(float));

   templatesFilePosition=ftell(fp);

   fseek(fp, templatesFilePosition, SEEK_SET);

   // MEMORY ALLOCATION
   ccmap_LR = (float *)calloc(LR.nx*LR.ny, sizeof(float));
   RPint[2] = (LR.nz-1)/2;

   for(int i=0; i<mhdr.nvol; i++)
   {

      if( fread(rp_template+i*mhdr.RPtemplatesize*mhdr.nangles, sizeof(float), mhdr.RPtemplatesize*mhdr.nangles, fp) != (size_t)(mhdr.RPtemplatesize*mhdr.nangles) )
      {
         printf("Failed to read \"rp_template\" variable.\n");
      }

      fseek(fp, sizeof(float)*mhdr.ACtemplatesize*mhdr.nangles, SEEK_CUR);
      fseek(fp, sizeof(float)*mhdr.PCtemplatesize*mhdr.nangles, SEEK_CUR);

      if(bigEndian())
      {
         for(int j=0; j<mhdr.RPtemplatesize*mhdr.nangles; j++)
            swapByteOrder( (char *)(rp_template+i*mhdr.RPtemplatesize*mhdr.nangles+j), sizeof(float));
      }
   }

   {
		int ns;
			
		ns = mhdr.RPtemplatesize*mhdr.nangles;

		for(int s=0; s<ns; s++)
		{
			mean_rp_template[s]=0.0;

			for(int i=0; i<mhdr.nvol; i++)
				mean_rp_template[s] += rp_template[ i*ns + s];

			mean_rp_template[s] /= mhdr.nvol;
		}
   }

	for(int x=mhdr.RPtemplateradius; x<LR.nx-mhdr.RPtemplateradius; x++)
	for(int y=mhdr.RPtemplateradius; y<LR.ny-mhdr.RPtemplateradius; y++)
	if(mask_LR[ RPint[2]*npLR + y*LR.nx + x])
	{
		extractArray(volumeMSP_LR, LR.nx, LR.ny, LR.nz, x, y, RPint[2], xRP,yRP,zRP,mhdr.RPtemplatesize, rp_test_vec);

		removeVectorMean( rp_test_vec, mhdr.RPtemplatesize);
		normalizeVector(rp_test_vec, mhdr.RPtemplatesize);

		ccmax=0.0;
		for(int j=0; j<mhdr.nangles; j++)
		{
			cc=dot(rp_test_vec, mean_rp_template + j*mhdr.RPtemplatesize, mhdr.RPtemplatesize);

            if(opt_T2) cc *= -1.0;

			if(cc>ccmax) ccmax=cc;
		}
		
		ccmap_LR[y*LR.nx + x] = ccmax;
	}

	ccmax=0.0;
	for(int x=mhdr.RPtemplateradius; x<LR.nx-mhdr.RPtemplateradius; x++)
	for(int y=mhdr.RPtemplateradius; y<LR.ny-mhdr.RPtemplateradius; y++)
	{
		cc = ccmap_LR[y*LR.nx + x];
		if(cc>ccmax) 
		{
			ccmax=cc;
			RPint[0]=x; RPint[1]=y;
		}
	}

	RP1[0]=RPint[0]; RP1[1]=RPint[1]; RP1[2]=RPint[2]; RP1[3]=1.0;
	ijk2xyz(I2X, LR.nx, LR.ny, LR.nz, LR.dx, LR.dy, LR.dz);
   	multi(I2X,4,4, RP1, 4, 1, RP1);

	for(int x=RPint[0]-mhdr.RPtemplateradius; x<RPint[0]+mhdr.RPtemplateradius; x++)
	for(int y=RPint[1]-mhdr.RPtemplateradius; y<RPint[1]+mhdr.RPtemplateradius; y++)
	{
		if( (y*LR.nx + x)>=0 && (y*LR.nx + x)<npLR)
			ccmap_LR[y*LR.nx + x]=0.0;
	}

	ccmax=0.0;
	for(int x=mhdr.RPtemplateradius; x<LR.nx-mhdr.RPtemplateradius; x++)
	for(int y=mhdr.RPtemplateradius; y<LR.ny-mhdr.RPtemplateradius; y++)
	{
		cc = ccmap_LR[y*LR.nx + x];
		if(cc>ccmax) 
		{
			ccmax=cc;
			RPint[0]=x; RPint[1]=y;
		}
	}

	RP2[0]=RPint[0]; RP2[1]=RPint[1]; RP2[2]=RPint[2]; RP2[3]=1.0;
	ijk2xyz(I2X, LR.nx, LR.ny, LR.nz, LR.dx, LR.dy, LR.dz);
   	multi(I2X,4,4, RP2, 4, 1, RP2);

	fclose(fp);

	delete ccmap_LR;
	delete rp_template;
	delete mean_rp_template;
	delete rp_test_vec;
}

void defineTemplate(int r, int h, short *x, short *y, short *z)
{
   int c=0;
   int r2;
   int h_2;

   h_2 = (h-1)/2;
   r2=r*r;

   for(int i=-r; i<=r; i++)
   for(int j=-r; j<=r; j++)
   for(int k=-h_2; k<=h_2; k++)
   {
      if( (i*i+j*j)  <= r2 )
      {
         x[c] = i;
         y[c] = j;
         z[c] = k;
         c++;
      }
   }
}

char *expandMask(short *mask_HR, DIM HR, float *RPmean, double RPsr)
{
   char *RPregion;
   double r1, r2, r;
   int npHR;

   npHR = HR.nx * HR.ny;

   RPregion = (char *)calloc(HR.nx*HR.ny, sizeof(char));

   for(int i=0; i<HR.nx; i++)
   for(int j=0; j<HR.ny; j++)
   for(int k=0; k<HR.nz; k++)
   {
      r1 = (i - (HR.nx-1.0)/2.0)*HR.dx - RPmean[0];
      r1 *= r1;
      r2 = (j - (HR.ny-1.0)/2.0)*HR.dy - RPmean[1];
      r2 *= r2;
      r = sqrt(r1 + r2);

      if(r >= RPsr) 
      {	
         mask_HR[k*npHR + j*HR.nx + i]=0;
         RPregion[j*HR.nx + i]=0;
      }
      else
      {	
         RPregion[j*HR.nx + i]=1;
      }
   }

   return(RPregion);
}


// Thresholds the input image `im' using Otsu's method of automatic threshold selection.
// The number of suprathreshold voxels is returned in `nbv'.
// The function returns and image `msk' where suprathreshold voxels have a value of 1, and
// subthreshold voxels have a value of 0.
short *thresholdImageOtsu(short *im, int nv, int *nbv)
{
	short *msk;
	int low, high, bw, k, thresh;
	double *h=NULL;

	msk = (short *)calloc(nv, sizeof(short));
	if(msk==NULL) return(NULL);

	setLowHigh(im, nv, &low, &high);

	h=findHistogram(im, nv, 256, low, high, &bw);
	if(h==NULL) { free(msk); return(NULL); }

	k=otsu(h, 256);
	if(k==-1) { free(msk); free(h); return(NULL); }

	thresh = low + k*bw;

	*nbv = 0;
	for(int i=0; i<nv; i++) 
	if(im[i]<=thresh) 
		msk[i]=0; 
	else 
	{ 
		msk[i]=1; 
		(*nbv)++; 
	}

	free(h);

	return(msk);
}

int detect_AC_PC_MSP(const char *imagefilename, char *orientation, char *modelfile,
float *AC, float *PC, float *RP, float *Tmsp, int opt_v, int opt_T2)
{
   char modelfilepath[DEFAULT_STRING_LENGTH];

   // (x, y, z) image orientation vectors (row, column, and slice orientation)
   float xvec[3], yvec[3], zvec[3];
   float msp_normal_xyz[3]; // unit normal to the MSP in xyz coordinates
   float msp_normal_RAS[3]; // unit normal to the MSP in NIFTI RAS system
   float msp_normal_LPS[3]; // unit normal to the MSP in DICOM LPS system
   float ACPCvec_xyz[3]; // vector from AC to PC in xyz coordinates
   float ACPCvec_RAS[3]; // vector from AC to PC in NIFTI RAS system
   float ACPCvec_LPS[3]; // vector from AC to PC in DICOM LPS system
   float d;

   model_file_hdr mhdr;
   model_file_tail mtail;
   DIM HR; 
   DIM LR; 
   DIM Orig;
   short *volOrig; // original input volume from the training set
   short *mask_HR;
   short *volumeMSP_HR; // MSP aligned volume
   short *volumeMSP_LR; // MSP aligned volume
   char *PCregion1, *PCregion2, *PCregion;
   char *ACregion1, *ACregion2, *ACregion;
   char *RPregion;
   float *invT; 
   int nbv;

   if(ARTHOME==NULL) getARTHOME();

   // ensure that the user has specified an image
   if(imagefilename[0]=='\0')
   {
      errorMessage("No input image filename in detect_AC_PC_MSP().");
   }

   // ensure that the specified image has either a .hdr or a .nii extension 
   if( !checkNiftiFileExtension((const char *)imagefilename) )
   {
      errorMessage("The image filename in detect_AC_PC_MSP() must have a `.hdr' or `.nii' extension.");
   }

   if(orientation[0]=='\0')
   {
      getNiftiImageOrientation(imagefilename, orientation);
   }

   if(orientation[0]=='\0')
   {
      errorMessage("Image orientation cannot be determined in detect_AC_PC_MSP().");
   }

   if ( isOrientationCodeValid(orientation) == 0)
   {
      printf("\nInput image orientation: %s\n",orientation);
      errorMessage("Invalid orientation code in detect_AC_PC_MSP(). The code is not one of the 48 legal ones.");
   }

   // If a specific model is not specified, 
   // then use the standard model (T1acpc.mdl) in ARTHOME directory. 

   if(modelfile[0]=='\0')
   {
      snprintf(modelfilepath,sizeof(modelfilepath),"%s/T1acpc.mdl",ARTHOME);
   }
   else
   {
      snprintf(modelfilepath,sizeof(modelfilepath),"%s/%s",ARTHOME,modelfile);
   }

   {
      // read information from the model file and initialize some variables

      // file pointer for opening the model file
      FILE *fp; 

      fp = fopen(modelfilepath,"r"); // open setup file for reading
      if( fp == NULL )
      {
         printf("\nI cannot open the model file: %s.\n\n",modelfilepath);
         exit(1);
      }

      if( fread(&mhdr, sizeof(mhdr), 1, fp) != 1 ) 
      {
         printf("Failed to read \"mhdr\" variable.\n");
      }

      if(bigEndian()) 
      {
         swap_model_file_hdr(&mhdr);
      }

      HR.nx=HR.ny=mhdr.nxHR;
      LR.nz=HR.nz=mhdr.nzHR;
      LR.dz=HR.dz=HR.dy=HR.dx=mhdr.dxHR; 
      LR.nx=LR.ny=mhdr.nxLR;
      LR.dy=LR.dx=2.0*HR.dx;

      fseek(fp, sizeof(float)*mhdr.RPtemplatesize*mhdr.nangles*mhdr.nvol, SEEK_CUR);
      fseek(fp, sizeof(float)*mhdr.ACtemplatesize*mhdr.nangles*mhdr.nvol, SEEK_CUR);
      fseek(fp, sizeof(float)*mhdr.PCtemplatesize*mhdr.nangles*mhdr.nvol, SEEK_CUR);

      if( fread(&mtail, sizeof(mtail), 1, fp) != 1 ) 
      {
         printf("Failed to read \"mtail\" variable.\n");
      }

      if(bigEndian()) 
      {
         swap_model_file_tail(&mtail);
      }

      fclose(fp);
   }

   ////////////////////////////////////////////////////////////////////////////

   volOrig = readNiftiImage( (const char *)imagefilename, &Orig, 0);

   if(volOrig == NULL)
   {
      errorMessage("I could not read the input image in detect_AC_PC_MSP().");
   }

   if(opt_MSP)
   {
      computeTmsp(orientation, volOrig, Orig, Tmsp);
   }

   if(!opt_AC)
   {
      float I2X[16];

      ijk2xyz(I2X, Orig.nx, Orig.ny, Orig.nz, Orig.dx, Orig.dy, Orig.dz);
      multi(I2X, 4, 4,  AC, 4,  1, AC);
      multi(Tmsp, 4, 4,  AC, 4,  1, AC);
   }

   if(!opt_PC)
   {
      float I2X[16];

      ijk2xyz(I2X, Orig.nx, Orig.ny, Orig.nz, Orig.dx, Orig.dy, Orig.dz);
      multi(I2X, 4, 4,  PC, 4,  1, PC);
      multi(Tmsp, 4, 4,  PC, 4,  1, PC);
   }

   if(!opt_RP)
   {
      float I2X[16];

      ijk2xyz(I2X, Orig.nx, Orig.ny, Orig.nz, Orig.dx, Orig.dy, Orig.dz);
      multi(I2X, 4, 4,  RP, 4,  1, RP);
      multi(Tmsp, 4, 4,  RP, 4,  1, RP);
   }

   ////////////////////////////////////////////////////////////////////////////

   invT = inv4(Tmsp);
   volumeMSP_HR=resliceImage(volOrig,Orig,HR,invT,LIN);
   free(invT);
   volumeMSP_LR = resizeXYZ(volumeMSP_HR, HR, LR);

   {
      short *maskOrig;
      maskOrig=thresholdImageOtsu(volOrig, Orig.nx*Orig.ny*Orig.nz, &nbv);

      invT = inv4(Tmsp);
      mask_HR = resliceImage(maskOrig,Orig,HR,invT,LIN);
      free(invT);

      delete maskOrig;
   }

   RPregion = expandMask(mask_HR, HR, mtail.RPmean, searchradius[0]);

   /////////////////////////////////////////////////////////////////////////////////////////////////

   {
      float ACccmax[2];
      float PCccmax[2];
      float AC1[4], AC2[4]; // variables storing AC coordinates
      float PC1[4], PC2[4]; // variables storing PC coordinates
      float RP1[4], RP2[4]; // variables storing RP coordinates
      short *xRP, *yRP, *zRP;
      short *xAC, *yAC, *zAC;
      short *xPC, *yPC, *zPC;
      short *mask_LR;

      xRP = (short *)calloc(mhdr.RPtemplatesize, sizeof(short));
      yRP = (short *)calloc(mhdr.RPtemplatesize, sizeof(short));
      zRP = (short *)calloc(mhdr.RPtemplatesize, sizeof(short));
      defineTemplate(mhdr.RPtemplateradius,mhdr.RPtemplateheight, xRP, yRP, zRP);

      xAC = (short *)calloc(mhdr.ACtemplatesize, sizeof(short));
      yAC = (short *)calloc(mhdr.ACtemplatesize, sizeof(short));
      zAC = (short *)calloc(mhdr.ACtemplatesize, sizeof(short));
      defineTemplate(mhdr.ACtemplateradius,mhdr.ACtemplateheight, xAC, yAC, zAC);

      xPC = (short *)calloc(mhdr.PCtemplatesize, sizeof(short));
      yPC = (short *)calloc(mhdr.PCtemplatesize, sizeof(short));
      zPC = (short *)calloc(mhdr.PCtemplatesize, sizeof(short));
      defineTemplate(mhdr.PCtemplateradius,mhdr.PCtemplateheight, xPC, yPC, zPC);

      mask_LR = resizeXYZ(mask_HR, HR, LR);
      if(opt_RP)
      {
         detectRP(RP1, RP2, modelfilepath, volumeMSP_LR, mask_LR, xRP, yRP, zRP, opt_T2);
      }
      else
      {
         for(int i=0; i<4; i++) RP1[i]=RP2[i]=RP[i];
      }

      PCregion1=definePCregion(HR, RP1, mtail.RPPCmean, searchradius[2]); 
      if(opt_PC)
      {
         PCccmax[0] = detectPC(PC1, modelfilepath, volumeMSP_HR, PCregion1, xPC, yPC, zPC, opt_T2);
      }
      else
      {
         for(int i=0; i<4; i++) PC1[i]=PC[i];
         PCccmax[0]=0.0;
      }

      ACregion1=defineACregion(HR, RP1, PC1, mtail.parcomMean, mtail.percomMean, searchradius[1]);
      if(opt_AC)
      {
         ACccmax[0] = detectAC(AC1, modelfilepath, volumeMSP_HR, ACregion1, xAC, yAC, zAC, opt_T2);
      }
      else
      {
         for(int i=0; i<4; i++) AC1[i]=AC[i];
         ACccmax[0]=0.0;
      }

      PCregion2=definePCregion(HR, RP2, mtail.RPPCmean, searchradius[2]); 
      if(opt_PC)
      {
         PCccmax[1] = detectPC(PC2, modelfilepath, volumeMSP_HR, PCregion2, xPC, yPC, zPC, opt_T2);
      }
      else
      {
         for(int i=0; i<4; i++) PC2[i]=PC[i];
         PCccmax[1]=0.0;
      }

      ACregion2=defineACregion(HR, RP2, PC2, mtail.parcomMean, mtail.percomMean, searchradius[1]);
      if(opt_AC)
      {
         ACccmax[1] = detectAC(AC2, modelfilepath, volumeMSP_HR, ACregion2, xAC, yAC, zAC, opt_T2);
      }
      else
      {
         for(int i=0; i<4; i++) AC2[i]=AC[i];
         ACccmax[1]=0.0;
      }

      if((ACccmax[0]+PCccmax[0])>(ACccmax[1]+PCccmax[1])) 
      { 
         PC[0]=PC1[0]; PC[1]=PC1[1]; PC[2]=PC1[2]; PC[3]=PC1[3]; 
         AC[0]=AC1[0]; AC[1]=AC1[1]; AC[2]=AC1[2]; AC[3]=AC1[3]; 
         RP[0]=RP1[0]; RP[1]=RP1[1]; RP[2]=RP1[2]; RP[3]=RP1[3]; 
         ACregion=ACregion1; PCregion=PCregion1;
      }
      else 
      { 
         PC[0]=PC2[0]; PC[1]=PC2[1]; PC[2]=PC2[2]; PC[3]=PC2[3]; 
         AC[0]=AC2[0]; AC[1]=AC2[1]; AC[2]=AC2[2]; AC[3]=AC2[3]; 
         RP[0]=RP2[0]; RP[1]=RP2[1]; RP[2]=RP2[2]; RP[3]=RP2[3]; 
         ACregion=ACregion2; PCregion=PCregion2;
      }

      delete xRP; delete yRP; delete zRP;
      delete xAC; delete yAC; delete zAC;
      delete xPC; delete yPC; delete zPC;
      delete mask_LR;
   }

   /////////////////////////////////////////////////////////////////////////////////////////////////
	
   // save the results 
   updateTmsp(imagefilename, Tmsp, RP, AC, PC);

   saveACPClocation(imagefilename, Tmsp, Orig, AC, PC, RP, opt_v);

   saveACPCimages(imagefilename, ACregion, PCregion, RPregion, AC, PC, RP, HR, Orig, volOrig, Tmsp);

   {
      // the code in this block changes the AC and PC cooridnates from the (x,y,z) system in MSP aligned PIL
      // orientation to (x,y,z) system in the original image orientation

      float *invT; 

      invT = inv4(Tmsp);

      multi(invT, 4, 4,  AC, 4,  1, AC);
      multi(invT, 4, 4,  PC, 4,  1, PC);
      multi(invT, 4, 4,  RP, 4,  1, RP);

      free(invT);
   }

   readOrientationVectorsFromFile(imagefilename, xvec, yvec, zvec);
   //printf("\nxvec = %f %f %f\n",xvec[0], xvec[1], xvec[2]);
   //printf("\nyvec = %f %f %f\n",yvec[0], yvec[1], yvec[2]);
   //printf("\nzvec = %f %f %f\n",zvec[0], zvec[1], zvec[2]);

   ACPCvec_xyz[0] = PC[0] - AC[0];
   ACPCvec_xyz[1] = PC[1] - AC[1];
   ACPCvec_xyz[2] = PC[2] - AC[2];

   ACPCvec_RAS[0] = ACPCvec_xyz[0] * xvec[0]  +  ACPCvec_xyz[1] * yvec[0] + ACPCvec_xyz[2] * zvec[0];
   ACPCvec_RAS[1] = ACPCvec_xyz[0] * xvec[1]  +  ACPCvec_xyz[1] * yvec[1] + ACPCvec_xyz[2] * zvec[1];
   ACPCvec_RAS[2] = ACPCvec_xyz[0] * xvec[2]  +  ACPCvec_xyz[1] * yvec[2] + ACPCvec_xyz[2] * zvec[2];

   ACPCvec_LPS[0] = -ACPCvec_RAS[0];
   ACPCvec_LPS[1] = -ACPCvec_RAS[1];
   ACPCvec_LPS[2] = ACPCvec_RAS[2];

   normalizeVector(ACPCvec_LPS,3);
   //printf("\nACPCvec_LPS = %f %f %f\n",ACPCvec_LPS[0], ACPCvec_LPS[1], ACPCvec_LPS[2]);
   //printf("\n%f\n",  ACPCvec_LPS[0]*ACPCvec_LPS[0] + ACPCvec_LPS[1]*ACPCvec_LPS[1] + ACPCvec_LPS[2]*ACPCvec_LPS[2] );

   compute_MSP_parameters_from_Tmsp(Tmsp, msp_normal_xyz, &d);

   msp_normal_RAS[0] = msp_normal_xyz[0] * xvec[0]  +  msp_normal_xyz[1] * yvec[0] + msp_normal_xyz[2] * zvec[0];
   msp_normal_RAS[1] = msp_normal_xyz[0] * xvec[1]  +  msp_normal_xyz[1] * yvec[1] + msp_normal_xyz[2] * zvec[1];
   msp_normal_RAS[2] = msp_normal_xyz[0] * xvec[2]  +  msp_normal_xyz[1] * yvec[2] + msp_normal_xyz[2] * zvec[2];

   msp_normal_LPS[0] = -msp_normal_RAS[0];
   msp_normal_LPS[1] = -msp_normal_RAS[1];
   msp_normal_LPS[2] = msp_normal_RAS[2];

   //printf("\nmsp_normal_LPS = %f %f %f\n",msp_normal_LPS[0], msp_normal_LPS[1], msp_normal_LPS[2]);

   {
      float angle[3];

      backwardTCSAP(msp_normal_LPS, ACPCvec_LPS, angle);

/*
      if(opt_v)
      {
         printf("\n**********************************************************************\n");
         printf("Trio FOV parameters for axial AC-PC aligned images (T>C>S A>>P mode):\n");
         printf("angles (deg.): T>C=%0.1f  >S=%0.1f  A>>P=%0.1f\n", angle[0], angle[1], angle[2]);
         printf("center (mm): L=%0.1f  P=%0.1f  H=%0.1f\n", AC[0], AC[1], AC[2]);
         printf("**********************************************************************\n");
      }
*/
   }


   {
      // the code in this block changes the AC and PC cooridnates from the (x,y,z) system
      // to (i,j,k) system in the original image orientation

      float X2I[16]; // 4x4 matrix that transforms a vector from xyz-coordinates to ijk-coordinates

      xyz2ijk(X2I, Orig);

      multi(X2I, 4, 4,  AC, 4,  1, AC);
      multi(X2I, 4, 4,  PC, 4,  1, PC);
      multi(X2I, 4, 4,  RP, 4,  1, RP);
   }

// The following block is just for testing to show that AC, PC, and RP are on the MSP
/*
{
   float ac[4], pc[4], rp[4];
   float I2X[16];
   ijk2xyz(I2X, Orig.nx, Orig.ny, Orig.nz, Orig.dx, Orig.dy, Orig.dz);
   float v[4]; // ac to rp vector
   float u[4]; // ac to pc vector
   float n[4]; // n = vxu / |vxu|

   multi(I2X, 4, 4,  RP, 4,  1, rp);
   multi(I2X, 4, 4,  AC, 4,  1, ac);
   multi(I2X, 4, 4,  PC, 4,  1, pc);

   for(int i=0; i<3; i++)
   {
      v[i] = rp[i] - ac[i];
      u[i] = pc[i] - ac[i];
   }

   n[0] = v[2]*u[1]-v[1]*u[2]; 
   n[1] = v[0]*u[2]-v[2]*u[0]; 
   n[2] = v[1]*u[0]-v[0]*u[1];

   normalizeVector(n,3);

   printf("\n%f %f %f %f\n", n[0], n[1], n[2], n[3]);
   printf("\n%f\n", n[0]*rp[0] + n[1]*rp[1] + n[2]*rp[2]);
   printf("\n%f\n", n[0]*ac[0] + n[1]*ac[1] + n[2]*ac[2]);
   printf("\n%f\n", n[0]*pc[0] + n[1]*pc[1] + n[2]*pc[2]);
}
*/


   //////////MEMORY RELEASE ////////////////////////////////////////////////////////////////////////

   delete ACregion1;
   delete ACregion2;
   delete PCregion1;
   delete PCregion2;
   delete RPregion;

   delete volOrig;

   delete volumeMSP_LR;
   delete volumeMSP_HR;

   delete mask_HR;

   return(0);
}

float reflectVertex(int pmax, float fac)
{
	//int j;
	float fac1,fac2;
	float newV;

	fac1=(1.0-fac)/3;
	fac2=fac1-fac;

	VertexNew[0]=vertexSum[0]*fac1-vertex[pmax][0]*fac2;
	VertexNew[1]=vertexSum[1]*fac1-vertex[pmax][1]*fac2;
	VertexNew[2]=vertexSum[2]*fac1-vertex[pmax][2]*fac2;

	newV=-reflection_cross_correlation2(gimage,gdim,VertexNew[0],VertexNew[1],
	VertexNew[2]);

	if (newV < value[pmax]) 
	{
		value[pmax]=newV;

		vertexSum[0] += VertexNew[0]-vertex[pmax][0]; 
		vertex[pmax][0]=VertexNew[0];

		vertexSum[1] += VertexNew[1]-vertex[pmax][1]; 
		vertex[pmax][1]=VertexNew[1];

		vertexSum[2] += VertexNew[2]-vertex[pmax][2]; 
		vertex[pmax][2]=VertexNew[2];
	}

	return(newV);
}

// Downhill Simplex Method of optimization
int dsm(void)
{
   int iter;

   int i,pmax,pmin,pnmax,j;
   float tol;
   float sum,ysave;
   float newV;

   for (j=0;j<3;j++) 
   {
      sum=0.0;
      for (i=0;i<4;i++) sum += vertex[i][j];
      vertexSum[j]=sum;
   }

   iter=1;
   while( iter < MAXITER )
   {
      iter++;

		if (value[0]>value[1] )
		{
			pmax=0; pnmax=1;
		}
		else
		{
			pmax=1; pnmax=0;
		} 

		for (i=0;i<4;i++) 
		if (value[i] > value[pmax]) 
		{
			pnmax=pmax;
			pmax=i;
		} 
		else if(value[i] > value[pnmax] && i != pmax)
		{
		  	pnmax=i;
		}

		pmin=0;
		if (value[1] < value[pmin]) pmin=1;
		if (value[2] < value[pmin]) pmin=2;
		if (value[3] < value[pmin]) pmin=3;

		tol=2.0*fabs(value[pmax]-value[pmin])/
		(fabs(value[pmax])+fabs(value[pmin]));

		if (tol < 1.0e-6 ) 
			break;

		newV=reflectVertex(pmax,-1.0);

		if (newV <= value[pmin])
		{
			newV=reflectVertex(pmax,2.0);
		}
		else if (newV >= value[pnmax]) 
		{
			ysave=value[pmax];
			newV=reflectVertex(pmax,0.5);
			if (newV >= ysave) 
			{
				for (i=0;i<4;i++) 
				{
					if (i != pmin) 
					{

						vertex[i][0]=vertexSum[0]=0.5*(vertex[i][0]+vertex[pmin][0]);
						vertex[i][1]=vertexSum[1]=0.5*(vertex[i][1]+vertex[pmin][1]);
						vertex[i][2]=vertexSum[2]=0.5*(vertex[i][2]+vertex[pmin][2]);

						value[i]=-reflection_cross_correlation2(gimage,gdim,
						vertexSum[0],vertexSum[1],vertexSum[2]);
					}
				}

   			for (j=0;j<3;j++) 
				{
					sum=0.0;
   				for (i=0;i<4;i++) sum += vertex[i][j];
   				vertexSum[j]=sum;
				}
			}
		} 
	}

	return(pmin);
}

/* Denote the input image by Q. If every point in Q is reflected with respect 
to the plane a*x+b*y+c*z=d,  a new image, R, will be obtained. This function
returns the cross-correlation between Q and R 
given by: C(Q,R) = Q.R/sqrt(Q.Q*R.R).  (a,b,c) is a unit vector perpendicular 
to the plane. d is the perpendicular distance between the origin and the
plane  */

float reflection_cross_correlation2(short *image, DIM dim, float A, float B, float C)
{
   float a,b,c,d;
   float dum;

   dum = A*A+B*B+C*C;

   if(dum<0.0)
   {
      return(0.0);
   }

   dum=sqrtf(dum);
   a=A/dum; b=B/dum; c=C/dum;
   d=1/dum;

   return( reflection_cross_correlation(image,dim,a,b,c,d) ) ;
}

float optimizeNormalVector(short *image,DIM dim, float *A, float *B, float *C)
{
   //int i;  /* loop indices */
   int pmin;

   //float dum;  /* general purpose dummy variable */
   float cc;
   float x[3];

   float T[16];
   float *invT;

   float z0;

   z0 = (1- (*C))/ (*C);

   T[0]=1.0; T[1]=0.0; T[2]=0.0; T[3]=0.0;
   T[4]=0.0; T[5]=1.0; T[6]=0.0; T[7]=0.0;
   T[8]=0.0; T[9]=0.0; T[10]=1.0; T[11]=-z0;
   T[12]=0.0; T[13]=0.0; T[14]=0.0; T[15]=1.0;

   *A /= (*C); 
   *B /= (*C); 
   *C = 1.0;

   invT = inv4(T);
   gimage=resliceImage(image, dim, dim, invT, LIN);
   free(invT);

   gdim=dim;

   x[0]=vertex[0][0]=(*A);
   x[1]=vertex[0][1]=(*B);
   x[2]=vertex[0][2]=(*C);
   value[0]=-reflection_cross_correlation2(gimage,gdim,x[0],x[1],x[2]);

   x[0]=vertex[1][0]=2.0/(dim.nx*dim.dx);
   x[1]=vertex[1][1]=0.0;
   x[2]=vertex[1][2]=0.0;
   value[1]=-reflection_cross_correlation2(gimage,gdim,x[0],x[1],x[2]);

   x[0]=vertex[2][0]=0.0;
   x[1]=vertex[2][1]=2.0/(dim.ny*dim.dy);
   x[2]=vertex[2][2]=0.0;
   value[2]=-reflection_cross_correlation2(gimage,gdim,x[0],x[1],x[2]);

   x[0]=vertex[3][0]=0.0;
   x[1]=vertex[3][1]=0.0;
   x[2]=vertex[3][2]=2.0/(dim.nz*dim.dz);
   value[3]=-reflection_cross_correlation2(gimage,gdim,x[0],x[1],x[2]);

   pmin=dsm();

   x[0]=vertex[pmin][0];
   x[1]=vertex[pmin][1];
   x[2]=vertex[pmin][2];

   cc=reflection_cross_correlation2(gimage,gdim,x[0],x[1],x[2]);

   *A= x[0]/(1+x[2]*z0);
   *B= x[1]/(1+x[2]*z0);
   *C= x[2]/(1+x[2]*z0);

   free(gimage);

   return(cc);
}

float reflection_cross_correlation(short *image, DIM dim, float a, float b, float c, float d)
{
   //static int dumcount=0;
   int   i,j,k;
   float x,y,z;
   float dp;        
   float dum1,dum2,dum3;
   int q;
   int np;
   //int nv;
   int N;

   float ic,jc,kc;

   float a1,b1,c1;

   float f,g;
   float Sfg,Sgg;
   float Sg;

   ic = (dim.nx-1.0)/2.0;
   jc = (dim.ny-1.0)/2.0;
   kc = (dim.nz-1.0)/2.0;

   np=dim.nx*dim.ny;
   //nv=dim.nz*np;

   a1 = a * dim.dx;
   b1 = b * dim.dy;
   c1 = c * dim.dz;

   Sf=Sg=0.0;
   Sff=Sfg=Sgg=0.0;
   N=q=0; 
   for(k=0;k<dim.nz;k++)
   {
      dum1=(k-kc)*c1;
      for(j=0;j<dim.ny;j++)
      {
         dum2=dum1+(j-jc)*b1;
         for(i=0;i<dim.nx;i++)
         {
            f=image[q];

            dum3=dum2 + (i-ic)*a1;

            /* if(f>thresh && (d-dum3)<0 ) */
            if(f!=0.0 && (d-dum3)<0 )
            {

               dp=2.0*(d-dum3);

               x=i+a*dp/dim.dx;
               y=j+b*dp/dim.dy;
               z=k+c*dp/dim.dz;

	       		g=linearInterpolator(x,y,z,image,dim.nx,dim.ny,dim.nz,np);

               /* if(g>thresh) */
               if(g!=0.0)
               {
                  Sfg += f*g;
                  Sgg += g*g;
                  Sg += g;
                  Sff += f*f;
                  Sf += f;

                  N++;
               }
            }

            q++;
 
         }
      }
   }
   
   if( N==0 || ((Sff - Sf*Sf/N)*(Sgg - Sg*Sg/N))==0.0 )
      return(0.0);
   else
      return( (Sfg-Sg*Sf/N)/sqrtf((Sff-Sf*Sf/N)*(Sgg-Sg*Sg/N)) );
}

// have to limit the search considering the fact that the input image will be almost PIL
void findInitialNormalVector(short *image, DIM dim, float *A, float *B,float *C)
{
   float dum;        /* dummy variable */
   float a,b,c;      /* direction cosines  ax+by+cz=d */
   float d;          /* the d in ax+by+cz=d */
   float A1,B1,C1;   /* Ax+By+Cz=1 */
   float cc;
   float ccmax;

   /* Coordinates of the input image "center of gravity" in mm. 
   Origin is taken to be the center of the image volume. */
   float x_cm,y_cm,z_cm;

   double pi;
   double phi0=20.0;
   double delphi=2.0;
   int nrings; // number of rings
   double ringlength; // length of a ring
   double *cumulativelength;
   double *theta;
   double *phi;
   double totallength; // sum of all ring lengths
   double arclength;
   int N; // number of samples

   pi = 4.0*atan(1.0);

   // convert from degrees to radians
   phi0 = pi*phi0/180.0;
   delphi = pi*delphi/180.0;

   nrings = (int)ceil(phi0/delphi) + 1;

   cumulativelength = (double *)calloc(nrings, sizeof(double));

   totallength=0.0;
   for(int i=0; i<nrings; i++)
   {
      ringlength = 2*pi*sin(i*delphi);

      totallength += ringlength;

      cumulativelength[i] = totallength;
   }

   N = (int)floor( totallength/delphi ) + 1;

   theta = (double *)calloc(N, sizeof(double));
   phi = (double *)calloc(N, sizeof(double));

   theta[0]=0.0;
   phi[0]=0.0;

   for(int i=1; i<N; i++)
   {
      arclength = i*delphi;

      for(int j=1; j<nrings; j++)
      {
         if( arclength<=cumulativelength[j] && arclength>cumulativelength[j-1])
         {
            arclength -= cumulativelength[j-1];
            phi[i] = j*delphi;
            theta[i] = arclength/sin(phi[i]);
            break;
         }
      }
   }

   // compute x_cm,y_cm,z_cm the coordinates of the image "center of gravity" 
   // in (mm) with respect to the image volume center as origin
   compute_cm(image, dim.nx, dim.ny, dim.nz, dim.dx, dim.dy, dim.dz, &x_cm, &y_cm, &z_cm);

//printf("\n******x_cm=%7.3f y_cm=%7.3f z_cm=%7.3f (mm)\n",x_cm,y_cm,z_cm); 
//printf("i = %f\n", (x_cm + dim.dx*(dim.nx-1.0)/2.0)/1.5 );
//printf("j = %f\n", (y_cm + dim.dy*(dim.ny-1.0)/2.0)/0.859375 );
//printf("k = %f\n", (z_cm + dim.dz*(dim.nz-1.0)/2.0)/0.859375 );
	
	// Sff must be computed before calling reflection_cross_correlation()
	// or reflection_cross_correlation2()
	Sff=0.0;
	Sf=0.0;
	for(int i=0;i<dim.nx*dim.ny*dim.nz;i++)
	{
		dum=image[i];
		if(dum!=0) 
		{
			Sff += (dum*dum);
			Sf += dum;
		}
   }

   ccmax=0.0;

   for(int i=0;i<N;i++) 
   {
      /* The samples theta and z define a direction in space. Find the 
      unit vector (a,b,c) in that direction. */
      a=(float)(sin(phi[i])*cos(theta[i]));
      b=(float)(sin(phi[i])*sin(theta[i]));
      c=(float)cos(phi[i]);

      d = a*x_cm + b*y_cm + c*z_cm; 

      /* make sure d is non-negative */
      if(d<0.0)  
      {
         a *= -1.0; 
         b *= -1.0; 
         c *= -1.0; 
         d *= -1.0;
      }

      /* find the cross-correlation between image and its reflection 
      about the plane ax+by+cz=d */
      cc=reflection_cross_correlation(image,dim,a,b,c,d);

      if(cc>ccmax)
      {
         ccmax=cc;
         A1=a/d;
         B1=b/d;
         C1=c/d;
      }
   }

   dum=(float)sqrt((double)A1*A1 + B1*B1 + C1*C1 ); 
   a=A1/dum; b=B1/dum; c=C1/dum; d=1./dum;
//printf("\nInitial guess:");
//printf("\nplane of symmetry: (%7.3f,%7.3f,%7.3f).(x,y,z) = %7.3f", a,b,c,d);
//printf("\ncross correlation = %6.4f\n",ccmax);

   cc=optimizeNormalVector(image,dim,&A1,&B1,&C1);
  
//dum=(float)sqrt((double)A1*A1 + B1*B1 + C1*C1 ); 
//a=A1/dum; b=B1/dum; c=C1/dum; d=1./dum;
//printf("\nRefined initial guess:");
//printf("\nplane of symmetry: (%7.3f,%7.3f,%7.3f).(x,y,z) = %7.3f", a,b,c,d);
//printf("\ncross correlation = %6.4f\n",cc);

   *A=A1; *B=B1; *C=C1;

   free(cumulativelength);
   free(phi);
   free(theta);
}

// (A,B,C)=(a,b,c)/d  (Ax+By+Cz=1)  (ax+by+cz=d)
float msp(short *im_in, int nx, int ny, int nz, float dx, float dy, float dz, float *A, float *B, float *C) 
{
   int low,high;
   int nv;
   DIM dim[3];
   short *image[3]; 
   float newvs;      /* new voxel size in mm */
   float cc;

   nv = nx*ny*nz;

   dim[0].nx=nx;
   dim[0].ny=ny;
   dim[0].nz=nz;
   dim[0].dx=dx;
   dim[0].dy=dy;
   dim[0].dz=dz;

   setLowHigh(im_in,nv,&low,&high);

   image[1]=image[2]=NULL;

   image[0]=(short *)calloc(nv,sizeof(short));
   if(image[0]==NULL) {
      printf("Error: memory allocation, aborting ...\n\n");
      exit(0);
   }

   for(int i=0;i<nv;i++)
   if(im_in[i] <= low || im_in[i] >= high) 
      image[0][i]=0;
   else
      image[0][i]=im_in[i];

   dim[2].nx=32;
   newvs = (dim[0].dx) * (dim[0].nx - 1)/(dim[2].nx -1);
   dim[2].ny= (int)( (dim[0].dy) * (dim[0].ny - 1)/newvs  + 1 + 0.5 );
   dim[2].nz= (int)( (dim[0].dz) * (dim[0].nz - 1)/newvs  + 1 + 0.5 );
   if(dim[2].nx!=1) dim[2].dx=dim[0].dx*(dim[0].nx-1)/(dim[2].nx-1); else dim[2].dx=dim[0].dx;
   if(dim[2].ny!=1) dim[2].dy=dim[0].dy*(dim[0].ny-1)/(dim[2].ny-1); else dim[2].dy=dim[0].dy;
   if(dim[2].nz!=1) dim[2].dz=dim[0].dz*(dim[0].nz-1)/(dim[2].nz-1); else dim[2].dz=dim[0].dz;
   image[2]=resizeXYZ(image[0], dim[0].nx,dim[0].ny,dim[0].nz, dim[0].dx,dim[0].dy,dim[0].dz,
   dim[2].nx,dim[2].ny,dim[2].nz, dim[2].dx,dim[2].dy,dim[2].dz);

   dim[1].nx=64;
   newvs = (dim[0].dx) * (dim[0].nx - 1)/(dim[1].nx -1);
   dim[1].ny= (int)( (dim[0].dy) * (dim[0].ny - 1)/newvs  + 1 + 0.5 );
   dim[1].nz= (int)( (dim[0].dz) * (dim[0].nz - 1)/newvs  + 1 + 0.5 );
   if(dim[1].nx!=1) dim[1].dx=dim[0].dx*(dim[0].nx-1)/(dim[1].nx-1); else dim[1].dx=dim[0].dx;
   if(dim[1].ny!=1) dim[1].dy=dim[0].dy*(dim[0].ny-1)/(dim[1].ny-1); else dim[1].dy=dim[0].dy;
   if(dim[1].nz!=1) dim[1].dz=dim[0].dz*(dim[0].nz-1)/(dim[1].nz-1); else dim[1].dz=dim[0].dz;
   image[1]=resizeXYZ(image[0], dim[0].nx,dim[0].ny,dim[0].nz, dim[0].dx,dim[0].dy,dim[0].dz,
   dim[1].nx,dim[1].ny,dim[1].nz, dim[1].dx,dim[1].dy,dim[1].dz);

   findInitialNormalVector(image[2],dim[2], A, B, C);
   cc=optimizeNormalVector(image[1],dim[1], A, B, C); 
//cc=optimizeNormalVector(image[0],dim[0], A, B, C); 

   free(image[0]);
   free(image[1]);
   free(image[2]);

   return(cc);
}

int save_as_ppm(const char *filename, int nx, int ny, unsigned char *R, unsigned char *G, unsigned char *B)
{
  FILE *fp;
  int np;

  np = nx*ny;

  fp = fopen(filename,"w");
  if(fp == NULL) return(1);

  fprintf(fp,"P6\n");
  fprintf(fp,"# Created by Automatic Registration Toolbox \n");
  fprintf(fp,"%d %d\n",nx, ny);
  fprintf(fp,"255\n");

  for(int i=0; i<np; i++)
  {
    fwrite(R+i, 1, 1, fp);
    fwrite(G+i, 1, 1, fp);
    fwrite(B+i, 1, 1, fp);
  }

  fclose(fp);

  if(opt_png)
  {
    char *pngfilename;
    char cmnd[DEFAULT_STRING_LENGTH];
    int L;

    L = strlen(filename);
    // there was bug here: I had length L instead of L+1 for pngfilename
    // this was simple to fix but very hard to find
    // as a result free(pngfilename) was failing
    pngfilename = (char *)calloc(L+1,sizeof(char)); 

    stpcpy(pngfilename, filename);
    pngfilename[L-1]='g';
    pngfilename[L-2]='n';
    pngfilename[L-3]='p';

    snprintf(cmnd,sizeof(cmnd),"pnmtopng %s > %s",filename,pngfilename); 
    //if(opt_png) system(cmnd);
    if(opt_png) (void)system(cmnd);

    free(pngfilename);
  }

  if(opt_ppm == NO ) remove(filename);

  return(0);
}

void computeTmsp(char *orientation, short *volOrig, DIM dim, float *Tmsp)
{
   //float cc; // a variable to store correlation coefficient values
   float R[16], T[16];
   float TPIL[16]; // Transformation from original to PIL orientation
   float dum; // dummy floating point variable
   float A,B,C,a,b,c,d; // (A,B,C) or (a,b,c) define the MSP as: Ax+By+Cz=1 or ax+by+cz=d
   float dxPIL, dyPIL, dzPIL; // voxel dimensions in PIL orientation
   int nxPIL, nyPIL, nzPIL; // matrix dimensions in PIL orientation
   short *volumePIL; // original volume after reorientation to PIL

   // Compute the transformation from original to PIL orientation.
   PILtransform(orientation, TPIL);

   // reorient the original volume to PIL orientation
   volumePIL=reorientVolume(volOrig,dim.nx,dim.ny,dim.nz,dim.dx,dim.dy,dim.dz,TPIL,nxPIL,nyPIL,nzPIL,dxPIL,dyPIL,dzPIL);

   // determine the MSP from the PIL orineted volume
   //cc=msp(volumePIL, nxPIL, nyPIL, nzPIL, dxPIL, dyPIL, dzPIL, &A, &B, &C);
   (void)msp(volumePIL, nxPIL, nyPIL, nzPIL, dxPIL, dyPIL, dzPIL, &A, &B, &C);
   delete volumePIL;

   // determine (a,b,c) from (A,B,C)
   dum=(float)sqrt( (double)(A*A + B*B + C*C) );
   //if(dum>0.0) // jus in case
   {
      a= A/dum; b=B/dum; c=C/dum; d=1./dum;
   }

   // ensure c>0.0, this ensures that the unit normal to the MSP is approximately in +z direction
   if(c<0.0) { a = -a; b = -b; c = -c; d = -d; A=-A; B=-B; C=-C;}

/////////////////////////////////
// Manual Intervention to correct OASIS subject 0221
//a=0.04; b=0.006; c=0.999; d=4.168;
//A=0.04/d; B=0.006/d; C=0.999/d; 
/////////////////////////////////

//
//printf("\n****Interhemispheric correlation = %6.4f\n",cc);
//printf("\n****Estimated mid-sagittal plane: (%7.3fx) + (%7.3fy) + (%7.3fz) = %7.3f (mm)\n", a,b,c,d);

   // ensures MSP on volumePIL pass through the center of the FOV
   T[0]=1.0;  T[1]=0.0;  T[2]=0.0;  T[3]=-a*d;
   T[4]=0.0;  T[5]=1.0;  T[6]=0.0;  T[7]=-b*d;
   T[8]=0.0;  T[9]=0.0;  T[10]=1.0; T[11]=-c*d;
   T[12]=0.0; T[13]=0.0; T[14]=0.0; T[15]=1.0;

   if(c>1.0) c=1.0; // just in case to prevent acos from getting into trouble
   dum = (float)acos((double)c);

   // ensures that MSP on volumePIL is parallel to the x-y plane
   rotate(R, dum, B, -A, 0.0);

   // Tmsp=R*T*TPIL transforms the volOrig to MSP aligned PIL
   multi(T, 4, 4,  TPIL, 4,  4, Tmsp);
   multi(R, 4, 4,  Tmsp, 4,  4, Tmsp);
}

void combine_warps_and_trans(int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T)
{
   float  x,y,z;   
   float  xx,yy,zz;   
   float xc,yc,zc;
   float *invT;		
   int v;
   int np;

   np = nx*ny;

   invT=inv4(T);

   xc=dx*(nx-1)/2.0;     /* +---+---+ */
   yc=dy*(ny-1)/2.0;
   zc=dz*(nz-1)/2.0;

   for(int k=0;k<nz;k++) 
   for(int j=0;j<ny;j++) 
   for(int i=0;i<nx;i++) 
   {
      v = k*np + j*nx + i;
      // (i*dx-xc) converts from image coordinates (i,j,z) to (x,y,z) coordinates
      xx = (i*dx - xc) + Xwarp[v];
      yy = (j*dy - yc) + Ywarp[v];
      zz = (k*dz - zc) + Zwarp[v];

      x = ( invT[0]*xx +invT[1]*yy +invT[2]*zz  +invT[3]  );
      y = ( invT[4]*xx +invT[5]*yy +invT[6]*zz  +invT[7]  );
      z = ( invT[8]*xx +invT[9]*yy +invT[10]*zz +invT[11] );

      Xwarp[v] = x - i*dx + xc;
      Ywarp[v] = y - j*dy + yc;
      Zwarp[v] = z - k*dz + zc;
   }

   free(invT);
}

// xvec, yvec, and zvec should be interpreted in DICOM LPS system
// Note: DICOM's LPS system is  is a right-handed system  LxP=S
// angles passed as input should be in degrees
// translations passed as input should be in mm
void forwardTCSAP(float *xvec, float *yvec, float *zvec, float *TLHC, float *angle, float *translation, DIM dim)
{
   float pi;
   float TLHC_ijk[3]; // a given position in ijk system
   float TLHC_xyz[3]; // a given position in xyz system
   float alpha, beta, gamma;

   pi = 4.0*atanf(1.0);

   alpha = angle[0]*pi/180.0;
   beta  = angle[1]*pi/180.0;
   gamma = angle[2]*pi/180.0;

   xvec[0] =  cosf(beta)  * cosf(gamma);
   xvec[1] = -cosf(alpha) * sinf(gamma) - sinf(alpha) * sinf(beta) * cosf(gamma);
   xvec[2] = -sinf(alpha) * sinf(gamma) + cosf(alpha) * sinf(beta) * cosf(gamma);

   yvec[0] =  cosf(beta)  * sinf(gamma);
   yvec[1] =  cosf(alpha) * cosf(gamma) - sinf(alpha) * sinf(beta) * sinf(gamma);
   yvec[2] =  sinf(alpha) * cosf(gamma) + cosf(alpha) * sinf(beta) * sinf(gamma);

   normalizeVector(xvec, 3);
   normalizeVector(yvec, 3);

   // right-handed cross product between xvec and yvec seems to give -zvec!
   zvec[0] = -(xvec[1]*yvec[2] - xvec[2]*yvec[1]);
   zvec[1] = -(xvec[2]*yvec[0] - xvec[0]*yvec[2]);
   zvec[2] = -(xvec[0]*yvec[1] - xvec[1]*yvec[0]);

   // it is IMPORTANT to take the integer part so as to drop the fractional part, that way the
   // in-plane FOV center is always on a voxel center.  This would not be the case when a voxel dimension
   // nx, ny, or nz is even if we didn't drop the fractional part.
   // discovered that this does not need to be the case for the slice direction, at least in 2D images
   TLHC_ijk[0] = -(int)(dim.nx/2.0);
   TLHC_ijk[1] = -(int)(dim.ny/2.0);
   TLHC_ijk[2] = -(dim.nz-1.0)/2.0;  // taking (int) is not necessary here

   // converting to units of mm in xyz system
   TLHC_xyz[0] = TLHC_ijk[0]*dim.dx; 
   TLHC_xyz[1] = TLHC_ijk[1]*dim.dy;
   TLHC_xyz[2] = TLHC_ijk[2]*dim.dz;

   TLHC[0] = xvec[0]*TLHC_xyz[0] + yvec[0]*TLHC_xyz[1] + zvec[0]*TLHC_xyz[2] + translation[0];
   TLHC[1] = xvec[1]*TLHC_xyz[0] + yvec[1]*TLHC_xyz[1] + zvec[1]*TLHC_xyz[2] + translation[1];
   TLHC[2] = xvec[2]*TLHC_xyz[0] + yvec[2]*TLHC_xyz[1] + zvec[2]*TLHC_xyz[2] + translation[2];

   return;
}

void  backwardTCSAP(float *xvec, float *yvec, float *angle)
{
   float pi;
   float pi_over_2;
   float alpha, beta, gamma;

   pi = 4.0*atanf(1.0);
   pi_over_2 = pi/2.0;

   normalizeVector(xvec,3);
   normalizeVector(yvec,3);

   if (xvec[0] == 0.0)
      gamma = pi_over_2;
   else
      gamma = atanf(yvec[0]/xvec[0]);

   if ( (int)(1000.0*gamma) == (int)(1000.0*pi_over_2) )
      alpha=0.0;
   else
      alpha= asinf( (yvec[2]-xvec[2]*tanf(gamma)) / (cosf(gamma)+sinf(gamma)*tanf(gamma)) );

   if ( (int)(1000.0*gamma) == (int)(1000.0*pi_over_2)  || alpha == 0.0)
      beta=0.0;
   else 
      beta = asinf( (-xvec[1]-cosf(alpha)*sinf(gamma)) / (cosf(gamma)*sinf(alpha)) );

   angle[0] = alpha*180.0/pi;
   angle[1] = beta*180.0/pi;
   angle[2] = gamma*180.0/pi;

   return;
}

void update_qsform(nifti_1_header &hdr, const char *neworient)
{
   char oldorient[4];
   float T[16];
   mat44 R;

   getNiftiImageOrientation(hdr, oldorient);

   float T_neworient_to_PIL[16];
   float T_PIL_to_oldorient[16];
   float T_neworient_to_oldorient[16];
   PILtransform(neworient, T_neworient_to_PIL);
   inversePILtransform(oldorient, T_PIL_to_oldorient);
   multi(T_PIL_to_oldorient, 4, 4, T_neworient_to_PIL, 4, 4, T_neworient_to_oldorient);

   hdr.sform_code = NIFTI_XFORM_ALIGNED_ANAT;
   T[0]=hdr.srow_x[0]; T[1]=hdr.srow_x[1]; T[2]=hdr.srow_x[2]; T[3]=hdr.srow_x[3];
   T[4]=hdr.srow_y[0]; T[5]=hdr.srow_y[1]; T[6]=hdr.srow_y[2]; T[7]=hdr.srow_y[3];
   T[8]=hdr.srow_z[0]; T[9]=hdr.srow_z[1]; T[10]=hdr.srow_z[2]; T[11]=hdr.srow_z[3];
   T[12]=T[13]=T[14]=0.0; T[15]=1.0;

   multi( T, 4, 4, T_neworient_to_oldorient, 4, 4, T);
   hdr.srow_x[0]=T[0]; hdr.srow_x[1]=T[1]; hdr.srow_x[2]=T[2]; hdr.srow_x[3]=T[3];
   hdr.srow_y[0]=T[4]; hdr.srow_y[1]=T[5]; hdr.srow_y[2]=T[6]; hdr.srow_y[3]=T[7];
   hdr.srow_z[0]=T[8]; hdr.srow_z[1]=T[9]; hdr.srow_z[2]=T[10]; hdr.srow_z[3]=T[11];

   hdr.qform_code = NIFTI_XFORM_ALIGNED_ANAT;
   R.m[0][0]=T[0];  R.m[0][1]=T[1];  R.m[0][2]=T[2];  R.m[0][3]=T[3];
   R.m[1][0]=T[4];  R.m[1][1]=T[5];  R.m[1][2]=T[6];  R.m[1][3]=T[7];
   R.m[2][0]=T[8];  R.m[2][1]=T[9];  R.m[2][2]=T[10]; R.m[2][3]=T[11];
   R.m[3][0]=T[12]; R.m[3][1]=T[13]; R.m[3][2]=T[14]; R.m[3][3]=T[15];

   nifti_mat44_to_quatern(R,  &(hdr.quatern_b), &(hdr.quatern_c), &(hdr.quatern_d),
   &(hdr.qoffset_x), &(hdr.qoffset_y), &(hdr.qoffset_z), 
   &(hdr.pixdim[1]), &(hdr.pixdim[2]), &(hdr.pixdim[3]), &(hdr.pixdim[0]));
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Find the transformation (pilT) that takes the image to standard PIL orientation
///////////////////////////////////////////////////////////////////////////////////////////////
void find_pil_transformation(char *imfile, DIM dim, float *pilT, float *AC, float *PC, float *VSPS)
{
   float ac[4], pc[4];  
   char orientation[4]="";
   char modelfile[1024]="";

   float Tmsp[16]; // transforms image to MSP aligned PIL orientation

   detect_AC_PC_MSP(imfile, orientation, modelfile, AC, PC, VSPS, Tmsp, 0, 0);

   // convert the AC/PC from (i,j,k) in original space to (x,y,z) in PIL space
   for(int i=0; i<4; i++) ac[i] = AC[i];
   for(int i=0; i<4; i++) pc[i] = PC[i];
   orig_ijk_to_pil_xyz(Tmsp, dim, ac, pc);

   ACPCtransform(pilT, Tmsp, ac, pc, 0);  // 0 is equivalent to opt_M=YES
}
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
// Find the transformation (pilT) that takes the image to standard PIL orientation
///////////////////////////////////////////////////////////////////////////////////////////////
void find_pil_transformation(char *imfile, DIM dim, float *pilT)
{
   float ac[4], pc[4];  
   char orientation[4]="";
   char modelfile[1024]="";

   float AC[4]={0.0, 0.0, 0.0, 1.0};
   float PC[4]={0.0, 0.0, 0.0, 1.0};
   float VSPS[4]={0.0, 0.0, 0.0, 1.0};

   float Tmsp[16]; // transforms image to MSP aligned PIL orientation

   detect_AC_PC_MSP(imfile, orientation, modelfile, AC, PC, VSPS, Tmsp, 0, 0);

   // convert the AC/PC from (i,j,k) in original space to (x,y,z) in PIL space
   for(int i=0; i<4; i++) ac[i] = AC[i];
   for(int i=0; i<4; i++) pc[i] = PC[i];
   orig_ijk_to_pil_xyz(Tmsp, dim, ac, pc);

   ACPCtransform(pilT, Tmsp, ac, pc, 0);  // 0 is equivalent to opt_M=YES
}
///////////////////////////////////////////////////////////////////////////////////////////////
