///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// Copyright (C) 2024 Babak A. Ardekani, PhD - All Rights Reserved.  //
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>  //      required by time(), stat()
#include <time.h>       //      required by time()
#include <sys/stat.h>   //      required by stat() 
#include <unistd.h>
#include <volume.h>
#include <spm_analyze.h>
#include <nifti1_io.h>
#include <babak_lib.h>
#include <smooth.h>
#include <minmax.h>
#include <stats.h>
#include <interpolator.h>

#define YES 1
#define NO 0

int NX=512;
int NY=512;
int NP;

// these are read from the atlas file which is in NIFTI format
int UPPER_LEFT_i;
int UPPER_LEFT_j;
int LOWER_RIGHT_i;
int LOWER_RIGHT_j;

//////////////////////////////////////////////////////////////////////////////////////////////////
float VSPS[4]={0.0, 0.0, 0.0, 1.0};
float AC[4]={0.0, 0.0, 0.0, 1.0};
float PC[4]={0.0, 0.0, 0.0, 1.0};

//int np;
float dx, dy, dz;
int N2;	// N*N

float *Xwarp, *Ywarp, *Zwarp;
int niter=4;

float *ARobj;	// array extracted from object and target images
float *ARtrg;

int Wx, Wy;
int Lx, Ly;

nifti_1_header sub_hdr;
nifti_1_header output_hdr;
short *subj_volume=NULL;
int Snx, Sny, Snz;
float Sdx, Sdy, Sdz;

float Tacpc[16];

float hampel_origin[2];
float hampel_axis[2];
int inferior_genu[2];
int inferior_splenium[2];
int anterior_point[2];
int posterior_point[2];
int posterior_genu[2];
int rostrum[2];
int cc_length_fifth;
int cc_length_third;
int cc_length_half;
int cc_length;

int *cci, *ccj;

int *medi, *medj; // (i,j) components of the medial axis pixels 0, 1, ..., nm-1
int nm=0;  // number of pixels on the medial axis
short *dist;

// normal unit vectors to the medial axis
// program is designed such that the normal vectors to the medial axis always point towards 
// the upper boundary of the CC.
float *normi, *normj; 

int *ubi, *ubj, *lbi, *lbj;
int nb; // number of border pixels
int nub; // number of upper boundary pixels
int nlb; // number of lower boundary pixels

float *lbindx, *ubindx;

int acpoint;
int pcpoint;

float ACx;
float PCx;

float ACi;
float PCi;
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
int opt;

static struct option options[] =
{
  {"-i", 1,  'i'},
  {"-v", 0,  'v'},
  {"-verbose", 0,  'v'},
  {"-V", 0,  'V'},
  {"-version", 0,  'V'},
  {"-Version", 0,  'V'},
  {"-h", 0,  'h'},
  {"-help", 0,  'h'},
  {"-Help", 0,  'h'},
  {"-n", 1,  'n'},
  {"-o", 1, 'o'},
  {"-csv", 1,  'c'},
  {"-H", 0,  'H'},
  {"-Hampel", 0,  'H'},
  {"-W", 0,  'W'},
  {"-Witelson", 0,  'W'},
  {"-png", 0,  'p'},
  {"-lm",1,'l'},  // landmark  file
  {"-A", 1,  'L'},
  {"-T", 1,  'T'},
  {"-mrx", 0,  'm'},
  {"-box", 0,  'b'},
  {"-secret", 0,  's'},
  {"-border", 0, 'B'},
  {"-cc", 1, 'C'},
  {"-t",1,'t'},
  {"-thresh",1,'t'},
  {"-threshold",1,'t'},
  {"-Atlas", 1,  'a'},
  {"-atlas",1,'a'},
  {"-a",1,'a'},
  {0, 0,  0}
};

int opt_cc=NO;
int opt_box=NO;
int opt_W=NO;
int opt_H=NO;
int opt_border=NO;
//////////////////////////////////////////////////////////////////////////////////////////////////

void print_help();

void invert_symmetric_matrix(float *X, int p)
{
   float *D;
   float *Ut;
   float *DUt;
   float *Xinv;

   // ensure memory was allocated for X and G
   if(X==NULL)
   {
      printf("\n\ninvert_symmetric_matrix(): NULL input matrix encountered.\n\n");
      exit(-1);
   }

   // make sure p is positive
   if(p<=0)
   {
      printf("\n\ninvert_symmetric_matrix(): Non-positive matrix dimension encountered.\n\n");
      exit(-1);
   }

   // Diagnolized X, such that X = U diag(D) U'
   // This function returns U' in place of X
   D=diagATA_float(X, p, 'L');
   Ut=X;

   //printMatrix(D,p,1,"D",NULL);

   if(D==NULL)
   {
      printf("\n\nninvert_symmetric_matrix(): Could not diagnolize X.\n\n");
      exit(-1);
   }

   DUt=(float *)calloc(p*p,sizeof(float));
   if(DUt==NULL) 
   {
      printf("\n\ninvert_symmetric_matrix(): Memory allocation problem.\n\n");
      exit(-1);
   }

   // compute diag(1/D) * U'
   for(int i=0; i<p; i++)
   {
      for(int j=0; j<p; j++) DUt[i*p + j]=Ut[i*p +j]/D[i];
   }

   //mat_trans_mat(Ut, p, p, DUt, p, DUt);
   //printMatrix(DUt,p,p,"UDUt",NULL);

   Xinv=(float *)calloc(p*p,sizeof(float));
   if(Xinv==NULL) 
   {
      printf("\n\ninvert_symmetric_matrix(): Memory allocation problem.\n\n");
      exit(-1);
   }

   // Compute Xinv = U * diag(1/D) * U'
   mat_trans_mat(Ut, p, p, DUt, p, Xinv);
   //printMatrix(Xinv,p,p,"Xinv",NULL);

   for(int i=0; i<p*p; i++) X[i]=Xinv[i];

   free(DUt);
   free(D);
   free(Xinv);
}

void computeWarpField(short *obj, short *trg, float sd, int bbnp, int bbnx, int bbny)
{
   int xopt, yopt;

   float Vt, Vo;
   float Sx, Sx2, Sy, Sy2, Sxy;
   float num, den;
   float CC, CCMAX;
   float *Xw, *Yw;
   float *Xtmp, *Ytmp;

   short *HRtrg;	// high res. target and object images
   short *HRobj;

   float HRdx, HRdy;
   int HRnx, HRny;

   for(int n=0; n<bbnp; n++) Xwarp[n]=Ywarp[n]=0.0;

   for(int iter=0; iter<niter; iter++)
   {
      // printf("\n\tIteration %d\n",iter);

      HRdx = (float)(dx*pow(2.0,(niter-iter-1.0)));
      HRdy = (float)(dy*pow(2.0,(niter-iter-1.0)));

      HRnx = (int)(bbnx/pow(2.0,(niter-iter-1.0)) + 0.5 );
      HRny = (int)(bbny/pow(2.0,(niter-iter-1.0)) + 0.5 );

      //printf("\tMatrix size = %d x %d (voxels)\n", HRnx, HRny);
      //printf("\tVoxel size = %8.6f x %8.6f (mm3)\n", HRdx,HRdy);

      Xw=resizeXY(Xwarp, bbnx, bbny, dx, dy, HRnx, HRny, HRdx, HRdy);
      Yw=resizeXY(Ywarp, bbnx, bbny, dx, dy, HRnx, HRny, HRdx, HRdy);

      {
         float *tmp;
         float StandDev;

         StandDev = (0.5/log(2.0)) * ( HRdx*HRdx - dx*dx );

         if(StandDev>0.0)
         {
            if( ( HRdx*HRdx - dx*dx ) > 0.0 )
            {
               StandDev=(float)( sqrt( (0.5/log(2.0)) * ( HRdx*HRdx - dx*dx ) )/dx );
            }

            tmp = smoothXY(obj,bbnx,bbny,StandDev);
            HRobj=computeReslicedImage(tmp, bbnx, bbny, dx, dy, HRnx, HRny, HRdx, HRdy, Xw, Yw);
            free(tmp);
         }
         else
         {
            HRobj=computeReslicedImage(obj, bbnx, bbny, dx, dy, HRnx, HRny, HRdx, HRdy, Xw, Yw);
         }
      }

      HRtrg=resizeXY(trg, bbnx ,bbny, dx, dy, HRnx, HRny, HRdx, HRdy);

      for(int j=0; j<HRny; j++)
      for(int i=0; i<HRnx; i++)
      {
			extractArray(HRtrg, HRnx, HRny, i, j, Lx, Ly, ARtrg);

			Sy=Sy2=0.0;
			for(int n=0; n<N2; n++)
			{
				Vt = ARtrg[n];
				Sy += Vt;
				Sy2 += (Vt*Vt);
			}

			if( Sy == 0.0 )
			{
				Xw[j*HRnx + i] = 0.0;
				Yw[j*HRnx + i] = 0.0;
				continue;
			}

			CCMAX=0.0; 	// IMPORTANT: we are not interested in -tive correlations
					// if CMAX is set to -1, program given unexpected results

			xopt=yopt=0;

			for(int x=-Wx; x<=Wx; x++)
			for(int y=-Wy; y<=Wy; y++)
			{
				extractArray(HRobj, HRnx, HRny, i+x, j+y, Lx, Ly, ARobj);
	
				Sx=Sx2=Sxy=0.0;
				for(int n=0; n<N2; n++)
				{
					Vo = ARobj[n];
					Sx += Vo;
					Sx2 += (Vo*Vo);
					Sxy += (Vo*ARtrg[n]);
				}

				num = Sxy-Sx*Sy/N2;
				den = (float)sqrt( (double)(Sx2-Sx*Sx/N2)*(Sy2-Sy*Sy/N2) );
	
				if(den==0.0) continue;
	
				CC = num/den;
		
				if( CC>CCMAX ) { CCMAX=CC; xopt=x; yopt=y; }
			}
	
			Xw[j*HRnx + i] = xopt*HRdx;
			Yw[j*HRnx + i] = yopt*HRdy;
      }
	
      Xtmp=smoothXY(Xw, HRnx, HRny, sd);
      Ytmp=smoothXY(Yw, HRnx, HRny, sd);
	
      free(Xw); free(Yw);
	
      Xw=resizeXY(Xtmp, HRnx, HRny, HRdx, HRdy, bbnx, bbny, dx, dy);
      Yw=resizeXY(Ytmp, HRnx, HRny, HRdx, HRdy, bbnx, bbny, dx, dy);
		
      free(Xtmp); free(Ytmp);
		
      for(int n=0; n<bbnp; n++)
      {
         Xwarp[n] += Xw[n];
         Ywarp[n] += Yw[n];
      }
		
      free(Xw); free(Yw);
	
      free(HRobj);
      free(HRtrg);
   }
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void update_qsform( const char *imagefilename , float *matrix)
{
   FILE *fp;
   nifti_1_header hdr; // 348 bytes
   nifti1_extender ext; // 4 bytes
   char *extension=NULL;
   int extension_size=0;
   char *data=NULL;
   int data_size=0;
   char swapflg=0;
   mat44 R;
   size_t sizeread;

   fp = fopen(imagefilename,"r");
   if(fp==NULL) file_open_error(imagefilename);
   sizeread = fread(&hdr, sizeof(nifti_1_header), 1, fp);  
   if(sizeread != 1 ) errorMessage("Read error, aborting ...");

   if(hdr.dim[0]<1 || hdr.dim[0]>7)
   {
      swapniftiheader(&hdr);
      swapflg=1;
   }

   //if( opt_sform)
   {
      hdr.sform_code = NIFTI_XFORM_TALAIRACH;
      hdr.srow_x[0]=matrix[0]; hdr.srow_x[1]=matrix[1]; hdr.srow_x[2]=matrix[2]; hdr.srow_x[3]=matrix[3];
      hdr.srow_y[0]=matrix[4]; hdr.srow_y[1]=matrix[5]; hdr.srow_y[2]=matrix[6]; hdr.srow_y[3]=matrix[7];
      hdr.srow_z[0]=matrix[8]; hdr.srow_z[1]=matrix[9]; hdr.srow_z[2]=matrix[10]; hdr.srow_z[3]=matrix[11];
   }

   //if( opt_qform)
   {
      hdr.qform_code = NIFTI_XFORM_TALAIRACH;
      R.m[0][0]=matrix[0];  R.m[0][1]=matrix[1];  R.m[0][2]=matrix[2];  R.m[0][3]=matrix[3];
      R.m[1][0]=matrix[4];  R.m[1][1]=matrix[5];  R.m[1][2]=matrix[6];  R.m[1][3]=matrix[7];
      R.m[2][0]=matrix[8];  R.m[2][1]=matrix[9];  R.m[2][2]=matrix[10]; R.m[2][3]=matrix[11];
      R.m[3][0]=matrix[12]; R.m[3][1]=matrix[13]; R.m[3][2]=matrix[14]; R.m[3][3]=matrix[15];

      nifti_mat44_to_quatern( R,  &(hdr.quatern_b), &(hdr.quatern_c), &(hdr.quatern_d),
      &(hdr.qoffset_x), &(hdr.qoffset_y), &(hdr.qoffset_z), 
      &(hdr.pixdim[1]), &(hdr.pixdim[2]), &(hdr.pixdim[3]), &(hdr.pixdim[0]));
   }

   if( hdr.magic[0]=='n' && hdr.magic[1]=='+' && hdr.magic[2]=='1' )
   {
      sizeread=fread(&ext, sizeof(nifti1_extender), 1, fp);
      if(sizeread != 1 ) errorMessage("Read error, aborting ...");

      extension_size = (int)(hdr.vox_offset)-352;

      if( extension_size > 0 )
      {
         extension = (char *)calloc(extension_size, 1);
         sizeread=fread(extension, 1, extension_size, fp);
         if(sizeread != (size_t)extension_size) errorMessage("Read error, aborting ...");
      }

      data_size = 1;
      for(int i=1; i<=hdr.dim[0]; i++)
      {
         data_size *= hdr.dim[i];
      }
      data_size *= (hdr.bitpix/8);

      if( data_size > 0 )
      {
         data = (char *)calloc(data_size, 1);
         sizeread = fread(data, 1, data_size, fp);
         if(sizeread != (size_t)data_size) errorMessage("Read error, aborting ...");
      }
   }

   fclose(fp);

   fp = fopen(imagefilename,"w");
   if(fp==NULL) file_open_error(imagefilename);

   if(swapflg)
   {
      swapniftiheader(&hdr);
   }

   fwrite(&hdr, sizeof(nifti_1_header), 1, fp);

   if( hdr.magic[0]=='n' && hdr.magic[1]=='+' && hdr.magic[2]=='1' )
   {
      fwrite(&ext, sizeof(nifti1_extender), 1, fp);
      if( extension_size > 0 )
      {
         fwrite(extension, 1, extension_size, fp);
         delete extension;
      }
      if( data_size > 0 )
      {
         fwrite(data, 1, data_size, fp);
         delete data;
      }
   }
   fclose(fp);
}
//////////////////////////////////////////////////////////////////////////////////////////////////

short *find_subject_msp(char *imagefilename, char *prefix, char *lmfile)
{
   DIM input_dim, output_dim;
   short *msp;

   float *invT;

   char outputfilename[512];
   FILE *fp;

   char iporient[4]="";

   new_PIL_transform(imagefilename, lmfile , iporient, Tacpc, 0);

   input_dim.nx = Snx;
   input_dim.ny = Sny;
   input_dim.nz = Snz;
   input_dim.dx = Sdx;
   input_dim.dy = Sdy;
   input_dim.dz = Sdz;

   output_dim.nx = NX;
   output_dim.ny = NY;
   output_dim.nz = 1;
   output_dim.nt = 1;
   output_dim.dx = dx;
   output_dim.dy = dy;
   output_dim.dz = dz;
   output_dim.dt = 0.0;

   invT = inv4(Tacpc);
   msp = resliceImage(subj_volume,input_dim, output_dim,invT,LIN);
   free(invT);
   
   {
      sprintf(outputfilename,"%s_msp.mrx",prefix);
      fp = fopen(outputfilename,"w");

      if(fp != NULL)
      {
         printMatrix(Tacpc, 4, 4, "ART acpcdetect tilt correction matrix:", fp);
         fclose(fp);
      }
      else
      {
         printf("Could not write to %s.\n", outputfilename);
      }
   }

   {
      float T_ijk2xyz[16];
      float PIL2RAS[16];

      output_hdr = read_NIFTI_hdr(imagefilename);
      output_hdr.pixdim[1]=output_dim.dx; 
      output_hdr.pixdim[2]=output_dim.dy; 
      output_hdr.pixdim[3]=output_dim.dz;
      output_hdr.dim[1]=output_dim.nx; 
      output_hdr.dim[2]=output_dim.ny; 
      output_hdr.dim[3]=output_dim.nz;
      output_hdr.magic[0]='n'; output_hdr.magic[1]='+'; output_hdr.magic[2]='1';
      sprintf(output_hdr.descrip,"Created by ART yuki");

      sprintf(outputfilename,"%s_msp.nii",prefix);
      save_nifti_image(outputfilename, msp, &output_hdr);

      //////////////////////////////////////////////////////////////////////////////////
      // This part of the code adjusts the SFORM matrix of the output image

      inversePILtransform("RAS", PIL2RAS);

      ijk2xyz(T_ijk2xyz, output_dim.nx, output_dim.ny, output_dim.nz, output_dim.dx, output_dim.dy, output_dim.dz);
      multi(PIL2RAS, 4, 4,  T_ijk2xyz, 4,  4, Tacpc);

      update_qsform( (const char *)outputfilename, Tacpc );
      //////////////////////////////////////////////////////////////////////////////////
   }

   return(msp);
}

short *find_subject_msp_using_transformation(char *imagefilename, char *prefix, char *msp_transformation_file)
{
   DIM input_dim, output_dim;
   short *msp;
   char orientation[4]="";
   char modelfile[1024]="";
   int opt_T2=NO;

   float Tmsp[16]; // transforms volOrig to MSP aligned PIL orientation

   float ac[4], pc[4];  

   float *invT;

   char outputfilename[512];

   detect_AC_PC_MSP(imagefilename, orientation, modelfile, AC, PC, VSPS, Tmsp, 0, opt_T2);

   input_dim.nx = Snx;
   input_dim.ny = Sny;
   input_dim.nz = Snz;
   input_dim.dx = Sdx;
   input_dim.dy = Sdy;
   input_dim.dz = Sdz;

   output_dim.nx = NX;
   output_dim.ny = NY;
   output_dim.nz = 1;
   output_dim.nt = 1;
   output_dim.dx = dx;
   output_dim.dy = dy;
   output_dim.dz = dz;
   output_dim.dt = 0.0;

   // convert the AC/PC from (i,j,k) in original space to (x,y,z) in PIL space
   for(int i=0; i<4; i++) ac[i] = AC[i];
   for(int i=0; i<4; i++) pc[i] = PC[i];

   // convert the AC/PC from (i,j,k) in original space to (x,y,z) in PIL space
   orig_ijk_to_pil_xyz(Tmsp, input_dim, ac, pc);
   ACPCtransform(Tacpc, Tmsp, ac, pc, 0);

   // this part of the code locates the ac and pc locations on the MSP-AC-PC aligned sagittal slice in (i,j) coordinates.
   // the actuall ac location is actually approximately acpoint+0.5.  also pc is approx. pcpoint+0.5 
   // on the output PPM image, we mark points acpoint, acpoint+1, pcpoint and pcpoint+1 on rows 255 and 256.
   {
      float I2X[16];
      float X2I[16];

      for(int i=0; i<4; i++) ac[i] = AC[i];
      for(int i=0; i<4; i++) pc[i] = PC[i];

      ijk2xyz(I2X, input_dim.nx, input_dim.ny, input_dim.nz, input_dim.dx, input_dim.dy, input_dim.dz);

      multi(I2X, 4, 4,  ac, 4,  1, ac);
      multi(I2X, 4, 4,  pc, 4,  1, pc);

      multi(Tacpc, 4, 4,  ac, 4,  1, ac);
      multi(Tacpc, 4, 4,  pc, 4,  1, pc);

      ACx=ac[0]; 
      PCx=pc[0];

      xyz2ijk(X2I, output_dim.nx, output_dim.ny, output_dim.nz, output_dim.dx, output_dim.dy, output_dim.dz);

      multi(X2I, 4, 4,  ac, 4,  1, ac);
      multi(X2I, 4, 4,  pc, 4,  1, pc);

      acpoint = (int)ac[0];
      pcpoint = (int)pc[0];

      ACi=ac[0];
      PCi=pc[0];
   }

   loadTransformation(msp_transformation_file, Tacpc);

   invT = inv4(Tacpc);
   msp = resliceImage(subj_volume,input_dim, output_dim,invT,LIN);
   free(invT);

   {
      float T_ijk2xyz[16];
      float PIL2RAS[16];

      output_hdr = read_NIFTI_hdr(imagefilename);
      output_hdr.pixdim[1]=output_dim.dx;
      output_hdr.pixdim[2]=output_dim.dy; 
      output_hdr.pixdim[3]=output_dim.dz;
      output_hdr.dim[1]=output_dim.nx; 
      output_hdr.dim[2]=output_dim.ny; 
      output_hdr.dim[3]=output_dim.nz;
      output_hdr.magic[0]='n'; output_hdr.magic[1]='+'; output_hdr.magic[2]='1';
      sprintf(output_hdr.descrip,"Created by ART yuki");

      sprintf(outputfilename,"%s_msp.nii",prefix);
      save_nifti_image(outputfilename, msp, &output_hdr);

      //////////////////////////////////////////////////////////////////////////////////
      // This part of the code adjusts the SFORM matrix of the output image

      inversePILtransform("RAS", PIL2RAS);

      ijk2xyz(T_ijk2xyz, output_dim.nx, output_dim.ny, output_dim.nz, output_dim.dx, output_dim.dy, output_dim.dz);
      multi(PIL2RAS, 4, 4,  T_ijk2xyz, 4,  4, Tacpc);

      update_qsform( (const char *)outputfilename, Tacpc );
      //////////////////////////////////////////////////////////////////////////////////
   }

   return(msp);
}

//////////////////////////////////////////////////////////////////////////////////

void output_ppm(short *trg, short *cc_est, const char *prefix) 
{
   char outputfile[1024]="";
   short min=0, max=0;

   ////////////////////////////////////////////////////////////
   unsigned char *R, *G, *B;

   R = (unsigned char *)calloc(NP, 1);
   if(R == NULL)
   {
      memory_allocation_error("R");
   }

   G = (unsigned char *)calloc(NP, 1);
   if(G == NULL)
   {
      memory_allocation_error("G");
   }

   B = (unsigned char *)calloc(NP, 1);
   if(B == NULL)
   {
      memory_allocation_error("B");
   }

   minmax(trg, NP, min, max);

   if(max==0) max=1; // to avoid division by zero

   //////////////////////////////////////////////////////////////////////////////////
   // outputs visualization of Witelson's subdivisions
   if(opt_W)
   {
      for(int i=0; i<NX; i++)
      for(int j=0; j<NY; j++)
      {
         int v;

         v = j*NX + i;

         if(cc_est[v]>0)
         {
            if(i >= (posterior_point[0]-cc_length_fifth) ) // W7 
            {
               R[v] = 143;
               G[v] = 0;
               B[v] = 255;
            }
            else if(i >= (posterior_point[0]-cc_length_third) ) // W6 
            {
               R[v] = 75;
               G[v] = 0;
               B[v] = 130;
            }
            else if(i >= (posterior_point[0]-cc_length_half) ) // W5 
            {
               R[v] = 0;
               G[v] = 0;
               B[v] = 255;
            }
            else if(i >= (anterior_point[0]+cc_length_third) ) // W4 
            {
               R[v] = 0;
               G[v] = 255;
               B[v] = 0;
            } 
            else if(i <= posterior_genu[0] ) // W2 
            {
               R[v] = 255;
               G[v] = 127;
               B[v] = 0;
            } 
            else if(j >= posterior_genu[1] ) // W1 
            {
               R[v] = 255;
               G[v] = 0;
               B[v] = 0;
            } 
            else // w3 
            {
               R[v] = 255;
               G[v] = 255;
               B[v] = 0;
            }
         }
         else
         {
            R[v] = (unsigned char)(trg[v]*255.0/max);
            G[v] = (unsigned char)(trg[v]*255.0/max);
            B[v] = (unsigned char)(trg[v]*255.0/max);
         }
      }

      sprintf(outputfile,"%s_cc_witelson.ppm",prefix);
      save_as_ppm((const char *)outputfile, NX, NY, R, G, B);
   }
   //////////////////////////////////////////////////////////////////////////////////
   
   //////////////////////////////////////////////////////////////////////////////////
   // outputs visualization of Hampel's subdivisions
   if(opt_H)
   {
      float d;
      float pi;
      float i0, j0;
      float theta, costheta;
      pi = 4.0*atanf(1.0);

      i0 = hampel_origin[0];
      j0 = hampel_origin[1];

      for(int i=0; i<NX; i++)
      for(int j=0; j<NY; j++)
      {
         int v;

         v = j*NX + i;

         if(cc_est[v]>0)
         {

            d = sqrtf ( (i-i0)*(i-i0) + (j-j0)*(j-j0) );
            costheta = (i-i0)*hampel_axis[0]/d + (j-j0)*hampel_axis[1]/d ;
            if(costheta > 1.0 ) costheta=1.0;
            if(costheta < -1.0 ) costheta=-1.0;
            theta = acosf( costheta );

            if( theta <= pi/5.0) // C5 (color from line 37 of rainbow)
            {
               R[v] = 0;
               G[v] = 0;
               B[v] = 255;
            }
            else if(theta <= 2.0*pi/5.0 ) // C4 (color from line 74 of rainbow)
            {
               R[v] = 0;
               G[v] = 255;
               B[v] = 0;
            }
            else if(theta <= 3.0*pi/5.0 ) // C3 (color from line 74 of rainbow)
            {
               R[v] = 255;
               G[v] = 255;
               B[v] = 0;
            }
            else if(theta <= 4.0*pi/5.0 ) // C2 (color from line 74 of rainbow)
            {
               R[v] = 255;
               G[v] = 127;
               B[v] = 0;
            } 
            else if(theta <= pi)  // C1 (color from line 185 of rainbow)
            {
               R[v] = 255;
               G[v] = 0;
               B[v] = 0;
            }
            else
            {
               R[v] = 255;
               G[v] = 255;
               B[v] = 255;
            }
         }
         else
         {
            R[v] = (unsigned char)(trg[v]*255.0/max);
            G[v] = (unsigned char)(trg[v]*255.0/max);
            B[v] = (unsigned char)(trg[v]*255.0/max);
         }
      }

      int O[2];
      O[0] = (int)rintf( hampel_origin[0] );
      O[1] = (int)rintf( hampel_origin[1] );
      for(int i=O[0]-5; i<=O[0]+5; i++)
      {
         R[O[1]*NX+i] = 0;
         G[O[1]*NX+i] = 0;
         B[O[1]*NX+i] = 255;
      }

      for(int j=O[1]-5; j<=O[1]+5; j++)
      {
         R[j*NX+O[0]] = 0;
         G[j*NX+O[0]] = 0;
         B[j*NX+O[0]] = 255;
      }

      sprintf(outputfile,"%s_cc_hampel.ppm",prefix);
      save_as_ppm((const char *)outputfile, NX, NY, R, G, B);
   }
   //////////////////////////////////////////////////////////////////////////////////

   for(int v=0; v<NP; v++)
   {
      R[v] = (unsigned char)(trg[v]*255.0/max);
      G[v] = (unsigned char)(trg[v]*255.0/max);
      B[v] = (unsigned char)(trg[v]*255.0/max);
   }

   for(int i=1; i<NX-1; i++)
   for(int j=1; j<NY-1; j++)
   {
      if(cc_est[j*NX + i]>0 && (cc_est[j*NX+i-1]==0 || cc_est[j*NX+i+1]==0 || cc_est[(j-1)*NX+i]==0 || cc_est[(j+1)*NX+i]==0)  )
      {
         R[j*NX+i] = 255;
         G[j*NX+i] = 0;
         B[j*NX+i] = 0;
      }
   }

//
// The following two lines if enabled saves the MSP image with a red CC border
//
   if(opt_border)
   {
      //////////////////////////////////////////////////////////////////////////////
      // draws the boudning box in *cc.ppm image
      //////////////////////////////////////////////////////////////////////////////
      if(opt_box)
      {
         for(int i=UPPER_LEFT_i; i<=LOWER_RIGHT_i; i++)
         {
            for(int j=UPPER_LEFT_j; j<=UPPER_LEFT_j+1; j++)
            {
               R[j*NX+i] = 255;
               G[j*NX+i] = 255;
               B[j*NX+i] = 255;
            }

            for(int j=LOWER_RIGHT_j-1; j<=LOWER_RIGHT_j; j++)
            {
               R[j*NX+i] = 255;
               G[j*NX+i] = 255;
               B[j*NX+i] = 255;
            }
         }

         for(int j=UPPER_LEFT_j; j<=LOWER_RIGHT_j; j++)
         {
            for(int i=UPPER_LEFT_i; i<=UPPER_LEFT_i+1; i++)
            {
               R[j*NX+i] = 255;
               G[j*NX+i] = 255;
               B[j*NX+i] = 255;
            }

            for(int i=LOWER_RIGHT_i-1; i<=LOWER_RIGHT_i; i++)
            {
               R[j*NX+i] = 255;
               G[j*NX+i] = 255;
               B[j*NX+i] = 255;
            }
         }
      }
      /////////////////////////////////////////////////////////////////////
      
      sprintf(outputfile,"%s_cc_border.ppm",prefix);
      //save_as_ppm((const char *)outputfile, NX, NY, R, G, B);
   }

   free(R);
   free(G);
   free(B);

   return;
}

//////////////////////////////////////////////////////////////////////////////////

void output_bounding_box_ppm(short *trg, const char *prefix) 
{
   char outputfile[1024]="";
   int min=0, max=0;

   ////////////////////////////////////////////////////////////
   unsigned char *R, *G, *B;

   R = (unsigned char *)calloc(NP, 1);
   if(R == NULL)
   {
      memory_allocation_error("R");
   }

   G = (unsigned char *)calloc(NP, 1);
   if(G == NULL)
   {
      memory_allocation_error("G");
   }

   B = (unsigned char *)calloc(NP, 1);
   if(B == NULL)
   {
      memory_allocation_error("B");
   }

   setLowHigh(trg, NP, &min, &max);

   if(max==0) max=1; // to avoid division by zero

   for(int v=0; v<NP; v++)
   {
      if(trg[v]>max) 
      {
         R[v] = (unsigned char)(255.0);
         G[v] = (unsigned char)(255.0);
         B[v] = (unsigned char)(255.0);
      }
      else
      {
         R[v] = (unsigned char)(trg[v]*255.0/max);
         G[v] = (unsigned char)(trg[v]*255.0/max);
         B[v] = (unsigned char)(trg[v]*255.0/max);
      }
   }

   //////////////////////////////////////////////////////////////////////////////
   // draws the boudning box in *cc.ppm image
   //////////////////////////////////////////////////////////////////////////////
   if(opt_box)
   {
      for(int i=UPPER_LEFT_i; i<=LOWER_RIGHT_i; i++)
      {
         for(int j=UPPER_LEFT_j; j<=UPPER_LEFT_j+1; j++)
         {
            R[j*NX+i] = 255;
            G[j*NX+i] = 255;
            B[j*NX+i] = 255;
         }

         for(int j=LOWER_RIGHT_j-1; j<=LOWER_RIGHT_j; j++)
         {
            R[j*NX+i] = 255;
            G[j*NX+i] = 255;
            B[j*NX+i] = 255;
         }
      }

      for(int j=UPPER_LEFT_j; j<=LOWER_RIGHT_j; j++)
      {
         for(int i=UPPER_LEFT_i; i<=UPPER_LEFT_i+1; i++)
         {
            R[j*NX+i] = 255;
            G[j*NX+i] = 255;
            B[j*NX+i] = 255;
         }

         for(int i=LOWER_RIGHT_i-1; i<=LOWER_RIGHT_i; i++)
         {
            R[j*NX+i] = 255;
            G[j*NX+i] = 255;
            B[j*NX+i] = 255;
         }
      }
   }
   /////////////////////////////////////////////////////////////////////

   if( acpoint<5 || acpoint>(NX-5)) acpoint = NX/2;
   if( pcpoint<5 || pcpoint>(NX-5)) pcpoint = NX/2;

   // mark the AC and PC locations
   /*
   for(int j=(NY/2-1)-5; j<=(NY/2)+5; j++)
   {

      if(j==(NY/2-1) || j==(NY/2))
      {
         for(int i=acpoint-5; i<=acpoint+1+5; i++)
         {
            R[j*NX+i] = 0;
            G[j*NX+i] = 255;
            B[j*NX+i] = 0;
         }

         for(int i=pcpoint-5; i<=pcpoint+1+5; i++)
         {
            R[j*NX+i] = 255;
            G[j*NX+i] = 0;
            B[j*NX+i] = 0;
         }
      }
      else
      {
         for(int i=acpoint; i<=acpoint+1; i++)
         {
            R[j*NX+i] = 0;
            G[j*NX+i] = 255;
            B[j*NX+i] = 0;
         }

         for(int i=pcpoint; i<=pcpoint+1; i++)
         {
            R[j*NX+i] = 255;
            G[j*NX+i] = 0;
            B[j*NX+i] = 0;
         }
      }
   }
   */
 
   /////////////////////////////////////////////////////////////////////////////////////
   // mark the CC border 
   // upper and lower boundaries are colored differently
   /////////////////////////////////////////////////////////////////////////////////////
   for(int b=0; b<nub; b++)
   {
      //R[ ubj[b]*NX+ubi[b]] = 255;
      //G[ ubj[b]*NX+ubi[b]] = 255;
      //B[ ubj[b]*NX+ubi[b]] = 0;
      R[ ubj[b]*NX+ubi[b]] = 255;
      G[ ubj[b]*NX+ubi[b]] = 0;
      B[ ubj[b]*NX+ubi[b]] = 0;
      
      if(opt_border)
      {
         R[ ubj[b]*NX+ubi[b]] = 255;
         G[ ubj[b]*NX+ubi[b]] = 0;
         B[ ubj[b]*NX+ubi[b]] = 0;
      }
   }
   for(int b=0; b<nlb; b++)
   {
      //R[ lbj[b]*NX+lbi[b]] = 0;
      //G[ lbj[b]*NX+lbi[b]] = 255;
      //B[ lbj[b]*NX+lbi[b]] = 255;
      R[ lbj[b]*NX+lbi[b]] = 255;
      G[ lbj[b]*NX+lbi[b]] = 0;
      B[ lbj[b]*NX+lbi[b]] = 0;

      if(opt_border)
      {
         R[ lbj[b]*NX+lbi[b]] = 255;
         G[ lbj[b]*NX+lbi[b]] = 0;
         B[ lbj[b]*NX+lbi[b]] = 0;
      }
   }

   /////////////////////////////////////////////////////////////////////////////////////
   // Mark the medial axis as red

   /*
   if(!opt_border)
   for(int v=0; v<nm; v++)
   {
      R[ medj[v]*NX + medi[v]] = 255;
      G[ medj[v]*NX + medi[v]] = 0;
      B[ medj[v]*NX + medi[v]] = 0;
   }
   */

   //////////////////////////////////////////////////////////////////////
   // mark rostrum, most anterior and most posterior points etc.
   
   /*
   for(int i=rostrum[0]-3; i<=rostrum[0]+3; i++)
   if(i<NX && i>=0)
   {
      R[rostrum[1]*NX+i] = 0;
      G[rostrum[1]*NX+i] = 0;
      B[rostrum[1]*NX+i] = 255;
   }

   for(int j=rostrum[1]-3; j<=rostrum[1]+3; j++)
   if(j<NY && j>=0)
   {
      R[j*NX+rostrum[0]] = 0;
      G[j*NX+rostrum[0]] = 0;
      B[j*NX+rostrum[0]] = 255;
   }
   
   for(int i=anterior_point[0]-3; i<=anterior_point[0]+3; i++)
   if(i<NX && i>=0)
   {
      R[anterior_point[1]*NX+i] = 0;
      G[anterior_point[1]*NX+i] = 0;
      B[anterior_point[1]*NX+i] = 255;
   }

   for(int j=anterior_point[1]-3; j<=anterior_point[1]+3; j++)
   if(j<NY && j>=0)
   {
      R[j*NX+anterior_point[0]] = 0;
      G[j*NX+anterior_point[0]] = 0;
      B[j*NX+anterior_point[0]] = 255;
   }

   for(int i=posterior_point[0]-3; i<=posterior_point[0]+3; i++)
   if(i<NX && i>=0)
   {
      R[posterior_point[1]*NX+i] = 0;
      G[posterior_point[1]*NX+i] = 0;
      B[posterior_point[1]*NX+i] = 255;
   }

   for(int j=posterior_point[1]-3; j<=posterior_point[1]+3; j++)
   if(j<NY && j>=0)
   {
      R[j*NX+posterior_point[0]] = 0;
      G[j*NX+posterior_point[0]] = 0;
      B[j*NX+posterior_point[0]] = 255;
   }

   for(int i=inferior_genu[0]-3; i<=inferior_genu[0]+3; i++)
   if(i<NX && i>=0)
   {
      R[inferior_genu[1]*NX+i] = 0;
      G[inferior_genu[1]*NX+i] = 0;
      B[inferior_genu[1]*NX+i] = 255;
   }

   for(int j=inferior_genu[1]-3; j<=inferior_genu[1]+3; j++)
   if(j<NY && j>=0)
   {
      R[j*NX+inferior_genu[0]] = 0;
      G[j*NX+inferior_genu[0]] = 0;
      B[j*NX+inferior_genu[0]] = 255;
   }

   for(int i=inferior_splenium[0]-3; i<=inferior_splenium[0]+3; i++)
   if(i<NX && i>=0)
   {
      R[inferior_splenium[1]*NX+i] = 0;
      G[inferior_splenium[1]*NX+i] = 0;
      B[inferior_splenium[1]*NX+i] = 255;
   }

   for(int j=inferior_splenium[1]-3; j<=inferior_splenium[1]+3; j++)
   if(j<NY && j>=0)
   {
      R[j*NX+inferior_splenium[0]] = 0;
      G[j*NX+inferior_splenium[0]] = 0;
      B[j*NX+inferior_splenium[0]] = 255;
   }

   for(int i=posterior_genu[0]-3; i<=posterior_genu[0]+3; i++)
   if(i<NX && i>=0)
   {
      R[posterior_genu[1]*NX+i] = 0;
      G[posterior_genu[1]*NX+i] = 0;
      B[posterior_genu[1]*NX+i] = 255;
   }

   for(int j=posterior_genu[1]-3; j<=posterior_genu[1]+3; j++)
   if(j<NY && j>=0)
   {
      R[j*NX+posterior_genu[0]] = 0;
      G[j*NX+posterior_genu[0]] = 0;
      B[j*NX+posterior_genu[0]] = 255;
   }
   */
   //////////////////////////////////////////////////////////////////////

   sprintf(outputfile,"%s_cc.ppm",prefix);
   save_as_ppm((const char *)outputfile, NX, NY, R, G, B);

   free(R);
   free(G);
   free(B);

   return;
}

//////////////////////////////////////////////////////////////////////////////////

// NOTE: The commented out lines reflect secret options
void print_help()
{
   printf("\nUsage: yuki [optional arguments] -i <input volume>.nii\n\n"
 
   "Required argument:\n\n"

   "-i <input volume>.nii\n\t3D MRI volume on which the corpus callosum is to be segmented.\n"
   "\tThis image must be in NIFTI format of type \"short int\"\n\n"

   "Optional arguments:\n\n"

   "-o <output prefix>\n\tPrefix for naming output files (default: input volume prefix)\n\n"

   "-verbose (-v)\n\tEnables verbose mode\n\n"

   "-version (-V)\n\tReports software version\n\n"

   "-help (-h)\n\tPrints help message\n\n"

   "-n <integer>\n\tSpecifies the number of atlases to be used (default=49)\n\n"

   "-threshold (-t) <float>\n\tThreshold used for label fusion (default=50.0)\n\n"

   "-csv <csvfile>\n\tCC measurements (area, perimeter, etc.) will be appended to this file\n"
   "\tin comma-separated values (CSV) format (default: <output-prefix>.csv)\n\n"

   "-Hampel (-H)\n\tSegments the CC according to Hampel's method and outputs the 5 sub-areas\n"
   "\tas well as <output-prefix>_cc_hampel.ppm and <output-prefix>_cc_hampel.nii images\n\n"

   "-Witelson (-W)\n\tSegments the CC according to Witelson's method and outputs the 7 sub-areas\n"
   "\tas well as <output-prefix>_cc_witelson.ppm and <output-prefix>_cc_witelson.nii images\n\n"

   "-png\n\tOutputs *.png images in addition to the *.ppm images\n\n"

   "-lm <filename>\n\tManually specifies AC/PC/VSPS landmarks for <input-filename>.nii\n\n"

   "-atlas (-a) <filename>\n\tSpecifies the atlas to be used (default: amir464).\n"
   "\tThe other option is: babak628.\n\n"

   "-A <filename>\n\tUses preselected set of atlases specified in <filename> instead of\n"
   "\tautomated atlas selection. <filename> is always the output of a previous yuki run.\n\n"

   "-cc <corrected_cc.nii>\n\tThis option is used when the out binary CC image is corrected\n"
   "\tmanually and we need to recalculate the CC related measurements (area, circularity, etc.)\n"
   "\tfor the corrected image.\n\n"

   "-T <filename.mrx>\n\tApplies the transformation matrix in <filename.mrx> to reorient\n"
   "\tthe <input-filename>.nii volume in preparation for CC detection. Thus, automatic\n"
   "\treorientation is disabled.\n\n");

   return;
}

void print_secret_help()
{
   printf("\nUsage: yuki [-version -h -o <output-prefix> -csv <csvfile>] -i <subject volume>\n\n"
   "[-v -mrx -ppm -box -W -border -n <# atlases>]\n"
   "-mrx : saves the transformation matrix that makes the input volume MSP/AC-PC aligned\n\n"
   "-box : draws the CC search window on <prefix>_cc.ppm\n\n"
   "-border : outputs <prefix>_cc_border.ppm\n\n");
   return;
}

//////////////////////////////////////////////////////////////////////////////////
int vertical_runs(short *cc, int nx, int ny, int i0)
{
   int sign_change=0;
   int runs=0;

   for(int j=0; j<ny-1; j++)
   {
      if( cc[j*nx + i0] != cc[(j+1)*nx + i0] ) sign_change++;
   }

   runs = sign_change/2;

   return( runs );
}

//////////////////////////////////////////////////////////////////////////////////
void cw_eight_neighbor(short *im, int nx, int ny, int del_i, int del_j, int &bi, int &bj, int &ci, int &cj)
{
   if( im[ (cj+del_j)*nx + (ci+del_i) ] == 1 )
   {
      ci += del_i;
      cj += del_j;
      //printf("ci=%d cj=%d\n", ci, cj);
      return;
   }
   else
   {
      bi = ci + del_i;
      bj = cj + del_j;
      // printf("bi=%d bj=%d\n", bi, bj);
   }


   if( del_i==1 && del_j==0)
   {
      cw_eight_neighbor(im, nx, ny, 1, 1, bi, bj, ci, cj);
      return;
   }
   else if( del_i==1 && del_j==1)
   {
      cw_eight_neighbor(im, nx, ny, 0, 1, bi, bj, ci, cj);
      return;
   }
   else if( del_i==0 && del_j==1)
   {
      cw_eight_neighbor(im, nx, ny, -1, 1, bi, bj, ci, cj);
      return;
   }
   else if( del_i==-1 && del_j==1)
   {
      cw_eight_neighbor(im, nx, ny, -1, 0, bi, bj, ci, cj);
      return;
   }
   else if( del_i==-1 && del_j==0)
   {
      cw_eight_neighbor(im, nx, ny, -1, -1, bi, bj, ci, cj);
      return;
   }
   else if( del_i==-1 && del_j==-1)
   {
      cw_eight_neighbor(im, nx, ny, 0, -1, bi, bj, ci, cj);
      return;
   }
   else if( del_i==0 && del_j==-1)
   {
      cw_eight_neighbor(im, nx, ny, 1, -1, bi, bj, ci, cj);
      return;
   }
   else if( del_i==1 && del_j==-1)
   {
      cw_eight_neighbor(im, nx, ny, 1, 0, bi, bj, ci, cj);
      return;
   }

   return;
}

//////////////////////////////////////////////////////////////////////////////////

float estimate_circumference(short *msk, int nx, int ny, float dx, float dy, int *borderi, int *borderj, int &nb)
{
   int si, sj;
   int b;
   int indx;
   int ci, cj;
   int ni[8];
   int nj[8];
   int zero_image_flg=1;

   // zeros the first and last rows
   for(int i=0; i<nx; i++)
   {
      msk[  0*nx    + i ]=0;
      msk[(ny-1)*nx + i ]=0;
   }

   // zeros the first and last columns
   for(int j=0; j<ny; j++)
   {
      msk[ j*nx +   0    ]=0;
      msk[ j*nx + (nx-1) ]=0;
   }

   // find s (Step 1)
   for(int n=0; n<nx*ny; n++)
   {
      if( msk[n]>0 ) 
      {
         si = n%nx;
         sj = n/nx;
         zero_image_flg=0;
         break;
      }
   }

   if( zero_image_flg )
   {
      nb=0;
      return(0.0);
   }

   // initialize c (Step 2)
   ci=si;
   cj=sj;

   // initialize b (Step 2)
   b = 0;

   nb=0;
   do {

      borderi[nb]=ci;
      borderj[nb]=cj;
      nb++;

      // set n (Step 3)
      ni[0]=ci-1; nj[0]=cj;
      ni[1]=ci-1; nj[1]=cj-1;
      ni[2]=ci;   nj[2]=cj-1;
      ni[3]=ci+1; nj[3]=cj-1;
      ni[4]=ci+1; nj[4]=cj;
      ni[5]=ci+1; nj[5]=cj+1;
      ni[6]=ci;   nj[6]=cj+1;
      ni[7]=ci-1; nj[7]=cj+1;

      // initialize next c (Step 3)
      for(int k=0; k<8; k++)
      {
         indx = (b+k)%8;   
         if( msk[ nj[indx]*nx + ni[indx] ] > 0)
         {
            ci = ni[indx];
            cj = nj[indx];

            // very important derivation takes:
            // index =0 or 1 to b=6
            // index =2 or 3 to b=0 
            // index =4 or 5 to b=2 
            // index =6 or 7 to b=4
            b = ((indx/2 - 1)*2 + 8)%8;
            break;
         }
      }

   } while( ci!=si || cj!=sj);

   float circumference=0.0;
   float d1, d2;

   for(int k=1; k<nb; k++)
   {
      d1 = (borderi[k]-borderi[k-1])*dx;
      d2 = (borderj[k]-borderj[k-1])*dy;
      circumference += sqrtf( d1*d1 + d2*d2 );
   }
   d1 = (borderi[0]-borderi[nb-1])*dx;
   d2 = (borderj[0]-borderj[nb-1])*dy;
   circumference += sqrtf( d1*d1 + d2*d2 );

   return(circumference);
}
////////////////////////////////////////////////////////////////////

void estimate_witelson(short *msk, int nx, int ny, float dx, float dy, float *W)
{
   int v;

   for(int i=0; i<=7; i++)
   {
      W[i]=0;
   }

   for(int i=0; i<nx; i++)
   for(int j=0; j<ny; j++)
   {
      v = j*nx + i;
      if(msk[v]>0 && i>=(posterior_point[0]-cc_length_fifth) ) // W7
      {
         W[7] += dx*dy;
         msk[v]=7;
      }
      else if(msk[v]>0 && i>=(posterior_point[0]-cc_length_third) ) // W6
      {
         W[6] += dx*dy;
         msk[v]=6;
      }
      else if(msk[v]>0 && i>=(posterior_point[0]-cc_length_half) ) // W5
      {
         W[5] += dx*dy;
         msk[v]=5;
      }
      else if(msk[v]>0 && i>=(anterior_point[0]+cc_length_third) ) // W4
      {
         W[4] += dx*dy;
         msk[v]=4;
      }
      else if(msk[v]>0 && i<=posterior_genu[0] ) // W2
      {
         W[2] += dx*dy;
         msk[v]=2;
      }
      else if(msk[v]>0 && j>=posterior_genu[1] ) // W1
      {
         W[1] += dx*dy;
         msk[v]=1;
      }
      else if(msk[v]>0) // W3
      {
         W[3] += dx*dy;
         msk[v]=3;
      }
   }
}

void estimate_hampel(short *msk, int nx, int ny, float dx, float dy, float *H)
{
   float d;
   float pi;
   float i0, j0;
   float theta;
   float costheta=0.0;

   pi = 4.0*atanf(1.0);

   i0 = hampel_origin[0];
   j0 = hampel_origin[1];

   int v;

   for(int i=0; i<=5; i++)
   {
      H[i]=0;
   }

   for(int i=0; i<nx; i++)
   for(int j=0; j<ny; j++)
   {
      v = j*nx + i;

      d = sqrtf ( (i-i0)*(i-i0) + (j-j0)*(j-j0) );

      costheta = (i-i0)*hampel_axis[0]/d + (j-j0)*hampel_axis[1]/d ;
      if(costheta > 1.0 ) costheta=1.0;
      if(costheta < -1.0 ) costheta=-1.0;
      theta = acosf( costheta );

      if(msk[v]>0 && theta <= pi/5.0) // H5
      {
         H[5] += dx*dy;
         msk[v]=5;
      }
      else if(msk[v]>0 && theta <= 2.0*pi/5.0 ) // H4
      {
         H[4] += dx*dy;
         msk[v]=4;
      }
      else if(msk[v]>0 && theta <= 3.0*pi/5.0 ) // H3
      {
         H[3] += dx*dy;
         msk[v]=3;
      }
      else if(msk[v]>0 && theta <= 4.0*pi/5.0 ) // H2
      {
         H[2] += dx*dy;
         msk[v]=2;
      }
      else if(msk[v]>0 && theta <= pi ) // H1
      {
         H[1] += dx*dy;
         msk[v]=1;
      }
      else if(msk[v]>0) // H0
      {
         H[0] += dx*dy;
         msk[v]=0;
      }
   }
}
//////////////////////////////////////////////////////////////////////////////////

float estimate_area(short *msk, int nx, int ny, float dx, float dy)
{
   int ncc;
   int np;

   np = nx*ny;

   ncc=0;
   for(int i=0; i<np; i++)
   {
      if(msk[i]>0)
      {
         ncc++;
      }
   }

   return( ncc*dx*dy );
}
//////////////////////////////////////////////////////////////////////////////////

float compute_circularity(float area, float perimeter)
{
   float pi;
   float circularity=0.0;

   pi = 4.0*atanf(1.0);

   if(perimeter > 0.0)
   {
      circularity = 4.0*pi*area/(perimeter*perimeter);
   }

   return(circularity);
}

//////////////////////////////////////////////////////////////////////////////////

void find_posterior_genu(short *cc_est)
{
   int j0=0, j1=0;
   posterior_genu[0]=anterior_point[0];
   posterior_genu[1]=anterior_point[1];

   // starting a AC and going towards to posterior directions (i.e., decreasing i)
   for(int i=anterior_point[0]; i<=posterior_point[0]; i++)
   {
      // find the first i-position where there are two vertical runs
      if( vertical_runs(cc_est, NX, NY, i) == 2 )
      {
         posterior_genu[0] = i-1;

         for(int j=0; j<NY-1; j++)
         {
            if(cc_est[j*NX + i]>0 && cc_est[(j+1)*NX + i]==0) 
            {
               j0 = j;
               break;
            }
         }

         for(int j=NY-1; j>0; j--)
         {
            if(cc_est[j*NX + i]>0 && cc_est[(j-1)*NX + i]==0) 
            {
               j1 = j;
               break;
            }
         }

         break;
      }
   }

   posterior_genu[1] = (int)((j0+j1)/2.0);
}

//////////////////////////////////////////////////////////////////////////////////
void find_rostrum(short *cc_est)
{
   for(int i=NX/2; i>=0; i--)
   {

      // find the first i-position where there are two vertical runs
      if( vertical_runs(cc_est, NX, NY, i) == 2 )
      {
         int sign_changes=0;
         int j0=0, j1=0;
         for( int j=NY-1; j>0; j-- )
         {
            if( cc_est[j*NX + i] != cc_est[(j-1)*NX + i] )
            {
               sign_changes++;

               if(sign_changes == 1) 
               { 
                  j1 = j;
               }
               else if(sign_changes == 2) 
               { 
                  j0 = j-1;
                  break;
               }
            }
         }

         rostrum[0] = i;
         rostrum[1] = j1 - (j1-j0)/2;
         break;
      }
   }
}

//////////////////////////////////////////////////////////////////////////////////
void find_inferior_genu(short *cc_est)
{
   int i0=0, j0=0, v=0;

   for(int i=0; i<NX/2; i++)
   for(int j=0; j<NY; j++)
   {
      v = j*NX + i;

      if( cc_est[v] > 0 )
      {
         if(j>j0) { j0=j; i0=i;}
      } 
   }

   inferior_genu[0] = i0;
   inferior_genu[1] = j0;
}

void find_inferior_splenium(short *cc_est)
{
   int i0=0, j0=0, v=0;

   for(int i=NX/2; i<NX; i++)
   for(int j=0; j<NY; j++)
   {
      v = j*NX + i;

      if( cc_est[v] > 0 )
      {
         if(j>j0) { j0=j; i0=i;}
      } 
   }

   /////////////////////////////////////////////////////////////////////
   // Ensure if there is a run of CC pixels at the bottom of the 
   // splenium the middle pixel is chosen.
   /////////////////////////////////////////////////////////////////////
   int i1, i2;

   i1=i2=i0;

   for(int i=NX/2; i<NX-1; i++)
   {
      v = j0*NX + i;

      if( cc_est[v] > 0 && cc_est[v-1] == 0) 
      {
         i1 = i;
         break;
      }
   }

   for(int i=NX-1; i>NX/2; i--)
   {
      v = j0*NX + i;

      if( cc_est[v] > 0 && cc_est[v+1] == 0) 
      {
         i2 = i;
         break;
      }
   }

   i0 = (i1 + i2)/2;
   /////////////////////////////////////////////////////////////////////
   
   inferior_splenium[0] = i0;
   inferior_splenium[1] = j0;
}

void find_hampel_coordinates()
{
   float ui=0.0, uj=0.0, d=0.0;
   float p0i=0.0, p0j=0.0; 
   float p1i=0.0, p1j=0.0; 
   float ri=0.0, rj=0.0;

   p0i = inferior_genu[0];
   p0j = inferior_genu[1];
   p1i = inferior_splenium[0];
   p1j = inferior_splenium[1];

   d = sqrtf( (p0i-p1i)*(p0i-p1i) + (p0j-p1j)*(p0j-p1j) );

   if(d!=0.0)
   {
      ui = (p1i-p0i)/d;
      uj = (p1j-p0j)/d;
   }

   ri = posterior_point[0]-p1i;
   rj = posterior_point[1]-p1j;
   p1i += (ri*ui + rj*uj)*ui;
   p1j += (ri*ui + rj*uj)*uj;

   ri = anterior_point[0]-p0i;
   rj = anterior_point[1]-p0j;
   p0i += (ri*ui + rj*uj)*ui;
   p0j += (ri*ui + rj*uj)*uj;

   hampel_origin[0] = p0i + (p1i-p0i)/2.0;
   hampel_origin[1] = p0j + (p1j-p0j)/2.0;

   hampel_axis[0] = ui;
   hampel_axis[1] = uj;
}

//////////////////////////////////////////////////////////////////////////////////
// This function splits the CC boundary (cci, ccj) into two parts (ubi, ubj) and
// (lbi, lbj) representing the upper and lower boundaries if the CC.  
// The first point of the low boundary coincides with the last point of the upper boundary.
// The first point of the upper boundary coindicides with the last point of the low boudary.
//////////////////////////////////////////////////////////////////////////////////
void separate_upper_and_lower_boundaries(int *&ubi, int *&ubj, int *&lbi, int *&lbj, int *cci, int *ccj, int nb)
{
   int b0=0, b1=0;

   // b0 is the index of the boundry point corresponding to the inferior splenium 
   // b1 is the index of the boundry point corresponding to the rostrum 
   for(int b=0; b<nb; b++)
   {
      if(ccj[b]==rostrum[1] && cci[b]==rostrum[0])
         b1=b;

      if(ccj[b]==inferior_splenium[1] && cci[b]==inferior_splenium[0])
         b0=b;
   }

   // nlb: number of points on the lower boundary of the CC defined as going
   // from b0 to b1 clockwise
   nlb = b1-b0+1;

   // nub: number of points on the upper boundary of the CC defined as going
   // from b1 to b0 clockwise
   nub = (b0 + 1) + (nb - b1);

   // (ubi, ubj) : coordinates of points on the upper boundary (0, 1, ..., nub-1)
   ubi = (int *)calloc(nub, sizeof(int));
   ubj = (int *)calloc(nub, sizeof(int));

   // (lbi, lbj) : coordinates of points on the lower boundary (0, 1, ..., nlb-1)
   lbi = (int *)calloc(nlb, sizeof(int));
   lbj = (int *)calloc(nlb, sizeof(int));

   // fill in (lbi, lbj)_b (b=0, ..., nlb-1) 
   // from (cci, ccj)_b0 to (cci, ccj)_b1
   for(int b=0; b<nlb; b++)
   {
      lbi[b] = cci[b+b0];
      lbj[b] = ccj[b+b0];
   }

   // fill in (lbi, lbj)_b (b=0, ..., nb-b1-1) 
   // from (cci, ccj)_b1 to (cci, ccj)_(nb-1)
   for(int b=0; b<(nb-b1); b++)
   {
      ubi[b] = cci[b+b1];
      ubj[b] = ccj[b+b1];
   }

   // fill in (lbi, lbj)_b (b=nb-b1, ..., nub-1) 
   // from (cci, ccj)_0 to (cci, ccj)_b0
   for(int b=(nb-b1); b<nub; b++)
   {
      ubi[b] = cci[b-nb+b1];
      ubj[b] = ccj[b-nb+b1];
   }

   return;
}
//////////////////////////////////////////////////////////////////////////////////
void find_medial_point(int i, int j, short *cc)
{
   short mindist=32767;
   int ii=-1;
   int jj=-1;

   if(i<0 && j<0) return;

   medi[nm]=i;
   medj[nm]=j;
   nm++;

   cc[ j*NX + i ] = 2;

   if( i>0 && cc[ j*NX + i-1 ] == 1 )
   {
      cc[ j*NX + i-1 ] = 2;
      mindist = dist[ j*NX + i-1 ];
      ii=i-1; jj=j;
   }

   if( i>0 && j<NY-1 && cc[ (j+1)*NX + i-1 ] == 1 )
   {
      cc[ (j+1)*NX + i-1 ] = 2;
      if( dist[ (j+1)*NX + i-1 ] < mindist )
      {
         mindist = dist[ (j+1)*NX + i-1 ];
         ii=i-1; jj=j+1;
      }
   }

   if( j<NY-1 && cc[ (j+1)*NX + i ] == 1 )
   {
      cc[ (j+1)*NX + i ] = 2;
      if( dist[ (j+1)*NX + i ] < mindist )
      {
         mindist = dist[ (j+1)*NX + i ];
         ii=i; jj=j+1;
      }
   }

   if( i<NX-1 && j<NY-1 && cc[ (j+1)*NX + i+1 ] == 1 )
   {
      cc[ (j+1)*NX + i+1 ] = 2;
      if( dist[ (j+1)*NX + i+1 ] < mindist )
      {
         mindist = dist[ (j+1)*NX + i+1 ];
         ii=i+1; jj=j+1;
      }
   }

   if( i<NX-1 && cc[ j*NX + i+1 ] == 1 )
   {
      cc[ j*NX + i+1 ] = 2;
      if( dist[ j*NX + i+1 ] < mindist )
      {
         mindist = dist[ j*NX + i+1 ];
         ii=i+1; jj=j;
      }
   }

   if( i<NX-1 && j>0 && cc[ (j-1)*NX + i+1 ] == 1 )
   {
      cc[ (j-1)*NX + i+1 ] = 2;
      if( dist[ (j-1)*NX + i+1 ] < mindist )
      {
         mindist = dist[ (j-1)*NX + i+1 ];
         ii=i+1; jj=j-1;
      }
   }

   if( j>0 && cc[ (j-1)*NX + i ] == 1 )
   {
      cc[ (j-1)*NX + i ] = 2;
      if( dist[ (j-1)*NX + i ] < mindist )
      {
         mindist = dist[ (j-1)*NX + i ];
         ii=i; jj=j-1;
      }
   }

   if( i>0 && j>0 && cc[ (j-1)*NX + (i-1) ] == 1 )
   {
      cc[ (j-1)*NX + (i-1) ] = 2;
      if( dist[ (j-1)*NX + (i-1) ] < mindist )
      {
         mindist = dist[ (j-1)*NX + (i-1) ];
         ii=i-1; jj=j-1;
      }
   }

   find_medial_point(ii, jj, cc);

   return;
}

void find_medial_axis(short *cc)
{
   double dub, dlb;
   double dubmin, dlbmin;

   // (medi, medj) are coordinates of points on the medial axis
   medi= (int *)calloc(nub, sizeof(int));
   medj= (int *)calloc(nub, sizeof(int));

   // "distance" map pixel values equal the absolute value of the difference between
   // mininum distances of the pixel and the lower and upper CC boundaries
   dist = (short *)calloc(NX*NY, sizeof(short));

   // binarize cc (just in case)
   for(int v=0; v<NX*NY; v++) 
   if(cc[v]>0) { cc[v]=1; } else { cc[v]=0; }

   for(int i=0; i<NX; i++) 
   for(int j=0; j<NY; j++) 
   if( cc[j*NX + i] == 1)
   {
      dubmin = NX*dx+NY*dy; // can't be larger than this!
      for(int k=0; k<nub; k++)
      {
         dub = (ubi[k]-i)*dx*(ubi[k]-i)*dx + (ubj[k]-j)*dy*(ubj[k]-j)*dy;
         dub = sqrt(dub);
         if(dub<dubmin) dubmin=dub;
      }

      dlbmin = NX*dx+NY*dy; // can't be larger than this!
      for(int k=0; k<nlb; k++)
      {
         dlb = (lbi[k]-i)*dx*(lbi[k]-i)*dx + (lbj[k]-j)*dy*(lbj[k]-j)*dy;
         dlb = sqrt(dlb);
         if(dlb<dlbmin) dlbmin=dlb;
      }

      dist[j*NX + i]= (short)(100*fabs(dlbmin-dubmin));
   }
   else
   {
      dist[j*NX + i]= 0;
   }

   //uncomment to save "distance map"
   //save_nifti_image("distance.nii", dist, &output_hdr);

   //find_medial_point(inferior_splenium[0], inferior_splenium[1], cc);
   find_medial_point(rostrum[0], rostrum[1], cc);

   //////////////////////////////////////////////////////////////////////////////////
   // reverse (medi, medj) order.  This is because we now start from rostrum instead
   // of the inferior splenium.
   //////////////////////////////////////////////////////////////////////////////////
   int *dumi;
   int *dumj;
   dumi = (int *)calloc(nm, sizeof(int));
   dumj = (int *)calloc(nm, sizeof(int));

   for(int n=0; n<nm; n++)
   {
      dumi[n]=medi[n];
      dumj[n]=medj[n];
   }

   for(int n=0; n<nm; n++)
   {
      medi[n] = dumi[nm-1-n];
      medj[n] = dumj[nm-1-n];
   }
   free(dumi);
   free(dumj);
}

//////////////////////////////////////////////////////////////////////////////////
// nbp: number of boundary pixels. Equal nub, or nlb depending one whether we are
// interested in finding the upper or lower boundary crossing.
// (bi, bj) are boundary pixel locations
// (ui, uj) is unit normal poiting to the search direction.
// bp: index of the bp (in 1/2 increments)
// m: index of the medial axis pixel
void find_boundary_crossing(int i0, int j0, float &bp, short *cc, float ui, float uj, int *bi, int *bj, int nbp)
{
   float dpmax=-1.0;
   float dp;
   int indx;

   if(i0<1 || j0<1 || i0>=(NX-1) || j0>=(NY-1))
   {
      bp = -1.0;
      return;
   }

   if( cc[ j0*NX + i0 ] == 2 )
   {
      bp = -1.0;
      for(int b=0; b<nbp; b++)
      {
         if(bi[b]==i0 && bj[b]==j0)
         {
            bp = b;
            break;
         }
      }
      return;
   }
   else if( cc[ (j0+1)*NX + (i0-1) ]==1 && cc[ j0*NX + (i0-1) ]==2 && cc[ (j0+1)*NX + i0 ]==2 )
   {
      bp = -1.0;

      for(int b=0; b<nbp-1; b++)
      {
         if(bi[b]==i0 && bj[b]==(j0+1) )
         {
            bp = b;
            break;
         }
      }

      for(int b=0; b<nbp-1; b++)
      {
         if(bi[b]==(i0-1) && bj[b]==j0 )
         {
            bp += b;
            break;
         }
      }

      bp /= 2.0;
      return;
   }
   else if( cc[ (j0+1)*NX + (i0+1) ]==1 && cc[ j0*NX + (i0+1) ]==2 && cc[ (j0+1)*NX + i0 ]==2 )
   {
      bp = -1.0;

      for(int b=0; b<nbp-1; b++)
      {
         if(bi[b]==(i0+1) && bj[b]==j0 )
         {
            bp = b;
            break;
         }
      }

      for(int b=0; b<nbp-1; b++)
      {
         if(bi[b]==i0 && bj[b]==(j0+1) )
         {
            bp += b;
            break;
         }
      }

      bp /= 2.0;
      return;
   }
   else if( cc[ (j0-1)*NX + (i0+1) ]==1 && cc[ (j0-1)*NX + i0 ]==2 && cc[ j0*NX + (i0+1) ]==2 )
   {
      bp = -1.0;
      for(int b=0; b<nbp-1; b++)
      {
         if(bi[b]==i0 && bj[b]==(j0-1) )
         {
            bp = b;
            break;
         }
      }

      for(int b=0; b<nbp-1; b++)
      {
         if(bi[b]==(i0+1) && bj[b]==j0 )
         {
            bp += b;
            break;
         }
      }

      bp /= 2.0;
      return;
   }
   else if( cc[ (j0-1)*NX + (i0-1) ]==1 && cc[ (j0-1)*NX + i0 ]==2 && cc[ j0*NX + (i0-1) ]==2 )
   {
      bp = -1.0;
      for(int b=0; b<nbp-1; b++)
      {
         if(bi[b]==i0 && bj[b]==(j0-1) )
         {
            bp = b;
            break;
         }
      }

      for(int b=0; b<nbp-1; b++)
      {
         if(bi[b]==(i0-1) && bj[b]==j0 )
         {
            bp += b;
            break;
         }
      }

      bp /= 2.0;
      return;
   }

   // NW pixel
   dp = ui*(-0.7071068) + uj*(-0.7071068);
   dpmax = dp; indx=0;

   // N pixel
   dp = -uj;
   if( dp > dpmax )
   {
      dpmax = dp; indx=1;
   }

   // NE pixel
   dp = ui*(0.7071068) + uj*(-0.7071068);
   if( dp > dpmax )
   {
      dpmax = dp; indx=2;
   }

   // E pixel
   dp = ui;
   if( dp > dpmax )
   {
      dpmax = dp; indx=3;
   }

   // SE pixel
   dp = ui*(0.7071068) + uj*(0.7071068);
   if( dp > dpmax )
   {
      dpmax = dp; indx=4;
   }

   // S pixel
   dp = uj;
   if( dp > dpmax )
   {
      dpmax = dp; indx=5;
   }

   // SW pixel
   dp = ui*(-0.7071068) + uj*(0.7071068);
   if( dp > dpmax )
   {
      dpmax = dp; indx=6;
   }

   // W pixel
   dp = -ui;
   if( dp > dpmax )
   {
      dpmax = dp; indx=7;
   }

   switch (indx) {
      case 0:
         find_boundary_crossing(i0-1, j0-1, bp, cc, ui, uj, bi, bj, nbp);
         break;
      case 1:
         find_boundary_crossing(i0, j0-1, bp, cc, ui, uj, bi, bj, nbp);
         break;
      case 2:
         find_boundary_crossing(i0+1, j0-1, bp, cc, ui, uj, bi, bj, nbp);
         break;
      case 3:
         find_boundary_crossing(i0+1, j0, bp, cc, ui, uj, bi, bj, nbp);
         break;
      case 4:
         find_boundary_crossing(i0+1, j0+1, bp, cc, ui, uj, bi, bj, nbp);
         break;
      case 5:
         find_boundary_crossing(i0, j0+1, bp, cc, ui, uj, bi, bj, nbp);
         break;
      case 6:
         find_boundary_crossing(i0-1, j0+1, bp, cc, ui, uj, bi, bj, nbp);
         break;
      case 7:
         find_boundary_crossing(i0-1, j0, bp, cc, ui, uj, bi, bj, nbp);
         break;
   }

   return;
}

//////////////////////////////////////////////////////////////////////////////////

void find_thickness_profile(short *cc, const char *prefix)
{
   FILE *fp;
   char outputfile[1024]="";
   float dum=0.0;
   float *thickness_profile;
   short *arclength;
   float *tmp;

   ////////////////////////////////////////////////////////////////////
   // Find the normal vectors to the medial axis.
   // Program is designed such that the normal vectors always point 
   // towards the upper boundary of the CC.
   ////////////////////////////////////////////////////////////////////
   normi = (float *)calloc(nm, sizeof(float));
   normj = (float *)calloc(nm, sizeof(float));

   int m0, m1, dm=3;

   // find a normal vector for each pixel on the medial axis
   for(int m=0; m<nm; m++)
   {
      m0 = m - dm;
      m1 = m + dm;

      // take care of boundary conditions by reflection
      if( m0 < 0 ) m0 *= -1; 
      if( m1 >= nm ) m1 = 2*nm-2-m1;

      // at this point this is only the tangent line
      normi[m] = medi[m1] - medi[m0];
      normj[m] = medj[m1] - medj[m0];

      dum = sqrtf( normi[m]*normi[m] + normj[m]*normj[m] );

      if(dum>0.0)
      {
         normi[m] /= dum;
         normj[m] /= dum;
      }
      else
      {
         normi[m] = 0.0;
         normj[m] = 0.0;
      }

      // up to this point, normi and normj has been really the tanget line
      // here, it converted to the real normal line by a 90-deg rotation (i,j) --> (-j,i)
      dum = normi[m];
      normi[m] = -normj[m];
      normj[m] = dum;
   }
   ////////////////////////////////////////////////////////////////////

   thickness_profile = (float *)calloc(1001, sizeof(float));
   arclength = (short *)calloc(nm, sizeof(short));
   tmp = (float *)calloc(nm, sizeof(float));
   ////////////////////////////////////////////////////////////////////

   // indices of the nm boundary point where the nm normals to the
   // medial axis points cross the upper boundary
   // possible values are 0 to nub
   ubindx = (float *)calloc(nm, sizeof(float));

   // indices of the nm boundary point where the nm normals to the
   // medial axis points cross the lower boundary
   // possible values are 0 to nlb
   lbindx = (float *)calloc(nm, sizeof(float));

   // ensure that cc is 0 or 1
   for(int i=0; i<NX; i++)
   for(int j=0; j<NY; j++)
   if( cc[j*NX + i] > 0 )
   {
      cc[j*NX + i] = 1;
   }
   else
   {
      cc[j*NX + i] = 0;
   }

   // make the boundary points 2
   for(int b=0; b<nub; b++)
   {
      cc[ ubj[b]*NX + ubi[b] ] = 2;
   }

   for(int b=0; b<nlb; b++)
   {
      cc[ lbj[b]*NX + lbi[b] ] = 2;
   }

   for(int m=0; m<nm; m++)
   {
      find_boundary_crossing(medi[m], medj[m], ubindx[m], cc, normi[m], normj[m], ubi, ubj, nub);
      find_boundary_crossing(medi[m], medj[m], lbindx[m], cc, -normi[m], -normj[m], lbi, lbj, nlb);

      if( ubindx[m] < 0.0) ubindx[m]=0.0;
      if( lbindx[m] < 0.0) lbindx[m]=0.0;
   }

   float maxindx;
   float minindx;

   // ensure that ubindx values are in decreasing order
   // those points that violate this are set to zero
   maxindx=ubindx[0];
   for(int m=1; m<nm; m++)
   if( ubindx[m]>0.0 ) 
   {
      if( ubindx[m] <= maxindx )
      {
         maxindx = ubindx[m];
      }
      else
      {
         ubindx[m] = 0.0;
      }
   }

   // ensure that lbindx values are in increasing order
   // those points that violate this are set to zero
   minindx=lbindx[0];
   for(int m=1; m<nm; m++)
   if( lbindx[m]>0.0 ) 
   {
      if( lbindx[m] >= minindx )
      {
         minindx = lbindx[m];
      }
      else
      {
         lbindx[m] = 0.0;
      }
   }

   for(int m=1; m<(nm-1); m++)
   {
      if ( ubindx[m]==0.0 || lbindx[m]==0.0 )
      {
          ubindx[m]=lbindx[m]=0.0;
      }
   }

   // those points in ibindx and ubindx that were set to zero are filled
   // by interpolating between neighboring non-zero points
   // zero filing of lbindx and ubindx
   int count=0;
   for(int m=1; m<(nm-1); m++)
   {
      if ( ubindx[m]==0.0 && lbindx[m]==0.0 )
      {
         count=1;
         m0 = m;
         while( ubindx[m0 + count]==0.0 && lbindx[m0 + count]==0.0 )
         {
            count++;
         }
         m1 = m0+count;

         // interpolation
         for( int c=0; c<count; c++)
         {
            lbindx[m0+c] = lbindx[m0-1] + (lbindx[m1]-lbindx[m0-1])*(c+1)/(m1-m0+1);
            ubindx[m0+c] = ubindx[m0-1] + (ubindx[m1]-ubindx[m0-1])*(c+1)/(m1-m0+1);
         }

         // jump over the 00000 segment that was just processed
         m += count;
      }
   }

   tmp[0]=0.0;
   for(int m=1; m<nm; m++)
   {
      tmp[m] = tmp[m-1] + 
      sqrtf( (medi[m]-medi[m-1])*(medi[m]-medi[m-1])*dx*dx + (medj[m]-medj[m-1])*(medj[m]-medj[m-1])*dy*dy );
   }

   if( tmp[nm-1] > 0.0 )
   for(int m=0; m<nm; m++)
   {
      tmp[m] /= tmp[nm-1];
      arclength[m] = (short)(tmp[m]*10000.0 + 0.5);
   }

   float res;
   float lx, ly, ux, uy;
   for(int m=0; m<nm; m++)
   {
      res = lbindx[m] - (int)lbindx[m];

      if( (int)lbindx[m]< (nlb-1) )
      {
         lx = dx* ( lbi[ (int)lbindx[m] ]*(1.0-res) + lbi[ (int)lbindx[m] + 1]*res );
         ly = dy* ( lbj[ (int)lbindx[m] ]*(1.0-res) + lbj[ (int)lbindx[m] + 1]*res );
      }
      else
      {
         lx = dx*lbi[ (int)lbindx[m] ];
         ly = dy*lbj[ (int)lbindx[m] ];
      }

      res = ubindx[m] - (int)ubindx[m];

      if( (int)ubindx[m]< (nub-1) )
      {
         ux = dx * ( ubi[ (int)ubindx[m] ]*(1.0-res) + ubi[ (int)ubindx[m] + 1]*res );
         uy = dy * ( ubj[ (int)ubindx[m] ]*(1.0-res) + ubj[ (int)ubindx[m] + 1]*res );
      }
      else
      {
         ux = dx*ubi[ (int)ubindx[m] ];
         uy = dy*ubj[ (int)ubindx[m] ];
      }

      tmp[m] = sqrtf( (ux-lx)*(ux-lx) + (uy-ly)*(uy-ly) );
   }

   thickness_profile[0]=thickness_profile[1000]=0.0;
   for(int i=1; i<1000; i++)
   {
      for(int m=1; m<nm; m++)
      if( arclength[m] >= (i*10)  )
      {
         if( (arclength[m]- arclength[m-1]) > 0 )
            res = ((i*10.0) - arclength[m-1])/(arclength[m]- arclength[m-1]);
         else
            res = 0.0;
         thickness_profile[i] = tmp[m-1]*(1.0-res) + tmp[m]*res;
         break;
      }
   }

   free(tmp);

   tmp = smoothY(thickness_profile, 1, 1001, 1, 10.0);
   tmp[0]=tmp[1000]=0.0;

   sprintf(outputfile,"%s_TP.txt",prefix);
   fp = fopen(outputfile,"w");

   if(fp != NULL)
   {
      for(int i=0; i<1001; i+=10)
      {
         fprintf(fp,"%f\n",tmp[i]);
      }

      fclose(fp);
   }

   //////////////////////////////////////////////////////////////////////////////
   free(arclength);
   free(tmp);
   free(normj);
   free(normi);
   free(thickness_profile);
   free(ubindx);
   free(lbindx);
}

//////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  char ipimagepath[1024]="";  // important to initialize to "" 
  char ipimagedir[1024]="";     // important to initialize to ""

  int number_of_atlases_used=49;
  float max_t=50.0;

   char atlas_filename[1024]="amir464";
   char filename[1024]=""; // for keeping various filenames temporarily

   char lmfile[DEFAULT_STRING_LENGTH]="";

   FILE *fp;

   int kmin=0, kmax=100;

   float W[8]; // areas of Witelson's subdivisions are stored in W[1], ... W[7]; W[0] is unused
   float H[6]; // areas of Hampel's subdivisions are stored in H[1], ... H[5]; H[0] is unused
   float CCarea;
   float CCperimeter;
   float CCcircularity;

   float *avg_warped_cc=NULL;
   float *dumf=NULL;
   short *cc_est=NULL;
   short *warped_cc=NULL;

   int count=0;
   short *atlas_msp_ptr;

   char preselected_atlases_file[1024]="";
   char selected_atlases_file[1024]="";
   char csvfile[1024]="";
   char csvfilepath[1024]="";

   float sd;

   int N=11; // N=2*L+1

   char output_prefix[1024]="";
   char outputfile[1024]="";
   char msp_transformation_file[1024]="";

   int window_width=5;

   /////////////////////////////////////////////

   // by setting this option, the program will not output the *ACPC_axial.ppm and *ACPC_sagittal.ppm files
   // for the AC/PC detection program.
   opt_ppm=YES;
   opt_png=NO;
   opt_txt=YES;

   ////////////////////////////////////////////////////////////////////////
   if( argc==1 ) 
   {
      print_help();
      exit(0);
   }
   ////////////////////////////////////////////////////////////////////////

  while( (opt=getoption(argc, argv, options)) != -1)
  {
    switch (opt) 
    {
         case 's':
            print_secret_help();
            exit(0);
      case 'V':
        printf("Version 3.1 (Nov. 2024)\n");
        printf("Author: Babak A. Ardekani, Ph.D.\n");
        exit(0);
      case 'h':
        print_help();
        exit(0);
      case 'n':
        number_of_atlases_used=atoi(optarg);
        if(number_of_atlases_used<=0) number_of_atlases_used=49;
        break;
      case 'i':
        sprintf(ipimagepath,"%s",optarg);
        break;
      case 'C':
        opt_cc=YES;
        sprintf(ipimagepath,"%s",optarg);
        break;
      case 't':
        max_t  = atof(optarg);
        break;
      case 'o':
        sprintf(output_prefix,"%s",optarg);
        break;
      case 'c':
        sprintf(csvfile,"%s",optarg);
        break;
      case 'v':
        opt_v=YES;
        break;
      case 'W':
        opt_W=YES;
        break;
      case 'H':
        opt_H=YES;
        break;
      case 'p':
        opt_png=YES;
        break;
      case 'l':
        strcpy(lmfile,optarg);
        break;
      case 'a':
        sprintf(atlas_filename,"%s",optarg);
        break;
      case 'T':
        sprintf(msp_transformation_file,"%s",optarg);
        break;
      case 'L':
        sprintf(preselected_atlases_file,"%s",optarg);
        break;
      case 'b':
        opt_box=YES;
        break;
      case 'B':
        opt_border=YES;
        break;
      case '?':
        print_help();
        exit(0);
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  // Guard against nonsense values entered
  /////////////////////////////////////////////////////////////////////////////////////////////
  if(max_t > 100.0 || max_t < 0 ) max_t = 50.0;
  /////////////////////////////////////////////////////////////////////////////////////////////
  
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  // get the value of the ARTHOME environment variable
  // The getenv() function searches the environment list for a string that matches "ARTHOME".
  // It returns a pointer to the value in the environment, or NULL if there is no match.
  /////////////////////////////////////////////////////////////////////////////////////////////
  char *ARTHOME;  // full path of the directory of the ART software

  ARTHOME=getenv("ARTHOME");

  if(ARTHOME == NULL)
  {
    printf("The ARTHOME environment variable is not defined. Aborting ...\n");
    exit(0);
  }

  if(opt_v)
  {
    printf("ARTHOME = %s\n",ARTHOME);
  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////
  //Ensure than a subject filename has been specified at the command line
  /////////////////////////////////////////////////////////////////////////////////////////////
  if( ipimagepath[0]=='\0')
  {
    printf("Please specify an input volume using: -i <input volume>.nii\n");
    exit(0);
  }

  // determine input image directory
  getDirectoryName(ipimagepath, ipimagedir);

  if(opt_v)
  {
    printf("Input volume = %s\n",ipimagepath);
  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////
  //Determine the outut_prefix
  /////////////////////////////////////////////////////////////////////////////////////////////
  if( output_prefix[0]=='\0')
  {
    if( niftiFilename(output_prefix, ipimagepath)==0 ) { exit(1); }
  }

  if(opt_v)
  {
    printf("Output prefix = %s\n",output_prefix);
  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////
  //If using preselected atlases, read the number_of_atlases_used and max_t for now
  /////////////////////////////////////////////////////////////////////////////////////////////
  if(preselected_atlases_file[0]!='\0') 
  {
    int dum;

    if(opt_v) printf("Preselected atlases = %s\n", preselected_atlases_file);

    fp = fopen(preselected_atlases_file,"r");
    if(fp==NULL) file_open_error(preselected_atlases_file);
    if( fscanf(fp,"%d\n",&number_of_atlases_used) == EOF ) errorMessage("Read error, aborting ...");
    for(int i=0; i<number_of_atlases_used; i++) 
    {
      if( fscanf(fp,"%d\n",&dum) == EOF) errorMessage("Read error, aborting ...");
    }
    if( fscanf(fp,"%f\n",&max_t) == EOF ) errorMessage("Read error, aborting ...");
    fclose(fp);

    if(max_t<0.0 || max_t>100.0) max_t=50.0;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  //Read atlas_data from atlas_filename and set related variables
  /////////////////////////////////////////////////////////////////////////////////////////////
  short *atlas_data=NULL; // initially all contents of the atlas_filename are stored here
  short *atlas_msp=NULL;  // the gray scale MSP images contained in the atlas_filename
  short *atlas_cc=NULL;   // the binary CC segmentations contained in the atlas_filename

  nifti_1_header atlas_hdr;
  int bbnx, bbny, bbnp;
  int number_of_atlases_available;
  short *atlas_cc_ptr;
  short *bbtrg=NULL; // the bounding box MSP subimage of the subject image
  short *cctrg=NULL;

  sprintf(filename,"%s/%s.nii",ARTHOME, atlas_filename);

  if(opt_v)
  {
    printf("Atlas filename = %s\n", filename);
  }

  atlas_data = (short *)read_nifti_image(filename, &atlas_hdr);
  if(atlas_data == NULL)
  {
    printf("Error reading %s, aborting ...\n", filename);
    exit(0);
  }

  //This code section is for seting the bounding box when creating the atlas set
  //atlas_hdr.dim[4]=130;
  //atlas_hdr.dim[5]=141;
  //atlas_hdr.dim[6]=354;
  //atlas_hdr.dim[7]=279;
  //save_nifti_image("tt.nii", atlas_data, &atlas_hdr);
  //exit(0);

  UPPER_LEFT_i=atlas_hdr.dim[4];
  UPPER_LEFT_j=atlas_hdr.dim[5];
  LOWER_RIGHT_i=atlas_hdr.dim[6];
  LOWER_RIGHT_j=atlas_hdr.dim[7];

  bbnx =  LOWER_RIGHT_i - UPPER_LEFT_i + 1;
  bbny =  LOWER_RIGHT_j - UPPER_LEFT_j + 1;
  bbnp = bbnx*bbny;

  NP=NX*NY;
  number_of_atlases_available = atlas_hdr.dim[3]/2;
  dx=atlas_hdr.pixdim[1]; 
  dy=atlas_hdr.pixdim[2]; 
  dz=atlas_hdr.pixdim[3];

  bbtrg = (short *)calloc(bbnp, sizeof(short));
  cctrg = (short *)calloc(bbnp, sizeof(short));

  // ensures that the number of atlases used does no exceed the number of atlases available
  if(number_of_atlases_available < number_of_atlases_used)
  {
    number_of_atlases_used = number_of_atlases_available;
  }

  warped_cc = (short *)calloc(bbnp*number_of_atlases_available, sizeof(short));

  if(opt_v)
  {
    printf("Atlas matrix size = %d x %d (pixels)\n", NX, NY);
    printf("Atlas voxel size = %8.6f x %8.6f (mm^2)\n", dx, dy);
    printf("Number of atlases available = %d\n", number_of_atlases_available);
    printf("Number of atlases used = %d\n", number_of_atlases_used);
    printf("Bounding box size = %d x %d\n",bbnx,bbny);
    printf("Bounding box upper left corner = (%d, %d)\n",UPPER_LEFT_i,UPPER_LEFT_j);
    printf("Bounding box lower right corner = (%d, %d)\n",LOWER_RIGHT_i,LOWER_RIGHT_j);
  }

  atlas_msp = (short *)calloc(bbnp*number_of_atlases_available, sizeof(short));
  atlas_cc = (short *)calloc(bbnp*number_of_atlases_available, sizeof(short));

  // This is important!
  for(int a=0; a<number_of_atlases_available; a++)
  {
    atlas_msp_ptr  = atlas_data + (2*a)*bbnp;
    atlas_cc_ptr   = atlas_data + (2*a+1)*bbnp;

    for(int v=0; v<bbnp; v++)
    {
      atlas_msp[a*bbnp + v] = atlas_msp_ptr[v];

      //important: ensures that segmented CC values are 0 or 100
      if(atlas_cc_ptr[v] > 0) atlas_cc[a*bbnp + v] = 100; else atlas_cc[a*bbnp + v] = 0; 
    }
  }
  free(atlas_data);
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  // read the subject volume if specified at the command line
  /////////////////////////////////////////////////////////////////////////////////////////////
  if( ipimagepath[0] != '\0')
  {
    // check to see if ipimagepath appears to be a NIFTI image
    if( not_magical_nifti(ipimagepath) )
    {
      exit(0);
    }

    if( opt_v )
    {
      printf("Subject volume = %s\n",ipimagepath);
    }

    subj_volume=(short *)read_nifti_image(ipimagepath, &sub_hdr);

    if(subj_volume == NULL)
    {
      printf("Error reading %s, aborting ...\n", ipimagepath);
      exit(0);
    }

    Snx=sub_hdr.dim[1]; Sny=sub_hdr.dim[2]; Snz=sub_hdr.dim[3];
    Sdx=sub_hdr.pixdim[1]; Sdy=sub_hdr.pixdim[2]; Sdz=sub_hdr.pixdim[3];

    if(sub_hdr.datatype != DT_SIGNED_SHORT && sub_hdr.datatype != 512)
    {
      printf("\nSorry, this program only handles images of datatype\n"
      "DT_SIGNED_SHORT=4 or DT_UINT16=512. %s has datatype %d. Aborting ...\n\n", 
      ipimagepath, sub_hdr.datatype);
      free(subj_volume);
      exit(0);
    }

    if(opt_v)
    {
      printf("Subject volume matrix size = %d x %d x %d (voxels)\n", Snx, Sny, Snz);
      printf("Subject voxel size = %8.6f x %8.6f x %8.6f (mm3)\n", Sdx, Sdy, Sdz);
    }
  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  // subj_volume_msp image will be NX*NY*1 (i.e., 512*512*1)
  // This is the MSP of the subject image in PIL 
  /////////////////////////////////////////////////////////////////////////////////////////////
  short *subj_volume_msp = NULL;

  sprintf(outputfile,"%s/%s",ipimagedir, output_prefix);
  if( ipimagepath[0] != '\0')
  {
    if( msp_transformation_file[0]=='\0')
    {
      subj_volume_msp = find_subject_msp(ipimagepath, outputfile, lmfile);
    }
    else
    {
      subj_volume_msp = find_subject_msp_using_transformation(ipimagepath, outputfile, 
      msp_transformation_file);
    }

    count=0;
    for(int j=UPPER_LEFT_j; j<=LOWER_RIGHT_j; j++)
    for(int i=UPPER_LEFT_i; i<=LOWER_RIGHT_i; i++)
    {
      bbtrg[count] = subj_volume_msp[j*NX + i];
      count++;
    }

    //uncomment the next 6 lines to save bbtrg as a test
    //atlas_hdr.dim[1]=bbnx;
    //atlas_hdr.dim[2]=bbny;
    //atlas_hdr.dim[3]=1;
    //sprintf(outputfile,"bbtrg.nii");
    //save_nifti_image(outputfile, bbtrg, &atlas_hdr);
    //exit(0);
  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////
  //Altas selection
  /////////////////////////////////////////////////////////////////////////////////////////////
  float *corr=NULL;
  int *atlas_indx=NULL;

  if(!opt_cc)
  {
    corr=(float *)calloc(number_of_atlases_available,sizeof(float));
    atlas_indx = (int *)calloc(number_of_atlases_available, sizeof(int));

    for(int a=0; a<number_of_atlases_available; a++)
    {
      atlas_msp_ptr = atlas_msp + a*bbnp;

      atlas_indx[a] = a;
      corr[a] = pearsonCorrelation(bbtrg, atlas_msp_ptr, bbnp );
    }
    hpsort(number_of_atlases_available, corr, atlas_indx);

      // if atlases are preselected, read them into last elements of atlas_indx
      if(preselected_atlases_file[0] != '\0') 
      {
         fp = fopen(preselected_atlases_file,"r");
         if(fp==NULL) file_open_error(preselected_atlases_file);
         if( fscanf(fp,"%d\n",&number_of_atlases_used) == EOF ) errorMessage("Read error, aborting ...");
         for(int i=0; i<number_of_atlases_used; i++)
         {
            if( fscanf(fp,"%d\n",&atlas_indx[number_of_atlases_available-1-i]) == EOF ) 
	      errorMessage("Read error, aborting ...");
         }
         fclose(fp);
      }

    // saves the selected atlases
    sprintf(selected_atlases_file,"%s/%s_A.txt",ipimagedir, output_prefix);
    fp = fopen(selected_atlases_file, "w");
    if(fp==NULL) file_open_error(selected_atlases_file);
    fprintf(fp, "%d\n", number_of_atlases_used);
    for(int i=0; i<number_of_atlases_used; i++)
    {
      fprintf(fp, "%d\n", atlas_indx[number_of_atlases_available-1-i]);
    }
    fclose(fp);
  }  //if opt_cc
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  if(!opt_cc)
  {
    output_hdr = atlas_hdr;
    output_hdr.dim[3] = 1;
    output_hdr.dim[1] = NX;
    output_hdr.dim[2] = NY;
    sprintf(output_hdr.descrip,"Created by ART yuki");

    if(N<=0) N=11;
    if( (N%2)==0 ) N += 1; // ensuring it's odd
	
    if(window_width<=0) window_width=5;
    if( (window_width%2)==0 ) window_width+=1; // ensuring it's odd

    if( niter<=0 ) niter=4;
	
    // correlation window is (2*Lx+1)x(2*Ly+1)
    Lx = (N-1)/2;
    Ly = (int)(Lx*dx/dy + 0.5);
    N2 = (2*Lx+1)*(2*Ly+1);

    //Search window size is (2*Wx+1)*(2*Wy+1)
    Wx = (window_width-1)/2;
    Wy = (int)(Wx*dx/dy + 0.5);

    sd = 2.0*Wx;

    /////////////////////////////////////////////////////////////////////////////////////////
    // allocate memory for various arrays
    cc_est = (short *)calloc(NP, sizeof(short));
    dumf = (float *)calloc(bbnp, sizeof(float));
    Xwarp=(float *)calloc(bbnp,sizeof(float));
    Ywarp=(float *)calloc(bbnp,sizeof(float));
    Zwarp=(float *)calloc(bbnp,sizeof(float));
    ARobj=(float *)calloc(N2,sizeof(float));
    ARtrg=(float *)calloc(N2,sizeof(float));
    avg_warped_cc = (float *)calloc(bbnp, sizeof(float));

    /////////////////////////////////////////////////////////////////////////////////////////
    for(int i=0; i<number_of_atlases_used; i++)
    {
      int a;
      short *warped_cc_ptr;

      a = atlas_indx[number_of_atlases_available-1-i];
      atlas_msp_ptr = atlas_msp + a*bbnp;
      atlas_cc_ptr  = atlas_cc +  a*bbnp;

      computeWarpField(atlas_msp_ptr, bbtrg, sd, bbnp, bbnx, bbny);

      warped_cc_ptr=computeReslicedImage(atlas_cc_ptr, bbnx,bbny,1, dx,dy,dz, bbnx,bbny,1, dx,dy,dz, 
        Xwarp, Ywarp, Zwarp);

      // copy warped cc back to its original memory
      for(int v=0; v<bbnp; v++) 
      {
        avg_warped_cc[v] += warped_cc_ptr[v];
//        warped_cc[a*bbnp + v] = warped_cc_ptr[v];
        warped_cc[i*bbnp + v] = warped_cc_ptr[v];
      }

      free(warped_cc_ptr);

      if(opt_v)
      {
        printf("%d/%d: atlas #%03d corr=%6.4f\n", 
        i+1, number_of_atlases_used, a, corr[number_of_atlases_available-1-i]);
      }
    }

    for(int v=0; v<bbnp; v++) 
    {
      avg_warped_cc[v] /= number_of_atlases_used;
    }

    atlas_hdr.dim[1]=bbnx;
    atlas_hdr.dim[2]=bbny;
    atlas_hdr.dim[3]=number_of_atlases_used;
    sprintf(outputfile,"%s/%s_WACC.nii",ipimagedir,output_prefix);  //warped atlases corpora callosa
    save_nifti_image(outputfile, warped_cc, &atlas_hdr);

    ////////////////////////////////////////////////////////////
    // find kmin and kmax
    {
         float ccarea[101];
         int run[100];
         int vox;
         for(int k=0; k<=100; k++)
         {
            ccarea[k]=0.0;
            for(int i=1; i<bbnx-1; i++)
            for(int j=1; j<bbny-1; j++)
            {
               vox = j*bbnx + i;
               if( avg_warped_cc[vox] > (k*1.0) ) 
               {
                  ccarea[k] += dx*dy;
               }
            }
         }

         for(int k=0; k<100; k++)
         {
            //if( (ccarea[k]-ccarea[k+1]) < 10.0) run[k]=1; else run[k]=0;
            if( (ccarea[k]-ccarea[k+1]) < 5.0) run[k]=1; else run[k]=0;
         }

         int nruns=0;
         int run_start_index[100];
         int run_length[100];
         int max_run=0;

         for(int k=0; k<100; k++)
         if( run[k]==1 )
         {
            run_start_index[nruns]=k;
            run_length[nruns] = 0;
            while( run[k]==1 && k<100) 
            {  
               run_length[nruns]++;
               k++;
            }  
            nruns++;
         }

         for(int i=0; i<nruns; i++) 
         {
             if( run_length[i] > max_run )
             {
                max_run = run_length[i];
                kmin = run_start_index[i];
                kmax = run_start_index[i]+max_run-1;
             }
         }

         //if( opt_v )
         //   printf("Threshold selection interval: [%d, %d]\n",kmin, kmax);
    }

      
    {
      float mean1, mean2;
      float ssd;
      int n1, n2;
      float *fdr; // Fisher's discriminant ratio

      copyarray(avg_warped_cc, dumf, bbnp);

      fdr = (float *)calloc(101, sizeof(float));

         for(int k=kmin; k<=kmax; k++)
         {
            float t;
            t = k*1.0;

            n1=n2=0;
            mean1=mean2=0.0;
            ssd=0.0;
            for(int i=0; i<bbnx; i++)
            {
               for(int j=0; j<bbny; j++)
               {
                  int voxel = j*bbnx + i;

                  if( avg_warped_cc[voxel] > t && dumf[voxel] > 0.0) 
                  { 
                     n1++;
                     mean1 += bbtrg[voxel];
                  }
                  else if ( dumf[voxel] > 0.0)
                  {
                     n2++;
                     mean2 += bbtrg[voxel];
                  }
               } // j
            } // i

            if( n1>0 ) mean1 /= n1; 
            if( n2>0 ) mean2 /= n2; 
            if( n1==0 || n2==0 ) mean1=mean2=0.0;

            for(int i=0; i<bbnx; i++)
            {
               for(int j=0; j<bbny; j++)
               {
                  int voxel = j*bbnx + i;

                  if( avg_warped_cc[voxel] > t && dumf[voxel] > 0.0) 
                  { 
                     ssd += (mean1-bbtrg[voxel])*(mean1-bbtrg[voxel]);
                  }
                  else if ( dumf[voxel] > 0.0)
                  {
                     ssd += (mean2-bbtrg[voxel])*(mean2-bbtrg[voxel]);
                  }
               } // j
            } // i

         if( n1==0 || n2==0 ) fdr[k]=0.0;

         if( (n1+n2-2.0) > 0 ) ssd /= (n1+n2-2.0);

         fdr[k] = sqrtf( (mean1 - mean2)*(mean1 - mean2) / ssd);

        // printf("%f %d %d %d mean1 = %f mean2 = %f fdr=%f\n", t, n1, n2, n1+n2, mean1, mean2, fdr[k]);
      }

      int maxidx=0;
      float maxfdr=0.0;
      for(int k=kmin; k<=kmax; k++)
      {
        if(fdr[k] > maxfdr) 
        {
          maxidx=k;
          maxfdr = fdr[k];
        }
      }
      free(fdr);

      //if(opt_v)
      //{
      //   printf("maxfdr=%f maxidx=%d\n",maxfdr, maxidx);
      //}
      if(max_t == 0.0)
      {
        max_t = maxidx*1.0;
      }

      if( opt_v )
      {
         printf("Label fusion threshold = %3.1f\n",max_t);
      }

      // save the threshold used for label fusion
      fp = fopen(selected_atlases_file, "a");
      if(fp==NULL) file_open_error(selected_atlases_file);
      fprintf(fp, "%f\n", max_t);
      fclose(fp);

      for(int i=0; i<bbnx; i++)
      {
        for(int j=0; j<bbny; j++)
        {
          int voxel = j*bbnx + i;

          if( avg_warped_cc[voxel] > max_t)
          { 
            cc_est[ (j+UPPER_LEFT_j)*NX + (i+UPPER_LEFT_i)]=1;
            cctrg[voxel]=1;
          }
          else 
          {
            cc_est[ (j+UPPER_LEFT_j)*NX + (i+UPPER_LEFT_i)]=0;
            cctrg[voxel]=0;
          }
        } // j
      } // i

    }
  }

  if(opt_cc)
  {
    cc_est=(short *)read_nifti_image(ipimagepath, &output_hdr);
    dx=output_hdr.pixdim[1]; 
    dy=output_hdr.pixdim[1]; 
  }
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////
   {
      // locates the most anterior point of the CC
      for(int i=UPPER_LEFT_i; i<=LOWER_RIGHT_i; i++)
      {
         if( vertical_runs(cc_est, NX, NY, i) )
         {
            int j0, j1; 

            j1 = 0; j0 = NY;
            for(int j=UPPER_LEFT_j; j<=LOWER_RIGHT_j; j++)
            {
               if(cc_est[j*NX + i])
               {
                  if( j>j1 ) j1=j;
                  if( j<j0 ) j0=j;
               }
            }

            anterior_point[0]=i;
            anterior_point[1]=(j1+j0)/2;
         
            break;
         }
      }

      // locates the most posterior point of the CC
      for(int i=LOWER_RIGHT_i; i>=UPPER_LEFT_i; i--)
      {
         if( vertical_runs(cc_est, NX, NY, i) )
         {
            int j0, j1; 

            j1 = 0; j0 = NY;
            for(int j=UPPER_LEFT_j; j<=LOWER_RIGHT_j; j++)
            {
               if(cc_est[j*NX + i])
               {
                  if( j>j1 ) j1=j;
                  if( j<j0 ) j0=j;
               }
            }

            posterior_point[0]=i;
            posterior_point[1]=(j1+j0)/2;
         
            break;
         }
      }

      // locates the most posterior point of the genu (not rostrum)
      find_posterior_genu(cc_est); // required for Witelson

      find_inferior_genu(cc_est); // required for Hampel
      find_inferior_splenium(cc_est); // required for Hampel
      find_hampel_coordinates();
      find_rostrum(cc_est);

      // fix a rare degenerate case
      if (cc_est[ posterior_point[1]*NX + posterior_point[0] ]==0)
      {
         cc_est[ posterior_point[1]*NX + posterior_point[0] ]=1;
      }

      cc_length = posterior_point[0]-anterior_point[0];
      cc_length_fifth = (int)nearbyint( cc_length/5.0 );
      cc_length_third = (int)nearbyint( cc_length/3.0 );
      cc_length_half = (int)nearbyint( cc_length/2.0 );

/*
      if(opt_v)
      {
         printf("Most anterior point: %d, %d\n", anterior_point[0], anterior_point[1]);
         printf("Most posterior point: %d, %d\n", posterior_point[0], posterior_point[1]);
         printf("CC length = %f mm\n", cc_length*dx);
         printf("CC length / 5 = %d pixels \n", cc_length_fifth);
      }
*/
   }

   //////////////////////////////////////////////////////////////////////////
   {
      if( ipimagepath[0]!='\0')
      {
         sprintf(outputfile,"%s/%s",ipimagedir, output_prefix);
         output_ppm(subj_volume_msp, cc_est, (const char *)outputfile);
      }

      output_hdr.pixdim[4]=ACi;
      output_hdr.pixdim[5]=PCi;
      output_hdr.pixdim[6]=ACx;
      output_hdr.pixdim[7]=PCx;
      output_hdr.pixdim[1]=0.5;
      output_hdr.pixdim[2]=0.5;
      output_hdr.pixdim[3]=1.0;
      output_hdr.dim[1]=NX;
      output_hdr.dim[2]=NY;
      output_hdr.dim[3]=1;
      output_hdr.datatype=16;

      sprintf(outputfile,"%s/%s_cc.nii",ipimagedir,output_prefix);
      save_nifti_image(outputfile, cc_est, &output_hdr);

      update_qsform( (const char *)outputfile, Tacpc );
   }
   
   {
      CCarea = estimate_area(cc_est, NX, NY, dx, dy);

      cci = (int *)calloc(bbnp,sizeof(int));
      ccj = (int *)calloc(bbnp,sizeof(int));

      CCperimeter = estimate_circumference(cc_est, NX, NY, dx, dy, cci, ccj, nb);

      // separates the CC boundary (cci, ccj) into an upper part (ubi, ubj)
      // and a lower part (lbi, lbj)
      separate_upper_and_lower_boundaries(ubi, ubj, lbi, lbj, cci, ccj, nb);

      // This was failing for some cases so it is removed in version 3.1
      // finds (medi, medj) representing a CC medial axis
      // find_medial_axis(cc_est);
      
      CCcircularity = compute_circularity(CCarea, CCperimeter);

    sprintf(outputfile,"%s/%s",ipimagedir, output_prefix);
    find_thickness_profile(cc_est, outputfile);

    if( ipimagepath[0]!='\0')
    {
      output_bounding_box_ppm(subj_volume_msp, (const char *)outputfile);
    }

    if(opt_W)
    {
      estimate_witelson(cc_est, NX, NY, dx, dy, W);

      if( ipimagepath[0]!='\0')
      {
        output_hdr.pixdim[4]=ACi;
        output_hdr.pixdim[5]=PCi;
        output_hdr.pixdim[6]=ACx;
        output_hdr.pixdim[7]=PCx;
        output_hdr.pixdim[1]=0.5;
        output_hdr.pixdim[2]=0.5; 
        output_hdr.pixdim[3]=1.0;
        output_hdr.dim[1]=NX; 
        output_hdr.dim[2]=NY; 
        output_hdr.dim[3]=1;
        output_hdr.datatype=16;

        sprintf(outputfile,"%s/%s_cc_witelson.nii", ipimagedir, output_prefix);
        save_nifti_image(outputfile, cc_est, &output_hdr);

        update_qsform( (const char *)outputfile, Tacpc );
      }
    }

    if(opt_H)
    {
      estimate_hampel(cc_est, NX, NY, dx, dy, H);

      if( ipimagepath[0]!='\0')
      {
        output_hdr.pixdim[4]=ACi;
        output_hdr.pixdim[5]=PCi;
        output_hdr.pixdim[6]=ACx;
        output_hdr.pixdim[7]=PCx;
        output_hdr.pixdim[1]=0.5;
        output_hdr.pixdim[2]=0.5; 
        output_hdr.pixdim[3]=1.0;
        output_hdr.dim[1]=NX; 
        output_hdr.dim[2]=NY; 
        output_hdr.dim[3]=1;
        output_hdr.datatype=16;

        sprintf(outputfile,"%s/%s_cc_hampel.nii", ipimagedir, output_prefix);
        save_nifti_image(outputfile, cc_est, &output_hdr);

        update_qsform( (const char *)outputfile, Tacpc );
      }
    }

    if(csvfile[0]=='\0')
    {
      sprintf(csvfile,"%s.csv",output_prefix);
    }
    sprintf(csvfilepath,"%s/%s",ipimagedir,csvfile);

      if(opt_v)
      {
         printf("CC area = %6.2f mm^2\n",CCarea);
         printf("CC perimeter = %6.2f mm\n",CCperimeter);
         printf("CC circularity = %8.6f\n",CCcircularity);
         printf("CC length = %5.1f mm\n",cc_length*dx);
         if(opt_W)
         {
            printf("Witelson's subdivisoins:\n");
            printf("\tW1 = %6.2f mm^2\n",W[1]);
            printf("\tW2 = %6.2f mm^2\n",W[2]);
            printf("\tW3 = %6.2f mm^2\n",W[3]);
            printf("\tW4 = %6.2f mm^2\n",W[4]);
            printf("\tW5 = %6.2f mm^2\n",W[5]);
            printf("\tW6 = %6.2f mm^2\n",W[6]);
            printf("\tW7 = %6.2f mm^2\n",W[7]);
         }
         if(opt_H)
         {
            printf("Hampels's subdivisoins:\n");
            printf("\tC1 = %6.2f mm^2\n",H[1]);
            printf("\tC2 = %6.2f mm^2\n",H[2]);
            printf("\tC3 = %6.2f mm^2\n",H[3]);
            printf("\tC4 = %6.2f mm^2\n",H[4]);
            printf("\tC5 = %6.2f mm^2\n",H[5]);
         }
      }

      if( csvfilepath[0]!='\0')
      {
         FILE *fp;

         if (checkFileExistence(csvfilepath)==0)
         {
            fp = fopen(csvfilepath,"a");
            if(fp==NULL) file_open_error(csvfilepath);
            fprintf(fp,"\"ID\", CC_area, CC_perimeter, CC_circularity, CC_length");
            if(opt_W) 
            {
               fprintf(fp,", W1, W2, W3, W4, W5, W6, W7");
            }
            if(opt_H) 
            {
               fprintf(fp,", C1, C2, C3, C4, C5");
            }
            fprintf(fp,"\n");
            fclose(fp);
         }

         fp = fopen(csvfilepath,"a");
         if(fp==NULL) file_open_error(csvfilepath);
         if(opt_cc) fprintf(fp,"\"%s\", ",ipimagepath); 
         else fprintf(fp,"\"%s_cc.nii\", ",output_prefix);
         fprintf(fp,"%6.2f, ",CCarea);
         fprintf(fp,"%6.2f, ",CCperimeter);
         fprintf(fp,"%8.6f, ",CCcircularity);
         fprintf(fp,"%5.1f",cc_length*dx);
         if(opt_W)
         {
            fprintf(fp,", %6.2f, ",W[1]);
            fprintf(fp,"%6.2f, ",W[2]);
            fprintf(fp,"%6.2f, ",W[3]);
            fprintf(fp,"%6.2f, ",W[4]);
            fprintf(fp,"%6.2f, ",W[5]);
            fprintf(fp,"%6.2f, ",W[6]);
            fprintf(fp,"%6.2f",W[7]);
         }
         if(opt_H)
         {
            fprintf(fp,", %6.2f, ",H[1]);
            fprintf(fp,"%6.2f, ",H[2]);
            fprintf(fp,"%6.2f, ",H[3]);
            fprintf(fp,"%6.2f, ",H[4]);
            fprintf(fp,"%6.2f",H[5]);
         }
         fprintf(fp,"\n");
         fclose(fp);
      }
   }

   ////////////////////////////////////////////////////////////
   // Free Memory
   ////////////////////////////////////////////////////////////
   free(warped_cc);

   if(bbtrg != NULL) free(bbtrg);
   if(cctrg != NULL) free(cctrg);

   if( ipimagepath[0]!='\0')
   {
      if(subj_volume != NULL) free(subj_volume);
      if(subj_volume != NULL) free(subj_volume_msp);
   }

   if(!opt_cc)
   {
      if( corr != NULL) free(corr);
      if( atlas_indx != NULL) free(atlas_indx);
      if( dumf != NULL) free(dumf);

      if( Xwarp != NULL ) free(Xwarp);
      if( Ywarp != NULL ) free(Ywarp);
      if( Zwarp != NULL ) free(Zwarp);
      if( ARobj != NULL ) free(ARobj);
      if( ARtrg != NULL ) free(ARtrg);
      if( avg_warped_cc != NULL ) free(avg_warped_cc);
   }

   if(cc_est != NULL) free(cc_est);
   ////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////
}
