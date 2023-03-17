#define _PILTRANSFORM

#include <babak_lib.h>
#include <sph.h>
#include <ctype.h>
#include <landmarks.h>
#include <stats.h>
#include <string.h>

double searchradius[3]={50.0, 15.0, 15.0}; // in units of mm
int opt_CENTER_AC=NO;

void transform_P(float *P, int nl, float *T)
{
   float dum[4];

   for(int i=0; i<nl; i++)
   {
      dum[0]=P[0*nl + i];
      dum[1]=P[1*nl + i];
      dum[2]=P[2*nl + i];
      dum[3]=1.0;
      multi(T,4,4,dum,4,1,dum);
      P[0*nl + i]=dum[0];
      P[1*nl + i]=dum[1];
      P[2*nl + i]=dum[2];
   }

   return;
}

// This function makes a ppm file from the MSP and displays 'nl' detected landmarks.
//im: the input image 
//lmx and lmy: are pointers to nl size arrays including the (x,y) for the nl landmarks
//ppmfile: the output file name
void mspPPM(SHORTIM im, int *ii, int *jj, int nl, const char *ppmfile)
{
   int colourflag;
//   unsigned char colour[3]={0xFF,0xFF,0x00};
//   unsigned char colour[3]={0x00,0xFF,0x00};
   unsigned char colour[3]={0xFF,0xFF,0xFF};
   unsigned char *imgTemp;
   FILE *fp;
   int system_return_value;

   int kk=im.nz/2;
   
   /*size of the mark*/
   int d=4; 

   float temp;

   imgTemp=(unsigned char *)calloc(im.nv,sizeof(unsigned char));

   int low, high;
   setLowHigh(im.v+kk*im.np, im.nx*im.ny, &low, &high, 1.0);

   for(int i=0;i<im.nv;i++)
   {
      if( im.v[i] >= high ) temp=255.0;
      else temp=im.v[i]*255.0/high;

      imgTemp[i]=(unsigned char)temp;
   }

   fp=fopen(ppmfile,"w");
   if(fp==NULL) file_open_error(ppmfile);

   /*Write the header part of the PPM file*/
   fprintf(fp,"P6\n");
   fprintf(fp,"%d %d\n", im.nx, im.ny);
   fprintf(fp,"255\n");

   for (int j=0;j<im.ny;j++)
   {
      for(int i=0;i<im.nx;i++)
      {
         for(int m=0;m<nl;m++)
         {
//uncomment to save in different colours
//            if(m==3) { colour[0]=255; colour[1]=0; colour[2]=0; } 
//            else if(m==1) { colour[0]=0; colour[1]=255; colour[2]=0; } 
//            else if(m==2) { colour[0]=0; colour[1]=0; colour[2]=255; } 
//            else if(m==0) { colour[0]=0; colour[1]=255; colour[2]=191; } 
//            else if(m==4) { colour[0]=0; colour[1]=128; colour[2]=255; } 
//            else if(m==5) { colour[0]=255; colour[1]=0; colour[2]=255; } 
//            else if(m==6) { colour[0]=255; colour[1]=128; colour[2]=0; } 
//            else if(m==7) { colour[0]=255; colour[1]=255; colour[2]=0; } 
            
            if( (i==ii[m] && jj[m]-d<j && j<jj[m]+d) || (j==jj[m] && ii[m]-d<i && i<ii[m]+d) )
            {
               fwrite(colour,1,3,fp);
               colourflag=1;
            }
         }
         if(colourflag==0) 
         {
            fwrite(imgTemp+im.np*kk+im.nx*j+i,1,1,fp);
            fwrite(imgTemp+im.np*kk+im.nx*j+i,1,1,fp);
            fwrite(imgTemp+im.np*kk+im.nx*j+i,1,1,fp);
         }
         colourflag=0;
      }
   }

   fclose(fp);

  if(opt_png)
  {
    char *pngfilename;
    char cmnd[DEFAULT_STRING_LENGTH];
    int L;

    L = strlen(ppmfile);

    pngfilename = (char *)calloc(L+1,sizeof(char));
    stpcpy(pngfilename, ppmfile);
    pngfilename[L-1]='g';
    pngfilename[L-2]='n';

    snprintf(cmnd,sizeof(cmnd),"pnmtopng %s > %s",ppmfile,pngfilename); 
    if(opt_png) 
    {
        system_return_value = system(cmnd);
        if( system_return_value == -1 || system_return_value == 127 )
        {
           printf("Execution of command: %s failed\n",cmnd);
        }
    }

    free(pngfilename);
  }

  if(opt_ppm==NO) remove(ppmfile);

  delete imgTemp;
}

void convert_to_xyz(float *P, int n, SHORTIM im)
{
   /////////////////////////////////////////////////////////
   // convert P from (i,j,k) to (x,y,z) coordinates
   float T[16];
   float dum[4]={0.0,0.0,0.0,1.0};

   DIM dim;
   set_dim(dim,im);
   ijk2xyz(T, dim);

   for(int i=0; i<n; i++)
   {
      // note that P is (3 x n)
      dum[0]=P[0*n + i]; 
      dum[1]=P[1*n + i]; 
      dum[2]=P[2*n + i];
      dum[3]=1.0;
      multi(T,4,4,dum,4,1,dum);
      P[0*n + i]=dum[0]; 
      P[1*n + i]=dum[1]; 
      P[2*n + i]=dum[2];
   }
}

void convert_to_ijk(float *P, int n, SHORTIM im)
{
   /////////////////////////////////////////////////////////
   // convert P from (i,j,k) to (x,y,z) coordinates
   float T[16];
   float dum[4]={0.0,0.0,0.0,1.0};

   DIM dim;
   set_dim(dim,im);
   xyz2ijk(T, dim);

   for(int i=0; i<n; i++)
   {
      // note that P is (3 x n)
      dum[0]=P[0*n + i]; 
      dum[1]=P[1*n + i]; 
      dum[2]=P[2*n + i];
      dum[3]=1.0;
      multi(T,4,4,dum,4,1,dum);
      P[0*n + i]=dum[0]; 
      P[1*n + i]=dum[1]; 
      P[2*n + i]=dum[2];
   }
}

// matches P2=(x2,y2) as closely as possible to P1=(x1,y1) 
void point_match(float *x1, float *y1, float *x2, float *y2, int N, float *T)
{
   float params[3];
   float meanx1=0,meanx2=0;
   float meany1=0,meany2=0;
   float num=0,den=0;
   float *u1;
   float *v1;
   float *u2;
   float *v2;

   int i;

   set_to_I(T,4);

   for(i=0;i<N;i++)
   {
      meanx1+=x1[i];
      meany1+=y1[i];
      meanx2+=x2[i];
      meany2+=y2[i];
   }

   meanx1/=N;
   meany1/=N;
   meanx2/=N;
   meany2/=N;

   u1=(float *)calloc(N,sizeof(float));
   v1=(float *)calloc(N,sizeof(float));
   u2=(float *)calloc(N,sizeof(float));
   v2=(float *)calloc(N,sizeof(float));
   for(i=0;i<N;i++)
   {
      u1[i]=x1[i]-meanx1;
      v1[i]=y1[i]-meany1;
      u2[i]=x2[i]-meanx2;
      v2[i]=y2[i]-meany2;
   }

   for(i=0;i<N;i++)
   {
      num+=v1[i]*u2[i]-u1[i]*v2[i];
      den+=u1[i]*u2[i]+v1[i]*v2[i];
   }

   params[0]=atan2(num,den);
   params[1]=meanx1-cos(params[0])*meanx2+sin(params[0])*meany2;
   params[2]=meany1-sin(params[0])*meanx2-cos(params[0])*meany2;

   T[0] = cosf(params[0]);
   T[1] = -sinf(params[0]);
   T[4] = sinf(params[0]);
   T[5] = cosf(params[0]);
   T[3]=params[1];
   T[7]=params[2];

   //printMatrix(T,4,4,"T:",NULL);

   free(u1);
   free(v1);
   free(u2);
   free(v2);

   return;
}

void Procrustes(float *Q, float *Qavg, int n, float *P, float *Pavg, float *TLM)
{
   float Ht[9]; // H=P*Q' Ht=Q*P' (' means transpose in my notation)
   //float Ut[9], V[9], I[9];
   float Ut[9], V[9];
   float T[3]; // 3x1 translation vector
   float R[9]; // 3x3 rotation matrix 
   float S[3]; // 3x1 vector of singular values

   Pavg[0] = (float)removeVectorMean(P, n);
   Pavg[1] = (float)removeVectorMean(P+n, n);
   Pavg[2] = (float)removeVectorMean(P+2*n, n);

   Qavg[0] = (float)removeVectorMean(Q, n);
   Qavg[1] = (float)removeVectorMean(Q+n, n);
   Qavg[2] = (float)removeVectorMean(Q+2*n, n);

   mat_mat_trans(Q,3,n,P,3,Ht);

      svd(Ht, 3, 3, Ut, V, S);

      multi(V,3,3,Ut,3,3,R);  // Eq. (13) Arun et al. 1987

      // if(opt_v) printf("det(R) = %f\n", det3(R));

      if( det3(R) < 0.0 ) 
      {  
         // if(opt_v) printf("Negative determinant (reflection) detected\n");
         V[2] *= -1.0; 
         V[5] *= -1.0; 
         V[8] *= -1.0; 
         multi(V,3,3,Ut,3,3,R);  // Eq. (13) Arun et al. 1987

         // if(opt_v) printf("Corrected det(R) = %f\n", det3(R));
      }

   multi(R,3,3,Pavg,3,1,T);
   for(int i=0; i<3; i++) T[i] = Qavg[i]-T[i];  // Eq. (10) Arun et al. 1987

   // if(opt_v) printMatrix(T,3,1,"Translation:",NULL);

   set_to_I(TLM,4);

   TLM[0] = R[0];
   TLM[1] = R[1];
   TLM[2] = R[2];
   TLM[3] = T[0];

   TLM[4] = R[3];
   TLM[5] = R[4];
   TLM[6] = R[5];
   TLM[7] = T[1];

   TLM[8] = R[6];
   TLM[9] = R[7];
   TLM[10] = R[8];
   TLM[11] = T[2];

   return;
}

void Procrustes(float *Q, int n, float *P, float *TLM)
{
   float Pavg[3];
   float Qavg[3];
   float Ht[9]; // H=P*Q' Ht=Q*P' (' means transpose in my notation)
   //float Ut[9], V[9], I[9];
   float Ut[9], V[9];
   float T[3]; // 3x1 translation vector
   float R[9]; // 3x3 rotation matrix 
   float S[3]; // 3x1 vector of singular values

   Pavg[0] = (float)removeVectorMean(P, n);
   Pavg[1] = (float)removeVectorMean(P+n, n);
   Pavg[2] = (float)removeVectorMean(P+2*n, n);

   Qavg[0] = (float)removeVectorMean(Q, n);
   Qavg[1] = (float)removeVectorMean(Q+n, n);
   Qavg[2] = (float)removeVectorMean(Q+2*n, n);

   mat_mat_trans(Q,3,n,P,3,Ht);

      svd(Ht, 3, 3, Ut, V, S);

      multi(V,3,3,Ut,3,3,R);  // Eq. (13) Arun et al. 1987

      //if(opt_v) printf("det(R) = %f\n", det3(R));

      if( det3(R) < 0.0 ) 
      {  
         //if(opt_v) printf("Negative determinant (reflection) detected\n");
         V[2] *= -1.0; 
         V[5] *= -1.0; 
         V[8] *= -1.0; 
         multi(V,3,3,Ut,3,3,R);  // Eq. (13) Arun et al. 1987

         //if(opt_v) printf("Corrected det(R) = %f\n", det3(R));
      }

      multi(R,3,3,Pavg,3,1,T);
      for(int i=0; i<3; i++) T[i] = Qavg[i]-T[i];  // Eq. (10) Arun et al. 1987

      //if(opt_v) printMatrix(T,3,1,"Translation:",NULL);

   set_to_I(TLM,4);

   TLM[0] = R[0];
   TLM[1] = R[1];
   TLM[2] = R[2];
   TLM[3] = T[0];

   TLM[4] = R[3];
   TLM[5] = R[4];
   TLM[6] = R[5];
   TLM[7] = T[1];

   TLM[8] = R[6];
   TLM[9] = R[7];
   TLM[10] = R[8];
   TLM[11] = T[2];

   return;
}

// lmfile (landmarks file) - This is a text file with format:
// i1 j1 k1
// i2 j2 k2
// i3 j3 k3
// where 
// (i1, j1, k1) are the coordinates of the AC
// (i2, j2, k2) are the coordinates of the PC, and
// (i3, j3, k3) are the coordinates of the RP.
// subfile (subject file) - 3D T1W volume of type short in NIFTI format
// TPIL - output 4x4 rigid-body transformation matrix that would transform
// subfile into a standardized PIL orientation
void new_PIL_transform(const char *subfile,const char *lmfile,char *orient,float *TPIL, int SAVE_MRX_FLAG)
{
   float Qavg[3]; // average of rows of Q
   float Pavg[3]; // average of rows of P
   float TPIL0[16]; // transforms the original image to MSP/AC-PC aligned PIL orientation
   int n; // number of landmarks
   float *LM;  // (3 x n) matrix of detected landmarks
   float *invT;
   char filename[2048];
   SHORTIM subimPIL; 
   nifti_1_header PILbraincloud_hdr;
   // subfile without the directory structure and extension
   char subfile_prefix[DEFAULT_STRING_LENGTH]; 
   char imagedir[DEFAULT_STRING_LENGTH]; 
   char modelfile[DEFAULT_STRING_LENGTH]="";

   if( niftiFilename(subfile_prefix, subfile)==0 ) exit(0);
   getDirectoryName(subfile, imagedir);

   getARTHOME();

   // initial TPIL0 using old PIL transformation
   char opt_CENTER_AC_old;   
   opt_CENTER_AC_old = opt_CENTER_AC; // record opt_CENTER_AC
   opt_CENTER_AC = NO; // temporarily set opt_CENTER_AC to NO
   standard_PIL_transformation(subfile, lmfile, orient, 0, TPIL0);
   opt_CENTER_AC = opt_CENTER_AC_old; // restore opt_CENTER_AC

   /////////////////////////////////////////////////////////
   // Read input volume from subfile
   SHORTIM subim; 
   nifti_1_header subim_hdr;  

   subim.v = (short *)read_nifti_image(subfile, &subim_hdr);

   if(subim.v==NULL)
   {
      printf("Error reading %s, aborting ...\n", subfile);
      exit(1);
   }

   set_dim(subim, subim_hdr);
   /////////////////////////////////////////////////////////

   /////////////////////////////////////////////////////////
   // Reslice subim to PIL space 
   sprintf(filename,"%s/PILbrain.nii",ARTHOME);
   PILbraincloud_hdr=read_NIFTI_hdr(filename);

   set_dim(subimPIL, PILbraincloud_hdr);
   invT = inv4(TPIL0);
   resliceImage(subim, subimPIL, invT, LIN);
   free(invT);

   // Detect 8 landmarks on the MSP on the subimPIL volume
   snprintf(modelfile,sizeof(modelfile),"%s/orion.mdl",ARTHOME);

   LM=detect_landmarks(subimPIL, modelfile, n);

   convert_to_xyz(LM, n, subimPIL);
 
   // This block outputs the locations of the detected landmarks 
   // in (i,j,k) coordinates of the native space of the input volume
   if(opt_txt)
   {
      FILE *fp;
      float landmark[4];
      invT = inv4(TPIL0);
      sprintf(filename,"%s/%s_orion.txt",imagedir,subfile_prefix);
      fp=fopen(filename,"w");
      if(fp==NULL) file_open_error(filename);
      for(int i=0; i<n; i++)
      {
         landmark[0]=LM[0*n + i];
         landmark[1]=LM[1*n + i];
         landmark[2]=LM[2*n + i];
         landmark[3]=1;
         multi(invT,4,4,landmark,4,1,landmark);
         convert_to_ijk(landmark, 1, subim);
         fprintf(fp,"%7.3f %7.3f %7.3f\n",landmark[0], landmark[1], landmark[2]);
      }
      fclose(fp);
      free(invT);
   }

   float *P;
   float *Q; // Image using Q insted of P' in Eq. (1) of Arun et al. 1987
   float TLM[16];

   Q =  (float *)calloc(3*n,sizeof(float));
   P =  (float *)calloc(3*n,sizeof(float));
   for(int i=0; i<3*n; i++) P[i]=LM[i];

   delete subimPIL.v;
   /////////////////////////////////////////////////////////

   // read Q
   {
      FILE *fp;
      int r, r2;
      int cm[3];

      fp=fopen(modelfile, "r");
      if(fp==NULL) file_open_error(modelfile);

      if( fread(&n, sizeof(int), 1, fp) == 1) {};
      if( fread(&r, sizeof(int), 1, fp) == 1) {};
      if( fread(&r2, sizeof(int), 1, fp) == 1) {};
      SPH refsph(r);

      for(int i=0; i<n; i++)
      {
         if( fread(cm, sizeof(int), 3, fp) == 3) {};
         if( fread(refsph.v, sizeof(float), refsph.n, fp) == (size_t)refsph.n ) {};
         Q[0*n + i]=cm[0];
         Q[1*n + i]=cm[1];
         Q[2*n + i]=cm[2];
      }

      fclose(fp);

      convert_to_xyz(Q, n, subimPIL);
   }

   Procrustes(Q, Qavg, n, P, Pavg, TLM);

   multi(TLM,4,4,TPIL0,4,4,TPIL);

   if(opt_CENTER_AC)
   {
     float ac[4]; 
     float Ttmp[16];

     ac[0] = P[1] + Pavg[0];  // landmark #1 is AC
     ac[1] = P[n+1] + Pavg[1];
     ac[2] = P[2*n+1] + Pavg[2];
     ac[3] = 1.0;

     multi(TLM,4,4,ac,4,1,ac);

     Ttmp[0]=1.0;  Ttmp[1]=0.0;  Ttmp[2]=0.0;  Ttmp[3]=-ac[0];  // makes ac the center of the FOV
     Ttmp[4]=0.0;  Ttmp[5]=1.0;  Ttmp[6]=0.0;  Ttmp[7]=-ac[1];  // ac[1] should be equal to pc[1]
     Ttmp[8]=0.0;  Ttmp[9]=0.0;  Ttmp[10]=1.0; Ttmp[11]=0;
     Ttmp[12]=0.0; Ttmp[13]=0.0; Ttmp[14]=0.0; Ttmp[15]=1.0;

     multi(Ttmp,4,4,TPIL,4,4,TPIL);
   }

   // create the *LM.ppm image and save orion landmarks in xyz coordinates in PIL orientation
   {
      int *lmx, *lmy;
      float lm[4]; 
      FILE *fp;
      float Ttmp[16];

      invT = inv4(TPIL);
      resliceImage(subim, subimPIL, invT, LIN);
      free(invT);

      lmx = (int *)calloc(n,sizeof(int));
      lmy = (int *)calloc(n,sizeof(int));

      if(opt_CENTER_AC)
      {
        float ac[4]; 

        ac[0] = P[1] + Pavg[0];  // landmark #1 is AC
        ac[1] = P[n+1] + Pavg[1];
        ac[2] = P[2*n+1] + Pavg[2];
        ac[3] = 1.0;

        multi(TLM,4,4,ac,4,1,ac);

        Ttmp[0]=1.0;  Ttmp[1]=0.0;  Ttmp[2]=0.0;  Ttmp[3]=-ac[0];  // makes ac the center of the FOV
        Ttmp[4]=0.0;  Ttmp[5]=1.0;  Ttmp[6]=0.0;  Ttmp[7]=-ac[1];  // ac[1] should be equal to pc[1]
        Ttmp[8]=0.0;  Ttmp[9]=0.0;  Ttmp[10]=1.0; Ttmp[11]=0;
        Ttmp[12]=0.0; Ttmp[13]=0.0; Ttmp[14]=0.0; Ttmp[15]=1.0;
      }

      sprintf(filename,"%s/%s_orion_PIL.txt",imagedir,subfile_prefix);
      fp=fopen(filename,"w");
      if(fp==NULL) file_open_error(filename);

      for(int i=0; i<n; i++)
      {
         lm[0] = P[i] + Pavg[0];
         lm[1] = P[n+i] + Pavg[1];
         lm[2] = P[2*n+i] + Pavg[2];
         lm[3] = 1.0;

         multi(TLM,4,4,lm,4,1,lm);

         if(opt_CENTER_AC) multi(Ttmp,4,4,lm,4,1,lm);

         fprintf(fp,"%7.3f %7.3f %7.3f\n",lm[0], lm[1], lm[2]);
         convert_to_ijk(lm, 1, subimPIL);

         lmx[i]=(int)( lm[0] + 0.5 );
         lmy[i]=(int)( lm[1] + 0.5 );
      }
      fclose(fp);

      if(opt_ppm || opt_png)
      {
         sprintf(filename,"%s/%s_orion.ppm",imagedir, subfile_prefix);
         mspPPM(subimPIL, lmx, lmy, n, filename);
      }

      delete lmx;
      delete lmy;
      delete subimPIL.v;
   }

  // save the PIL transformation in <subfile_prefix>_PIL.mrx
  if(SAVE_MRX_FLAG == 1)
  {
    FILE *fp;
    snprintf(filename,sizeof(filename),"%s/%s_PIL.mrx",imagedir, subfile_prefix);
    fp=fopen(filename,"w");
    if(fp==NULL) file_open_error(filename);
    printMatrix(TPIL,4,4,"",fp);
    fclose(fp);
  }

   {
      float ssd1=0.0;
      float ssd3=0.0;
      float x[4], y[4];
      
      x[3]=y[3]=1.0;

      for(int i=0; i<n; i++)
      {
         x[0] = Q[i] + Qavg[0];
         x[1] = Q[i+n] + Qavg[1];
         x[2] = Q[i+2*n] + Qavg[2];

         y[0] = P[i]+Pavg[0];
         y[1] = P[i+n]+Pavg[1];
         y[2] = P[i+2*n]+Pavg[2];

         ssd1 += (x[0]-y[0])*(x[0]-y[0]); 
         ssd1 += (x[1]-y[1])*(x[1]-y[1]); 
         ssd1 += (x[2]-y[2])*(x[2]-y[2]); 

         multi(TLM,4,4,y,4,1,y);
         ssd3 += (x[0]-y[0])*(x[0]-y[0]); 
         ssd3 += (x[1]-y[1])*(x[1]-y[1]); 
         ssd3 += (x[2]-y[2])*(x[2]-y[2]); 
      }
      //if(opt_v) printf("SSD (MSP + AC/PC transformation) = %f\n",ssd1);
      //if(opt_v) printf("SSD (MSP + AC/PC + LM transformation) = %f\n",ssd3);
   }

   delete subim.v;
   free(P);
   free(Q);
}

// If is empty, i.e. orient="", that the orientation is read from imfile header
// otherwise, the orien overrides the header information
void standard_PIL_transformation(const char *imfile, const char *lmfile, char *orient, int verbose, float *TPIL)
{
   char modelfile[1024]="";

   DIM dim;
   float ac[4]={0.0, 0.0, 0.0, 1.0};
   float pc[4]={0.0, 0.0, 0.0, 1.0};
   float rp[4]={0.0, 0.0, 0.0, 1.0};

   //////////////////////////////////////////////////////////////////////////////////
   // Detect the MSP plane 
   //////////////////////////////////////////////////////////////////////////////////
   float Tmsp[16]; // takes the image to PIL orientation -- no AC/PC alignment

   findMSP(imfile,orient,lmfile,Tmsp,verbose,dim);

   ///////////////////////////////////////////////////////////////////////////////////////////////
   // AC/PC detection 
   ///////////////////////////////////////////////////////////////////////////////////////////////

   opt_AC=opt_PC=opt_RP=YES;  // by default find these automatically

   if( lmfile[0]!='\0' )
   {
      FILE *fp;
      
      fp = fopen(lmfile,"r");
      if(fp==NULL) file_open_error(lmfile);
      if( fscanf(fp,"%f %f %f\n", &ac[0], &ac[1], &ac[2])==3 ) {}; ac[3]=1;
      if( fscanf(fp,"%f %f %f\n", &pc[0], &pc[1], &pc[2])==3 ) {}; pc[3]=1;
      if( fscanf(fp,"%f %f %f\n", &rp[0], &rp[1], &rp[2])==3 ) {}; rp[3]=1;
      fclose(fp);

      opt_AC=opt_PC=opt_RP=NO;  // Do not find AC, PC, RP
   }

   opt_MSP=NO;
   detect_AC_PC_MSP(imfile,orient,modelfile,ac,pc,rp,Tmsp,verbose,0);
   //restore to default state, this is important for subsequent calls to detect_AC_PC_MSP()
   opt_MSP=YES; 

   // convert the AC/PC from (i,j,k) in original space to (x,y,z) in PIL space
   orig_ijk_to_pil_xyz(Tmsp, dim, ac, pc);

   // if opt_CENTER_AC=NO, puts FOV center at halfway btw AC and PC to match $ARTHOME/PILbrain.nii
   // if opt_CENTER_AC=YES, makes AC the FOV center
   ACPCtransform(TPIL, Tmsp, ac, pc, opt_CENTER_AC);

   return;
}

// Matrices of type T are called "signed permutation matrices"
// A permutation matrix is one in which there is exactly one 1 in 
// each row and column.  In a signed permutation matrix the non-zero
// entries can be 1 or -1.
void PILtransform(const char *inputOrientCode, float *T)
{
   char c;

   // Initialize T
   for(int i=0; i<15; i++) T[i]=0.0;
   T[15]=1.0;

   for(int j=0; j<3; j++)
   {
      // read the jth element of the inputOrientCode vector
      c = inputOrientCode[j];		

      c=toupper(c);	// convert c to upper case

      // set the jth column of the T
      if( c=='P' )		T[j]   =  1.0;
      else if( c=='A' )	T[j]   = -1.0;
      else if( c=='I' )	T[4+j] =  1.0;
      else if( c=='S' )	T[4+j] = -1.0;
      else if( c=='L' )	T[8+j] =  1.0;
      else if( c=='R' )	T[8+j] = -1.0;
   }

   return;
}

void inversePILtransform(const char *outputOrientCode, float *T)
{
   char c;

   // Initialize T
   for(int i=0; i<16; i++) T[i]=0.0;
   T[15]=1.0;

   for(int i=0; i<3; i++)
   {
      // read the ith element of the outputOrientCode vector
      c = outputOrientCode[i];		
   
      c=toupper(c);	// convert c to upper case

      // set the ith row of the T
      if( c=='P' )		T[4*i] 		=  1.0;
      else if( c=='A' )	T[4*i] 		= -1.0;
      else if( c=='I' )	T[4*i + 1]	=  1.0;
      else if( c=='S' )	T[4*i + 1]	= -1.0;
      else if( c=='L' )	T[4*i + 2]	=  1.0;
      else if( c=='R' )	T[4*i + 2]	= -1.0;
   }

   return;
}
