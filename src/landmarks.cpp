#include <babak_lib.h>
#include <sph.h>
#include <stats.h>
#include <nifti1.h>
#include <stdio.h>
#include <stdlib.h>

extern void read_nifti_image(const char *filename, unsigned char **im, nifti_1_header *hdr);
extern void read_nifti_image(const char *filename, short **im, nifti_1_header *hdr);
extern void standardize(float *x, int n);

//im: the input image 
//nx,ny,nz: size of the image
//lm[3]: coordinate of the landmark (satrting from 0)
//ppmfile: name of the ppmfile as the output
void makePPM(SHORTIM im, int *lm, const char *ppmfile)
{
   unsigned char yellow[3]={0xFF,0xFF,0x00};
   unsigned char *imgTemp;
   FILE *fp;

   int nv=im.nx*im.ny*im.nz;
   int np=im.nx*im.ny;

   int ii=lm[0];
   int jj=lm[1];
   int kk=lm[2];

   float temp;

   imgTemp=(unsigned char *)calloc(nv,sizeof(unsigned char));

   int low, high;
   setLowHigh(im.v, nv, &low, &high);

   for(int i=0;i<nv;i++)
   {
      if( im.v[i] >= high ) temp=255.0;
      else temp=im.v[i]*255.0/high;

      imgTemp[i]=(unsigned char)temp;
   }

   fp=fopen(ppmfile,"w");
   if(fp==NULL) file_open_error(ppmfile);

   /*Write the header part of the PPM file*/
   fprintf(fp,"P6\n");
   fprintf(fp,"# x=%d, y=%d, z=%d\n", ii, jj, kk);
   fprintf(fp,"%d %d\n", im.nx+im.nz+im.nz, im.ny);
   fprintf(fp,"255\n");

   /*Write the data part of the PPM file*/
   for (int j=0,ik=0;j<im.ny & ik<im.nx;j++,ik++)
   {
      for(int i=0;i<im.nx;i++)
      {
         if(i==ii | j==jj) fwrite(yellow,1,3,fp);
         else {
            fwrite(imgTemp+np*kk+im.nx*j+i,1,1,fp);
            fwrite(imgTemp+np*kk+im.nx*j+i,1,1,fp);
            fwrite(imgTemp+np*kk+im.nx*j+i,1,1,fp);
         }
      }

      for(int k=0;k<im.nz;k++)
      {
         if(j==jj | k==kk) fwrite(yellow,1,3,fp);
         else {
            fwrite(imgTemp+np*k+im.nx*j+ii,1,1,fp);
            fwrite(imgTemp+np*k+im.nx*j+ii,1,1,fp);
            fwrite(imgTemp+np*k+im.nx*j+ii,1,1,fp);
         }
      }

      for(int k=0;k<im.nz;k++)
      {
         if(ik==ii | k==kk) fwrite(yellow,1,3,fp);
         else {
            fwrite(imgTemp+np*k+im.nx*jj+ik,1,1,fp);
            fwrite(imgTemp+np*k+im.nx*jj+ik,1,1,fp);
            fwrite(imgTemp+np*k+im.nx*jj+ik,1,1,fp);
         }
      }
   }

   fclose(fp);

   delete imgTemp;
}

float detect_lm(SPH &searchsph, SPH &testsph, SHORTIM testim, int lmcm[], SPH &refsph, int lm[])
{
   //float ccmax=-1.0;
   float ccmax=0.0;

   lm[0]=lmcm[0];
   lm[1]=lmcm[1];
   lm[2]=lmcm[2];

   for(int n=0; n<searchsph.n; n++)
   {
      testsph.set(testim, lmcm[0]+searchsph.i[n], lmcm[1]+searchsph.j[n], lmcm[2]+searchsph.k[n]);
      //standardize(testsph.v, testsph.n);
      standardize(testsph.v, refsph.v, testsph.n);  //refsph is used as a mask
      searchsph.v[n] = dot(testsph.v, refsph.v, testsph.n);
      if( searchsph.v[n] > ccmax )
      {
         ccmax = searchsph.v[n];
         lm[0] = lmcm[0]+searchsph.i[n];
         lm[1] = lmcm[1]+searchsph.j[n];
         lm[2] = lmcm[2]+searchsph.k[n];
      }
   }

   return(ccmax);
}

//nl: number of landmarks (returned)
float *detect_landmarks(const char *subfile, const char *mdlfile, int &nl, char ppmflg)
{
   char prefix[1024];
   char filename[1024];
   SHORTIM subim; 
   nifti_1_header hdr;
   FILE *fpi;
   int r; // patch radius
   int R; // search region radius
   int lm[3], lmcm[3];
   float ccmax;
   float *P;

   if( niftiFilename(prefix, subfile)==0 )
   {
      exit(0);
   }

   subim.v = (short *)read_nifti_image(subfile, &hdr);
   subim.nx = hdr.dim[1];
   subim.ny = hdr.dim[2];
   subim.nz = hdr.dim[3];
   subim.nv = subim.nx*subim.ny*subim.nz;
   subim.np = subim.nx*subim.ny;

   fpi=fopen(mdlfile, "r");
   if(fpi==NULL) file_open_error(mdlfile);

   fread(&nl, sizeof(int), 1, fpi);
   fread(&r, sizeof(int), 1, fpi);
   fread(&R, sizeof(int), 1, fpi);
   SPH refsph(r);
   SPH subsph(r);
   SPH searchsph(R);

   P = (float *)calloc(3*nl, sizeof(float));

   for(int l=0; l<nl; l++)
   {
      fread(&lmcm[0], sizeof(int), 1, fpi);
      fread(&lmcm[1], sizeof(int), 1, fpi);
      fread(&lmcm[2], sizeof(int), 1, fpi);

      fread(refsph.v, sizeof(float), refsph.n, fpi);

      ccmax = detect_lm(searchsph, subsph, subim, lmcm, refsph, lm);

      P[0*nl + l]=lm[0];
      P[1*nl + l]=lm[1];
      P[2*nl + l]=lm[2];

      if(ppmflg)
      {
         sprintf(filename,"%s_lm_%03d.ppm",prefix,l);
         makePPM(subim, lm, filename);
      }
   }

   fclose(fpi);
 
   return(P);
}

float *detect_landmarks(SHORTIM subim, const char *mdlfile, int &nl)
{
   FILE *fpi;
   int r; // patch radius
   int R; // search region radius
   int lm[3], lmcm[3];
   float ccmax;
   float *P;

   fpi=fopen(mdlfile, "r");

   if(fpi==NULL)
   {
      printf("Could not find %s, aborting ...\n",mdlfile);
      exit(0);
   }

   fread(&nl, sizeof(int), 1, fpi);
   fread(&r, sizeof(int), 1, fpi);
   fread(&R, sizeof(int), 1, fpi);
   SPH refsph(r);
   SPH subsph(r);
   SPH searchsph(R);

   P = (float *)calloc(3*nl, sizeof(float));

   for(int l=0; l<nl; l++)
   {
      fread(&lmcm[0], sizeof(int), 1, fpi);
      fread(&lmcm[1], sizeof(int), 1, fpi);
      fread(&lmcm[2], sizeof(int), 1, fpi);
      fread(refsph.v, sizeof(float), refsph.n, fpi);

      ccmax = detect_lm(searchsph, subsph, subim, lmcm, refsph, lm);

      P[0*nl + l]=lm[0];
      P[1*nl + l]=lm[1];
      P[2*nl + l]=lm[2];
   }

   fclose(fpi);
 
   return(P);
}

float *read_landmark_centers(const char *mdlfile, int &nl)
{
   FILE *fpi;
   int r; // patch radius
   int R; // search region radius
   int lmcm[3];
   float *Q;

   fpi=fopen(mdlfile, "r");
   if(fpi==NULL) file_open_error(mdlfile);

   fread(&nl, sizeof(int), 1, fpi);
   fread(&r, sizeof(int), 1, fpi);
   fread(&R, sizeof(int), 1, fpi);
   SPH refsph(r);

   Q = (float *)calloc(3*nl, sizeof(float));

   for(int l=0; l<nl; l++)
   {
      fread(&lmcm[0], sizeof(int), 1, fpi);
      fread(&lmcm[1], sizeof(int), 1, fpi);
      fread(&lmcm[2], sizeof(int), 1, fpi);
      fread(refsph.v, sizeof(float), refsph.n, fpi);

      Q[0*nl + l]=lmcm[0];
      Q[1*nl + l]=lmcm[1];
      Q[2*nl + l]=lmcm[2];
   }

   fclose(fpi);
 
   return(Q);
}
