#include <babak_lib.h>
#include <math.h>

// Inputs:
// filename: NIFTI image filename on which the MSP is to be detected
// orient: if not "" overrides the orientation code in filename
// verbose: flags verbose mode
// Outputs:
// Tmsp: 4x4 matrix, rigid-body transformation takes the original image to PIL space without AC-PC alignment, that is, only the
// MSP of the original image is transformed to the x-y plane with x pointing posteriorly, y pointing inferiorly, and z pointing
// to the left.  However, the x axis is not necessarily aligned with the AC-PC line.
void findMSP(const char *filename, char *orient, const char *lmfile, float *Tmsp, int verbose, DIM &dim)
{
   float A, B, C; // parameters in Ax+By+Cz=1 in the original image space
   float tmpT[16]; // temporary matrix container 
   float nrml[3]; // unit vector normal to the MSP
   short *im;
   float L; // sqrtf(A*A + B*B + C*C)
   float d; // shortest distance from original to the plane
   nifti_1_header hdr;  // image NIFTI header

   ///////////////////////////////////////////////////////////////////////////////////////////////
   // determine the orientation code 
   ///////////////////////////////////////////////////////////////////////////////////////////////
   
   // if orient is not provided at the input, determine the orientation from the input image file
   if(orient[0]=='\0')
   {
      getNiftiImageOrientation(filename, orient);
   }

   if(orient[0]=='\0')
   {
      errorMessage("Image orientation could not be determined.");
   }

   if ( isOrientationCodeValid(orient) == 0)
   {
      printf("\nImage orientation: %s\n",orient);
      errorMessage("Image orientation code is not one of the 48 legal ones.");
   }

   if(verbose)
   {
      printf("Image orientation: %s\n",orient);
   }

   ///////////////////////////////////////////////////////////////////////////////////////////////
   // read image 
   ///////////////////////////////////////////////////////////////////////////////////////////////
   im = (short *)read_nifti_image(filename, &hdr);

   if(im==NULL)
   {
      printf("Error reading %s, aborting ...\n", filename);
      exit(1);
   }

   if(verbose)
   {
      printf("Image matrix size = %d x %d x %d (voxels)\n", hdr.dim[1], hdr.dim[2], hdr.dim[3]);
      printf("Image voxel size = %8.6f x %8.6f x %8.6f (mm^3)\n", hdr.pixdim[1], hdr.pixdim[2], hdr.pixdim[3]);
   }

   set_dim(dim, hdr);

   ///////////////////////////////////////////////////////////////////////////////////////////////
   // find the MSP for the baseline image 
   ///////////////////////////////////////////////////////////////////////////////////////////////

   // At this point, Tmsp only makes the orientation PIL without MSP alignment
   PILtransform(orient, Tmsp);

   if( lmfile[0]=='\0')
   {
      float cc; // a variable to store correlation coefficient values
      short *imPIL; // original volume after reorientation to PIL
      float dxPIL, dyPIL, dzPIL; // voxel dimensions in PIL orientation
      int nxPIL, nyPIL, nzPIL; // matrix dimensions in PIL orientation

      // reorient the original volume to PIL orientation
      imPIL=reorientVolume(im,hdr.dim[1],hdr.dim[2],hdr.dim[3],hdr.pixdim[1],hdr.pixdim[2],hdr.pixdim[3],
      Tmsp,nxPIL,nyPIL,nzPIL,dxPIL,dyPIL,dzPIL);

      cc=msp(imPIL,nxPIL,nyPIL,nzPIL,dxPIL,dyPIL,dzPIL,&A,&B,&C);
      delete imPIL;

      // normalize (A,B,C)
      L=sqrtf(A*A+B*B+C*C);
      nrml[0]= A/L;
      nrml[1]= B/L;
      nrml[2]= C/L;

      d = 1/L;

      if(verbose) printf("Interhemispheric correlation = %6.4f\n",cc);
   }
   else // read from input file
   {
      FILE *fp;
      float ac[4]={0.0, 0.0, 0.0, 1.0};
      float pc[4]={0.0, 0.0, 0.0, 1.0};
      float rp[4]={0.0, 0.0, 0.0, 1.0};
      float u[3]; // ac-rp
      float v[3]; // pc-rp

      // read landmarks, add code to detect errors here
      fp = fopen(lmfile,"r");
      if(fp==NULL) file_open_error(lmfile);
      if( fscanf(fp,"%f %f %f\n", &ac[0], &ac[1], &ac[2]) != 3)
      {
         printf("Failed to read \"ac\" variable.\n");
      }
      if( fscanf(fp,"%f %f %f\n", &pc[0], &pc[1], &pc[2]) != 3)
      {
         printf("Failed to read \"pc\" variable.\n");
      }
      if( fscanf(fp,"%f %f %f\n", &rp[0], &rp[1], &rp[2]) != 3)
      {
         printf("Failed to read \"rp\" variable.\n");
      }
      fclose(fp);

      // convert ac, pc, and rp from ijk to xyz coordinates (still in the original space)
      ijk2xyz(tmpT, hdr.dim[1], hdr.dim[2], hdr.dim[3], hdr.pixdim[1], hdr.pixdim[2], hdr.pixdim[3]);
      multi(tmpT, 4, 4,  ac, 4,  1, ac);
      multi(tmpT, 4, 4,  pc, 4,  1, pc);
      multi(tmpT, 4, 4,  rp, 4,  1, rp);

      // convert ac, pc, and rp from xyz coordinates in the original space to PIL space
      multi(Tmsp, 4, 4,  ac, 4,  1, ac);
      multi(Tmsp, 4, 4,  pc, 4,  1, pc);
      multi(Tmsp, 4, 4,  rp, 4,  1, rp);

      u[0]=ac[0]-rp[0];
      u[1]=ac[1]-rp[1];
      u[2]=ac[2]-rp[2];

      v[0]=pc[0]-rp[0];
      v[1]=pc[1]-rp[1];
      v[2]=pc[2]-rp[2];

      // (A,B,C) is the cross-product uXv
      nrml[0]=u[1]*v[2]-u[2]*v[1];
      nrml[1]=u[2]*v[0]-u[0]*v[2];
      nrml[2]=u[0]*v[1]-u[1]*v[0];

      normalizeVector(nrml, 3);

      d = nrml[0]*ac[0]+nrml[1]*ac[1]+nrml[2]*ac[2];

      // IMPORTANT
      if(d<0.0)
      {
         nrml[0] *= -1.0; 
         nrml[1] *= -1.0; 
         nrml[2] *= -1.0;
         d *= -1.0;
      }

      // makes the equation Ax+By+Cz=1
      A = nrml[0]/d;
      B = nrml[1]/d;
      C = nrml[2]/d;
   }

   // ensures that the MSP passes through the center of the FOV in PIL space
   tmpT[0]=1.0;  tmpT[1]=0.0;  tmpT[2]=0.0;  tmpT[3]=-d*nrml[0];
   tmpT[4]=0.0;  tmpT[5]=1.0;  tmpT[6]=0.0;  tmpT[7]=-d*nrml[1];
   tmpT[8]=0.0;  tmpT[9]=0.0;  tmpT[10]=1.0; tmpT[11]=-d*nrml[2];
   tmpT[12]=0.0; tmpT[13]=0.0; tmpT[14]=0.0; tmpT[15]=1.0; 
   multi(tmpT, 4, 4,  Tmsp, 4,  4, Tmsp);

   // ensures that MSP is parallel to the x-y plane
   if( nrml[2] > 0 )
   {
      float c, theta;
      c = nrml[2];
      if( c>1.0 ) c=1.0; // to prevent acos from getting into trouble
      if( c<0.0 ) c=0.0; // to prevent acos from getting into trouble
      theta = (float)acos((double)c);

      // ensures that MSP is parallel to the x-y plane
      rotate(tmpT, theta, nrml[1], -nrml[0], 0.0);
      multi(tmpT, 4, 4,  Tmsp, 4,  4, Tmsp);
   }
   else
   {
      float c, theta;
      c = -nrml[2];
      if( c>1.0 ) c=1.0; // to prevent acos from getting into trouble
      if( c<0.0 ) c=0.0; // to prevent acos from getting into trouble
      theta = (float)acos((double)c);

      // ensures that MSP is parallel to the x-y plane
      rotate(tmpT, theta, -nrml[1], nrml[0], 0.0);
      multi(tmpT, 4, 4,  Tmsp, 4,  4, Tmsp);
   }

   // switch (A,B,C) from PIL back to the original orientation
   {
      float dum[4]; 
      float iTPIL[16]; // Transformation from PIL to original orientation

      inversePILtransform(orient, iTPIL);
      dum[0]=A; dum[1]=B; dum[2]=C; dum[3]=1;
      multi(iTPIL,4,4,dum,4,1,dum);
      A=dum[0]; B=dum[1]; C=dum[2];
   }

   if(verbose)
   {
      printf("Estimated mid-sagittal plane: (%7.3fx) + (%7.3fy) + (%7.3fz) = 1\n", A,B,C);
   }

//{
//FILE *fp;
//fp = fopen("tt.mrx","w");
//printMatrix(Tmsp, 4, 4, "ART acpcdetect tilt correction matrix:", fp);
//fclose(fp);
//}

   ///////////////////////////////////////////////////////////////////////////////////////////////
   delete im;
}
