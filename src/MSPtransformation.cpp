#include <cfloat>
#include <math.h>
#include "babak_lib.h"
#include "bbk_linear_algebra.h"
#include "bbk.h"

bool MSPtransformation(
   const char *filename,
   const char *orient,
   const char *lmfile,
   float *Tmsp,
   DIM &dim
)
{
   // Validate arguments.
   if(filename == nullptr ||
      orient == nullptr ||
      lmfile == nullptr ||
      Tmsp == nullptr)
   {
      return false;
   }

   // Ensure that a nonempty orientation code contains exactly three characters.
   if(orient[0] != '\0' &&
      (orient[1] == '\0' || orient[2] == '\0' || orient[3] != '\0'))
   {
      return false;
   }

   char imorient[4];

   // If an orientation code is not provided as input, then
   // determine the orientation from the input image file.
   if(orient[0] == '\0')
   {
      if(!getNiftiImageOrientation(filename, imorient))
      {
         return false;
      }
   }
   else
   {
      imorient[0] = orient[0];
      imorient[1] = orient[1];
      imorient[2] = orient[2];
      imorient[3] = '\0';
   }

   // Ensure that the orientation code is a valid code.
   if(!valid_orientation_code(imorient))
   {
      return false;
   }

   // At this point, Tmsp only makes the orientation PIL without MSP alignment
   PILtransform(imorient, Tmsp);

   nifti_1_header hdr;

   if ( !read_nifti_hdr(filename, &hdr) )
   {
      return false;
   }

   set_dim(dim, hdr);

   float tmpT[16];
   float nrml[3];
   float d;
   FILE *fp;

   float ac[4] = {0.0f, 0.0f, 0.0f, 1.0f};
   float pc[4] = {0.0f, 0.0f, 0.0f, 1.0f};
   float rp[4] = {0.0f, 0.0f, 0.0f, 1.0f};

   float u[3]; // ac - rp
   float v[3]; // pc - rp

   // Read landmarks.
   // Add code to detect errors here.
   fp = fopen(lmfile, "r");

   if(fp == nullptr)
   {
      return false;
   }

   if(fscanf(fp, "%f %f %f\n", &ac[0], &ac[1], &ac[2]) != 3)
   {
      fclose(fp);
      return false;
   }

   if(fscanf(fp, "%f %f %f\n", &pc[0], &pc[1], &pc[2]) != 3)
   {
      fclose(fp);
      return false;
   }

   if(fscanf(fp, "%f %f %f\n", &rp[0], &rp[1], &rp[2]) != 3)
   {
      fclose(fp);
      return false;
   }

   fclose(fp);

   ijk2xyz(
      tmpT,
      hdr.dim[1],
      hdr.dim[2],
      hdr.dim[3],
      hdr.pixdim[1],
      hdr.pixdim[2],
      hdr.pixdim[3]
   );

   // Convert ac, pc, and rp from ijk to xyz coordinates
   // (still in the original space).
   multi(tmpT, 4, 4, ac, 4, 1, ac);
   multi(tmpT, 4, 4, pc, 4, 1, pc);
   multi(tmpT, 4, 4, rp, 4, 1, rp);

   // Convert ac, pc, and rp from xyz coordinates in the original space
   // to PIL space.
   multi(Tmsp, 4, 4, ac, 4, 1, ac);
   multi(Tmsp, 4, 4, pc, 4, 1, pc);
   multi(Tmsp, 4, 4, rp, 4, 1, rp);

   u[0] = ac[0] - rp[0];
   u[1] = ac[1] - rp[1];
   u[2] = ac[2] - rp[2];

   v[0] = pc[0] - rp[0];
   v[1] = pc[1] - rp[1];
   v[2] = pc[2] - rp[2];

   crossProduct(u, v, nrml);

   normalizeVector(nrml, 3);

   d = nrml[0] * ac[0] +
       nrml[1] * ac[1] +
       nrml[2] * ac[2];

   // IMPORTANT
   if(d < 0.0f)
   {
      nrml[0] *= -1.0f;
      nrml[1] *= -1.0f;
      nrml[2] *= -1.0f;

      d *= -1.0f;
   }

   // Ensures that the MSP passes through the center of the FOV in PIL space.
   tmpT[0]  = 1.0f;
   tmpT[1]  = 0.0f;
   tmpT[2]  = 0.0f;
   tmpT[3]  = -d * nrml[0];

   tmpT[4]  = 0.0f;
   tmpT[5]  = 1.0f;
   tmpT[6]  = 0.0f;
   tmpT[7]  = -d * nrml[1];

   tmpT[8]  = 0.0f;
   tmpT[9]  = 0.0f;
   tmpT[10] = 1.0f;
   tmpT[11] = -d * nrml[2];

   tmpT[12] = 0.0f;
   tmpT[13] = 0.0f;
   tmpT[14] = 0.0f;
   tmpT[15] = 1.0f;

   multi(tmpT, 4, 4, Tmsp, 4, 4, Tmsp);

   // Ensures that the MSP is parallel to the x-y plane.
   if(nrml[2] > 0.0f)
   {
      float c;
      float theta;

      c = nrml[2];

      // Clamp c to prevent acos from getting into trouble.
      if(c > 1.0f)
	 c = 1.0f;

      if(c < 0.0f)
	 c = 0.0f;

      theta = (float)acos((double)c);

      rotationMatrix(tmpT, theta, nrml[1], -nrml[0], 0.0f);

      multi(tmpT, 4, 4, Tmsp, 4, 4, Tmsp);
   }
   else
   {
      float c;
      float theta;

      c = -nrml[2];

      // Clamp c to prevent acos from getting into trouble.
      if(c > 1.0f)
	 c = 1.0f;

      if(c < 0.0f)
	 c = 0.0f;

      theta = (float)acos((double)c);

      rotationMatrix(tmpT, theta, -nrml[1], nrml[0], 0.0f);

      multi(tmpT, 4, 4, Tmsp, 4, 4, Tmsp);
   }

   return true;
}

// Inputs:
// filename: NIFTI-1 image filename on which the MSP is to be detected
// orient: if not "" overrides the orientation code in filename
// Outputs:
// Tmsp: 4x4 matrix, rigid-body transformation takes the original image 
// to PIL space without AC-PC alignment, that is, only the MSP of the original 
// image is transformed to the x-y plane with x pointing posteriorly, y pointing 
// inferiorly, and z pointing to the left.  However, the x axis is not necessarily 
// aligned with the AC-PC line.
bool MSPtransformation(const char *filename, const char *orient, float *Tmsp, DIM &dim)
{
   // Validate arguments.
   if(filename == nullptr ||
      orient == nullptr ||
      Tmsp == nullptr)
   {
      return false;
   }

   // Ensure that a nonempty orientation code contains exactly three characters.
   if(orient[0] != '\0' &&
      (orient[1] == '\0' || orient[2] == '\0' || orient[3] != '\0'))
   {
      return false;
   }

   char imorient[4];

   // If an orientation code is not provided as input, then
   // determine the orientation from the input image file.
   if(orient[0] == '\0')
   {
      if(!getNiftiImageOrientation(filename, imorient))
      {
         return false;
      }
   }
   else
   {
      imorient[0] = orient[0];
      imorient[1] = orient[1];
      imorient[2] = orient[2];
      imorient[3] = '\0';
   }

   // Ensure that the orientation code is a valid code.
   if(!valid_orientation_code(imorient))
   {
      return false;
   }

   // At this point, Tmsp only makes the orientation PIL without MSP alignment
   PILtransform(imorient, Tmsp);

   short *im;
   nifti_1_header hdr; // image NIFTI header
   im = (short *)read_nifti_image(filename, &hdr);

   if(im == nullptr)
   {
      return false;
   }

   set_dim(dim, hdr);

   short *imPIL; // original volume after reorientation to PIL
   float dxPIL, dyPIL, dzPIL; // voxel dimensions in PIL orientation
   int nxPIL, nyPIL, nzPIL; // matrix dimensions in PIL orientation

   // Reorient the original volume to PIL orientation.
   imPIL = reorientVolume(
      im,
      hdr.dim[1],
      hdr.dim[2],
      hdr.dim[3],
      hdr.pixdim[1],
      hdr.pixdim[2],
      hdr.pixdim[3],
      Tmsp,
      nxPIL,
      nyPIL,
      nzPIL,
      dxPIL,
      dyPIL,
      dzPIL
   );

   free(im);

   if(imPIL == nullptr)
   {
      return false;
   }

   float A, B, C; // parameters in Ax + By + Cz = 1 in PIL space

   msp(
      imPIL,
      nxPIL,
      nyPIL,
      nzPIL,
      dxPIL,
      dyPIL,
      dzPIL,
      &A,
      &B,
      &C
   );

   free(imPIL);

   float tmpT[16]; // temporary matrix container
   float nrml[3]; // unit vector normal to the MSP
   float L; // sqrtf(A * A + B * B + C * C)
   float d; // distance from the origin to the MSP

   // Normalize (A,B,C)
   L = sqrtf(A * A + B * B + C * C);

   if(L <= FLT_EPSILON)
   {
      return false;
   }

   nrml[0] = A / L;
   nrml[1] = B / L;
   nrml[2] = C / L;

   d = 1.0f / L;

   // Translate the MSP to pass through the origin in PIL space.
   tmpT[0] = 1.0f;
   tmpT[1] = 0.0f;
   tmpT[2] = 0.0f;
   tmpT[3] = -d * nrml[0];

   tmpT[4] = 0.0f;
   tmpT[5] = 1.0f;
   tmpT[6] = 0.0f;
   tmpT[7] = -d * nrml[1];

   tmpT[8] = 0.0f;
   tmpT[9] = 0.0f;
   tmpT[10] = 1.0f;
   tmpT[11] = -d * nrml[2];

   tmpT[12] = 0.0f;
   tmpT[13] = 0.0f;
   tmpT[14] = 0.0f;
   tmpT[15] = 1.0f;

   // Update Tmsp.
   multi(tmpT, 4, 4, Tmsp, 4, 4, Tmsp);

   // Ensures that the MSP is parallel to the x-y plane.
   if(nrml[2] > 0.0f)
   {
      float c, theta;

      c = nrml[2];

      // Clamp c to prevent acos from getting into trouble.
      if(c > 1.0f)
         c = 1.0f;

      if(c < 0.0f)
         c = 0.0f; 

      theta = (float)acos((double)c);

      // ensures that MSP is parallel to the x-y plane
      rotationMatrix(tmpT, theta, nrml[1], -nrml[0], 0.0f);

      // Update Tmsp. 
      multi(tmpT, 4, 4, Tmsp, 4, 4, Tmsp);
   }
   else
   {
      float c, theta;

      c = -nrml[2];

      // Clamp c to prevent acos from getting into trouble.
      if(c > 1.0f)
         c = 1.0f;

      if(c < 0.0f)
         c = 0.0f;

      theta = (float)acos((double)c);

      // ensures that MSP is parallel to the x-y plane
      rotationMatrix(tmpT, theta, -nrml[1], nrml[0], 0.0f);

      // Update Tmsp. 
      multi(tmpT, 4, 4, Tmsp, 4, 4, Tmsp);
   }

   return true;
}
