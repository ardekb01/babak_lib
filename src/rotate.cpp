#include <math.h>
#include <stdlib.h>
#include "babak_lib.h"
#include "babak_linear_algebra.h"

// Produces a 4x4 transformation matrix for rotating a point by an angle
// alpha about the (x,y,z) axis.
bool rotate(float *R, float alpha, float x, float y, float z)
{
   float CosAlpha, SinAlpha;
   float dum;
   float ax, ay, az;

   if (!set_to_I(R, 4))
      return false;

   if (x == 0.0 && y == 0.0 && z == 0.0)
   {
      return false;
   }

   CosAlpha = (float)cos((float8)alpha);
   SinAlpha = (float)sin((float8)alpha);

   /* Find the unit vector (ax, ay, az) in the direction of (x, y, z). */
   dum = (float)sqrt((float8)x * x + y * y + z * z);

   ax = x / dum;
   ay = y / dum;
   az = z / dum;

   R[0] = ax * ax + CosAlpha - CosAlpha * ax * ax;
   R[1] = ax * ay - CosAlpha * ax * ay - SinAlpha * az;
   R[2] = ax * az - CosAlpha * ax * az + SinAlpha * ay;

   R[4] = ay * ax - CosAlpha * ay * ax + SinAlpha * az;
   R[5] = ay * ay + CosAlpha - CosAlpha * ay * ay;
   R[6] = ay * az - CosAlpha * ay * az - SinAlpha * ax;

   R[8] = az * ax - CosAlpha * az * ax - SinAlpha * ay;
   R[9] = az * ay - CosAlpha * az * ay + SinAlpha * ax;
   R[10] = az * az + CosAlpha - CosAlpha * az * az;

   return true;
}

void SE3_to_se3(float *M, float *w, float *v, float &theta)
{
   float A[9];
   float R[9];
   float wh[9];
   float wwT[9];
   float *invA;
   float d[3];

   R[0] = M[0]; 
   R[1] = M[1]; 
   R[2] = M[2];
   R[3] = M[4]; 
   R[4] = M[5]; 
   R[5] = M[6];
   R[6] = M[8]; 
   R[7] = M[9]; 
   R[8] = M[10];

   // determines w and theta from R
   irodrigues_formula(R, w, theta);

   d[0] = M[3];
   d[1] = M[7];
   d[2] = M[11];

   // temporarily set A = I-R
   A[0]=1.0-R[0]; A[1]=-R[1];    A[2]=-R[2];
   A[3]=-R[3];    A[4]=1.0-R[4]; A[5]=-R[5];
   A[6]=-R[6];    A[7]=-R[7];    A[8]=1.0-R[8];

   wh[0]=0.0; wh[1]=-w[2]; wh[2]=w[1];
   wh[3]=w[2]; wh[4]=0.0; wh[5]=-w[0];
   wh[6]=-w[1]; wh[7]=w[0]; wh[8]=0.0;

   // set A = (I-R)*wh
   multi(A,3,3,wh,3,3,A);

   //compute w*wT
   wwT[0]=w[0]*w[0];
   wwT[1]=w[0]*w[1];
   wwT[2]=w[0]*w[2];
   wwT[3]=w[1]*w[0];
   wwT[4]=w[1]*w[1];
   wwT[5]=w[1]*w[2];
   wwT[6]=w[2]*w[0];
   wwT[7]=w[2]*w[1];
   wwT[8]=w[2]*w[2];

   for(int4 i=0; i<9; i++) A[i] = A[i] + theta*wwT[i];

   invA = inv3(A);

   // compute v

   multi(invA,3,3, d,3,1, v);

   free(invA);

   return;
}

// M: output 4x4 matrix where M is a member of SE(3)
// w: input 3x1 unit vector  (rotation is about this direction)
// v: input 3x1 vector (translation parallel to w will be w^T*v*theta)
// theta: rotation angle in radians
void se3_to_SE3(float *M, float *w, float *v, float theta)
{
   float R[9];
   float wxv[3];
   float wTv;
   float Rwxv[3];

   // set the last row of M
   M[12] = 0.0;
   M[13] = 0.0;
   M[14] = 0.0;
   M[15] = 1.0;

   // also ensure w is a unit vector
   rodrigues_formula(R, w, theta);  // determines R given w and theta

   // Set the upper-right 3x3 submatrix of M equal to R
   M[0]  = R[0];
   M[1]  = R[1];
   M[2]  = R[2];
   M[4]  = R[3];
   M[5]  = R[4];
   M[6]  = R[5];
   M[8]  = R[6];
   M[9]  = R[7];
   M[10] = R[8];

   // the remaining code is for computing M[3], M[7], and M[11]
   
   // compute the cross-product between w and v
   wxv[0] = w[1]*v[2] - w[2]*v[1];
   wxv[1] = w[2]*v[0] - w[0]*v[2];
   wxv[2] = w[0]*v[1] - w[1]*v[0];

   // compute the the dot-producted between w and v
   wTv = w[0]*v[0] + w[1]*v[1] + w[2]*v[2];

   // compute R*wxV
   multi(R,3,3,wxv,3,1,Rwxv);

   // Finally put everything together
   M[3]  =  wxv[0] - Rwxv[0] + w[0]*wTv*theta;
   M[7]  =  wxv[1] - Rwxv[1] + w[1]*wTv*theta;
   M[11] =  wxv[2] - Rwxv[2] + w[2]*wTv*theta;

   return;
}

void irodrigues_formula(float *R, float *w, float &theta)
{
   float tr; // trace of R
   float ctheta; // cos(theta)
   float sthetax2;  // 2 x sin(theta)

   tr = R[0] + R[4] + R[8];
   
   ctheta = (tr-1.0)/2.0;
   if(ctheta >  1.0)  ctheta=1.0;
   if(ctheta < -1.0) ctheta=-1.0;

   theta=acosf(ctheta);

   // for theta closte to 0.0 w is undefined
   if( theta<0.001) 
   {
      w[0]=1.0; w[2]=0.0; w[3]=0.0; 
      return;
   }

   sthetax2 = 2.0*sinf( theta );

   w[0] = (R[7]-R[5])/sthetax2;
   w[1] = (R[2]-R[6])/sthetax2;
   w[2] = (R[3]-R[1])/sthetax2;

   return;
}

// R: 3x3 output rotation matrix
// w: 3x1 unit vector
// theta: angle in radians
void rodrigues_formula(float *R, float *w, float theta)
{
   float ctheta, stheta, vtheta;
   float wn; // L2 norm of w

   // ensure w is a unit vector
   wn = sqrtf(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

   if( wn == 0.0 )
   {
      R[0] = 1.0;
      R[1] = 0.0;
      R[2] = 0.0;

      R[3] = 0.0;
      R[4] = 1.0;
      R[5] = 0.0;

      R[6] = 0.0;
      R[7] = 0.0;
      R[8] = 0.1;

      return;
   }

   w[0]/=wn;
   w[1]/=wn;
   w[2]/=wn;

   ctheta = cosf( theta );
   stheta = sinf( theta );
   vtheta = 1.0 - ctheta;

   R[0] = w[0]*w[0]*vtheta + ctheta; 
   R[1] = w[0]*w[1]*vtheta - w[2]*stheta; 
   R[2] = w[0]*w[2]*vtheta + w[1]*stheta; 

   R[3] = w[0]*w[1]*vtheta + w[2]*stheta; 
   R[4] = w[1]*w[1]*vtheta + ctheta; 
   R[5] = w[1]*w[2]*vtheta - w[0]*stheta; 

   R[6] = w[0]*w[2]*vtheta - w[1]*stheta; 
   R[7] = w[1]*w[2]*vtheta + w[0]*stheta; 
   R[8] = w[2]*w[2]*vtheta + ctheta; 

   return;
}

// R: 4x4 output homogenious rotation matrix
// w: 3x1 unit vector
// theta: angle in radians
void rodrigues_formula4x4(float *R, float *w, float theta)
{
   float ctheta, stheta, vtheta;
   float wn; // L2 norm of w

   // ensure w is a unit vector
   wn = sqrtf(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

   if( wn == 0.0 )
   {
      R[0] = 1.0;
      R[1] = 0.0;
      R[2] = 0.0;
      R[3] = 0.0;

      R[4] = 0.0;
      R[5] = 1.0;
      R[6] = 0.0;
      R[7] = 0.0;

      R[8]  = 0.0;
      R[9]  = 0.0;
      R[10] = 1.0;
      R[11] = 0.0;

      R[12] = 0.0;
      R[13] = 0.0;
      R[14] = 0.0;
      R[15] = 1.0;

      return;
   }

   w[0]/=wn;
   w[1]/=wn;
   w[2]/=wn;

   ctheta = cosf( theta );
   stheta = sinf( theta );
   vtheta = 1.0 - ctheta;

   R[0] = w[0]*w[0]*vtheta + ctheta; 
   R[1] = w[0]*w[1]*vtheta - w[2]*stheta; 
   R[2] = w[0]*w[2]*vtheta + w[1]*stheta; 
   R[3] = 0.0;

   R[4] = w[0]*w[1]*vtheta + w[2]*stheta; 
   R[5] = w[1]*w[1]*vtheta + ctheta; 
   R[6] = w[1]*w[2]*vtheta - w[0]*stheta; 
   R[7] = 0.0;

   R[8]  = w[0]*w[2]*vtheta - w[1]*stheta; 
   R[9]  = w[1]*w[2]*vtheta + w[0]*stheta; 
   R[10] = w[2]*w[2]*vtheta + ctheta; 
   R[11] = 0.0;

   R[12] = 0.0;
   R[13] = 0.0;
   R[14] = 0.0;
   R[15] = 1.0;

   return;
}


float *rotate(float alpha, float x, float y, float z)
{
   float CosAlpha,SinAlpha;
   float dum;
   float ax,ay,az;
   float *T;

   T=(float *)calloc(16,sizeof(float));
   T[15]=1.0;
   
   if(x==0.0 && y==0.0 && z==0.0)
   {
      T[0]=T[5]=T[10]=1.0;    
      return(T);
   }

   CosAlpha=(float)cos((float8)alpha);
   SinAlpha=(float)sin((float8)alpha);

   /* find the unit vector (ax,ay,az) in the direction of (x,y,z) */
   dum=(float)sqrt((float8)x*x+y*y+z*z);
   ax=x/dum; ay=y/dum; az=z/dum;
   
   T[0] = ax*ax + CosAlpha - CosAlpha*ax*ax; 
   T[1] = ax*ay - CosAlpha*ax*ay - SinAlpha*az; 
   T[2] = ax*az - CosAlpha*ax*az + SinAlpha*ay; 

   T[4] = ay*ax - CosAlpha*ay*ax + SinAlpha*az; 
   T[5] = ay*ay + CosAlpha - CosAlpha*ay*ay; 
   T[6] = ay*az - CosAlpha*ay*az - SinAlpha*ax; 

   T[8]  = az*ax - CosAlpha*az*ax - SinAlpha*ay; 
   T[9]  = az*ay - CosAlpha*az*ay + SinAlpha*ax; 
   T[10] = az*az + CosAlpha - CosAlpha*az*az; 

   return(T);
}

void rotate(float *R, float CosAlpha, float SinAlpha, float x, float y, float z)
{
   float dum;
   float ax,ay,az;

   R[0]=R[5]=R[10]=R[15]=1.0;
   R[1]=R[2]=R[3]=0.0;
   R[4]=R[6]=R[7]=0.0;
   R[8]=R[9]=R[11]=0.0;
   R[12]=R[13]=R[14]=0.0;
      
   if(x==0.0 && y==0.0 && z==0.0) return;
   
   /* find the unit vector (ax,ay,az) in the direction of (x,y,z) */
   dum=(float)sqrt((float8)x*x+y*y+z*z);
   ax=x/dum; ay=y/dum; az=z/dum;
      
   R[0] = ax*ax + CosAlpha - CosAlpha*ax*ax; 
   R[1] = ax*ay - CosAlpha*ax*ay - SinAlpha*az; 
   R[2] = ax*az - CosAlpha*ax*az + SinAlpha*ay; 
   
   R[4] = ay*ax - CosAlpha*ay*ax + SinAlpha*az; 
   R[5] = ay*ay + CosAlpha - CosAlpha*ay*ay; 
   R[6] = ay*az - CosAlpha*ay*az - SinAlpha*ax; 
   
   R[8]  = az*ax - CosAlpha*az*ax - SinAlpha*ay; 
   R[9]  = az*ay - CosAlpha*az*ay + SinAlpha*ax; 
   R[10] = az*az + CosAlpha - CosAlpha*az*az; 
   
   return;
}
