float *compute_L_inverse(float *x, float *y, int n)
{
   float *L;
   float *Li;
   
   int d;
   d = n+3;
   L = (float *)calloc(d*d, sizeof(float));
   Li = (float *)calloc(d*d, sizeof(float));

   ////////////////////////////////////////////////
   // set K
   ////////////////////////////////////////////////
   float *K;  // nxn

   K = (float *)calloc(n*n, sizeof(float));

   for(int i=0; i<n; i++)
   {
      for(int j=0; j<=i; j++)
      {
         K[i*n + j] = tpsU(x[i]-x[j],y[i]-y[j]);
      }
   }

   for(int j=0; j<n; j++)
   {
      for(int i=0; i<j; i++)
      {
         K[i*n + j] =  K[j*n + i];
      }
   }

   //printMatrix(K,n,n,"K:",NULL);  // for testing only
   ////////////////////////////////////////////////

   ////////////////////////////////////////////////
   // compute Ki=inverse(K)
   ////////////////////////////////////////////////
   float *Ki; // inverse(K)  nxn

   Ki = (float *)calloc(n*n, sizeof(float));

   for(int i=0; i<n*n; i++) 
   {
      Ki[i]=K[i]; // copy K in Ki temporarily
   }

   invert_symmetric_matrix(Ki, n);

   //multi(Ki,n,n,K,n,n,Ki);  // for testing only
   //printMatrix(Ki,n,n,"K*invserse(K):",NULL);  // for testing only
   ////////////////////////////////////////////////

   ////////////////////////////////////////////////
   // set P and Pt=transpose(P)
   ////////////////////////////////////////////////
   float *P;  // nx3
   float *Pt; // transpose(P) 3xn

   P = (float *)calloc(n*3, sizeof(float));

   for(int i=0; i<n; i++)
   {
      P[3*i] = 1.0;
      P[3*i + 1] = x[i];
      P[3*i + 2] = y[i];
   }

   Pt = trans(P, n, 3);

   //printMatrix(P,n,3,"P:",NULL);  // for testing only
   //printMatrix(Pt,3,n,"Pt:",NULL);  // for testing only
   ////////////////////////////////////////////////

   //////////////////////////////////////////////////////
   // compute B and Bt=transpose(B) where B=inverse(K)*P
   //////////////////////////////////////////////////////
   float *B; // B=inverse(K)*P
   float *Bt; // transpose(B)

   B = (float *)calloc(n*3, sizeof(float));

   multi(Ki,n,n,P,n,3,B); // compute B=inverse(K)*P 
   Bt = trans(B, n, 3);

   //printMatrix(B,n,3,"B:",NULL);  // for testing only
   //printMatrix(Bt,3,n,"Bt:",NULL);  // for testing only
   //////////////////////////////////////////////////////

   //////////////////////////////////////////////////////
   // compute A and Ai=inverse(A) 
   // where A=transpose(P)*inverse(K)*P
   //////////////////////////////////////////////////////
   float *A;
   float *Ai;

   A = (float *)calloc(3*3, sizeof(float));

   multi(Pt,3,n,B,n,3,A); // compute PtKiP=transpose(P)inverse(K)*P 
   Ai = inv3(A);

   //printMatrix(A,3,3,"A:",NULL);  // for testing only
   //printMatrix(Ai,3,3,"Ai:",NULL);  // for testing only
   //multi(Ai,3,3,A,3,3,Ai);  // for testing only
   //printMatrix(Ai,3,3,"A*invserse(A):",NULL);  // for testing only
   //////////////////////////////////////////////////////

   //////////////////////////////////////////////////////
   // compute BAi = B*inverse(A)
   // and AiBt = inverse(A)*transpose(B)
   //////////////////////////////////////////////////////
   float *BAi;
   float *AiBt;

   BAi = (float *)calloc(n*3, sizeof(float));
   multi(B,n,3,Ai,3,3,BAi);
   AiBt = trans(BAi, n, 3);

   //printMatrix(BAi,n,3,"BAi:",NULL);  // for testing only
   //printMatrix(AiBt,3,n,"AiBt:",NULL);  // for testing only
   //////////////////////////////////////////////////////

   //////////////////////////////////////////////////////
   // compute C
   //////////////////////////////////////////////////////
   float *C;
   C = (float *)calloc(n*n, sizeof(float));

   multi(BAi,n,3,Bt,3,n,C);

   //printMatrix(C,n,n,"C:",NULL);  // for testing only
   //////////////////////////////////////////////////////

   for(int i=0; i<n; i++)
   for(int j=0; j<n; j++)
   {
      L[i*d + j] = K[i*n + j];

      Li[i*d + j] = Ki[i*n + j] - C[i*n + j];
   }

   for(int i=0; i<n; i++)
   {
      L[i*d + n] = P[i*3];
      L[i*d + n + 1] = P[i*3 + 1];
      L[i*d + n + 2] = P[i*3 + 2];

      Li[i*d + n] = BAi[i*3];
      Li[i*d + n + 1] = BAi[i*3 + 1];
      Li[i*d + n + 2] = BAi[i*3 + 2];
   }

   for(int j=0; j<n; j++)
   {
      L[n*d + j] = Pt[j];
      L[(n+1)*d + j] = Pt[n + j];
      L[(n+2)*d + j] = Pt[2*n + j];

      Li[n*d + j] = AiBt[j];
      Li[(n+1)*d + j] = AiBt[n + j];
      Li[(n+2)*d + j] = AiBt[2*n + j];
   }

   for(int i=0; i<3; i++)
   for(int j=0; j<3; j++)
   {
      Li[(i+n)*d + j+n] = -Ai[i*3 + j];
   }

   //printMatrix(L,d,d,"L:",NULL);  // for testing only
   //printMatrix(Li,d,d,"Li:",NULL);  // for testing only
   //multi(L,d,d,Li,d,d,Li);
   //printMatrix(Li,d,d,"I:",NULL);  // for testing only

   free(K); free(Ki);
   free(P); free(Pt);
   free(B); free(Bt);
   free(A); free(Ai);
   free(BAi); free(AiBt);
   free(C);
   free(L);

   return(Li);
}
