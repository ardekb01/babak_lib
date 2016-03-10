/////////////////////////////////////////////////////////////////////////////////////////////
void intensity_norm(double *x0, double *x1, double *w, double *w_const, int n, double *u, double &d, double &x0a, double &x1a)
{

   double *x0c;  // x0 centered wrt x0a
   double *x1c;  // x1 centered wrt x1a
   double a11, a22, a12;
   double L; // the smaller eigenvalue of A

/*
   //////////////////////////////////////////////////////
   // Compute the weighted averages of x's and y's
   //If this part is commented the line will pass from the given point (x0a,x1a).   
   //////////////////////////////////////////////////////
   {
      double Sw;   // sum of weights

      Sw=x0a=x1a=0.0;
      for(int i=0; i<n; i++)
      {
         x0a += w[i]*w_const[i]*x0[i];
         x1a += w[i]*w_const[i]*x1[i];
         Sw += w[i]*w_const[i];
      }

      x0a /= Sw;
      x1a /= Sw;
   

      if(opt_v)   printf("x0a=%lf x1a=%lf Sw=%lf\n",x0a,x1a,Sw);
   }
*/

   // allocate memory for x0c and x1c
   x0c = (double *)calloc(n, sizeof(double));
   x1c = (double *)calloc(n, sizeof(double));

   // center x0 and x1 
   for(int i=0; i<n; i++)
   {
      x0c[i] = x0[i] - x0a;
      x1c[i] = x1[i] - x1a;
   }

   a11=a12=a22=0.0;
   for(int i=0; i<n; i++)
   {
      a11 += w[i]*w_const[i]*x0c[i]*x0c[i];
      a12 += w[i]*w_const[i]*x0c[i]*x1c[i];
      a22 += w[i]*w_const[i]*x1c[i]*x1c[i];
   }

   // to make the numbers more manageable 
   a11 /= n;
   a12 /= n;
   a22 /= n;

   {
      double a,b,c,D;
      double L1,L2;

      a=1;
      b=-(a11+a22);
      c=a11*a22 - a12*a12;
      D = b*b - 4*a*c;
      L1 = (-b + sqrt(D))/(2*a);
      L2 = (-b - sqrt(D))/(2*a);

      L=L2;
   }

//   if(opt_v)   printf("Smaller eigen-value = %lf\n",L);

   {
      double u1, u2;

      u1 = sqrt( (a22 - L)/(a11+a22-2*L) );
      u2 = sqrt( 1.0 - u1*u1);

      if( (int)(1000*L*u1)!= (int)(1000*a11*u1+a12*u2) || (int)(1000*L*u2)!=(int)(1000*a12*u1+a22*u2))
      {
         u2 = -u2;
      }

/* 
      if(opt_v)
      {
         printf("L*u^t=[%lf %lf]\n",L*u1, L*u2);
         printf("u^t*A=[%lf %lf]\n",a11*u1+a12*u2, a12*u1+a22*u2);
      }
*/
      d = u1*x0a + u2*x1a;

      u[0]=u1;
      u[1]=u2;

   }

   //////////////////////////////////////////////////////
   // update weights in preparation for the next iteration
   //////////////////////////////////////////////////////
   {
      double r;
      //Huber parameter
      for(int i=0; i<n; i++)
      {
         r = u[0]*x0c[i] + u[1]*x1c[i];
         r = fabs(r);

         // amounts to M-esitmator with Geman-McClure residuals
         w[i] = 1.0/((1.0+r*r)*(1.0+r*r));

      }
   }
   free(x0c);
   free(x1c);
}


//*********************************************************************
// "hist2D_line" plots the histogram of the baseline and follow-up image 
// intensities and finds the best fitting line.
// Input images should be registered in the PIL space.
// The PILbraincloud is used for weighting each voxel of the image.
// The program generates 3 files in the running directory in the following
// formats:  hist_%s_%s.mat 
//           hist_%s_%s.plt 
//           hist_%s_%s.png
//
// Inputs:
//    bfile: the baseline image file name
//    ffile: the follow-up image file name
//    PILbrainfile: the PIL brain cloud file name
// 
// Output:
//    slope: slope of the best fitting line
// 
double hist2D_line(const char *bfile, const char *ffile, const char *PILbrainfile)
{

   FILE *fp;
   char filename[1024]=""; // a generic filename for reading/writing stuff
   char command[1024]="";  // a generic filename for reading/writing stuff

   char bprefix[1024]="";  //baseline image prefix
   char fprefix[1024]="";  //follow-up image prefix

   double slope;           //The slope of the best fitting line
   double slope_tmp;


   if(opt_v)   printf("-----------------------------------\n");
   /////////////////////////////////////////////////////////////////////////////////////////////
   // Read PILbrainfile from the $ARTHOME directory
   short *PILbraincloud;
   DIM PILbraincloud_dim;
   nifti_1_header PILbraincloud_hdr;
   sprintf(filename,"%s/%s", ARTHOME, PILbrainfile);
   PILbraincloud = (short *)read_nifti_image(filename, &PILbraincloud_hdr);

   if(PILbraincloud==NULL)
   {
         printf("Error reading %s, aborting ...\n", filename);
         exit(1);
   }

   set_dim(PILbraincloud_dim, PILbraincloud_hdr);

   /////////////////////////////////////////////////////////////////////////////////////////////
   // Read the baseline and follow-up images
   if( niftiFilename(bprefix, bfile)==0 )  exit(0);
   if(opt_v)  printf("Baseline image prefix: %s\n",bprefix);

   if( niftiFilename(fprefix, ffile)==0 )  exit(0);
   if(opt_v)  printf("Follow-up image prefix: %s\n",fprefix);

   short* bim;          // baseline image
   short* fim;          // follow-up image

   DIM bdim;            // baseline image dimensions structure
   DIM fdim;            // follow-up image dimensions structure

   nifti_1_header bhdr; // baseline image NIFTI header
   nifti_1_header fhdr; // follow-up image NIFTI header

   bim = (short *)read_nifti_image(bfile, &bhdr);
   fim = (short *)read_nifti_image(ffile, &fhdr);

   if(bim==NULL)
   {
      printf("Error reading %s, aborting ...\n", bfile);
      exit(1);
   }
   if(fim==NULL)
   {
      printf("Error reading %s, aborting ...\n", ffile);
      exit(1);
   }

   set_dim(bdim, bhdr);
   set_dim(fdim, fhdr);
  
   
   /////////////////////////////////////////////////////////////////////////////////////////////
   // Check the image sizes to make sure they match the PILbraincloud
   if(bdim.nx!=PILbraincloud_dim.nx || bdim.ny!=PILbraincloud_dim.ny || bdim.nz!=PILbraincloud_dim.nz)
   {
      printf("Size of image %s does not match the PILbraincloud, aborting ...\n", bprefix);
      exit(1);
   }
   else if(fdim.nx!=PILbraincloud_dim.nx || fdim.ny!=PILbraincloud_dim.ny || fdim.nz!=PILbraincloud_dim.nz)
   {
      printf("Size of image %s does not match the PILbraincloud, aborting ...\n", fprefix);
      exit(1);
   }
   else
   {
      int highb=0, highf=0;
      int ntot,n;
      int *hist;

      double *x0, *x1, *w, *w_const;
      double x0a, x1a;         //The best fitting line is passing through [x0a,x1a] 
      double u[2], d;


      int high_square;
      int rangex,rangey;

      for(int i=0; i<PILbraincloud_dim.nv; i++)
      if(PILbraincloud[i]>0)
      {
         if(bim[i]>highb) highb=bim[i];
         if(fim[i]>highf) highf=fim[i];
      }

      hist=(int *)calloc((highb+1)*(highf+1), sizeof(int));

      for(int i=0; i<PILbraincloud_dim.nv; i++)
      if(PILbraincloud[i]>0)
      {
         hist[bim[i]+(highb+1)*fim[i]]+=PILbraincloud[i];
      }

      n=0;
      for(int i=0; i<(highb+1)*(highf+1); i++)
      if(hist[i]>100)   n++;
      ntot=n;

      for(int k=0; k<=highf+highb; k++)
      {
         n=0;
         for(int i=0; i<=k; i++)
         {
            for(int j=0; j<=k; j++)
            {
               if(i<highb && j<highf && hist[i+(highb+1)*j]>100) n++;
            }
         }

         if(n>(998*ntot)/1000)
         {
            high_square=k-1;
            break;
         }
      }

      rangex=rangey=high_square;
      if (high_square>highb) rangex=highb;
      if (high_square>highf) rangey=highf;

      if(opt_v)
      {
         printf("\n");
         printf("           Intensity_range  Display_range\n");
         printf("Baseline:     [0,%4d]         [0,%4d]\n",highb, rangex);
         printf("Follow-up:    [0,%4d]         [0,%4d]\n",highf, rangey);
      }

      n=(rangex+1)*(rangey+1);
      x0 = (double *)calloc(n, sizeof(double));
      x1 = (double *)calloc(n, sizeof(double));
      w = (double *)calloc(n, sizeof(double));
      w_const = (double *)calloc(n, sizeof(double));

      sprintf(filename,"hist_%s_%s.mat",bprefix,fprefix);
      fp = fopen(filename,"w");
      n=0;
      for(int j=0; j<=rangey; j++)
      {
         for(int i=0; i<=rangex; i++)
         {
           fprintf(fp,"%f ",(float)hist[i+(highb+1)*j]/100);
           x0[n]=(double)i;
           x1[n]=(double)j;
           w_const[n]=(double)hist[i+(highb+1)*j]/100;
           w[n]=1;
         
           n++;
         }
         fprintf(fp,"\n");
      }
      fclose(fp);

      //////////////////////////////////////////////////////////////////////
      //Calculating the best fitting line
      x0a=0;        // to pass through the origin
      x1a=0;

      for(int iter=0; iter<2000; iter++)
      {
         intensity_norm(x0, x1, w, w_const, n, u, d, x0a, x1a);
         slope=(-u[0]/u[1]);            //The slope of the best fitting line

         if(fabs(slope_tmp-slope)< 1.0e-5)       
         {
            if(opt_v) printf("\nNumber of iteration=%d  Slope=%lf\n",iter, slope);
            break;
         }
         slope_tmp=slope;
      }

      //////////////////////////////////////////////////////////////////////
      //Writing the plot file
      {
         sprintf(filename,"hist_%s_%s.plt",bprefix,fprefix);
         fp = fopen(filename,"w");
   
         fprintf(fp,"unset key\n");
         fprintf(fp,"plot 'hist_%s_%s.mat' matrix with image\n",bprefix,fprefix);

         fprintf(fp,"set lmargin screen 0.15\n");
         fprintf(fp,"set rmargin screen 0.85\n");
         fprintf(fp,"set tmargin screen 0.9\n");
         fprintf(fp,"set bmargin screen 0.15\n");

         fprintf(fp,"set xlabel \"Baseline intensity\"\n");
         fprintf(fp,"set ylabel \"Follow-up intensity\" \n");

         fprintf(fp,"if (GPVAL_DATA_Y_MAX+1 > GPVAL_DATA_X_MAX+1)\\\n");
         fprintf(fp,"    MAX=GPVAL_DATA_Y_MAX; \\\n");
         fprintf(fp,"else \\\n");
         fprintf(fp,"    MAX=GPVAL_DATA_X_MAX\n");
   
         if(rangex>rangey)
         {
            fprintf(fp,"set xrange [0:GPVAL_DATA_X_MAX]\n");
            fprintf(fp,"set yrange [0:GPVAL_DATA_X_MAX]\n");
         }
         else
         {
            fprintf(fp,"set xrange [0:GPVAL_DATA_Y_MAX]\n");
            fprintf(fp,"set yrange [0:GPVAL_DATA_Y_MAX]\n");
         }
         
         fprintf(fp,"set termoption dashed\n");
         fprintf(fp,"set style line 1 linecolor rgb '#8B008B' lt 1 lw 1\n");
         fprintf(fp,"set style line 2 linecolor rgb '#191970' lt 3 lw 1\n");
         fprintf(fp,"set style line 3 linecolor rgb '#00FF00' lt 3 lw 1\n");
         fprintf(fp,"set style line 4 linecolor rgb '#006400' lt 3 lw 1\n");

         fprintf(fp,"f(x) = (%lf - (%lf*x))/(%lf)\n",d,u[0],u[1]);
         fprintf(fp,"g(x) = x\n");
         
         fprintf(fp,"plot 'hist_%s_%s.mat' matrix with image, f(x) ls 1, g(x) ls 2\n",bprefix,fprefix);
  
         fprintf(fp,"set grid\n");
         fprintf(fp,"set style line 11 lc rgb '#252525' lt 1\n");
         fprintf(fp,"set tics nomirror out scale 0.75\n");
         fprintf(fp,"set palette defined ( 0 '#F7FBFF',\\\n");
         fprintf(fp,"                      1   '#3288BD',\\\n");
         fprintf(fp,"                      211 '#66C2A5',\\\n");
         fprintf(fp,"                      311 '#ABDDA4',\\\n");
         fprintf(fp,"                      411 '#E6F598',\\\n");
         fprintf(fp,"                      511 '#FEE08B',\\\n");
         fprintf(fp,"                      611 '#FDAE61',\\\n");
         fprintf(fp,"                      711 '#F46D43',\\\n");
         fprintf(fp,"                      811 '#D53E4F')\n");

         fprintf(fp,"set size ratio -1\n");
         fprintf(fp,"set terminal png font 'arial,28' enhanced nocrop size 2000,1600\n");
         //fprintf(fp,"set terminal png enhanced nocrop size 2000,1600\n");
         
         fprintf(fp,"set output 'hist_%s_%s.png'\n",bprefix,fprefix);
         fprintf(fp,"replot\n");
         fclose(fp);
         
         //////////////////////////////////////////////
         //Checking the .plt file
         sprintf(command,"test -f hist_%s_%s.plt && echo || echo \"Error in writing the plt file.\"",bprefix,fprefix);
         system(command);

         sprintf(command,"test -f hist_%s_%s.plt && echo \"\\\"hist_%s_%s.plt\\\" file generated.\" ",bprefix,fprefix,bprefix,fprefix);
         if(opt_v)   system(command);

         ////////////////////////////////////////////////
         //Running the gnuplot to generate the .png file.
         sprintf(command,"gnuplot hist_%s_%s.plt",bprefix,fprefix);
         system(command);

         ///////////////////////////////////////////////
         //Checking the .png file
         sprintf(command,"test -f hist_%s_%s.png && echo || echo \"Error in running the gnuplot.\"",bprefix,fprefix);
         system(command);

         sprintf(command,"test -f hist_%s_%s.png && echo \"\\\"hist_%s_%s.png\\\" file generated.\" ",bprefix,fprefix,bprefix,fprefix);
         if(opt_v)   system(command);
      }
    
      free(x0);
      free(x1);
      free(w);
      free(w_const);
      free(hist);

   }

   if(opt_v)   printf("-----------------------------------\n");

   return slope;
}

