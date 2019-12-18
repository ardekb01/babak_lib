#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>  //      required by time(), stat()
#include <time.h>       //      required by time()
#include <sys/stat.h>   //      required by stat() 
#include <unistd.h>
#include <spm_analyze.h>
#include <babak_lib.h>

#define YES 1
#define NO 0

int opt;

static struct option options[] =
{
        {"-scale", 1, 'S'},
        {"-threshold", 1, 't'},
        {"-n", 1, 'n'},
        {"-o", 1, 'o'},
        {"-s", 0, 's'},
        {"-b", 0, 'b'},
        {"-h", 0, 'h'},
        {"-help", 0, 'h'},
        {0, 0, 0}
};

int opt_o=NO;
int opt_b=NO;

void print_help_and_exit()
{
   printf("\nUsage: avgImage [-h -help -b -scale <scale> -s -threshold <threshold>] -o <output image name> <input images ...>\n"
   "-b: Binarizes (0 or 1) input images before averaging\n"
   "-h or -help: Prints help message\n" 
   "-o <output image name>: Specifies the filename for the outputted average image\n\n"
   "-scale <scale>: Multiplies the average image by <scale>\n\n"
   "-s: Save the output image as type short\n\n"
   "-threshold <threshold>: Binarizes the average images using <threhsold> level\n\n"
   "\n\n");

   exit(0);
}

float *avg(int N, char **imagefile)
{
   nifti_1_header hdr;

   int nx, ny, nz;
   float dx, dy, dz;
   int nv;
   int type;
   char *image;
   float *avg_image;

   if(N == 0) return(NULL);

   image = read_nifti_image(imagefile[0], &hdr);
   if(image==NULL) return(NULL);
   nx=hdr.dim[1]; ny=hdr.dim[2]; nz=hdr.dim[3];
   dx=hdr.pixdim[1]; dy=hdr.pixdim[2]; dz=hdr.pixdim[3];

   nv = nx*ny*nz;
   avg_image = (float *)calloc(nv, sizeof(float));

   type=hdr.datatype;
   switch(type) 
   {
		case 2:
            if(opt_b)
            {
			   for(int i=0; i<nv; i++) 
               {
                  if( ((unsigned char *)image)[i] > 0 )
                     ((unsigned char *)image)[i] = 1;
                  else
                     ((unsigned char *)image)[i] = 0;
               }
            }
			for(int i=0; i<nv; i++) avg_image[i] = ((unsigned char *)image)[i];
			break;
		case 4:
            if(opt_b)
            {
			   for(int i=0; i<nv; i++) 
               {
                  if( ((short *)image)[i] > 0 )
                     ((short *)image)[i] = 1;
                  else
                     ((short *)image)[i] = 0;
               }
            }
			for(int i=0; i<nv; i++) avg_image[i] = ((short *)image)[i];
			break;
		case 8:
            if(opt_b)
            {
			   for(int i=0; i<nv; i++) 
               {
                  if( ((int *)image)[i] > 0 )
                     ((int *)image)[i] = 1;
                  else
                     ((int *)image)[i] = 0;
               }
            }
			for(int i=0; i<nv; i++) avg_image[i] = ((int *)image)[i];
			break;
		case 16:
            if(opt_b)
            {
			   for(int i=0; i<nv; i++) 
               {
                  if( ((float *)image)[i] > 0 )
                     ((float *)image)[i] = 1;
                  else
                     ((float *)image)[i] = 0;
               }
            }
			for(int i=0; i<nv; i++) avg_image[i] = ((float *)image)[i];
			break;
		case 64:
            if(opt_b)
            {
			   for(int i=0; i<nv; i++) 
               {
                  if( ((double *)image)[i] > 0 )
                     ((double *)image)[i] = 1;
                  else
                     ((double *)image)[i] = 0;
               }
            }
			for(int i=0; i<nv; i++) avg_image[i] = ((double *)image)[i];
			break;
		case 512:
            if(opt_b)
            {
			   for(int i=0; i<nv; i++) 
               {
                  if( ((unsigned short *)image)[i] > 0 )
                     ((unsigned short *)image)[i] = 1;
                  else
                     ((unsigned short *)image)[i] = 0;
               }
            }
			for(int i=0; i<nv; i++) avg_image[i] = ((unsigned short *)image)[i];
			break;
		case 256:
            if(opt_b)
            {
			   for(int i=0; i<nv; i++) 
               {
                  if( ((char *)image)[i] > 0 )
                     ((char *)image)[i] = 1;
                  else
                     ((char *)image)[i] = 0;
               }
            }
			for(int i=0; i<nv; i++) avg_image[i] = ((char *)image)[i];
			break;
   }

   delete image;

   for(int i=1; i<N; i++)
   {
      image = read_nifti_image(imagefile[i], &hdr);
      if(image==NULL) return(NULL);

      type=hdr.datatype; // this statement was missing; major bug fixed Aug. 8, 2012
      switch(type) {
         case 2:
            if(opt_b)
            {
			      for(int j=0; j<nv; j++) 
                  {
                     if( ((unsigned char *)image)[j] > 0 )
                        ((unsigned char *)image)[j] = 1;
                     else
                        ((unsigned char *)image)[j] = 0;
                  }
            }
            for(int j=0; j<nv; j++) avg_image[j] += ((unsigned char *)image)[j];
            break;
         case 4:
            if(opt_b)
            {
			      for(int j=0; j<nv; j++) 
                  {
                     if( ((short *)image)[j] > 0 )
                        ((short *)image)[j] = 1;
                     else
                        ((short *)image)[j] = 0;
                  }
            }
            for(int j=0; j<nv; j++) avg_image[j] += ((short *)image)[j];
            break;
         case 8:
            if(opt_b)
            {
			      for(int j=0; j<nv; j++) 
                  {
                     if( ((int *)image)[j] > 0 )
                        ((int *)image)[j] = 1;
                     else
                        ((int *)image)[j] = 0;
                  }
            }
            for(int j=0; j<nv; j++) avg_image[j] += ((int *)image)[j];
            break;
         case 16:
            if(opt_b)
            {
			      for(int j=0; j<nv; j++) 
                  {
                     if( ((float *)image)[j] > 0 )
                        ((float *)image)[j] = 1;
                     else
                        ((float *)image)[j] = 0;
                  }
            }
            for(int j=0; j<nv; j++) avg_image[j] += ((float *)image)[j];
            break;
          case 64:
            if(opt_b)
            {
			      for(int j=0; j<nv; j++) 
                  {
                     if( ((double *)image)[j] > 0 )
                        ((double *)image)[j] = 1;
                     else
                        ((double *)image)[j] = 0;
                  }
            }
            for(int j=0; j<nv; j++) avg_image[j] += ((double *)image)[j];
            break;
         case 512:
            if(opt_b)
            {
			      for(int j=0; j<nv; j++) 
                  {
                     if( ((unsigned short *)image)[j] > 0 )
                        ((unsigned short *)image)[j] = 1;
                     else
                        ((unsigned short *)image)[j] = 0;
                  }
            }
            for(int j=0; j<nv; j++) avg_image[j] += ((unsigned short *)image)[j];
            break;
         case 256:
            if(opt_b)
            {
			      for(int j=0; j<nv; j++) 
                  {
                     if( ((char *)image)[j] > 0 )
                        ((char *)image)[j] = 1;
                     else
                        ((char *)image)[j] = 0;
                  }
            }
            for(int j=0; j<nv; j++) avg_image[j] += ((char *)image)[j];
            break;
      }

      delete image;
   }

   for(int i=0; i<nv; i++) avg_image[i] /= N;

   return(avg_image);
}

float *avg4d(char *imagefile, int n)
{
   nifti_1_header hdr;

   int nx, ny, nz, nt;
   float dx, dy, dz;
   int nv;
   char *image;
   float *avg_image;

   image = read_nifti_image(imagefile, &hdr);
   if(image==NULL) 
   {
      return(NULL);
   }

   nx=hdr.dim[1]; ny=hdr.dim[2]; nz=hdr.dim[3]; nt=hdr.dim[4];
   dx=hdr.pixdim[1]; dy=hdr.pixdim[2]; dz=hdr.pixdim[3];

   if(n<=0) n=nt;

   printf("Matrix size = %d x %d x %d x %d\n", nx, ny, nz, nt);
   printf("Voxel size = %f x %f x %f\n", dx, dy, dz);
   printf("Averaging %d images ...\n",n);

   if(nt == 0)
   {
      printf("nt cannot be equal to 0, aborting ...\n");
      exit(0);
   }

   nv = nx*ny*nz;
   avg_image = (float *)calloc(nv, sizeof(float));
   if(avg_image == NULL)
   {
      memory_allocation_error("\"avg_image\" in avg4d()");
   }

   for(int v=0; v<nv; v++) avg_image[v]=0.0;

   for(int i=0; i<n; i++)
   {
      switch(hdr.datatype) {
         case 2:
            for(int v=0; v<nv; v++) avg_image[v] += ((unsigned char *)image)[i*nv + v];
            break;
         case 4:
            for(int v=0; v<nv; v++) avg_image[v] += ((short *)image)[i*nv + v];
            break;
         case 8:
            for(int v=0; v<nv; v++) avg_image[v] += ((int *)image)[i*nv + v];
            break;
         case 16:
            for(int v=0; v<nv; v++) avg_image[v] += ((float *)image)[i*nv + v];
            break;
         case 64:
            for(int v=0; v<nv; v++) avg_image[v] += ((double *)image)[i*nv + v];
            break;
         case 512:
            for(int v=0; v<nv; v++) avg_image[v] += ((unsigned short *)image)[i*nv + v];
            break;
         case 256:
            for(int v=0; v<nv; v++) avg_image[v] += ((char *)image)[i*nv + v];
            break;
      }
   }

   delete image;

   for(int v=0; v<nv; v++) avg_image[v] /= n;

   return(avg_image);
}

// returns 1 if all images have the same dimensions nx, ny, and nz, 0 otherwise
int checkDimension_avgImage(int N, char **imagefile, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz)
{
   nifti_1_header hdr;
   short dataType;

   if(N==0) return(1);

   printf("Image %d: %s\n",1,imagefile[0]);

   hdr = read_NIFTI_hdr(imagefile[0]);
   *nx = hdr.dim[1];
   *ny = hdr.dim[2];
   *nz = hdr.dim[3];
   *dx = hdr.pixdim[1];
   *dy = hdr.pixdim[2];
   *dz = hdr.pixdim[3];
   dataType = hdr.datatype;

   for(int i=1; i<N; i++)
   {
      printf("Image %d: %s\n",i+1,imagefile[i]);
      hdr = read_NIFTI_hdr(imagefile[i]);

      if( *nx != hdr.dim[1] ||  *ny != hdr.dim[2] ||  *nz != hdr.dim[3]) 
      {
            return(0);
      }
   }

   return(1);
}

int main(int argc, char **argv)
{
   int n=0; // number of images from the top to be averaged if the input image is 4D
   int opt_s=NO;
   int opt_threshold=NO;

   float scale=1.0;
   float threshold=0.0;

   nifti_1_header hdr;
   int nx, ny, nz, nv;
   float dx, dy, dz;
   int number_of_images; // number of images to be averaged
   float *avg_image;
   char outputfile[1024];

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) 
      {
         case 't':
            threshold = atof(optarg);
            opt_threshold=YES;
            break;
         case 'S':
            scale = atof(optarg);
            break;
         case 'o':
            sprintf(outputfile,"%s",optarg);
            opt_o=YES;
            break;
         case 'h':
            print_help_and_exit();
            break;
         case 'b':
            opt_b=YES;
            break;
         case 's':
            opt_s=YES;
            break;
         case 'n':
            n = atoi(optarg);
            break;
         case '?':
            print_help_and_exit();
      }
   }

   if(argc==1) print_help_and_exit();

   if(!opt_o)
   {
      printf("Please specify an output file using the -o argument.\n");
      exit(0);
   }

   if(opt_b)
   {
      printf("Binarlization = ON\n");
   }

   number_of_images = argc-optind;

   printf("Number of input images = %d\n",number_of_images);

   if(number_of_images==0) exit(0);

   if( !checkDimension_avgImage(number_of_images, argv+optind, &nx, &ny, &nz, &dx, &dy, &dz) )
   {
      printf("\n\nAll input images must be the same size. Aborting ...");
      printf("\n\n");
      exit(0);
   }

   hdr = read_NIFTI_hdr( (argv+optind)[0]);

   if( hdr.dim[0] == 4 && number_of_images == 1)
   {
      avg_image = avg4d( (argv+optind)[0], n);
      hdr.dim[0]=3;
   }
   else
   {
      avg_image = avg(number_of_images, argv+optind);
   }

   nv = nx*ny*nz;
   if(scale != 1.0)
   {
      printf("Scale = %f\n",scale);
      for(int v=0; v<nv; v++) avg_image[v] *= scale;
   }

   if(opt_threshold)
   {
      printf("\nThreshold = %f\n", threshold);
      for(int v=0; v<nv; v++) 
      {
         if( avg_image[v]>threshold ) 
         {
            avg_image[v]=1.0; 
         }
         else 
         {
            avg_image[v]=0.0;
         }
      }
   }
 
   if(opt_s == YES)
   {
      short *tmp;
      tmp = (short *)calloc(nv, sizeof(short));
      for(int i=0; i<nv; i++) tmp[i] = (short)(avg_image[i] + 0.5);
      save_nifti_image(outputfile, tmp, &hdr);
      free(tmp);
   }
   else
   {
      save_nifti_image(outputfile, avg_image, &hdr);
   }

   return 0;
}
