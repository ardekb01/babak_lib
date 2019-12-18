//NOTES: subject image is expect to be 255*255*189, 1*1*1 mm^3 PIL
#include <babak_lib.h>
#include <sph.h>
#include <landmarks.h>

//////////////////////////////////////////////////////////////////////////////////////////////////
int opt;

static struct option options[] =
{
   {"-m",1,'m'},
   {"-i",1,'i'},
   {"-o",1,'o'},
   {"-v",0,'v'},
   {0, 0,  0}
};

int main(int argc, char **argv)
{
   FILE *fp;
   int nl;
   char mdlfile[256]=""; 
   char subfile[256]=""; 
   char outputfile[256]=""; 
   float *P;

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) {
         case 'v':
            opt_v=YES;
            break;
         case 'm':
            sprintf(mdlfile,"%s",optarg);
            break;
         case 'i':
            sprintf(subfile,"%s",optarg);
            break;
         case 'o':
            sprintf(outputfile,"%s",optarg);
            break;
         case '?':
            exit(0);
      }
   }

   ///////////////////////////////////////////////////////////////////////////////////////

   if( subfile[0]=='\0' )
   {
      printf("Please specify a subject image using -i argument.\n");
      exit(0);
   }

   if( outputfile[0]=='\0' )
   {
      printf("Please specify an output filename using -o argument.\n");
      exit(0);
   }

   if( mdlfile[0]=='\0' )
   {
      printf("Please specify a landmark model using -m argument.\n");
      exit(0);
   }

   ///////////////////////////////////////////////////////////////////////////////////////
   if(opt_v)
   {
      printf("Subject image: %s\n",subfile);
      printf("Output file: %s\n",outputfile);
      printf("Landmark model file: %s\n",mdlfile);
   }

   P=detect_landmarks( (const char *)subfile, (const char *)mdlfile, nl, 0);

   if(opt_v)
   {
      printf("Number of landmarks = %d\n",nl);
   }

   fp = fopen(outputfile,"w");
   fprintf(fp,"%d\n",nl);
   for(int L=0; L<nl; L++)
   {
      fprintf(fp,"%d %d %d\n", (int)P[0*nl + L], (int)P[1*nl + L], (int)P[2*nl + L]);
      //printf("Landmark %d detected at: %5.1f %5.1f %5.1f\n", L+1, P[0*nl + L], P[1*nl + L], P[2*nl + L]);
   }
   fclose(fp);

   free(P);
}
