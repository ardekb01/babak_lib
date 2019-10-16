#include <babak_lib.h>
#include <sph.h>
#include <landmarks.h>

//////////////////////////////////////////////////////////////////////////////////////////////////
int opt;

static struct option options[] =
{
   {"-m",1,'m'},
   {"-i",1,'i'},
   {"-v",0,'v'},
   {0, 0,  0}
};

int main(int argc, char **argv)
{
   int nl;
   char mdlfile[256]=""; 
   char subfile[256]=""; 
   float *P;

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) {
         case 'm':
            sprintf(mdlfile,"%s",optarg);
            break;
         case 'i':
            sprintf(subfile,"%s",optarg);
            break;
         case '?':
            exit(0);
      }
   }

   ///////////////////////////////////////////////////////////////////////////////////////

   if( subfile[0]=='\0' )
   {
      printf("Please specify a subject image using -i option.\n");
      exit(0);
   }

   if( mdlfile[0]=='\0' )
   {
      printf("Please specify a landmark model using -m option.\n");
      exit(0);
   }

   ///////////////////////////////////////////////////////////////////////////////////////
   printf("Subject image: %s\n",subfile);
   printf("Landmark model file: %s\n",mdlfile);

   P=detect_landmarks( (const char *)subfile, (const char *)mdlfile, nl);
   printf("Number of landmarks = %d\n",nl);

   for(int l=0; l<nl; l++)
   {
      printf("Landmark %d detected at: %5.1f %5.1f %5.1f\n", l+1, P[0*nl + l], P[1*nl + l], P[2*nl + l]);
   }
}
