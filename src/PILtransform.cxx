#include <babak_lib.h>
#include <sph.h>
#include <landmarks.h>

//////////////////////////////////////////////////////////////////////////////////////////////////
int opt;

static struct CmdOption options[] =
{
   {"-lm",1,'m'},
   {"-i",1,'i'},
   {"-v",0,'v'},
   {0, 0,  0}
};

void print_help_and_exit()
{
   printf("\nThis program computes a PIL transformation matrix.\n");
   printf("\nUsage:\n"
   "\tPILtransform [-v] [-lm <landmarks>] -i <image>.nii\n\n");
   exit(0);
}

int main(int argc, char **argv)
{
   char orient[4]=""; 
   char lmfile[256]=""; 
   char subfile[256]=""; 
   float TPIL[16];

   if(argc==1) print_help_and_exit();

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) {
         case 'm':
            snprintf(lmfile,sizeof(lmfile),"%s",optarg);
            break;
         case 'i':
            snprintf(subfile,sizeof(subfile),"%s",optarg);
            break;
         case 'v':
            opt_v=YES;
            break;
         case '?':
            exit(0);
      }
   }

   if( subfile[0] == '\0')
   {
      printf("\nPlease specify an inptu image using -i <filename.nii> ...\n\n");
      exit(0);
   }

   if(opt_v)
   {
      printf("Input image file: %s\n",subfile);
   }

   if(lmfile[0] != '\0' && opt_v) printf("Landmarks are read from %s\n",lmfile);

   orient[0] = '\0';
   new_PIL_transform(subfile,lmfile,orient,TPIL,1);
}
