#include <stdlib.h>
//#include <iostream.h>
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
   {"-from", 1, 'f'},
   {"-to", 1, 't'},
   {"-h", 0, 'h'},
   {"-pixdim4", 1, '4'},
   {"-pixdim5", 1, '5'},
   {"-pixdim6", 1, '6'},
   {"-pixdim7", 1, '7'},
   {0, 0, 0}
};

int opt_o=NO;
int opt_pixdim4=NO;
int opt_pixdim5=NO;
int opt_pixdim6=NO;
int opt_pixdim7=NO;

void print_help_and_exit()
{
   printf("Usage: cphdr -from <nifti image> -to <nifti image>\n");
   exit(0);
}

int main(int argc, char **argv)
{
   char outputfile[1024]="";
   char inputfile[1024]="";
   nifti_1_header input_hdr;
   nifti_1_header output_hdr;
   float pixdim4=0.0, pixdim5=0.0, pixdim6=0.0, pixdim7=0.0;

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) 
      {
         case 't':
            sprintf(outputfile,"%s",optarg);
            break;
         case 'f':
            sprintf(inputfile,"%s",optarg);
            break;
         case 'h':
            print_help_and_exit();
            break;
         case '4':
            pixdim4 = atof(optarg);
            opt_pixdim4=YES;
            break;
         case '5':
            pixdim5 = atof(optarg);
            opt_pixdim5=YES;
            break;
         case '6':
            pixdim6 = atof(optarg);
            opt_pixdim6=YES;
            break;
         case '7':
            pixdim7 = atof(optarg);
            opt_pixdim7=YES;
            break;
         case '?':
            print_help_and_exit();
      }
   }

   if(argc==1) print_help_and_exit();

   if(inputfile[0]=='\0')
   {
      printf("Please specify an input file using -from argument.\n");
      exit(0);
   }

   if(outputfile[0]=='\0')
   {
      printf("Please specify an output file using -to argument.\n");
      exit(0);
   }

   char *input_hdr_segment;
   char *output_hdr_segment;
   char *output_data_segment;

   int input_hdr_size;
   int output_hdr_size, output_data_size;

   FILE *input_fp;
   FILE *output_fp;

   struct stat output_buf; // variable - file stats returned by stat()

   nifti_1_header *input_hdr_ptr;
   nifti_1_header *output_hdr_ptr;

   printf("Input file: %s\n", inputfile);
   printf("Output file: %s\n", outputfile);

   input_hdr = read_NIFTI_hdr(inputfile);
   output_hdr = read_NIFTI_hdr(outputfile);

   printf("ACi = %f\n",input_hdr.pixdim[4]);
   printf("PCi = %f\n",input_hdr.pixdim[5]);
   printf("ACx = %f\n",input_hdr.pixdim[6]);
   printf("PCx = %f\n",input_hdr.pixdim[7]);

   input_hdr_size = (int)input_hdr.vox_offset;
   output_hdr_size = (int)output_hdr.vox_offset;

   printf("Input header size = %d\n",input_hdr_size);
   printf("Output header size = %d\n",output_hdr_size);

   input_hdr_segment = (char *)calloc(input_hdr_size, 1);
   input_fp = fopen(inputfile,"r");
   if(input_fp==NULL) file_open_error(inputfile);
   fread(input_hdr_segment, 1, input_hdr_size, input_fp);
   fclose(input_fp);
   input_hdr_ptr = (nifti_1_header *)input_hdr_segment;

   if(opt_pixdim4)
   {
      input_hdr_ptr->pixdim[4]=pixdim4;
   }

   if(opt_pixdim5)
   {
      input_hdr_ptr->pixdim[5]=pixdim5;
   }

   if(opt_pixdim6)
   {
      input_hdr_ptr->pixdim[6]=pixdim6;
   }

   if(opt_pixdim7)
   {
      input_hdr_ptr->pixdim[7]=pixdim7;
   }

   stat(outputfile, &output_buf);
   output_data_size = output_buf.st_size - output_hdr_size;
   output_hdr_segment = (char *)calloc(output_hdr_size, 1);
   output_data_segment = (char *)calloc(output_data_size, 1);

   printf("Output data size = %d\n",output_data_size);

   output_fp = fopen(outputfile,"r");
   if(output_fp==NULL) file_open_error(outputfile);
   fread(output_hdr_segment, 1, output_hdr_size, output_fp);
   fread(output_data_segment, 1, output_data_size, output_fp);
   fclose(output_fp);
   output_hdr_ptr = (nifti_1_header *)output_hdr_segment;

   output_fp = fopen(outputfile,"w");
   if(output_fp==NULL) file_open_error(outputfile);
   fwrite(input_hdr_segment, 1, input_hdr_size, output_fp);
   fwrite(output_data_segment, 1, output_data_size, output_fp);
   fclose(output_fp);

   return 0;
}

