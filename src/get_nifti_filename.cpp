#include "check_nifti_file_extension.h"
#include "check_nifti1_magic.h"
#include <cstdio> 
#include <cstring> 

// extracts the "filename" from the full "path" string
// Example: If path="/home/babak/images/test.nii", then filename="test".
// Note that the extension is not included in the output filename
// Returns 0 on failture and 1 on success
bool get_nifti_filename(char *filename, const char *path)
{
   int i;
   int len;	// length of the path string
   int pos;	// position of the filename

   // ensure that the specified image has .nii extension
   if( check_nifti_file_extension(path) == false )
   {
      printf("%s does not have `.nii' extension\n",path);
      return false;
   }

   if( check_nifti1_magic(path) == false )
   {
      return false;
   }

   len=(int)strlen(path);

   if(len<=0)
   {
      printf("Error: unexpected string length for the NIFTI image path, aborting ...\n");
      return false;
   }

   // finds the position of the first '/' character from right if any
   // returns 0 of no '/' charater is there
   // Examples: path=/sss/yyy would give pos=5
   // path=sss would give pos=0
   // path=/x/ would give pos=3
   i=len-1;
   while( i>=0 && path[i] != '/' )
   {
      i--;
   }
   pos=i+1;

   // copy everything to the right of the first '/' character from right
   // into filename
   // Examples: path=/sss/yyy would give filename=yyy  len=3
   // path=sss would give filename=sss  len=3
   // path=/x/ would give filename=""  len=0
   strcpy(filename,path+pos);
   len=(int)strlen(filename);

   if(len<=0)
   {
      printf("Error: unexpected string length for the NIFTI image filename, aborting ...\n");
      return false;
   }

   if( len>=2 && filename[len-2]=='g' && filename[len-1]=='z' )
   {
      printf("Sorry but this program currently does not handle gzipped images, aborting ...\n");
      return false;
   }

   // finds the position of the first '.' character from right if any
   // returns 0 of no '.' charater is there
   // Examples: path=/sss.yyy would give pos=5
   // path=sss would give pos=0
   // path=/x. would give pos=3
   i=len-1;
   while( i>=0 && filename[i] != '.' )
   {
      i--;
   }
   pos=i+1;

   if(pos>0) filename[pos-1]='\0';

   len=(int)strlen(filename);

   if(len<=0)
   {
      printf("Error: unexpected string length for the NIFTI image prefix, aborting ...\n");
      return false;
   }

   return true;
}
