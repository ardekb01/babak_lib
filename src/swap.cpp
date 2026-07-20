#include <stdlib.h>
#include <stdio.h>
#include <babak_lib.h>

// Returns true if the computer is big-endian (e.g., SUN).
// Returns false if the computer is little-endian (e.g., IBM PC).
bool bigEndian()
{
   short s;
   unsigned char *cp;

   // We will set s = 1. If the computer is big-endian, the first byte
   // in s will be 0 and the second byte will be 1. Therefore, this
   // function returns true if the second byte is 1.

   s = 1;

   cp = (unsigned char *)(&s);

   // For s = 1:
   // Big-endian:    bytes are 00 01 -> cp[1] == 1 -> true.
   // Little-endian: bytes are 01 00 -> cp[1] == 0 -> false.

   return cp[1] == 1;
}

void swapByteOrder(char *in, int4 N)
{
   char *dum;

   if(N<=0)
   {
      printf("\n\nWarning swapByteOrder(): Non-positive array dimension argument.\n\n");
      return;
   }

   if(in==NULL)
   {
      printf("\n\nWarning swapByteOrder(): NULL array argument.\n\n");
      return;
   }

   dum = (char *)malloc(N);

   if(dum==NULL)
   {
      printf("\n\nWarning swapByteOrder(): Memory allocation failure.\n\n");
      return;
   }

   for( int4 i=0; i<N; i++)
   {
      dum[i]=in[N-1-i];
   }

   for( int4 i=0; i<N; i++)
   {
      in[i]=dum[i];
   }

   free(dum);
}       
////////////////////////////////////////////////////////////////////////////////

void swapN(char *in, int4 N)
{
	char dum[2]; 

	for(int4 i=0; i<N/2; i++)
	{    
		dum[0]=in[i*2+1];
		dum[1]=in[i*2];
        
		in[i*2]=dum[0];
		in[i*2+1]=dum[1];
	}
}       

void swap_double_array( float8 *x, int4 n)
{
   for(int4 i=0; i<n; i++)
   {    
      swapByteOrder( (char *)(x+i), sizeof(float8) );
   }
}       

void swap_float_array( float4 *x, int4 n)
{
   for(int4 i=0; i<n; i++)
   {    
      swapByteOrder( (char *)(x+i), sizeof(float4) );
   }
}       

void swap_int_array( int4 *x, int4 n)
{
   for(int4 i=0; i<n; i++)
   {    
      swapByteOrder( (char *)(x+i), sizeof(int4) );
   }
}       

void swap_model_file_hdr(model_file_hdr *hdr)
{
   swapByteOrder( (char *)(&hdr->nxHR), sizeof(int4));
   swapByteOrder( (char *)(&hdr->nzHR), sizeof(int4));
   swapByteOrder( (char *)(&hdr->dxHR), sizeof(float4));
   swapByteOrder( (char *)(&hdr->nxLR), sizeof(int4));
   swapByteOrder( (char *)(&hdr->nvol), sizeof(int4));
   swapByteOrder( (char *)(&hdr->RPtemplateradius), sizeof(int4));
   swapByteOrder( (char *)(&hdr->RPtemplateheight), sizeof(int4));
   swapByteOrder( (char *)(&hdr->RPtemplatesize), sizeof(int4));
   swapByteOrder( (char *)(&hdr->ACtemplateradius), sizeof(int4));
   swapByteOrder( (char *)(&hdr->ACtemplateheight), sizeof(int4));
   swapByteOrder( (char *)(&hdr->ACtemplatesize), sizeof(int4));
   swapByteOrder( (char *)(&hdr->PCtemplateradius), sizeof(int4));
   swapByteOrder( (char *)(&hdr->PCtemplateheight), sizeof(int4));
   swapByteOrder( (char *)(&hdr->PCtemplatesize), sizeof(int4));
   swapByteOrder( (char *)(&hdr->nangles), sizeof(int4));
}

void swap_model_file_tail(model_file_tail *tail)
{
   swapByteOrder( (char *)(&tail->RPPCmean[0]), sizeof(float4));
   swapByteOrder( (char *)(&tail->RPPCmean[1]), sizeof(float4));
   swapByteOrder( (char *)(&tail->parcomMean), sizeof(float4));
   swapByteOrder( (char *)(&tail->percomMean), sizeof(float4));
   swapByteOrder( (char *)(&tail->RPmean[0]), sizeof(float4));
   swapByteOrder( (char *)(&tail->RPmean[1]), sizeof(float4));
}
