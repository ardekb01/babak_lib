#include <stdlib.h>
#include <stdio.h>
#include "babak_lib.h"
#include "swap.h"

// x[0] -> byte-reversed
// x[1] -> byte-reversed
// ...
bool swap_int_array(int *x, size_t n)
{
   if (x == nullptr || n == 0)
   {
      return false;
   }

   for (size_t i = 0; i < n; i++)
   {
      if (!swapByteOrder((char *)(x + i), sizeof(int)))
      {
         return false;
      }
   }

   return true;
}

bool swap_double_array(double *x, size_t n)
{
   if (x == nullptr || n == 0)
   {
      return false;
   }

   for (size_t i = 0; i < n; i++)
   {
      if (!swapByteOrder((char *)(x + i), sizeof(double)))
      {
         return false;
      }
   }

   return true;
}

bool swap_float_array(float *x, size_t n)
{
   if (x == nullptr || n == 0)
   {
      return false;
   }

   for (size_t i = 0; i < n; i++)
   {
      if (!swapByteOrder((char *)(x + i), sizeof(float)))
      {
         return false;
      }
   }

   return true;
}

// AB CD EF GH -> BA DC FE HG
bool swapN(char *in, size_t N)
{
   char dum[2];

   // N should be even and positive.
   if (N == 0 || N % 2 != 0)
   {
      return false;
   }

   if (in == nullptr)
   {
      return false;
   }

   for (size_t i = 0; i < N / 2; i++)
   {
      dum[0] = in[i * 2 + 1];
      dum[1] = in[i * 2];

      in[i * 2] = dum[0];
      in[i * 2 + 1] = dum[1];
   }

   return true;
}

// [A B C D] -> [D C B A]
bool swapByteOrder(char *in, size_t N)
{
   char *dum;

   if (in == nullptr || N == 0)
   {
      return false;
   }

   dum = (char *)malloc(N);

   if (dum == nullptr)
   {
      return false;
   }

   for (size_t i = 0; i < N; i++)
   {
      dum[i] = in[N - 1 - i];
   }

   for (size_t i = 0; i < N; i++)
   {
      in[i] = dum[i];
   }

   free(dum);

   return true;
}

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
