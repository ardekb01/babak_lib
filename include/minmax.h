///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// Copyright (C) 2024 Babak A. Ardekani, PhD - All Rights Reserved.  //
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

#ifndef _minmax_h
#include <cstddef>

//Input: 'array' of size 'n'
//Ouputs: The maximum values of 'array' is returned in 'max'
template<class TYPE> void arraymax(TYPE *array, int n, TYPE &max)
{
   if(array==NULL) return;

   max=array[0];

   for(int i=0; i<n; i++)
   {
      if(array[i]>max) 
      {
         max=array[i];
      }
   }

   return;
}

//Input: 'array' of size 'n'
//Ouputs: The minimum and maximum values of 'array' are returned 
//in 'min' and 'max' variables.
template<class TYPE> void minmax(TYPE *array, int n, TYPE &min, TYPE &max)
{
   if(array==NULL) return;

   min=max=array[0];

   for(int i=0; i<n; i++)
   {
      if(array[i]<min) 
      {
         min=array[i];
      }
      else if(array[i]>max) 
      {
         max=array[i];
      }
   }

   return;
}

// This version returns the minimum and maximum values within the masked region specified by 'msk'
// Returns the number of voxels in 'msk'
template<class TYPE1, class TYPE2> int minmax(TYPE1 *array, TYPE2 *msk, int n, TYPE1 &min, TYPE1 &max)
{
   int nbv=0;

   if(array==NULL || msk==NULL ) return(nbv);

   for(int i=0; i<n; i++)
   {
      if( msk[i]>0.0)
      {
         min=max=array[i];
         break;
      }
   }

   for(int i=0; i<n; i++)
   {
      if( msk[i]>0.0)
      {
         nbv++;

         if(array[i]<min) 
         {
            min=array[i];
         }

         if(array[i]>max) 
         {
            max=array[i];
         }
      }
   }

   return(nbv);
}

#define _minmax_h
#endif
