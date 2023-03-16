#define _binomial

#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

/***************************************************************************
*  L O G   G A M M A                                                       *
*                                                                          *
*  Computes the natural log of the gamma function using the Lanczos        *
*  approximation formula.  Gamma is defined by                             *
*                                                                          *
*                                 ( z - 1 )   -t                           *
*         gamma( z ) = Integral[ t           e    dt ]                     *
*                                                                          *
*                                                                          *
*  where the integral ranges from 0 to infinity.  The gamma function       *
*  satisfies                                                               *
*                    gamma( n + 1 ) = n!                                   *
*                                                                          *
*  This algorithm has been adapted from "Numerical Recipes", p. 157.       *
*                                                                          *
***************************************************************************/
double LogGamma( double x )
    {
    static const double 
        coeff0 =  7.61800917300E+1,
        coeff1 = -8.65053203300E+1,
        coeff2 =  2.40140982200E+1,
        coeff3 = -1.23173951600E+0,
        coeff4 =  1.20858003000E-3,
        coeff5 = -5.36382000000E-6,
        stp    =  2.50662827465E+0,
        half   =  5.00000000000E-1,
        fourpf =  4.50000000000E+0,
        one    =  1.00000000000E+0,
        two    =  2.00000000000E+0, 
        three  =  3.00000000000E+0,
        four   =  4.00000000000E+0, 
        five   =  5.00000000000E+0;
    double r = coeff0 / ( x        ) + coeff1 / ( x + one   ) +
               coeff2 / ( x + two  ) + coeff3 / ( x + three ) +
               coeff4 / ( x + four ) + coeff5 / ( x + five  ) ;
    double s = x + fourpf;
    double t = ( x - half ) * log( s ) - s;
    return t + log( stp * ( r + one ) );
    }




/***************************************************************************
*  L O G   F A C T                                                         *
*                                                                          *
*  Returns the natural logarithm of n factorial.  For efficiency, some     *
*  of the values are cached, so they need be computed only once.           *
*                                                                          *
***************************************************************************/
double LogFact( int n )
    {
    static const int Cache_Size = 100;
    static double c[ Cache_Size ] = { 0.0 }; // Cache some of the values.
    if( n <= 1 ) return 0.0;
    if( n < Cache_Size )
        {
        if( c[n] == 0.0 ) c[n] = LogGamma((double)(n+1));
        return c[n];
        }
    return LogGamma((double)(n+1)); // gamma(n+1) == n!
    }


/***************************************************************************
*  B I N O M I A L    C O E F F                                            *
*                                                                          *
*  Compute a given binomial coefficient.  Several rows of Pascal's         *
*  triangle are stored for efficiently computing the small coefficients.   *
*  Higher-order terms are computed using LogFact.                          *
*                                                                          *
***************************************************************************/
double BinomialCoeff( int n, int k )
    {
    double b;
    int    p = n - k;
    if( k <= 1 || p <= 1 )  // Check for errors and special cases.
        {
        if( k == 0 || p == 0 ) return 1;
        if( k == 1 || p == 1 ) return n;
        //cerr << form( "BinomialCoeff(%d,%d) is undefined", n, k );
        return 0;
        }
    static const int  // Store part of Pascal's triange for small coeffs.
        n0[] = { 1 },
        n1[] = { 1, 1 },
        n2[] = { 1, 2, 1 },
        n3[] = { 1, 3, 3, 1 },
        n4[] = { 1, 4, 6, 4, 1 },
        n5[] = { 1, 5, 10, 10, 5, 1 },
        n6[] = { 1, 6, 15, 20, 15, 6, 1 },
        n7[] = { 1, 7, 21, 35, 35, 21, 7, 1 },
        n8[] = { 1, 8, 28, 56, 70, 56, 28, 8, 1 },
        n9[] = { 1, 9, 36, 84, 126, 126, 84, 36, 9, 1 };
    switch( n )
        {
        case 0 : b = n0[k]; break;
        case 1 : b = n1[k]; break;
        case 2 : b = n2[k]; break;
        case 3 : b = n3[k]; break;
        case 4 : b = n4[k]; break;
        case 5 : b = n5[k]; break;
        case 6 : b = n6[k]; break;
        case 7 : b = n7[k]; break;
        case 8 : b = n8[k]; break;
        case 9 : b = n9[k]; break;
        default:
            {
            double x = LogFact( n ) - LogFact( p ) - LogFact( k );
            b = floor( exp( x ) + 0.5 );
            }
        }
    return b;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//******************************************************************************

void k_subset ( int k, double rank, int a[] )

//******************************************************************************
//
//  Purpose:
//
//    Returns the k-subset of a given rank.
//
//  Discussion:
//
//    The routine is given a rank and returns the corresponding subset of K
//    elements of a set of N elements.
//
//    It uses the same ranking that KSUB_NEXT2 uses to generate all the subsets
//    one at a time.
//
//    Note that the value of N itself is not input, nor is it needed.
//
//  Modified:
//
//    13 June 2004
//
//  Reference:
//
//    Albert Nijenhuis and Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int K, the number of elements in the subset.
//
//    Input, int RANK, the rank of the desired subset.
//    There are ( N*(N-1)*...*(N+K-1)) / ( K*(K-1)*...*2*1) such
//    subsets, so RANK must be between 1 and that value.
//
//    Output, int A[K], K distinct integers in order between
//    1 and N, which define the subset.
//
{
  int i;
  //int ii;
  int ip;
  double iprod;
  double jrank;
//
  jrank = rank - 1;
  for ( i = k; 1 <= i; i-- )
  {
    ip = i - 1;
    iprod = 1;
    for ( ; ; )
    {
      ip = ip + 1;

      if ( ip != i )
      {
        iprod = ( ip * iprod ) / ( ip - i );
      }

      if ( jrank < iprod )
      {
        break;
      }

    }
    if ( ip != i )
    {
      iprod = ( ( ip - i ) * iprod ) / ip;
    }

    jrank = jrank - iprod;
    a[i-1] = ip;

  }

  return;
}

