//<license>
/*
   This file is part of Behmeth.
   Author: Behram Kapadia, wiowou@hotmail.com

    Behmeth is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Behmeth is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Behmeth.  If not, see <http://www.gnu.org/licenses/>.
*/
//</license>

#ifndef _MATH_SOLVER_h
#define _MATH_SOLVER_h

#include "dot.h"
#include "sum.h"
#include "vect.h"
#include "matrix.h"

namespace math{

#ifdef MYDEBUG
  class GaussTest;
#endif //MYDEBUG

//! Base class for all other solver classes
template< template <typename> class Storage = Vect, typename T = double >
class Solver 
{
public:

  void BackSolveUp( Matrix< Storage, T >& A, Vect<T>& x, Vect<T>& b ) //tested
  { 
    if ( A.IsTransposed() )
    {
      BackSolveUpTransposed( A, x, b);
      return;
    }
    unsigned long long i = A.Cols();
    x[i-1] = b[i-1] / A[i-1][i-1];
    --i;
    vops::Dot<T> dot;
    for ( ; i > 0; --i )
    {
      x[i-1] = ( b[i-1] - dot( x, A[i-1], i ) ) / A[i-1][i-1];
    }
  }
  
  void BackSolveDown( Matrix< Storage, T >& A, Vect<T>& x, Vect<T>& b ) //tested
  {
    unsigned long long i = 0;
    x[i] = b[i] / A[i][i];
    ++i;
    vops::Dot<T> dot;
    for ( ; i < A.Cols(); ++i )
    {
      x[i] = ( b[i] - dot( i, x, A[i] ) ) / A[i][i];
    }
  }

  void InitialGuess( Matrix< Storage, T >& A, Vect<T>& x )
  {
    vops::Sum<T> sum;
    for ( unsigned long long i = 0; i < A.Cols(); ++i )
    {
      x[i] = sum( A[i] ) / A[i].NNZ(); //use the average for the row
    }
  }
  
  void Zeros( Matrix<Vect,T> &A ) //tested
  {
    for ( unsigned long long i = 0; i < A.Rows(); ++i )
    {
      for ( unsigned long long j = 0; j < A.Cols(); ++j )
      {
        A(i,j) = T();
      }
    }
  }

  void Zeros( Matrix<SparseVect,T> &A ) //tested
  {
    for ( unsigned long long i = 0; i < A.Rows(); ++i )
    {
      A[i].Clear();
    }
  }
  
  void Identity( Matrix<Storage,T> &A ) //tested
  {
    unsigned long long size = A.Cols() < A.Rows() ? A.Cols() : A.Rows();
    Zeros(A);   
    for ( unsigned long long i = 0; i < size; ++i )
    {
      A[i].Set(i, T(1.0) );
    }
  }
  
protected:
  //! Starting with A as a lower triangular matrix with the transposed flag.
  //! Start on the last column of A and work to the first column
  void BackSolveUpTransposed( Matrix< Storage, T >& A, Vect<T>& x, Vect<T>& b ) 
  {
    unsigned long long i = A.Cols();
    x[i-1] = b[i-1] / A[i-1][i-1];
    --i;
    for ( ; i > 0; --i ) // iterating through rows of A transpose, cols of A
    {
      T sum = T();
      for ( unsigned long long j = i; j < A.Cols(); ++j ) //iterating through cols of A transpose, rows of A
      {
        if ( A[j].Pos(0) > i ) // not tested
        {
          break; //means we've reached the end of a band. A banded matrix is always assumed, no holes
        }
        sum += x[j] * A[j][i-1];
      }
      x[i-1] = ( b[i-1] - sum ) / A[i-1][i-1];
    }
  }
};

}/*math*/ 

#endif /*_MATH_SOLVER_h */
