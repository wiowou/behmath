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

#ifndef _MATH_VOPS_DOT_h
#define _MATH_VOPS_DOT_h

#include "sparseVect.h"
#include "vect.h"
#include "mult.h"
#include "sum.h"
#include "matrix.h"


namespace math{
namespace vops{

#ifdef MYDEBUG
  class DotTest;
#endif //MYDEBUG

//! Dot product of two vectors: Z = X dot Y
//!     Operations occur up to but not including parameter "end"
//!     Operations begin at parameter "start"
template<class T = double>
class Dot
{

public:
  //! This is a matrix * diag matrix * matrix, regular matrix multiplication. The rows are dotted, so the algorithm
  //! ignores whether the m_isTransposed flag is set. No resizing occurs.
  template< template <typename> class Vec, template <typename> class Vec2, template <typename> class Vec3, template <typename> class Vec4 >
  void operator()( Matrix< Vec, T > &lhs, Vec2<T>& diag, Matrix< Vec3, T > &rhs, Matrix< Vec4, T > &result ) 
  {
    for ( unsigned long long i = 0; i < result.PRows(); ++i )
    {
      for ( unsigned long long j = 0; j < result.PCols(); ++j )
      {
        result[i][j] = operator()( lhs[i], diag, rhs[j] );
      }
    }
  }

  //! This is a matrix * diag matrix * matrix, regular matrix multiplication, but lhs is multiplied with itself. 
  //! Now, however, there's a diagonal matrix inbetween the lhs and itself.
  //! The rows are dotted, so the algorithm ignores whether the m_isTransposed flag is set. No resizing occurs.
  template< template <typename> class Vec, template <typename> class Vec2, template <typename> class Vec3 >
  void operator()( Matrix< Vec, T > &lhs, Vec2<T>& diag, Matrix< Vec3, T > &result ) 
  {
    for ( unsigned long long i = 0; i < result.PRows(); ++i )
    {
      for ( unsigned long long j = 0; j < result.PCols(); ++j )
      {
        result[i][j] = operator()( lhs[i], diag, lhs[j] );
      }
    }
  }

  //! This is a matrix * matrix, regular matrix multiplication. The rows are dotted, so the algorithm
  //! ignores whether the m_isTransposed flag is set. No resizing occurs.
  template< template <typename> class Vec, template <typename> class Vec2, template <typename> class Vec3 >
  void operator()( Matrix< Vec, T > &lhs, Matrix< Vec2, T > &rhs, Matrix< Vec3, T > &result )
  {
    for ( unsigned long long i = 0; i < result.PRows(); ++i )
    {
      for ( unsigned long long j = 0; j < result.PCols(); ++j )
      {
        result[i][j] = operator()( lhs[i], rhs[j] );
      }
    }
  }

  //! This is a matrix * matrix, regular matrix multiplication, but lhs is multiplied with itself. 
  //! The rows are dotted, so the algorithm ignores whether the m_isTransposed flag is set. No resizing occurs.
  template< template <typename> class Vec, template <typename> class Vec2 >
  void operator()( Matrix< Vec, T > &lhs, Matrix< Vec2, T > &result )
  {
    for ( unsigned long long i = 0; i < result.PRows(); ++i )
    {
      for ( unsigned long long j = 0; j < result.PCols(); ++j )
      {
        result[i][j] = operator()( lhs[i], lhs[j] );
      }
    }
  }
  
  //! This is a matrix * vector. The rows are dotted, so the algorithm
  //! ignores whether the m_isTransposed flag is set. No resizing occurs.  
  template< template <typename> class Vec, template <typename> class Vec2 >
  void operator()( Matrix< Vec, T > &matrix, Vec2<T>& vector, Vect<T>& result ) //tested
  {
    for ( unsigned long long i = 0; i < matrix.PRows(); ++i )
    {
      result.Set( i, operator()( matrix[i], vector ) );
    }
  }

  template< class Vec1, class Vec2, class Vec3 >
  T operator()( Vec1 &W, Vec2 &X, Vec3 &Y, unsigned long long start = 0 )
  {
    Mult<T> mult;
    SparseVect<T> Z;
    mult( W, X, Y, Z, start );
    Sum<T> sum;
    return sum( Z );
  }
  
  template< class Vec1, class Vec2, class Vec3 >
  T operator()( unsigned long long end, Vec1 &W, Vec2 &X, Vec3 &Y )
  {
    Mult<T> mult;
    SparseVect<T> Z;
    mult( end, W, X, Y, Z );
    Sum<T> sum;
    return sum( Z );
  }

  T operator()( Vect<T> &W, Vect<T> &X, Vect<T> &Y, unsigned long long start = 0 )
  {
    Mult<T> mult;
    Vect<T> Z;
    Z.Resize( X.Size() );
    mult( W, X, Y, Z, start );
    Sum<T> sum;
    return sum( Z, start );
  }

  T operator()( unsigned long long end, Vect<T> &W, Vect<T> &X, Vect<T> &Y )
  {
    Mult<T> mult;
    Vect<T> Z;
    Z.Resize( X.Size() );
    mult( end, W, X, Y, Z );
    Sum<T> sum;
    return sum( end, Z );
  }
  
  template< class Vec1, class Vec2 >
  T operator()( Vec1 &X, Vec2 &Y, unsigned long long start = 0 ) //tested
  {
    Mult<T> mult;
    SparseVect<T> Z;
    mult( X, Y, Z, start );
    Sum<T> sum;
    return sum( Z );
  }

  template< class Vec1, class Vec2 >
  T operator()( unsigned long long end, Vec1 &X, Vec2 &Y )  //tested
  {
    Mult<T> mult;
    SparseVect<T> Z;
    mult( end, X, Y, Z );
    Sum<T> sum;
    return sum( end, Z );
  }
  
  T operator()( Vect<T> &X, Vect<T> &Y, unsigned long long start = 0 ) //tested
  {
    Mult<T> mult;
    Vect<T> Z;
    Z.Resize( X.Size() );
    mult( X, Y, Z, start );
    Sum<T> sum;
    return sum( Z, start );
  }

  T operator()( unsigned long long end, Vect<T> &X, Vect<T> &Y ) //tested
  {
    Mult<T> mult;
    Vect<T> Z;
    Z.Resize( X.Size() );
    mult( end, X, Y, Z );
    Sum<T> sum;
    return sum( end, Z );
  }
  
protected:

private:


#ifdef MYDEBUG
  friend class DotTest;
#endif //MYDEBUG
};

}/*vops*/ }/*math*/ 

#endif /*_MATH_VOPS_DOT_h */
