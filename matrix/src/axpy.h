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

#ifndef _MATH_VOPS_AXPY_h
#define _MATH_VOPS_AXPY_h

#include "vect.h"
#include "sparseVect.h"
#include "typedefs.h"

namespace math{
namespace vops{

#ifdef MYDEBUG
  class AxpyTest;
#endif //MYDEBUG

//! Y = A*X+Y 
//! Z = A*X+Y
//!     Operations occur up to but not including parameter "end"
//!     Operations begin at parameter "start"
template<class T = double>
class Axpy
{

public:
  void operator()( T A, SparseVect<T> &X, SparseVect<T> &Y, ULong start = 0 ) //tested
  {
    ULong i = 0;
    ULong j = 0;
    if ( start != 0 )
    {
      i = X.UpperBoundIdx(start);
      j = Y.UpperBoundIdx(start);
    }
    while ( i < X.NNZ() && j < Y.NNZ() )
    {
      if ( X.m_pos[i] == Y.m_pos[j] )
      {
        Y.m_data[j] += A * X.m_data[i];
        ++i; ++j;
      }
      else if ( X.m_pos[i] < Y.m_pos[j] )
      {
        Y.Set( X.m_pos[i], A * X.m_data[i] );
        ++j;
        ++i;
      }
      else
      {
        ++j;
      }
    }
    for ( ; i < X.NNZ(); ++i )
    {
      Y.PushBack( X.m_pos[i], A * X.m_data[i] );
    }
  }

  void operator()( ULong end, T A, SparseVect<T> &X, SparseVect<T> &Y ) //tested
  {
    ULong i = 0;
    ULong j = 0;
    ULong iend = X.NNZ();
    ULong jend = Y.NNZ();
    if ( X.m_pos[iend - 1] >= end )
    {
      iend = X.UpperBoundIdx(end);
    }
    if ( Y.m_pos[jend - 1] >= end )
    {
      jend = Y.UpperBoundIdx(end);
    }
    
    while ( i < iend && j < jend )
    {
      if ( X.m_pos[i] == Y.m_pos[j] )
      {
        Y.m_data[j] += A * X.m_data[i];
        ++i; ++j;
      }
      else if ( X.m_pos[i] < Y.m_pos[j] )
      {
        Y.Set( X.m_pos[i], A * X.m_data[i] );
        ++jend;
        ++j;
        ++i;
      }
      else
      {
        ++j;
      }
    }
    for ( ; i < iend; ++i )
    {
      Y.PushBack( X.m_pos[i], A * X.m_data[i] );
    }
  }
  
  void operator()( T A, SparseVect<T> &X, SparseVect<T> &Y, SparseVect<T> &Z, ULong start = 0 ) //tested
  {
    ULong i = 0;
    ULong j = 0;
    if ( start != 0 )
    {
      i = X.UpperBoundIdx(start);
      j = Y.UpperBoundIdx(start);
    }
    while ( i < X.NNZ() && j < Y.NNZ() )
    {
      if ( X.m_pos[i] == Y.m_pos[j] )
      {
        Z.PushBack( X.m_pos[i], A * X.m_data[i] + Y.m_data[j] );
        ++i; ++j;
      }
      else if ( X.m_pos[i] < Y.m_pos[j] )
      {
        Z.PushBack( X.m_pos[i], A * X.m_data[i] );
        ++i;
      }
      else
      {
        Z.PushBack( Y.m_pos[j], Y.m_data[j] );
        ++j;
      }
    }
    for ( ; i < X.NNZ(); ++i )
    {
      Z.PushBack( X.m_pos[i], A * X.m_data[i] );
    }
    for ( ; j < Y.NNZ(); ++j )
    {
      Z.PushBack( Y.m_pos[j], Y.m_data[j] );
    }
  }

  void operator()( ULong end, T A, SparseVect<T> &X, SparseVect<T> &Y, SparseVect<T> &Z ) //tested
  {
    ULong i = 0;
    ULong j = 0;
    ULong iend = X.NNZ();
    ULong jend = Y.NNZ();
    if ( X.m_pos[iend - 1] >= end )
    {
      iend = X.UpperBoundIdx(end);
    }
    if ( Y.m_pos[jend - 1] >= end )
    {
      jend = Y.UpperBoundIdx(end);
    }
    
    while ( i < iend && j < jend )
    {
      if ( X.m_pos[i] == Y.m_pos[j] )
      {
        Z.PushBack( X.m_pos[i], A * X.m_data[i] + Y.m_data[j] );
        ++i; ++j;
      }
      else if ( X.m_pos[i] < Y.m_pos[j] )
      {
        Z.PushBack( X.m_pos[i], A * X.m_data[i] );
        ++i;
      }
      else
      {
        Z.PushBack( Y.m_pos[j], Y.m_data[j] );
        ++j;
      }
    }
    for ( ; i < iend; ++i )
    {
      Z.PushBack( X.m_pos[i], A * X.m_data[i] );
    }
    for ( ; j < jend; ++j )
    {
      Z.PushBack( Y.m_pos[j], Y.m_data[j] );
    }
  }
  
  void operator()( T A, SparseVect<T> &X, Vect<T> &Y, ULong start = 0 ) //tested
  {
    ULong i = 0;
    if ( start != 0 )
    {
      i = X.UpperBoundIdx(start);
    }
    while ( i < X.NNZ() && X.m_pos[i] < Y.Size() )
    {
      Y.m_data[ X.m_pos[i] ] += A * X.m_data[i];
      ++i;
    }
  }

  void operator()( T A, Vect<T> &X, SparseVect<T> &Y, ULong start = 0 ) //tested
  {
    operator()( A, X, Y, start );
  }  
  
  void operator()( T A, Vect<T> &X, Vect<T> &Y, ULong start = 0 )
  {
    for ( ULong i = start; i < Y.Size(); ++i )
    {
      Y.m_data[i] += A * X.m_data[i];
    }
  }

  void operator()( ULong end, T A, Vect<T> &X, Vect<T> &Y )
  {
    //end = Y.Size() < end ? Y.Size() : end;
    for ( ULong i = 0; i < end; ++i )
    {
      Y.m_data[i] += A * X.m_data[i];
    }
  }
protected:

private:


#ifdef MYDEBUG
  friend class AxpyTest;
#endif //MYDEBUG
};

}/*vops*/ }/*math*/ 

#endif /*_MATH_VOPS_AXPY_h */
