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

#ifndef _MATH_VOPS_ADD_h
#define _MATH_VOPS_ADD_h

#include "vect.h"
#include "sparseVect.h"


namespace math{
namespace vops{

#ifdef MYDEBUG
  class AddTest;
#endif //MYDEBUG

//!Z = X + Y
//!     Operations occur up to but not including parameter "end"
//!     Operations begin at parameter "start"
template<class T = double>
class Add
{

public:
  void operator()( SparseVect<T> &X, SparseVect<T> &Y, SparseVect<T> &Z, unsigned long long start = 0 ) //tested
  {
    unsigned long long i = 0;
    unsigned long long j = 0;
    if ( start != 0 )
    {
      i = X.UpperBoundIdx(start);
      j = Y.UpperBoundIdx(start);
    }
    while ( i < X.NNZ() && j < Y.NNZ() )
    {
      if ( X.m_pos[i] == Y.m_pos[j] )
      {
        Z.PushBack( X.m_pos[i], X.m_data[i] + Y.m_data[j] );
        ++i; ++j;
      }
      else if ( X.m_pos[i] < Y.m_pos[j] )
      {
        Z.PushBack( X.m_pos[i], X.m_data[i] );
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
      Z.PushBack( X.m_pos[i], X.m_data[i] );
    }
    for ( ; j < Y.NNZ(); ++j )
    {
      Z.PushBack( Y.m_pos[j], Y.m_data[j] );
    }
  }

  void operator()( unsigned long long end, SparseVect<T> &X, SparseVect<T> &Y, SparseVect<T> &Z ) //tested
  {
    unsigned long long i = 0;
    unsigned long long j = 0;
    unsigned long long iend = X.NNZ();
    unsigned long long jend = Y.NNZ();
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
        Z.PushBack( X.m_pos[i], X.m_data[i] + Y.m_data[j] );
        ++i; ++j;
      }
      else if ( X.m_pos[i] < Y.m_pos[j] )
      {
        Z.PushBack( X.m_pos[i], X.m_data[i] );
        ++i;
      }
      else
      {
        Z.PushBack( Y.m_pos[j], Y.m_data[j] );
        ++j;
      }
    }
    while ( i < iend && i < X.NNZ() )
    {
      Z.PushBack( X.m_pos[i], X.m_data[i] );
      ++i;
    }
    while ( j < jend  && j < Y.NNZ() )
    {
      Z.PushBack( Y.m_pos[j], Y.m_data[j] );
      ++j;
    }
  }
  
  void operator()( Vect<T> &X, SparseVect<T> &Y, Vect<T> &Z, unsigned long long start = 0  ) //tested
  {
    unsigned long long j = 0;
    if ( start != 0 )
    {
      j = Y.UpperBoundIdx(start);
    }
    for( unsigned long long i = start; i < Z.Size(); ++i )
    {
      if ( i == Y.m_pos[j] )
      {
        Z[i] = X.m_data[i] + Y.m_data[j];
        if ( j < Y.NNZ() )
        {
          ++j;
        }
      }
      else
      {
        Z[i] = X[i];
      }
    }
  }
  
  void operator()( SparseVect<T> &X, Vect<T> &Y, Vect<T> &Z, unsigned long long start = 0  ) //tested
  {
    operator()( Y, X, Z, start );
  }
  
  void operator()( unsigned long long end, Vect<T> &X, SparseVect<T> &Y, Vect<T> &Z ) //tested
  {
    //end = Z.Size() < end ? Z.Size() : end;
    unsigned long long j = 0;
    unsigned long long jend = Y.NNZ();
    if ( Y.m_pos[jend - 1] >= end )
    {
      jend = Y.UpperBoundIdx(end);
    }
    for( unsigned long long i = 0; i < end; ++i )
    {
      if ( i == Y.m_pos[j] )
      {
        Z[i] = X.m_data[i] + Y.m_data[j];
        if ( j < Y.NNZ() )
        {
          ++j;
        }
      }
      else
      {
        Z[i] = X[i];
      }
    }
  }

  void operator()( unsigned long long end, SparseVect<T> &X, Vect<T> &Y, Vect<T> &Z  ) //tested
  {
    operator()( end, Y, X, Z );
  }
  
  void operator()( Vect<T> &X, Vect<T> &Y, Vect<T> &Z, unsigned long long start = 0  ) //tested
  {
    for ( unsigned long long i = start; i < X.Size(); ++i )
    {
      Z[i] = X[i] + Y[i];
    }
  }

  void operator()( unsigned long long end, Vect<T> &X, Vect<T> &Y, Vect<T> &Z  ) //tested
  {
    //end = Z.Size() < end ? Z.Size() : end;
    for ( unsigned long long i = 0; i < end; ++i )
    {
      Z[i] = X[i] + Y[i];
    }
  }
protected:

private:


#ifdef MYDEBUG
  friend class AddTest;
#endif //MYDEBUG
};

}/*vops*/ }/*math*/ 

#endif /*_MATH_VOPS_ADD_h */
