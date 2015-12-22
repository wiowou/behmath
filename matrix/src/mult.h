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

#ifndef _MATH_VOPS_MULT_h
#define _MATH_VOPS_MULT_h

#include "vect.h"
#include "sparseVect.h"
#include "typedefs.h"

namespace math{
namespace vops{

#ifdef MYDEBUG
  class MultTest;
#endif //MYDEBUG

//!Z = X * Y
//!     Operations occur up to but not including parameter "end"
//!     Operations begin at parameter "start"
//!     Z != X, Z != Y for SparseVect's
template<class T = double>
class Mult
{

public:
  //! Z = X * Y * W
  template< template <typename> class Vec1, template <typename> class Vec2, template <typename> class Vec3, template <typename> class Vec4 >
  void operator()( Vec1<T> &W, Vec2<T> &X, Vec3<T> &Y, Vec4<T> &Z, ULong start = 0  )
  {
    operator()( W, X, Z, start );
    operator()( Z, Y, start );
  }

  //! Z = X * Y * W
  template< template <typename> class Vec1, template <typename> class Vec2, template <typename> class Vec3, template <typename> class Vec4 >
  void operator()( ULong end, Vec1<T> &W, Vec2<T> &X, Vec3<T> &Y, Vec4<T> &Z )
  {
    operator()( end, W, X, Z );
    operator()( end, Z, Y );
  }  
  
  //! X = X * Y
  void operator()( SparseVect<T> &X, SparseVect<T> &Y, ULong start = 0  ) //tested
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
        X.m_data[i] *= Y.m_data[j];
        ++i; ++j;
      }
      else if ( X.m_pos[i] < Y.m_pos[j] )
      {
        // resizing X so that it still consists of actual non-zero entries
        --X.m_nnz;
        for ( ULong k = i; k < X.NNZ(); ++k )
        {
          X.m_pos[k] = X.m_pos[k+1];
          X.m_data[k] = X.m_data[k+1];
        }
      }
      else
      {
        ++j;
      }
    } 
  }
  
  //! X = X * Y
  void operator()( ULong end, SparseVect<T> &X, SparseVect<T> &Y ) //tested
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
        X.m_data[i] *= Y.m_data[j];
        ++i; ++j;
      }
      else if ( X.m_pos[i] < Y.m_pos[j] )
      {
        // resizing X so that it still consists of actual non-zero entries
        --X.m_nnz;
        for ( ULong k = i; k < X.NNZ(); ++k )
        {
          X.m_pos[k] = X.m_pos[k+1];
          X.m_data[k] = X.m_data[k+1];
        }
        --iend;
      }
      else
      {
        ++j;
      }
    } 
  }  
  
  //! X = X * Y
  void operator()( SparseVect<T> &X, Vect<T> &Y, ULong start = 0  ) //tested
  {
    ULong j = 0;
    if ( start != 0 )
    {
      j = X.UpperBoundIdx(start);
    }
    while ( j < X.NNZ() && X.m_pos[j] < Y.Size() )
    {
      X.m_data[j] *= Y[ X.m_pos[j] ];
      ++j;
    }
  }

  //! X = X * Y
  void operator()( ULong end, SparseVect<T> &X, Vect<T> &Y ) //tested
  {
    ULong j = 0;
    ULong jend = X.NNZ();
    if ( X.m_pos[jend - 1] >= end )
    {
      jend = X.UpperBoundIdx(end);
    }
    while ( j < jend && X.m_pos[j] < Y.Size() )
    {
      X.m_data[j] *= Y[ X.m_pos[j] ];
      ++j;
    }
  }

  void operator()( Vect<T> &X, SparseVect<T> &Y, ULong start = 0 ) //tested
  {
    ULong j = 0;
    if ( start != 0 )
    {
      j = Y.UpperBoundIdx(start);
    }
    for ( ULong i = start; i < X.Size(); ++i )
    {
      if ( i == Y.m_pos[j] )
      {
        X[i] *= Y.m_data[j++];
      }
      else
      {
        X[i] = T();
      }
    }
  }

  void operator()( ULong end, Vect<T> &X, SparseVect<T> &Y ) //tested
  {
    ULong j = 0;
    //end = Z.Size() < end ? Z.Size() : end;
    for ( ULong i = 0; i < end; ++i )
    {
      if ( i == Y.m_pos[j] )
      {
        X[i] *= Y.m_data[j++];
      }
      else
      {
        X[i] = T();
      }
    }
  }
  
  void operator()( Vect<T> &X, Vect<T> &Y, ULong start = 0  )
  {
    for ( ULong i = start; i < X.Size(); ++i )
    {
      X[i] *= Y[i];
    }
  }

  void operator()( ULong end, Vect<T> &X, Vect<T> &Y  )
  {
    //end = Z.Size() < end ? Z.Size() : end;
    for ( ULong i = 0; i < end; ++i )
    {
      X[i] *= Y[i];
    }
  }
  
  void operator()( SparseVect<T> &X, SparseVect<T> &Y, SparseVect<T> &Z, ULong start = 0  ) //tested
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
        Z.PushBack( X.m_pos[i], X.m_data[i] * Y.m_data[j] );
        ++i; ++j;
      }
      else if ( X.m_pos[i] < Y.m_pos[j] )
      {
        ++i;
      }
      else
      {
        ++j;
      }
    } 
  }

  void operator()( ULong end, SparseVect<T> &X, SparseVect<T> &Y, SparseVect<T> &Z  ) //tested
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
        Z.PushBack( X.m_pos[i], X.m_data[i] * Y.m_data[j] );
        ++i; ++j;
      }
      else if ( X.m_pos[i] < Y.m_pos[j] )
      {
        ++i;
      }
      else
      {
        ++j;
      }
    } 
  }
  
  void operator()( SparseVect<T> &X, SparseVect<T> &Y, Vect<T> &Z, ULong start = 0  ) //tested
  {
    ULong i = 0;
    ULong j = 0;
    if ( start != 0 )
    {
      i = X.UpperBoundIdx(start);
      j = Y.UpperBoundIdx(start);
    }
    
    ULong startZeros;
    if ( Y.m_pos[j] > X.m_pos[i] )
    {
      startZeros = Y.m_pos[j] < Z.Size() ? Y.m_pos[j] : Z.Size();
      i = X.UpperBoundIdx( startZeros );
    }
    else
    {
      startZeros = X.m_pos[i] < Z.Size() ? X.m_pos[i] : Z.Size();
      j = Y.UpperBoundIdx( startZeros );
    }

    ULong k;
    for ( k = start; k < startZeros; ++k )
    {
      Z[k] = T();
    }
    for ( ; k < Z.Size(); ++k )
    {
      if ( X.m_pos[i] == Y.m_pos[j] && X.m_pos[i] == k )
      {
        Z[k] = X.m_data[i++] * Y.m_data[j++];
      }
      else if ( X.m_pos[i] < Y.m_pos[j] )
      {
        Z[k] = T();
        ++i;
      }
      else if ( X.m_pos[i] > Y.m_pos[j] )
      {
        Z[k] = T();
        ++j;
      }
      else
      {
        Z[k] = T();
      }
    }
  }

  void operator()( ULong end, SparseVect<T> &X, SparseVect<T> &Y, Vect<T> &Z ) //tested
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
    
    ULong startZeros;
    if ( Y.m_pos[j] > X.m_pos[i] )
    {
      startZeros = Y.m_pos[j] < Z.Size() ? Y.m_pos[j] : end;
      i = X.UpperBoundIdx( startZeros );
    }
    else
    {
      startZeros = X.m_pos[i] < Z.Size() ? X.m_pos[i] : end;
      j = Y.UpperBoundIdx( startZeros );
    }

    ULong k;
    for ( k = 0; k < startZeros; ++k )
    {
      Z[k] = T();
    }
    for ( ; k < end; ++k )
    {
      if ( X.m_pos[i] == Y.m_pos[j] && X.m_pos[i] == k )
      {
        Z[k] = X.m_data[i++] * Y.m_data[j++];
      }
      else if ( X.m_pos[i] < Y.m_pos[j] )
      {
        Z[k] = T();
        ++i;
      }
      else if ( X.m_pos[i] > Y.m_pos[j] )
      {
        Z[k] = T();
        ++j;
      }
      else
      {
        Z[k] = T();
      }
    }
  }
  
  void operator()( Vect<T> &X, SparseVect<T> &Y, SparseVect<T> &Z, ULong start = 0  ) //tested
  {
    ULong j = 0;
    if ( start != 0 )
    {
      j = Y.UpperBoundIdx(start);
    }
    while ( j < Y.NNZ() && Y.m_pos[j] < X.Size() )
    {
      Z.PushBack( Y.m_pos[j], Y.m_data[j] * X.m_data[ Y.m_pos[j] ] );
      ++j;
    }
  }
  
  void operator()( SparseVect<T> &X, Vect<T> &Y, SparseVect<T> &Z, ULong start = 0  ) //tested
  {
    operator()( Y, X, Z, start );
  }

  void operator()( ULong end, Vect<T> &X, SparseVect<T> &Y, SparseVect<T> &Z ) //tested
  {
    ULong j = 0;
    ULong jend = Y.NNZ();
    if ( Y.m_pos[jend - 1] >= end )
    {
      jend = Y.UpperBoundIdx(end);
    }
    while ( j < jend && Y.m_pos[j] < X.Size() )
    {
      Z.PushBack( Y.m_pos[j], Y.m_data[j] * X.m_data[ Y.m_pos[j] ] );
      ++j;
    }
  }
  
  void operator()( ULong end, SparseVect<T> &X, Vect<T> &Y, SparseVect<T> &Z  ) //tested
  {
    operator()( end, Y, X, Z );
  }
  
  void operator()( Vect<T> &X, SparseVect<T> &Y, Vect<T> &Z, ULong start = 0 )  //tested
  {
    ULong j = 0;
    if ( start != 0 )
    {
      j = Y.UpperBoundIdx(start);
    }
    for ( ULong i = start; i < Z.Size(); ++i )
    {
      if ( i == Y.m_pos[j] )
      {
        Z[i] = X[i] * Y.m_data[j++];
      }
      else
      {
        Z[i] = T();
      }
    }
  }
  
  void operator()( SparseVect<T> &X, Vect<T> &Y, Vect<T> &Z, ULong start = 0  )  //tested
  {
    operator()( Y, X, Z, start );
  }

  void operator()( ULong end, Vect<T> &X, SparseVect<T> &Y, Vect<T> &Z )  //tested
  {
    ULong j = 0;
    //end = Z.Size() < end ? Z.Size() : end;
    for ( ULong i = 0; i < end; ++i )
    {
      if ( i == Y.m_pos[j] )
      {
        Z[i] = X[i] * Y.m_data[j++];
      }
      else
      {
        Z[i] = T();
      }
    }
  }

  void operator()( ULong end, SparseVect<T> &X, Vect<T> &Y, Vect<T> &Z )  //tested
  {
    operator()( end, Y, X, Z );
  }
  
  void operator()( Vect<T> &X, Vect<T> &Y, Vect<T> &Z, ULong start = 0  )
  {
    for ( ULong i = start; i < X.Size(); ++i )
    {
      Z[i] = X[i] * Y[i];
    }
  }

  void operator()( ULong end, Vect<T> &X, Vect<T> &Y, Vect<T> &Z  )
  {
    //end = Z.Size() < end ? Z.Size() : end;
    for ( ULong i = 0; i < end; ++i )
    {
      Z[i] = X[i] * Y[i];
    }
  }
  
protected:

private:


#ifdef MYDEBUG
  friend class MultTest;
#endif //MYDEBUG
};

}/*vops*/ }/*math*/ 

#endif /*_MATH_VOPS_MULT_h */
