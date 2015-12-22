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

#ifndef _MATH_VOPS_SUM_h
#define _MATH_VOPS_SUM_h

#include "vect.h"
#include "sparseVect.h"
#include "typedefs.h"

namespace math{
namespace vops{

#ifdef MYDEBUG
  class SumTest;
#endif //MYDEBUG

//! Reduction sum
//!     Operations occur up to but not including parameter "end"
//!     Operations begin at parameter "start"
template<class T = double>
class Sum
{

public:
  T operator()( SparseVect<T> &X, ULong start = 0 ) //tested
  {
    T res = T();
    for ( ULong i = start; i < X.NNZ(); ++i )
    {
      res += X.m_data[i];
    }
    return res;
  }

  T operator()( ULong end, SparseVect<T> &X )  //tested
  {
    T res = T();
    ULong i = 0;
    ULong iend = X.NNZ();
    if ( X.m_pos[iend - 1] >= end )
    {
      iend = X.UpperBoundIdx(end);
    }

    for ( ULong i = 0; i < iend; ++i )
    {
      res += X.m_data[i];
    }
    return res;
  }
  
  T operator()( Vect<T> &X, ULong start = 0  )  //tested
  {
    T res = T(); //make sure this operation works for T
    for ( ULong i = start; i < X.Size(); ++i )
    {
      res += X.m_data[i];
    }
    return res;
  }

  T operator()( ULong end, Vect<T> &X ) //tested
  {
    T res = T(); //make sure this operation works for T
    for ( ULong i = 0; i < end; ++i )
    {
      res += X.m_data[i];
    }
    return res;
  }
protected:

private:


#ifdef MYDEBUG
  friend class SumTest;
#endif //MYDEBUG
};

}/*vops*/ }/*math*/ 

#endif /*_MATH_VOPS_SUM_h */
