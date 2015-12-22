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

#ifndef _MATH_SOLVER_ITERATIVE_h
#define _MATH_SOLVER_ITERATIVE_h

#include "typedefs.h"
#include "impl/matrixConfig.h"
#include "matrix.h"
#include "sum.h"
#include <cmath>

namespace math{
namespace solver{

#ifdef MYDEBUG
  class IterativeTest;
#endif //MYDEBUG

template< template <typename> class Storage = SparseVect, typename T = double >
class Iterative
{

public:
  Iterative()
  {
    Init();
  }
  
  void Clear()
  {
    Init();
  }
  
  void A( Matrix< Storage, T >& matrix )
  {
    m_A = &matrix;
  }
  
  void X( Vect<T>& v )
  {
    m_x = &v;
    m_xNew.Resize( m_x->Size() );
    m_notWithinTol.Resize( m_x->Size() );
  }
  
  void B( Vect<T>& v )
  {
    m_b = &v;
  }
  
  T Tolerance()
  {
    return m_tolerance;
  }
  
  void Tolerance( T tol )
  {
    m_tolerance = tol;
  }

  void MaxIter( ULong imax )
  {
    m_maxIter = imax;
  }
  
  ULong MaxIter()
  {
    return m_maxIter;
  }
  
  ULong NIter()
  {
    return m_nIter;
  }
  
protected:
  void InitialGuess()
  {
    Matrix< Storage, T >& A = *m_A;
    Vect<T>& x = *m_x;
    vops::Sum<T> sum;
    for ( ULong i = 0; i < A.Rows(); ++i )
    {
      x[i] = sum( A[i] ) / A[i].NNZ(); //use the average for the row
    }
  }
  
  inline bool Continue()
  {
    //this can be parallelized, it's a simple redux sum
    for ( ULong i = 0; i < m_xNew.Size(); ++i )
    {
      if ( m_notWithinTol[i] )
      {
        return true;
      }
    }
    return false;
  }
  
  void Init()
  {
    m_A = NULL;
    m_b = NULL;
    m_x = NULL;
    m_tolerance = 1e-5;
    m_maxIter = 0;
    m_nIter = 0;
  }
  
protected:
  Matrix< Storage, T >* m_A;
  Vect<T>* m_b;
  Vect<T>* m_x;
  Vect<bool> m_notWithinTol;
  Vect<T> m_xNew;
  T m_tolerance;
  ULong m_maxIter;
  ULong m_nIter;

private:


#ifdef MYDEBUG
  friend class IterativeTest;
#endif //MYDEBUG
};

}/*solver*/ }/*math*/ 

#endif /*_MATH_SOLVER_ITERATIVE_h */
