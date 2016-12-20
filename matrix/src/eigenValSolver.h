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

#ifndef _MATH_SOLVER_EIGENVALSOLVER_h
#define _MATH_SOLVER_EIGENVALSOLVER_h

#include "matrix.h"
#include "sparseVect.h"
#include "vect.h"
#include <cmath> //std::abs

namespace math{
namespace solver{

#ifdef MYDEBUG
  class EigenValSolverTest;
#endif //MYDEBUG

template< template <typename> class Storage = SparseVect, typename T = double >
class EigenValSolver
{

public:
  
  Matrix<Vect,T>* EigVect()
  {
    if ( !m_EigVect.IsTransposed() )
    {
      m_EigVect.Transpose();
    }
    return &m_EigVect;
  }

  Matrix<SparseVect,T>* EigVal()
  {
    return &m_EigVal;
  }
  
protected:
  bool ExceedsOffDiagTol( Matrix< Vect, T >& A, T tol )
  {
    for ( unsigned long long i = 0; i < A.Rows(); ++i )
    {
      for ( unsigned long long j = 0; j < A[i].NNZ(); ++j )
      {
        if ( A[i].Pos(j) >= i )
        {
          break; // focus only on zeroing the lower half
        }
        if ( std::abs( A[i](j) ) > tol  )
        {
          return true;
        }
      }          
    }
    return false;
  }
  
protected:
  Matrix< Vect, T > m_EigVect;
  Matrix<SparseVect,T> m_EigVal;


#ifdef MYDEBUG
  friend class EigenValSolverTest;
#endif //MYDEBUG
};

}/*solver*/ }/*math*/ 

#endif /*_MATH_SOLVER_EIGENVALSOLVER_h */
