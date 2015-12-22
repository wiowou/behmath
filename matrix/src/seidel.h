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

#ifndef _MATH_SOLVER_SEIDEL_h
#define _MATH_SOLVER_SEIDEL_h

#include "iterative.h"
#include "solver.h"
#include "subtr.h"
#include "dot.h"
#include <cmath> // std::abs

namespace math{
namespace solver{

#ifdef MYDEBUG
  class SeidelTest;
#endif //MYDEBUG

template< template <typename> class Storage = SparseVect, typename T = double >
class Seidel : public Iterative<Storage,T>
{
  using Iterative<Storage,T>::m_A;
  using Iterative<Storage,T>::m_x;
  using Iterative<Storage,T>::m_b;
  using Iterative<Storage,T>::m_notWithinTol;
  using Iterative<Storage,T>::m_xNew;
  using Iterative<Storage,T>::m_tolerance;
  using Iterative<Storage,T>::m_maxIter;
  using Iterative<Storage,T>::m_nIter;
  using Iterative<Storage,T>::Continue;
  using Iterative<Storage,T>::InitialGuess;
  
public:
  void Solve()
  {
    Matrix< Storage, T >& L = *m_A;
    Matrix< Storage, T >& U = *m_A;
    Vect<T>& x = *m_x;
    Vect<T>& b = *m_b;
    Vect<T>& xCheck = m_xNew;
    vops::Dot<T> dot;
    vops::Subtr<T> subtr;
    Solver<Storage,T> solver;
    
    InitialGuess();
    m_nIter = 0;
    while ( m_nIter < m_maxIter && Continue() )
    {
      // [x] = [U][x]
      for ( ULong i = 0; i < U.Rows(); ++i )
      {
        x[i] = dot( U[i], x, i + 1 );
      }
      
      // [x] = [b] - [x]
      subtr( b, x, x );
      
      // [xNew] = L^(-1)[x]
      solver.BackSolveDown( L, x, x );
      
      //update tolerance array
      UpdateNotWithinTol( xCheck, x );
      
      //We don't need xCheck or m_notWithinTol if a residual is used instead of an element-by-element tol check
      //Can save some space
      xCheck = x;
      
      ++m_nIter;
    }
  }
  
protected:
  void UpdateNotWithinTol( Vect<T> &v1, Vect<T> &v2 )
  {
    // can be parallelized
    for ( ULong i = 0; i < m_notWithinTol.Size(); ++i )
    {
      m_notWithinTol[i] = std::abs( v1[i] - v2[i] ) > m_tolerance;
    }
  }
  
private:


#ifdef MYDEBUG
  friend class SeidelTest;
#endif //MYDEBUG
};

}/*solver*/ }/*math*/ 

#endif /*_MATH_SOLVER_SEIDEL_h */
