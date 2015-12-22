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

#ifndef _MATH_ROOT_LINFRAC_h
#define _MATH_ROOT_LINFRAC_h

#include "rootFinder.h"
#include "matrix/elimination.h"
#include "matrix/matrix.h"
#include "matrix/vect.h"

namespace math{
namespace root{

#ifdef MYDEBUG
  class LinFracTest;
#endif //MYDEBUG

template<class F, class T = double>
class LinFrac : public RootFinder<F,T>
{
  using RootFinder<F,T>::m_f;
  using RootFinder<F,T>::m_tolerance;
  using RootFinder<F,T>::m_x;
  using RootFinder<F,T>::m_root;
  using RootFinder<F,T>::m_maxIter;
  using RootFinder<F,T>::m_iter;
  using RootFinder<F,T>::m_runCalc;
  using RootFinder<F,T>::m_rootFound;
  
public:
  LinFrac()
  {
    m_A.Resize(3);
    m_A[0][0] = 1.0;
    m_A[1][0] = 1.0;
    m_A[2][0] = 1.0;
    m_x.Resize(3);
    m_b.Resize(3);
  }
  
  T Root()
  {
    if ( m_runCalc )
    {
      Calc();
    } 
    return m_root;
  }
  
protected:
  void Calc()
  {
    Matrix<Vect,T>& A = m_A;
    Vect<T>& x = m_x;
    Vect<T>& b = m_b;

    A[0][2] = -m_f(x[0]);
    A[1][2] = -m_f(x[1]);
    A[2][2] = -m_f(x[2]);
    
    A[0][1] = -A[0][2] * x[0];
    A[1][1] = -A[1][2] * x[1];
    A[2][1] = -A[2][2] * x[2];
    
    solver::Elimination<Vect,T> sol;
    sol.Tolerance(1e-5);
    sol.PartialPivoting(true);
    
    Matrix<Vect,T> Acpy = A;
    sol.A(Acpy);
    Vect<T> xcpy = x;
    sol.B(xcpy); //the vector of outputs, b, is the known input x
    sol.X(b); //we are solving for the coefficients, b, of the parabola
    sol.Solve();
    m_root = b[0];
    
    m_iter = 0;
    while ( m_iter < m_maxIter && std::abs( m_f( m_root ) ) > m_tolerance )
    {
      x[0] = x[1];
      x[1] = x[2];
      x[2] = b[0];
      
      A[0][2] = A[1][2];
      A[1][2] = A[2][2];
      A[2][2] = -m_f(x[2]);

      A[0][1] = A[1][1];
      A[1][1] = A[2][1];
      A[2][1] = -A[2][2] * x[2];
      
      Acpy = A;
      xcpy = x;
      sol.Solve();
      m_root = b[0];
    
      ++m_iter;
    }
    m_runCalc = false;
    if ( m_iter < m_maxIter )
    {
      m_rootFound = true;
    }
  }
  
protected:
  Vect<T> m_b;
  Matrix<Vect,T> m_A;

private:


#ifdef MYDEBUG
  friend class LinFracTest;
#endif //MYDEBUG
};

}/*root*/ }/*math*/ 

#endif /*_MATH_ROOT_LINFRAC_h */
