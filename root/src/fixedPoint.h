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

#ifndef _MATH_ROOT_FIXEDPOINT_h
#define _MATH_ROOT_FIXEDPOINT_h

#include "rootFinder.h"

namespace math{
namespace root{

#ifdef MYDEBUG
  class FixedPointTest;
#endif //MYDEBUG

//! Uses the fixed-point method to find the root of f(x)
template<class F, class T = double>
class FixedPoint : public RootFinder<F,T>
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
    m_iter = 0;
    m_root = G( m_x[0] );
    while ( m_iter < m_maxIter && std::abs( m_f(m_root) ) > m_tolerance )
    {
      m_root = G(m_root);
      ++m_iter;
    }
    m_runCalc = false;
    if ( m_iter < m_maxIter )
    {
      m_rootFound = true;
    }
  }
  
  //! The original function is modified so that x = x + f(x) in order
  //! that the user need specify only f(x).
  T G( T x )
  {
    return x + m_f(x);
  }
private:


#ifdef MYDEBUG
  friend class FixedPointTest;
#endif //MYDEBUG
};

}/*root*/ }/*math*/ 

#endif /*_MATH_ROOT_FIXEDPOINT_h */
