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

#ifndef _MATH_ROOT_ROOTFINDER_h
#define _MATH_ROOT_ROOTFINDER_h

#include "impl/rootConfig.h"
#include "matrix/vect.h"
#include "typedefs.h"
#include <cmath> //std::abs

namespace math{
namespace root{

#ifdef MYDEBUG
  class RootFinderTest;
#endif //MYDEBUG

template<class F, class T = double>
class RootFinder
{

public:
  RootFinder()
  {
    m_maxIter = m_iter = 0;
    m_x.Resize(3);
    m_low = m_high = T(0);
    m_tolerance = 1e-5;
    m_runCalc = true;
    m_rootFound = false;
  }
  
  void Tolerance( T tol )
  {
    m_tolerance = tol;
    m_runCalc = true;
    m_rootFound = false;
  }
  
  T Tolerance()
  {
    return m_tolerance;
    m_rootFound = false;
  }
  
  //! Specify the range that the root must lie within.
  //! Not applicable for all methods.
  void Range( T low, T high )
  {
    m_low = low;
    m_high = high;
    m_runCalc = true;
    m_rootFound = false;
  }
  
  void X0( T x )
  {
    m_x[0] = x;
    m_runCalc = true;
    m_rootFound = false;
  }

  void X1( T x )
  {
    m_x[1] = x;
    m_runCalc = true;
    m_rootFound = false;
  }

  void X2( T x )
  {
    m_x[2] = x;
    m_runCalc = true;
    m_rootFound = false;
  }
  
  void Function( F &f )
  {
    m_f = f;
    m_runCalc = true;
    m_rootFound = false;
  }
  
  void MaxIter( ULong iter )
  {
    m_maxIter = iter;
    m_runCalc = true;
    m_rootFound = false;
  }

  ULong Iter()
  {
    return m_iter;
  }
  
  bool RootFound()
  {
    return m_rootFound;
  }
  
  void ExpandRange()
  {
    T dist = 0.5 * ( m_high - m_low );
    m_high += dist;
    m_low -= dist;
  }
  
  bool RootInRange()
  {
    return !( m_f(m_low) * m_f(m_high) > T(0) );
  }
  
protected:
  F m_f;
  T m_tolerance;
  T m_low;
  T m_high;
  T m_root;
  ULong m_maxIter;
  ULong m_iter;
  Vect<T> m_x;
  bool m_runCalc;
  bool m_rootFound;
private:


#ifdef MYDEBUG
  friend class RootFinderTest;
#endif //MYDEBUG
};

}/*root*/ }/*math*/ 

#endif /*_MATH_ROOT_ROOTFINDER_h */
