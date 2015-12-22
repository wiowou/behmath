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

#ifndef _MATH_ROOT_BISECTION_h
#define _MATH_ROOT_BISECTION_h

#include "rootFinder.h"

namespace math{
namespace root{

#ifdef MYDEBUG
  class BisectionTest;
#endif //MYDEBUG

template<class F, class T = double>
class BisectionMidpoint : public RootFinder<F,T>
{
  using RootFinder<F,T>::m_x;
public:
  T Root()
  {
    return m_x[0];
  }
protected:
};

template<class F, template< typename F_, typename T_> class MP = BisectionMidpoint, class T = double>
class Bisection : public RootFinder<F,T>
{
  using RootFinder<F,T>::m_f;
  using RootFinder<F,T>::m_tolerance;
  using RootFinder<F,T>::m_low;
  using RootFinder<F,T>::m_high;
  using RootFinder<F,T>::m_root;
  using RootFinder<F,T>::m_maxIter;
  using RootFinder<F,T>::m_iter;
  using RootFinder<F,T>::m_runCalc;
  using RootFinder<F,T>::m_rootFound;
  

public:
  using RootFinder<F,T>::RootInRange;
  
  Bisection()
  {
    m_mid.MaxIter(1);
    m_useHybrid = true;
    m_switch = 3;
    m_fracTarget = 1.0;
    for ( ULong i = 0; i < m_switch; ++i )
    {
      m_fracTarget /= 2.0;
    }
  }
  
  //! Specify the function for which the root will be found
  void Function( F &f )
  {
    m_mid.Function(f);
    RootFinder<F,T>::Function(f);
  }
  
  //! Return the root. Calculate the root if it hasn't been
  //! calculated already.
  T Root()
  {
    if ( m_runCalc )
    {
      Calc();
    } 
    return m_root;
  }

protected:
  //The main routine to find the root
  void Calc()
  {
    T& low = m_low;
    T& high = m_high;
    T& mid = m_root; 
    
    if ( !RootInRange() )
    {
      return;
    }
    
    mid = ( high + low ) / 2.0;
    
    //added to allow for a better midpoint
    T mid2;
    m_mid.X0(mid);
    m_mid.X1(high);
    m_mid.X2(low);
    mid2 = m_mid.Root();
    mid = mid2 > high || mid2 < low ? mid : mid2;
    //
    T oldRange = ( high - low ) / m_fracTarget * 2.0;
    
    m_iter = 0;
    while ( m_iter < m_maxIter && std::abs( m_f(mid) ) > m_tolerance )
    {
      if ( m_f(mid) < T(0) )
      {
        low = mid;
      }
      else
      {
        high = mid;
      }
      mid = ( high + low ) / 2.0;
    
      // added to make sure the other method is keeping up with the bisection method.
      // Checks every m_switch iterations and gives hybrid method another chance every m_switch iterations
      if ( m_iter % m_switch == 0 )
      {
        if ( !m_useHybrid )
        {
          m_useHybrid = true;
        }  
        else if ( ( high - low ) > oldRange * m_fracTarget )
        {
          m_useHybrid = false;
        }
        oldRange = high - low;
      }
      
      if ( m_useHybrid )
      {
        //added to allow for a better midpoint
        T mid2;
        m_mid.X0(mid);
        m_mid.X1(high);
        m_mid.X2(low);
        mid2 = m_mid.Root();
        mid = mid2 > high || mid2 < low ? mid : mid2;
        //        
      }

      ++m_iter;
    }
    m_runCalc = false;
    if ( m_iter < m_maxIter )
    {
      m_rootFound = true;
    }
  }
protected:
  MP<F,T> m_mid;
  ULong m_switch;
  bool m_useHybrid;
  T m_fracTarget;
private:


#ifdef MYDEBUG
  friend class BisectionTest;
#endif //MYDEBUG
};

}/*root*/ }/*math*/ 

#endif /*_MATH_ROOT_BISECTION_h */
