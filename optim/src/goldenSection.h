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

#ifndef _MATH_OPTIM_GOLDENSECTION_h
#define _MATH_OPTIM_GOLDENSECTION_h

#include "minFinder.h"

namespace math{
namespace optim{

#ifdef MYDEBUG
  class GoldenSectionTest;
#endif //MYDEBUG

template<class F, class T = double>
class GoldenSection : public MinFinder<F,T>
{

  using MinFinder<F,T>::m_f;
  using MinFinder<F,T>::m_tolerance;
  using MinFinder<F,T>::m_min;
  using MinFinder<F,T>::m_maxIter;
  using MinFinder<F,T>::m_iter;
  using MinFinder<F,T>::m_runCalc;
  using MinFinder<F,T>::m_minFound;
  using MinFinder<F,T>::m_low;
  using MinFinder<F,T>::m_high;
  
public:
  T MinX()
  {
    if ( m_runCalc )
    {
      Calc();
    } 
    return m_min;
  }
  
protected:
  void Calc()
  {
    const double tau = 0.381966;
    const double tau1 = 0.618034;
    double a = m_low;
    double b = m_high;
    double range = b - a;
    double x1 = a + tau * range;
    double x2 = a + tau1 * range;
    double f1 = (*m_f)(x1);
    double f2 = (*m_f)(x2);
    
    m_iter = 0;
    while ( range > m_tolerance && m_iter < m_maxIter )
    {
      if ( f1 > f2 )
      {
        a = x1;
        x1 = x2;
        f1 = f2;
        range = b - a;
        x2 = a + tau1 * range;
        f2 = (*m_f)(x2);
      }
      else
      {
        b = x2;
        x2 = x1;
        f2 = f1;
        range = b - a;
        x1 = a + tau * range;
        f1 = (*m_f)(x1);
      }
      ++m_iter;
    }
    m_min = 0.5 * ( a + b );
    m_runCalc = false;
    if ( m_iter < m_maxIter )
    {
      m_minFound = true;
    }
  }
private:


#ifdef MYDEBUG
  friend class GoldenSectionTest;
#endif //MYDEBUG
};

}/*optim*/ }/*math*/ 

#endif /*_MATH_OPTIM_GOLDENSECTION_h */
