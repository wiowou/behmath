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

#ifndef _MATH_OPTIM_MINFINDER_h
#define _MATH_OPTIM_MINFINDER_h

#include "impl/optimConfig.h"
#include "typedefs.h"

namespace math{
namespace optim{

#ifdef MYDEBUG
  class MinFinderTest;
#endif //MYDEBUG

template<class F, class T = double>
class MinFinder
{

public:
  MinFinder()
  {
    Init();
  }
  
  void Clear()
  {
    Init();
  }

  void Tolerance( double tol )
  {
    m_tolerance = tol;
    m_runCalc = true;
    m_minFound = false;
  }

  double Tolerance()
  {
    return m_tolerance;
    m_minFound = false;
  }

  void Range( T &low, T &high )
  {
    m_low = low;
    m_high = high;
    m_runCalc = true;
    m_minFound = false;
  }

  void Function( F &f )
  {
    m_f = &f;
    m_runCalc = true;
    m_minFound = false;
  }
  
  void MaxIter( ULong iter )
  {
    m_maxIter = iter;
    m_runCalc = true;
    m_minFound = false;
  }

  ULong Iter()
  {
    return m_iter;
  }
  
  bool RootFound()
  {
    return m_minFound;
  }
  
  void Gain( T &gain )
  {
    m_gain = gain;
  }

  void Guess( T &guess )
  {
    m_guess = guess;
  }
  
protected:
  void Init()
  {
    m_maxIter = 0;
    m_tolerance = 1e-5;
    m_runCalc = true;
    m_minFound = false;    
  }
  
protected:
  F* m_f;
  double m_tolerance;
  T m_low;
  T m_high;
  T m_min;
  T m_gain;
  T m_guess;
  ULong m_maxIter;
  ULong m_iter;
  bool m_runCalc;
  bool m_minFound;
private:


#ifdef MYDEBUG
  friend class MinFinderTest;
#endif //MYDEBUG
};

} /*optim*/ }/*math*/ 

#endif /*_MATH_OPTIM_MINFINDER_h */
