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

#ifndef _MATH_OPTIM_BIGD_h
#define _MATH_OPTIM_BIGD_h

#include "minFinder.h"
#include "matrix/vect.h"
#include <cmath> //std::abs

namespace math{
namespace optim{

#ifdef MYDEBUG
  class BigdTest;
#endif //MYDEBUG

//! Bounded Inertial Gradient Descender
/*!
      - Reverse direction on contact with boundary
      - Set the initial time step size
      - stop when your gradient vector is below tolerance
*/
template<class F>
class Bigd : public MinFinder<F,Vect<double> >
{

  using MinFinder<F,Vect<double> >::m_f;
  using MinFinder<F,Vect<double> >::m_tolerance;
  using MinFinder<F,Vect<double> >::m_min;
  using MinFinder<F,Vect<double> >::m_maxIter;
  using MinFinder<F,Vect<double> >::m_iter;
  using MinFinder<F,Vect<double> >::m_runCalc;
  using MinFinder<F,Vect<double> >::m_minFound;
  using MinFinder<F,Vect<double> >::m_low;
  using MinFinder<F,Vect<double> >::m_high;
  using MinFinder<F,Vect<double> >::m_guess;
  
public:
  Bigd()
  {
    Init();
  }
  
  //! Specify the size of the step taken. Smaller steps
  //! are generally more stable but slower.
  void Dt( double d )
  {
    m_dt = d;
  }
  
  double Dt()
  {
    return m_dt;
  }  

  //! Factor placed on acceleration. Initialized to 1.0.
  void Cf( double d )
  {
    m_cf = d;
  }
  
  double Cf()
  {
    return m_cf;
  }
  
  //! The minimum function location within the Range is returned
  Vect<double>* MinX()
  {
    if ( m_runCalc )
    {
      Calc();
    } 
    return &m_min;
  }

  //! Set the lower and upper bounds of the range so roots are not
  //! found outside this range.
  void Range( Vect<double> &low, Vect<double> &high )
  {
    MinFinder<F,Vect<double> >::Range( low, high );
    m_guess.Resize( low.Size() );
    m_vel.Resize( low.Size() );
    m_acc.Resize( low.Size() );
    for ( ULong i = 0; i < m_guess.Size(); ++i )
    {
      m_guess[i] = 0.5 * ( high[i] + low[i] );
      m_dt = m_dt > 0.1 * ( high[i] - low[i] ) ? 0.1 * ( high[i] - low[i] ) : m_dt;
    }
    m_minVal = (*m_f)( m_guess );
  }
  
protected:
  void Calc()
  {
    Vect<double> &pos = m_guess;
    Vect<double> &vel = m_vel;
    Vect<double> &acc = m_acc;
    double dt = m_dt;
    
    m_f->Grad( pos, acc );
    m_iter = 0;
    while ( Continue( acc, m_tolerance ) && m_iter < m_maxIter )
    {
      for ( ULong i = 0; i < pos.Size(); ++i )
      {
        pos[i] -= acc[i] * dt + vel[i] * dt;
        vel[i] += m_cf * acc[i] * dt;
      }
      Bounce( pos, vel );
      
      double val = (*m_f)( pos );
      if ( val < m_minVal )
      {
        m_min = pos;
        m_minVal = val;
      }
      
      m_f->Grad( pos, acc );
      ++m_iter;
    }
    
    m_runCalc = false;
    if ( m_iter < m_maxIter )
    {
      m_minFound = true;
    }
  }
  
  //! Keeps iterating while the acc=grad is large enough
  bool Continue( Vect<double>& acc, double tol )
  {
    for ( ULong i = 0; i < acc.Size(); ++i )
    {
      if ( std::abs( acc[i] ) > tol )
      {
        return true;
      }
    }
    return false;
  }
  
  //! Will bounce off the walls of the Range to remain within the Range.
  void Bounce( Vect<double>& pos, Vect<double>& vel )
  {
    for ( ULong i = 0; i < pos.Size(); ++i )
    {
      if ( pos[i] < m_low[i] )
      {
        pos[i] = m_low[i];
        vel[i] *= -1.0;
      }
      else if ( pos[i] > m_high[i] )
      {
        pos[i] = m_high[i];
        vel[i] *= -1.0;
      }
    }
  }
  
  void Init()
  {
    m_dt = 1.0;
    m_cf = 1.0;
    m_tolerance = -1.0;
  }
private:
  double m_minVal;
  double m_dt;
  double m_cf;
  Vect<double> m_acc;
  Vect<double> m_vel;
  
#ifdef MYDEBUG
  friend class BigdTest;
#endif //MYDEBUG
};



}/*optim*/ }/*math*/ 

#endif /*_MATH_OPTIM_BIGD_h */
