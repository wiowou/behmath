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

#ifndef _MATH_ODE_EULEREXPL_h
#define _MATH_ODE_EULEREXPL_h

#include "impl/odeConfig.h"
#include "surface/surface.h"
#include "solver.h"
#include "matrix/vect.h"
#include "matrix/dot.h"
#include <cmath> //std::abs

namespace math{
namespace ode{

#ifdef MYDEBUG
  class EulerExplTest;
#endif //MYDEBUG

template<class T>
class EulerExpl : Solver<T>
{

public:

protected:

private:


#ifdef MYDEBUG
  friend class EulerExplTest;
#endif //MYDEBUG
};

/*! 
This is the template to use for ODE's that are a function of time and x parameter:
dx/dt = F( x(t), t ). One example is dx/dt = 2x + cos(t)
*/
template<>
class EulerExpl<Surface> : public Solver<double>
{
  using Solver<double>::m_initialCondition;
  using Solver<double>::m_tolerance;
  using Solver<double>::m_tInitial;
  using Solver<double>::m_tFinal;
  using Solver<double>::m_steps;
  using Solver<double>::m_maxSteps;
  using Solver<double>::m_fLow;
  using Solver<double>::m_fUp;
  
public:
  EulerExpl()
  {
    Init();
  }
  
  void F( Surface* s )
  {
    m_f = s;
  }
  
  //! Uses the Huen-Euler butcher tableau for a 1st and 2nd order Runge Kutta and compares the two. 
  //! Adapts the time step by 2x or 1/2x according to RK1 and RK2 difference comparison.
  void Solve()
  {
    vops::Dot<double> dot;
    double resultRK1;
    double resultRK2 =  m_initialCondition;
    double lastResult =  m_initialCondition;
    double dt = m_tFinal - m_tInitial;
    double dtDouble = 0.1 * m_tolerance;
    int dtDoubleCounter = 0;
    m_steps = 0;
    Vect<double> x(2);
    double epsilon = ( m_tFinal - m_tInitial ) / static_cast<double>( m_maxSteps * 100 );
    Vect<double> k(2);

    Vect<double> bRK1(2);
    Vect<double> bRK2(2);
    bRK1[0] = 1.0; bRK1[1] = 0.0;
    bRK2[0] = 0.5; bRK2[1] = 0.5;   
    
    m_time = m_tFinal;
    while ( m_time < m_tFinal + epsilon && m_steps < m_maxSteps )
    {
      x[1] = m_time;
      
      x[0] = resultRK2;
      k[0] = (*m_f)(x);
      
      x[0] = resultRK2 + k[0] * dt;
      k[1] = (*m_f)(x);
      
      resultRK1 = resultRK2 + dot( k, bRK1 ) * dt;
      lastResult = resultRK2;
      resultRK2 = resultRK2 + dot( k, bRK2 ) * dt;
      
      double error = std::abs( resultRK2 - resultRK1 );
      //Adaptive time-stepping
      if ( error < m_tolerance )
      {
        if ( !ResultInRange(resultRK2) )
        {
          m_result = lastResult;
          resultRK2 = lastResult;
          m_time -= dt;
          dt *= 0.5;
          if ( std::abs( m_result - m_fLow ) < m_rangeTol || std::abs( m_result - m_fUp ) < m_rangeTol )
          {
            return;
          }  
        }
        //The dtDoubleCounter prevents successive doubling and halving of the time step
        if ( error < dtDouble && dtDoubleCounter > 16 )
        {
          dt *= 2.0;
          dtDoubleCounter = 0;
        }
      }
      else
      {
        m_time -= dt;
        resultRK2 = lastResult;
        dt *= 0.5;
      }
      
      if ( m_time + dt > m_tFinal )
      {
        dt = m_tFinal - m_time < epsilon ? epsilon : m_tFinal - m_time;
      }
      m_time += dt;
      ++dtDoubleCounter;
      ++m_steps;
    }
    m_result = resultRK2;
  }

protected:
  bool ResultInRange( double f )
  {
    return ( f > m_fLow && f < m_fUp );
  }
  void Init()
  {
    m_f = NULL;
    Solver<double>::Init();
    double m_dtChange = 0.1 * m_tolerance;
    m_fLow = -1e30;
    m_fUp = 1e30;
  }
protected:
  Surface* m_f;
private:


#ifdef MYDEBUG
  friend class FehlbergTest;
#endif //MYDEBUG
};
}/*ode*/ }/*math*/ 

#endif /*_MATH_ODE_EULEREXPL_h */
