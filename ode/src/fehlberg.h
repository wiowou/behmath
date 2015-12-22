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

#ifndef _MATH_ODE_FEHLBERG_h
#define _MATH_ODE_FEHLBERG_h

#include "impl/odeConfig.h"
#include "surface/surface.h"
#include "solver.h"
#include "matrix/vect.h"
#include "matrix/dot.h"
#include <cmath> //std::abs

namespace math{
namespace ode{

#ifdef MYDEBUG
  class FehlbergTest;
#endif //MYDEBUG

template<class T>
class Fehlberg
{

public:

protected:

private:


#ifdef MYDEBUG
  friend class FehlbergTest;
#endif //MYDEBUG
};

/*! 
This is the template to use for ODE's that are a function of a time and x parameter:
dx/dt = F( x(t), t ). One example is dx/dt = 2x + cos(t)
*/
template<>
class Fehlberg<double> : Solver<double>
{
  using Solver<double>::m_initialCondition;
  using Solver<double>::m_tolerance;
  using Solver<double>::m_tInitial;
  using Solver<double>::m_tFinal;
  using Solver<double>::m_steps;
  using Solver<double>::m_maxSteps;
  using Solver<double>::m_fLow;
  using Solver<double>::m_fUp;
  using Solver<double>::m_result;
  
public:
  Fehlberg()
  {
    Init();
  }
  
  void F( Surface* s )
  {
    m_f = s;
  }
  
  //! Uses the Fehlberg butcher tableau for a 4th and 5th order Runge Kutta and compares the two. 
  //! Adapts the time step by 2x or 1/2x according to RK4 and RK5 difference comparison.
  void Solve()
  {
    vops::Dot<double> dot;
    double resultRK4;
    double resultRK5 =  m_initialCondition;
    double lastResult =  m_initialCondition;
    double dt = m_tFinal - m_tInitial;
    double dtDouble = 0.1 * m_tolerance;
    int dtDoubleCounter = 0;
    m_steps = 0;
    Vect<double> x(2);
    double epsilon = ( m_tFinal - m_tInitial ) / static_cast<double>( m_maxSteps * 100 );
    Vect<double> k(6);
    
    Vect<double> a[5];
    for ( int i = 0; i < 5; ++i )
    {
      a[i].Resize(i+1);
    }
    a[0][0] = 0.25;
    a[1][0] = 0.09375; a[1][1] = 0.28125;
    a[2][0] = 0.879380974056; a[2][1] = -3.277196176604; a[2][2] = 3.320892125626;
    a[3][0] = 2.032407407407; a[3][1] = -8.0; a[3][2] = 7.173489278752; a[3][3] = -0.205896686160;
    a[4][0] = -0.296296296296; a[4][1] = 2.0; a[4][2] = -1.381676413255; a[4][3] = 0.452972709552; a[4][4] = -0.275;
    
    Vect<double> bRK4(6);
    Vect<double> bRK5(6);
    bRK4[0] = 0.115740741; bRK4[1] = 0.0; bRK4[2] = 0.548927875; bRK4[3] = 0.535331384; bRK4[4] = -0.2; bRK4[5] = 0.0;
    bRK5[0] = 0.118518519; bRK5[1] = 0.0; bRK5[2] = 0.518986355; bRK5[3] = 0.50613149; bRK5[4] = -0.18; bRK5[5] = 0.036363636;   
    
    m_time = m_tFinal;
    while ( m_time < m_tFinal + epsilon && m_steps < m_maxSteps && m_f->Continue() )
    {
      x[1] = m_time;
      
      x[0] = resultRK5;
      k[0] = (*m_f)(x);
      
      for ( int i = 0; i < 5; ++i )
      {
        x[0] = resultRK5 + dot( i + 1, a[i], k) * dt;
        k[i + 1] = (*m_f)(x);
      }
      resultRK4 = resultRK5 + dot( k, bRK4 ) * dt;
      lastResult = resultRK5;
      resultRK5 = resultRK5 + dot( k, bRK5 ) * dt;
      
      double error = std::abs( resultRK5 - resultRK4 );
      //Adaptive time-stepping
      if ( error < m_tolerance )
      {
        if ( !ResultInRange(resultRK5) )
        {
          m_result = lastResult;
          resultRK5 = lastResult;
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
        resultRK5 = lastResult;
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
    m_result = resultRK5;
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

#endif /*_MATH_ODE_FEHLBERG_h */
