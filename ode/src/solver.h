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

#ifndef _MATH_ODE_SOLVER_h
#define _MATH_ODE_SOLVER_h

#include "impl/odeConfig.h"
#include "typedefs.h"

namespace math{
namespace ode{

#ifdef MYDEBUG
  class SolverTest;
#endif //MYDEBUG

template<class T>
class Solver
{

public:
  Solver()
  {
    Init();
  }
  
  void Clear()
  {
    Init();
  }
  
  //! Provide an initial condition to the ode
  void InitialCondition( T &ic )
  {
    m_initialCondition = ic;
  }
  
  //! Adaptive time-stepping tolerance
  void Tolerance( double tol)
  {
    m_tolerance = tol;
  }
  
  //! Set the bounds for the integration
  void Range( double tInitial, double tFinal )
  {
    m_tInitial = tInitial;
    m_tFinal = tFinal;
  }
  
  //! Set the bounds for the result of the integration.
  //! If the integration exceeds the bounds, the process will stop.
  void ResultRange( T low, T up, T rangeTol )
  {
    m_fLow = low;
    m_fUp = up;
    m_rangeTol = rangeTol;
  }
  
  ULong Steps()
  {
    return m_steps;
  }
  
  //! Set the maximum number of steps to take during integration
  void MaxSteps( ULong maxSteps )
  {
    m_maxSteps = maxSteps;
  }
  
  //! Return the current time step
  double Time()
  {
    return m_time;
  }
  
  T Result()
  {
    return m_result;
  }
protected:
  void Init()
  {
    m_tolerance = 1e-4;
    m_tInitial = 0.0;
    m_tFinal = 0.0;
    m_steps = 0;
    m_maxSteps = 1;
  }
protected:
  T m_initialCondition;
  T m_result;
  T m_fLow;
  T m_fUp;
  T m_rangeTol;
  ULong m_steps;
  ULong m_maxSteps;
  double m_tolerance;
  double m_tInitial;
  double m_tFinal;
  double m_time;

private:


#ifdef MYDEBUG
  friend class SolverTest;
#endif //MYDEBUG
};

}/*ode*/ }/*math*/ 

#endif /*_MATH_ODE_SOLVER_h */
