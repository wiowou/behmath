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

#ifndef _MATH_SURFACE_h
#define _MATH_SURFACE_h

#include "impl/surfaceConfig.h"
#include "matrix/matrix.h"
#include "matrix/vect.h"

namespace math{

#ifdef MYDEBUG
  class SurfaceTest;
#endif //MYDEBUG

class Surface
{

public:
  Surface()
  {
    m_deltaX = 1e-4;
  }
  
  void DeltaX( double dx )
  {
    m_deltaX = dx;
  }
  
  //! The value of DeltaX is used as the finite difference step size
  //! to calcualte the gradient. Default value is 1e-4.
  double DeltaX()
  {
    return m_deltaX;
  }
  
  //! The surface returns a double for each input vector, x
  virtual double operator()( Vect<double>& x )
  {
    return 0.0;
  }
  
  virtual bool Continue()
  {
    return true;
  }
  
  //! Gradient vector of the surface. Should be overloaded if an
  //! explicit expression for the gradient vector is known.
  //! Uses a central difference scheme, 2nd order accurate
  virtual void Grad( Vect<double> &x, Vect<double> &grad )
  {
    grad.Resize( x.Size() );
    double u[2];
    for ( ULong i = 0; i < x.Size(); ++i )
    {
      x[i] += m_deltaX;
      u[1] = operator()( x );
      x[i] -= m_deltaX;
      x[i] -= m_deltaX;
      u[0] = operator()( x );
      grad[i] = ( u[1] - u[0] ) / m_deltaX;
      grad[i] *= 0.5;
      x[i] += m_deltaX;
    }
  }
  
  //! Hessian of the surface. Hessian is the second mixed partials
  //! Uses a central difference scheme, 2nd order accurate
  virtual void Hessian( Vect<double> &x, Matrix<Vect> &H )
  {
    H.Resize( x.Size() );
    double u[2][2];
    double f = operator()( x );
    
    for ( ULong i = 0; i < H.Rows(); ++i )
    {
      for ( ULong j = i + 1; j < H.Cols(); ++j )
      {
        x[i] += m_deltaX;
        x[j] += m_deltaX;
        u[1][1] = operator()( x );
        x[j] -= m_deltaX;
        x[j] -= m_deltaX;
        u[1][0] = operator()( x );
        x[i] -= m_deltaX;
        x[i] -= m_deltaX;
        u[0][0] = operator()( x );
        x[j] += m_deltaX;
        x[j] += m_deltaX;
        u[0][1] = operator()( x );
        x[j] -= m_deltaX;
        x[i] += m_deltaX;
        H[i][j] = u[1][1] - u[1][0] - u[0][1] + u[0][0];
        H[i][j] /= m_deltaX;
        H[i][j] /= m_deltaX;
        H[i][j] *= 0.25;
        H[j][i] = H[i][j];
      }
      x[i] += m_deltaX;
      u[1][1] = operator()( x );
      x[i] -= m_deltaX;
      x[i] -= m_deltaX;
      u[0][0] = operator()( x );
      H[i][i] = u[1][1] - f - f + u[0][0];
      H[i][i] /= m_deltaX;
      H[i][i] /= m_deltaX;
      x[i] += m_deltaX;
    }
  }
protected:
  double m_deltaX;
private:


#ifdef MYDEBUG
  friend class SurfaceTest;
#endif //MYDEBUG
};

}/*math*/ 

#endif /*_MATH_SURFACE_h */
