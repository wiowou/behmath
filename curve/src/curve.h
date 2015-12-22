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

#ifndef _MATH_CURVE_h
#define _MATH_CURVE_h

#include "impl/curveConfig.h"
#include "matrix/QR.h"
#include "matrix/rectangular.h"
#include "matrix/matrix.h"
#include "matrix/vect.h"
#include "matrix/dot.h"
#include "typedefs.h"
#include <math.h> //log, exp, sqrt functions
#include <cmath> //std::abs for printing
#include <iostream> //for printing

namespace math{

#ifdef MYDEBUG
  class CurveTest;
#endif //MYDEBUG
class Curve
{
public:
  Curve()
  {
    Init(2);
  }
  
  void Clear()
  {
    Init( m_coeff.Size() );
  }

	//! Uses the Sol template class, which is a matrix solver like QR or Elimination, to fit a curve to x and y
  template< template< template <typename> class Storage, typename T > class Sol >
  void Fit( Vect<double>& x, Vect<double>& y )
  {
    Sol<Vect,double> s;
    Fit( x, y, s );
  }
  
  //! Get a coefficient at an index
  double Coeff( const ULong idx ) const
  {
    return m_coeff.Get(idx);
  }
  
  //! Calculate an error estimate between the generated curve and the actual data
  double Error()
  {
    return m_error;
  }
  
  //! First derivative, central difference
  virtual double Deriv1( double x )
  {
    return 0.5 * ( operator()( x + m_deltaX ) - operator()( x - m_deltaX ) ) / m_deltaX;
  }
  
  //! Second derivative, central difference
  virtual double Deriv2( double x )
  {
    double d = operator()( x + m_deltaX ) - operator()( x ) - operator()( x ) + operator()( x - m_deltaX );
    d /= m_deltaX;
    return d / m_deltaX;
  }
  
  //! use this to return y(x)
  virtual double operator()( double x )
  {
    return 0.0;
  }
protected:
  //! This method does sometimes change x and y in calculating the fit
  template< template< template <typename> class Storage, typename T > class Sol >
  void Fit( Vect<double>& x, Vect<double>& y, Sol<Vect,double> s )
  {
    Matrix<Vect,double> AT( m_coeff.Size(), x.Size() );
    //creating A transpose
    CreateAT(x, AT);
    
    //specifying the input parameters to the rectangular solver
    solver::Rectangular<Sol,Vect,double> solver;
    solver.ATranspose(AT);
    solver.B(y);
    solver.X(m_coeff);
    solver.Solve();
  }

	//! Solves for the coefficients of the various curves.
  void Fit( Vect<double>& x, Vect<double>& y, solver::QR<Vect,double> s )
  {
    Matrix<Vect,double> AT( m_coeff.Size(), x.Size() );
    //creating A transpose
    CreateAT(x, AT);
    
    s.ATranspose(AT);
    s.B(y);
    s.X(m_coeff);
    s.Solve();
  }
  
  //! create A transpose from x
  virtual void CreateAT( Vect<double>& x, Matrix<Vect,double>& AT )
  {
    return;
  }
  
  //! calculation:  sqrt( dot( curve(x) - y ) ) / #pts which calculates the Euclidian norm between 
  //! the curve and the actual points and divides by the number of points
  void ErrorEuclidian( Vect<double>& x, Vect<double>& y )
  {
    vops::Dot<double> dot;
    Vect<double> err( y.Size() );
    for ( ULong i = 0; i < err.Size(); ++i )
    {
      err[i] = operator()( x[i] ) - y[i];
    }
    m_error = sqrt( dot( err, err ) ) / static_cast<double>( err.Size() );
  }

  void Init( ULong size = 2 )
  {
    m_coeff.Clear();
    m_coeff.Resize(size);
    m_error = 0.0;
    m_deltaX = 1e-4;
  }
protected:
  Vect<double> m_coeff;
  double m_error;
  double m_deltaX;
private:


#ifdef MYDEBUG
  friend class CurveTest;
#endif //MYDEBUG
};

}/*math*/ 

#endif /*_MATH_CURVE_h */
