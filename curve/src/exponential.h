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

#ifndef _MATH_CURVE_EXPONENTIAL_h
#define _MATH_CURVE_EXPONENTIAL_h

#include "curve.h"

namespace math{
namespace curve{

#ifdef MYDEBUG
  class ExponentialTest;
#endif //MYDEBUG

//y=a*exp(b*x). Relies on the Curve class to create a size of 2 for m_coeff
class Exponential : public Curve
{
public:
  //! Get the value of the function
  double operator()( double x )
  {
    return m_coeff[0] * exp( m_coeff[1] * x );
  }

  //! Will fit an exponential to the data, where a = m_coeff[0] and b = m_coeff[1]
  template< template< template <typename> class Storage, typename T > class Sol >
  void Fit( Vect<double>& x, Vect<double>& y )
  {
    for ( ULong i = 0; i < y.Size(); ++i )
    {
      y[i] = log( y[i] );
    }
    Curve::Fit<Sol>( x, y );
    m_coeff[0] = exp( m_coeff[0] );
    ErrorEuclidian( x, y );
  }
  
  friend std::ostream& operator<<( std::ostream& os, const Exponential& p )
  {
    os << " y(x) = ";
    os << p.Coeff(0) << " e^(" << p.Coeff(1) << "x)" << std::endl;
    return os;
  }
protected:
  void CreateAT( Vect<double>& x, Matrix<Vect,double>& AT )
  {
    for ( ULong i = 0; i < x.Size(); ++i )
    {
      AT[0][i] = log( x[i] );
      AT[1][i] = x[i];
    }
  }
private:
#ifdef MYDEBUG
  friend class CurveTest;
#endif //MYDEBUG
};

}/*curve*/ }/*math*/ 

#endif /*_MATH_CURVE_EXPONENTIAL_h */
