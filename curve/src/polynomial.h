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

#ifndef _MATH_CURVE_POLYNOMIAL_h
#define _MATH_CURVE_POLYNOMIAL_h

#include "curve.h"

namespace math{
namespace curve{

#ifdef MYDEBUG
  class PolynomialTest;
#endif //MYDEBUG

//y=a0+a1*x1+a2*x2^2...an*xn^n
class Polynomial : public Curve
{
public:
  Polynomial() {}
  explicit Polynomial( ULong order )
  {
    Order( order );
  }
  
  //! Specify the order of the polynomial, such as 2 for a parabola, 3 for a cubic, etc
  void Order( ULong order )
  {
    m_coeff.Resize( order + 1 );
  }
  
  //! Returns the order of the polynomial
  ULong Order() const
  {
    return m_coeff.Size();
  }
  
	//! Get the value of the function
  double operator()( double x )
  {
    double sum = m_coeff[0];
    for ( ULong i = 1; i < m_coeff.Size(); ++i )
    {
      sum += m_coeff[i] * pow( x, i );
    }
    return sum;
  }
  
  //! First derivative
  double Deriv1( double x )
  {
    double sum = m_coeff[1];
    for ( ULong i = 2; i < m_coeff.Size(); ++i )
    {
      sum += static_cast<double>(i) * m_coeff[i] * pow( x, i - 1 );
    }
    return sum;
  }
  
  //! Second derivative
  double Deriv2( double x )
  {
    double sum = m_coeff[2];
    for ( ULong i = 3; i < m_coeff.Size(); ++i )
    {
      sum += static_cast<double>(i) * static_cast<double>(i - 1) * m_coeff[i] * pow( x, i - 2 );
    }
    return sum;
  }

	//! Fit a set of x and y points to the function
  template< template< template <typename> class Storage, typename T > class Sol >
  void Fit( Vect<double>& x, Vect<double>& y )
  {
    Curve::Fit<Sol>( x, y );
    ErrorEuclidian( x, y );
  }
  
	//! Print the function
  friend std::ostream& operator<<( std::ostream& os, const Polynomial& p )
  {
    os << " y(x) = ";
    ULong i = p.Order() - 1;
    os << p.Coeff(i) << " * x^" << i;
    --i;
    for ( ; i > 1; --i )
    {
      if ( p.Coeff(i) > 0.0 )
      {
        os << " + ";
      }
      else
      {
        os << " - ";
      }
      os << std::abs( p.Coeff(i) ) << " * x^" << i;
    }
    if ( p.Coeff(i) > 0.0 )
    {
      os << " + ";
    }
    else
    {
      os << " - ";
    }
    os << std::abs( p.Coeff(i) ) << " * x";
    --i;
    if ( p.Coeff(i) > 0.0 )
    {
      os << " + ";
    }
    else
    {
      os << " - ";
    }
    os << std::abs( p.Coeff(i) );
    return os;
  }
protected:
  void CreateAT( Vect<double>& x, Matrix<Vect,double>& AT )
  {
    /*
    for ( ULong j = 0; j < m_coeff.Size(); ++j )
    {
      for ( ULong i = 0; i < x.Size(); ++i )
      {
        AT[j][i] = pow( x[i], j );
      }
    }
    */
    for ( ULong i = 0; i < x.Size(); ++i )
    {
      AT[0][i] = 1.0;
      for ( ULong j = 1; j < m_coeff.Size(); ++j )
      {
        AT[j][i] *= x[i];
      }
    }
  }
  

private:
#ifdef MYDEBUG
  friend class CurveTest;
#endif //MYDEBUG
};

}/*curve*/ }/*math*/ 

#endif /*_MATH_CURVE_POLYNOMIAL_h */
