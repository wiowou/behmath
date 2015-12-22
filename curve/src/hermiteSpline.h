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

#ifndef _MATH_CURVE_HERMITESPLINE_h
#define _MATH_CURVE_HERMITESPLINE_h

#include "curve.h"

namespace math{
namespace curve{

#ifdef MYDEBUG
  class HermiteSplineTest;
#endif //MYDEBUG
class HermiteSpline : public Curve
{

public:
  HermiteSpline() : m_x(2), m_y(2) 
  {}
	
	//! Get the value of the function
  double operator()( double x )
  {
    ULong idx = UpperBoundIdx( x );
    if ( idx == 0 )
    {
      return m_y[0];
    }
    if ( idx == m_x.Size() )
    {
      return m_y[ m_x.Size() - 1 ];
    }
    double t =  ( x - m_x[idx-1] ) / ( m_x[idx] - m_x[idx-1] );
    
    Vect<double>& slope = m_coeff;
    double d = 
       ( 1.0 - 3.0 * t * t + 2.0 * t * t * t ) * m_y[idx-1]
     + ( t - 2.0 * t * t + t * t * t ) * slope[idx-1]
     + ( 3.0 * t * t - 2.0 * t * t * t ) * m_y[idx]
     + ( t * t * t - t * t ) * slope[idx];
     return d;
  }
  
  //! First derivative
  double Deriv1( double x )
  {
    ULong idx = UpperBoundIdx( x );
    if ( idx == 0 )
    {
      return m_y[0];
    }
    if ( idx == m_x.Size() )
    {
      return m_y[ m_x.Size() - 1 ];
    }
    double t =  ( x - m_x[idx-1] ) / ( m_x[idx] - m_x[idx-1] );
    
    Vect<double>& slope = m_coeff;
    double d = 
       ( - 6.0 * t + 6.0 * t * t ) * m_y[idx-1]
     + ( 1.0 - 4.0 * t + 3.0 * t * t ) * slope[idx-1]
     + ( 6.0 * t - 6.0 * t * t ) * m_y[idx]
     + ( 3.0 * t * t - 2.0 * t ) * slope[idx];
     return d;
  }
  
	//! Will fit a cubic Hermite spline to the data. Hermite splines pass through all the points.
  void Fit( Vect<double>& x, Vect<double>& y )
  {
    m_x = x;
    m_y = y;
    
    if ( m_x.Size() < 2 || m_x.Size() != m_y.Size() )
    {
      return;
    }
    Vect<double>& slope = m_coeff;
    slope.Resize( m_x.Size() );
    
    slope[0] = ( m_y[1] - m_y[0] ) / ( m_x[1] - m_x[0] );
    
    ULong size = m_x.Size();
    slope[ size - 1 ] = ( m_y[ size - 1 ] - m_y[ size - 2 ] ) / ( m_x[ size - 1 ] - m_x[ size - 2 ]  );
    for ( ULong i = 1; i < size - 1; ++i )
    {
      slope[i] = ( m_y[i + 1] - m_y[i - 1] ) / ( m_x[i + 1] - m_x[i - 1] );
    }
  }
protected:

	//! Locates the upper bound of the index that is greater than val in the array of y points
  ULong UpperBoundIdx( const double val )
  {
    ULong beg = 0;
    ULong end = m_x.Size() - 1;
    
    if ( val <= m_x[beg] )
    {
      return beg;
    }
    if ( val > m_x[end] )
    {
      return m_x.Size();
    }
    
    while ( ( end - beg ) > 2 )
    {
      ULong pos = ( end + beg ) >> 1; // >> 1 is division by 2
      if ( m_x[pos] > val )
      {
        end = pos;
      }
      else
      {
        beg = pos;
      }
    }

    for ( ULong i = beg; i <= end; ++i )
    {
      if ( m_x[i] >= val )
      {
        return i;
      }   
    }
  }
protected:
  Vect<double> m_x;
  Vect<double> m_y;
private:


#ifdef MYDEBUG
  friend class HermiteSplineTest;
#endif //MYDEBUG
};

}/*curve*/ }/*math*/ 

#endif /*_MATH_CURVE_HERMITESPLINE_h */
