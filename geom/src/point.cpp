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

#include "point.h"

namespace math{
namespace geom{

const Point origin(0.,0.,0.);
const Point* porigin = &origin;

Point::Point()
{}

Point::Point( double x, double y, double z ) : Point()
{
  Set(x,y,z);
}

void Point::Set( double x, double y, double z )
{
  m_crd[0] = x;
  m_crd[1] = y;
  m_crd[2] = z;
}

double* Point::Crd()
{
  return m_crd;
}

double& Point::operator[]( int i )
{
  return m_crd[i];
}

double Point::X() const
{
  return m_crd[0];
}
void Point::X( double d )
{
  m_crd[0] = d;
}

double Point::Y() const
{
  return m_crd[1];
}
void Point::Y( double d )
{
  m_crd[1] = d;
}

double Point::Z() const
{
  return m_crd[2];
}
void Point::Z( double d )
{
  m_crd[2] = d;
}

Point& Point::operator-=( const Point &rhs )
{
  for ( int i = 0; i < 3; ++i )
  {
    m_crd[i] -= rhs.m_crd[i];
  }
  return *this;
}

Point& Point::operator+=( const Point &rhs )
{
  for ( int i = 0; i < 3; ++i )
  {
    m_crd[i] += rhs.m_crd[i];
  }
  return *this;
}

Point operator-( const Point &lhs, const Point &rhs )
{
  Point p;
  for ( int i = 0; i < 3; ++i )
  {
    p.m_crd[i] = lhs.m_crd[i] - rhs.m_crd[i];
  }
  return p;
}

Point operator+( const Point &lhs, const Point &rhs )
{
  Point p;
  for ( int i = 0; i < 3; ++i )
  {
    p.m_crd[i] = lhs.m_crd[i] + rhs.m_crd[i];
  }
  return p;
}

}/*geom*/ }/*math*/ 
