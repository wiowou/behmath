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

#include "vect.h"
#include <math.h>

namespace math{
namespace geom{

Vect::Vect() : Point()
{}

Vect::Vect( double x, double y, double z ) : Point(x,y,z)
{}

Vect::Vect( const Point* head )
{
  Set(head);
}

Vect::Vect( const Point& head )
{
  Set(head);
}

Vect::Vect( const Point* head, const Point* tail )
{
  Set(head,tail);
}

Vect::Vect( const Point& head, const Point& tail )
{
  Set(head,tail);
}

Vect& Vect::operator+=( const Vect& rhs )
{
  for ( int i = 0; i < 3; ++i )
  {
    m_crd[i] += rhs.m_crd[i];
  }
  return *this;
}

Vect& Vect::operator-=( const Vect& rhs )
{
  for ( int i = 0; i < 3; ++i )
  {
    m_crd[i] -= rhs.m_crd[i];
  }
  return *this;
}

Vect& Vect::operator*=( const double rhs )
{
  for ( int i = 0; i < 3; ++i )
  {
    m_crd[i] *= rhs;
  }
  return *this;
}

Vect& Vect::operator/=( const double rhs )
{
  for ( int i = 0; i < 3; ++i )
  {
    m_crd[i] /= rhs;
  }
  return *this;
}

Vect& Vect::operator=( const Point& rhs )
{
  m_crd[0] = rhs.X();
  m_crd[1] = rhs.Y();
  m_crd[2] = rhs.Z();
  return *this;
}

void Vect::Set( const Point* head, const Point* tail )
{
  operator=( (*head) - (*tail) );
}

void Vect::Set( const Point& head, const Point& tail )
{
  operator=( head - tail );
}

void Vect::Set( const Point* head )
{
  operator=( *head );
}

void Vect::Set( const Point& head )
{
  operator=( head );
}

double Vect::Mag() const
{
  double mag = 0.0;
  for ( int i = 0; i < 3; ++i )
  {
    mag += m_crd[i] * m_crd[i];
  }
  mag = sqrt(mag);
  return mag;
}

double Vect::PerpDistance(Point* p) const
{
  PerpDistance(*p);
}

double Vect::PerpDistance(Point& p) const
{
  Vect v(p);
  Vect bhat = VectorProj(v, *this);
  Vect d(v,bhat);
  return d.Mag();
}
  
Vect Cross( const Vect &a, const Vect &b )
{
  Vect v;
  v.m_crd[0] = a.m_crd[1] * b.m_crd[2] - a.m_crd[2] * b.m_crd[1];
  v.m_crd[1] = a.m_crd[2] * b.m_crd[0] - a.m_crd[0] * b.m_crd[2];
  v.m_crd[2] = a.m_crd[0] * b.m_crd[1] - a.m_crd[1] * b.m_crd[0];
  return v;
}

double Dot( const Vect &a, const Vect &b )
{
  double d = 0.0;
  d = a.m_crd[0] * b.m_crd[0] + a.m_crd[1] * b.m_crd[1] + a.m_crd[2] * b.m_crd[2];
  return d;
}

double TripleProd( const Vect &a, const Vect &b, const Vect &c )
{
  return Dot(a, Cross(b, c));
}

Vect operator+( const Vect &lhs, const Vect &rhs )
{
  Vect v;
  for ( int i = 0; i < 3; ++i )
  {
    v.m_crd[i] = lhs.m_crd[i] + rhs.m_crd[i];
  }
  return v;
}

Vect operator-( const Vect &lhs, const Vect &rhs )
{
  Vect v;
  for ( int i = 0; i < 3; ++i )
  {
    v.m_crd[i] = lhs.m_crd[i] - rhs.m_crd[i];
  }
  return v;
}

Vect operator*( const double lhs, const Vect &rhs )
{
  Vect v;
  for ( int i = 0; i < 3; ++i )
  {
    v.m_crd[i] = lhs * rhs.m_crd[i];
  }
  return v;
}

Vect operator*( const Vect &lhs, const double rhs )
{
  Vect v;
  for ( int i = 0; i < 3; ++i )
  {
    v.m_crd[i] = lhs.m_crd[i] * rhs;
  }
  return v;
}

Vect operator/( const Vect &lhs, const double rhs )
{
  Vect v;
  for ( int i = 0; i < 3; ++i )
  {
    v.m_crd[i] = lhs.m_crd[i] / rhs;
  }
  return v;
}

double ScalarProj( const Vect &a, const Vect &b )
{
  return Dot(a,b) / b.Mag();
}

Vect VectorProj( const Vect &a, const Vect &b )
{
  Vect c = b;
  double scaleFactor = Dot(a,b) / Dot(b,b);
  c *= scaleFactor;
  return c;
}

}/*geom*/ }/*math*/ 
