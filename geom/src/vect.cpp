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

Vect::Vect( const Point* head )
{
  Set(head);
}

Vect::Vect( const Point head )
{
  Set(&head);
}

Vect::Vect( const Point* head, const Point* tail )
{
  Set(head,tail);
}

Vect::Vect( const Point head, const Point tail )
{
  Set(&head,&tail);
}

Vect& Vect::operator=( const Point& rhs )
{
  m_crd[0] = rhs.X();
  m_crd[1] = rhs.Y();
  m_crd[2] = rhs.Z();
}

void Vect::Set( const Point* head, const Point* tail )
{
  operator=( (*head) - (*tail) );
}

double Vect::Mag()
{
  double mag = 0.0;
  for ( int i = 0; i < 3; ++i )
  {
    mag += m_crd[i] * m_crd[i];
  }
  mag = sqrt(mag);
  return mag;
}

Vect Cross( const Vect &a, const Vect &b )
{
  Vect v;
  v.m_crd[0] = a.m_crd[1] * b.m_crd[2] - a.m_crd[2] * b.m_crd[1];
  v.m_crd[1] = a.m_crd[2] * b.m_crd[0] - a.m_crd[0] * b.m_crd[2];
  v.m_crd[2] = a.m_crd[0] * b.m_crd[1] - a.m_crd[1] * b.m_crd[0];
  return v;
}

Vect Cross( Vect &a, Vect &b )
{
  Vect v;
  v.m_crd[0] = a.m_crd[1] * b.m_crd[2] - a.m_crd[2] * b.m_crd[1];
  v.m_crd[1] = a.m_crd[2] * b.m_crd[0] - a.m_crd[0] * b.m_crd[2];
  v.m_crd[2] = a.m_crd[0] * b.m_crd[1] - a.m_crd[1] * b.m_crd[0];
  return v;
}

}/*geom*/ }/*math*/ 
