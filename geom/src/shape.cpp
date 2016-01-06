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

#include "shape.h"

namespace math{
namespace geom{

//! Point the m_pts array to entry 0
void Shape::Set( Point** p )
{
  m_pts = p;
}

bool Shape::Empty()
{
  return m_pts == nullptr;
}

Point Shape::Average( ULong size )
{
  Point p(0.,0.,0.);
  for ( ULong j = 0; j < size; ++j )
  {
    for ( int i = 0; i < 3; ++i )
    {
      p[i] += (*m_pts[j])[i];
    }
  }
  double denom = 1.0 / static_cast<double>( size );
  for ( int i = 0; i < 3; ++i )
  {
    p[i] *= denom;
  }
  return p;
}

Vect Shape::Perp()
{
  return Cross( *m_pts[1] - *m_pts[0], *m_pts[2] - *m_pts[0] );
}

double Shape::Area( ULong size )
{
  Vect perp(0.,0.,0.);
  
  if ( size < 3 )
  {
    return 0.0;
  }
  
  for ( int i = 0; i < size - 1; ++i )
  {
    Vect v = Cross( m_pts[i], m_pts[i+1] );
    perp += v;
  }
  Point* p = m_pts[0];
  Vect v = Cross( m_pts[size-1], m_pts[0] );
  perp += v;
    
  for ( int i = 0; i < 3; ++i )
  {
    perp[i] *= 0.5;
  }

  return perp.Mag();
}

double Shape::AreaTri()
{
  return 0.5 * Perp().Mag();
}

double Shape::VolHex()
{
  double vol = 0.0;
  Point** pts = m_pts;
  vol = TripleProd( *pts[6] - *pts[1] + *pts[7] - *pts[0], *pts[6] - *pts[3], *pts[2] - *pts[0] );
  vol += TripleProd( *pts[7] - *pts[0], *pts[6] - *pts[3] + *pts[5] - *pts[0], *pts[6] - *pts[4] );
  vol += TripleProd( *pts[6] - *pts[1], *pts[5] - *pts[0], *pts[6] - *pts[4] +* pts[2] - *pts[0] );
  return vol * 0.0833333333333333334;
}

double Shape::VolTet()
{
  double vol = 0.0;
  Point** pts = m_pts;
  vol = TripleProd( *pts[1] - *pts[0], *pts[2] - *pts[0], *pts[3] - *pts[0] );
  return vol * 0.166666666666666667;
}

}/*geom*/ }/*math*/ 
