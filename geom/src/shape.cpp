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

Vect Shape::PerpTri()
{
  return Cross( *m_pts[1] - *m_pts[0], *m_pts[2] - *m_pts[0] );
}

double Shape::AreaTri()
{
  return 0.5 * PerpTri().Mag();
}

Vect Shape::PerpQuad()
{
  return Cross( *m_pts[2] - *m_pts[0], *m_pts[3] - *m_pts[1] );
}

double Shape::AreaQuad()
{
  return 0.5 * PerpQuad().Mag();
}

double Shape::VolHex()
{
  double vol = 0.0;
  Point** pts = m_pts;
                                  Vect r1c2 = *pts[6] - *pts[3];  Vect r1c3 = *pts[2] - *pts[0];
  Vect r2c1 = *pts[7] - *pts[0];                                  Vect r2c3 = *pts[6] - *pts[4];
  Vect r3c1 = *pts[6] - *pts[1];  Vect r3c2 = *pts[5] - *pts[0];
  
  vol = TripleProd(  r2c1 + r3c1, r1c2,        r1c3 );
  vol += TripleProd( r2c1,        r1c2 + r3c2, r2c3 );
  vol += TripleProd( r3c1,        r3c2,        r1c3 + r2c3 );
  return vol * 0.083333333333333333;
}

double Shape::VolPrism()
{
  double vol = 0.0;
  Point** pts = m_pts;
                                  Vect r1c2 = *pts[5] - *pts[2];  Vect r1c3 = *pts[2] - *pts[0];
  Vect r2c1 = *pts[5] - *pts[0];                                  Vect r2c3 = *pts[5] - *pts[3];
  Vect r3c1 = *pts[5] - *pts[1];  Vect r3c2 = *pts[4] - *pts[0];
  
  vol = TripleProd(  r2c1 + r3c1, r1c2,        r1c3 );
  vol += TripleProd( r2c1,        r1c2 + r3c2, r2c3 );
  vol += TripleProd( r3c1,        r3c2,        r1c3 + r2c3 );
  return vol * 0.083333333333333333;
}

double Shape::VolPyra()
{
  Vect perp = PerpQuad();
  double base = 0.5 * perp.Mag();
  Vect a( *m_pts[4], *m_pts[0] );
  double height = ScalarProj( a, perp);
  return 0.3333333333333333333 * base * height;
}

double Shape::VolTet()
{
  double vol = 0.0;
  Point** pts = m_pts;
  vol = TripleProd( *pts[1] - *pts[0], *pts[2] - *pts[0], *pts[3] - *pts[0] );
  return vol * 0.166666666666666667;
}

}/*geom*/ }/*math*/ 
