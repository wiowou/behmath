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

#include "polygon.h"

namespace math{
namespace geom{

void Polygon::Add( Point* p )
{
  m_pts.push_back(p);
}

void Polygon::Add( Point &p )
{
  m_pts.push_back(&p);
}

void Polygon::Remove( Point* p )
{
  for ( int i = 0; i < m_pts.size(); ++i )
  {
    if ( m_pts[i] == p )
    {
      m_pts.erase( m_pts.begin() + i );
      return;
    }
  }
}

void Polygon::Remove( Point& p )
{
  for ( int i = 0; i < m_pts.size(); ++i )
  {
    if ( m_pts[i] == &p )
    {
      m_pts.erase( m_pts.begin() + i );
      return;
    }
  }
}

Vect Polygon::Perp()
{
  Vect perp(0.,0.,0.);
  
  if ( m_pts.size() == 0 )
  {
    return perp;
  }
  
  m_pts.push_back( m_pts[0] );
  
  for ( int i = 0; i < m_pts.size() - 1; ++i )
  {
    Vect v = Cross( m_pts[i], m_pts[i+1] );
    perp += v;
  }
  for ( int i = 0; i < 3; ++i )
  {
    perp[i] *= 0.5;
  }
  m_pts.pop_back();
  return perp;
}

Point Polygon::Centroid()
{
  Point p(0.,0.,0.);
  for ( int j = 0; j < m_pts.size(); ++j )
  {
    for ( int i = 0; i < 3; ++i )
    {
      p[i] += (*m_pts[j])[i];
    }
  }
  double denom = 1.0 / static_cast<double>( m_pts.size() );
  for ( int i = 0; i < 3; ++i )
  {
    p[i] *= denom;
  }
  return p;
}
  
}/*geom*/ }/*math*/ 
