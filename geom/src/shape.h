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

#ifndef _MATH_GEOM_SHAPE_h
#define _MATH_GEOM_SHAPE_h

#include "typedefs.h"
#include "point.h"
#include "vect.h"

namespace math{
namespace geom{

#ifdef MYDEBUG
  class ShapeTest;
#endif //MYDEBUG
class Shape
{

public:
  Shape() : m_pts(nullptr)
  {}
  //! Point the m_pts array to an array of pointers
  void Set( Point** p );
  //! Does the shape have any points?
  bool Empty();
  //! Average of the first size points
  Point Average( ULong size );
  //! Normal vector to the plane
  Vect Perp();
  //! The area of a planar set of points with FEA ordering
  double Area( ULong size );
  //! Area of general triangle
  double AreaTri();
  //! Area of general quadrilateral
  double AreaQuad();
  //! Volume of an irregular hexahedron, assuming points have FEA ordering
  double VolHex();
  //! Volume of an irregular tetrahedron, assuming points have FEA ordering
  double VolTet();
protected:
  Point** m_pts;
  
private:


#ifdef MYDEBUG
  friend class ShapeTest;
#endif //MYDEBUG
};

}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_SHAPE_h */
