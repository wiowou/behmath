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

#ifndef _MATH_GEOM_LINE_h
#define _MATH_GEOM_LINE_h

#include "curve.h"
#include "csys.h"

namespace math{
namespace geom{

#ifdef MYDEBUG
  class LineTest;
#endif //MYDEBUG

//! Lines can be of zero length and mesh properly
class Line : public Curve
{

public:
  Line();
  void Clear();
  void Mesh();
  //! keypoints are in global cs
  void Endpoint(KeyPoint* start, KeyPoint* end);
  void SetCsys(Csys* csys);
  Csys* GetCsys();
  //! Returns a point on the curve that splits curve into curve of length t and 1-t
	void PointWithRatio(double ratio, Point &point);//ratio ranges from 0 to 1.
	void PointWithRatio(std::vector<double> &ratio, std::vector<Point> &point);

private:
  Csys* m_csys;
  //! points are in local coordinate system
  Point m_point[2];

#ifdef MYDEBUG
  friend class LineTest;
#endif //MYDEBUG
};

}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_LINE_h */
