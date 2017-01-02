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

#ifndef _MATH_GEOM_HERMITESPLINE_h
#define _MATH_GEOM_HERMITESPLINE_h

#include <vector>
#include <set>

#include "curve.h"
#include "surfaceI.h"
#include "keyPoint.h"
#include "point.h"
#include "node.h"
#include "vect.h"

namespace math{
namespace geom{

#ifdef MYDEBUG
  class HermiteSplineTest;
#endif //MYDEBUG

class HermiteSpline : public Curve
{
public:
	HermiteSpline() = default;
  HermiteSpline(const HermiteSpline &other) = delete;
  HermiteSpline& operator=(HermiteSpline other) = delete;
  void Clear();
  void Mesh();
  //! Returns a point on the curve that splits curve into curve of length t and 1-t
	void PointWithRatio(double ratio, Point &point); //ratio ranges from 0 to 1.
	void PointWithRatio(std::vector<double> &ratio, std::vector<Point> &point);
	//! returns true on error, false otherwise
  bool Fit(std::vector<Point*> pt);
  
private:
  
	void PointAtTargetLengthOnSegment(unsigned long long idx, double targLength, Point &p);
	void PointOnSegment(unsigned long long idx, double t, Point &p);
	void DerivOnSegment(unsigned long long idx, double t, Vect &v);
	double LengthOnSegment(unsigned long long idx, double t);
  //! slope and ratio are recalculated
	void Update(std::vector<double> &deltat);
  
	std::vector<Point*> m_point;
	std::vector<Point> m_slope;
	std::vector<double> m_ratio;
  
#ifdef MYDEBUG
  friend class HermiteSplineTest;
#endif //MYDEBUG
};

}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_HERMITESPLINE_h */
