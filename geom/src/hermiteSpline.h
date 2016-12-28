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

#include "point.h"
#include "vect.h"

namespace math{
namespace geom{

class KeyPoint;

#ifdef MYDEBUG
  class HermiteSplineTest;
#endif //MYDEBUG

class HermiteSpline
{
public:
	HermiteSpline();
	//! Returns a point on the curve that splits curve into curve of length t and 1-t
	void PointWithRatio(double ratio, Point &point); //ratio ranges from 0 to 1.
	void PointWithRatio(std::vector<double> &ratio, std::vector<Point> &point);
	//! returns true on error, false otherwise
  bool Fit(std::vector<KeyPoint*> pt);
	void Clear();
	//! slope and ratio are recalculated
	void Update(std::vector<double> &deltat);
	bool Empty();
	double Length();
  friend bool operator<(const HermiteSpline &lhs, const HermiteSpline &rhs);
  friend bool operator==(const HermiteSpline &lhs, const HermiteSpline &rhs);
private:
	void PointAtTargetLengthOnSegment(unsigned long long idx, double targLength, Point &p);
	void PointOnSegment(unsigned long long idx, double t, Point &p);
	void DerivOnSegment(unsigned long long idx, double t, Vect &v);
	double LengthOnSegment(unsigned long long idx, double t);
	
	std::vector<KeyPoint*> m_point;
	std::vector<Point> m_slope;
	std::vector<double> m_ratio;
	double m_length;
  
  KeyPoint* m_id[2];
#ifdef MYDEBUG
  friend class HermiteSplineTest;
#endif //MYDEBUG
};

extern bool operator<(const HermiteSpline &lhs, const HermiteSpline &rhs);
extern bool operator==(const HermiteSpline &lhs, const HermiteSpline &rhs);

}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_HERMITESPLINE_h */
