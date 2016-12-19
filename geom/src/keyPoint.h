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

#ifndef _MATH_GEOM_KEYPOINT_h
#define _MATH_GEOM_KEYPOINT_h

#include <set>

#include "point.h"
#include "UID.h"


namespace math{
namespace geom{

class HermiteSpline;

#ifdef MYDEBUG
  class KeyPointTest;
#endif //MYDEBUG

class KeyPoint : public Point
{

public:
	unsigned long long ID() const
	{
		return m_id.ID();
	}
	
	void ID( unsigned long long id )
	{
		m_id.ID(id);
	}
	
	void Associate(HermiteSpline* p)
	{
		m_spline.insert(p);
	}
	
	void Disassociate(HermiteSpline* p)
	{
		m_spline.erase(p);
	}
	
protected:

private:
	UID<KeyPoint> m_id;
	std::set<HermiteSpline*> m_spline;

#ifdef MYDEBUG
  friend class KeyPointTest;
#endif //MYDEBUG
};

}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_KEYPOINT_h */
