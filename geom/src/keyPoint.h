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
#include "node.h"
#include "UID.h"


namespace math{
namespace geom{

class Curve;

#ifdef MYDEBUG
  class KeyPointTest;
#endif //MYDEBUG

class KeyPoint : public Point
{

public:
  KeyPoint();
  KeyPoint(double x, double y, double z);
  KeyPoint(const KeyPoint &other);
  KeyPoint& operator=(KeyPoint other);
  KeyPoint& operator=(Point other);
  ~KeyPoint();
  
  void Swap(KeyPoint &other);
	unsigned long long ID() const;
	void ID( unsigned long long id );
	void Associate(Curve* p);
	void Disassociate(Curve* p);
  std::set<Curve*> AssociatedCurve();
  bool IsMeshed();
  void Mesh();
  void UnMesh();
  Node* GetNode();
  friend bool operator<(const KeyPoint &lhs, const KeyPoint &rhs)
  {
    return lhs.m_id < rhs.m_id;
  }
	
protected:

private:
	UID<KeyPoint> m_id;
	std::set<Curve*> m_curve;
  Node* m_node;
  
#ifdef MYDEBUG
  friend class KeyPointTest;
#endif //MYDEBUG
};

}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_KEYPOINT_h */
