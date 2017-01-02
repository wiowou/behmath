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

#ifndef _MATH_GEOM_NODE_h
#define _MATH_GEOM_NODE_h

#include "point.h"
#include "UID.h"

namespace math{
namespace geom{

#ifdef MYDEBUG
  class NodeTest;
#endif //MYDEBUG

class Node : public Point
{

public:
  Node() : Point() {}
  Node(Point& p) : Point(p) {}
  Node(double x, double y, double z) : Point(x,y,z){}
	Node& operator=(Point other)
  {
    X(other.X());
    Y(other.Y());
    Z(other.Z());
    return *this;
  }
  unsigned long long ID() const
	{
		return m_id.ID();
	}
	void ID( unsigned long long id )
	{
		m_id.ID(id);
	}
  static void Compress()
  {
    UID<Node>::Compress();
  }
  static unsigned long long Count()
  {
    return UID<Node>::Count();
  }
  friend bool operator<(const Node &lhs, const Node &rhs)
  {
    return lhs.m_id < rhs.m_id;
  }
protected:

private:
  UID<Node> m_id;

#ifdef MYDEBUG
  friend class NodeTest;
#endif //MYDEBUG
};

}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_NODE_h */
