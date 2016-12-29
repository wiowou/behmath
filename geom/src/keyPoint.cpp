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

#include "keyPoint.h"

namespace math{
namespace geom{

KeyPoint::KeyPoint() : Point(), m_node(nullptr) 
{
  
}

KeyPoint::KeyPoint(double x, double y, double z) : Point(x,y,z), m_node(nullptr)
{
  
}

KeyPoint::KeyPoint(const KeyPoint &other)
{
  m_node = nullptr;
  *this = static_cast<Point>(other);
}

KeyPoint& KeyPoint::operator=(KeyPoint other)
{
  *this = static_cast<Point>(other);
  return *this;
}

KeyPoint& KeyPoint::operator=(Point other)
{
  X(other.X());
  Y(other.Y());
  Z(other.Z());
  return *this;
}

KeyPoint::~KeyPoint()
{
  UnMesh();
}

void KeyPoint::Swap(KeyPoint &other)
{
  unsigned long long tmp_ID = ID();
  Point tmp_p = *this;
  
  ID(other.ID() );
  X(other.X());
  Y(other.Y());
  Z(other.Z());
  
  other.ID(tmp_ID);
  other.X(tmp_p.X());
  other.Y(tmp_p.Y());
  other.Z(tmp_p.Z());
}
unsigned long long KeyPoint::ID() const
{
  return m_id.ID();
}

void KeyPoint::ID( unsigned long long id )
{
  m_id.ID(id);
}

void KeyPoint::Associate(Curve* p)
{
  m_curve.insert(p);
}

void KeyPoint::Disassociate(Curve* p)
{
  m_curve.erase(p);
}

std::set<Curve*> KeyPoint::AssociatedCurve()
{
  return m_curve;
}

bool KeyPoint::IsMeshed()
{
  return m_node != nullptr;
}

void KeyPoint::Mesh()
{
  m_node = new Node();
}

void KeyPoint::UnMesh()
{
  if (IsMeshed()) delete m_node;
  m_node = nullptr;
}

Node* KeyPoint::GetNode()
{
  return m_node;
}

}/*geom*/ }/*math*/ 
