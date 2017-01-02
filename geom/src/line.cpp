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

#include "line.h"
#include "vect.h"

namespace math{
namespace geom{

Line::Line() : Curve()
{
  m_csys = &g_csys[0];
}

Line::Line(Point* start, Point* end, const Csys* csys)
{
  m_csys = csys;
  m_point[0] = *start;
  m_point[1] = *end;
  if (m_csys->ID() != 0)
  {
    m_csys->ToLocal(m_point[0]);
    m_csys->ToLocal(m_point[1]);
  }
  m_length = Vect(end, start).Mag();
}

void Line::Clear()
{
  m_csys = 0;
  Curve::Clear();
}

void Line::Mesh()
{
  if (m_meshRatio == nullptr) return;
  m_node.clear();
  if (m_endpoint[0] == m_endpoint[1]) //zero length, collapsed line
  {
    for (unsigned long long i = 0; i < m_meshRatio->size() + 2; ++i)
    {
      m_node.push_back(m_endpoint[0]->GetNode());
    }
    return;
  }
  if (!m_endpoint[0]->IsMeshed()) m_endpoint[0]->Mesh();
  m_node.push_back(m_endpoint[0]->GetNode());
  std::vector<Point> point;
  PointWithRatio(*m_meshRatio, point);
  for (Point& p : point)
  {
    Node* n = new Node(p);
    m_node.push_back(n);
  }
  if (!m_endpoint[1]->IsMeshed()) m_endpoint[1]->Mesh();
  m_node.push_back(m_endpoint[1]->GetNode());
}

void Line::Endpoint(KeyPoint* start, KeyPoint* end)
{
  m_point[0] = *start;
  m_point[1] = *end;
  if (m_csys->ID() != 0)
  {
    m_csys->ToLocal(m_point[0]);
    m_csys->ToLocal(m_point[1]);
  }
  m_length = Vect(end, start).Mag();
  Curve::Endpoint(start, end);
}

void Line::SetCsys(Csys* csys)
{
  m_csys = csys;
}

unsigned long long Line::CsysID() const
{
  return m_csys->ID();
}
  
void Line::PointWithRatio(double ratio, Point &point)
{
  if (ratio < 0.0) ratio = 0.0;
  if (ratio > 1.0) ratio = 1.0;
  for (int i = 0; i < 3; ++i)
  {
    point[i] = (1.0 - ratio) * m_point[0][i] + ratio * m_point[1][i];
  }
  m_csys->ToGlobal(point);
}

void Line::PointWithRatio(std::vector<double> &ratio, std::vector<Point> &point)
{
  if (ratio.empty()) return;
  point.resize(ratio.size());
  for (unsigned long long i = 0; i < ratio.size(); ++i)
  {
    PointWithRatio(ratio[i], point[i]);
  }
}

}/*geom*/ }/*math*/ 
