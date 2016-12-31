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

#include "curve.h"

namespace math{
namespace geom{

void Curve::Clear()
{
	if (m_endpoint[0] != nullptr)
  {
    m_endpoint[0]->Disassociate(this);
  }
  if (m_endpoint[1] != nullptr)
  {
    m_endpoint[1]->Disassociate(this);
  }
	m_length = 0.0;
  m_surface.clear();
  m_id[0] = m_id[1] = m_endpoint[0] = m_endpoint[1] = nullptr;
  m_meshRatio = nullptr;
  UnMesh(true);
}

bool Curve::IsMeshed()
{
  return !m_node.empty();
}

void Curve::UnMesh(bool below)
{
  if (m_node.empty() ) return;
  if (below) UnMeshEndPoints();
  for (unsigned long long i = 1; i < m_node.size() - 1; ++i)
  {
    if (m_node[i] != nullptr) delete m_node[i];
    m_node[i] = nullptr;
  }
}

void Curve::UnMeshEndPoints()
{
  //unmesh the keypoints if possible
  for (int i = 0; i < 2; ++i)
  {
    if (m_endpoint[i] == nullptr) continue;
    unsigned long long meshCount = 0;
    for (CurveI* p : m_endpoint[i]->AssociatedCurve())
    {
      if (p->IsMeshed()) ++meshCount;
    }
    if (meshCount == 0) m_endpoint[i]->UnMesh();
  }  
}

void Curve::MeshRatio(std::vector<double> *meshRatio)
{
  m_meshRatio = meshRatio;
}

bool Curve::Empty()
{
	return (m_endpoint[0] == nullptr || m_endpoint[1] == nullptr);
}

KeyPoint* Curve::Endpoint(int i, bool reverse)
{
  if (reverse) return m_endpoint[1-i];
  return m_endpoint[i];
}

void Curve::Endpoint(KeyPoint* start, KeyPoint* end)
{
  m_endpoint[0] = start;
  m_endpoint[1] = end;
  if (m_endpoint[0] < m_endpoint[1])
  {
    m_id[0] = m_endpoint[0];
    m_id[1] = m_endpoint[1];
  }
  else
  {
    m_id[0] = m_endpoint[1];
    m_id[1] = m_endpoint[0];
  }
  m_endpoint[0]->Associate(this);
  m_endpoint[1]->Associate(this);
}

double Curve::Length()
{
	return m_length;
}

void Curve::Associate(SurfaceI* p)
{
  m_surface.insert(p);
}

void Curve::Disassociate(SurfaceI* p)
{
  m_surface.erase(p);
}

Node* Curve::GetNode(unsigned long long i, bool reverse)
{
  if (reverse)
  {
    unsigned long long iend = m_node.size() - 1;
    return m_node[iend-i];
  }
  return m_node[i];
}

double Curve::GetMeshRatio(unsigned long long i, bool reverse)
{
  unsigned long long iend = m_meshRatio->size() + 1;
  if (reverse)
  {
    if (i == 0) return 1.0;
    if (i == iend) return 0.0;
    return 1.0 - m_meshRatio->at(iend-i-1);
  }
  if (i == 0) return 0.0;
  if (i == iend) return 1.0;
  return m_meshRatio->at(i-1);
}

unsigned long long Curve::NumNode()
{
  return m_node.size();
}

bool operator<(const Curve &lhs, const Curve &rhs)
{
  for (int i = 0; i < 2; ++i)
  {
    if (lhs.m_id[i] == rhs.m_id[i]) continue;
    else return lhs.m_id[i] < rhs.m_id[i];
  }
  return false;
}

bool operator==(const Curve &lhs, const Curve &rhs)
{
  for (int i = 0; i < 2; ++i)
  {
    if (lhs.m_id[i] != rhs.m_id[i]) return false;
  }
  return true;
}

}/*geom*/ }/*math*/ 
