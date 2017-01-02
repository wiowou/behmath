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

#include "surface.h"

namespace math{
namespace geom{

Surface::Surface()
{
  Initialize();
}

Surface::~Surface()
{
  UnMesh();
}

void Surface::Associate(VolumeI* p)
{
  for (int i = 0; i < 2; ++i)
  {
    if (m_volume[i] == nullptr) m_volume[i] = p;
  }
}

void Surface::Disassociate(VolumeI* p)
{
  for (int i = 0; i < 2; ++i)
  {
    if (m_volume[i] == p) m_volume[i] = nullptr;
  }
}

VolumeI* Surface::AssociatedVolume(int i)
{
  return m_volume[i];
}

std::vector<CurveI*> Surface::AssociatedCurve()
{
  return m_curve;
}

void Surface::UnMesh(bool below)
{
  if (m_node.empty() || m_face.empty()) return;
  if (below)
  {
    UnMeshLines(below);
  }
  for (Node* p : m_node)
  {
    if (p != nullptr) delete p;
  }
  m_node.clear();
  for (Face* p : m_face)
  {
    if (p != nullptr) delete p;
  }
  m_face.clear();
  for (Node* p : m_ratioNode)
  {
    if (p != nullptr) delete p;
  }
  m_ratioNode.clear();
  for (Face* p : m_ratioFace)
  {
    if (p != nullptr) delete p;
  }
  m_ratioFace.clear();
}

bool Surface::IsMeshed()
{
  return !m_node.empty();
}

void Surface::Clear()
{
  Initialize();
  m_curve.clear();
  m_id.clear();
  m_direction.clear();
  UnMesh(true);
}

bool Surface::Empty()
{
  return m_id.empty();
}

void Surface::UnMeshLines(bool below)
{
  for (CurveI* p : m_curve)
  {
    unsigned long long meshCount = 0;
    for (SurfaceI* s : p->AssociatedSurface() )
    {
      if (s != nullptr && s->IsMeshed()) 
      {
        ++meshCount;
      }
    }
    if (meshCount == 0) p->UnMesh(below);
  }
}

void Surface::Initialize()
{
  m_volume[0] = nullptr;
  m_volume[1] = nullptr;
}

bool operator<(const Surface &lhs, const Surface &rhs)
{
  return lhs.m_id < rhs.m_id;
}

bool operator==(const Surface &lhs, const Surface &rhs)
{
  return lhs.m_id == rhs.m_id;
}

}/*geom*/ }/*math*/ 
