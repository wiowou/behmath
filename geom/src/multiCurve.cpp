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

#include "multiCurve.h"
#include "keyPoint.h"

namespace math{
namespace geom{

void MultiCurve::Clear()
{
  m_cvMeshRatio.clear();
  m_curve.clear();
  m_node.clear();
  m_direction.clear();
  Curve::Clear();
}

bool MultiCurve::IsMeshed()
{
  for (CurveI* cv : m_curve)
  {
    if (cv->IsMeshed() ) return true;
  }
  return false;
}

void MultiCurve::Mesh()
{
  for (CurveI* cv : m_curve)
  {
    if (!cv->IsMeshed() ) cv->Mesh();
  }
}

void MultiCurve::UnMesh(bool below)
{
  for (CurveI* cv : m_curve)
  {
    cv->UnMesh(below);
  }
}

bool MultiCurve::Empty()
{
  return m_curve.empty();
}

void MultiCurve::Define(std::vector<CurveI*>& curve)
{
  if (curve.size() < 2) return; 
  m_curve = curve;
  m_length = 0.0;
  for (CurveI* cv : m_curve)
  {
    m_length += cv->Length();
  }
  CurveDirection();
}

void MultiCurve::UpdateNodes()
{
  m_node.clear();
  m_cvMeshRatio.clear();
  unsigned long long i = 0;
  double prevSegLength = 0.0;
  bool rev = !m_direction[i];
  m_node.push_back(m_curve.front()->GetNode(0, rev));
  //! don't include 0.0 and 1.0 in m_cvMeshRatio
  for (CurveI* cv : m_curve)
  {
    rev = !m_direction[i];
    for (unsigned long long j = 1; j < cv->NumNode(); ++j)
    {
      m_node.push_back(cv->GetNode(j, rev));
      double lenOnSeg = cv->Length() * cv->GetMeshRatio(j, rev);
      m_cvMeshRatio.push_back( (prevSegLength + lenOnSeg) / m_length);
    }
    ++i;
    prevSegLength += cv->Length();
  }
  m_cvMeshRatio.pop_back();
}

void MultiCurve::CurveDirection()
{
  //! set to true for a curve if its second endpoint is shared with
  //! one of the endpoints from the following curve 
  KeyPoint* ep0[2];
  KeyPoint* ep1[2];
  unsigned long long nc = m_curve.size();
  for (unsigned long long i = 0; i < nc - 1; ++i)
  {
    ep0[0] = m_curve[i]->Endpoint(0);
    ep0[1] = m_curve[i]->Endpoint(1);
    ep1[0] = m_curve[i+1]->Endpoint(0);
    ep1[1] = m_curve[i+1]->Endpoint(1);
    if (ep0[1] == ep1[0] || ep0[1] == ep1[1])
    {
      //curve runs in same direction as curve list 
      m_direction.push_back(true);
    }
    else
    {
      //curve runs opposite to curve list
      m_direction.push_back(false);
    }
  }
  ep0[0] = m_curve[nc-1]->Endpoint(0); //last curve
  ep0[1] = m_curve[nc-1]->Endpoint(1); //last curve
  ep1[0] = m_curve[nc-2]->Endpoint(0); //sec to last curve
  ep1[1] = m_curve[nc-2]->Endpoint(1); //sec to last curve
  if (ep0[0] == ep1[0] || ep0[0] == ep1[1])
  {
    //curve runs in same direction as curve list
    m_direction.push_back(true);
  }
  else
  {
    //curve runs opposite to curve list
    m_direction.push_back(false);
  }
}

}/*geom*/ }/*math*/ 
