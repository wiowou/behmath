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
  m_length = 0.0;
  m_meshRatio.clear();
  m_curve.clear();
  m_node.clear();
  m_id[0] = m_id[1] = m_endpoint[0] = m_endpoint[1] = nullptr;
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

KeyPoint* MultiCurve::Endpoint(int i)
{
  return m_endpoint[i];
}

double MultiCurve::Length()
{
	return m_length;
}

void MultiCurve::Associate(SurfaceI* p)
{
  m_surface.insert(p);
}

void MultiCurve::Disassociate(SurfaceI* p)
{
  m_surface.erase(p);
}

Node* MultiCurve::GetNode(unsigned long long i)
{
  return m_node[i];
}

double MultiCurve::GetMeshRatio(unsigned long long i)
{
  return m_meshRatio[i];
}

unsigned long long MultiCurve::NumNode()
{
  return m_node.size();
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
  KeyPoint* ep[2];
  ep[0] = m_curve.front()->Endpoint(0);
  ep[1] = m_curve.front()->Endpoint(1);
  if (m_direction.front() ) m_endpoint[0] = ep[0];
  else m_endpoint[0] = ep[1];
  
  ep[0] = m_curve.back()->Endpoint(0);
  ep[1] = m_curve.back()->Endpoint(1);
  if (m_direction.back() ) m_endpoint[1] = ep[0];
  else m_endpoint[1] = ep[1];
  
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
}

void MultiCurve::UpdateNodes()
{
  m_node.clear();
  m_meshRatio.clear();
  unsigned long long i = 0;
  double prevSegLength = 0.0;
  for (CurveI* cv : m_curve)
  {
    if (m_direction[i])
    {
      if (i != 0) 
      {
        m_node.push_back(cv->Endpoint(0)->GetNode());
        m_meshRatio.push_back( (prevSegLength) / m_length);
      }
      for (unsigned long long j = 0; j < cv->NumNode(); ++j)
      {
        m_node.push_back(cv->GetNode(j));
        double lenOnSeg = cv->Length() * cv->GetMeshRatio(j);
        m_meshRatio.push_back( (prevSegLength + lenOnSeg) / m_length);
      }
      if (i != m_curve.size()-1) 
      {
        m_node.push_back(cv->Endpoint(1)->GetNode());
      }
    }
    else
    {
      if (i != 0) 
      {
        m_node.push_back(cv->Endpoint(1)->GetNode());
        m_meshRatio.push_back( (prevSegLength) / m_length);
      }
      for (unsigned long long j = cv->NumNode() - 1; j > 0; --j)
      {
        m_node.push_back(cv->GetNode(j));
        double lenOnSeg = cv->Length() * (1.0 - cv->GetMeshRatio(j));
        m_meshRatio.push_back( (prevSegLength + lenOnSeg) / m_length);
      }
      m_node.push_back(cv->GetNode(0));
      double lenOnSeg = cv->Length() * (1.0 - cv->GetMeshRatio(0));
      m_meshRatio.push_back( (prevSegLength + lenOnSeg) / m_length);
      if (i != m_curve.size()-1) 
      {
        m_node.push_back(cv->Endpoint(0)->GetNode());
      }
    }
    ++i;
    prevSegLength += cv->Length();
  }
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
