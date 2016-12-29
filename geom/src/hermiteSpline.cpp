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
#include <math.h>
#include <cmath>

#include <algorithm>

#include "hermiteSpline.h"

namespace math{
namespace geom{

const double quadPt[] = 
{
	1./2. - sqrt(5.0 + 2.0 * sqrt(10.0/7.0))/6.0,
	1./2. - sqrt(5.0 - 2.0 * sqrt(10.0/7.0))/6.0, 
	1./2.,
	1./2. + sqrt(5.0 - 2.0 * sqrt(10.0/7.0))/6.0,
	1./2. + sqrt(5.0 + 2.0 * sqrt(10.0/7.0))/6.0
};

const double quadWt[] = 
{
	(322. - 13. * sqrt(70.)) / 900.0,
	(322. + 13. * sqrt(70.)) / 900.0, 
	128. / 225.,
	(322. + 13. * sqrt(70.)) / 900.0,
	(322. - 13. * sqrt(70.)) / 900.0
};

HermiteSpline::HermiteSpline()
{
	Clear();
}

HermiteSpline::~HermiteSpline()
{
  UnMesh();
}

void HermiteSpline::Clear()
{
	for (unsigned long long i = 0; i < m_point.size(); ++i)
	{
		m_point[i]->Disassociate(this);
	}
	m_length = 0.0;
	m_ratio.clear();
	m_point.clear();
	m_slope.clear();
  
  m_id[0] = m_id[1] = m_endpoint[0] = m_endpoint[1] = nullptr;
  m_meshRatio = nullptr;
}

bool HermiteSpline::IsMeshed()
{
  return !m_node.empty();
}

void HermiteSpline::Mesh()
{
  m_node.clear();
  if (!m_endpoint[0]->IsMeshed()) m_endpoint[0]->Mesh();
  m_node.push_back(m_endpoint[0]->GetNode());
  if (m_meshRatio != nullptr)
  {
    std::vector<Point> point;
    PointWithRatio(*m_meshRatio, point);
    for (Point& p : point)
    {
      Node* n = new Node(p);
      m_node.push_back(n);
    }
  }
  if (!m_endpoint[1]->IsMeshed()) m_endpoint[1]->Mesh();
  m_node.push_back(m_endpoint[1]->GetNode());
}

void HermiteSpline::UnMesh(bool below)
{
  if (below) UnMeshEndPoints();
  for (unsigned long long i = 1; i < m_node.size() - 1; ++i)
  {
    delete m_node[i];
  }
}

void HermiteSpline::UnMeshEndPoints()
{
  //unmesh the keypoints if possible
  for (int i = 0; i < 2; ++i)
  {
    unsigned long long meshCount = 0;
    for (Curve* p : m_endpoint[i]->AssociatedCurve())
    {
      if (p->IsMeshed()) ++meshCount;
    }
    if (meshCount == 0) m_endpoint[i]->UnMesh();
  }  
}

void HermiteSpline::PointWithRatio(std::vector<double> &ratio, std::vector<Point> &point)
{
	point.clear();
	if (ratio.empty() ) return;
	std::vector<double>::iterator it = std::upper_bound(m_ratio.begin(), m_ratio.end(), ratio[0]);
	if (it == m_ratio.begin() )
	{
		point.push_back(*m_point[0]);
	}
	unsigned long long idx = it - m_ratio.begin();
	while (point.size() != ratio.size() )
	{
		if (it == m_ratio.end())
		{
			point.push_back(*m_point.back());
			continue;
		}
		if (ratio[point.size()] >= *it) 
		{
			//++it;
			it = std::upper_bound(it, m_ratio.end(), ratio[point.size()]);
			continue;
		}
		unsigned long long idx = it - m_ratio.begin();
		double targLength = (ratio[point.size()] - m_ratio[idx-1]) * m_length;
		Point p;
		PointAtTargetLengthOnSegment(idx, targLength, p);
		point.push_back(p);
	}
	return;
}

void HermiteSpline::PointWithRatio(double ratio, Point &point)
{
	std::vector<double>::iterator it = std::upper_bound(m_ratio.begin(), m_ratio.end(), ratio);
	if (it == m_ratio.begin() )
	{
		point = *m_point[0];
		return;
	}
	else if (it == m_ratio.end() )
	{
		point = *m_point.back();
		return;
	}
	unsigned long long idx = it - m_ratio.begin();
	double targLength = (ratio - m_ratio[idx-1]) * m_length;
	PointAtTargetLengthOnSegment(idx, targLength, point);
	return;
}

void HermiteSpline::PointAtTargetLengthOnSegment(unsigned long long idx, double targLength, Point &p)
{
	double tlow = 0.0;
	double tup = 1.0;
	double t = 0.0;
	double tol = 0.0001 * targLength;
	double length = 0.0;
	do
	{
		t = 0.5 * (tlow + tup);
		length = LengthOnSegment(idx, t);
		if (length < targLength)
		{
			tlow = t;
		}
		else 
		{
			tup = t;
		}
	} 
	while (std::abs(length - targLength) > tol);
	t = 0.5 * (tlow + tup);
	PointOnSegment(idx, t, p);
	return;
}

bool HermiteSpline::Fit(std::vector<KeyPoint*> pt)
{
	unsigned long long npts = pt.size();
  if (npts < 1) return true;
  m_point = pt;
  if (pt.back() == pt.front() ) return true; //it's a loop
  if (pt.back() < pt.front() )
  {
    m_id[0] = m_endpoint[0] = pt.back();
    m_id[1] = m_endpoint[1] = pt.front();
  }
  else
  {
    m_id[1] = m_endpoint[1] = pt.back();
    m_id[0] = m_endpoint[0] = pt.front();
  }
  std::vector<double> deltat;
  //Vect vlenNorm(m_point.back(), m_point.front());
  
  Vect vfirst(m_point[1], m_point[0]);
  deltat.push_back(vfirst.Mag());
  m_point[0]->Associate(this);
  for (unsigned long long i = 1; i < npts - 1; ++i)
	{
    Vect v(m_point[i+1], m_point[i-1]);
    deltat.push_back(v.Mag());
    m_point[i]->Associate(this);
	}
  Vect vlast(m_point[npts-1], m_point[npts-2]);
  deltat.push_back(vfirst.Mag());
  m_point[npts-1]->Associate(this);
  
	Update(deltat);
  return false;
}

void HermiteSpline::Update(std::vector<double> &deltat)
{
	unsigned long long npts = m_point.size();
	if (npts < 1) return;
	m_ratio.resize(npts);
	m_slope.resize(npts);
	m_ratio[0] = 0.0;
	if (npts == 2)
	{	
		m_ratio[1] = 1.0;
		for (int i = 0; i < 3; ++i)
		{
			m_slope[0][i] = (*m_point[1])[i] - (*m_point[0])[i];
			m_slope[1][i] = m_slope[0][i];
		}
		return;
	}
	//calculate slopes of first and last point
	for (int i = 0; i < 3; ++i)
	{
		m_slope[0][i] = (*m_point[1])[i] - (*m_point[0])[i];
    m_slope[0][i] /= deltat[0];
		m_slope[npts-1][i] = (*m_point[npts-1])[i] - (*m_point[npts-2])[i];
    m_slope[npts-1][i] /= deltat[npts-1];
	}
	//calculate slopes for points in between
	for (unsigned long long n = 1; n < npts - 1; ++n)
	{
		for (int i = 0; i < 3; ++i)
		{
			m_slope[n][i] = ((*m_point[n+1])[i] - (*m_point[n-1])[i]) / deltat[n];
		}
	}
	//calculate segment distances
	for (unsigned long long n = 1; n < npts; ++n)
	{
		m_ratio[n] = m_ratio[n-1] + LengthOnSegment(n, 1.0);
	}
	m_length = m_ratio[npts-1];
	//turning distances into ratios
	for (unsigned long long n = 1; n < npts; ++n)
	{
		m_ratio[n] /= m_length;
	}
	return;
}

double HermiteSpline::LengthOnSegment(unsigned long long idx, double t)
{
	double length = 0.0;
	for (int i = 0; i < 5; ++i)
	{
		Vect v;
		DerivOnSegment(idx, quadPt[i] * t, v);
		length += quadWt[i] * v.Mag();
	}
	length *= t / 2.0;
	return length;
}

void HermiteSpline::PointOnSegment(unsigned long long idx, double t, Point &p)
{
	for (int i = 0; i < 3; ++i)
	{
		p[i] = 
		 ( 1.0 - 3.0 * t * t + 2.0 * t * t * t ) * (*m_point[idx-1])[i]
		+ ( t - 2.0 * t * t + t * t * t ) * m_slope[idx-1][i]
		+ ( 3.0 * t * t - 2.0 * t * t * t ) * (*m_point[idx])[i]
		+ ( t * t * t - t * t ) * m_slope[idx][i];
	}
	return;
}

void HermiteSpline::DerivOnSegment(unsigned long long idx, double t, Vect &v)
{
	for (int i = 0; i < 3; ++i)
	{
		v[i] = 
		 ( - 6.0 * t + 6.0 * t * t ) * (*m_point[idx-1])[i]
		+ ( 1.0 - 4.0 * t + 3.0 * t * t ) * m_slope[idx-1][i]
		+ ( 6.0 * t - 6.0 * t * t ) * (*m_point[idx])[i]
		+ ( 3.0 * t * t - 2.0 * t ) * m_slope[idx][i];
	}
	return;
}

void HermiteSpline::MeshRatio(std::vector<double> *meshRatio)
{
  m_meshRatio = meshRatio;
}

bool HermiteSpline::Empty()
{
	return m_point.size() > 1;
}

KeyPoint** HermiteSpline::Endpoint()
{
  return m_endpoint;
}

double HermiteSpline::Length()
{
	return m_length;
}

void HermiteSpline::Associate(Surface* p)
{
  m_surface.insert(p);
}

void HermiteSpline::Disassociate(Surface* p)
{
  m_surface.erase(p);
}

std::vector<Node*>& HermiteSpline::GetNode()
{
  return m_node;
}

bool operator<(const HermiteSpline &lhs, const HermiteSpline &rhs)
{
  for (int i = 0; i < 2; ++i)
  {
    if (lhs.m_id[i] == rhs.m_id[i]) continue;
    else return lhs.m_id[i] < rhs.m_id[i];
  }
  return false;
}

bool operator==(const HermiteSpline &lhs, const HermiteSpline &rhs)
{
  for (int i = 0; i < 2; ++i)
  {
    if (lhs.m_id[i] != rhs.m_id[i]) return false;
  }
  return true;
}

}/*geom*/ }/*math*/ 
