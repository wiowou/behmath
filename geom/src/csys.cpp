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

#include "csys.h"
#include <math.h>
#include <cmath>
#include "vect.h"

namespace math{
namespace geom{

const long double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148L;

Csys::Csys()
{
  m_theta[0] = m_theta[1] = m_theta[2] = 0.0;
  m_origin = origin;
  m_isCyl = false;
}

Csys::Csys(int i)
{
	switch(i)
	{
		case 0:
		{
			m_theta[0] = m_theta[1] = m_theta[2] = 0.0;
			m_origin = origin;
			m_isCyl = false;
			break;
		}
		case 1:
		{
			m_theta[0] = m_theta[1] = m_theta[2] = 0.0;
			m_origin = origin;
			m_isCyl = true;
			break;
		}
		case 5:
		{
			Point o(0,0,0);
			Point x(1,0,0);
			Point y(0,0,-1);
			Initialize(o,x,y,true);
			break;
		}
		case 6:
		{
			Point o(0,0,0);
			Point x(0,0,-1);
			Point y(0,1,0);
			Initialize(o,x,y,true);
			break;
		}
	}
}

Csys::Csys(Point &o, Point &x, Point &y, bool isCyl)
{
  Initialize(o,x,y,isCyl);
}

Csys::Csys(Point *o, Point *x, Point *y, bool isCyl)
{
  Initialize(*o,*x,*y,isCyl);
}

Csys::Csys(Point o, Point x, Point y, bool isCyl)
{
  Initialize(o,x,y,isCyl);
}

void Csys::Initialize(Point &o, Point &x, Point &y, bool isCyl)
{
  m_origin = o;
  Vect vx(x,o);
  Vect vy(y,o);
  Vect vz = Cross(vx,vy); //because y won't necessarily be 90 degrees from x
  m_theta[2] = atan2(vx[1],vx[0]);
  m_theta[1] = atan2(-vx[2],vx[0]);
  //m_theta[0] = atan2(vy[2],vy[1]);
	m_theta[0] = atan2(-vz[1],vz[2]);
  m_isCyl = isCyl;
}

bool Csys::IsCyl()
{
  return m_isCyl;
}

void Csys::IsCyl(bool icyl)
{
  m_isCyl = icyl;
}
  
void Csys::ToGlobal(Point &p)
{
  if (m_isCyl)    //double x = p[0] * cos(theta);
  {
		double theta = p[1] * PI / 180.0L;
    p[1] = p[0] * sin(theta);
    p[0] *= cos(theta);
  }
  RX(p);
  RY(p);
  RZ(p);
  p += m_origin;
}

void Csys::ToGlobal(Point *p)
{
  ToGlobal(*p);
}

void Csys::ToLocal(Point &p)
{
  p -= m_origin;
  ARX(p);
  ARY(p);
  ARZ(p);
  if (m_isCyl)
  {
    double theta = atan2(p[1], p[0]) * 180.0L / PI;
    p[0] = std::hypot(p[0], p[1]);
    p[1] = theta;
  }
}

void Csys::ToLocal(Point *p)
{
  ToLocal(*p);
}

Csys& Csys::RotateX(double deg)
{
  double rad = deg * PI / 180.0;
  m_theta[0] += rad;
  return *this;
}

Csys& Csys::RotateY(double deg)
{
  double rad = deg * PI / 180.0;
  m_theta[1] += rad;
  return *this;
}

Csys& Csys::RotateZ(double deg) 
{
  double rad = deg * PI / 180.0;
  m_theta[2] += rad;
  return *this;
}

Csys& Csys::OffsetX(double d)
{
	Point p(d,0,0);
	ToGlobal(p);
	m_origin = p;
	return *this;
}

Csys& Csys::OffsetY(double d)
{
	Point p(0,d,0);
	ToGlobal(p);
	m_origin = p;
	return *this;
}

Csys& Csys::OffsetZ(double d)
{
	Point p(0,0,d);
	ToGlobal(p);
	m_origin = p;
	return *this;
}

Csys& Csys::Offset(Point *p)
{
  m_origin = *p;
  return *this;
}

Csys& Csys::Offset(Point p)
{
  m_origin = p;
  return *this;
}
  
Csys& Csys::RX(Point &p)
{
  double ct = cos(m_theta[0]);
  double st = sin(m_theta[0]);
  double y = ct * p[1] - st * p[2];
  p[2] = st * p[1] + ct * p[2];
  p[1] = y;
  return *this;
}

Csys& Csys::ARX(Point &p)
{
  double ct = cos(m_theta[0]);
  double st = sin(m_theta[0]);
  double y = ct * p[1] + st * p[2];
  p[2] = -st * p[1] + ct * p[2];
  p[1] = y;
  return *this;
}

Csys& Csys::RY(Point &p)
{
  double ct = cos(m_theta[1]);
  double st = sin(m_theta[1]);
  double x = ct * p[0] + st * p[2];
  p[2] = -st * p[0] + ct * p[2];
  p[0] = x;
  return *this;
}

Csys& Csys::ARY(Point &p)
{
  double ct = cos(m_theta[1]);
  double st = sin(m_theta[1]);
  double x = ct * p[0] - st * p[2];
  p[2] = st * p[0] + ct * p[2];
  p[0] = x;
  return *this;
}

Csys& Csys::RZ(Point &p)
{
  double ct = cos(m_theta[2]);
  double st = sin(m_theta[2]);
  double x = ct * p[0] - st * p[1];
  p[1] = st * p[0] + ct * p[1];
  p[0] = x;
  return *this;
}

Csys& Csys::ARZ(Point &p)
{
  double ct = cos(m_theta[2]);
  double st = sin(m_theta[2]);
  double x = ct * p[0] + st * p[1];
  p[1] = -st * p[0] + ct * p[1];
  p[0] = x;
  return *this;
}
  
}/*geom*/ }/*math*/ 
