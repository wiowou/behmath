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

#ifndef _MATH_GEOM_CURVE_h
#define _MATH_GEOM_CURVE_h

#include <vector>
#include <set>

#include "curveI.h"
#include "surfaceI.h"
#include "keyPoint.h"
#include "point.h"
#include "node.h"

namespace math{
namespace geom{

class KeyPoint;

#ifdef MYDEBUG
  class CurveTest;
#endif //MYDEBUG

class Curve : public CurveI
{

public:

  virtual void Clear();
  virtual bool IsMeshed();

  virtual void MeshRatio(std::vector<double> *meshRatio);
  virtual void UnMesh(bool below = false);
	virtual bool Empty();
  virtual KeyPoint* Endpoint(int i, bool reverse = false);
  virtual void Endpoint(KeyPoint* start, KeyPoint* end);
	virtual double Length();
  virtual void Associate(SurfaceI* p);
	virtual void Disassociate(SurfaceI* p);
  virtual Node* GetNode(unsigned long long i, bool reverse = false);
  virtual double GetMeshRatio(unsigned long long i, bool reverse = false);
  virtual unsigned long long NumNode();
  
  friend bool operator<(const Curve &lhs, const Curve &rhs);
  friend bool operator==(const Curve &lhs, const Curve &rhs);
  
protected:
  virtual void UnMeshEndPoints();
  
	//! includes endpoints
  std::vector<Node*> m_node;
  //! provides a way to store the desired mesh points along the spline
  std::vector<double> *m_meshRatio;
	double m_length;
  std::set<SurfaceI*> m_surface;
  KeyPoint* m_id[2];
  KeyPoint* m_endpoint[2];

#ifdef MYDEBUG
  friend class CurveTest;
#endif //MYDEBUG

};

extern bool operator<(const Curve &lhs, const Curve &rhs);
extern bool operator==(const Curve &lhs, const Curve &rhs);
  
}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_CURVE_h */
