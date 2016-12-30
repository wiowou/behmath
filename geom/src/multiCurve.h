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

#ifndef _MATH_GEOM_MULTICURVE_h
#define _MATH_GEOM_MULTICURVE_h

#include <set>

#include "curveI.h"
#include "node.h"

namespace math{
namespace geom{

#ifdef MYDEBUG
  class MultiCurveTest;
#endif //MYDEBUG

class MultiCurve : public CurveI
{

public:
  void Clear();
  bool IsMeshed();
  void Mesh();
  void UnMesh(bool below = false);
	bool Empty();
  KeyPoint* Endpoint(int i);
	double Length();
  void Associate(SurfaceI* p);
	void Disassociate(SurfaceI* p);
  Node* GetNode(unsigned long long i);
  double GetMeshRatio(unsigned long long i);
  unsigned long long NumNode();
  void Define(std::vector<CurveI*>& curve);
  
private:
  void CurveDirection();
  void UpdateNodes();

  std::vector<CurveI*> m_curve;
  std::vector<double> m_meshRatio;
  //! doesn't include endpoints
  std::vector<Node*> m_node;
  //! set to true for a curve if its second endpoint is shared with
  //! one of the endpoints from the following curve 
  std::vector<bool> m_direction;
  double m_length;
  std::set<SurfaceI*> m_surface;
  KeyPoint* m_id[2];
  KeyPoint* m_endpoint[2];
  
#ifdef MYDEBUG
  friend class MultiCurveTest;
#endif //MYDEBUG
};

}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_MULTICURVE_h */
