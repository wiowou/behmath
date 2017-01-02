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

#ifndef _MATH_GEOM_SURFACE_h
#define _MATH_GEOM_SURFACE_h

#include <set>
#include <vector>

#include "surfaceI.h"
#include "curveI.h"
#include "node.h"
#include "face.h"

namespace math{
namespace geom{

#ifdef MYDEBUG
  class SurfaceTest;
#endif //MYDEBUG

class Surface : public SurfaceI
{

public:
  Surface();
  virtual ~Surface();
  virtual void PointAt(double u, double v) {};
  virtual void Associate(VolumeI* p);
	virtual void Disassociate(VolumeI* p);
  virtual VolumeI* AssociatedVolume(int i);
  virtual std::vector<CurveI*> AssociatedCurve();
  virtual void UnMesh(bool below = true);
  virtual bool IsMeshed();
  virtual void Clear();
  virtual bool Empty();
  
  friend bool operator<(const Surface &lhs, const Surface &rhs);
  friend bool operator==(const Surface &lhs, const Surface &rhs);
  
protected:
  void UnMeshLines(bool below);
  void Initialize();
  std::vector<Node*> m_node;
  std::vector<Face*> m_face;
  //! m_ratioNode and m_ratioFace provide a template of u, v coords 
  //! to apply to intermediate faces of a swept mesh
  std::vector<Node*> m_ratioNode;
  std::vector<Face*> m_ratioFace;
  VolumeI* m_volume[2];
  //! the curves should create a loop when read in order
  std::vector<CurveI*> m_curve;
  std::set<CurveI*> m_id;
  //! m_direction set up so that following the loop, thumb points
  //! towards center of m_volume[0] using right hand rule
  std::vector<bool> m_direction;
  
private:


#ifdef MYDEBUG
  friend class SurfaceTest;
#endif //MYDEBUG
};

extern bool operator<(const Surface &lhs, const Surface &rhs);
extern bool operator==(const Surface &lhs, const Surface &rhs);
  
}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_SURFACE_h */
