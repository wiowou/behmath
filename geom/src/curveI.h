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

#include "surfaceI.h"

namespace math{
namespace geom{

class KeyPoint;
class Node;

#ifdef MYDEBUG
  class CurveTest;
#endif //MYDEBUG

class CurveI
{

public:
  virtual void Clear() = 0;
  virtual bool IsMeshed() = 0;
  virtual void Mesh() = 0;
  //virtual void MeshRatio(std::vector<double> *meshRatio) = 0;
  virtual void UnMesh(bool below = false) = 0;
	virtual bool Empty() = 0;
  virtual KeyPoint* Endpoint(int i) = 0;
	virtual double Length() = 0;
  virtual void Associate(SurfaceI* p) = 0;
	virtual void Disassociate(SurfaceI* p) = 0;
  virtual Node* GetNode(unsigned long long i) = 0;
  virtual double GetMeshRatio(unsigned long long i) = 0;
  virtual unsigned long long NumNode() = 0;
protected:

private:


#ifdef MYDEBUG
  friend class CurveTest;
#endif //MYDEBUG
};

}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_CURVE_h */
