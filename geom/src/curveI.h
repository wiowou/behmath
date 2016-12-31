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

#ifndef _MATH_GEOM_CURVEI_h
#define _MATH_GEOM_CURVEI_h

#include <vector>

#include "point.h"
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
  virtual void Mesh() {};
  virtual void MeshRatio(std::vector<double> *meshRatio) = 0;
  virtual void UnMesh(bool below = false) = 0;
	virtual bool Empty() = 0;
  virtual KeyPoint* Endpoint(int i, bool reverse = false) = 0;
  virtual void Endpoint(KeyPoint* start, KeyPoint* end) = 0;
	virtual double Length() = 0;
  virtual void Associate(SurfaceI* p) = 0;
	virtual void Disassociate(SurfaceI* p) = 0;
  virtual Node* GetNode(unsigned long long i, bool reverse = false) = 0;
  virtual double GetMeshRatio(unsigned long long i, bool reverse = false) = 0;
  virtual unsigned long long NumNode() = 0;
  //! Returns a point on the curve that splits curve into curve of length t and 1-t
	virtual void PointWithRatio(double ratio, Point &point) {} //ratio ranges from 0 to 1.
	virtual void PointWithRatio(std::vector<double> &ratio, std::vector<Point> &point) {};
protected:

private:


#ifdef MYDEBUG
  friend class CurveTest;
#endif //MYDEBUG
};

}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_CURVEI_h */
