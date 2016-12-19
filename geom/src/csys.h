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

#ifndef _MATH_GEOM_CSYS_h
#define _MATH_GEOM_CSYS_h

#include "point.h"

namespace math{
namespace geom{

#ifdef MYDEBUG
  class CsysTest;
#endif //MYDEBUG
class Csys
{

public:
  Csys();
	Csys(int i);
  Csys(Point &o, Point &x, Point &y, bool isCyl = false);
  Csys(Point *o, Point *x, Point *y, bool isCyl = false);
  Csys(Point o, Point x, Point y, bool isCyl = false);
  void ToGlobal(Point &p); //converts p to csys 0 while assuming it's in local csys
  void ToGlobal(Point *p); //converts p to csys 0 while assuming it's in local csys
  void ToLocal(Point &p); //converts p to local csys while assuming it's in global csys
  void ToLocal(Point *p); //converts p to local csys while assuming it's in global csys
  Csys& RotateX(double deg); //rotates the csys deg about local x axis
  Csys& RotateY(double deg); //rotates the csys deg about local y axis
  Csys& RotateZ(double deg); //rotates the csys deg about local z axis 
	Csys& OffsetX(double d); //offsets the csys d distance along its x direction
	Csys& OffsetY(double d); //offsets the csys d distance along its y direction
	Csys& OffsetZ(double d); //offsets the csys d distance along its z direction
	
  bool IsCyl();
  void IsCyl(bool icyl);
protected:
  void Initialize(Point &o, Point &x, Point &y, bool isCyl);
  Csys& RX(Point &p); //rotate p in csys about X
  Csys& RY(Point &p); //rotate p in csys about Y
  Csys& RZ(Point &p); //rotate p in csys about Z
  Csys& ARX(Point &p); //anti-rotate p in csys about X
  Csys& ARY(Point &p); //anti-rotate p in csys about Y
  Csys& ARZ(Point &p); //anti-rotate p in csys about Z
private:
  Point m_origin;
  double m_theta[3]; //theta x,y,z in radians
  bool m_isCyl;

#ifdef MYDEBUG
  friend class CsysTest;
#endif //MYDEBUG
};

extern const long double PI;

}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_CSYS_h */
