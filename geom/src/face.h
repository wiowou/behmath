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

#ifndef _MATH_GEOM_FACE_h
#define _MATH_GEOM_FACE_h

#include "vect.h"

namespace math{
namespace geom{

class Cell;
class Attribute;
#ifdef MYDEBUG
  class FaceTest;
#endif //MYDEBUG

class Face
{

public:
  Face()
  {
    pID = 0;
    for (int i = 0; i < 4; ++i) point[i] = nullptr;
    cell[0] = cell[1] = nullptr;
    area = 0.0;
  }
  
  Point* point[4];
  Point ctr;
  Vect norm;
  double area;
  mutable Cell* cell[2];
  mutable unsigned long long pID;
  mutable Attribute* att;
protected:

private:


#ifdef MYDEBUG
  friend class FaceTest;
#endif //MYDEBUG
};

}/*model*/ }/*math*/ 

#endif /*_MATH_MODEL_FACE_h */
