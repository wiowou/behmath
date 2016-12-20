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

#ifndef _MATH_MODEL_MESHFACE_h
#define _MATH_MODEL_MESHFACE_h

#include <set>
#include "face.h"

namespace math{
namespace model{

#ifdef MYDEBUG
  class MeshFaceTest;
#endif //MYDEBUG

class MeshFace : public Face
{

public:
  MeshFace()
  {
    for (int i = 0; i < 4; ++i) point[i] = nullptr;
    for (int i = 0; i < 2; ++i) cell[i] = nullptr;
  }
  
  std::set<unsigned long long> pointID;
  
  friend bool operator<( const MeshFace& lhs, const MeshFace& rhs )
  {
    return lhs.pointID < rhs.pointID;
  }
protected:

private:


#ifdef MYDEBUG
  friend class MeshFaceTest;
#endif //MYDEBUG
};

}/*model*/ }/*math*/ 

#endif /*_MATH_MODEL_MESHFACE_h */
