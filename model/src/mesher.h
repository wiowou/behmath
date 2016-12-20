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

#ifndef _MATH_MODEL_MESHER_h
#define _MATH_MODEL_MESHER_h

#include "domain.h"
#include "meshFace.h"
#include <set>
#include <vector>

namespace math{
namespace model{

#ifdef MYDEBUG
  class MesherTest;
#endif //MYDEBUG
class Mesher : public Domain
{

public:
  MeshFace MakeMFace(unsigned long long n0, unsigned long long n1, unsigned long long n2);
  MeshFace MakeMFace(unsigned long long n0, unsigned long long n1, unsigned long long n2, unsigned long long n3);
  void CreateCellsFaces();
  
  std::set<MeshFace> mFace;
  unsigned long long el0;
  std::vector<MeshFace> mFaceBC;
protected:

private:


#ifdef MYDEBUG
  friend class MesherTest;
#endif //MYDEBUG
};

}/*model*/ }/*math*/ 

#endif /*_MATH_MODEL_MESHER_h */
