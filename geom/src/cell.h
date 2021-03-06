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

#ifndef _MATH_GEOM_CELL_h
#define _MATH_GEOM_CELL_h

#include "point.h"
#include "eng/material.h"
#include "neighbor.h"
#include <vector>

namespace math{
namespace geom{

class Attribute;
#ifdef MYDEBUG
  class CellTest;
#endif //MYDEBUG

class Cell
{

public:
  Cell()
  {
    pID = 0;
    mat = nullptr;
    vol = 0.0;
  }
  unsigned long long pID;
  double vol;
  Point ctr;
  std::vector<Neighbor> nbor;
  std::vector<unsigned long long> pointID;
  eng::Material* mat;
  Attribute* att;
  
protected:

private:


#ifdef MYDEBUG
  friend class CellTest;
#endif //MYDEBUG
};

}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_CELL_h */
