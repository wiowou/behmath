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

#ifndef _MATH_GEOM_NEIGHBOR_h
#define _MATH_GEOM_NEIGHBOR_h

namespace math{
namespace geom{

class Cell;
class Face;

#ifdef MYDEBUG
  class NeighborTest;
#endif //MYDEBUG

class Neighbor
{

public:
  Neighbor()
  {
    cell = nullptr;
    face = nullptr;
  }
  Cell* cell;
  Face* face;
protected:

private:


#ifdef MYDEBUG
  friend class NeighborTest;
#endif //MYDEBUG
};

}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_NEIGHBOR_h */
