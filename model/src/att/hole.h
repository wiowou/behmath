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

#ifndef _MATH_MODEL_ATT_HOLE_h
#define _MATH_MODEL_ATT_HOLE_h

#include "attribute.h"

namespace math{
namespace model{
namespace att{

#ifdef MYDEBUG
  class HoleTest;
#endif //MYDEBUG

class Hole : public Attribute
{

public:
  double diameter;
  double cd;
protected:

private:


#ifdef MYDEBUG
  friend class HoleTest;
#endif //MYDEBUG
};

}/*att*/ }/*model*/ }/*math*/ 

#endif /*_MATH_MODEL_ATT_HOLE_h */
