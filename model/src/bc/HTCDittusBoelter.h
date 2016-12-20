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

#ifndef _MATH_MODEL_HTCDITTUSBOELTER_h
#define _MATH_MODEL_HTCDITTUSBOELTER_h

#include "boundaryCondition.h"
#include "interp/orthoLinTerp.h"
#include "eng/material.h"
#include "domain.h"


namespace math{
namespace model{
namespace bc{

#ifdef MYDEBUG
  class HTCDittusBoelterTest;
#endif //MYDEBUG

class HTCDittusBoelter : public BoundaryCondition
{

public:
  HTCDittusBoelter()
  {
    TFluidCell = nullptr;
  }
  void Apply( Face* f1);
  OrthoLinTerp<4> TFluid;
  Cell* TFluidCell;
  double velAvg;
  double dhyd;
  eng::Material* fluid;
protected:

private:


#ifdef MYDEBUG
  friend class HTCDittusBoelterTest;
#endif //MYDEBUG
};

}/*bc*/ }/*model*/ }/*math*/ 

#endif /*_MATH_MODEL_BC_HTCDITTUSBOELTER_h */
