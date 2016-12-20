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

#ifndef _MATH_MODEL_BC_HTC_h
#define _MATH_MODEL_BC_HTC_h

#include "boundaryCondition.h"
#include "interp/orthoLinTerp.h"
#include "eng/material.h"
#include "domain.h"

namespace math{
namespace model{
namespace bc{

#ifdef MYDEBUG
  class HTCTest;
#endif //MYDEBUG

template<std::size_t ndim>
class HTC : public BoundaryCondition
{

public:
  HTC()
  {
    TambCell = nullptr;
  }
  
  void Apply( Face* f1)
  {
    if (pID != f1->pID ) return;
    unsigned long long fidx = dom->Index(f1);
    unsigned long long cidx = dom->Index(f1->cell[0]);
    
    double h = htc(f1->ctr[0], f1->ctr[1], f1->ctr[2], time);
    double hA = h * f1->area;
    double coeff_m = dom->A_ts(cidx, cidx);
    coeff_m -= hA;
    dom->A_ts.Set(cidx, cidx, coeff_m);
    
    double tamb;
    if (TambCell == nullptr) 
    {
      tamb = Tamb(f1->ctr[0], f1->ctr[1], f1->ctr[2], time);
      //Positive heat flux is into the volume, on the left side of Ax = b.
      //This positive heat flux becomes negative on the right side of Ax = b.
      dom->Q_b[cidx] -= hA * tamb;
    }
    else
    {
      unsigned long long bcidx = dom->Index(TambCell);
      double coeff_n = dom->A_ts(bcidx, bcidx);
      coeff_n += hA;
      dom->A_ts.Set(bcidx, bcidx, coeff_n);
    }
  }
  
  OrthoLinTerp<ndim> Tamb;
  OrthoLinTerp<ndim> htc;
  Cell* TambCell;
protected:

private:


#ifdef MYDEBUG
  friend class HTCTest;
#endif //MYDEBUG
};

}/*bc*/ }/*model*/ }/*math*/ 

#endif /*_MATH_MODEL_BC_HTC_h */
