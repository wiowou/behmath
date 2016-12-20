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

#include "HTCDittusBoelter.h"
#include <cmath>

namespace math{
namespace model{
namespace bc{

void HTCDittusBoelter::Apply( Face* f1)
{
  if (pID != f1->pID ) return;
  unsigned long long fidx = dom->Index(f1);
  unsigned long long cidx = dom->Index(f1->cell[0]);
  
  double tsurf = dom->ts[cidx]; //get surface temp from face's cell
  double tfluid;
  if (TFluidCell == nullptr) 
  {
    tfluid = TFluid(f1->ctr[0], f1->ctr[1], f1->ctr[2], time);
  }
  else
  {
    unsigned long long bcidx = dom->Index(TFluidCell);
    tfluid = dom->ts[bcidx];
  }
  double tfilm = 0.5 * (tsurf + tfluid);
  double Re = velAvg * dhyd / fluid->KVisc(tfilm);
  double Pr = fluid->Pr(tfilm); //this will not be correct for compressible flow
  
  double Nu;
  if (tsurf > tfluid) Nu = 0.023 * pow(Re, 0.8) * pow(Pr, 0.4);
  else Nu = 0.023 * pow(Re, 0.8) * pow(Pr, 0.3);
  double hA = Nu * fluid->cond(tfilm) / dhyd * f1->area;
  dom->Q_a[cidx] -= hA;
  dom->Q_b[cidx] += hA * tfluid;
}

}/*bc*/} /*model*/ }/*math*/ 
