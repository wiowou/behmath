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

#include "domain.h"

namespace math{
namespace model{

Domain::Domain()
{
}
  
unsigned long long Domain::Index(Cell* p)
{
  return p - &cell[0];
}

unsigned long long Domain::Index(Face* p)
{
  return p - &face[0];
}
  
void Domain::Clear()
{
  cell.clear();
  cellBC.clear();
  face.clear();
  faceBC.clear();
  point.clear();
}

void Domain::Swap(Domain& other)
{
  cell.swap(other.cell);
  cellBC.swap(other.cellBC);
  face.swap(other.face);
  faceBC.swap(other.faceBC);
  point.swap(other.point);
}

void Domain::AssignMaterial(std::map<unsigned long long, eng::Material*> &matMap)
{
  for (unsigned long long i = 0; i < cell.size(); ++i)
  {
    unsigned long long pID = cell[i].pID;
    if (matMap.find(pID) == matMap.end() )
    {
      cell[i].mat = matMap.begin()->second;
    }
    else cell[i].mat = matMap[pID];
  }
}

}/*model*/ }/*math*/ 
