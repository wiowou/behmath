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

#ifndef _MATH_MODEL_DOMAIN_h
#define _MATH_MODEL_DOMAIN_h

#include "geom/point.h"
#include "cell.h"
#include "face.h"
#include <vector>
#include "matrix/vect.h"
#include <string>
#include "matrix/matrix.h"
#include <map>

namespace math{
namespace model{

class BoundaryCondition;

#ifdef MYDEBUG
  class DomainTest;
#endif //MYDEBUG

class Domain
{

public:
  Domain();
  void Clear();
  void Swap( Domain& other);
  unsigned long long Index(Cell* p);
  unsigned long long Index(Face* p);
  
  void AssignMaterial(std::map<unsigned long long, eng::Material*> &matMap);
  virtual void SolveFVM( double endTime = -1.0){}
  virtual void SolveFEM( double endTime = -1.0){}
  
  std::vector<Cell> cell;
  std::vector<Cell*> cellBC; //cells that have a boundary condition on them, like heat gen
  std::vector<Face> face;
  std::vector<Face*> faceBC; //faces that have a boundary condition on them, like a htc
  std::vector<geom::Point> point;
    
  //cell DOFs
  Vect<double> ts;
  //face DOFs
  // y = ax + b
  Vect<double> Q;
  Vect<double> Q_a;
  Vect<double> Q_b;
  Matrix<SparseVect> A_ts; 

protected:

  
  
private:
  Domain& operator=(const Domain& rhs)
  {
    return *this;
  }
  
  Domain(const Domain& dom){}

#ifdef MYDEBUG
  friend class DomainTest;
#endif //MYDEBUG
};

}/*model*/ }/*math*/ 

#endif /*_MATH_MODEL_DOMAIN_h */
