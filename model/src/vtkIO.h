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

#ifndef _MATH_MODEL_VTKIO_h
#define _MATH_MODEL_VTKIO_h

#include "domain.h"
#include <string>
#include <vector>

namespace math{
namespace model{

#ifdef MYDEBUG
  class VtkIOTest;
#endif //MYDEBUG

class VtkIO
{

public:
  void AddCellData( Vect<double>& data, std::string desc);
  void AddPointData( Vect<double>& data, std::string desc);
  void Clear();
  void Write(std::string fname, std::string title = "3D Unstructured data");
  Domain* dom;
protected:  
  std::vector<Vect<double>*> cellData;
  std::vector<std::string> cellDataDesc;
  std::vector<Vect<double>*> pointData;
  std::vector<std::string> pointDataDesc;
private:


#ifdef MYDEBUG
  friend class VtkIOTest;
#endif //MYDEBUG
};

}/*model*/ }/*math*/ 

#endif /*_MATH_MODEL_VTKIO_h */
