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

#ifndef _MATH_MODEL_GMSHIO_h
#define _MATH_MODEL_GMSHIO_h

#include "mesher.h"
#include <string>
#include <vector>
#include <set>

namespace math{
namespace model{

#ifdef MYDEBUG
  class GmshIOTest;
#endif //MYDEBUG

class GmshIO : public Mesher
{

public:
  GmshIO()
  {
    el0 = 0;
    el0--;
  }
  void Clear();
  void ReadFVM( std::string fname, bool D3 = true);
  
protected:
  void AllocateStorage( std::string fname, bool D3 = true);
  void ReadNode( std::string& line);
  void ReadTet( std::string& line);
  void ReadHex( std::string& line);
  void ReadPrism( std::string& line);
  void ReadPyra( std::string& line);
  void ReadTriBC( std::string& line);
  void ReadQuadBC( std::string& line);

private:
  

#ifdef MYDEBUG
  friend class GmshIOTest;
#endif //MYDEBUG
};

}/*model*/ }/*math*/ 

#endif /*_MATH_MODEL_GMSHIO_h */
