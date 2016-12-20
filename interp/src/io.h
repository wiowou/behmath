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

#ifndef _IO_h
#define _IO_h

#include <vector>
#include "json/json.h"

namespace math{

#ifdef MYDEBUG
    class IOTest;
#endif //MYDEBUG
class IO
{

public:

  //! Initializes Relation with \param jsonArray.
  /*! \param jsonArray is a string containing the array. For example:
      std::string jsonArray = "[ [1.1, 1.5, 3.8], [2.1, 5.4, 6.2] ]";
  */
  bool ConvertJSON( const std::string &jsonStr, std::vector<std::vector<double> > &dataTable )
  {
    json::Value jsonValArray = json::Deserialize( jsonStr );
    return Convert( jsonValArray, dataTable );
  }
  
  bool Convert( const json::Value &jsonValArray, std::vector<std::vector<double> > &dataTable )
  {
    if ( jsonValArray.GetType() == json::NULLVal ) return true;
    std::size_t ndimp1 = jsonValArray.size(); //# of dims plus 1
    dataTable.resize( ndimp1 );
    for ( int i = 0; i < ndimp1; ++i )
    {
      std::size_t dimSize = jsonValArray[i].size();
      dataTable[i].resize( dimSize );
      for ( int j = 0; j < dimSize; ++j )
      {
        dataTable[i][j] = jsonValArray[i][j].ToDouble();
      }
    }
    return false;
  }

  bool ConvertJSON( const std::string &jsonStr, std::vector<double> &dataTable )
  {
    json::Value jsonValArray = json::Deserialize( jsonStr );
    return Convert( jsonValArray, dataTable );
  }
  
  bool Convert( const json::Value &jsonValArray, std::vector<double> &dataTable )
  {
    if ( jsonValArray.GetType() == json::NULLVal ) return true;
    dataTable.clear();
    std::size_t ndimp1 = jsonValArray.size(); //# of dims plus 1
    dataTable.push_back(ndimp1 - 0.99);
    for ( int i = 0; i < ndimp1 - 1; ++i )
    {
      dataTable.push_back( jsonValArray[i].size() + 0.01 );
    }
    for ( int i = 0; i < ndimp1; ++i )
    {
      std::size_t dimSize = jsonValArray[i].size();
      for ( int j = 0; j < dimSize; ++j )
      {
        dataTable.push_back( jsonValArray[i][j].ToDouble() );
      }
    }
    return false;
  }


protected:

private:


#ifdef MYDEBUG
    friend class IOTest;
#endif //MYDEBUG
};

}/*math*/ 

#endif /*_IO_h */
