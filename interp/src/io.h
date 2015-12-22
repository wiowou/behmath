//<license>

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
