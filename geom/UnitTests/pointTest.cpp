#include "../src/point.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace geom{

class PointTest
{
public:
  math::geom::Point point;
  std::string testDir;
  
  int Exec()
  {
#ifdef _WIN32
    testDir = "/cygdrive/f/Documents/prog/proj/${Project}/UnitTests";
#else 
    testDir = "/media/wirelessHD/Documents/prog/proj/${Project}/UnitTests";
#endif //WIN32
    int i = 0;
    
    return i;
  }
  
  int FileMatch( std::string outFileName, std::string expFileName )
  {
    std::ifstream out, expected;
    std::string lineOut, lineExpected;
    out.open( outFileName );
    if ( out.fail() )
    {
      std::cout << "Failed open" << outFileName << std::endl;
      return 1;
    }
    expected.open( expFileName );
    if ( expected.fail() )
    {
      std::cout << "Failed open" << expFileName << std::endl;
      return 1;
    }
    while ( std::getline( expected, lineExpected ) )
    {
      std::getline( out, lineOut );
      if ( lineOut != lineExpected )
      {
        std::cout << "Failed " << outFileName << std::endl;
        return 1;
      }
    }
    if ( std::getline( out, lineOut ) )
    {
      std::cout << "Failed truncate " << outFileName << std::endl;
      return 1;
    }
    return 0;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::PointTest test;
  return test.Exec();
}
