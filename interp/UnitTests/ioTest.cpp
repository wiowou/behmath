#include "../src/io.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace math{

class IOTest
{
public:
    math::IO io;
    std::string testDir;
    
    int Exec()
    {
#ifdef _WIN32
        testDir = "/cygdrive/f/Documents/prog/proj/${Project}/UnitTests";
#else 
        testDir = "/media/wirelessHD/Documents/prog/proj/${Project}/UnitTests";
#endif //WIN32
        int i = 0;
        i = JSONStringToRelation1D();
        return i;
    }
    
    int JSONStringToRelation1D()
    {
      std::vector<std::vector<double> > dataTable;
      std::string s;
      s = "[ [ 1.2, 3.6, 7.0], [ -2.5, 3.5, 1.75] ]";
      io.ConvertJSON( s, dataTable);

      if ( std::abs( dataTable[0][0] - 1.2 ) > 0.0001 )
      {
        std::cout << "failed JSONStringToRelation 1" << std::endl;
        return 1;
      }

      if ( std::abs( dataTable[0][2] - 7.0 ) > 0.0001 )
      {
        std::cout << "failed JSONStringToRelation 2" << std::endl;
        return 1;
      }

      size_t ind = 0;
      if ( std::abs( dataTable[1][0] + 2.5 ) > 0.0001 )
      {
        std::cout << "failed JSONStringToRelation 3" << std::endl;
        return 1;
      }
 
      ind = 2;
      if ( std::abs( dataTable[1][2] - 1.75 ) > 0.0001 )
      {
        std::cout << "failed JSONStringToRelation 4" << std::endl;
        return 1;
      }
      
      std::vector<double> dataTable2;
      io.ConvertJSON( s, dataTable2 );
      
      if ( std::abs( dataTable2[0] - 1.01 ) > 0.0001 )
      {
        std::cout << "failed JSONStringToRelation 5" << std::endl;
        return 1;
      }
      
      if ( std::abs( dataTable2[1] - 3.01 ) > 0.0001 )
      {
        std::cout << "failed JSONStringToRelation 6" << std::endl;
        return 1;
      }
      
      if ( std::abs( dataTable2[3] - 3.6 ) > 0.0001 )
      {
        std::cout << "failed JSONStringToRelation 7" << std::endl;
        return 1;
      }
      
      if ( std::abs( dataTable2[7] - 1.75 ) > 0.0001 )
      {
        std::cout << "failed JSONStringToRelation 8" << std::endl;
        return 1;
      }
      
      if ( dataTable2.size() != 8 )
      {
        std::cout << "failed JSONStringToRelation 9" << std::endl;
        return 1;
      }
      
      s = "[ [ 1.0, 2.0, 3.0 ], [ 1.0, 2.0, 3.0, 4.0 ], [ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 ] ]";
      std::vector<double> dataTable3;
      io.ConvertJSON( s, dataTable3 );
      
      if ( dataTable3.size() != 22 )
      {
        std::cout << "failed JSONStringToRelation 9" << std::endl;
        return 1;
      }
      
      if ( std::abs( dataTable3[0] - 2.01 ) > 0.0001 )
      {
        std::cout << "failed JSONStringToRelation 10" << std::endl;
        return 1;
      }
      
      if ( std::abs( dataTable3[1] - 3.01 ) > 0.0001 )
      {
        std::cout << "failed JSONStringToRelation 11" << std::endl;
        return 1;
      }
      
      if ( std::abs( dataTable3[2] - 4.01 ) > 0.0001 )
      {
        std::cout << "failed JSONStringToRelation 12" << std::endl;
        return 1;
      }

      return 0;
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

}/*math*/ 

int main()
{
    math::IOTest test;
    return test.Exec();
}
