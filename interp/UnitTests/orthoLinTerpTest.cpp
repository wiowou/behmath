#include "../src/orthoLinTerp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace math{
  
class OrthoLinTerpTest
{
public:
    OrthoLinTerp<1> lt1;
    OrthoLinTerp<2> lt2;
    OrthoLinTerp<0> lt3;
    std::string testDir;
    
    int Exec()
    {
#ifdef _WIN32
        testDir = "/cygdrive/f/Documents/prog/proj/math/UnitTests";
#else 
        testDir = "/media/wirelessHD/Documents/prog/proj/math/UnitTests";
#endif //WIN32
        int i = 0;
        std::vector<double> r, x, y, z;
        x.push_back(1.2);
        x.push_back(3.6);
        x.push_back(7.0);
        
        r.push_back(-2.5);
        r.push_back(3.5);
        r.push_back(1.75);
        
        math::VecVec dataTable;
        dataTable.resize(2);
        dataTable[0].swap(x);
        dataTable[1].swap(r);
        
        lt1.DataTable(dataTable);
        //for(int j = 0; j < 10000000; ++j)
        {
          i += Interp1(1.8, -1.0);
          i += Interp1(1.5, -1.75);
          i += Interp1(1.0, -2.5);
          i += Interp1(7.1, 1.75);
          i += Interp1(4.5, 3.036765);
        }
        
        x.clear(); y.clear(); z.clear(); r.clear();
        
        x.push_back(1.2);
        x.push_back(3.6);
        x.push_back(7.0);
        

        y.push_back(-6.0);
        y.push_back(2.0);
        y.push_back(4.0);

        r.push_back(-2.5);
        r.push_back(3.5);
        r.push_back(1.75);
        r.push_back(8.0);
        r.push_back(0.0);
        r.push_back(5.0);
        r.push_back(3.0);
        r.push_back(-4.0);
        r.push_back(6.0);
        
        dataTable.resize(3);
        dataTable[0].swap(x);
        dataTable[1].swap(y);
        dataTable[2].swap(r);

        lt2.DataTable(dataTable);

        i += Interp2(1.8, 0.0, 2.0);
        
        lt3.val = 5.0;
        
        double diff = std::abs(lt3(16.0) - 5.0 );
        i += diff > 0.001;
        return i;
    }

    int Interp1( double d, double ans )
    {
      if ( std::abs( lt1(d) - ans ) > 0.0001 )
      {
        std::cout << "failed Interp1 " << d << std::endl;
        return 1;
      }
      return 0;
    }

    int Interp2( double x1, double x2, double ans )
    {
      if ( std::abs( lt2( x1, x2 ) - ans ) > 0.0001 )
      {
        std::cout << "failed Interp2 " << x1 << " " << x2 << std::endl;
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
    math::OrthoLinTerpTest test;
    return test.Exec();
}
