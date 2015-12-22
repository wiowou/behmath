#include "../src/table.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace math{

class TableTest
{
public:
  math::Table table;
  std::string testDir;
  
  int Exec()
  {
#ifdef _WIN32
    testDir = "/cygdrive/f/Documents/prog/proj/${Project}/UnitTests";
#else 
    testDir = "/media/wirelessHD/Documents/prog/proj/${Project}/UnitTests";
#endif //WIN32
    int i = 0;
    
    
    i += UpperBound( 22.0, 3 );
    i += UpperBound( 2.0, 1 );
    i += UpperBound( 0.5, 0 );
    i += UpperBound( 56.3, -4 );
    i += UpperBound( 36.0, 4 );
    i += UpperBound( 3.5, -1 );
    i += UpperBound( 35.8, -3 );
    
    i += Interp1(1.8, -1.0);
    i += Interp1(1.5, -1.75);
    i += Interp1(1.0, -2.5);
    i += Interp1(4.5, 3.036765);
    
    i += table.Size() != 12;
    
    i += Interp2(2.5, 2.5, 8.5);
    i += Interp2(0.5, 0.5, 1.0);
    i += Interp2(3.0, 2.0, 10.0);
    i += Interp2(0.0, 1.5, 1.5);
    i += Interp2(0.0, 3.5, 3.5);
    i += Interp2(5.5, 1.5, 17.5);
    i += Interp2(5.5, 3.5, 19.5);
    i += Interp2(6.0, 5.0, 20.0);
    i += Interp2(2.5, 0.5, 7.0);
    i += Interp2(4.75, 0.5, 16.0);
    i += Interp2(2.5, 4.5, 10.0);
    i += Interp2(4.75, 4.5, 19.0);
    i += Interp2(4.0, 3.0, 15.0);
    i += Interp2(3.75, 2.75, 13.75);
    i += Interp2(0.5, 3.0, 3.0);
    i += Interp2(0.5, 1.0, 1.0);
    i += Interp2(0.5, 4.0, 4.0);
    i += Interp2(3.0, 1.5, 9.5);
    i += Interp2(3.0, 3.75, 11.75);
    i += Interp2(5.0, 1.5, 17.5);
    i += Interp2(5.0, 3.75, 19.75);
    i += Interp2(5.0, 4.0, 20.0);
    
    i += Interp3(3.0, 4.0, 2.0, 54.0);
    i += Interp3(4.0, 5.0, 4.0, 80.0);
    i += Interp3(1.0, 1.0, 1.0, 1.0);
    i += Interp3(2.5, 3.0, 2.0, 40.0);
    i += Interp3(2.75, 3.0, 2.0, 45.0);
    i += Interp3(2.5, 3.5, 2.0, 42.0);
    i += Interp3(2.5, 3.0, 4.0, 42.0);
    i += Interp3(2.75, 3.0, 4.0, 47.0);
    i += Interp3(2.5, 3.5, 4.0, 44.0);
    i += Interp3(2.5, 3.5, 3.0, 43.0);
    i += Interp3(2.5, 3.5, 2.5, 42.5);
    
    
    return i;
  }
  
  int UpperBound( double val, unsigned long long idxExp )
  {
    double arr[] = { 1.1, 3.5, 21.6, 35.8, 56.2 };
    long long idx;
    idx = table.UpperBoundIdx( arr, arr + 5, val );
    if ( idx != idxExp )
    {
      std::cout << "failed Upper Bound " << val << std::endl;
      return 1;
    }
    return 0;
  }
  
  int Interp1( double val, double valExp )
  {
    double arr[] = { 1.01, 5.01, 
      1.2, 3.6, 7.0, 9.2, 10.0,
     -2.5, 3.5, 1.75, 11.0, -5.25 };
    table.Data(arr);
    double ans = table(val);
    if ( std::abs( (ans - valExp) / valExp ) > 0.001 )
    {
      std::cout << "failed Interp1 " << val << std::endl;
      return 1;
    }
    if ( table.Size() != 12 )
    {
      std::cout << "failed Interp1 Size" << std::endl;
      return 1;
    }
    return 0;
  }
  
  int Interp2( double val1, double val2, double valExp )
  {
    double arr[] = { 2.01, 5.01, 4.01,
      1.0, 2.0, 3.0, 4.0, 5.0,
      1.0, 2.0, 3.0, 4.0,
      1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20. };
    table.Data(arr);
    double ans = table(val1, val2);
    if ( std::abs( (ans - valExp) / valExp ) > 0.001 )
    {
      std::cout << "failed Interp2 " << val1 << " " << val2 << std::endl;
      return 1;
    }
    if ( table.Size() != 32 )
    {
      std::cout << "failed Interp2 Size" << std::endl;
      return 1;
    }
    return 0;
  }
  
  int Interp3( double val1, double val2, double val3, double valExp )
  {
    double arr[] = { 3.01, 4.01, 5.01, 4.01,
      1.0, 2.0, 3.0, 4.0, 
      1.0, 2.0, 3.0, 4.0, 5.0,
      1.0, 2.0, 3.0, 4.0,
      1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 
      21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 
      41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50., 51., 52., 53., 54., 55., 56., 57., 58., 59., 60., 
      61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70., 71., 72., 73., 74., 75., 76., 77., 78., 79., 80.
      };
    table.Data(arr);
    double ans = table(val1, val2, val3);
    if ( std::abs( (ans - valExp) / valExp ) > 0.001 )
    {
      std::cout << "failed Interp3 " << val1 << " " << val2 << " " << val3 << std::endl;
      return 1;
    }
    if ( table.Size() != 97 )
    {
      std::cout << "failed Interp3 Size" << std::endl;
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
  math::TableTest test;
  return test.Exec();
}
