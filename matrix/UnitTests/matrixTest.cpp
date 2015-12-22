#include "../src/matrix.h"
#include "../src/vect.h"
#include "../src/sparseVect.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace math{

class MatrixTest
{
public:
  
  std::string testDir;
  
  int Exec()
  {
#ifdef _WIN32
    testDir = "/cygdrive/f/Documents/prog/proj/${Project}/UnitTests";
#else 
    testDir = "/media/wirelessHD/Documents/prog/proj/${Project}/UnitTests";
#endif //WIN32
    int i = 0;
    Matrix<SparseVect> matrix(3,5);
    
    std::cout << matrix << std::endl;
    
    matrix.Transpose();
    
    std::cout << matrix << std::endl;
    
    Matrix<Vect> A(5);
    
    double dataA[] = 
    {
      2, 4, 2, 4, 2,
      -1, -2, 6, -4, 1,
      1, -1, 5, 1, 1,
      1, 4, -1, 2, -1,
      5, -2, -1, 2, -3
    };
    A.Data(dataA);
    
    std::cout << A << std::endl;
    
    i += CheckVal(A, 1, 2, 6.0);
    
    A.Transpose();
    i += CheckVal(A, 2, 1, 6.0);
    
    std::cout << A << std::endl;
    
    A.Set(2,1, 6.6);
    i += CheckVal(A, 2, 1, 6.6);
    std::cout << A << std::endl;
    
    A.Transpose();
    A.RowSwap(3,1);
    i += CheckVal(A, 1, 1, 4.0);
    std::cout << A << std::endl;
    
    int maxRow = A.MaxIdxBelow(0);
    
    i += maxRow != 4;
    
    A(0,3) = 5.0;
    maxRow = A.MaxIdxAbove(3);
    i += maxRow != 0;
    
    
    return i;
  }
  
  int CheckVal( Matrix<Vect> &A, int m, int n, double d )
  {
    if ( std::abs( (A(m,n) - d) / d ) > 0.0001 )
    {
      std::cout << "Failed " << m << " " << n << ": " << d << " != " << A(m,n) << std::endl;
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
  math::MatrixTest test;
  return test.Exec();
}
