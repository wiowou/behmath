#include "../src/rectangular.h"
#include "../src/cholesky.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace math{
namespace solver{

class RectangularTest
{
public:
  
  
  int Exec()
  {
    int i = 0;
    // specify a Cholesky solver for the rectangular system
    math::solver::Rectangular<Cholesky, Vect> rect;
    
    // create the inputs
    double dataAT[] = 
    {
      1,2,3,4,5,6,
      .1,.1,.2,.4,.7,.8,
      .05,.05,.05,.05,.1,.2,
      .1,.2,.5,.2,.1,.05
    };
    Matrix<Vect> AT(4,6);
    AT.Data( dataAT );
    
    std::cout << "A transpose:" << AT << std::endl;
    
    double dataB[] = 
    {
      21,15.2,17,.5,4,-.2
    };
    Vect<double> b(6);
    b.Data(dataB);
    
    std::cout << "B:" << b << std::endl;
    
    double dataF[] = 
    {
      50,-100,75,63
    };
    Vect<double> f(4);
    f.Data(dataF);
    std::cout << "f:" << f << std::endl;
    
    double dataC[] = 
    {
      2,5,6,4,3,7
    };
    Vect<double> c(6);
    c.Data(dataC);
    std::cout << "Diagonal of C:" << c << std::endl;
    
    Vect<double> x(4);
    
    rect.ATranspose(AT);
    rect.B(b);
    rect.F(f);
    rect.C(c);
    rect.X(x);
    
    rect.Solve();
    
    std::cout << "X:" << x << std::endl;
    
    //check the answer
    i += CheckPair( x, 0, -80.7546 );
    i += CheckPair( x, 1, 1385.6138 );
    i += CheckPair( x, 2, -3612.217 );
    i += CheckPair( x, 3, 314.812 );  
    
    std::cout << " problem 3, hw3: " << std::endl;
    
    // specify a Cholesky solver for the rectangular system
    rect.Clear();
    //specify the new A transpose

    double dataAT2[] = 
    {
      0.1,0.1,0.2,0.4,0.7,0.8,1,0.9,0.7,0.4,0.2,0.1,0.1,0.1,0.05,0.05,0.05,0.05,0.05,0.05,
      0.05,0.05,0.05,0.05,0.1,0.2,0.3,0.35,0.5,0.7,0.9,1,0.9,0.7,0.5,0.2,0.1,0.05,0.05,0.05,
      0.1,0.2,0.5,0.2,0.1,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.1,0.4,0.7,0.4,0.2,0.1
    };
    AT.Resize(3,20);
    AT.Data( dataAT2 );
    
    std::cout << "A Transpose: " << AT << std::endl;

    double dataB2[] = 
    {
      0.0775,0.1125,0.2375,0.1725,0.22,0.2675,0.3525,0.355,0.3825,0.4125,0.4625,0.4875,0.4425,0.3525,0.27,0.24,0.3,0.1725,0.1025,0.0675
    };
    b.Clear();
    b.Resize(20);
    b.Data(dataB2);
    
    std::cout << "B: " << b << std::endl;
    
    x.Clear();
    x.Resize(3);
    
    rect.ATranspose(AT);
    rect.B(b);
    rect.X(x);

    rect.Solve();
    
    std::cout << "X:" << x << std::endl;
    //check the answer
    i += CheckPair( x, 0, 0.2 );
    i += CheckPair( x, 1, 0.45 );
    i += CheckPair( x, 2, 0.35 );
    
    return i;
  }

  //helper method
  int CheckPair( SparseVect<double> &sv, ULong pos, double d )
  {
    if ( std::abs(sv[pos] - d) > 1e-4 )
    {
      std::cout << pos << " bad data " << sv[pos] << "!=" << d << std::endl;
      return 1;
    }
    return 0;
  }
  //helper method
  int CheckPair( Vect<double> &v, ULong pos, double d )
  {
    if ( std::abs(d) < 1e-4 )
    {
      if ( v[pos] > 1e-4 )
      {
        std::cout << pos << " bad data " << v[pos] << "!=" << d << std::endl;
        return 1;
      }
      return 0;
    }
    if ( std::abs(v[pos] - d) > 1e-4 )
    {
      std::cout << pos << " bad data " << v[pos] << "!=" << d << std::endl;
      return 1;
    }
    return 0;
  }
  
};

}/*solver*/ }/*math*/ 

int main()
{
  math::solver::RectangularTest test;
  return test.Exec();
}
