#include "../src/hermiteSpline.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

namespace math{
namespace curve{

class HermiteSplineTest
{
public:
  
  int Exec()
  {
    int i = 0;
    Vect<double> x(9);
    Vect<double> y(9);
    
    for ( int i = 0; i < 9; ++i )
    {
      x[i] = i;
      y[i] = i * i;
    }
    
    math::curve::HermiteSpline hs;
    hs.Fit(x,y);
    
    std::cout << hs( 7.5 ) << std::endl;
    std::cout << hs.Deriv1( 7.5 ) << std::endl;
    
    math::Vect<double> ax(10);
    math::Vect<double> dia(10);
    double pi = 3.14159265358979;
    
    for ( int j = 0; j < 10; j++ )
    {
      double d = pi / 10.0 * static_cast<double>(j);
      ax[j] = d;
      dia[j] = 0.01 * ( cos(d) + 1.5 );
    }
    
    hs.Clear();
    hs.Fit( ax, dia );
    std::cout << hs( pi / 10.0 ) << std::endl;
    
    return i;
  }

};

}/*curve*/ }/*math*/ 

int main()
{
  math::curve::HermiteSplineTest test;
  return test.Exec();
}
