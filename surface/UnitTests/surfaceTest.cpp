#include "../src/surface.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath> //std::abs

namespace math{

struct Surf1 : public Surface
{
  double operator()( Vect<double>& x )
  {
    return 16.0 * x[0] * x[0] * x[0] * x[0] + 16.0 * x[1] * x[1] * x[1] * x[1] + x[2] * x[2] * x[2] * x[2];
  }
};

struct Surf2 : public Surface
{
  double operator()( Vect<double>& x )
  {
    return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  }
};

struct Surf3 : public Surface
{
  double operator()( Vect<double>& x )
  {
    return x[0] * x[0] * x[0] - x[1];
  }
};

class SurfaceTest
{
public:
  
  int Exec()
  {
    int i = 0;
    Vect<double> x(3);
    Vect<double> ddx(3);
    
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    
    double res;
    Surf1 s1;
    Surf2 s2;
    Surf3 s3;
    
    res = s1( x );
    i += ( std::abs(res - 353.0) ) > 1e-4;
    
    res = s2( x );
    i += ( std::abs(res - 14.0) ) > 1e-4;
    
    res = s3( x );
    i += ( std::abs(res + 1.0) ) > 1e-4;
    
    s1.Grad(x, ddx);
    i += ( std::abs(ddx[0] - 64.0) ) > 1e-4;
    i += ( std::abs(ddx[1] - 512.0) ) > 1e-4;
    i += ( std::abs(ddx[2] - 108.0) ) > 1e-4;
    std::cout << ddx << std::endl;
    
    s2.Grad(x, ddx);
    i += ( std::abs( ddx[0] - 2.0 ) ) > 1e-4;
    i += ( std::abs( ddx[1] - 4.0 ) ) > 1e-4;
    i += ( std::abs( ddx[2] - 6.0 ) ) > 1e-4;
    std::cout << ddx << std::endl;
    
    s3.Grad( x, ddx );
    i += ( std::abs( ddx[0] - 3.0 ) ) > 1e-4;
    i += ( std::abs( ddx[1] + 1.0 ) ) > 1e-4;
    i += ( std::abs( ddx[2] - 0.0 ) ) > 1e-4;
    std::cout << ddx << std::endl;
    
    Matrix<Vect> H;
    s1.Hessian(x, H);
    i += ( std::abs( H[0][0] - 192.0 ) ) > 1e-4;
    i += ( std::abs( H[1][1] - 768.0 ) ) > 1e-4;
    i += ( std::abs( H[2][2] - 108.0 ) ) > 1e-4;
    i += ( std::abs( H[1][0] ) ) > 1e-4;
    i += ( std::abs( H[0][1] ) ) > 1e-4;
    std::cout << H << std::endl;
    
    return i;
  }

};

}/*math*/ 

int main()
{
  math::SurfaceTest test;
  return test.Exec();
}
