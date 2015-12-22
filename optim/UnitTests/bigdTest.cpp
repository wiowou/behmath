#include "../src/bigd.h"
#include "surface/surface.h"
#include <iostream>
#include <fstream>
#include <string>

#include <cmath> //std::abs

namespace math{
namespace optim{

class F5 : public Surface
{
  double operator()( Vect<double>&x )
  {
    return 2.0 * x[0] * x[0] * x[0] - 4.0 * x[0] * x[0] - 6.0 * x[0] * x[1] * ( x[0] - x[1] - 1.0 );
  }
};

class F6 : public Surface
{
  double operator()( Vect<double>&x )
  {
    double d = -5.0 * x[0] * x[1] * x[3] + 4.0 * x[0] * pow( x[1], 1.4 ) + 0.75 * pow( x[2], 0.6 )
      + 1.0 * std::abs( x[0] * x[3] - 8.4 * x[1] * x[2] * ( 1.0 - x[3] ) * ( 1.0 - x[3] ) );
    return d;
  }
};

class BigdTest
{
public:
  
  
  int Exec()
  {
    int i = 0;
    
    math::optim::Bigd<Surface> bigd;
    F5 f5;
    F6 f6;
    Vect<double>* minx;
    
    std::cout << " function 5: " << std::endl;
    
    Vect<double> low;
    Vect<double> high;
    
    low.Resize(2);
    high.Resize(2);
    
    low[0] = 1.0;
    low[1] = 0.0;
    high[0] = 5.0;
    high[1] = 5.0;
    
    bigd.Function(f5);
    bigd.Range(low, high);
    bigd.MaxIter(1000);
    bigd.Dt(0.05);
    bigd.Tolerance(1e-4);
    minx = bigd.MinX();
    std::cout << "iterations: " << bigd.Iter() << std::endl;
    std::cout << "min point: " << *minx << std::endl;
    
    i += std::abs( (*minx)[0] - 1.8685 ) > 1e-4;
    i += std::abs( (*minx)[1] - 0.4343 ) > 1e-4;
    i += bigd.Iter() > 175;
    
    low.Resize(4);
    high.Resize(4);
    
    low[0] = 0.0;
    low[1] = 0.0;
    low[2] = 0.0;
    low[3] = 0.0;
    
    high[0] = 100.0;
    high[1] = 100.0;
    high[2] = 100.0;
    high[3] = 1.0;
    
    bigd.Function(f6);
    bigd.Range(low, high);
    bigd.MaxIter(10000);
    bigd.Dt(0.5);
    bigd.Tolerance(1e-4);
    minx = bigd.MinX();
    std::cout << "iterations: " << bigd.Iter() << std::endl;
    std::cout << "min point: " << *minx << std::endl;
    
    
    
    return i;
  }

};

}/*optim*/ }/*math*/ 

int main()
{
  math::optim::BigdTest test;
  return test.Exec();
}
