#include "../src/fixedPoint.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <cmath> //std::abs

namespace math{
namespace root{

class F2
{
public:
  double operator()( double x )
  {
    return 1.0 - x * x / 5.0;
  }
};

class Kepler
{
public:
  double operator()( double x )
  {
    return 1.0 + 0.5 * sin(x) - x;
  }
};

class FixedPointTest
{
public:
  math::root::FixedPoint<F2> rf;
  
  int Exec()
  {
    int i = 0;
    rf.X0( 1.0 );
    rf.Tolerance( 1e-4 );
    rf.MaxIter( 100 );
    
    std::cout << "root: " << rf.Root() << std::endl;
    std::cout << "Iterations: " << rf.Iter() << std::endl;
    
    i += std::abs( rf.Root() - 2.23606 ) > 1e-4;
    i += rf.Iter() > 6; 
    
    FixedPoint<Kepler> rf2;
    rf2.X0(1.0);
    rf2.Tolerance(1e-4);
    rf2.MaxIter( 100 );
    
    std::cout << std::endl << " Kepler: " << std::endl;
    std::cout << "root: " << rf2.Root() << std::endl;
    std::cout << "Iterations: " << rf2.Iter() << std::endl; 
    
    i += std::abs( rf2.Root() - 1.4987 ) > 1e-4;
    i += rf2.Iter() > 4;
    
    return i;
  }

};

}/*root*/ }/*math*/ 

int main()
{
  math::root::FixedPointTest test;
  return test.Exec();
}
