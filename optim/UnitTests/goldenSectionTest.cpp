#include "../src/goldenSection.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include <cmath> //std::abs

namespace math{
namespace optim{

struct FBase
{
  virtual double operator()( double x )
  {
    return 0.0;
  }
};

class Fa : public FBase
{
public:
  double operator()( double x )
  {
    return x * x * x * x - 14.0 * x * x * x + 60.0 * x * x - 70.0 * x;
  }
};

class Fb : public FBase
{
public:
  double operator()( double x )
  {
    return 0.5 * x * x - sin(x);
  }
};

class Fc : public FBase
{
public:
  double operator()( double x )
  {
    return x * x + 4.0 * cos(x);
  }
};

class Fd : public FBase
{
public:
  double operator()( double x )
  {
    double sum = 0.0;
    
    for ( double inc = 1e-6; inc < 100.1; inc *= 100.0 )
    {
      double f0 = f(inc, x);
      double f1;
      for ( double t = 2.0 * inc; t < 100.1 * inc; t += inc )
      {
        f1 = f(t,x);
        sum += 0.5 * ( f0 + f1 ) * inc;
        f0 = f1;
      }      
    }
    return sum;
  }
  
  double f( double t, double x )
  {
    return pow( t, x - 1.0 ) * exp( -t );
  }
};

class GoldenSectionTest
{
public:
  
  std::string testDir;
  
  int Exec()
  {
    GoldenSection<FBase> gs1;
    int i = 0;
    
    Fa fa;
    Fb fb;
    Fc fc;
    Fd fd;
    
    double minx;
    double low = 0.0;
    double high = 3.0;
    std::cout << "minimums in interval 0 - 3" << std::endl;
    
    std::cout << "function a: " << std::endl;
    gs1.Clear();
    gs1.Function(fa);
    gs1.Range( low, high );
    gs1.MaxIter(100);
    minx = gs1.MinX();
    std::cout << minx << std::endl;
    
    i += std::abs( minx - .780886 ) > 1e-4;
    
    std::cout << "function b: " << std::endl;
    gs1.Clear();
    gs1.Function(fb);
    gs1.Range( low, high );
    gs1.MaxIter(100);
    minx = gs1.MinX();
    std::cout << minx << std::endl;
    
    i += std::abs( minx - .739086 ) > 1e-4;
    
    std::cout << "function c: " << std::endl;
    gs1.Clear();
    gs1.Function(fc);
    gs1.Range( low, high );
    gs1.MaxIter(100);
    minx = gs1.MinX();
    std::cout << minx << std::endl;
    
    i += std::abs( minx - 1.89549 ) > 1e-4;
    
    std::cout << "function d: " << std::endl;
    gs1.Clear();
    gs1.Function(fd);
    gs1.Range( low, high );
    gs1.MaxIter(100);
    minx = gs1.MinX();
    std::cout << minx << std::endl;
    
    i += std::abs( minx - 1.49611 ) > 1e-4;
    
    return i;
  }

};

}/*optim*/ }/*math*/ 

int main()
{
  math::optim::GoldenSectionTest test;
  return test.Exec();
}
