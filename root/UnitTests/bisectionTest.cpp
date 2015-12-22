#include "../src/bisection.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <cmath> //std::abs

#include "../src/newton.h"
#include "../src/secant.h"

namespace math{
namespace root{

class F2
{
public:
  double operator()( double x )
  {
    return x * x - 5.0;
  }
};

class Neutron
{
public:
  double operator()( double x )
  {
    return ( x * x - 1.0 ) / ( 2 * x ) - 1.0 / tan(x);
  }
};

class Beam
{
public:
  double operator()( double x )
  {
    return tan(x) * tanh(x) + 1.0;
  }
};

class BisectionTest
{
public:
  
  
  int Exec()
  {
    int j = 0;
    
    Bisection<F2> rf;
    rf.Range( 0.0, 5.0 );
    rf.Tolerance( 1e-4 );
    rf.MaxIter( 100 );
    
    std::cout << "Bisection " << std::endl;
    std::cout << "root: " << rf.Root() << std::endl;
    std::cout << "Iterations: " << rf.Iter() << std::endl;
    
    j += std::abs( rf.Root() - 2.23607 ) > 1e-3;
    j += rf.Iter() > 17;

    Bisection<F2,Secant> rf_hy;
    rf_hy.Range( 0.0, 5.0 );
    rf_hy.Tolerance( 1e-4 );
    rf_hy.MaxIter( 100 );
    
    std::cout << "Bisection/Secant hybrid" << std::endl;
    std::cout << "root: " << rf_hy.Root() << std::endl;
    std::cout << "Iterations: " << rf_hy.Iter() << std::endl;
    
    j += std::abs( rf_hy.Root() - 2.23607 ) > 1e-3;
    j += rf_hy.Iter() > 4;
    
    Bisection<Neutron> rf2;
    rf2.Tolerance( 1e-4 );
    rf2.MaxIter( 100 );

    Bisection<Neutron,Secant> rf2_hy;
    rf2_hy.Tolerance( 1e-4 );
    rf2_hy.MaxIter( 100 );    
    
    double pi = 3.14159265;
    for ( int i = 0; i < 20; ++i )
    {
      rf2.Range( i * pi + 1e-3 , ( i + 1 ) * pi - 1e-3 );
      rf2.Root();
      rf2_hy.Range( i * pi + 1e-3 , ( i + 1 ) * pi - 1e-3 );
      rf2_hy.Root();
      if ( rf2.RootFound() && rf2_hy.RootFound() )
      {
        std::cout << "root: " << rf2.Root() << "  ";
        std::cout << "Iterations: " << rf2.Iter() << "  Hybrid Iterations: " << rf2_hy.Iter() << std::endl;        
      }
      else if ( rf2.RootFound() )
      {
        std::cout << "root: " << rf2.Root() << "  ";
        std::cout << "Iterations: " << rf2.Iter() << "  Hybrid Iterations: unconverged" << std::endl;         
      }
      else if ( rf2_hy.RootFound() )
      {
        std::cout << "root: " << rf2.Root() << "  ";
        std::cout << "Iterations: unconverged" << "  Hybrid Iterations: " << rf2_hy.Iter() << std::endl;        
      }      
    }

    std::cout << "Beam " << std::endl;
    
    Bisection<Beam> rf3;
    rf3.Tolerance( 1e-4 );
    rf3.MaxIter( 100 );
    
    for ( int i = 0; i < 20; ++i )
    {
      rf3.Range( i * pi / 2.0 + 1e-3 , ( i + 1 ) * pi / 2.0 - 1e-3 );
      rf3.Root();
      if ( rf3.RootFound() )
      {
        std::cout << "root: " << rf3.Root() << "  ";
        std::cout << "Iterations: " << rf3.Iter() << std::endl;        
      } 
    }    
    
    return j;
  }

};

}/*root*/ }/*math*/ 

int main()
{
  math::root::BisectionTest test;
  return test.Exec();
}
