#include "../src/secant.h"
#include <iostream>
#include <fstream>
#include <string>

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

class SecantTest
{
public:
  math::root::Secant<F2> rf;
  
  int Exec()
  {
    int i = 0;
    rf.X0( 1.0 );
    rf.X1( 1.1 );
    rf.Tolerance( 1e-4 );
    rf.MaxIter( 200 );
    
    std::cout << "root: " << rf.Root() << std::endl;
    std::cout << "Iterations: " << rf.Iter() << std::endl;    
    
    i += std::abs( rf.Root() - 2.23606 ) > 1e-4;
    i += rf.Iter() > 5;  

    Secant<Neutron> rf2;
    rf2.Tolerance( 1e-4 );
    rf2.MaxIter( 100 );
    
    double pi = 3.14159265;
    for ( int i = 0; i < 20; ++i )
    {
      rf2.X0( i * pi + pi / 4.0 );
      rf2.X1( i * pi + pi / 4.0 + 0.01 );
      rf2.Root();
      if ( rf2.RootFound() )
      {
        std::cout << "root: " << rf2.Root() << "  ";
        std::cout << "Iterations: " << rf2.Iter() << std::endl;        
      } 
    }
    
    std::cout << "Beam " << std::endl;
    
    Secant<Beam> rf3;
    rf3.Tolerance( 1e-4 );
    rf3.MaxIter( 100 );
    
    for ( int i = 0; i < 20; ++i )
    {
      rf3.X0( i * pi + pi / 2.5 );
      rf3.X1( i * pi + pi / 3.0 );
      rf3.Root();
      if ( rf3.RootFound() )
      {
        std::cout << "root: " << rf3.Root() << "  ";
        std::cout << "Iterations: " << rf3.Iter() << std::endl;        
      } 
    }
    
    return i;
  }

};

}/*root*/ }/*math*/ 

int main()
{
  math::root::SecantTest test;
  return test.Exec();
}
