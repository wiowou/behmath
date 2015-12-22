#include "../src/invQuad.h"
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

class InvQuadTest
{
public:
  math::root::InvQuad<F2> rf;
  
  int Exec()
  {
    int i = 0;
    rf.X0( 1.0 );
    rf.X1( 1.1 );
    rf.X2( 1.2 );
    rf.Tolerance( 1e-4 );
    rf.MaxIter( 200 );
    
    std::cout << "root: " << rf.Root() << std::endl;
    std::cout << "Iterations: " << rf.Iter() << std::endl;
    
    i += std::abs( rf.Root() - 2.23606 ) > 1e-4;
    i += rf.Iter() > 5; 
              
    return i;
  }


};

}/*root*/ }/*math*/ 

int main()
{
  math::root::InvQuadTest test;
  return test.Exec();
}
