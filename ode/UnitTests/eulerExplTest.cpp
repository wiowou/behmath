#include "../src/eulerExpl.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace ode{

class F1 : public Surface
{
  double operator()( Vect<double> &x )
  {
    return -2.0 * x[0];
  }
};

class EulerExplTest
{
public:
  
  
  int Exec()
  {
    int i = 0;
    math::ode::EulerExpl<Surface> ex;
    F1 f;
    ex.F( &f );
    double ic = 0.1;
    ex.InitialCondition(ic);
    ex.Range( 0.0, 0.5 );
    ex.MaxSteps( 1000 );
    //ex.ResultRange( 0.05, 100.0, 0.005 );
    ex.Solve();
    std::cout << ex.Steps() << std::endl;
    std::cout << ex.Time() << std::endl;
    std::cout << ex.Result() << std::endl;
    i += std::abs( ex.Result() - 0.03679 ) > 1e-4;
    return i;
  }

};

}/*ode*/ }/*math*/ 

int main()
{
  math::ode::EulerExplTest test;
  return test.Exec();
}
