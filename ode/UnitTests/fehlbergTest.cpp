#include "../src/fehlberg.h"
#include <iostream>
#include <fstream>
#include <string>

#include <cmath>

namespace math{
namespace ode{

class F1 : public Surface
{
  double operator()( Vect<double> &x )
  {
    m_x = x;
    return -2.0 * x[0];
  }
  
  /*
  bool Continue()
  {
    return ( m_x[0] > 10.0 );
  */
  
protected:
  Vect<double> m_x;
};

class FehlbergTest
{
public:
  
  int Exec()
  {
    int i = 0;
    math::ode::Fehlberg<double> fb;
    F1 f;
    fb.F( &f );
    double ic = 0.1;
    fb.InitialCondition(ic);
    fb.Range( 0.0, 0.5 );
    fb.MaxSteps( 1000 );
    //fb.ResultRange( 0.05, 100.0, 0.005 );
    fb.Solve();
    std::cout << fb.Steps() << std::endl;
    std::cout << fb.Time() << std::endl;
    std::cout << fb.Result() << std::endl;
    i += std::abs( fb.Result() - 0.03678 ) > 1e-4;
    return i;
  }

};

}/*ode*/ }/*math*/ 

int main()
{
  math::ode::FehlbergTest test;
  return test.Exec();
}
