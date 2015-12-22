#include "../src/solver.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace ode{

class SolverTest
{
public:
  math::ode::Solver<double> solver;
  
  int Exec()
  {
    int i = 0;
    
    return i;
  }

};

}/*ode*/ }/*math*/ 

int main()
{
  math::ode::SolverTest test;
  return test.Exec();
}
