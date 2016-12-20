#define MYDEBUG

#include "../src/solver.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{

class SolverTest
{
public:
  math::Solver<> solver;
  std::string testDir;
  
  int Exec()
  {
    int i = 0;
    
    return i;
  }

};

}/*math*/ 

int main()
{
  math::SolverTest test;
  return test.Exec();
}
