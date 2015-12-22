#include "../src/minFinder.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace optim{

class MinFinderTest
{
public:
  math::optim::MinFinder<int> minFinder;
  std::string testDir;
  
  int Exec()
  {
    int i = 0;
    
    return i;
  }
 

};

} /*optim*/ }/*math*/ 

int main()
{
  math::optim::MinFinderTest test;
  return test.Exec();
}
