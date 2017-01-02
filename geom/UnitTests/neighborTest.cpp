#define MYDEBUG

#include "../src/neighbor.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace geom{

class NeighborTest
{
public:
  Neighbor neighbor;
  
  int Exec()
  {
    int i = 0;
    
    return i;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::NeighborTest test;
  return test.Exec();
}
