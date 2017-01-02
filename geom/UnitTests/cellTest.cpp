#define MYDEBUG

#include "../src/cell.h"
#include "../src/point.cpp"
#include "../src/vect.cpp"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace geom{

class CellTest
{
public:
  Cell cell;
  
  int Exec()
  {
    int i = 0;
    
    return i;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::CellTest test;
  return test.Exec();
}
