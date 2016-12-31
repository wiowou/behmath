#define MYDEBUG

#include "../src/line.h"
#include "../src/keyPoint.cpp"
#include "../src/point.cpp"
#include "../src/vect.cpp"
#include "../src/csys.cpp"
#include "../src/curve.cpp"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace geom{

class LineTest
{
public:

  int Exec()
  {
    Line line;
    int i = 0;
    
    return i;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::LineTest test;
  return test.Exec();
}
