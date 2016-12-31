#define MYDEBUG

#include "../src/curve.h"
#include "../src/keyPoint.cpp"
#include "../src/point.cpp"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace geom{

class CurveTest
{
public:

  int Exec()
  {
    Curve curve;
    int i = 0;
    
    return i;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::CurveTest test;
  return test.Exec();
}
