#define MYDEBUG

#include "../src/multiCurve.h"
#include "../src/keyPoint.cpp"
#include "../src/point.cpp"
#include "../src/vect.cpp"
#include "../src/curve.cpp"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace geom{

class MultiCurveTest
{
public:

  int Exec()
  {
    MultiCurve multiCurve;
    double d = multiCurve.Length();
    int i = 0;
    
    return i;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::MultiCurveTest test;
  return test.Exec();
}
