#include "../src/point.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace geom{

class PointTest
{
public:
  math::geom::Point point;
  std::string testDir;
  
  int Exec()
  {
    int i = 0;
    Point p1(1.0, 2.0, 3.0);
    Point p2(5.0, 8.5, 9.3);
    
    double d = Dist(p1,p2);
    i += !( (d - 9.89646) > -0.001 && (d - 9.89646) < 0.001);
    return i;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::PointTest test;
  return test.Exec();
}
