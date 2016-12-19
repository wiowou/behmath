#include "../src/csys.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace math{
namespace geom{

class CsysTest
{
public:
  
  int Exec()
  {
		int i;
    Csys csys6(Point(0,0,0), Point(0,0,-1), Point(0,1,0));
    Csys csys11(Point(2,-5,-3), Point(3,-5,-3), Point(2,-4,-3));
    Point p(1,2,3);
    csys11.ToGlobal(p);
    i = Verify(p,3,-3,0);
		csys11.ToLocal(p);
		i = Verify(p,1,2,3);
		csys11 = 6;
		p[0] = 1.0; p[1] = 30.0; p[2] = 0.0;
		csys11.ToGlobal(p);
		i = Verify(p,0.0,0.5,-0.866025);
		csys11.ToLocal(p);
		i = Verify(p,1.0,30.0,0.0);
		csys11 = 0;
		csys11.RotateZ(90.0);
		p[0] = 1.0; p[1] = 0.0; p[2] = 0.0;
		csys11.ToGlobal(p);
		Verify(p,0,1,0);
    return i;
  }
  
  bool Verify(Point p, double x, double y, double z)
  {
    if ( std::abs(p.X() - x) > 0.0001) return true;
    if ( std::abs(p.Y() - y) > 0.0001) return true;
    if ( std::abs(p.Z() - z) > 0.0001) return true;
    return false;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::CsysTest test;
  return test.Exec();
}
