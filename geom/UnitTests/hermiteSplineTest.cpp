#include "../src/hermiteSpline.h"
#include "../src/keyPoint.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

namespace math{
namespace geom{

class HermiteSplineTest
{
public:
  
  
  int Exec()
  {
    /*
		std::vector<double> v {-0.1, 0.2, 0.5, 0.8, 0.8, 1.2, 3.5, 6.2};
    std::vector<double>::iterator low, up;
		low = std::lower_bound(v.begin(), v.end(), -0.1);
		up = std::upper_bound(v.begin(), v.end(), -0.1);
		if (low == v.begin() ) std::cout << "before first low" << std::endl;
		if (up == v.begin() ) std::cout << "before first up" << std::endl;
		if (low == v.end() ) std::cout << "after last low" << std::endl;
		if (up == v.end() ) std::cout << "after last up" << std::endl;
		std::cout << low - v.begin() << std::endl;
		std::cout << up - v.begin() << std::endl;
    */
    int i = 0;
    KeyPoint p1(0.0,2.0,0);
    KeyPoint p2(0.4,1.959591794,0);
    KeyPoint p3(0.8,1.833030278,0);
    KeyPoint p4(1.2,1.6,0);
    KeyPoint p5(1.6,1.2,0);
    KeyPoint p6(2.0,0.0,0);
    
    std::vector<KeyPoint*> pt;
    pt.push_back(&p1);
    pt.push_back(&p2);
    pt.push_back(&p3);
    pt.push_back(&p4);
    pt.push_back(&p5);
    pt.push_back(&p6);
    
    HermiteSpline hs;
    hs.Fit(pt);
    Point p;
    hs.PointWithRatio(0.5, p);
    return i;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::HermiteSplineTest test;
  return test.Exec();
}
