#include "../src/hermiteSpline.h"
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
  math::geom::HermiteSpline hermiteSpline;
  
  int Exec()
  {
    int i = 0;
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
    return i;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::HermiteSplineTest test;
  return test.Exec();
}
