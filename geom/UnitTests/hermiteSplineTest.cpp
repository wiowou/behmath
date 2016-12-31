#define MYDEBUG

#include "../src/hermiteSpline.h"
#include "../src/point.cpp"
#include "../src/vect.cpp"
#include "../src/UID.cpp"
#include "../src/keyPoint.cpp"
#include "../src/curve.cpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <math.h>

namespace math{
namespace geom{

class HermiteSplineTest
{
public:
  HermiteSpline hs;
  
  int Exec()
  {
    int i = 0;
    KeyPoint p1(0.0,2.0,0);
    KeyPoint p2(0.4,1.959591794,0);
    KeyPoint p3(0.8,1.833030278,0);
    KeyPoint p4(1.2,1.6,0);
    KeyPoint p5(1.6,1.2,0);
    KeyPoint p6(2.0,0.0,0);
    
    std::vector<Point*> pt;
    pt.push_back(&p1);
    pt.push_back(&p2);
    pt.push_back(&p3);
    pt.push_back(&p4);
    pt.push_back(&p5);
    pt.push_back(&p6);
    
    hs.Endpoint(&p1, &p6);
    hs.Fit(pt);
    
    i += PointWithRatio();
    i += Mesh();
    
    
    return i;
  }
  
  int PointWithRatio()
  {
    int i = 0;
    
    Point p;
    hs.PointWithRatio(0.5, p);
    i += ComparePoints(p, Point(sqrt(2.0),sqrt(2.0),0));
    
    std::vector<double> ratio = {0.0, 0.25, 0.5, 0.75, 1.0};
    std::vector<Point> pv;
    hs.PointWithRatio(ratio, pv);
    std::vector<Point> actual;
    actual.push_back(Point(0.0,2.0,0.0));
    actual.push_back(Point(0.765,1.847,0.0));
    actual.push_back(Point(1.414,1.414,0.0));
    actual.push_back(Point(1.847,0.765,0.0));
    actual.push_back(Point(2.0,0.0,0.0));
    
    for (int j = 0; j < actual.size(); ++j) i += ComparePoints(pv[j], actual[j]);
    KeyPoint p7;
    return i;
  }
  
  int Mesh()
  {
    int i = 0;
    std::vector<double> ratio = {0.25, 0.5, 0.75};
    hs.MeshRatio(&ratio);
    hs.Mesh();
    std::vector<Point> actual;
    actual.push_back(Point(0.0,2.0,0.0));
    actual.push_back(Point(0.765,1.847,0.0));
    actual.push_back(Point(1.414,1.414,0.0));
    actual.push_back(Point(1.847,0.765,0.0));
    actual.push_back(Point(2.0,0.0,0.0));
    
    for (int j = 0; j < actual.size(); ++j) i += ComparePoints(*hs.GetNode(j), actual[j]);
    
    return i;
  }
  
  int ComparePoints(Point val, Point actual)
  {
    double tol = 0.05;
    for (int i = 0; i < 3; ++i)
    {
      if (std::abs((val[i] - actual[i])) > std::abs(tol * actual[i]) + 1e-8)
      {
        return 1;
      }
    }
    return 0;
  }
  
};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::HermiteSplineTest test;
  return test.Exec();
}
