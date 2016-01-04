#include "../src/polygon.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace math{
namespace geom{

class PolygonTest
{
public:
  
  int Exec()
  {
    int i = 0;
    math::geom::Polygon poly;
    Point p1(4.167, 7.563, 15.395);
    Point p2(-5.03, 4.38, 13.09);
    Point p3(-5.158, -1.231, 21.37);
    Point p4(4.04, 1.95, 23.67);
    
    poly.Add(p1);
    poly.Add(p2);
    poly.Add(p3);
    poly.Add(p4);
    
    std::cout << "perpendicular x,y,z" << std::endl;
    std::cout << poly.Perp()[0] << std::endl;
    std::cout << poly.Perp()[1] << std::endl;
    std::cout << poly.Perp()[2] << std::endl;
    double d = poly.Perp().Mag();
    std::cout << "perpendicular magnitude" << std::endl;
    std::cout << d << std::endl;
    
    i += !( d > 100.023 && d < 100.025 );
    
    std::cout << "centroid x,y,z" << std::endl;
    d = poly.Centroid()[0];
    std::cout << d << std::endl;
    i += !( d > -0.49526 && d < -0.49524 );
    
    std::cout << "centroid x,y,z" << std::endl;
    d = poly.Centroid()[1];
    std::cout << d << std::endl;
    i += !( d > 3.1654 && d < 3.1656 );
    
    std::cout << "centroid x,y,z" << std::endl;
    d = poly.Centroid()[2];
    std::cout << d << std::endl;
    i += !( d > 18.3812 && d < 18.3814 );
    
    return i;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::PolygonTest test;
  return test.Exec();
}
