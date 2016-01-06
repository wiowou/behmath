#include "../src/shape.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace math{
namespace geom{

class ShapeTest
{
public:
  math::geom::Shape shape;
  std::string testDir;
  
  int Exec()
  {
    int i = 0;
    Point p1(4.167, 7.563, 15.395);
    Point p2(-5.03, 4.38, 13.09);
    Point p3(-5.158, -1.231, 21.37);
    Point p4(4.04, 1.95, 23.67);
    
    Point* pts[8];
    
    i += !shape.Empty();
    
    shape.Set( &pts[0] );
    
    pts[0] = &p1;
    pts[1] = &p2;
    pts[2] = &p3;
    pts[3] = &p4;
    
    double d = shape.Area(4);
    std::cout << d << std::endl;
    
    i += !( d > 100.023 && d < 100.025 );
    
    Point p = shape.Average(4);
    std::cout << "average x" << std::endl;
    i += !( p[0] > -0.49526 && p[0] < -0.49524 );
    
    std::cout << "average y" << std::endl;
    i += !( p[1] > 3.1654 && p[1] < 3.1656 );
    
    std::cout << "average z" << std::endl;
    i += !( p[2] > 18.3812 && p[2] < 18.3814 );
    
    p1.Set(0., 0., 0.);
    p2.Set(2., 0., 0.);
    p3.Set(2., 1., 0.);
    p4.Set(0., 1., 0.);
    
    Point p5(0., 0., 3.);
    Point p6(2., 0., 3.);
    Point p7(2., 1., 3.);
    Point p8(0., 1., 3.);
    
    pts[4] = &p5;
    pts[5] = &p6;
    pts[6] = &p7;
    pts[7] = &p8;
    
    d = shape.VolHex();
    i += !( d > 5.99 && d < 6.01 );
    
    return i;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::ShapeTest test;
  return test.Exec();
}