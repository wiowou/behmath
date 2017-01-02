#define MYDEBUG

#include "../src/face.h"
#include "../src/point.cpp"
#include "../src/vect.cpp"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace geom{

class FaceTest
{
public:
  Face face;
  
  int Exec()
  {
    int i = 0;
    
    return i;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::FaceTest test;
  return test.Exec();
}
