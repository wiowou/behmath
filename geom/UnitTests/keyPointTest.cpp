#include "../src/keyPoint.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace geom{

class KeyPointTest
{
public:
  
  int Exec()
  {
		KeyPoint kp;
    int i = 0;
    KeyPoint kp2 = kp;
		
		kp2.X(1.2);
		kp.X(5.0);
		
		std::cout << kp.ID() << std::endl;
		std::cout << kp2.ID() << std::endl;
		
    return i;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::KeyPointTest test;
  return test.Exec();
}
