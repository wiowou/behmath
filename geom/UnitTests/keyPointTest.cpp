#define MYDEBUG

#include "../src/keyPoint.h"
#include <iostream>
#include <set>
#include <list>

namespace math{
namespace geom{

class KeyPointTest
{
public:
  
  int Exec()
  {
    int i = 0;
		KeyPoint kp;
    KeyPoint kp2;
    KeyPoint kp3;
    kp3.X(2.3);
    KeyPoint kp4 = kp3;
    kp = kp3;
		
		kp2.X(1.2);
		kp.X(5.0);
		
		std::cout << kp.ID() << std::endl;
		std::cout << kp2.ID() << std::endl;
    std::cout << kp3.ID() << std::endl;
    
    std::list<KeyPoint> lkp;
		lkp.push_back(kp);
    std::list<KeyPoint>::iterator lit = lkp.begin();
    
    std::cout << kp.ID() << std::endl;
    std::cout << lit->ID() << std::endl;
    return i;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::KeyPointTest test;
  return test.Exec();
}
