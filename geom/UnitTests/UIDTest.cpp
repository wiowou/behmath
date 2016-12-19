#include "../src/UID.h"
#include <iostream>
#include <fstream>
#include <string>

#include <memory>
#include <vector>

namespace math{
namespace geom{

class UIDTest
{
public:
  
  int Exec()
  {
		UID<> u1;
		std::cout << "u1 " << u1.ID() << std::endl;
		UID<> u2 = u1;
		std::cout << "u2 " << u2.ID() << std::endl;
		
		UID<> u3;
		std::cout << "u3 " << u3.ID() << std::endl;
		
		u3 = u1;
		std::cout << "u1 " << u1.ID() << std::endl;
		std::cout << "u3 " << u3.ID() << std::endl;
		
    int i = 0;
    
    return i;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::UIDTest test;
  return test.Exec();
}
