#define MYDEBUG

#include "../src/node.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace geom{

class NodeTest
{
public:

  int Exec()
  {
    Node node;
    int i = 0;
    
    return i;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::NodeTest test;
  return test.Exec();
}
