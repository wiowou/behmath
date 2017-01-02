#define MYDEBUG

#include "../src/volumeI.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace geom{

class VolumeITest
{
public:

  int Exec()
  {
    VolumeI volumeI;
    int i = 0;
    
    return i;
  }

};

}/*geom*/ }/*math*/ 

int main()
{
  math::geom::VolumeITest test;
  return test.Exec();
}
