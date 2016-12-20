#define MYDEBUG

#include "../src/sum.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace math{
namespace vops{

class SumTest
{
public:
  math::vops::Sum<> sum;
  std::string testDir;
  
  int Exec()
  {

    int i = 0;
    SparseVect<double> sv1;
    sv1.Set(0, -5.2);
    sv1.Set(1, -8.5);
    sv1.Set(2, 8.2);
    sv1.Set(3, 4.4);
    sv1.Set(5, 5.1);
    sv1.Set(8, 10.5);
    sv1.Set(12, 16.3);
    
    double d = sum(sv1, 3);
    if ( std::abs( d - 36.3 ) > 0.001 )
    {
      i++;
      std::cout << " sparse prob " << d << "!=" << 36.3 << std::endl;
    }
    
    d = sum(5, sv1);
    if ( std::abs( d + 1.1 ) > 0.001 )
    {
      i++;
      std::cout << " sparse prob " << d << "!=" << -1.1 << std::endl;
    }
    
    Vect<double> v1;
    v1.Resize(10);
    for ( int j = 0; j < 10; ++j )
    {
      v1[j] = j;
    }
    d = sum(v1, 5);
    if ( std::abs( d - 35.0 ) > 0.001 )
    {
      i++;
      std::cout << " dense prob " << d << "!=" << 35.0 << std::endl;
    }
    
    d = sum(5, v1);
    if ( std::abs( d - 10.0 ) > 0.001 )
    {
      i++;
      std::cout << " dense prob " << d << "!=" << 10.0 << std::endl;
    }
    
    return i;
  }

};

}/*vops*/ }/*math*/ 

int main()
{
  math::vops::SumTest test;
  return test.Exec();
}
