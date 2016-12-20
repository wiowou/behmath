#define MYDEBUG

#include "../src/mult.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace math{
namespace vops{

class MultTest
{
public:
  math::vops::Mult<> mult;
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

    SparseVect<double> sv2;
    sv2.Set(0, -5.2);
    sv2.Set(3, 4.4); 
    sv2.Set(5, 5.1);
    sv2.Set(11, 11.0);
    sv2.Set(13, 13.0);

    SparseVect<double> sv3;
    
    mult( sv1, sv2, sv3, 3);
    i += CheckPair( sv3, 0, 0.0);
    i += CheckPair( sv3, 2, 0.0);
    i += CheckPair( sv3, 3, 19.36);
    i += CheckPair( sv3, 4, 0.0);
    i += CheckPair( sv3, 5, 26.01);
    i += CheckPair( sv3, 6, 0.0);
    
    i += sv3.NNZ() != 2;
    
    sv3.Clear();
    mult( 5, sv1, sv2, sv3 );
    i += CheckPair( sv3, 0, 27.04);
    i += CheckPair( sv3, 1, 0.0);
    i += CheckPair( sv3, 3, 19.36);
    i += CheckPair( sv3, 4, 0.0);
    i += CheckPair( sv3, 5, 0.0);
    
    i += sv3.NNZ() != 2;
    
    Vect<double> v1;
    v1.Resize(15);
    for (int j = 0; j < 15; ++j )
    {
      v1[j] = j;
    }
    mult( sv1, sv2, v1, 3 );
    i += CheckPair( v1, 0, 0.0);
    i += CheckPair( v1, 1, 1.0);
    i += CheckPair( v1, 2, 2.0);
    i += CheckPair( v1, 3, 19.36);
    i += CheckPair( v1, 4, 0.0);
    i += CheckPair( v1, 5, 26.01);
    i += CheckPair( v1, 6, 0.0);
    i += CheckPair( v1, 8, 0.0);
    i += CheckPair( v1, 11, 0.0);
    i += CheckPair( v1, 12, 0.0);
    
    for (int j = 0; j < 15; ++j )
    {
      v1[j] = j;
    }    
    mult( 5, sv1, sv2, v1 );
    i += CheckPair( v1, 0, 27.04);
    i += CheckPair( v1, 1, 0.0);
    i += CheckPair( v1, 3, 19.36);
    i += CheckPair( v1, 4, 0.0);
    i += CheckPair( v1, 5, 5.0);

    for (int j = 0; j < 15; ++j )
    {
      v1[j] = j;
    }

    sv3.Clear();
    mult( sv2, v1, sv3, 5);
    i += CheckPair( sv3, 0, 0.0);
    i += CheckPair( sv3, 3, 0.0);
    i += CheckPair( sv3, 5, 25.5);
    i += CheckPair( sv3, 6, 0.0);
    i += CheckPair( sv3, 10, 0.0);
    i += CheckPair( sv3, 11, 121.0);
    
    sv3.Clear();
    mult( 11, sv2, v1, sv3 );
    i += CheckPair( sv3, 0, 0.0);
    i += CheckPair( sv3, 3, 13.2);
    i += CheckPair( sv3, 5, 25.5);
    i += CheckPair( sv3, 6, 0.0);
    i += CheckPair( sv3, 10, 0.0);
    i += CheckPair( sv3, 11, 0.0);    
    
    Vect<double> v2;
    v2.Resize(15);
    for (int j = 0; j < 15; ++j )
    {
      v2[j] = j;
      v1[j] = j;
    }
    mult( sv2, v1, v2, 5);
    i += CheckPair( v2, 0, 0.0);
    i += CheckPair( v2, 3, 3.0);
    i += CheckPair( v2, 5, 25.5);
    i += CheckPair( v2, 6, 0.0);
    i += CheckPair( v2, 10, 0.0);
    i += CheckPair( v2, 11, 121.0);
    
    for (int j = 0; j < 15; ++j )
    {
      v2[j] = j;
      v1[j] = j;
    }
    mult( 5, sv2, v1, v2 );
    i += CheckPair( v2, 0, 0.0);
    i += CheckPair( v2, 3, 13.2);
    i += CheckPair( v2, 4, 0.0);
    i += CheckPair( v2, 5, 5.0);
    i += CheckPair( v2, 6, 6.0);
    i += CheckPair( v2, 10, 10.0);
    i += CheckPair( v2, 11, 11.0);

    sv1.Clear();
    sv1.Set(0, -5.2);
    sv1.Set(1, -8.5);
    sv1.Set(2, 8.2);
    sv1.Set(3, 4.4);
    sv1.Set(5, 5.1);
    sv1.Set(8, 10.5);
    sv1.Set(12, 16.3);

    sv2.Clear();
    sv2.Set(0, -5.2);
    sv2.Set(3, 4.4); 
    sv2.Set(5, 5.1);
    sv2.Set(11, 11.0);
    sv2.Set(13, 13.0);
    
    mult( sv1, sv2, 3);
    i += CheckPair( sv1, 0, -5.2);
    i += CheckPair( sv1, 2, 8.2);
    i += CheckPair( sv1, 3, 19.36);
    i += CheckPair( sv1, 4, 0.0);
    i += CheckPair( sv1, 5, 26.01);
    i += CheckPair( sv1, 8, 0.0);
    i += CheckPair( sv1, 12, 0.0);
 
     sv1.Clear();
    sv1.Set(0, -5.2);
    sv1.Set(1, -8.5);
    sv1.Set(2, 8.2);
    sv1.Set(3, 4.4);
    sv1.Set(5, 5.1);
    sv1.Set(8, 10.5);
    sv1.Set(12, 16.3);

    sv2.Clear();
    sv2.Set(0, -5.2);
    sv2.Set(3, 4.4); 
    sv2.Set(5, 5.1);
    sv2.Set(11, 11.0);
    sv2.Set(13, 13.0);
    
    mult( 8, sv1, sv2 );
    i += CheckPair( sv1, 0, 27.04);
    i += CheckPair( sv1, 2, 0.0);
    i += CheckPair( sv1, 3, 19.36);
    i += CheckPair( sv1, 4, 0.0);
    i += CheckPair( sv1, 5, 26.01);
    i += CheckPair( sv1, 8, 10.5);
    i += CheckPair( sv1, 12, 16.3);

    sv2.Clear();
    sv2.Set(0, -5.2);
    sv2.Set(3, 4.4); 
    sv2.Set(5, 5.1);
    sv2.Set(11, 11.0);
    sv2.Set(13, 13.0);

    for (int j = 0; j < 15; ++j )
    {
      v2[j] = j;
    }
    
    mult( sv2, v2, 3);
    i += CheckPair( sv2, 0, -5.2 );
    i += CheckPair( sv2, 3, 13.2 );
    i += CheckPair( sv2, 5, 25.5 );
    i += CheckPair( sv2, 11, 121.0 );
    i += CheckPair( sv2, 13, 169.0 );

    sv2.Clear();
    sv2.Set(0, -5.2);
    sv2.Set(3, 4.4); 
    sv2.Set(5, 5.1);
    sv2.Set(11, 11.0);
    sv2.Set(13, 13.0);
    
    mult( 8, sv2, v2);
    i += CheckPair( sv2, 0, 0.0 );
    i += CheckPair( sv2, 3, 13.2 );
    i += CheckPair( sv2, 5, 25.5 );
    i += CheckPair( sv2, 11, 11.0 );
    i += CheckPair( sv2, 13, 13.0 );
    
    sv2.Clear();
    sv2.Set(0, -5.2);
    sv2.Set(3, 4.4); 
    sv2.Set(5, 5.1);
    sv2.Set(11, 11.0);
    sv2.Set(13, 13.0);

    for (int j = 0; j < 15; ++j )
    {
      v2[j] = j;
    }

    mult( v2, sv2, 4);
    i += CheckPair( v2, 0, 0.0 );
    i += CheckPair( v2, 3, 3.0 );
    i += CheckPair( v2, 4, 0.0 );
    i += CheckPair( v2, 5, 25.5 );
    i += CheckPair( v2, 11, 121.0 );
    i += CheckPair( v2, 12, 0.0 );
    i += CheckPair( v2, 13, 169.0 );
    i += CheckPair( v2, 14, 0.0 );
    
    for (int j = 0; j < 15; ++j )
    {
      v2[j] = j;
    }

    mult( 8, v2, sv2 );
    i += CheckPair( v2, 0, 0.0 );
    i += CheckPair( v2, 3, 13.2 );
    i += CheckPair( v2, 4, 0.0 );
    i += CheckPair( v2, 5, 25.5 );
    i += CheckPair( v2, 11, 11.0 );
    i += CheckPair( v2, 12, 12.0 );
    i += CheckPair( v2, 13, 13.0 );
    i += CheckPair( v2, 14, 14.0 );
    return i;
  }
  
  int CheckPair( SparseVect<double> &sv, unsigned long long pos, double d )
  {
    if ( std::abs( (sv[pos] - d) / d) > 0.001 )
    {
      std::cout << pos << " bad data " << sv[pos] << "!=" << d << std::endl;
      return 1;
    }
    return 0;
  }
  int CheckPair( Vect<double> &v, unsigned long long pos, double d )
  {
    if ( std::abs(d) < 1e-10 )
    {
      if ( v[pos] > 1e-10 )
      {
        std::cout << pos << " bad data " << v[pos] << "!=" << d << std::endl;
        return 1;
      }
      return 0;
    }
    if ( std::abs(v[pos] - d) / d > 0.0001 )
    {
      std::cout << pos << " bad data " << v[pos] << "!=" << d << std::endl;
      return 1;
    }
    return 0;
  } 

};

}/*vops*/ }/*math*/ 

int main()
{
  math::vops::MultTest test;
  return test.Exec();
}
