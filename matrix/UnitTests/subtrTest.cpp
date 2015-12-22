#include "../src/subtr.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace math{
namespace vops{

class SubtrTest
{
public:
  math::vops::Subtr<> subtr;
  SparseVect<double> sv1, sv2, sv3;
  Vect<double> v1, v2;
  std::string testDir;
  
  int Exec()
  {
    int i = 0;

    sv1.Set(0, -5.2);
    sv1.Set(1, -8.5);
    sv1.Set(2, 8.2);
    sv1.Set(3, 4.4);
    sv1.Set(5, 5.1);
    sv1.Set(8, 10.5);
    sv1.Set(12, 16.3);
    
    
    sv2.Set(0, -5.2);
    sv2.Set(3, 4.4); 
    sv2.Set(5, 5.1);
    sv2.Set(11, 11.0);
    sv2.Set(13, 13.0);
    
    subtr(sv1, sv2, sv3);
    
    i += CheckPair( 0, 0.0 );
    i += CheckPair( 1, -8.5 );
    i += CheckPair( 2, 8.2 );
    i += CheckPair( 3, 0.0 );
    i += CheckPair( 5, 0.0 );
    i += CheckPair( 8, 10.5 );
    i += CheckPair( 11, -11.0 );
    i += CheckPair( 12, 16.3 );
    i += CheckPair( 13, -13.0 );
    
    i += sv3.NNZ() != 9;
    
    sv3.Clear();
    subtr(sv1, sv2, sv3, 2);
    i += CheckPair( 2, 8.2 );
    i += CheckPair( 3, 0.0 );
    i += CheckPair( 5, 0.0 );
    i += CheckPair( 8, 10.5 );
    i += CheckPair( 11, -11.0 );
    i += CheckPair( 12, 16.3 );
    i += CheckPair( 13, -13.0 );
    
    i += sv3.NNZ() != 7;
    
    sv3.Clear();
    subtr(11, sv1, sv2, sv3 );
    i += CheckPair( 0, 0.0 );
    i += CheckPair( 1, -8.5 );
    i += CheckPair( 2, 8.2 );
    i += CheckPair( 3, 0.0 );
    i += CheckPair( 5, 0.0 );
    i += CheckPair( 8, 10.5 );

    i += sv3.NNZ() != 6;
    
    v1.Resize(15);
    for ( int j = 0; j < 15; ++j )
    {
      v1[j] = 0.0;
    }
    v1[1] = 5.0;
    v1[11] = 3.0;
    v1[14] = 7.0;
    
    v2.Resize(15);
    subtr(sv2, v1, v2);
    
    i += v2.Size() != 15;
    i += CheckPair2(0, -5.2);
    i += CheckPair2(1, -5.0);
    i += CheckPair2(3, 4.4); 
    i += CheckPair2(5, 5.1);
    i += CheckPair2(11, 8.0);
    i += CheckPair2(13, 13.0);
    i += CheckPair2(14, -7.0);
    
    for ( int i = 0; i < v2.Size(); ++i )
    {
      v2[i] = 0.0;
    }
    subtr(sv2, v1, v2, 11);
    i += CheckPair2(3, 0.0);
    i += CheckPair2(5, 0.0);
    i += CheckPair2(11, 8.0);
    i += CheckPair2(13, 13.0);
    i += CheckPair2(14, -7.0);  

    for ( int i = 0; i < v2.Size(); ++i )
    {
      v2[i] = 0.0;
    }
    subtr(11, sv2, v1, v2 );
    i += CheckPair2(0, -5.2);
    i += CheckPair2(1, -5.0);
    i += CheckPair2(3, 4.4); 
    i += CheckPair2(5, 5.1);
    i += CheckPair2(11, 0.0);
    i += CheckPair2(13, 0.0);
    i += CheckPair2(14, 0.0);
    
    Vect<double> v3;
    v3.Resize(v1.Size());
    
    for ( int i = 0; i < v3.Size(); ++i )
    {
      v2[i] = v1[1] = 0.0;
      v3[i] = i;
    }
    v1[1] = 5.0;
    v1[11] = 3.0;
    v1[14] = 7.0;
    
    subtr(3, v1, v3, v2);
    i += CheckPair2(0, 0.0);
    i += CheckPair2(1, 4.0);
    i += CheckPair2(2, -2.0);
    i += CheckPair2(3, 0.0);

    for ( int i = 0; i < v3.Size(); ++i )
    {
      v2[i] = v1[1] = 0.0;
      v3[i] = i;
    }
    v1[1] = 5.0;
    v1[11] = 3.0;
    v1[14] = 7.0;
    subtr( v1, v3, v2, 7);
    i += CheckPair2(6, 0.0);
    i += CheckPair2(7, -7.0);
    i += CheckPair2(10, -10.0);
    i += CheckPair2(11, -8.0);
    i += CheckPair2(14, -7.0);
    
    return i;
  }
  
  int CheckPair( ULong pos, double d )
  {
    if ( std::abs( ( sv3[pos] - d) ) / d > 0.0001 )
    {
      std::cout << pos << " bad data " << sv3[pos] << "!=" << d << std::endl;
      return 1;
    }
    return 0;
  }
  
  int CheckPair2( ULong pos, double d )
  {
    if ( std::abs(d) < 1e-10 )
    {
      if ( v2[pos] > 1e-10 )
      {
        std::cout << pos << " bad data " << v2[pos] << "!=" << d << std::endl;
        return 1;
      }
      return 0;
    }
    if ( std::abs( ( v2[pos] - d) / d ) > 0.0001 )
    {
      std::cout << pos << " bad data " << v2[pos] << "!=" << d << std::endl;
      return 1;
    }
    return 0;
  }

};

}/*vops*/ }/*math*/ 

int main()
{
  math::vops::SubtrTest test;
  return test.Exec();
}
