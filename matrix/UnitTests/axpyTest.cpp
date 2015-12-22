#include "../src/axpy.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace math{
namespace vops{

class AxpyTest
{
public:
  math::vops::Axpy<> axpy;
  SparseVect<double> sv1, sv2, sv3;
  Vect<double> v1, v2;
  
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
    
    axpy(2.0, sv1, sv2);
    
    i += CheckPair( 0, -15.6 );
    i += CheckPair( 1, -17.0 );
    i += CheckPair( 2, 16.4 );
    i += CheckPair( 3, 13.2 );
    i += CheckPair( 5, 15.3 );
    i += CheckPair( 8, 21.0 );
    i += CheckPair( 11, 11.0 );
    i += CheckPair( 12, 32.6 );
    i += CheckPair( 13, 13.0 );
    
    i += sv2.NNZ() != 9;
    
    v1.Resize(15);
    for ( int j = 0; j < 15; ++j )
    {
      v1[j] = 0.0;
    }
    v1[1] = 5.0;
    v1[11] = 3.0;
    v1[14] = 7.0;
    
    axpy(-1.0, v1, sv2, 2);
    
    i += v1.Size() != 15;
    
    i += CheckPair( v1, 0, 0.0 );
    i += CheckPair( v1, 1, 0.0 );
    i += CheckPair( v1, 2, -16.4 );
    i += CheckPair( v1, 3, -13.2 );
    i += CheckPair( v1, 5, -15.3 );
    i += CheckPair( v1, 8, -21.0 );
    i += CheckPair( v1, 11, -8.0 );
    i += CheckPair( v1, 12, -32.6 );
    i += CheckPair( v1, 13, -13.0 );
    i += CheckPair( v1, 14, 7.0 );
    
    sv3.Clear();
    sv2.Clear();
    sv2.Set(0, -5.2);
    sv2.Set(3, 4.4); 
    sv2.Set(5, 5.1);
    sv2.Set(11, 11.0);
    sv2.Set(13, 13.0);
    axpy(2.0, sv1, sv2, sv3);
    
    i += CheckPair( sv3, 0, -15.6 );
    i += CheckPair( sv3,  1, -17.0 );
    i += CheckPair( sv3,  2, 16.4 );
    i += CheckPair( sv3,  3, 13.2 );
    i += CheckPair( sv3,  5, 15.3 );
    i += CheckPair( sv3,  8, 21.0 );
    i += CheckPair( sv3,  11, 11.0 );
    i += CheckPair( sv3,  12, 32.6 );
    i += CheckPair( sv3,  13, 13.0 );

    sv3.Clear();
    axpy(2.0, sv1, sv2, sv3, 3);
    i += CheckPair( sv3, 0, 0.0 );
    i += CheckPair( sv3,  1, 0.0 );
    i += CheckPair( sv3,  2, 0.0 );
    i += CheckPair( sv3,  3, 13.2 );
    i += CheckPair( sv3,  5, 15.3 );
    i += CheckPair( sv3,  8, 21.0 );
    i += CheckPair( sv3,  11, 11.0 );
    i += CheckPair( sv3,  12, 32.6 );
    i += CheckPair( sv3,  13, 13.0 );

    sv3.Clear();
    axpy(11, 2.0, sv1, sv2, sv3 );
    i += CheckPair( sv3, 0, -15.6 );
    i += CheckPair( sv3,  1, -17.0 );
    i += CheckPair( sv3,  2, 16.4 );
    i += CheckPair( sv3,  3, 13.2 );
    i += CheckPair( sv3,  5, 15.3 );
    i += CheckPair( sv3,  8, 21.0 );
    i += CheckPair( sv3,  11, 0.0 );
    i += CheckPair( sv3,  12, 0.0 );
    i += CheckPair( sv3,  13, 0.0 );
    
    sv1.Clear();
    sv2.Clear();
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
    
    axpy(8, 2.0, sv1, sv2);
    
    i += CheckPair( 0, -15.6 );
    i += CheckPair( 1, -17.0 );
    i += CheckPair( 2, 16.4 );
    i += CheckPair( 3, 13.2 );
    i += CheckPair( 5, 15.3 );
    i += CheckPair( 8, 0.0 );
    i += CheckPair( 11, 11.0 );
    i += CheckPair( 12, 0.0 );
    i += CheckPair( 13, 13.0 );
    
    v1.Resize(15);
    for ( int j = 0; j < 15; ++j )
    {
      v1[j] = 0.0;
    }
    v1[1] = 5.0;
    v1[11] = 3.0;
    v1[14] = 7.0;
    
    return i;
  }
  
  int CheckPair( ULong pos, double d )
  {
    if ( std::abs( ( sv2[pos] - d) / d ) > 0.0001 )
    {
      std::cout << pos << " bad data " << sv2[pos] << "!=" << d << std::endl;
      return 1;
    }
    return 0;
  }
  int CheckPair( SparseVect<double> &sv, ULong pos, double d )
  {
    if ( std::abs( ( sv[pos] - d) / d ) > 0.0001 )
    {
      std::cout << pos << " bad data " << sv[pos] << "!=" << d << std::endl;
      return 1;
    }
    return 0;
  }
  int CheckPair( Vect<double> &v, ULong pos, double d )
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
    if ( std::abs( ( v[pos] - d) / d ) > 0.0001 )
    {
      std::cout << pos << " bad data " << v[pos] << "!=" << d << std::endl;
      return 1;
    }
    return 0;
  }    
  int CheckPair2( ULong pos, double d )
  {
    if ( std::abs( ( v1[pos] - d) / d ) > 0.0001 )
    {
      std::cout << pos << " bad data " << v1[pos] << "!=" << d << std::endl;
      return 1;
    }
    return 0;
  }

};

}/*vops*/ }/*math*/ 

int main()
{
  math::vops::AxpyTest test;
  return test.Exec();
}
