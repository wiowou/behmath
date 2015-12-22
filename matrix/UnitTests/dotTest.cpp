#include "../src/dot.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "matrix.h"

namespace math{
namespace vops{

class DotTest
{
public:
  math::vops::Dot<> dot;
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
    sv2.Set(1, -8.5);
    sv2.Set(2, 8.2);
    sv2.Set(3, 4.4);
    sv2.Set(5, 5.1);
    sv2.Set(8, 10.5);
    sv2.Set(12, 16.3);
    
    double d = dot(sv1,sv2, 3);
    if ( std::abs( d - 421.31 ) > 0.001 )
    {
      i++;
      std::cout << " sparse prob " << d << "!=" << 421.31 << std::endl;
    }
    
    d = dot(5,sv1,sv2);
    if ( std::abs( d - 185.89 ) > 0.01 )
    {
      i++;
      std::cout << " sparse prob 2 " << d << "!=" << 185.89 << std::endl;
    }
    
    Vect<double> v1;
    v1.Resize(13);
    for ( int j = 0; j < 13; ++j )
    {
      v1[j] = j;
    }
    
    d = dot(sv1,v1, 3);
    if ( std::abs( d - 318.3 ) > 0.001 )
    {
      i++;
      std::cout << " sparse-dense prob " << d << "!=" << 318.3 << std::endl;
    }
    
    d = dot(3,sv1,v1);
    if ( std::abs( d - 7.9 ) > 0.01 )
    {
      i++;
      std::cout << " sparse-dense prob 2 " << d << "!=" << 7.9 << std::endl;
    }
    
    double dataA[] = 
    {
      2, 4, 2, 4, 2,
      -1, -2, 6, -4, 1,
      1, -1, 5, 1, 1,
      1, 4, -1, 2, -1,
      5, -2, -1, 2, -3
    };
    Matrix<Vect> A(5);
    A.Data(dataA);

    double dataB[] = { 1, -1, 2, 3, 0 };
    Vect<double> B(5);
    B.Data(dataB);
    
    Vect<double> v5(5);
    dot(A,B,v5);
    
    i += CheckPair( v5, 0, 14.0);
    i += CheckPair( v5, 1, 1.0);
    i += CheckPair( v5, 2, 15.0);
    i += CheckPair( v5, 3, 1.0);
    i += CheckPair( v5, 4, 11.0);
   
    Matrix<SparseVect> ASparse(5);
    ASparse.Set(0,0, 2.0);
    ASparse.Set(0,2, 2.0);
    ASparse.Set(0,3, 4.0);
    
    ASparse.Set(1,0, -1.0);
    ASparse.Set(1,1, -2.0);
    ASparse.Set(1,3, -4.0);
    
    ASparse.Set(2,0, 1.0);
    ASparse.Set(2,2, 5.0);
    ASparse.Set(2,3, 1.0);
    ASparse.Set(2,4, 1.0);
    
    ASparse.Set(3,0, 1.0);
    ASparse.Set(3,2, -1.0);
    
    ASparse.Set(4,3, 2.0);
    ASparse.Set(4,4, -3.0);

    dot(ASparse, B, v5 );
     
    i += CheckPair( v5, 0, 18.0);
    i += CheckPair( v5, 1, -11.0);
    i += CheckPair( v5, 2, 14.0);
    i += CheckPair( v5, 3, -1.0);
    i += CheckPair( v5, 4, 6.0);

    return i;
  }

  int CheckPair( SparseVect<double> &sv, ULong pos, double d )
  {
    if ( std::abs ( (sv[pos] - d) / d ) > 0.0001 )
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
    if ( std::abs(( v[pos] - d) / d ) > 0.0001 )
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
  math::vops::DotTest test;
  return test.Exec();
}
