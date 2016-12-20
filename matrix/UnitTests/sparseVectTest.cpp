#define MYDEBUG

#include "../src/sparseVect.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace math{

class SparseVectTest
{
public:
  math::SparseVect<double> sparseVect;
  std::string testDir;
  
  int Exec()
  {
    int i = 0;
    
    sparseVect.Set(3, 4.4);
    sparseVect.Set(2, 8.2);
    sparseVect.Set(5, 5.1);
    sparseVect.Set(1, -6.2);
    sparseVect.Set(8, 10.5);
    sparseVect.Set(12, 16.3);
    sparseVect.Set(0, -5.2);
    
    sparseVect.Set(1, -8.5);
    
    for ( int i = 0; i < sparseVect.NNZ(); ++i )
    {
      std::cout << i << "  \t " << sparseVect.m_pos[i] << " \t " << sparseVect.m_data[i] << std::endl;
    }
    
    i += CheckPair( 0, 0, -5.2 );
    i += CheckPair( 1, 1, -8.5 );
    i += CheckPair( 2, 2, 8.2 );
    i += CheckPair( 3, 3, 4.4 );
    i += CheckPair( 4, 5, 5.1 );
    i += CheckPair( 5, 8, 10.5 );
    i += CheckPair( 6, 12, 16.3 );
    
    if ( sparseVect.Find(5) != 4 )
    {
      i++;
    }
    std::cout << sparseVect.Find(5) << std::endl;
    
    i+= sparseVect.Find(9) != sparseVect.NNZ();
    std::cout << sparseVect.Find(9) << std::endl;
    
    i+= sparseVect.Find(0) != 0;
    i+= sparseVect.Find(3) != 3;
    i+= sparseVect.Find(8) != 5;
    i+= sparseVect.Find(2) != 2;
    i+= sparseVect.Find(5) != 4;
    
    unsigned long long max;
    max = sparseVect.MaxIdx();
    i += ( max != 6 );
    std::cout << max << std::endl;
    
    sparseVect.PushBack(15, 17.5);
    sparseVect.PushBack(17, 19.5);
    sparseVect.PushBack(19, 20.5);
    
    i += CheckPair( 7, 15, 17.5 );
    i += CheckPair( 8, 17, 19.5 );
    i += CheckPair( 9, 19, 20.5 );
    
    //check that the resizing worked
    i += CheckPair( 0, 0, -5.2 );
    i += CheckPair( 1, 1, -8.5 );
    i += CheckPair( 2, 2, 8.2 );
    i += CheckPair( 3, 3, 4.4 );
    i += CheckPair( 4, 5, 5.1 );
    i += CheckPair( 5, 8, 10.5 );
    i += CheckPair( 6, 12, 16.3 );
    i += CheckPair( 7, 15, 17.5 );
    i += CheckPair( 8, 17, 19.5 );
    i += CheckPair( 9, 19, 20.5 );
    
    sparseVect.RowReduce( 2.0, 6, 19 );
    i += sparseVect.Find(6) != sparseVect.NNZ();
    i += CheckPair( 9, 19, 20.5 );
    
    sparseVect.RowReduce( 2.0, 8, 19 );
    i += CheckPair( 9, 19, 41.5 );
    
    sparseVect.RowReduce( 2.0, 19, 6 );
    i += CheckPair( 5, 6, 83.0 );
    
    math::SparseVect<double> spVect( sparseVect ); //test copy constructor
    i += CheckPair2( spVect, 5, 6, 83.0 );
    
    spVect.Set(6, 57.0);
    i += CheckPair2( spVect, 5, 6, 57.0 );
    
    spVect.Swap(sparseVect);
    i += CheckPair2( spVect, 5, 6, 83.0 );
    i += CheckPair( 5, 6, 57.0 );
    
    spVect = sparseVect;
    i += CheckPair2( spVect, 5, 6, 57.0 );
    
    std::cout << spVect << std::endl;
    
    
    return i;
  }
  
  int CheckPair( unsigned long long i, unsigned long long pos, double d )
  {
    if ( sparseVect.m_pos[i] != pos )
    {
      std::cout << "bad pos " << sparseVect.m_pos[i] << "!=" << pos << std::endl;
      return 1;
    }
    if ( std::abs( ( sparseVect.m_data[i] - d) / d ) > 0.0001 )
    {
      std::cout << "bad data " << sparseVect.m_data[i] << "!=" << d << std::endl;
      return 1;
    }
    return 0;
  }
  
  int CheckPair2( math::SparseVect<double> &spVect, unsigned long long i, unsigned long long pos, double d )
  {
    if ( spVect.m_pos[i] != pos )
    {
      std::cout << "bad pos " << spVect.m_pos[i] << "!=" << pos << std::endl;
      return 1;
    }
    if ( std::abs(spVect.m_data[i] - d) / d > 0.0001 )
    {
      std::cout << "bad data " << spVect.m_data[i] << "!=" << d << std::endl;
      return 1;
    }
    return 0;
  }

};

}/*math*/ 

int main()
{
  math::SparseVectTest test;
  return test.Exec();
}
