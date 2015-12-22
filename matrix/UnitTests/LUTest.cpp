#include "../src/LU.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace solver{

class LUTest
{
public:
  
  int Exec()
  {
    int i = 0;
    
    std::cout << " LU solution 4:" << std::endl;
    
    //specify that we want an LU solver that stores a matrix as dense vectors
    solver::LU<Vect> LU;
    
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
    
    std::cout << "A:" << A << std::endl;
    //set the pivot tolerance
    LU.Tolerance(0.001);
    LU.A(A);
    LU.Factorize();
    std::cout << "A:" << A << std::endl;
    Matrix<SparseVect>* P;
    //provide the method with a blank sparse matrix so that the P matrix can be printed.
    P = LU.PMatrix();
    std::cout << " P:" << *P << std::endl;
    
    std::cout << " LU solution 5:" << std::endl;
    
    //create the input vector B
    double dataB[] = { 1, -1, 2, 3, 0 };
    Vect<double> B(5);
    B.Data(dataB);
    
    //create the input vector X
    Vect<double> X(5);
    
    //give the LU solver the vectors
    LU.B(B);
    LU.X(X); 
    
    //use the matrices in the LU solver from solution 4 as the means to solve the system
    LU.Solve();
    
    std::cout << "X:" << X << std::endl; 
    
    std::cout << " LU solution 6:" << std::endl;
    
    //use data from A in problems 4 and 5
    
    
    //Create a new matrix, Ainv and feed it to the LU solver. 
    Matrix<Vect> Ainv;
    LU.InverseTranspose(Ainv);
    Ainv.Transpose(); //only necessary for printing purposes
    std::cout << "Ainv:" << Ainv << std::endl;

    //check the answer
    vops::Dot<double> dot;
    
    Matrix<Vect> Achk(5);
    Achk.Data(dataA);
    Matrix<Vect> Ident(5);
    dot( Achk, Ainv, Ident);
    
    std::cout << "A * Ainv:" << Ident << std::endl;
    for ( int ii = 0; ii < 5; ++ii )
    {
      for ( int j = 0; j < 5; ++j )
      {
        if ( ii == j )
        {
          i += CheckPair( Ident[ii], ii, 1.0 );
        }
        else
        {
          i += CheckPair( Ident[ii], j, 0.0 );
        }
      } 
    }
    
    

    return i;
  }

  int CheckPair( SparseVect<double> &sv, ULong pos, double d )
  {
    if ( std::abs(sv[pos] - d) > 1e-4 )
    {
      std::cout << pos << " bad data " << sv[pos] << "!=" << d << std::endl;
      return 1;
    }
    return 0;
  }
  //helper method
  int CheckPair( Vect<double> &v, ULong pos, double d )
  {
    if ( std::abs(d) < 1e-4 )
    {
      if ( v[pos] > 1e-4 )
      {
        std::cout << pos << " bad data " << v[pos] << "!=" << d << std::endl;
        return 1;
      }
      return 0;
    }
    if ( std::abs(v[pos] - d) > 1e-4 )
    {
      std::cout << pos << " bad data " << v[pos] << "!=" << d << std::endl;
      return 1;
    }
    return 0;
  }  

};

}/*solver*/ }/*math*/ 

int main()
{
  math::solver::LUTest test;
  return test.Exec();
}
