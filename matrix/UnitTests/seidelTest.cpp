#define MYDEBUG

#include "../src/seidel.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace solver{

class SeidelTest
{
public:
  int Exec()
  {

    int i = 0;
    //specify solver
    solver::Seidel<Vect> seidel;
    //Create the matrix
    double dataA[] = 
    {
      10, 4, 2, 4, 2,
      -1, -12, 6, -4, 1,
      1, -1, 15, 1, 1,
      1, 4, -1, 12, -1,
      5, -2, -1, 2, -13
    };
    Matrix<Vect> A(5);
    A.Data(dataA);
    
    //Create the input vector
    double dataB[] = { 1, -1, 2, 3, 0 };
    Vect<double> B(5);
    B.Data(dataB);
    
    //Create the output vector
    Vect<double> X(5);
    
    //Specify the inputs to the solver
    seidel.Tolerance(1e-5);
    seidel.MaxIter(100); // max num of iterations before quitting
    seidel.A(A);
    seidel.X(X);
    seidel.B(B);
    
    //call solve 
    seidel.Solve();
    std::cout << " Gauss-Seidel method solves in " << seidel.NIter() << " iterations." << std::endl;
    std::cout << "A:" << A << std::endl;
    std::cout << "B:" << B << std::endl;
    std::cout << "X:" << X << std::endl;  
    
    //checking output for automated unit testing
    
    //Create an object to calculate dot product
    vops::Dot<double> dot;
    //calc dot product of X with A and check against B
    Vect<double> Bcheck(5);
    dot(A,X,Bcheck);
    std::cout << "B checked:" << Bcheck << std::endl;
    i += CheckPair( Bcheck, 0, 1.0);
    i += CheckPair( Bcheck, 1, -1.0);
    i += CheckPair( Bcheck, 2, 2.0);
    i += CheckPair( Bcheck, 3, 3.0);
    i += CheckPair( Bcheck, 4, 0.0);
    
    return i;
  }
  //helper method
  int CheckPair( SparseVect<double> &sv, unsigned long long pos, double d )
  {
    if ( std::abs(sv[pos] - d) > 1e-4 )
    {
      std::cout << pos << " bad data " << sv[pos] << "!=" << d << std::endl;
      return 1;
    }
    return 0;
  }
  //helper method
  int CheckPair( Vect<double> &v, unsigned long long pos, double d )
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
  math::solver::SeidelTest test;
  return test.Exec();
}
