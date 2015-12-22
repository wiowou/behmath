#include "../src/jacobi.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace solver{

class JacobiTest
{
public:
  
  
  int Exec()
  {
    std::cout << "Jacobi solution 8: "  << std::endl;
    int i = 0;
    math::solver::Jacobi<Vect> itSolver;
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
    itSolver.Tolerance(1e-5);
    itSolver.MaxIter(100); // max num of iterations before quitting
    itSolver.A(A);
    itSolver.X(X);
    itSolver.B(B);
    
    //call solve 
    itSolver.Solve();
    std::cout << " Jacobi method solves in " << itSolver.NIter() << " iterations." << std::endl;
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
  math::solver::JacobiTest test;
  return test.Exec();
}
