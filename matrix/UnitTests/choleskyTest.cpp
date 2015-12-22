#include "../src/cholesky.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace math{
namespace solver{

class CholeskyTest
{
public:
  
  int Exec()
  {
    int i = 0;
    std::cout << " Cholesky solution 7:" << std::endl;
    
    //We can take advantage of symmetry by creating a sparse matrix with only the lower half defined.
    Matrix<SparseVect> C(5);
    C[0].PushBack(0, 15.0 );
    C[1].PushBack(0, 4.0 ); C[1].PushBack(1, 10.0 );
    C[2].PushBack(0, 2.0 ); C[2].PushBack(1, 6.0 ); C[2].PushBack(2, 8.0 );
    C[3].PushBack(0, 4.0 ); C[3].PushBack(1, -4.0 ); C[3].PushBack(2, 1.0 ); C[3].PushBack(3, 10.0 );
    C[4].PushBack(0, 5.0 ); C[4].PushBack(1, 1.0 ); C[4].PushBack(2, 1.0 ); C[4].PushBack(3, -1.0 ); C[4].PushBack(4, 9.0 );    
    
    std::cout << "Lower Half of C: " << C << std::endl;

    //create the cholesky solver
    math::solver::Cholesky<SparseVect> chol;
    
    //provide the matrix as the L argument to the solver
    chol.L( C );
    
    //ask the solver to factorize matrix C
    chol.Factorize();
    
    std::cout << "C factorized into L^T*L: " << C << std::endl;
    
    //check the answer
    i+= CheckPair( C[0], 0, 3.873 );
    i+= CheckPair( C[2], 0, 0.5164 );
    i+= CheckPair( C[2], 1, 1.829 );
    i+= CheckPair( C[4], 2, 0.2565 );
    i+= CheckPair( C[4], 4, 2.1175 );

    std::cout << "Solve a system using Cholesky: " << std::endl;

    //create the B vector
    double dataB[] = { 1, -1, 2, 3, 0 };
    Vect<double> B(5);
    B.Data(dataB);
    
    std::cout << "B:" << B << std::endl;
    //create the X vector
    Vect<double> X(5);
    
    //pass to solver
    chol.B(B);
    chol.X(X);
    
    chol.Solve();
    
    std::cout << "X:" << X << std::endl;

    //check the answer
    i += CheckPair( X, 0, .2816 );
    i += CheckPair( X, 1, -.75980 );
    i += CheckPair( X, 2, .7994 );
    i += CheckPair( X, 3, -.2150 );
    i += CheckPair( X, 4, -.1847 );
    
    //We can do banded matrices as well
    Matrix<SparseVect> C2(5);
    C2[0].PushBack(0, 15.0 );
    C2[1].PushBack(0, 4.0 ); C2[1].PushBack(1, 10.0 );
    /*C2[2].PushBack(0, 2.0 );*/ C2[2].PushBack(1, 6.0 ); C2[2].PushBack(2, 8.0 );
    /*C2[3].PushBack(0, 4.0 );*/ C2[3].PushBack(1, -4.0 ); C2[3].PushBack(2, 1.0 ); C2[3].PushBack(3, 10.0 );
    /*C2[4].PushBack(0, 5.0 ); C2[4].PushBack(1, 1.0 );*/ C2[4].PushBack(2, 1.0 ); C2[4].PushBack(3, -1.0 ); C2[4].PushBack(4, 9.0 );
    
    //provide the matrix as the L argument to the solver
    chol.L( C2 );
    //ask the solver to factorize matrix C2
    chol.Factorize();
    
    std::cout << "C2 factorized into L^T*L: " << C2 << std::endl;
    
    std::cout << " Cholesky hw3, solution 1:" << std::endl;
    double dataA[] = 
    {
      15.0, 4.0, 2.0, 4.0, 5.0,
      4.0, 10.0, 6.0, -4.0, 1.0,
      2.0, 6.0, 8.0, 1.0, 1.0,
      4.0, -4.0, 1.0, 10.0, -1.0,
      5.0, 1.0, 1.0, -1.0, 9.0
    };
    
    //create a dense matrix and fill it with the symmetric matrix values
    Matrix<Vect> A(5);
    A.Data(dataA);
    std::cout << "A:" << A << std::endl;
    
    //create a new solver and ask the solver to factorize A
    solver::Cholesky<Vect> chol2;
    chol2.L(A);
    chol2.Factorize();
    std::cout << "A factorized into L^T*L: " << A << std::endl;
    std::cout << "B:" << B << std::endl;
    //clear X
    X.Clear();
    X.Resize(5);
    //pass to solver
    B.Data(dataB);
    chol2.B(B);
    chol2.X(X);
    
    chol2.Solve();
    
    std::cout << "X:" << X << std::endl;

    //check the answer
    i += CheckPair( X, 0, .2816 );
    i += CheckPair( X, 1, -.75980 );
    i += CheckPair( X, 2, .7994 );
    i += CheckPair( X, 3, -.2150 );
    i += CheckPair( X, 4, -.1847 );
    
    
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
  math::solver::CholeskyTest test;
  return test.Exec();
}
