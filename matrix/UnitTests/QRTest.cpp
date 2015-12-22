#include "../src/QR.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace math{
namespace solver{

class QRTest
{
public:
  
  
  int Exec()
  {
    int i = 0;
    
    //create the solver object and specify dense matrix storage. 
    math::solver::QR<Vect> QR;
    
    //specify the input data in transposed form, for calculation efficiency
    double dataAT[] = 
    {
      2, -1, 1, 1, 5,
      4, -2, -1, -3, -2,
      2, 6, 5, -1, 1
    };
    Matrix<Vect> AT(3,5);
    AT.Data(dataAT);
    
    std::cout << "A Transpose:" << AT << std::endl;
    
    //provide rectangular A to the QR object and factorize it
    QR.ATranspose(AT);
    QR.Factorize();
    
    //create the input vector B
    double dataB[] = { 1, 2, 3, 4, 5 };
    Vect<double> B(5);
    B.Data(dataB);

    //create the solution vector X
    Vect<double> X(3);
    
    //pass X and B to the QR solver
    QR.X(X);
    QR.B(B);
    
    QR.Solve();
    std::cout << "A after Householder: " << AT << std::endl;
    //std::cout << "alphas: " << QR.m_alpha << std::endl;
    std::cout << "X:" << X << std::endl;
    
    Matrix<Vect>* Q = QR.Q();
    Matrix<Vect>* R = QR.R();
    
    std::cout << "Q: " << *Q << std::endl;
    std::cout << "R: " << *R << std::endl;
    
    i += CheckPair( (*Q)[0], 0, .3536 );
    i += CheckPair( (*Q)[0], 1, -.1768 );
    i += CheckPair( (*Q)[0], 2, .1768 );
    i += CheckPair( (*Q)[0], 3, .1768 );
    i += CheckPair( (*Q)[0], 4, .8839 );
    
    i += CheckPair( (*Q)[1], 0, .7343 );
    i += CheckPair( (*Q)[2], 0, .3084 );
    
    i += CheckPair( (*Q)[1], 1, -.3671 );
    i += CheckPair( (*Q)[2], 2, .5746 );
    
    std::cout << "Using QR for hw3, problem 3: " << std::endl;
    // create the inputs
    double dataAT2[] = 
    {
      0.1,0.1,0.2,0.4,0.7,0.8,1,0.9,0.7,0.4,0.2,0.1,0.1,0.1,0.05,0.05,0.05,0.05,0.05,0.05,
      0.05,0.05,0.05,0.05,0.1,0.2,0.3,0.35,0.5,0.7,0.9,1,0.9,0.7,0.5,0.2,0.1,0.05,0.05,0.05,
      0.1,0.2,0.5,0.2,0.1,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.1,0.4,0.7,0.4,0.2,0.1
    };
    AT.Resize(3,20);
    AT.Data( dataAT2 );
    
    std::cout << "A Transpose: " << AT << std::endl;

    double dataB2[] = 
    {
      0.0775,0.1125,0.2375,0.1725,0.22,0.2675,0.3525,0.355,0.3825,0.4125,0.4625,0.4875,0.4425,0.3525,0.27,0.24,0.3,0.1725,0.1025,0.0675
    };
    B.Clear();
    B.Resize(20);
    B.Data(dataB2);
    
    std::cout << "B: " << B << std::endl;
    
    X.Clear();
    X.Resize(3);
    
    QR.Clear();
    QR.ATranspose(AT);
    QR.B(B);
    QR.X(X);

    QR.Solve();
    
    std::cout << "X:" << X << std::endl;
    
    i += CheckPair( X, 0, .2);
    i += CheckPair( X, 1, .45);
    i += CheckPair( X, 2, .35);
    
    std::cout << "Eigen extraction using QR: " << std::endl;
    QR.Clear();
    //specify the input data in transposed form, for calculation efficiency. 
    //This time it's symmetric, so it doesn't matter.
    double dataA[] = 
    {
      4,-2,2,1,
      -2,6,-1,1,
      2,-1,3,-2,
      1,1,-2,5
    };
    Matrix<Vect> A(4);
    A.Data(dataA);
    
    std::cout << "A^T: " << A << std::endl;
    
    QR.ATranspose(A);
    ULong iter = QR.Eig( 50, 1e-5, true ); //50 iterations max, tolerance of 1e-5 for off diagonal entries, calculate eigenvectors too
    
    std::cout << "iterations: " << iter << std::endl;
    Matrix<Vect>* EVect = QR.EigVect();
    Matrix<SparseVect>* EVal = QR.EigVal();   
    
    std::cout << "Eigenvalues: " << *EVal << std::endl;
    std::cout << "Eigenvectors: " << *EVect << std::endl;
    
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
  math::solver::QRTest test;
  return test.Exec();
}
