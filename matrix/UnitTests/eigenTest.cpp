#include "../src/eigen.h"
#include "../src/QR.h"
#include "../src/eigJacobi.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace solver{

class EigenTest
{
public:
  
  
  int Exec()
  {
    int i = 0;
    //create the eigen solver that uses the QR method for eigenvalue extraction
    Eigen<QR,Vect> eigenQR;
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
    
    eigenQR.ATranspose( A );
    eigenQR.MaxIter(100);
    eigenQR.Tolerance(1e-5);
    Matrix<SparseVect>* EVal = eigenQR.EigVal();
    std::cout << "iterations : " << eigenQR.Iter() << std::endl;
    std::cout << "eigenvalues using QR: " << *EVal << std::endl;
    
    std::cout << "A: " << A << std::endl;
    
    Matrix<Vect>* EVect = eigenQR.EigVect();
    std::cout << "eigenvectors using QR: " << *EVect << std::endl;
    
    //create the eigen solver that uses the Jacobi method for eigenvalue extraction
    Eigen<EigJacobi,Vect> eigenJac;
    A.Data(dataA);
    eigenJac.ATranspose( A );
    eigenJac.MaxIter(100);
    eigenJac.Tolerance(1e-5);
    EVal = eigenJac.EigVal();
    std::cout << "iterations : " << eigenJac.Iter() << std::endl;
    std::cout << "eigenvalues using Jacobi: " << *EVal << std::endl;
    
    EVect = eigenJac.EigVect();
    std::cout << "eigenvectors using Jacobi: " << *EVect << std::endl;

    std::cout << " calculate SVD using Jacobi: " << std::endl;
    double dataATsvd[] = 
    {
      4,4,3,1,
      2,6,-1,1,
      1,-2,2,-2,
      1,5,-2,5
    };
    Matrix<Vect> ATsvd(4);
    ATsvd.Data(dataATsvd);
    //ATsvd.Transpose();
    
    std::cout << " A transpose: " << ATsvd << std::endl;
    eigenJac.ATranspose(ATsvd); //provide A^T to eigenJac for efficiency purposes
    eigenJac.SVD();
    Matrix<Vect>* U = eigenJac.U();
    Matrix<SparseVect>* Sigma = eigenJac.S();
    Matrix<Vect>* VT = eigenJac.VT();
    
    std::cout << "U: " << *U << std::endl;
    std::cout << "Sigma: " << *Sigma << std::endl;
    std::cout << "V^T: " << *VT << std::endl;
    
    return i;
  }

};

}/*solver*/ }/*math*/ 

int main()
{
  math::solver::EigenTest test;
  return test.Exec();
}
