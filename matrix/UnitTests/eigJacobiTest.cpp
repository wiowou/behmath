#define MYDEBUG

#include "../src/eigJacobi.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace solver{

class EigJacobiTest
{
public:
  math::solver::EigJacobi<Vect> eigJacobi;
  
  int Exec()
  {
    int i = 0;
    
    //specify the data as-is, not transposed, though transposition makes no difference to eigenvals
    double dataA[] = 
    {
      4,-2,2,1,
      -2,6,-1,1,
      2,-1,3,-2,
      1,1,-2,5
    };
    Matrix<Vect> A(4);
    A.Data(dataA);

    eigJacobi.A(A);
    eigJacobi.Eig( 20, 1e-5, true ); //20 iterations max, tolerance of 1e-5 for off diagonal entries, calculate eigenvectors too
    Matrix<Vect>* EVect = eigJacobi.EigVect();
    Matrix<SparseVect>* EVal = eigJacobi.EigVal();
    
    std::cout << "A: " << A << std::endl;
    std::cout << "Eigenvalues: " << *EVal << std::endl;
    std::cout << "Eigenvectors: " << *EVect << std::endl;
    
    return i;
  }

};

}/*solver*/ }/*math*/ 

int main()
{
  math::solver::EigJacobiTest test;
  return test.Exec();
}
