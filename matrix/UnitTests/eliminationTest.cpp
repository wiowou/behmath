#define MYDEBUG

#include "../src/elimination.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace solver{

class EliminationTest
{
public:
  
  int Exec()
  {
    int i = 0;
    
    //specify that we want a gaussian elim solver that stores a matrix as dense vectors
    solver::Elimination<Vect> gauss;
    
    //specify the inputs to the solver
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
    
    Vect<double> X(5);
    
    gauss.A(A);
    gauss.B(B);
    gauss.X(X);
    gauss.Tolerance(0.001);
    
    gauss.Solve();
    std::cout << "A:" << A << std::endl; 
    
    //Create an object to calculate dot product
    vops::Dot<double> dot;
    //give A and B back their original values
    A.Data(dataA);
    B.Data(dataB);
    std::cout << " Gauss Elim solution 1:" << std::endl;
    std::cout << "A:" << A << std::endl;
    std::cout << "B:" << B << std::endl;
    std::cout << "X:" << X << std::endl; 
    //calc dot product of X with A and check against B
    Vect<double> Bcheck(5);
    dot(A,X,Bcheck);
    std::cout << "B checked:" << Bcheck << std::endl;
    i += CheckPair( Bcheck, 0, 1.0);
    i += CheckPair( Bcheck, 1, -1.0);
    i += CheckPair( Bcheck, 2, 2.0);
    i += CheckPair( Bcheck, 3, 3.0);
    i += CheckPair( Bcheck, 4, 0.0);
    
    std::cout << "Gauss Elim solution 2" << std::endl;
    //specify that we want a gaussian elim solver that stores a matrix as sparse vectors
    solver::Elimination<SparseVect> gauss2;
    
    //specify the inputs to the solver
    
    Matrix<SparseVect> A2(6);
    double row1[] = 
    {
      2, 3, -1 
    };
    unsigned long long row1pos[] =
    {
      0, 1, 2
    };
    A2[0].Data(row1, row1pos, 3);
    
    double row2[] = 
    {
      3, 6, -3, 3 
    };
    unsigned long long row2pos[] =
    {
      0, 1, 2, 3
    };
    A2[1].Data(row2, row2pos, 4);
    
    double row3[] = 
    {
      -1, -3, 9, 4, -2 
    };
    unsigned long long row3pos[] =
    {
      0, 1, 2, 3, 4
    };
    A2[2].Data(row3, row3pos, 5);
    
    double row4[] = 
    {
      3, 4, -3, 5, 2 
    };
    unsigned long long row4pos[] =
    {
      1, 2, 3, 4, 5
    };
    A2[3].Data(row4, row4pos, 5);
    
    double row5[] = 
    {
      -2, 5, 3, -3 
    };
    unsigned long long row5pos[] =
    {
      2, 3, 4, 5
    };
    A2[4].Data(row5, row5pos, 4);
    
    double row6[] = 
    {
      2, -3, 4 
    };
    unsigned long long row6pos[] =
    {
      3, 4, 5
    };
    A2[5].Data(row6, row6pos, 3);
    
    double dataB2[] = { 1, -1, 2, 3, 0, 1 };
    Vect<double> B2(6);
    B2.Data(dataB2);
    
    Vect<double> X2(6);
    
    std::cout << "A:" << A2 << std::endl;
    std::cout << "B:" << B2 << std::endl;
   
    gauss2.A(A2);
    gauss2.B(B2);
    gauss2.X(X2);
    gauss2.Solve();
    std::cout << "X:" << X2 << std::endl; 
    std::cout << "A:" << A2 << std::endl;

    //give A and B back their original values
    B2.Data(dataB2);
    A2[0].Data(row1, row1pos, 3);
    A2[1].Data(row2, row2pos, 4);
    A2[2].Data(row3, row3pos, 5);
    A2[3].Data(row4, row4pos, 5);
    A2[4].Data(row5, row5pos, 4);
    A2[5].Data(row6, row6pos, 3);
    
    //calc dot product of X with A and check against B
    Vect<double> Bcheck2(6);
    dot(A2,X2,Bcheck2);
    std::cout << "B checked:" << Bcheck << std::endl;
    i += CheckPair( Bcheck2, 0, 1.0);
    i += CheckPair( Bcheck2, 1, -1.0);
    i += CheckPair( Bcheck2, 2, 2.0);
    i += CheckPair( Bcheck2, 3, 3.0);
    i += CheckPair( Bcheck2, 4, 0.0);
    i += CheckPair( Bcheck2, 5, 1.0);
    
    std::cout << " Gauss Elim solution 3:" << std::endl;
    
    //use gauss solver 1 and data from A in problem 1
    
    
    //Create a new matrix, Ainv and initialize it to the identity matrix.
    //Ainv will be fed row by row to the gauss solver, since my matrices are in row major storage format
    //It will then be transposed so that Ainv can be recognized by humans.
    Matrix<Vect> Ainv(5);
    for ( int i = 0; i < 5; ++i )
    {
      for ( int j = 0; j < 5; ++j )
      {
        Ainv(i,j) = 0.0;
      }
      Ainv(i,i) = 1.0;
    }
    std::cout << "Ainv initial:" << Ainv << std::endl;
    
    X.Clear();
    X.Resize(5);
    gauss.X(X);
    for ( int i = 0; i < 5; ++i )
    {
      A.Data(dataA);
      gauss.B( Ainv[i] );
      gauss.Solve();
      Ainv[i].Swap(X); //the swap is a fast exchange of pointers and places the solution into Ainv
    }
    Ainv.Transpose(); //only necessary for printing purposes
    std::cout << "Ainv:" << Ainv << std::endl;
    
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
  math::solver::EliminationTest test;
  return test.Exec();
}
