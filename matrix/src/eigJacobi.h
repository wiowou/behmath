//<license>
/*
   This file is part of Behmeth.
   Author: Behram Kapadia, wiowou@hotmail.com

    Behmeth is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Behmeth is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Behmeth.  If not, see <http://www.gnu.org/licenses/>.
*/

//</license>

#ifndef _MATH_SOLVER_EIGJACOBI_h
#define _MATH_SOLVER_EIGJACOBI_h

#include "eigenValSolver.h"
#include "matrix.h"
#include "vect.h"
#include "sparseVect.h"
#include "solver.h"

namespace math{
namespace solver{

#ifdef MYDEBUG
  class EigJacobiTest;
#endif //MYDEBUG

template< template <typename> class Storage = SparseVect, typename T = double >
class EigJacobi : public EigenValSolver<Storage,T>
{
  using EigenValSolver<Storage,T>::m_EigVect;
  using EigenValSolver<Storage,T>::m_EigVal;
  using EigenValSolver<Storage,T>::ExceedsOffDiagTol;
public:
  
  void A( Matrix<Storage,T>& matrix )
  {
    m_A = &matrix;
  }
  
  void ATranspose( Matrix<Storage,T>& matrix )
  {
    A(matrix);
  }
  
  //! calculates the the eigenvalues and eigenvectors of A (same as A Transpose)
  //! for a non-singular square matrix
  ULong Eig( ULong maxIter, T tol = 1e-5, bool calcEigVect = false )
  {
    Matrix<Storage,T>& A = *m_A;
    m_EigVal.Clear();
    
    if ( calcEigVect )
    {
      Solver<Storage,T> s;
      m_EigVect.Resize( A.Cols() );
      s.Identity(m_EigVect);      
    }
    
    Matrix<Storage,T> res( A.Rows() );
    ULong iter = 0;    
    while ( iter < maxIter && ExceedsOffDiagTol( A, tol) )
    {
      for ( ULong i = 0; i < A.Rows(); ++i )
      {
        for ( ULong j = 0; j < A[i].NNZ(); ++j )
        {
          if ( A[i].Pos(j) >= i )
          {
            break; // focus only on zeroing the lower half
          }
          if ( std::abs( A[i](j) ) <= tol )
          {
            continue;
          }
          Matrix<SparseVect,T> G( A.Rows() );
          CreateGivensRotMatrix( i, j, A, G );
          vops::Dot<T> dot;
          dot( G, A, res );
          dot( res, G, A );
          if ( calcEigVect )
          {
            TransposeGivens( i, j, G );
            dot( m_EigVect, G, res );
            m_EigVect.Swap(res);
          }
        }
      }
      ++iter;
    }
    if( iter < maxIter )
    {
      m_EigVal.Resize( A.Rows() );
      for ( ULong i = 0; i < A.Rows(); ++i )
      {
        m_EigVal[i].PushBack( i, A[i][i] );
      }      
    }
    return iter;
  }
protected:
  void CreateGivensRotMatrix( ULong i, ULong j, Matrix<Storage,T>& A, Matrix<SparseVect,T> &G )
  {
    T Aij = A[i](j);
    T Aij2 = Aij * Aij;
    T y = 0.5 * ( A[i](i) - A[j](j) );
    T d = std::abs(y) + sqrt( Aij2 + y * y );
    T r = sqrt( Aij2 + d * d );
    T c = d / r;
    T s = ( Aij / r ) *  static_cast<T>( ( T(0) < y ) - ( y < T(0) ) ) ; //last expression is signum function
    
    Solver<SparseVect,T> solver;
    solver.Identity(G);
    G[j].Set(j, c); G[j].Set(i, -s);
    G[i].Set(j, s); G[i].Set(i, c); 
  }
  void TransposeGivens( ULong i, ULong j, Matrix<SparseVect,T> &G )
  {
    G[i][j] *= -1.0;
    G[j][i] *= -1.0;    
  }
protected:
  
  Matrix<Storage,T>* m_A;
private:


#ifdef MYDEBUG
  friend class EigJacobiTest;
#endif //MYDEBUG
};

}/*solver*/ }/*math*/ 

#endif /*_MATH_SOLVER_EIGJACOBI_h */
