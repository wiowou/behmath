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

#ifndef _MATH_SOLVER_QR_h
#define _MATH_SOLVER_QR_h

#include "eigenValSolver.h"
#include "solver.h"
#include "matrix.h"
#include "vect.h"
#include "sparseVect.h"
#include "dot.h"
#include <math.h> //for sqrt
#include <cmath> //for std::abs

namespace math{
namespace solver{

#ifdef MYDEBUG
  class QRTest;
#endif //MYDEBUG

//! Generates the QR decomposition of [A], which can be used in linear
//! least squares problems and can also be used to find eigenvectors and eigenvalues of [A]
template< template <typename> class Storage = Vect, typename T = double >
class QR
{
public:
  
protected:

private:


#ifdef MYDEBUG
  friend class QRTest;
#endif //MYDEBUG
};

template< typename T >
class QR<Vect,T> : public EigenValSolver<Vect,T>
{
  using EigenValSolver<Vect,T>::m_EigVect;
  using EigenValSolver<Vect,T>::m_EigVal;
  using EigenValSolver<Vect,T>::ExceedsOffDiagTol;
public:
  QR()
  {
    Init();
  }
  
  void Clear()
  {
    Init();
  }
  
  void ATranspose( Matrix<Vect,T>& matrix )
  {
    if ( !matrix.IsTransposed() )
    {
      matrix.Transpose();
    }
    A(matrix);
    m_factorized = false;
  }
  
  void X( Vect<T>& v )
  {
    m_x = &v;
  }
  
  void B( Vect<T>& v )
  {
    m_b = &v;
  }
  
  Matrix<Vect,T>* Q()
  {
    if ( !m_factorized )
    {
      CreateQTr( m_QTr, true );
    }
    return &m_QTr;
  }

  Matrix<Vect,T>* R()
  {
    if ( !m_factorized )
    {
      CreateQTr( m_QTr, true );
    }
    return m_A;
  }
  
  void Factorize()
  {
    Matrix< Vect, T >& A = *m_A;
    
    for ( unsigned long long n = 0; n < A.Cols(); ++n )
    {
      StoreAlpha( n, A );
      for ( unsigned long long j = n + 1; j < A.Cols(); ++j )
      {
        HouseHolderQn_x_( A[j], n );
      }
    }
    CreateQTr( m_QTr );
    m_factorized = true;
  }
  
  void Solve()
  {
    // Ax = b
    // (Q*R)x = b, Q = Q^T = Q^-1
    // R*x = Q^T * b
    // R is upper triangular, back solve up for x
    Matrix< Vect, T >& R = *m_A;
    //Vect<T>& b = *m_b;
    Vect<T>& x = *m_x;
    
    if ( !m_factorized )
    {
      Factorize();
    }
    
    //represents Q^T * b
    vops::Dot<T> dot;
    Vect<T> b( m_b->Size() );
    dot( m_QTr, *m_b, b );

    Solver<Vect,T> solver;
    solver.BackSolveUp( R, x, b ); //backsolving up for x
  }
  
  //! calculates the the eigenvalues and eigenvectors of A (same as A Transpose)
  //! for a non-singular square matrix. If convergence is not met, the eigenvalue
  //! matrix returned will be 2 x 1 of zeros.
  unsigned long long Eig( unsigned long long maxIter, T tol = 1e-5, bool calcEigVect = false )
  {
    Matrix< Vect, T >& R = *m_A;
    Matrix< Vect, T >& A = *m_A;
    Matrix< Vect, T >& evect = m_EigVect;
    Matrix< Vect, T >& Q = m_QTr;
    Matrix<SparseVect,T>& eval = m_EigVal; 
    
    if ( !m_factorized )
    {
      Factorize();
    }
    
    eval.Clear();
    if ( calcEigVect )
    {
      evect.Resize( Q.Rows() );
      //Transpose the eigenvector matrix for easy multiplication
      for ( unsigned long long i = 0; i < Q.Rows(); ++i )
      {
        for ( unsigned long long j = 0; j < Q.Rows(); ++j )
        {
          evect[i][j] = Q[j][i];
        }
      }      
    }

    vops::Dot<T> dot;
    Matrix< Vect, T > AP( Q.Rows() );
    dot( R, Q, AP );
    
    unsigned long long iter = 0;    
    while ( iter < maxIter && ExceedsOffDiagTol( AP, tol) )
    {
      A.Swap( AP );
      Factorize();
      if ( calcEigVect )
      {
        dot( evect, Q, AP ); //use A' as temporary storage of the matrix result
        evect.Swap(AP); //swap result into the eigenvector matrix        
      }
      dot( R, Q, AP );
      ++iter;
    }
    
    if( iter < maxIter )
    {
      eval.Resize( AP.Rows() );
      for ( unsigned long long i = 0; i < AP.Rows(); ++i )
      {
        eval[i].PushBack( i, AP[i][i] );
      }      
    }
    return iter;
  }

protected:
  void A( Matrix<Vect,T>& matrix )
  {
    m_A = &matrix;
    m_alpha.Clear();
    m_alpha.Resize( matrix.Cols() );
    m_factorized = false;
  }
  
  void StoreAlpha( unsigned long long col, Matrix< Vect, T >& A )
  {
    vops::Dot<T> dot;
    m_alpha[col] = sqrt( dot( A[col], A[col], col ) );
    A[col][col] -= m_alpha[col];
  }
  
  void HouseHolderQn_x_( Vect<T>& v, unsigned long long n )
  {
    Matrix< Vect, T >& A = *m_A;
    vops::Dot<T> dot;

    T qu = dot( A[n], v, n );
    T uuRecip = -2.0 / dot( A[n], A[n], n ); //multiply by the -2.0 and combine with reciprocal
    for ( unsigned long long i = n; i < v.Size(); ++i )
    {
      v[i] += qu * A[n][i] * uuRecip;
    }
  }
  
  void SwapAlphaWithRPivots( Matrix< Vect, T >& R, Vect<T>& alpha )
  {
    T tmp;
    for ( unsigned long long i = 0; i < R.Cols(); ++i )
    {
      tmp = R[i][i];
      R[i][i] = alpha[i];
      alpha[i] = tmp;
    }     
  }

  void CreateQTr( Matrix< Vect, T >& QTr, bool prettyR = false )
  {
    Matrix< Vect, T >& A = *m_A;
    //QTr.Clear();
    QTr.Resize( A.Cols(), A.Rows() );
    Solver<Vect,T> solver;
    solver.Identity(QTr);
    
    //creating Q transpose
    for ( unsigned long long i = 0; i < A.Cols(); ++i )
    {
      for ( unsigned long long j = A.Cols(); j > 0; --j ) //cols?
      {
        HouseHolderQn_x_( QTr[i], j - 1 );
      }      
    }
    
    //creating R
    CreateR();
    if ( !QTr.IsTransposed() )
    {
      QTr.Transpose(); //so that we can print Q
    }  
  }
  
  void CreateR()
  {
    Matrix< Vect, T >& R = *m_A;
    //creating R
    if ( R.IsTransposed() )
    {
      R.Transpose(); //just undo the transpose flag
    }
    for ( unsigned long long i = 0; i < R.Rows(); ++i )
    {
      R[i][i] = m_alpha[i];
      for ( unsigned long long j = i + 1; j < R.Rows(); ++j )
      {
        R[i][j] = R[j][i];
        R[j][i] = T();
      }        

    }
    R.Resize( R.Rows() );
    //m_alpha.Clear(); //this is slightly wasteful but plays havoc with Eig if cleared.
  }
  
  void Init()
  {
    m_A = NULL;
    m_b = NULL;
    m_x = NULL;
    m_factorized = false;
  }
  
protected:
  Matrix< Vect, T > m_QTr;
  Matrix< Vect, T >* m_A;
  Vect<T> m_alpha;
  Vect<T>* m_b;
  Vect<T>* m_x;
  bool m_factorized;
private:


#ifdef MYDEBUG
  friend class QRTest;
#endif //MYDEBUG
};

}/*solver*/ }/*math*/ 

#endif /*_MATH_SOLVER_QR_h */
