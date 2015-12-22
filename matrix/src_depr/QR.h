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

/*
This implementation is unique in that it turns Q into a square matrix while
keeping R as rectangular. It's the same Q, but with extra columns. Same R,
but with zeros as extra rows.
*/

#ifndef _MATH_SOLVER_QR_h
#define _MATH_SOLVER_QR_h

#include "solver.h"
#include "matrix.h"
#include "dot.h"
#include <math.h> //for sqrt

namespace math{
namespace solver{

#ifdef MYDEBUG
  class QRTest;
#endif //MYDEBUG

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
class QR<Vect,T>
{

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
    matrix.Transpose();
    A(matrix);
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
    if ( !m_QRcalculated )
    {
      CreateQ( m_Q );
    }
    return &m_Q;
  }

  Matrix<Vect,T>* R()
  {
    if ( !m_QRcalculated )
    {
      CreateQ( m_Q );
    }
    return m_A;
  }
  
  void Factorize()
  {
    Matrix< Vect, T >& A = *m_A;
    
    for ( ULong n = 0; n < A.Cols(); ++n )
    {
      StoreAlpha( n, A );
      for ( ULong j = n + 1; j < A.Cols(); ++j )
      {
        HouseHolderQn_x_( A[j], n );
      }
    }
    m_factorized = true;
  }
  
  void Solve()
  {
    // Ax = b
    // (Q*R)x = b, Q = Q^T = Q^-1
    // R*x = Q^T * b
    // R is upper triangular, back solve up for x
    Matrix< Vect, T >& R = *m_A;
    Vect<T>& b = *m_b;
    Vect<T>& x = *m_x;
    
    if ( !m_factorized )
    {
      Factorize();
    }
    
    //multiplying the HouseHolder Qn's by b in reverse order
    //represents Q^T * b
    for ( ULong i = R.Cols(); i > 0; --i )
    {
      HouseHolderQn_x_( b, i - 1 );
    }

    SwapAlphaWithRPivots( R, m_alpha ); //swap once
    Solver<Vect,T> solver;
    solver.BackSolveUp( R, x, b ); //backsolving up for x
    SwapAlphaWithRPivots( R, m_alpha ); //swap back
  }

protected:
  void A( Matrix<Vect,T>& matrix )
  {
    m_A = &matrix;
    m_alpha.Resize( matrix.Cols() );
    m_QRcalculated = false;
    m_factorized = false;
  }
  
  void StoreAlpha( ULong col, Matrix< Vect, T >& A )
  {
    vops::Dot<T> dot;
    m_alpha[col] = sqrt( dot( A[col], A[col], col ) );
    A[col][col] -= m_alpha[col];
  }
  
  void HouseHolderQn_x_( Vect<T>& v, ULong n )
  {
    Matrix< Vect, T >& A = *m_A;
    vops::Dot<T> dot;

    T qu = dot( A[n], v, n );
    T uuRecip = -2.0 / dot( A[n], A[n], n ); //multiply by the -2.0 and combine with reciprocal
    for ( ULong i = n; i < v.Size(); ++i )
    {
      v[i] += qu * A[n][i] * uuRecip;
    }
  }
  
  void SwapAlphaWithRPivots( Matrix< Vect, T >& R, Vect<T>& alpha )
  {
    T tmp;
    for ( ULong i = 0; i < R.Cols(); ++i )
    {
      tmp = R[i][i];
      R[i][i] = alpha[i];
      alpha[i] = tmp;
    }     
  }

  void CreateQ( Matrix< Vect, T >& Q )
  {
    Matrix< Vect, T >& A = *m_A;
    
    Q.Resize( A.Rows() );
    Solver<Vect,T> solver;
    solver.Identity(Q);
    
    //creating Q transpose
    for ( ULong i = 0; i < A.Rows(); ++i )
    {
      for ( ULong j = A.Cols(); j > 0; --j )
      {
        HouseHolderQn_x_( Q[i], j - 1);
      }      
    }
    
    //creating R
    for ( ULong i = 0; i < A.Cols(); ++i )
    {
      A[i][i] = m_alpha[i];
      for ( ULong j = i + 1; j < A.Rows(); ++j )
      {
        A[i][j] = T();
      }
    }
    Q.Transpose();
    m_QRcalculated = true;
  }
  
  void Init()
  {
    m_A = NULL;
    m_b = NULL;
    m_x = NULL;
    m_QRcalculated = false;
    m_factorized = false;
  }
  
protected:
  Matrix< Vect, T > m_Q;
  Matrix< Vect, T >* m_A;
  Vect<T> m_alpha;
  Vect<T>* m_b;
  Vect<T>* m_x;
  bool m_QRcalculated;
  bool m_factorized;
private:


#ifdef MYDEBUG
  friend class QRTest;
#endif //MYDEBUG
};

}/*solver*/ }/*math*/ 

#endif /*_MATH_SOLVER_QR_h */
