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

#ifndef _MATH_SOLVER_LU_h
#define _MATH_SOLVER_LU_h

#include "elimination.h"
#include "solver.h"

namespace math{
namespace solver{

#ifdef MYDEBUG
  class LUTest;
#endif //MYDEBUG

template< template <typename> class Storage = Vect, typename T = double >
class LU
{

public:
  
protected:

private:


#ifdef MYDEBUG
  friend class LUTest;
#endif //MYDEBUG
};


template< typename T >
class LU<Vect,T> : public Elimination<Vect,T >
{
  using Elimination<Vect,T>::m_A;
  using Elimination<Vect,T>::m_b;
  using Elimination<Vect,T>::m_x;
  using Elimination<Vect,T>::m_factorized;
  using Elimination<Vect,T>::UpperTriangular;
  using Elimination<Vect,T>::ReOrder;
  
public:
  void Factorize()
  {
    UpperTriangular();
    m_factorized = true;
    return;
  }
  
  //! Creates the Inverse of [A], Transposed. LU is used to do a row by row
  //! solve on the identity matrix. [A] must have been provided to the class prior to 
  //! using this method.
  void InverseTranspose( Matrix<Vect,T> &InvT )
  {
    InvT.Resize( m_A->Cols() );
    Solver<Vect,T> solver;
    solver.Identity(InvT);
    Vect<T> x( m_A->Cols() );
    m_x = &x;
    if ( !m_factorized )
    {
      Factorize();
    }
    for ( ULong i = 0; i < m_A->Cols(); ++i )
    {
      m_b = &InvT[i];
      Solve();
      x.Swap( InvT[i] );
    }
  }

  void Solve()
  {
    //LU*x = b
    //L(U*x) = b, U*x = c
    //L*c = b <- solve for c first
    //U*x = c <- solve for x next  
    if ( !m_factorized )
    {
      Factorize();
    }
    Solver<Vect,T> solver;
    //This re-naming should clarify the application of the algorithm
    Matrix< Vect, T >& L = *m_A;
    Matrix< Vect, T >& U = *m_A;
    Vect<T>& b = *m_b;
    Vect<T>& c = *m_b;
    Vect<T>& x = *m_x;
    
    ReOrder( b ); //use permutation matrix to re order b
    LInvB( L, b ); //multiply L^-1 * b and store in b
    
    solver.BackSolveUp( U, x, c );  //solving for x
    m_x->Swap( x ); //putting x where the user expects to find it
  }  
protected:
  void LInvB( Matrix< Vect, T >& L, Vect<T>& b )
  {
    //go up the column, then back a column
    for ( ULong n = 0; n < L.Cols(); ++n )
    {
      for ( ULong m = n + 1; m < L.Rows(); ++m )
      {
        b[m] += L[m][n] * b[n];
      }
    }
  }
private:


#ifdef MYDEBUG
  friend class LUTest;
#endif //MYDEBUG
};

template< typename T >
class LU<SparseVect,T> : public Elimination<SparseVect,T >
{
  using Elimination<SparseVect,T>::m_A;
  using Elimination<SparseVect,T>::m_L;
  using Elimination<SparseVect,T>::m_b;
  using Elimination<SparseVect,T>::m_x;
  using Elimination<SparseVect,T>::m_factorized;
  using Elimination<SparseVect,T>::UpperTriangular;
  using Elimination<SparseVect,T>::ReOrder;
  
public:
  void Factorize( Matrix< SparseVect, T >& L, Matrix< SparseVect, T >& U )
  {
    L = U;
    m_L = &L;
    Elimination<SparseVect,T>::U( U );
    UpperTriangular();
    Matrix< SparseVect, T >& A = *m_A;
    
    for ( ULong i = 0; i < L.Rows(); ++i )
    {
      L[i].Set(i, 1.0 );
      for ( ULong j = i + 1; j < L.Cols(); ++j )
      {
         L[i][j] = T();
      }
    }
    Vect<T> colU( L.Rows() );
    vops::Dot<T> dot;
    for ( ULong i = 0; i < L.Cols(); ++i )
    {
      for ( ULong k = 0; k < i + 1; ++k )
      {
        colU[k] = A[k][i];
      }
      for ( ULong j = i + 1; j <  L.Rows(); ++j )
      {
        T d = dot( i, colU, L[j] );
        L[j][i] = ( L[j][i] - d ) / colU[i];
      }      
    }
    m_factorized = true;
    return;
  }
  
  //! Creates the Inverse of [A], Transposed. LU is used to do a row by row
  //! solve on the identity matrix. [A] must have been provided to the class prior to 
  //! using this method.
  void InverseTranspose( Matrix<SparseVect,T> &InvT )
  {
    InvT.Resize( m_A->Cols() );
    Solver<SparseVect,T> solver;
    solver.Identity(InvT);
    Vect<T> x( m_A->Cols() );
    m_x = &x;
    if ( !m_factorized )
    {
      Factorize( *m_L, *m_A );
    }
    for ( ULong i = 0; i < m_A->Cols(); ++i )
    {
      m_b = &InvT[i];
      Solve();
      x.Swap( InvT[i] );
    }
  }

  void Solve()
  {
    //LU*x = b
    //L(U*x) = b, U*x = c
    //L*c = b <- solve for c first
    //U*x = c <- solve for x next  
    if ( !m_factorized )
    {
      Factorize( *m_L, *m_A );
    }
    Solver<SparseVect,T> solver;
    //This re-naming should clarify the application of the algorithm
    Matrix< SparseVect, T >& L = *m_L;
    Matrix< SparseVect, T >& U = *m_A;
    Vect<T>& b = *m_b;
    Vect<T>& c = *m_x;
    
    ReOrder( b ); //use permutation matrix to re order b
    solver.BackSolveDown( L, c, b ); //solving for c
    Vect<T>& x = *m_b; //re-naming so we can reuse memory rather than having to reallocate
    solver.BackSolveUp( U, x, c );  //solving for x
    m_x->Swap( x ); //putting x where the user expects to find it
  }  
protected:

private:


#ifdef MYDEBUG
  friend class LUTest;
#endif //MYDEBUG
};

}/*solver*/ }/*math*/ 

#endif /*_MATH_SOLVER_LU_h */
