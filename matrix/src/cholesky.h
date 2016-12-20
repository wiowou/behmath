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

#ifndef _MATH_SOLVER_CHOLESKY_h
#define _MATH_SOLVER_CHOLESKY_h

#include "solver.h"
#include "vect.h"
#include "sparseVect.h"
#include "dot.h"
#include <math.h> //for sqrt

namespace math{
namespace solver{

#ifdef MYDEBUG
  class CholeskyTest;
#endif //MYDEBUG

//! Cholesky decomposition and solve.
template< template <typename> class Storage = Vect, typename T = double >
class CholeskyBase
{

public:
  CholeskyBase()
  {
    Init();
  }
  
  void Clear()
  {
    Init();
  }
  
  void L( Matrix<Storage,T>& matrix )
  {
    m_L = &matrix;
    m_factorized = false;
  }

  void A( Matrix<Storage,T>& matrix )
  {
    L(matrix);
  }
  
  void X( Vect<T>& v )
  {
    m_x = &v;
  }
  
  void B( Vect<T>& v )
  {
    m_b = &v;
  }

  void Solve()
  {
    //L(LT)*x = b
    //L(LT*x) = b, LT*x = c
    //L*c = b <- solve for c first
    //LT*x = c <- solve for x next  
    
    Solver<Storage,T> solver;
    //This re-naming should clarify the application of the algorithm
    Matrix<Storage, T >& L = *m_L;
    Vect<T>& b = *m_b;
    Vect<T>& c = *m_x;
    
    solver.BackSolveDown( L, c, b ); //solving for c
    Vect<T>& x = *m_b; //re-naming so we can reuse memory rather than having to reallocate
    if ( !L.IsTransposed() )
    {
      L.Transpose();
    }
    solver.BackSolveUp( L, x, c ); //solving for x
    m_x->Swap( x ); //putting x where the user expects to find it
  }
protected:
  void Init()
  {
    m_L = NULL;
    m_b = NULL;
    m_x = NULL;
    m_factorized = false;
  }
protected:
  Matrix< Storage, T >* m_L;
  Vect<T>* m_b;
  Vect<T>* m_x;
  bool m_factorized;
private:


#ifdef MYDEBUG
  friend class CholeskyTest;
#endif //MYDEBUG
};

//! Cholesky decomposition and solve.
template< template <typename> class Storage = Vect, typename T = double >
class Cholesky
{

public:
  
protected:

private:


#ifdef MYDEBUG
  friend class CholeskyTest;
#endif //MYDEBUG
};

//! Cholesky decomposition and solve. The SparseVect matrix only needs to define the lower half
//! of the symmetric matrix. If the matrix is banded, only the lower half of the bands need
//! to be specified. 
template< typename T >
class Cholesky<SparseVect, T> : public CholeskyBase<SparseVect,T>
{
  using CholeskyBase<SparseVect,T>::m_L;
  using CholeskyBase<SparseVect,T>::m_b;
  using CholeskyBase<SparseVect,T>::m_x;
  using CholeskyBase<SparseVect,T>::m_factorized;

public:
  using CholeskyBase<SparseVect,T>::L;
  using CholeskyBase<SparseVect,T>::X;
  using CholeskyBase<SparseVect,T>::A;
  using CholeskyBase<SparseVect,T>::B;
  //using CholeskyBase<SparseVect,T>::Solve;

  void Solve()
  {
    if ( !m_factorized )
    {
      Factorize();
    }
    CholeskyBase<SparseVect,T>::Solve();
  }
  
  void Factorize()
  {
    Matrix< SparseVect, T >& L = *m_L;
    vops::Dot<T> dot;
    
    L[0].Back() = sqrt( L[0].Back() );
    for ( unsigned long long i = 1; i < L.Rows(); ++i )
    {
      for ( unsigned long long j = L[i].Pos(0), k = 0; j < i + 1; ++j )
      {
        if ( k == L[i].NNZ() )
        {
          break;
        }
        if ( L[i].Pos(k) != j )
        {
          continue;
        }
        if ( j == i )
        {
          T d = dot( j, L[j], L[j] );
          L[j].Back() = sqrt( L[j].Back() - d );
        }
        else if ( L[i].Pos(k) == 0 )
        {
          L[i].Front() = L[i].Front() / L[0].Front();
        }
        else
        {
          T d = dot( j, L[i], L[j] );
          L[i](k) = ( L[i](k) - d )/ L[j].Back();
        }
        ++k;
      }
    }
    m_factorized = true;
  }
  

protected:

protected:


private:


#ifdef MYDEBUG
  friend class CholeskyTest;
#endif //MYDEBUG
};

//! Cholesky decomposition and solve.
//! The dense matrix need only specify the bottom half of the symmetric matrix
template< typename T >
class Cholesky<Vect, T> : public CholeskyBase<Vect,T>
{

  using CholeskyBase<Vect,T>::m_L;
  using CholeskyBase<Vect,T>::m_b;
  using CholeskyBase<Vect,T>::m_x;
  using CholeskyBase<Vect,T>::m_factorized;
public:
  using CholeskyBase<Vect,T>::L;
  using CholeskyBase<Vect,T>::X;
  using CholeskyBase<Vect,T>::A;
  using CholeskyBase<Vect,T>::B;
  //using CholeskyBase<Vect,T>::Solve;
  
  void Solve()
  {
    if ( !m_factorized )
    {
      Factorize();
    }
    CholeskyBase<Vect,T>::Solve();
  }
  
  void Factorize() //tested
  {
    Matrix< Vect, T >& L = *m_L;
    vops::Dot<T> dot;
    
    L[0][0] = sqrt( L[0][0] );
    for ( unsigned long long i = 1; i < L.Rows(); ++i )
    {
      for ( unsigned long long j = 0; j < i + 1; ++j )
      {
        if ( j == i )
        {
          T d = dot( j, L[j], L[j] );
          L[j][j] = sqrt( L[j][j] - d );
        }
        else if ( j == 0 )
        {
          L[i][0] = L[i][0] / L[0][0];
        }
        else
        {
          T d = dot( j, L[i], L[j] );
          L[i][j] = ( L[i][j] - d )/ L[j][j];
        }
      }
    }
    m_factorized = true;
  }
  

protected:

protected:


private:


#ifdef MYDEBUG
  friend class CholeskyTest;
#endif //MYDEBUG
};

}/*solver*/ }/*math*/ 

#endif /*_MATH_SOLVER_CHOLESKY_h */
