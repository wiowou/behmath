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

#ifndef _MATH_SOLVER_ELIMINATION_h
#define _MATH_SOLVER_ELIMINATION_h

#include "solver.h"
#include <cmath> //std::abs
#include "axpy.h"
#include "exception.h"
#include "sparseVect.h"
#include "vect.h"

namespace math{
namespace solver{

#ifdef MYDEBUG
  class EliminationTest;
#endif //MYDEBUG

//! Base class for all gaussian elimination based solvers
template< template <typename> class Storage = Vect, typename T = double >
class Elimination : public Solver<Storage,T>
{
  using Solver<Storage,T>::BackSolveUp;
  using Solver<Storage,T>::BackSolveDown;

public:
  Elimination()
  {
    Init();
  }

  void Clear()
  {
    Init();
  }
  
  void A( Matrix< Storage, T >& matrix )
  {
    m_A = &matrix;
    m_PMatrix.Resize( m_A->Rows(), m_A->Cols() );
    for ( ULong i = 0; i < m_A->Rows(); ++i )
    {
      m_PMatrix[i].PushBack(i, 1.0);
    }
    m_factorized = false;
  }
  
  void U( Matrix< Storage, T >& matrix )
  {
    A(matrix);
  }
  
  void L( Matrix< Storage, T >& matrix )
  {
    m_L = &matrix;
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
  
  T Tolerance()
  {
    return m_tolerance;
  }
  
  void Tolerance( T tol )
  {
    m_tolerance = tol;
  }
  
  void PartialPivoting( bool p )
  {
    m_partialPivoting = p;
    m_partialPivotingMax = !p;
  }
  
  bool PartialPivoting()
  {
    return m_partialPivoting;
  }

  void PartialPivotingMax( bool p )
  {
    m_partialPivotingMax = p;
    m_partialPivoting = !p;
  }
  
  bool PartialPivotingMax()
  {
    return m_partialPivotingMax;
  }

  //! Returns a pointer to the Sparse permutation matrix
  Matrix< SparseVect, T >* PMatrix( ) //tested
  {
    return &m_PMatrix;
  }
  
  void UpperTriangular() //tested
  {
    RowReduceDown( *m_A );
  }
  
  void LowerTriangular() //tested
  {
    RowReduceUp( *m_A, m_A->Rows() - 1 );
  }
  
  void Solve()
  {
    UpperTriangular();
    BackSolveUp( *m_A, *m_x, *m_b );
  }
  
protected:
  inline void Reduce( Matrix< Storage, T > &A, const ULong &pivotIdx, T &pivotVal, ULong i ) //tested
  {
    T factor = -A[i][pivotIdx] / pivotVal;
    A[i][pivotIdx] = factor; //store the elimination matrix in A
    // Can use GPU to speed up AXPY
    vops::Axpy<T> AXPY;
    //a start idx was added to AXPY 9/28/2015. This will cause printed matrices
    //to look as though they haven't been reduced. The start idx will save computational time
    //because AXPY need only affect the columns ahead of the pivotIdx, since cells behind it should be 0
    AXPY( factor, A[pivotIdx], A[i], pivotIdx + 1 );
    //AXPY( factor, A[pivotIdx], A[i] );

    if ( m_b != NULL )
    {
      m_b->RowReduce( factor, pivotIdx, i );
    }
  }
 
  void RowReduceDown( Matrix< Storage, T > &A, ULong pivotIdx = 0 ) //tested
  {
    if ( pivotIdx == A.Rows() )
    {
      return;
    }
    
    T pivotVal = A[pivotIdx][pivotIdx];
    ULong newPivotIdx = PivotDown( A, pivotIdx, pivotVal ); //pivots according to options set
    RowSwap( pivotIdx, newPivotIdx );
    
    // This loop can be parallelized
    for ( ULong i = pivotIdx + 1; i < A.Rows(); ++i )
    {
      Reduce( A, pivotIdx, pivotVal, i );
    }
    RowReduceDown( A, pivotIdx + 1 ); //proper tail recursion
  }

  void RowReduceUp( Matrix< Storage, T > &A, ULong pivotIdx ) //tested
  {
    if ( pivotIdx == 0 )
    {
      return;
    }
    
    T pivotVal;
    ULong newPivotIdx = PivotUp( A, pivotIdx, pivotVal ); //pivots according to options set
    RowSwap( pivotIdx, newPivotIdx );

    vops::Axpy<T> AXPY;
    T factor;
    
    // This loop can be parallelized
    for ( ULong i = pivotIdx; i > 0; --i )
    {
      Reduce( A, pivotIdx, pivotVal, i - 1 );
    }
    RowReduceUp( A, pivotIdx - 1 ); //proper tail recursion
  }
  
  ULong PivotDown( Matrix< Storage, T > &A, const ULong pivotIdx, T &pivotVal )
  {
    if ( m_partialPivotingMax )
    {
      return PartialPivotingMaxDown( A, pivotIdx, pivotVal );
    }
    else if ( m_partialPivoting )
    {
      return PartialPivotingDown( A, pivotIdx, pivotVal );
    }   
  }

  ULong PivotUp( Matrix< Storage, T > &A, const ULong pivotIdx, T &pivotVal )
  {
    if ( m_partialPivotingMax )
    {
      return PartialPivotingMaxUp( A, pivotIdx, pivotVal );
    }
    else if ( m_partialPivoting )
    {
      return PartialPivotingUp( A, pivotIdx, pivotVal );
    }    
  }
  
  ULong PartialPivotingDown( Matrix< Storage, T > &A, const ULong pivotIdx, T &pivotVal )
  { 
    ULong newIdx = pivotIdx;
    while ( newIdx < A.Rows() && std::abs(pivotVal) < m_tolerance )
    {
      ++newIdx;
      pivotVal = A[newIdx][pivotIdx];
    }
    if ( newIdx == A.Rows() )
    {
      throw Singular();
    }
    return newIdx;
  }

  ULong PartialPivotingUp( Matrix< Storage, T > &A, const ULong pivotIdx, T &pivotVal )
  { 
    ULong newIdx = pivotIdx;
    while ( newIdx > 0 && std::abs(pivotVal) < m_tolerance )
    {
      --newIdx;
      pivotVal = A[newIdx][pivotIdx];
    }
    if ( newIdx == 0 )
    {
      if ( std::abs(pivotVal) < m_tolerance )
      {
        throw Singular();
      }
      pivotVal = A[newIdx][pivotIdx];
    }
    return newIdx;
  }
  
  ULong PartialPivotingMaxDown( Matrix< Storage, T > &A, const ULong pivotIdx, T &pivotVal )
  {
    ULong maxIdx = A.MaxIdxBelow( pivotIdx );
    pivotVal = A[maxIdx][pivotIdx];
    if ( std::abs(pivotVal) < Tolerance() )
    {
      throw Singular();
    }
    return maxIdx;
  }

  ULong PartialPivotingMaxUp( Matrix< Storage, T > &A, const ULong pivotIdx, T &pivotVal )
  {
    ULong maxIdx = A.MaxIdxAbove( pivotIdx );
    pivotVal = A[maxIdx][pivotIdx];
    if ( std::abs(pivotVal) < Tolerance() )
    {
      throw Singular();
    }
    return maxIdx;
  }
  
  void RowSwap( ULong pivotIdx, ULong swapIdx ) //tested
  {
    m_A->RowSwap( pivotIdx, swapIdx );
    m_PMatrix.RowSwap( pivotIdx, swapIdx );
    if ( m_b != NULL )
    {
      m_b->RowSwap( pivotIdx, swapIdx );
    }
    if ( m_L != NULL )
    {
      m_L->RowSwap( pivotIdx, swapIdx );
    }
  }

  void ReOrder( Vect<T>& v )
  {
    Vect<T> vNew( v.Size() );
    for ( ULong i = 0; i < v.Size(); ++i )
    {
      vNew[i] = v[ m_PMatrix[i].Pos(0) ];
    }
    v.Swap(vNew);
  }
  
  void Init()
  {
    m_A = NULL;
    m_L = NULL;
    m_b = NULL;
    m_x = NULL;
    m_PMatrix.Clear();
    m_tolerance = 0.001;
    m_partialPivoting = true;
    m_partialPivotingMax = false;
    m_factorized = false;
  }  
  
protected:
  Matrix< Storage, T >* m_A;
  Matrix< Storage, T >* m_L;
  Vect<T>* m_b;
  Vect<T>* m_x;
  Matrix< SparseVect, T > m_PMatrix;
  T m_tolerance;
  bool m_partialPivoting;
  bool m_partialPivotingMax;
  bool m_factorized;
  
private:


#ifdef MYDEBUG
  friend class EliminationTest;
#endif //MYDEBUG
};

}/*solver*/ }/*math*/ 

#endif /*_MATH_SOLVER_ELIMINATION_h */
