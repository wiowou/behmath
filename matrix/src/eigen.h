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

#ifndef _MATH_SOLVER_EIGEN_h
#define _MATH_SOLVER_EIGEN_h

#include "matrix.h"
#include "vect.h"
#include "sparseVect.h"
#include "elimination.h"
#include "exception.h"

namespace math{
namespace solver{

#ifdef MYDEBUG
  class EigenTest;
#endif //MYDEBUG

template< template< template <typename> class Storage, typename T > class Sol, template <typename> class Storage = SparseVect, typename T = double >
class Eigen
{

public:
  Eigen()
  {
    Init();
  }
  
  void Clear()
  {
    Init();
  }
  
  void ATranspose( Matrix<Storage,T>& matrix )
  {
    m_A = &matrix;
    m_calcedEigVect = false;
    m_calcedEigVal = false;
    m_needToPivot.Resize( m_A->Rows() );
    for ( ULong i = 0; i < m_A->Rows(); ++i )
    {
      m_needToPivot[i] = true;
    }
  }
  
  void EigVal( Matrix<SparseVect,T>& eigVal )
  {
    m_EigVal = &eigVal;
    m_calcedEigVal = true;
    m_calcedEigVect = false;
  }
  
  void Tolerance( T tol )
  {
    m_tolerance = tol;
  }
  
  void MaxIter( ULong maxIter )
  {
    m_maxIter = maxIter;
  }
  
  Matrix<Vect,T>* EigVect()
  {
    if ( !m_calcedEigVal )
    {
      CalcEigVal();
    }
    CalcEigVect();
    if ( !m_EigVect.IsTransposed() )
    {
      m_EigVect.Transpose();
    }
    return &m_EigVect;
  }

  Matrix<SparseVect,T>* EigVal()
  {
    CalcEigVal();
    return &m_EigVal;
  }
  
  Matrix<Vect,T>* VT()
  {
    return &m_VT;
  }

  Matrix<Vect,T>* U()
  {
    if ( !m_U.IsTransposed() )
    {
      m_U.Transpose();
    }
    return &m_U;
  }
  
  Matrix<SparseVect,T>* S()
  {
    return &m_EigVal;
  }
  
  ULong Iter()
  {
    return m_iter;
  }
  
  void SVD()
  {
    //undo any transposes in EigVect
    if ( m_EigVect.IsTransposed() )
    {
      m_EigVect.Transpose();
    }
    
    // m_A represents A transpose
    Matrix<Storage,T>& AT = *m_A;
    
    //use the physical dimensions since not sure whether transpose flag is set
    ULong nrow = AT.PRows();
    ULong ncol = AT.PCols();
    
    //save the location of m_A so it can be restored
    Matrix<Storage,T>* pA = m_A;
    vops::Dot<T> dot;
    
    //create A*A^T by creating A first then dotting it
    Matrix<Vect,T> A(ncol, nrow);  //mem alloc

    for ( ULong i = 0; i < nrow; ++i )
    {
      for ( ULong j = 0; j < ncol; ++j )
      {
        A[j][i] = AT[i][j];
      }
    }
    
    Matrix<Vect,T> AAT(ncol);  //mem alloc
    dot( A, A, AAT );
    
    ATranspose(AAT); //feed AAT to solver
    CalcEigVal(); //ask solver to calculate eigen values
    CalcEigVect(); //ask solver to calculate eigen vects
    m_U.Swap( m_EigVect ); //store the eig vects in U
    
    //create A^T*A
    Matrix<Vect,T>& ATA = AAT; //reuse memory from AAT
    ATA.Resize( nrow );
    dot( AT, AT, ATA );
    
    //find its eigenvals and eigenvects
    ATranspose( ATA );
    CalcEigVect();
    m_VT.Swap( m_EigVect ); //store the eig vects in VT

    //sqrt of eigvals of AAT or ATA = eigvals of A or AT
    for ( ULong i = 0; i < nrow; ++i )
    {
      m_EigVal[i][i] = sqrt( m_EigVal[i][i] );
    }

    ATranspose( *pA ); //restoring some of the internal state

  }
  
protected:
  void CalcEigVal()
  {
    if ( m_calcedEigVal )
    {
      return;
    }
    
    Sol<Storage,T> solver;
    Matrix<Storage,T> A = *m_A; //making a copy of A so that eigenvector extraction works later on. Price to pay for performance.
    solver.ATranspose(A);
    m_iter = solver.Eig( m_maxIter, m_tolerance, false ); //don't let the solver calculate the eigenvectors
    
    m_calcedEigVal = true;
    Matrix<SparseVect,T>* tmp = solver.EigVal();
    m_EigVal.Swap(*tmp);
  }
  
  //! The eigenvectors will come out in the physical rows of the matrix. Therefore,
  //! the m_EigVect / evect matrix is transposed for printing purposes.
  void CalcEigVect()
  {
    if ( m_calcedEigVect )
    {
      return;
    }
    
    Matrix<SparseVect,T>& eval = m_EigVal;
    Matrix<Vect,T>& evect = m_EigVect;
    ULong nrow = m_A->Rows();
    
    Vect<T> b( nrow - 1 ); //size b for A - last row

    evect.Resize( nrow );
    
    for ( ULong i = 0; i < nrow; ++i )
    {
      evect[i][nrow - 1] = T(1.0);
      Elimination<Storage,T> sol;
      Matrix<Storage,T> A = *m_A; //copy values from m_A since gauss solver changes A
      
      Shift( A, eval[i][i] ); //shift A
  
      ULong knownRow = nrow;
      while ( knownRow > 0 )
      {
        A.RowSwap( nrow - 1, knownRow - 1 );
        for ( ULong j = 0; j < b.Size(); ++j )
        {
          b[j] = -A[j][ nrow - 1 ];
        }
        
        bool breakLoop = true;
        try
        {
          Elimination<Storage,T> solver;
          A.Resize( nrow - 1 );
          evect[i].Resize( nrow - 1 );
          solver.A( A );
          solver.B( b );
          
          solver.X( evect[i] );
          solver.PartialPivoting(true);
          solver.Tolerance(1e-5);
          solver.Solve();
          A.Resize( nrow ); //make sure this doesn't delete the data in A
          evect[i].Resize( nrow );
        }
        catch ( Singular &err )
        {
          breakLoop = false; 
        }
        if ( breakLoop )
        {
          break;
        }
        --knownRow;
      }
      if ( knownRow == 0 )
      {
        for ( ULong j = 0; j < nrow - 1; ++j )
        {
          evect[i][j] = T();
        }
      }
    }
    m_calcedEigVect = true;
  }
  
  void Shift( Matrix<Storage,T>& A, T val )
  {
    for ( ULong i = 0; i < A.Rows(); ++i )
    {
      A[i][i] -= val;
    }
  }

  void Init()
  {
    m_EigVect.Clear();
    //m_EigVect.Transpose();
    m_A = NULL;
    m_EigVal.Clear();
    m_calcedEigVect = false;
    m_calcedEigVal = false;
    m_maxIter = 0;
    m_tolerance = 1e-5;
    m_iter = 0;
    m_VT.Clear();
    m_U.Clear();
  }
  
protected:
  Matrix<Storage,T>* m_A;
  Matrix<SparseVect,T> m_EigVal;
  Matrix<Vect,T> m_EigVect;
  ULong m_maxIter;
  T m_tolerance;
  bool m_calcedEigVect;
  bool m_calcedEigVal;
  Vect<bool> m_needToPivot;
  ULong m_iter;
  
  Matrix<Vect,T> m_VT;
  Matrix<Vect,T> m_U;
private:


#ifdef MYDEBUG
  friend class EigenTest;
#endif //MYDEBUG
};

}/*solver*/ }/*math*/ 

#endif /*_MATH_SOLVER_EIGEN_h */
