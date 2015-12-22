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

#ifndef _MATH_RMAJ_MATRIX_h
#define _MATH_RMAJ_MATRIX_h

#include "impl/matrixConfig.h"
#include "vect.h"
#include "sparseVect.h"
#include <iostream>
#include <iomanip>
#include "typedefs.h"

namespace math{

#ifdef MYDEBUG
  class MatrixTest;
#endif //MYDEBUG


//! A mathematical 2-variable matrix. Will always have dimension of at least 2 rows,
//! even after a Clear().
template< template <typename> class Storage = Vect, typename T = double >
class Matrix
{
public:
  Matrix()
  {
    Init(2,1);
  }
  
  explicit Matrix( ULong size )
  {
    Init(size);
  }
  
  Matrix( ULong mrow, ULong ncol  )
  {
    Init(mrow, ncol);
  }
  
  ~Matrix()
  {
    Destroy();
  }

  //! copy constructor
  Matrix( const Matrix &other )
  {
    m_size = other.m_size;
    m_isTransposed = other.m_isTransposed;
    m_row = new Storage<T> [m_size];
    for ( ULong i = 0; i < m_size; ++i )
    {
      m_row[i] = other.m_row[i];
    }
  }
  
  //! assignment operator
  Matrix& operator=( Matrix other )
  {
    Swap(other);
    return *this;
  }
  
  //! Returns number of rows of the matrix, taking account of transposition
  ULong Rows()
  {
    if ( m_isTransposed )
    {
      return m_row[0].Size();
    }
    return m_size;
  }

  //! Returns number of physical rows of the matrix, not taking account of transposition
  ULong PRows()
  {
    return m_size;
  }
  
  //! Returns number of cols of the matrix, taking account of transposition
  ULong Cols()
  {
    if ( m_isTransposed )
    {
      return m_size;
    }
    return m_row[0].Size();
  }

  //! Returns number of cols of the matrix, not taking account of transposition
  ULong PCols()
  {
    return m_row[0].Size();
  }
  
  void Swap( Matrix &other )
  {
    bool tmp_isTransposed = other.m_isTransposed;
    ULong tmp_size = other.m_size;
    Storage<T>* tmp_row = other.m_row;
    
    other.m_isTransposed = m_isTransposed;
    other.m_size = m_size;
    other.m_row = m_row;
    
    m_isTransposed = tmp_isTransposed;
    m_size = tmp_size;
    m_row = tmp_row;
  }
  
  //! Re-dimensions the matrix and its child row vectors. Deletes previous data contained within.
  //! This method assumes a square matrix
  void Resize( ULong size )
  {
    if ( size <= m_size )
    {
      m_size = size;
      for ( ULong i = 0; i < m_size; ++i )
      {
        m_row[i].Resize(m_size);
      }
      return;
    }    
    Clear();
    Init(size);
  }
  
  //! Re-dimensions the matrix and its child row vectors. Deletes previous data contained within.
  void Resize( ULong mrow, ULong ncol )
  {
    m_isTransposed = false;
    if ( mrow <= m_size )
    {
      m_size = mrow;
      for ( ULong i = 0; i < m_size; ++i )
      {
        m_row[i].Resize(ncol);
      }
      return;
    }

    Storage<T>* row = new Storage<T> [mrow];  
    for ( ULong i = 0; i < m_size; ++i )
    {
      row[i] = m_row[i];
    }
    m_size = mrow;
    Destroy();
    m_row = row;
    
    for ( ULong i = 0; i < m_size; ++i )
    {
      m_row[i].Resize(ncol);
    }
  }
  
  void Clear()
  {
    Destroy();
    Init(2,1);
  }
  
  //! Sets data for a dense matrix. Data is read until nrow*ncol is reached.
  //! If matrix is sparse, fill data in using the vector's Data method: matrix[i].Data()
  void Data( T* data ) //tested
  {
    ULong nrow = m_size;
    ULong ncol = m_row[0].Size();
    for ( ULong i = 0; i < nrow; ++i )
    {
      for ( ULong j = 0; j < ncol; ++j )
      {
        Set( i,j, *data++ );
      }
    }
  }
  
  //! Returns a row vector
  Storage<T>& operator[]( ULong m )
  {
    return m_row[m];
  }
  
  //! Gets/Sets the m,n element of the matrix. Recognizes if matrix is transposed.
  T& operator()( ULong m, ULong n ) //tested
  {
    if ( m_isTransposed )
    {
      return m_row[n][m];
    }
    return m_row[m][n];
  }
  
  //! Sets the m,n element of the matrix. Recognizes if matrix is transposed.
  void Set( ULong m, ULong n, T val ) //tested
  {
    if ( m_isTransposed )
    {
      m_row[n].Set(m, val);
      return;
    }
    m_row[m].Set(n, val);
    return;
  }
  
  void Transpose() //tested
  {
    m_isTransposed = !m_isTransposed;
  }
  
  bool IsTransposed()
  {
    return m_isTransposed;
  }
  
  void RowSwap( ULong rowNum1, ULong rowNum2 ) //tested
  {
    m_row[rowNum1].Swap( m_row[rowNum2] );
  }

  //! Returns the row number for the maximum pivot value below @param pivotIdx row
  ULong MaxIdxBelow( const ULong pivotIdx ) const //tested
  {
    ULong rowEnd = m_size;
    ULong max = pivotIdx;
    T maxVal = m_row[pivotIdx][pivotIdx];
    for ( ULong i = pivotIdx + 1; i < rowEnd; ++i )
    {   
      T val = m_row[i][pivotIdx];
      if ( val * val > maxVal * maxVal  )
      {
        maxVal = val;
        max = i;
      }
    }
    return max;
  }
  
  //! Returns the row number for the maximum pivot value above @param pivotIdx row
  ULong MaxIdxAbove( const ULong pivotIdx ) const //tested
  {
    ULong rowEnd = 0;
    ULong max = pivotIdx;
    T maxVal = m_row[pivotIdx][pivotIdx];
    for ( ULong i = pivotIdx; i > rowEnd; --i )
    {   
      T val = m_row[i - 1][pivotIdx];
      if ( val * val > maxVal * maxVal )
      {
        maxVal = val;
        max = i - 1;
      }
    }
    return max;
  }
  
  //! Print the matrix
  friend std::ostream& operator<<( std::ostream& os, Matrix& A )
  {
    ULong nrow = A.m_size;
    ULong ncol = A[0].Size();
    if ( A.IsTransposed() )
    {
      nrow = A[0].Size();
      ncol = A.m_size;
    }
    std::ios::fmtflags fla(os.flags() );
    os << "Matrix( ";
    os << nrow << "x" << ncol << " )" << " = " << std::endl;
    
    for (ULong i = 0; i < nrow; ++i )
    {
      os << "  [";
      for (ULong j = 0; j < ncol; ++j )
      {
        os << std::fixed << std::setw(10) << std::setprecision(4);
        if ( A.IsTransposed() )
        {
          os << A[j][i] << " ";
        }
        else
        {
          os << A[i][j] << " ";
        }
        os.flags(fla);
      }
      os << " ]" << std::endl;
    }
    return os;
  }

protected: //methods
  
  void Init( ULong size )
  {
    m_isTransposed = false;
    size = size < 2 ? 2 : size;
    m_row = new Storage<T> [size];
    m_size = size;
    
    for ( ULong i = 0; i < m_size; ++i )
    {
      m_row[i].Resize(m_size);
    }
  }
  
  void Init( ULong mrow, ULong ncol )
  { 
    m_isTransposed = false;
    mrow = mrow < 2 ? 2 : mrow;
    m_row = new Storage<T> [mrow];
    m_size = mrow;
    
    for ( ULong i = 0; i < m_size; ++i )
    {
      m_row[i].Resize(ncol);
    }
  }
  
  void Destroy()
  {
    delete[] m_row;
  }
  
protected: //members
  bool m_isTransposed;
  
  //! number of rows in the matrix
  ULong m_size;
  
  //! pointer to a row vector
  Storage<T>* m_row;

private:

#ifdef MYDEBUG
  friend class MatrixTest;
#endif //MYDEBUG
};

}/*math*/ 

#endif /*_MATH_RMAJ_MATRIX_h */
