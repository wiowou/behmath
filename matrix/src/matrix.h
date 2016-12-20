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

#include "vect.h"
#include "sparseVect.h"
#include <iostream>
#include <iomanip>


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
  
  explicit Matrix( unsigned long long size )
  {
    Init(size);
  }
  
  Matrix( unsigned long long mrow, unsigned long long ncol  )
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
    for ( unsigned long long i = 0; i < m_size; ++i )
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
  unsigned long long Rows()
  {
    if ( m_isTransposed )
    {
      return m_row[0].Size();
    }
    return m_size;
  }

  //! Returns number of physical rows of the matrix, not taking account of transposition
  unsigned long long PRows()
  {
    return m_size;
  }
  
  //! Returns number of cols of the matrix, taking account of transposition
  unsigned long long Cols()
  {
    if ( m_isTransposed )
    {
      return m_size;
    }
    return m_row[0].Size();
  }

  //! Returns number of cols of the matrix, not taking account of transposition
  unsigned long long PCols()
  {
    return m_row[0].Size();
  }
  
  void Swap( Matrix &other )
  {
    bool tmp_isTransposed = other.m_isTransposed;
    unsigned long long tmp_size = other.m_size;
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
  void Resize( unsigned long long size )
  {
    if ( size <= m_size )
    {
      m_size = size;
      for ( unsigned long long i = 0; i < m_size; ++i )
      {
        m_row[i].Resize(m_size);
      }
      return;
    }    
    Clear();
    Init(size);
  }
  
  //! Re-dimensions the matrix and its child row vectors. Deletes previous data contained within.
  void Resize( unsigned long long mrow, unsigned long long ncol )
  {
    m_isTransposed = false;
    if ( mrow <= m_size )
    {
      m_size = mrow;
      for ( unsigned long long i = 0; i < m_size; ++i )
      {
        m_row[i].Resize(ncol);
      }
      return;
    }

    Storage<T>* row = new Storage<T> [mrow];  
    for ( unsigned long long i = 0; i < m_size; ++i )
    {
      row[i] = m_row[i];
    }
    m_size = mrow;
    Destroy();
    m_row = row;
    
    for ( unsigned long long i = 0; i < m_size; ++i )
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
    unsigned long long nrow = m_size;
    unsigned long long ncol = m_row[0].Size();
    for ( unsigned long long i = 0; i < nrow; ++i )
    {
      for ( unsigned long long j = 0; j < ncol; ++j )
      {
        Set( i,j, *data++ );
      }
    }
  }
  
  //! Returns a row vector
  Storage<T>& operator[]( unsigned long long m )
  {
    return m_row[m];
  }
  
  //! Gets/Sets the m,n element of the matrix. Recognizes if matrix is transposed.
  T& operator()( unsigned long long m, unsigned long long n ) //tested
  {
    if ( m_isTransposed )
    {
      return m_row[n][m];
    }
    return m_row[m][n];
  }
  
  //! Sets the m,n element of the matrix. Recognizes if matrix is transposed.
  void Set( unsigned long long m, unsigned long long n, T val ) //tested
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
  
  void RowSwap( unsigned long long rowNum1, unsigned long long rowNum2 ) //tested
  {
    m_row[rowNum1].Swap( m_row[rowNum2] );
  }

  //! Returns the row number for the maximum pivot value below @param pivotIdx row
  unsigned long long MaxIdxBelow( const unsigned long long pivotIdx ) const //tested
  {
    unsigned long long rowEnd = m_size;
    unsigned long long max = pivotIdx;
    T maxVal = m_row[pivotIdx][pivotIdx];
    for ( unsigned long long i = pivotIdx + 1; i < rowEnd; ++i )
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
  unsigned long long MaxIdxAbove( const unsigned long long pivotIdx ) const //tested
  {
    unsigned long long rowEnd = 0;
    unsigned long long max = pivotIdx;
    T maxVal = m_row[pivotIdx][pivotIdx];
    for ( unsigned long long i = pivotIdx; i > rowEnd; --i )
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
    unsigned long long nrow = A.m_size;
    unsigned long long ncol = A[0].Size();
    if ( A.IsTransposed() )
    {
      nrow = A[0].Size();
      ncol = A.m_size;
    }
    std::ios::fmtflags fla(os.flags() );
    os << "Matrix( ";
    os << nrow << "x" << ncol << " )" << " = " << std::endl;
    
    for (unsigned long long i = 0; i < nrow; ++i )
    {
      os << "  [";
      for (unsigned long long j = 0; j < ncol; ++j )
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
  
  void Init( unsigned long long size )
  {
    m_isTransposed = false;
    size = size < 2 ? 2 : size;
    m_row = new Storage<T> [size];
    m_size = size;
    
    for ( unsigned long long i = 0; i < m_size; ++i )
    {
      m_row[i].Resize(m_size);
    }
  }
  
  void Init( unsigned long long mrow, unsigned long long ncol )
  { 
    m_isTransposed = false;
    mrow = mrow < 2 ? 2 : mrow;
    m_row = new Storage<T> [mrow];
    m_size = mrow;
    
    for ( unsigned long long i = 0; i < m_size; ++i )
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
  unsigned long long m_size;
  
  //! pointer to a row vector
  Storage<T>* m_row;

private:

#ifdef MYDEBUG
  friend class MatrixTest;
#endif //MYDEBUG
};

}/*math*/ 

#endif /*_MATH_RMAJ_MATRIX_h */
