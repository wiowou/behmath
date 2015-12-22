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

#ifndef _MATH_SPARSEVECT_h
#define _MATH_SPARSEVECT_h

#include "impl/matrixConfig.h"
#include <iostream>
#include <iomanip>
#include "typedefs.h"

namespace math{

#ifdef MYDEBUG
  class SparseVectTest;
#endif //MYDEBUG

namespace vops
{
  template<class T> class Sum;
  template<class T> class Add;
  template<class T> class Subtr;
  template<class T> class Mult;
  template<class T> class Axpy;
  template<class T> class Dot;
}

template<class T>
class SparseVect
{
public:
  SparseVect() //tested
  {
    Init(2);
  }
  
  explicit SparseVect( ULong allocSize ) //tested
  {
    Init( allocSize );
  }
  
  ~SparseVect() //tested
  {
    Destroy();
  }

  //! copy constructor
  SparseVect( const SparseVect &other )
  {
    m_nnz = other.m_nnz;
    m_size = other.m_size;
    m_allocSize = other.m_allocSize;
    m_pos = new ULong [m_allocSize];
    m_data = new T [m_allocSize];
    for ( ULong i = 0; i < m_nnz; ++i )
    {
      m_pos[i] = other.m_pos[i];
      m_data[i] = other.m_data[i];
    }
  }
  
  //! assignment operator
  SparseVect& operator=( SparseVect other )
  {
    Swap(other);
    return *this;
  }
  
  void Swap( SparseVect &other )
  {
    ULong tmp_nnz = other.m_nnz;
    ULong tmp_size = other.m_size;
    ULong tmp_allocSize = other.m_allocSize;
    ULong* tmp_pos = other.m_pos;
    T* tmp_data = other.m_data;
    
    other.m_nnz = m_nnz;
    other.m_size = m_size;
    other.m_allocSize = m_allocSize;
    other.m_pos = m_pos;
    other.m_data = m_data;
    
    m_nnz = tmp_nnz;
    m_size = tmp_size;
    m_allocSize = tmp_allocSize;
    m_pos = tmp_pos;
    m_data = tmp_data;
  }
  
  void Clear()
  {
    Destroy();
    Init(2);
  }
  
  //! Input the data. If @param pos is not specified, the data is read until the size of the
  //! vector is reached. If @param pos is specified, it should be an array of size @param nnz.
  //! The @param pos should indicate the index, or position, of the corresponding entry in @param data.
  //! The @param nnz indicates the size of both data and pos. Vector is resized to allow @param nnz entries.
  void Data( T* data, ULong* pos = NULL, ULong nnz = 0 )
  {
    if ( pos == NULL || nnz == 0 )
    {
      return;
    }
    Clear();
    Allocate(nnz);
    for ( ULong i = 0; i < nnz; ++i )
    {
      PushBack( pos[i], data[i] );
    }
  }
  
  T& Front()
  {
    return m_data[0];
  }
  
  T& Back()
  {
    return m_data[m_nnz - 1];
  }
  
  //! Returns the location in m_data where the max occurs.
  ULong MaxIdx() //tested
  {
    ULong max = 0;
    for ( ULong i = 1; i < m_nnz; ++i )
    {
      if ( m_data[i] * m_data[i] > m_data[max] * m_data[max] )
      {
        max = i;
      }
    }
    return max;
  }
  
  //! Returns the index located in m_pos at location @param i
  ULong Pos( ULong i ) const
  {
    return m_pos[i];
  }
  
  //! Returns the index located in m_data at location @param i.
  //! Use this in conjunction with MaxIdx()
  T& operator()( ULong i )
  {
    return m_data[i];
  }
  
  //! Returns the number of entries not logically considered to be zero.
  //! ie, returns the size of m_data
  ULong NNZ() const
  {
    return m_nnz;
  }
  
  //! Returns the dimensionality that the vector is considered to be
  ULong Size() const
  {
    return m_size;
  }
  
  //! Changes the size of the number of entries logically considered to be zero
  void Resize( ULong size )
  {
    m_size = size;
    if ( m_nnz > 0 && size < m_pos[ m_nnz - 1] )
    {
      m_nnz = UpperBoundIdx(size);
    }
  }
  
  //! Allocates memory in m_data of @param size. Copies old data to new vector
  //! if necessary
  void Allocate( ULong size ) //tested
  {
    if ( size <= m_allocSize )
    {
      return;
    }
    m_size = m_allocSize = size;
    ULong* pos = new ULong [m_allocSize];
    T* data = new T [m_allocSize];
    for ( ULong i = 0; i < m_nnz; ++i)
    {
      pos[i] = m_pos[i];
      data[i] = m_data[i];
    }
    ULong* tmp = m_pos;
    T* tmp2 = m_data;
    m_pos = pos;
    m_data = data;
    delete[] tmp;
    delete[] tmp2;
  }
  
  //! Returns the value at index idx. If idx is logically considered to be zero,
  //! returns default constructor for template type T
  T& operator[]( const ULong idx ) //tested
  {
    ULong ub = UpperBoundIdx(idx);
    if ( ub == m_nnz || m_pos[ub] != idx )
    {
      m_default = T();
      return m_default;
    }
    return m_data[ub];
  }
  
  //! Sets and overwrites any previous value at idx. Idx is the logical index to the vector, not the physical one.
  //! Will resize the allocated data if necessary.
  void Set( ULong idx, T val ) //tested
  {
    ULong ub = UpperBoundIdx(idx);
    Set( idx, val, ub );
  }

  /*
  void SetData( ULong i, T val )
  {
    m_data[i] = val;
  }
  */
  
  //! Retrieves data value at idx. Idx is the logical index to the vector, not the physical one.
  T Get( const ULong idx ) const
  {
    ULong ub = UpperBoundIdx(idx);
    if ( ub == m_nnz || m_pos[ub] != idx )
    {
      m_default = T();
      return m_default;
    }
    return m_data[ub];
  }
  
  //! Adds a logical non-zero entry to the end of the vector
  void PushBack( ULong idx, T val ) //tested
  {
    if ( m_nnz == m_allocSize - 1 )
    {
      Allocate( 2 * m_allocSize );
    }
    m_data[m_nnz] = val;
    m_pos[m_nnz] = idx;
    ++m_nnz;
    if ( m_size < m_nnz )
    {
      m_size = m_nnz;
    }
  }
  
  //! Returns the physical index associated with logical index @param idx. If @param idx
  //! is not found, returns m_nnz;
  ULong Find( const ULong idx ) const //tested
  {
    ULong ub = UpperBoundIdx(idx);
    if ( ub == m_nnz || m_pos[ub] != idx )
    {
      return m_nnz;
    }
    return ub;
  }
  
  //! Performs a row reduction, assuming vector is a column vector.
  //! vector[reducedRow] += factor * vector[reducingRow]. Does the right thing if either
  //! or both indexes are not found
  inline void RowReduce(  T factor, ULong reducingRow, ULong reducedRow ) //tested
  {
    ULong reducingRowIdx = Find( reducingRow );
    if ( reducingRowIdx == m_nnz )
    {
      return;
    }
    ULong reducedRowIdx = Find( reducedRow );
    if ( reducedRowIdx == m_nnz )
    {
      Set( reducedRow, factor * m_data[reducingRowIdx] );
      return;
    }
    m_data[reducedRowIdx] += factor * m_data[reducingRowIdx];
  }
  
  //! Swaps entries, assuming vector is a column vector
  inline void RowSwap( ULong rowNum1, ULong rowNum2 )
  {
    if ( rowNum1 == rowNum2 )
    {
      return;
    }
    ULong rowNum1idx = Find(rowNum1);
    ULong rowNum2idx = Find(rowNum2);
    
    if ( rowNum1idx == m_nnz && rowNum2idx == m_nnz )
    {
      return;
    }
    
    if ( rowNum1idx == m_nnz )
    {
      Insert( rowNum1, m_data[rowNum2idx] );
      m_data[rowNum2idx] = T();
      return;
    }
    
    if ( rowNum2idx == m_nnz )
    {
      Insert( rowNum2, m_data[rowNum1idx] );
      m_data[rowNum1idx] = T();
      return;
    }
    
    T tmp = m_data[rowNum1idx];
    m_data[rowNum1idx] = m_data[rowNum2idx];
    m_data[rowNum2idx] = tmp;
  }
  
  //! Print the vector
  friend std::ostream& operator<<( std::ostream& os, SparseVect& A )
  {
    std::ios::fmtflags fla(os.flags() );
    os << "Vector( ";
    os << A.Size() << " )" << " = " << std::endl;
    for ( ULong i = 0; i < A.Size(); ++i )
    {
      os << "  [";
      os << std::fixed << std::setw(10) << std::setprecision(4);
      os << A[i];
      os.flags(fla); 
      os << " ]" << std::endl;       
    }
    return os;
  }
  
protected:
  void Init( ULong allocSize ) //tested
  {
    m_nnz = 0;
    m_size = allocSize;
    m_allocSize = allocSize;
    m_data = new T [m_allocSize];
    m_pos = new ULong [m_allocSize];
  }
  
  void Destroy() //tested
  {
    delete[] m_data;
    delete[] m_pos;
  }
  
  //! Returns the upper bound physical index in m_data of @param val, which is a logical index.
  //! If val > the end, returns m_nnz. If val == an logical index, returns the corresponding
  //! physical index
  ULong UpperBoundIdx( const ULong val ) const //tested
  {
    ULong const* beg = m_pos;
    ULong const* end = &m_pos[m_nnz];
    if ( val <= *beg )
    {
      return 0;
    }
    --end;
    if ( val > *end )
    {
      return ( end - beg ) + 1;
    }
    ULong const* zeroPos = beg;
    while ( ( end - beg ) > 2 )
    {
      //ULong const* mid = &beg[ ( end - beg ) >> 1 ]; // >> 1 is division by 2
      double pos = ( static_cast<double>(val) - static_cast<double>(*beg) ) / 
        ( static_cast<double>(*end) - static_cast<double>(*beg) ) * static_cast<double>(end - beg);
      double fraction = pos - static_cast<ULong>(pos);
      if ( fraction < 0.001 )
      {
        pos += 0.01;
      }
      ULong const* mid = &beg[ static_cast<ULong>(pos) ];
      if ( *mid == val )
      {
        return (mid - zeroPos);
      }
      if ( mid == beg )
      {
        ++mid;
      }
      else if ( mid == end )
      {
        --mid;
      }
      if ( *mid > val )
      {
        end = mid;
      }
      else
      {
        beg = mid;
      }
    }
    ULong dist = end - beg;
    for ( ULong i = 0; i <= dist; ++i )
    {
      if ( *beg >= val )
      {
        return (beg - zeroPos);
      }
      ++beg;    
    }
  }

  //! Sets and overwrites any previous value at idx. Idx is the logical index to the vector, not the physical one.
  //! Will resize the allocated data if necessary.  
  void Set( ULong idx, T val, ULong ub ) //tested
  {
    if ( m_nnz == 0 )
    {
      m_pos[0] = idx;
      m_data[0] = val;
      m_nnz++;
      if ( m_size < m_nnz )
      {
        m_size = m_nnz;
      }
      return;
    }
    if ( ub < m_nnz && m_pos[ub] == idx )
    {
      m_data[ub] = val;
      return;
    }
    if ( m_nnz == m_allocSize - 1 )
    {
      m_allocSize = 2 * m_allocSize;
      ULong* pos = new ULong [m_allocSize];
      T* data = new T [m_allocSize];
      for ( ULong i = 0; i < ub; ++i)
      {
        pos[i] = m_pos[i];
        data[i] = m_data[i];
      }
      pos[ub] = idx;
      data[ub] = val;
      for ( ULong i = ub + 1; i < m_nnz + 1; ++i )
      {
        pos[i] = m_pos[i - 1];
        data[i] = m_data[i - 1];      
      }
      ULong* tmp = m_pos;
      T* tmp2 = m_data;
      m_pos = pos;
      m_data = data;
      delete[] tmp;
      delete[] tmp2;
      
      m_nnz++;
      if ( m_size < m_nnz )
      {
        m_size = m_nnz;
      }
      return;
    }
    for ( ULong i = m_nnz; i > ub; --i )
    {
      m_pos[i] = m_pos[i - 1];
      m_data[i] = m_data[i - 1];
    }
    m_pos[ub] = idx;
    m_data[ub] = val;
    
    m_nnz++;
    if ( m_size < m_nnz )
    {
      m_size = m_nnz;
    }
    return;
  }

protected: //members
  //! number of non-zero entries
  ULong m_nnz;
  
  //! allocation size
  ULong m_allocSize;
  
  //! logical size of the vector
  ULong m_size;
  
  //! the logical indexes corresponding to the physical indexes of m_data
  ULong* m_pos;
  
  //! the array of data 
  T* m_data;
  
  //! A default T value to return when the logical index isn't found
  T m_default;

  friend class vops::Sum<T>;
  friend class vops::Add<T>;
  friend class vops::Subtr<T>;
  friend class vops::Mult<T>;
  friend class vops::Axpy<T>;
  friend class vops::Dot<T>;
  
  
private:


#ifdef MYDEBUG
  friend class SparseVectTest;
#endif //MYDEBUG
};

}/*math*/ 

#endif /*_MATH_SPARSEVECT_h */
