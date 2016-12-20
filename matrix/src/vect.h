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

#ifndef _MATH_VECT_h
#define _MATH_VECT_h

#include <iostream>
#include <iomanip>


namespace math{

#ifdef MYDEBUG
  class VectTest;
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
class Vect
{
public:
  Vect()
  {
    Init(2);    
  }
  
  explicit Vect( unsigned long long size )
  {
    Init(size);
  }
  
  ~Vect()
  {
    Destroy();
  }

  //! copy constructor
  Vect( const Vect &other )
  {
    m_size = other.m_size;
    m_allocSize = other.m_allocSize;
    m_data = new T [m_allocSize];
    for ( unsigned long long i = 0; i < m_size; ++i )
    {
      m_data[i] = other.m_data[i];
    }
  }
  
  //! assignment operator
  Vect& operator=( Vect other )
  {
    Swap(other);
    return *this;
  }
  
  void Swap( Vect &other )
  {
    unsigned long long tmp_size = other.m_size;
    unsigned long long tmp_allocSize = other.m_allocSize;
    T* tmp_data = other.m_data;
    
    other.m_size = m_size;
    other.m_allocSize = m_allocSize;
    other.m_data = m_data;
    
    m_size = tmp_size;
    m_allocSize = tmp_allocSize;
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
  void Data( T* data, unsigned long long* pos = NULL, unsigned long long nnz = 0 )
  {
    if ( pos != NULL )
    {
      if ( nnz == 0 )
      {
        return;
      }
      Resize( pos[nnz - 1] );
      for ( unsigned long long i = 0, j = 0; i < m_size; ++i )
      {
        if ( i == pos[j] )
        {
          Set( i, *data++ );
          ++j;
        }
        else
        {
          Set( i, T() );
        }
      }
      return;
    }
    
    for ( unsigned long long i = 0; i < m_size; ++i )
    {
      Set( i, *data++ );
    }
  }
  
  T& Front()
  {
    return m_data[0];
  }
  
  T& Back()
  {
    return m_data[m_size - 1];
  }
  
  //! Returns the location in m_data where the max occurs.
  unsigned long long MaxIdx( unsigned long long start = 0 )
  {
    unsigned long long max = start;
    for ( unsigned long long i = start + 1; i < m_size; ++i )
    {
      if ( m_data[i] * m_data[i] > m_data[max] * m_data[max] )
      {
        max = i;
      }
    }
    return max;
  }

  //! Returns i. For compatibility with SparseVect
  unsigned long long Pos( unsigned long long i ) const
  {
    return i;
  }
  
  //! Returns the index located in m_data at location @param i.
  T& operator()( unsigned long long i )
  {
    return m_data[i];
  }
  
  //! Returns the number of entries not logically considered to be zero.
  //! ie, returns the size of m_data  
  unsigned long long NNZ() const
  {
    return m_size;
  }
  
  //! Returns the dimensionality that the vector is considered to be
  unsigned long long Size() const
  {
    return m_size;
  }
  
  //! Changes the size of the number of entries logically considered to be zero
  void Resize( unsigned long long size )
  {
    if ( size <= m_allocSize )
    {
      m_size = size;
      return;
    }
    T* tmp = m_data;
    m_data = new T [size];
    for ( unsigned long long i = 0; i < m_size; ++i )
    {
      m_data[i] = tmp[i];
    }
    delete[] tmp;
    m_size = m_allocSize = size;
  }

  //! Allocates memory in m_data of @param size. Copies old data to new vector
  //! if necessary  
  void Allocate( unsigned long long size )
  {
    Resize(size);
  }

  //! Returns the value at index idx. If idx is logically considered to be zero,
  //! returns default constructor for template type T  
  T& operator[]( const unsigned long long idx )
  {
    return m_data[ idx ];
  }
  
  //! Sets and overwrites any previous value at idx. Idx is the logical index to the vector, not the physical one.
  //! Will resize the allocated data if necessary.
  void Set( unsigned long long idx, T val )
  {
    if ( idx >= m_size )
    {
      Resize( idx + 1);
    }
    m_data[ idx ] = val;
    return;
  }
  
  //! Retrieves data value at idx. Idx is the logical index to the vector, not the physical one.
  T Get( const unsigned long long idx ) const
  {
    return m_data[idx];
  }
  
  
  //! Adds a logical non-zero entry to the end of the vector
  void PushBack( unsigned long long idx, T val )
  {
    Set( idx, val );
  }

  //! Performs a row reduction, assuming vector is a column vector.
  //! vector[reducedRow] += factor * vector[reducingRow]. Does the right thing if either
  //! or both indexes are not found  
  inline void RowReduce(  T factor, unsigned long long reducingRow, unsigned long long reducedRow )
  {
    m_data[reducedRow] += factor * m_data[reducingRow];
  }
  
  //! Swaps entries, assuming vector is a column vector
  inline void RowSwap( unsigned long long rowNum1, unsigned long long rowNum2 )
  {
    if ( rowNum1 == rowNum2 )
    {
      return;
    }
    T tmp = m_data[rowNum1];
    m_data[rowNum1] = m_data[rowNum2];
    m_data[rowNum2] = tmp;
  }

  //! Print the vector
  friend std::ostream& operator<<( std::ostream& os, Vect& A )
  {
    std::ios::fmtflags fla(os.flags() );
    os << "Vector( ";
    os << A.Size() << " )" << " = " << std::endl;
    for ( unsigned long long i = 0; i < A.Size(); ++i )
    {
      os << "  [";
      os << std::fixed << std::setw(10) << std::setprecision(4);
      os << A[i];
      os.flags(fla); 
      os << " ]" << std::endl;       
    }
    return os;
  }
  
protected: //methods
  void Init( unsigned long long size )
  {
    m_size = size;
    m_allocSize = size;
    m_data = new T [m_allocSize];
  }
  
  void Destroy()
  {
    delete[] m_data;
  }
  
protected: //members
  //! the logical and physical size of vector
  unsigned long long m_size;
  
  //! The memory allocation size
  unsigned long long m_allocSize;
  
  //! the array of data
  T* m_data;

  friend class vops::Sum<T>;
  friend class vops::Add<T>;
  friend class vops::Subtr<T>;
  friend class vops::Mult<T>;
  friend class vops::Axpy<T>;
  friend class vops::Dot<T>;
  
private:


#ifdef MYDEBUG
  friend class VectTest;
#endif //MYDEBUG
};

}/*math*/ 

#endif /*_MATH_VECT_h */
