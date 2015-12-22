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

#ifndef _MATH_SOLVER_RECTANGULAR_h
#define _MATH_SOLVER_RECTANGULAR_h

#include "matrix.h"
#include "dot.h"
#include "mult.h"
#include "subtr.h"
#include <math.h> //for sqrt
//#include "exception.h"

namespace math{
namespace solver{

#ifdef MYDEBUG
  class RectangularTest;
#endif //MYDEBUG

//! Solves rectangular m x n systems where m > n. The solver is a template parameter.
//! Solves the system [A^T][C]( [b] - [A][x] ) = [f] -> [A^T][C][A][x] = [A^T][C][b] - [f]
template< template< template <typename> class Storage, typename T > class Sol, template <typename> class Storage = SparseVect, typename T = double >
class Rectangular
{

public:
  Rectangular()
  {
    Init();
  }
  
  void Clear()
  {
    Init();
  }
  
  //! This is the transpose of the incidence matrix. n x m
  void ATranspose( Matrix< Storage, T >& AT )
  {
    m_AT = &AT;
    m_updateATCA = true;
    m_updataATCBminusF = true;
  }
  
  //! The batteries or specified potentials. m x 1
  void B( Storage<T>& b )
  {
    m_b = &b;
    m_updataATCBminusF = true;
  }
  
  //! The drop in potential across each edge. n x 1
  void X( Vect<T>& x )
  {
    m_x = &x;
  }

  //! The weights for each edge. m x m diagonal matrix, but specified as an m x 1 vector
  //! for the purposes of efficiency.
  void C( Vect<T>& c )
  {
    m_c = &c;
    m_updateATCA = true;
    m_updataATCBminusF = true;
  }

  //! The specified flows along each edge. n x 1
  void F( Storage<T>& f )
  {
    m_f = &f;
    m_updataATCBminusF = true;
  }
  
  bool DimProblem()
  {
    Matrix< Storage, T >& AT = *m_AT;
    Storage<T>& b = *m_b;
    Vect<T>& x = *m_x;
    Vect<T>& c = *m_c;
    Storage<T>& f = *m_f;
    
    if ( AT.Cols() < AT.Rows() )
    {
      return true;
    }
    if ( DimCheck( b, AT.Cols() ) == true  )
    {
      return true;
    }
    if ( DimCheck( f, AT.Rows() ) == true  )
    {
      return true;
    }
    if ( DimCheck( x, AT.Rows() ) == true  )
    {
      return true;
    }
    if ( DimCheck( c, AT.Cols() ) == true  )
    {
      return true;
    }
    
    return false;
  }
  
  void Solve()
  {
    if ( m_updataATCBminusF == true )
    {
      UpdateATCBminusF();
    }

    if ( m_updateATCA == true )
    {
      UpdateATCA();
    }
    std::cout << m_ATCA << std::endl;
    std::cout << m_ATCBminusF << std::endl;
    Sol<Storage,T> solver;
    solver.A( m_ATCA );
    solver.B( m_ATCBminusF );
    solver.X( *m_x );
    solver.Solve();
  }

protected:
  void UpdateATCA()
  {
    Matrix< Storage, T >& AT = *m_AT;
    Vect<T>& c = *m_c;

    m_ATCA.Resize( AT.Rows() );
    vops::Dot<T> dot;
    
    if ( m_c != NULL )
    {
      /*
      for ( ULong i = 0; i < c.NNZ(); ++i )
      {
        c(i) = sqrt( c(i) ); //weighting factors must be positive
      }
      vops::Mult<T> mult;
      for ( ULong i = 0; i < AT.Rows(); ++i )
      {
        mult( AT[i], c );
      }
      */
      dot( AT, c, m_ATCA ); //AT dotted with itself with c in between, output to m_ATCA
      m_updateATCA = false;
      return;
    }
    
    dot( AT, m_ATCA ); //AT dotted with itself, output to m_ATCA
    m_updateATCA = false;
  }
  
  void UpdateATCBminusF()
  {
    Matrix< Storage, T >& AT = *m_AT;
    Vect<T>& c = *m_c;
    Storage<T>& b = *m_b;
    Storage<T>& f = *m_f;
    
    if ( m_c != NULL )
    {
      vops::Mult<T> mult;
      mult( b, c );
    }
    
    m_ATCBminusF.Resize( AT.Rows() );
    Storage<T> tmp( AT.Rows() );
    vops::Dot<T> dot;
    dot( AT, b, tmp );
    
    m_updataATCBminusF = false;
    if ( m_f == NULL )
    {
      m_ATCBminusF.Swap(tmp);
      return;
    }
    
    vops::Subtr<T> subtr;
    subtr( tmp, f, m_ATCBminusF );
  }
  
  bool DimCheck( SparseVect<T> &v, ULong size )
  {
    return v.NNZ() > size;
  }
  
  bool DimCheck( Vect<T> &v, ULong size )
  {
    return v.Size() != size;
  }
  
  void Init()
  {
    m_AT = NULL;
    m_b = NULL;
    m_x = NULL;
    m_c = NULL;
    m_f = NULL;
    m_updateATCA = true;
    m_updataATCBminusF = true;
  }
protected:
  Matrix< Storage, T >* m_AT;
  Storage<T>* m_b;
  Vect<T>* m_x;
  Vect<T>* m_c;
  Storage<T>* m_f;
  
  //! [A^T][C][A], n x n
  Matrix< Storage, T > m_ATCA;
  //! [A^T][C][b] - [f], n x 1
  Storage<T> m_ATCBminusF;
  bool m_updateATCA;
  bool m_updataATCBminusF;
private:


#ifdef MYDEBUG
  friend class RectangularTest;
#endif //MYDEBUG
};

}/*solver*/ }/*math*/ 

#endif /*_MATH_SOLVER_RECTANGULAR_h */
