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

#ifndef _MATH_SOLVER_ROWSOLVER_h
#define _MATH_SOLVER_ROWSOLVER_h

#include "typedefs.h"
#include "impl/matrixConfig.h"
#include "dot.h"

namespace math{
namespace solver{

#ifdef MYDEBUG
  class RowSolverTest;
#endif //MYDEBUG

template< template <typename> class Storage = SparseVect, typename T = double >
class RowSolver
{
public:
  RowSolver()
  {
    m_row = NULL;
    m_x = NULL;
    m_xNew = NULL;
    m_b = NULL;
    m_notWithinTol = NULL;
    m_idx = 0;
    m_coeff = 1.0;
    m_tolerance = -1.0; //this means absDelta for this idx is always     
  }
  
  void SetRow( Storage<T>* row )
  {
    m_row = row;
  }

  void SetX( Vect<T>* x )
  {
    m_x = x;
  }

  void SetXNew( Vect<T>* xNew )
  {
    m_xNew = xNew;
  }


  void SetB( Vect<T>* b )
  {
    m_b = b;
  }

  void SetNotWithinTol( Vect<bool>* notWithinTol )
  {
    m_notWithinTol = notWithinTol;
  }

  void SetIdx( ULong idx )
  {
    m_idx = idx;
  }

  void SetCoeff()
  {
    m_coeff = m_row->Get(m_idx);
  }

  void SetTolerance( T tolerance )
  {
    m_tolerance = tolerance;
  }

protected:
  Storage<T>* m_row;
  Vect<T>* m_x;
  Vect<T>* m_xNew;
  Vect<T>* m_b;
  Vect<bool>* m_notWithinTol;
  ULong m_idx;
  T m_coeff;
  T m_tolerance;
  vops::Dot<T> dot;
private:


#ifdef MYDEBUG
  friend class RowSolverTest;
#endif //MYDEBUG
};

}/*solver*/ }/*math*/ 

#endif /*_MATH_SOLVER_ROWSOLVER_h */
