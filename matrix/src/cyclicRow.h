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

#ifndef _MATH_SOLVER_CYCLICROW_h
#define _MATH_SOLVER_CYCLICROW_h

#include "rowSolver.h"
#include <cmath> //std::abs

namespace math{
namespace solver{

#ifdef MYDEBUG
  class CyclicRowTest;
#endif //MYDEBUG

template< template <typename> class Storage = SparseVect, typename T = double >
class CyclicRow : public RowSolver<Storage, T>
{
  using RowSolver<Storage,T>::m_row;
  using RowSolver<Storage,T>::m_x;
  using RowSolver<Storage,T>::m_xNew;
  using RowSolver<Storage,T>::m_b;
  using RowSolver<Storage,T>::m_notWithinTol;
  using RowSolver<Storage,T>::m_idx;
  using RowSolver<Storage,T>::m_coeff;
  using RowSolver<Storage,T>::dot;
  using RowSolver<Storage,T>::m_tolerance;
public:
  void Exec()
  {
    T res = dot( m_idx, *m_xNew, *m_row ) + dot( *m_x, *m_row, m_idx + 1 );
    (*m_xNew)[m_idx] = ( (*m_b)[m_idx] - res ) / m_coeff;
    (*m_notWithinTol)[m_idx] = std::abs( (*m_xNew)[m_idx] - (*m_x)[m_idx] ) > m_tolerance;
  }
protected:

private:


#ifdef MYDEBUG
  friend class CyclicRowTest;
#endif //MYDEBUG
};

}/*solver*/ }/*math*/ 

#endif /*_MATH_SOLVER_CYCLICROW_h */
