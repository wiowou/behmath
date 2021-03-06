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

#ifndef _MATH_SOLVER_CYCLIC_h
#define _MATH_SOLVER_CYCLIC_h

#include "iterative.h"
#include "cyclicRow.h"

namespace math{
namespace solver{

#ifdef MYDEBUG
  class CyclicTest;
#endif //MYDEBUG

template< template <typename> class Storage = SparseVect, typename T = double >
class Cyclic : public Iterative<Storage,T>
{
  using Iterative<Storage,T>::m_A;
  using Iterative<Storage,T>::m_x;
  using Iterative<Storage,T>::m_b;
  using Iterative<Storage,T>::m_notWithinTol;
  using Iterative<Storage,T>::m_xNew;
  using Iterative<Storage,T>::m_tolerance;
  using Iterative<Storage,T>::m_maxIter;
  using Iterative<Storage,T>::m_nIter;
  using Iterative<Storage,T>::Continue;
  using Iterative<Storage,T>::InitialGuess;
  
public:
  void Solve()
  {
    InitialGuess();
    
    Vect<CyclicRow<Storage,T> > rowSolver( m_x->Size() );
    for ( unsigned long long i = 0; i < rowSolver.Size(); ++i )
    {
      rowSolver[i].SetRow( &(*m_A)[i] );
      rowSolver[i].SetX( m_x );
      rowSolver[i].SetXNew( &m_xNew );
      rowSolver[i].SetB( m_b );
      rowSolver[i].SetNotWithinTol( &m_notWithinTol );
      rowSolver[i].SetIdx(i);
      rowSolver[i].SetTolerance( m_tolerance );
      rowSolver[i].SetCoeff();
    }
    
    m_nIter = 0;
    while ( m_nIter < m_maxIter && Continue() )
    {
      //this loop cannot be parallelized
      for ( unsigned long long i = 0; i < rowSolver.Size(); ++i )
      {
        rowSolver[i].Exec();
      }
      m_xNew.Swap(*m_x);
      ++m_nIter;
    }
    return;
  }
protected:

private:


#ifdef MYDEBUG
  friend class CyclicTest;
#endif //MYDEBUG
};

}/*solver*/ }/*math*/ 

#endif /*_MATH_SOLVER_CYCLIC_h */
