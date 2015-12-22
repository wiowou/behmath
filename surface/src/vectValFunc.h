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

#ifndef _MATH_VECTVALFUNC_h
#define _MATH_VECTVALFUNC_h

#include "surface.h"
#include "matrix/matrix.h"
#include "matrix/subtr.h"
#include "matrix/elimination.h"
#include <cmath> //std::abs

namespace math{

#ifdef MYDEBUG
  class VectValFuncTest;
#endif //MYDEBUG

//! Represents a vector-valued function that takes a vector as an argument and returns
//! an output vector of the same size.
class VectValFunc
{

public:
  VectValFunc()
  {
    Init();
  }
  
  //! Explicitly sets the size of the output and Surface component vector
  explicit VectValFunc( ULong size )
  {
    Init( size );
  }
  
  void Clear()
  {
    Init();
  }
  
  //! The input vector
  void X( Vect<double>& x )
  {
    m_x = &x;
  }
  
  //! The output vector
  void B( Vect<double>& b )
  {
    m_b = &b;
  }
  
  //! Tolerance required for finding the x input for a given b output.
  void Tolerance( double tol )
  {
    m_tolerance = tol;
  }
  
  //! Maximum number of Newton iterations to try for finding the x input
  //! for a given b output.
  void MaxIter( ULong iter )
  {
    m_maxIter = iter;
  }

  ULong MaxIter()
  {
    return m_maxIter;
  }
  
  //! The number of interations required for finding the x input for a 
  //! given b output.
  ULong Iter()
  {
    return m_iter;
  }
  
  //! Set one of the components of the output vector to the function represented by Surface
  void Comp( ULong i, Surface* s )
  {
    m_comp[i] = s;
  }
  
  //! Sets the size of the output vector and size of the Surface component vector
  void Resize( ULong size )
  {
    m_comp.Resize(size);
  }
  
  //! Returns the size of the output and Surface component vector
  ULong Size()
  {
    return m_comp.Size();
  }
  
  virtual bool Continue()
  {
    return true;
  }

  //! x is the input vector, b is the output vector. The vector b must have been
  //! supplied for the method to provide output. Proivides for flexibility in that
  //! x and b do not have to be the same size.
  Vect<double>* operator()( Vect<double>& x )
  {
    Vect<double>& b = *m_b;
    
    //b.Resize( x.Size() );
    for ( ULong i = 0; i < b.Size(); ++i )
    {
      b[i] = (*m_comp[i])(x);
    }
    
    return m_b;
  }

  //! Find the input vector, x, for the given output vector b
  //! Uses Newton's method: {x1} = {x0} - {xmod}, {xmod} = [J]^-1( f{x0} - {b} )
  //! x and b must be of the same size for [J] to be square and invertible.
  void Solve()
  {
    ULong size = m_comp.Size();
    
    Vect<double> b = *m_b; //copying the original b
    Vect<double>& fx = *m_b;
    Vect<double>& x = *m_x;
    x.Resize( size );
    
    Vect<double> xmod( size ); //the modification to x
    Matrix<Vect> J( size );
    
    for ( ULong i = 0; i < size; ++i )
    {
      xmod[i] = m_tolerance * 2.0; //to ensure while loop is entered at least once
    }
    
    vops::Subtr<double> subtr;
    solver::Elimination<Vect> solver;
    solver.Tolerance(m_tolerance);
    solver.PartialPivoting(true);
    
    m_iter = 0;
    while ( m_iter < m_maxIter && Continue(xmod) )
    {
      operator()(x); //finding fx
      
      subtr( fx, b, xmod );
      
      xmod.Swap(fx);
      //create the Jacobian
      for ( ULong i = 0; i < size; ++i )
      {
        m_comp[i]->Grad( x, J[i] );
      }
      
      solver.A(J);
      solver.B(fx);
      solver.X(xmod);
      solver.Solve();
      
      subtr( x, xmod, x );
      
      ++m_iter;
    }
    
    m_b->Swap(b); //replacing b
  }  
protected:
  inline bool Continue( Vect<double>& xmod )
  {
    //this can be parallelized
    for ( ULong i = 0; i < xmod.Size(); ++i )
    {
      if ( std::abs( xmod[i] ) > m_tolerance )
      {
        return true;
      }
    }
    return false;
  }

  void Init( ULong size = 2 )
  {
    m_tolerance = 1e-5;
    m_maxIter = 0;
    m_iter = 0;
    m_comp.Resize(size);
    m_x = NULL;
    m_b = NULL;
  }
protected:
  double m_tolerance;
  ULong m_maxIter;
  ULong m_iter;
  
  //! various components of the vect val function. The i component is defined by the first surface: m_comp[0].
  //! The j component is defined by the second surface, m_comp[1] and so forth.
  Vect<Surface*> m_comp;
  
  //! The input vector when the vector b is the output vector
  Vect<double>* m_x;
  
  //! The output vector when the vector x is the input vector
  Vect<double>* m_b;
private:


#ifdef MYDEBUG
  friend class VectValFuncTest;
#endif //MYDEBUG
};

}/*math*/ 

#endif /*_MATH_VECTVALFUNC_h */
