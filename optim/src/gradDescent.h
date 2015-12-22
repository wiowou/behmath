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

#ifndef _MATH_OPTIM_GRADDESCENT_h
#define _MATH_OPTIM_GRADDESCENT_h

#include "minFinder.h"
#include "matrix/vect.h"
#include "matrix/mult.h"
#include "matrix/subtr.h"
#include <cmath> //std::abs

namespace math{
namespace optim{

#ifdef MYDEBUG
  class GradDescentTest;
#endif //MYDEBUG

template<class F, class T>
class GradDescent : public MinFinder<F,T>
{

public:
  
protected:

private:


#ifdef MYDEBUG
  friend class GradDescentTest;
#endif //MYDEBUG
};

template<class F>
class GradDescent<F,double> : public MinFinder<F,double>
{

  using MinFinder<F,double>::m_f;
  using MinFinder<F,double>::m_tolerance;
  using MinFinder<F,double>::m_min;
  using MinFinder<F,double>::m_maxIter;
  using MinFinder<F,double>::m_iter;
  using MinFinder<F,double>::m_runCalc;
  using MinFinder<F,double>::m_minFound;
  using MinFinder<F,double>::m_low;
  using MinFinder<F,double>::m_high;
  using MinFinder<F,double>::m_gain;
  using MinFinder<F,double>::m_guess;
  
public:
  GradDescent()
  {
    Init();
  }
  
  void Clear()
  {
    Init();
    MinFinder<F,double>::Clear();
  }
  
  double MinX()
  {
    if ( m_runCalc )
    {
      Calc();
    } 
    return m_min;
  }
  
protected:
  void Calc()
  {
    double x = m_guess;
    double fprime = m_f->Deriv1( x );
    
    m_iter = 0;
    while ( std::abs( fprime ) > m_tolerance && m_iter < m_maxIter )
    {
      x -= m_gain * fprime;
      fprime = m_f->Deriv1( x );
      ++m_iter;
    }
    m_min = x;
    
    m_runCalc = false;
    if ( m_iter < m_maxIter )
    {
      m_minFound = true;
    }
  }
  
  void Init()
  {
    m_guess = 1.0;
    m_gain = 1.0;    
  }
private:


#ifdef MYDEBUG
  friend class GradDescentTest;
#endif //MYDEBUG
};

template<class F>
class GradDescent<F,Vect<double> > : public MinFinder<F,Vect<double> >
{

  using MinFinder<F,Vect<double> >::m_f;
  using MinFinder<F,Vect<double> >::m_tolerance;
  using MinFinder<F,Vect<double> >::m_min;
  using MinFinder<F,Vect<double> >::m_maxIter;
  using MinFinder<F,Vect<double> >::m_iter;
  using MinFinder<F,Vect<double> >::m_runCalc;
  using MinFinder<F,Vect<double> >::m_minFound;
  using MinFinder<F,Vect<double> >::m_low;
  using MinFinder<F,Vect<double> >::m_high;
  using MinFinder<F,Vect<double> >::m_gain;
  using MinFinder<F,Vect<double> >::m_guess;
  
public:
  
  Vect<double>* MinX()
  {
    if ( m_runCalc )
    {
      Calc();
    } 
    return &m_min;
  }
  
protected:
  void Calc()
  {
    Vect<double> &x = m_min;
    x = m_guess;
    vops::Mult<double> mult;
    vops::Subtr<double> subtr;
    
    Vect<double> s( x.Size() );
    
    m_f->Grad( x, s ); 
    
    m_iter = 0;
    while ( Continue(s, m_tolerance ) && m_iter < m_maxIter )
    {
      mult( s, m_gain );
      subtr( x, s, x );
      m_f->Grad( x, s );
      ++m_iter;
    }
    
    m_runCalc = false;
    if ( m_iter < m_maxIter )
    {
      m_minFound = true;
    }
  }
  
  bool Continue( Vect<double> &s, double tol )
  {
    for ( ULong i = 0; i < s.Size(); ++i )
    {
      if ( std::abs( s[i] ) > tol )
      {
        return true;
      }
    }
    return false;
  }
private:


#ifdef MYDEBUG
  friend class GradDescentTest;
#endif //MYDEBUG
};

}/*optim*/ }/*math*/ 

#endif /*_MATH_OPTIM_GRADDESCENT_h */
