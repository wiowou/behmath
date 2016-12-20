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

#ifndef _ORTHOLINTERP_h
#define _ORTHOLINTERP_h

#include <vector>
#include <algorithm>
#include <cmath>

namespace math{

#ifdef MYDEBUG
    class OrthoLinTerpTest;
#endif //MYDEBUG

//! \typedef VecVec is a std::vector< std::vector<double> >
typedef std::vector< std::vector<double> > VecVec;

//! This class provides a way to interpolate on a regular, orthogonal grid of independent variables.
/*!
  /param ndim is the number of dimensions of the independent variable. It supports an n-dimensional
  interpolation. The interpolation method in 1D is specialized (about 4x faster than the general method).
  Specializations in 2D should be considered. Element shape functions using natural coordinates are used
  as the methodology behind the interpolation. A reference is provided in the doc folder:
  NatCoordShapeFunc.pdf, page 11-5, formula 11.3
*/
template < std::size_t ndim > 
class OrthoLinTerp
{

public:
  OrthoLinTerp()
  {
    m_pow2ndim = 1 << ndim; // = 2^ndim
    m_shapeFuncCoeff = 1.0 / static_cast<double>(m_pow2ndim); //use this coeff to multiply the result against at the end.
  }
  
  //! DataTable provides write access to \param m_dataTable.
  /*!
  \param dataTable is a vector<vector<double> > of size ndim + 1 whose 0->ndim columns are the independent variables
  and whose ndim column is the dependent variable. For example:
  
      x00 x10 x20 v00
      x01 x11 x21 v01
      x02 x12     v02
          x13     v03
                  v04
                  ...
                  v23

  where v is of length x0*x1*x2
  */
  void DataTable( const VecVec &dataTable )
  {
    m_dataTable = &dataTable;
    CalcJumpSizes( dataTable );
  }

  //! \param v is an input vector of the variables to interpolate on in each dimension.
  /*! \param v is of length ndim. It provides a general method for an n dimensional regular table interpolation.
  */
  double operator()( const std::vector<double> &v )
  {
    return Interp( v.data(), *m_dataTable, false );
  }
  
  //! provides a method that can be used in convenient form to do the interpolation for up to ndim = 5.
  double operator()( const double x1, const double x2 = 0.0, const double x3 = 0.0, const double x4 = 0.0, const double x5 = 0.0 )
  {
    switch(ndim)
    {
      case 1: return Interp1( x1, *m_dataTable );
      case 2: m_x[0] = x1; m_x[1] = x2; break;
      case 3: m_x[0] = x1; m_x[1] = x2; m_x[2] = x3; break;
      case 4: m_x[0] = x1; m_x[1] = x2; m_x[2] = x3; m_x[3] = x4; break;
      case 5: m_x[0] = x1; m_x[1] = x2; m_x[2] = x3; m_x[3] = x4; m_x[4] = x5; break;
    }
    return Interp( m_x, *m_dataTable, false );
  }
  
protected:

  //! The general interpolation routine based on element shape functions.
  /*! \param v[] is an array of the independent input variables that will be used to find
      the factors from \param vv the vector of vectors that includes the independent and 
      dependent variable. \param calcJumpSizes is set to true when the jump sizes, used to 
      navigate through the dependent variable vector, need to be calculated.
  */
  double Interp( const double v[], const VecVec &vv, bool calcJumpSizes )
  {
    //finding the indexes of the lower and upper bound for each independent paramater, v[i]
    //find the natural coordinate for each independent variable
    for ( int i = 0; i < ndim; ++i )
    {
      m_up  = std::upper_bound( vv[i].begin(), vv[i].end(), v[i] ); // <algorithm>
      if ( m_up == vv[i].begin() ) 
      {
        m_idx[0][i] = m_idx[1][i] = 0;
        m_eta[i] = -1.0;
      }
      else if ( m_up == vv[i].end() ) 
      {
        m_idx[0][i] = m_idx[1][i] = vv[i].size() - 1;
        m_eta[i] = 1.0;
      }
      else
      {
        m_idx[1][i] = m_up - vv[i].begin();
        m_idx[0][i] = m_idx[1][i] - 1;
        //multiplying by 2 and subtracting 1 converts the ratio to natural coordinates
        m_eta[i] = 2.0 * ( v[i] - vv[i][ m_idx[0][i] ] ) / ( vv[i][ m_idx[1][i] ] - vv[i][ m_idx[0][i] ] ) - 1.0; 
      }
    }
    // done finding bounds, natural coordinates
    
    //convenience calc to determine the size of the jumps needed to navigate the dependent variable vector, vv[ndim]
    if ( calcJumpSizes )
    {
      CalcJumpSizes(vv);
    }
    // done finding jump sizes
    
    //calculate the 2^ndim shape functions. Each function has ndim coordinates that are multiplied together
    // formula for the shape functions in n dimensions comes from NatCoordShapeFunc.pdf, page 11-5, formula 11.3
    double result = 0.0;
    for ( int i = 0; i < m_pow2ndim; ++ i )
    {
      std::size_t b[ndim];
      ToBitArray( i, b );
      std::size_t resIdx = 0;
      double shapeFuncVal = 1.0;
      for ( int k = 0; k < ndim; ++k )
      {
        resIdx += m_idx[ b[k] ][k] * m_jmp[k]; //b[k] should resolve to 0 for false and 1 for true if the bit is set.
        shapeFuncVal *= 1.0 + ( 2.0 * static_cast<double>(static_cast<int>( b[k] ) ) - 1.0 ) * m_eta[k];
      }
      shapeFuncVal *= vv[ndim][resIdx];
      result += shapeFuncVal;
    }
    return result * m_shapeFuncCoeff;
  }

  //! Calculates the jump sizes needed to navigate around the dependent variable vector
  /*! For the example \param dataTable in \fn DataTable, the jump sizes would be 
      8=4*2 for the first column, 2 for the second column, and 1 for the third.
  */
  void CalcJumpSizes( const VecVec &vv )
  {
    for ( std::size_t n = 1; n < ndim; ++n )
    {
      std::size_t m = n - 1;
      m_jmp[m] = 1;
      for ( std::size_t i = n; i < ndim; ++i )
      {
        m_jmp[m] *= vv[i].size();
      }
    }
    m_jmp[ndim - 1] = 1;
  }
  
  void ToBitArray( std::size_t ival, std::size_t b[] )
  {
    std::size_t mask = 1;
    for( std::size_t i = 0; i < ndim; ++i )
    {
      b[i] = static_cast<std::size_t>( static_cast<bool>(ival & mask) );
      mask  <<= 1;
    }
  }

  //! The common 1D interpolation routine, 4x faster than the general routine.
  /*! \param v is the independent variable. \param vv[0] is the vector of the independent variable
      \param v[1] is the vector of the dependent variable.
  */
  double Interp1( double v, const VecVec &vv )
  {
    m_up  = std::upper_bound( vv[0].begin(), vv[0].end(), v );
    if ( m_up == vv[0].begin() ) 
    {
      return vv[1].front();
    }
    if ( m_up == vv[0].end() ) 
    {
      return vv[1].back();
    }
    std::size_t ind = m_up - vv[0].begin();
    double  yr =  ( v - vv[0][ind - 1] ) / ( vv[0][ind] - vv[0][ind - 1] );
    return yr * ( vv[1][ind] - vv[1][ind - 1] ) + vv[1][ind - 1];  
  }
  
  const VecVec *m_dataTable;
  std::size_t m_jmp[ndim];
  //! Stores the pre and post indexes of the independent variable. First col is pre, second is post.
  std::size_t m_idx[2][ndim];
  //! The natural coordinate in each dimension.
  double m_eta[ndim];
  std::size_t m_pow2ndim;
  double m_shapeFuncCoeff;
  std::vector<double>::const_iterator m_up;

private:
  double m_x[5];

#ifdef MYDEBUG
    friend class OrthoLinTerpTest;
#endif //MYDEBUG
};

template <> 
class OrthoLinTerp<0>
{
public:
  double val;
  //! \param v is an input vector of the variables to interpolate on in each dimension.
  /*! \param v is of length ndim. It provides a general method for an n dimensional regular table interpolation.
  */
  double operator()( const std::vector<double> &v )
  {
    return val;
  }
  
  //! provides a method that can be used in convenient form to do the interpolation for up to ndim = 5.
  double operator()( const double x1, const double x2 = 0.0, const double x3 = 0.0, const double x4 = 0.0, const double x5 = 0.0 )
  {
    return val;
  }
};

}/*math*/ 

#endif /*_ORTHOLINTERP_h */
