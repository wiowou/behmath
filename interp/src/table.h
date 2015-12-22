//<license>

//</license>

#ifndef _MATH_TABLE_h
#define _MATH_TABLE_h

#include "impl/interpConfig.h"

namespace math{

typedef unsigned long long Uint;
typedef long long Sint;

#ifdef MYDEBUG
  class TableTest;
#endif //MYDEBUG

/*! Class to perform linear interpolations in up to 3 dimensions.
The class is not responsible for freeing or allocating any memory. The user
must do this.
*/
class Table
{

public:
  Table() = default;
  //! Constructor that sets the data table
  explicit Table( double* d );
  /*!Set the data table
  
            0: number of dimensions, ndim
            1: number of values in dimension 1, nv1
            2: number of values in dimension 2, nv2
            ...
            ndim: number of values in dimension ndim, nvndim
            ndim+1: 1st value for dimension 1
            ndim+nv1: last value for dimension 1
            ndim+nv1+1: 1st value for dimension 2
            ndim+nv1+nv2: last value for dimension 2
            ...
            ndim+nv1+nv2...+nvndim+1: start of the dependent values
  */
  void Data( double* d );
  //! Return the size of the entire data table
  Uint Size();
  //! 1D, linear, orthogonal interpolation
  double operator()( double d1 );
  //! 2D, linear, orthogonal interpolation
  double operator()( double d1, double d2);
  //! 3D, linear, orthogonal interpolation
  double operator()( double d1, double d2, double d3 );

protected:
  Sint UpperBoundIdx( double* beg, double* end, const double val );
  double Interp1( double size[], double valI1[], double valD[], double d1 );
  double Interp2( double size[], double valI1[], double valI2[], double valD[], double d1, double d2 );
  double Interp3( double size[], double valI1[], double valI2[], double valI3[], double valD[], double d1, double d2, double d3 );

private:

  double* m_d;


#ifdef MYDEBUG
  friend class TableTest;
#endif //MYDEBUG
};

}/*math*/ 

#endif /*_MATH_TABLE_h */
