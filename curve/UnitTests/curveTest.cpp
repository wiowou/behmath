#define MYDEBUG

#include "../src/curve.h"
#include "../src/polynomial.h"
#include "../src/power.h"
#include "../src/exponential.h"
#include "matrix/cholesky.h"
#include "matrix/QR.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{

class CurveTest
{
public:
  
  
  int Exec()
  {
    int i = 0;
    //create a polynomial object of order 2, a quadratic
    curve::Polynomial parabola(2);
    
    //set the x and y data for the sin function between 0 and pi
    const double pi = 3.141592654;
    const ULong npts = 21;
    Vect<double> x(npts), y(npts);
    for ( int i = 0; i < npts; ++i )
    {
      x[i] = pi / static_cast<double>(npts - 1) * static_cast<double>(i);
      y[i] = sin(x[i]);
    }
    //print x and y
    std::cout << "x: " << x << std::endl;
    std::cout << "y: " << y << std::endl;
    //pass x and y to the polynomial object's Fit method. Indicate that the Cholesky least squares solver be used
    //since a symmetric Vandermonde matrix is being inverted.
    parabola.Fit<solver::Cholesky>(x,y);
    
    //view the equation generated and the error
    std::cout << "equation of parabola: " << std::endl;
    std::cout << parabola << std::endl;
    std::cout << "error estimate: " << parabola.Error() << std::endl;
    //demonstration of how to get a value out of the curve object
    std::cout << " calculate the value of the curve at pi/3: " << std::endl;
    std::cout << parabola(pi/3.0) << std::endl;
    
    i += CheckPair( parabola.m_coeff, 0, -.03673 );
    i += CheckPair( parabola.m_coeff, 1, 1.2902 );
    i += CheckPair( parabola.m_coeff, 2, -.4107 );
    
    i += ( parabola(pi/3.0) < .863) || ( parabola(pi/3.0) > .865 );  
    
    std::cout << " problem 4, hw3, using the QR solver for problem 2: " << std::endl;
    
    // reset the x and y data
    for ( int i = 0; i < npts; ++i )
    {
      x[i] = pi / static_cast<double>(npts - 1) * static_cast<double>(i);
      y[i] = sin(x[i]);
    }
    //print x and y
    std::cout << "x: " << x << std::endl;
    std::cout << "y: " << y << std::endl;
    
    //reset the curve coefficients
    parabola.Clear();
    
    //use the QR solver
    parabola.Fit<solver::QR>(x,y);
    
    //view the equation generated and the error
    std::cout << "equation of parabola: " << std::endl;
    std::cout << parabola << std::endl;
    std::cout << "error estimate: " << parabola.Error() << std::endl;

    i += CheckPair( parabola.m_coeff, 0, -.03673 );
    i += CheckPair( parabola.m_coeff, 1, 1.2902 );
    i += CheckPair( parabola.m_coeff, 2, -.4107 );
    
    double deriv = parabola.Deriv1(1.0);
    
    i += ( std::abs( deriv - 0.4688 ) > 0.001 );
    
    deriv = parabola.Deriv2(1.0);
    
    i += ( std::abs( deriv + 0.8214 ) > 0.001 );
    
    return i;
  }

    int CheckPair( SparseVect<double> &sv, ULong pos, double d )
  {
    if ( std::abs( (sv[pos] - d) / d) > 0.001 )
    {
      std::cout << pos << " bad data " << sv[pos] << "!=" << d << std::endl;
      return 1;
    }
    return 0;
  }
  int CheckPair( Vect<double> &v, ULong pos, double d )
  {
    if ( std::abs(d) < 1e-10 )
    {
      if ( v[pos] > 1e-10 )
      {
        std::cout << pos << " bad data " << v[pos] << "!=" << d << std::endl;
        return 1;
      }
      return 0;
    }
    if ( std::abs(v[pos] - d) / d > 0.0001 )
    {
      std::cout << pos << " bad data " << v[pos] << "!=" << d << std::endl;
      return 1;
    }
    return 0;
  } 
};

}/*math*/ 

int main()
{
  math::CurveTest test;
  return test.Exec();
}
