#include "../src/conjGrad.h"
#include "surface/surface.h"
#include <iostream>
#include <fstream>
#include <string>

#include <cmath> //std::abs

namespace math{
namespace optim{

class F5 : public Surface
{
  double operator()( Vect<double>&x )
  {
    return 2.0 * x[0] * x[0] * x[0] - 4.0 * x[0] * x[0] - 6.0 * x[0] * x[1] * ( x[0] - x[1] - 1.0 );
  }
};

class F6 : public Surface
{
  double operator()( Vect<double>&x )
  {
    double d = -5.0 * x[0] * x[1] * x[3] + 4.0 * x[0] * pow( x[1], 1.4 ) + 0.75 * pow( x[2], 0.6 )
      + 10.0 * std::abs( x[0] * x[3] - 8.4 * x[1] * x[2] * ( 1.0 - x[3] ) * ( 1.0 - x[3] ) );
    
    for ( int i = 0; i < 4; ++i )
    {
      if ( x[i] < 0.0 )
      {
        x[i] = 0.0;
        return -10.0 * x[i];
      }
    }
    if ( x[3] > 1.0 )
    {
      x[3] = 1.0 * x[3];
      return 10.0;
    }
    return d;
  }
};

class ConjGradTest
{
public:
  
  int Exec()
  {
    int i = 0;
    
    math::optim::ConjGrad<Surface> gd2;
    Vect<double> gain2(2);
    Vect<double> guess2(2);
    Vect<double>* minx2;
    F5 f5;
    
    for ( ULong i = 0; i < 2; ++i )
    {
      gain2[i] = 0.02;
      guess2[i] = 0.5;
    }
    std::cout << " function 5: " << std::endl;
    gd2.Clear();
    gd2.Function(f5);
    gd2.Guess(guess2);
    gd2.Gain(gain2);
    gd2.MaxIter(1000);
    minx2 = gd2.MinX();
    std::cout << "iterations: " << gd2.Iter() << std::endl;
    std::cout << "min point: " << *minx2 << std::endl;

    i += std::abs( (*minx2)[0] - 1.8685 ) > 1e-4;
    i += std::abs( (*minx2)[1] - 0.4343 ) > 1e-4;
    i += gd2.Iter() > 240;
    
    F6 f6;
    f6.DeltaX(1e-8);
    Vect<double> gain4(4);
    Vect<double> guess4(4);
    Vect<double>* minx4;
    for ( ULong i = 0; i < 4; ++i )
    {
      gain4[i] = 0.005;
    }
    guess4[0] = 1.0;
    guess4[1] = 2.0;
    guess4[2] = 3.0;
    guess4[3] = 0.0;
    
    std::cout << " function 6: " << std::endl;
    gd2.Clear();
    gd2.Function(f6);
    gd2.Guess(guess4);
    gd2.Gain(gain4);
    gd2.MaxIter(1000);
    minx4 = gd2.MinX();
    std::cout << "iterations: " << gd2.Iter() << std::endl;
    std::cout << "min point: " << *minx4 << std::endl;
    
    return i;
  }

};

}/*optim*/ }/*math*/ 

int main()
{
  math::optim::ConjGradTest test;
  return test.Exec();
}
