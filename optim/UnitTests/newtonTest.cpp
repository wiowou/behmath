#include "../src/newton.h"
#include "curve/curve.h"
#include "surface/surface.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

namespace math{
namespace optim{

class Fa : public Curve
{
public:
  double operator()( double x )
  {
    return x * x * x * x - 14.0 * x * x * x + 60.0 * x * x - 70.0 * x;
  }
};

class Fb : public Curve
{
public:
  double operator()( double x )
  {
    return 0.5 * x * x - sin(x);
  }
};

class Fc : public Curve
{
public:
  double operator()( double x )
  {
    return x * x + 4.0 * cos(x);
  }
};

class Fd : public Curve
{
public:
  double operator()( double x )
  {
    double sum = 0.0;
    
    for ( double inc = 1e-6; inc < 100.1; inc *= 100.0 )
    {
      double f0 = f(inc, x);
      double f1;
      for ( double t = 2.0 * inc; t < 100.1 * inc; t += inc )
      {
        f1 = f(t,x);
        sum += 0.5 * ( f0 + f1 ) * inc;
        f0 = f1;
      }      
    }
    return sum;
  }
  
  double f( double t, double x )
  {
    return pow( t, x - 1.0 ) * exp( -t );
  }
};

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

class ExamSurf7 : public Surface
{
  double operator()( Vect<double>&x )
  {
    return 100.0 * ( x[1] - x[0] * x[0] ) * ( x[1] - x[0] * x[0] ) + ( 1.0 - 2.0 * x[0] ) * ( 1.0 - 2.0 * x[0] );
  }
};

class ExamSurf8 : public Surface
{
  double operator()( Vect<double>&x )
  {
    double d = 0.5 * x[0] * x[0] + x[1] * x[1]
      + std::abs( x[1] * x[1] - ( x[0] - 3.0 ) * ( x[0] - 3.0 ) * ( x[0] - 3.0 ) );
    if ( x[0] < -5.0 || x[1] < -5.0 )
    {
      x[0] = -5.0;
      return -d;
    }
    if ( x[0] > 5.0 || x[1] > 5.0 )
    {
      x[0] = 5.0;
      return d;
    }
    return d;
  }
};

class NewtonTest
{
public:
  
int Exec()
  {
    math::optim::Newton<Curve,double> nw1;
    int i = 0;

    Fa fa;
    Fb fb;
    Fc fc;
    Fd fd;

    double minx;
    double gain;
    double guess = 1.5;
    std::cout << "minimum around x=1.5 " << std::endl;
    
    std::cout << "function a: " << std::endl;
    nw1.Clear();
    nw1.Function(fa);
    nw1.Guess( guess );
    nw1.MaxIter(100);
    minx = nw1.MinX();
    std::cout << "iterations: " << nw1.Iter() << std::endl;
    std::cout << minx << std::endl;
    i += std::abs( minx - .78088 ) > 1e-4;

    std::cout << "function b: " << std::endl;
    nw1.Clear();
    nw1.Function(fb);
    nw1.Guess( guess );
    nw1.MaxIter(1000);
    minx = nw1.MinX();
    std::cout << "iterations: " << nw1.Iter() << std::endl;
    std::cout << minx << std::endl;   
    i += std::abs( minx - .739085 ) > 1e-4; 

    std::cout << "function c: " << std::endl;
    nw1.Clear();
    nw1.Function(fc);
    nw1.Guess( guess );
    nw1.MaxIter(1000);
    minx = nw1.MinX();
    std::cout << "iterations: " << nw1.Iter() << std::endl;
    std::cout << minx << std::endl;
    i += std::abs( minx - 1.89549 ) > 1e-4; 

    std::cout << "function d: " << std::endl;
    nw1.Clear();
    nw1.Function(fd);
    nw1.Guess( guess );
    nw1.MaxIter(1000);
    minx = nw1.MinX();
    std::cout << "iterations: " << nw1.Iter() << std::endl;
    std::cout << minx << std::endl;
    i += std::abs( minx - 1.49611 ) > 1e-4; 
    
    math::optim::Newton<Surface,Vect<double> > nw2;
    Vect<double> guess2(2);
    Vect<double>* minx2;
    F5 f5;
    
    guess2[0] = 1.5;
    guess2[1] = 0.5;
    std::cout << " function 5: " << std::endl;
    nw2.Clear();
    nw2.Function(f5);
    nw2.Guess(guess2);
    nw2.MaxIter(1000);
    minx2 = nw2.MinX();
    std::cout << "iterations: " << nw2.Iter() << std::endl;
    std::cout << "min point: " << *minx2 << std::endl;
    i += std::abs( (*minx2)[0] - 1.8685 ) > 1e-4;
    i += std::abs( (*minx2)[1] - 0.4343 ) > 1e-4;
    
    F6 f6;
    f6.DeltaX(1e-4);
    Vect<double> guess4(4);
    Vect<double>* minx4;

    guess4[0] = 1.0;
    guess4[1] = 1.0;
    guess4[2] = 1.0;
    guess4[3] = 0.2;
    
    std::cout << " function 6: " << std::endl;
    nw2.Clear();
    nw2.Function(f6);
    nw2.Guess(guess4);
    nw2.MaxIter(1000);
    minx4 = nw2.MinX();
    std::cout << "iterations: " << nw2.Iter() << std::endl;
    std::cout << "min point: " << *minx4 << std::endl;

    ExamSurf7 f7;
    f7.DeltaX(1e-4);
    Vect<double> guess5(2);
    Vect<double>* minx5;

    nw2.Clear();
    nw2.Function(f7);
    nw2.Guess(guess5);
    nw2.MaxIter(1000);
    
    std::cout << " exam problem 5, (-1,1): " << std::endl;
    guess5[0] = -1.0;
    guess5[1] = 1.0;    
    minx5 = nw2.MinX();
    std::cout << "iterations: " << nw2.Iter() << std::endl;
    std::cout << "min point: " << *minx5 << std::endl;
    i += std::abs( (*minx2)[0] - 0.5 ) > 1e-4;
    i += std::abs( (*minx2)[1] - 0.25 ) > 1e-4;

    nw2.Clear();
    nw2.Function(f7);
    nw2.Guess(guess5);
    nw2.MaxIter(1000);
    
    std::cout << " exam problem 5, (0,1): " << std::endl;
    guess5[0] = 0.0;
    guess5[1] = 1.0;    
    minx5 = nw2.MinX();
    std::cout << "iterations: " << nw2.Iter() << std::endl;
    std::cout << "min point: " << *minx5 << std::endl;
    i += std::abs( (*minx2)[0] - 0.5 ) > 1e-4;
    i += std::abs( (*minx2)[1] - 0.25 ) > 1e-4;

    nw2.Clear();
    nw2.Function(f7);
    nw2.Guess(guess5);
    nw2.MaxIter(1000);
    
    std::cout << " exam problem 5, (2,1): " << std::endl;
    guess5[0] = 2.0;
    guess5[1] = 1.0;    
    minx5 = nw2.MinX();
    std::cout << "iterations: " << nw2.Iter() << std::endl;
    std::cout << "min point: " << *minx5 << std::endl;
    i += std::abs( (*minx2)[0] - 0.5 ) > 1e-4;
    i += std::abs( (*minx2)[1] - 0.25 ) > 1e-4;
    
    ExamSurf8 f8;
    
    nw2.Clear();
    nw2.Function(f8);
    nw2.Guess(guess5);
    nw2.MaxIter(1000);
    
    std::cout << " exam problem 6, guess (0,0): " << std::endl;
    guess5[0] = 0.0;
    guess5[1] = 0.0;    
    minx5 = nw2.MinX();
    std::cout << "iterations: " << nw2.Iter() << std::endl;
    std::cout << "min point: " << *minx5 << std::endl;
    i += std::abs( (*minx2)[0] - 2.1529 ) > 1e-4;
    i += std::abs( (*minx2)[1] - 0.0 ) > 1e-4;
    
    return i;
  }

};

}/*optim*/ }/*math*/ 

int main()
{
  math::optim::NewtonTest test;
  return test.Exec();
}
