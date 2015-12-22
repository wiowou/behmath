#include "../src/vectValFunc.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h> //log, exp, sqrt functions
#include <cmath> //std:abs
#include <vector>

namespace math{

struct Surf1 : public Surface
{
  double operator()( Vect<double>& x )
  {
    return 16.0 * x[0] * x[0] * x[0] * x[0] + 16.0 * x[1] * x[1] * x[1] * x[1] + x[2] * x[2] * x[2] * x[2];
  }
};

struct Surf2 : public Surface
{
  double operator()( Vect<double>& x )
  {
    return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  }
};

struct Surf3 : public Surface
{
  double operator()( Vect<double>& x )
  {
    return x[0] * x[0] * x[0] - x[1];
  }
};

struct ExamSurf1 : public Surface
{
  double operator()( Vect<double>& x )
  {
    return 1.0 / tan( x[0] ) - 1.0 / ( x[0] * x[0] ) - 1.0 / x[1] - 1.0 / ( x[1] * x[1] );
  }
};

struct ExamSurf2 : public Surface
{
  double operator()( Vect<double>& x )
  {
    return x[0] * x[0] + x[1] * x[1];
  }
};

class VectValFuncTest
{
public:
  
  
  int Exec()
  {
    int i = 0;
    
    Vect<double> x(3);
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;

    Vect<double> b(3);
    b[0] = 16.0;
    b[1] = 3.0;
    b[2] = 0.0;

    Surf1 s1;
    Surf2 s2;
    Surf3 s3;
    
    std::cout << "showing how to find the vector output for input vector { 1, 2, 3}: " << std::endl;
    math::VectValFunc vvf(3);
    vvf.Comp(0, &s1);
    vvf.Comp(1, &s2);
    vvf.Comp(2, &s3);
    
    vvf.X(x);
    vvf.B(b);
    
    Vect<double>* out = vvf(x);
    std::cout << "f( {1, 2, 3} ) = " << std::endl << *out << std::endl;
    i += std::abs( (*out)[0] - 353.0 > 1e-3 );
    i += std::abs( (*out)[1] - 14.0 > 1e-3 );
    i += std::abs( (*out)[2] + 1.0 > 1e-3 );
  
    std::cout << " solving for f({x}) = { 16.0, 3.0, 0.0 } " << std::endl;
    //resetting b
    b[0] = 16.0;
    b[1] = 3.0;
    b[2] = 0.0;
    
    vvf.Tolerance(1e-5);
    vvf.MaxIter(100);
    vvf.Solve(); //solving for x, given b
    
    std::cout << "x: " << x << std::endl;
    i += std::abs( x[0] - 0.8780 > 1e-3 );
    i += std::abs( x[1] - 0.6768 > 1e-3 );
    i += std::abs( x[2] - 1.3309 > 1e-3 );
    
    //beginning of Exam 2, question 4
    vvf.Clear();
    vvf.Resize(2);
    b.Resize(2);
    //setting the surface equal to the given constant
    b[0] = 0.0;
    b[1] = 25.0;
    x.Resize(2);

    ExamSurf1 es1;
    ExamSurf2 es2;
    
    vvf.Comp(0, &es1);
    vvf.Comp(1, &es2);
    
    vvf.X(x);
    vvf.B(b);
    
    vvf.Tolerance(1e-5);
    vvf.MaxIter(50);
    
    Vect<double> xold(2);
    xold[0] = -10.0;
    xold[1] = -10.0;
    
    Vect<double> xguess(2);
    std::vector<Vect<double> > solutions;
    
    for ( xguess[0] = -4.99; xguess[0] < 5.01; xguess[0] += 0.05 )
    {
      for ( xguess[1] = -4.99; xguess[1] < 5.01; xguess[1] += 0.05 )
      {
        x = xguess;
        vvf.Solve(); //solving for x, given b
        if ( x[0] > -5.0 && x[0] < 5.0 &&  x[1] > -5.0 && x[1] < 5.0 && vvf.Iter() < vvf.MaxIter() )
        {
          bool match = false;
          for ( int i = 0; i < solutions.size(); ++i )
          {        
            if ( std::abs( x[0] - solutions[i][0] ) < 0.001 && std::abs( x[1] - solutions[i][1] ) < 0.001 )
            {
              match = true;
              break;
            }        
          }
          if ( !match )
          {
            solutions.push_back(x);
            std::cout << "x: " << x << std::endl;
          }          
        }
      }
    }
    
    
    
    
    return i;
  }

};

}/*math*/ 

int main()
{
  math::VectValFuncTest test;
  return test.Exec();
}
