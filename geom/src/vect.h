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

#ifndef _MATH_GEOM_VECT_h
#define _MATH_GEOM_VECT_h

#include "point.h"

namespace math{
namespace geom{

#ifdef MYDEBUG
  class VectTest;
#endif //MYDEBUG


class Vect : public Point
{

public:
  Vect();
  Vect( double x, double y, double z );
  Vect( const Point* head );
  Vect( const Point& head );
  Vect( const Point* head, const Point* tail );
  Vect( const Point& head, const Point& tail );
  
  void Set( const Point* head, const Point* tail );
  void Set( const Point* head );
  void Set( const Point& head, const Point& tail );
  void Set( const Point& head );
  
  double Mag() const;
  double PerpDistance(Point* p) const;
  double PerpDistance(Point& p) const;
  
  Vect& operator=( const Point& rhs );
  Vect& operator+=( const Vect& rhs );
  Vect& operator-=( const Vect& rhs );
  Vect& operator*=( const double rhs );
  Vect& operator/=( const double rhs );
  friend Vect operator+( const Vect &lhs, const Vect &rhs );
  friend Vect operator-( const Vect &lhs, const Vect &rhs );
  friend Vect operator*( const double lhs, const Vect &rhs );
  friend Vect operator*( const Vect &lhs, const double rhs );
  friend Vect operator/( const Vect &lhs, const double rhs );
  friend Vect Cross( const Vect &a, const Vect &b );
  friend double Dot( const Vect &a, const Vect &b );
  friend double TripleProd( const Vect &a, const Vect &b, const Vect &c );
  
protected:
  
private:


#ifdef MYDEBUG
  friend class VectTest;
#endif //MYDEBUG
};

extern Vect operator+( const Vect &lhs, const Vect &rhs );
extern Vect operator-( const Vect &lhs, const Vect &rhs );
extern Vect operator*( const double lhs, const Vect &rhs );
extern Vect operator*( const Vect &lhs, const double rhs );
extern Vect operator/( const Vect &lhs, const double rhs );
extern Vect Cross( const Vect &a, const Vect &b );
extern double Dot( const Vect &a, const Vect &b );
extern double TripleProd( const Vect &a, const Vect &b, const Vect &c );

//! The scalar projection of a onto b, same as |a| * cos (theta), theta is angle between a and b
extern double ScalarProj( const Vect &a, const Vect &b );

//! The vector projection of a onto b, same as |a| * cos (theta) / |b| * b
extern Vect VectorProj( const Vect &a, const Vect &b );
}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_VECT_h */
