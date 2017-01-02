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

#ifndef _MATH_GEOM_POINT_h
#define _MATH_GEOM_POINT_h

namespace math{
namespace geom{

#ifdef MYDEBUG
  class PointTest;
#endif //MYDEBUG

class Point
{

public:
  Point();
  
  Point( double x, double y, double z );
  
  void Set( double x, double y, double z );
  
  double* Crd();
  
  double& operator[]( int i );
  
  double X() const;
  void X( double d );
  
  double Y() const;
  void Y( double d );
  
  double Z() const;
  void Z( double d );
  
  Point& operator-=( const Point &rhs );
  Point& operator+=( const Point &rhs );
  Point& operator*=( const double rhs );
  Point& operator/=( const double rhs );
  friend Point operator-( const Point &lhs, const Point &rhs );
  friend Point operator+( const Point &lhs, const Point &rhs );
  friend Point operator*( const double lhs, const Point &rhs );
  friend Point operator*( const Point &lhs, const double rhs );
  friend Point operator/( const Point &lhs, const double rhs );
  //virtual ~Point()
  //{}
protected:
  double m_crd[3];
  
private:


#ifdef MYDEBUG
  friend class PointTest;
#endif //MYDEBUG
};

extern Point operator-( const Point &lhs, const Point &rhs );
extern Point operator+( const Point &lhs, const Point &rhs );
extern Point operator*( const double lhs, const Point &rhs );
extern Point operator*( const Point &lhs, const double rhs );
extern Point operator/( const Point &lhs, const double rhs );
extern double Dist( const Point &lhs, const Point &rhs );

extern const Point origin;
extern const Point* porigin;
  
}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_POINT_h */
