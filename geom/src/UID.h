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

#ifndef _MATH_GEOM_UID_h
#define _MATH_GEOM_UID_h

namespace math{
namespace geom{

#ifdef MYDEBUG
  class UIDTest;
#endif //MYDEBUG

template <class T = void>
class UID
{

public:
	UID()
	{
		m_id = m_count++;
	}
	
	UID(const UID &u)
	{
		m_id = m_count++;
	}
	
	UID& operator=(UID u)
	{
		return *this;
	}
	
	unsigned long long ID() const
	{
		return m_id;
	}
	
	void ID( unsigned long long id )
	{
		m_id = id;
	}

  
protected:
  unsigned long long m_id;
	static unsigned long long m_count;
	
private:


#ifdef MYDEBUG
  friend class UIDTest;
#endif //MYDEBUG
};

template <class T>
unsigned long long UID<T>::m_count = 0;

}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_UID_h */
