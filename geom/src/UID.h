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

#include <vector>

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
    Increment();
    //m_id = m_count++;
	}
	
	UID(const UID &u)
	{
    Increment();
		//m_id = m_count++;
	}
	
	UID& operator=(UID u)
	{
		return *this;
	}
	
  ~UID()
  {
    m_pobj[m_id] = nullptr;
    m_deleteCount++;
  }
  
  void swap(UID &u)
  {
    unsigned long long tmp_id = m_id;
    m_id = u.m_id;
    u.m_id = tmp_id;
  }
  
	unsigned long long ID() const
	{
		return m_id;
	}
	
	void ID( unsigned long long id )
	{
    m_id = id;
	}
    
  static void Compress()
  {
    std::vector<UID<T>*> pobj;
    unsigned long long size = m_pobj.size() - m_deleteCount;
    pobj.resize(size);
    unsigned long long j = 0;
    for (unsigned long long i = 0; i < m_pobj.size(); ++i)
    {
      if (m_pobj[i] != nullptr)
      {
        m_pobj[i]->m_id = j;
        pobj[j] = m_pobj[i];
        ++j;
      }
    }
    m_pobj.swap(pobj);
    m_count = size;
    m_deleteCount = 0;
  }

  static unsigned long long Count()
  {
    return m_count;
  }

  friend bool operator<(const UID &lhs, const UID &rhs)
  {
    return lhs.m_id < rhs.m_id;
  }
  
protected:
  //unsigned long long m_id;
  //an iterator to the element before the one in the object
	unsigned long long m_id;
  static unsigned long long m_count;
  static unsigned long long m_deleteCount;
  static std::vector<UID<T>*> m_pobj;
  
private:
  void Increment()
  {
    m_id = m_count++;
    m_pobj.push_back(this);
  }

#ifdef MYDEBUG
  friend class UIDTest;
#endif //MYDEBUG
};

template <class T>
unsigned long long UID<T>::m_count = 0;

template <class T>
unsigned long long UID<T>::m_deleteCount = 0;

template <class T>
std::vector<UID<T>*> UID<T>::m_pobj;

}/*geom*/ }/*math*/ 

#endif /*_MATH_GEOM_UID_h */
