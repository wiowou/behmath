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

#include "mesher.h"
#include "geom/shape.h"

namespace math{
namespace model{

MeshFace Mesher::MakeMFace(unsigned long long n0, unsigned long long n1, unsigned long long n2)
{
  MeshFace mf;
  mf.pointID.insert(n0);
  mf.point[0] = &point[n0];
  mf.pointID.insert(n1);
  mf.point[1] = &point[n1];
  mf.pointID.insert(n2);
  mf.point[2] = &point[n2];
  geom::Shape surf;
  surf.Set(mf.point);
  mf.norm = surf.PerpTri();
  mf.area = surf.AreaTri();
  mf.ctr = surf.Average(3);
  return mf;
}

MeshFace Mesher::MakeMFace(unsigned long long n0, unsigned long long n1, unsigned long long n2, unsigned long long n3)
{
  MeshFace mf;
  mf.pointID.insert(n0);
  mf.point[0] = &point[n0];
  mf.pointID.insert(n1);
  mf.point[1] = &point[n1];
  mf.pointID.insert(n2);
  mf.point[2] = &point[n2];
  mf.pointID.insert(n3);
  mf.point[3] = &point[n3];
  geom::Shape surf;
  surf.Set(mf.point);
  mf.norm = surf.PerpQuad();
  mf.area = surf.AreaQuad();
  mf.ctr = surf.Average(4);
  return mf;
}

void Mesher::CreateCellsFaces()
{
  faceBC.reserve(mFaceBC.size() );
  for ( const auto& mf : mFaceBC)
  {  
    auto it = mFace.find(mf);
    if ( it != mFace.end() ) it->pID = mf.pID;
  }
  
  face.resize( mFace.size() );
  unsigned long long i = 0;
  for ( const auto& mf : mFace )
  {
    face[i] = mf;
    if (face[i].pID != 0) faceBC.push_back(&face[i]);
    
    //found next blank entry in gf.cell[0]
    mf.cell[0]->nbor.emplace_back();
    if (mf.cell[1] != nullptr) //check if there is another cell attached to this face
    {
      //since there is another cell attached to this face, 
      //gf.cell[0] has that cell specified as a neighbor
      mf.cell[0]->nbor.back().cell = mf.cell[1];
    }
    //add the address of this face to gf.cell[0]'s faces
    mf.cell[0]->nbor.back().face = &face[i];
    
    if (mf.cell[1] != nullptr) //skip this gf if it is only attached to one cell
    {
      mf.cell[1]->nbor.emplace_back();
      mf.cell[1]->nbor.back().cell = mf.cell[0];
      mf.cell[1]->nbor.back().face = &face[i];
    }
    
    ++i;
  }

  mFaceBC.clear();
  mFace.clear();
}

}/*model*/ }/*math*/ 
