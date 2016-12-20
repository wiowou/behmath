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

#include "gmshIO.h"
#include "geom/shape.h"
#include <fstream>
#include <sstream>

namespace math{
namespace model{

void GmshIO::ReadFVM( std::string fname, bool D3)
{
  AllocateStorage(fname, D3);
  std::ifstream fin;
  fin.open( fname.c_str(), std::ios::in);
  std::string line;
  bool elemLine = false;
  bool nodeLine = false;
  
  while (std::getline(fin, line) )
  {
    if ( line == "$Nodes" )
    {
      std::getline(fin, line);
      std::getline(fin, line);
      elemLine = false;
      nodeLine = true;
    }
    else if ( line == "$EndNodes" || line == "$EndElements" )
    {
      elemLine = false;
      nodeLine = false;
    }
    else if ( line == "$Elements" )
    {
      std::getline(fin, line);
      std::getline(fin, line);
      elemLine = true;
      nodeLine = false;
    }
    
    if (nodeLine) ReadNode(line);
    else if (elemLine)
    {
      std::stringstream ss(line);
      std::string token;
      if (!std::getline(ss, token, ' ') ) continue;
      if (!std::getline(ss, token, ' ') ) continue;
      std::stringstream ss2(token);
      int etype;
      ss2 >> etype;
      if (D3)
      {
        switch(etype)
        {
          case 2: ReadTriBC(line); break;
          case 3: ReadQuadBC(line); break;
          case 4: ReadTet(line); break;
          case 5: ReadHex(line); break;
          case 6: ReadPrism(line); break;
          case 7: ReadPyra(line); break;
        }
      }
      else
      {
      }
    }
  }
  fin.close();
  CreateCellsFaces();
}

void GmshIO::ReadTet( std::string& line)
{
  const unsigned long long nFaces = 4;
  const unsigned long long nNodes = 4;
  
  std::stringstream ss(line);
  std::vector<unsigned long long> id;
  
  unsigned long long val;
  while (ss >> val )
  {
    val--;
    id.push_back(val);
  }
  if (id[0] < el0) el0 = id[0];
  cell[id[0] - el0 ].pID = id[3]+1;
  if (cell[id[0] - el0 ].pID != 0) cellBC.push_back(&cell[id[0] - el0 ]);
  int s = id.size() - nNodes;
  MeshFace mf[nFaces];
  mf[0] = MakeMFace(id[s], id[s+1], id[s+2]);
  mf[1] = MakeMFace(id[s+1], id[s+3], id[s+2]);
  mf[2] = MakeMFace(id[s+2], id[s+3], id[s]);
  mf[3] = MakeMFace(id[s], id[s+3], id[s+1]);
  
  geom::Point* pts[nNodes];
  geom::Shape vol;
  for (int i = 0; i < nNodes; ++i)
  {
    pts[i] = &point[id[s+i] ];
    cell[id[0] - el0].pointID.push_back(id[s+i]);
  }
  vol.Set(pts);
  cell[id[0] - el0 ].ctr = vol.Average(nNodes);
  cell[id[0] - el0 ].vol = vol.VolTet();
  for (int i = 0; i < nFaces; ++i)
  {
    mf[i].cell[0] = &cell[id[0] - el0 ];
    auto it = mFace.find(mf[i]);
    if (it != mFace.end() ) 
    {
      it->cell[1] = &cell[id[0] - el0];
    }
    else mFace.insert(mf[i]);
  }
}

void GmshIO::ReadHex( std::string& line)
{
  const unsigned long long nFaces = 6;
  const unsigned long long nNodes = 8;
  
  std::stringstream ss(line);
  std::vector<unsigned long long> id;
  
  unsigned long long val;
  while (ss >> val )
  {
    val--;
    id.push_back(val);
  }
  if (id[0] < el0) el0 = id[0];
  cell[id[0] - el0 ].pID = id[3]+1;
  if (cell[id[0] - el0 ].pID != 0) cellBC.push_back(&cell[id[0] - el0 ]);
  int s = id.size() - nNodes;
  MeshFace mf[nFaces];
  mf[0] = MakeMFace(id[s], id[s+1], id[s+2], id[s+3]);
  mf[1] = MakeMFace(id[s+1], id[s+5], id[s+6], id[s+2]);
  mf[2] = MakeMFace(id[s], id[s+3], id[s+7], id[s+4]);
  mf[3] = MakeMFace(id[s+4], id[s+7], id[s+6], id[s+5]);
  mf[4] = MakeMFace(id[s], id[s+4], id[s+5], id[s+1]);
  mf[5] = MakeMFace(id[s+2], id[s+6], id[s+7], id[s+3]);
  
  geom::Point* pts[nNodes];
  geom::Shape vol;
  for (int i = 0; i < nNodes; ++i)
  {
    pts[i] = &point[id[s+i] ];
    cell[id[0] - el0].pointID.push_back(id[s+i]);
  }
  vol.Set(pts);
  cell[id[0] - el0 ].ctr = vol.Average(nNodes);
  cell[id[0] - el0 ].vol = vol.VolHex();
  for (int i = 0; i < nFaces; ++i)
  {
    mf[i].cell[0] = &cell[id[0] - el0 ];
    auto it = mFace.find(mf[i]);
    if (it != mFace.end() ) 
    {
      it->cell[1] = &cell[id[0] - el0];
    }
    else mFace.insert(mf[i]);
  }
}

void GmshIO::ReadPrism( std::string& line)
{
  const unsigned long long nFaces = 5;
  const unsigned long long nNodes = 6;
  
  std::stringstream ss(line);
  std::vector<unsigned long long> id;
  
  unsigned long long val;
  while (ss >> val )
  {
    val--;
    id.push_back(val);
  }
  if (id[0] < el0) el0 = id[0];
  cell[id[0] - el0 ].pID = id[3]+1;
  if (cell[id[0] - el0 ].pID != 0) cellBC.push_back(&cell[id[0] - el0 ]);
  int s = id.size() - nNodes;
  MeshFace mf[nFaces];
  mf[0] = MakeMFace(id[s], id[s+1], id[s+2]);
  mf[1] = MakeMFace(id[s+3], id[s+5], id[s+4]);
  mf[2] = MakeMFace(id[s], id[s+2], id[s+5], id[s+3]);
  mf[3] = MakeMFace(id[s], id[s+3], id[s+4], id[s+1]);
  mf[4] = MakeMFace(id[s+1], id[s+4], id[s+5], id[s+2]);
  
  geom::Point* pts[nNodes];
  geom::Shape vol;
  for (int i = 0; i < nNodes; ++i)
  {
    pts[i] = &point[id[s+i] ];
    cell[id[0] - el0].pointID.push_back(id[s+i]);
  }
  vol.Set(pts);
  cell[id[0] - el0 ].ctr = vol.Average(nNodes);
  cell[id[0] - el0 ].vol = vol.VolPrism();
  for (int i = 0; i < nFaces; ++i)
  {
    mf[i].cell[0] = &cell[id[0] - el0 ];
    auto it = mFace.find(mf[i]);
    if (it != mFace.end() ) 
    {
      it->cell[1] = &cell[id[0] - el0];
    }
    else mFace.insert(mf[i]);
  }
}

void GmshIO::ReadPyra( std::string& line)
{
  const unsigned long long nFaces = 5;
  const unsigned long long nNodes = 5;
  
  std::stringstream ss(line);
  std::vector<unsigned long long> id;
  
  unsigned long long val;
  while (ss >> val )
  {
    val--;
    id.push_back(val);
  }
  if (id[0] < el0) el0 = id[0];
  cell[id[0] - el0 ].pID = id[3]+1;
  if (cell[id[0] - el0 ].pID != 0) cellBC.push_back(&cell[id[0] - el0 ]);
  int s = id.size() - nNodes;
  MeshFace mf[nFaces];
  mf[0] = MakeMFace(id[s], id[s+1], id[s+2], id[s+3]);
  mf[1] = MakeMFace(id[s+1], id[s+4], id[s+2]);
  mf[2] = MakeMFace(id[s+2], id[s+4], id[s+3]);
  mf[3] = MakeMFace(id[s], id[s+3], id[s+4]);
  mf[4] = MakeMFace(id[s], id[s+4], id[s+1]);
  
  geom::Point* pts[nNodes];
  geom::Shape vol;
  for (int i = 0; i < nNodes; ++i)
  {
    pts[i] = &point[id[s+i] ];
    cell[id[0] - el0].pointID.push_back(id[s+i]);
  }
  vol.Set(pts);
  cell[id[0] - el0 ].ctr = vol.Average(nNodes);
  cell[id[0] - el0 ].vol = vol.VolPyra();
  for (int i = 0; i < nFaces; ++i)
  {
    mf[i].cell[0] = &cell[id[0] - el0 ];
    auto it = mFace.find(mf[i]);
    if (it != mFace.end() ) 
    {
      it->cell[1] = &cell[id[0] - el0];
    }
    else mFace.insert(mf[i]);
  }
}

void GmshIO::ReadTriBC( std::string& line)
{
  unsigned long long none = 0;
  --none;
  const unsigned long long nNodes = 3;
  
  std::stringstream ss(line);
  std::vector<unsigned long long> id;
  
  unsigned long long val;
  while (ss >> val )
  {
    val--;
    id.push_back(val);
  }
  if (id[3] == none) return; //isn't a member of a physical entity group
  int s = id.size() - nNodes;
  MeshFace mf;
  mf = MakeMFace(id[s], id[s+1], id[s+2]);
  mf.pID = id[3]+1;
  mFaceBC.push_back(mf);
}

void GmshIO::ReadQuadBC( std::string& line)
{
  unsigned long long none = 0;
  --none;
  const unsigned long long nNodes = 4;
  
  std::stringstream ss(line);
  std::vector<unsigned long long> id;
  
  unsigned long long val;
  while (ss >> val )
  {
    val--;
    id.push_back(val);
  }
  if (id[3] == none) return; //isn't a member of a physical entity group
  int s = id.size() - nNodes;
  MeshFace mf;
  mf = MakeMFace(id[s], id[s+1], id[s+2], id[s+3]);
  mf.pID = id[3]+1;
  mFaceBC.push_back(mf);
}

void GmshIO::ReadNode( std::string& line)
{
  std::stringstream ss(line);
  unsigned long long id;
  ss >> id;
  --id;
  for (int i = 0; i < 3; ++i)
  {
    ss >> point[id][i];
  }
}

void GmshIO::Clear()
{
  Domain::Clear();
  el0 = 0;
  --el0;
  mFace.clear();
  mFaceBC.clear();
}


void GmshIO::AllocateStorage( std::string fname, bool D3)
{
  Clear();
  std::ifstream fin;
  fin.open( fname.c_str(), std::ios::in);
  std::string line;
  bool elemLine = false;
  bool nodeLine = false;
  unsigned long long nodeCount = 0;
  unsigned long long elemCount = 0;
  
  while (std::getline(fin, line) )
  {
    if ( line == "$Nodes" )
    {
      std::getline(fin, line);
      std::getline(fin, line);
      elemLine = false;
      nodeLine = true;
    }
    else if ( line == "$EndNodes" || line == "$EndElements" )
    {
      elemLine = false;
      nodeLine = false;
    }
    else if ( line == "$Elements" )
    {
      std::getline(fin, line);
      std::getline(fin, line);
      elemLine = true;
      nodeLine = false;
    }
    
    if (nodeLine) ++nodeCount;
    else if (elemLine)
    {
      std::stringstream ss(line);
      std::string token;
      if (!std::getline(ss, token, ' ') ) continue;
      if (!std::getline(ss, token, ' ') ) continue;
      std::stringstream ss2(token);
      int etype;
      ss2 >> etype;
      if (D3 && etype > 3 && etype < 8) ++elemCount;
      else if (!D3 && etype > 1 && etype < 4) ++elemCount;
    }
  }
  
  point.resize(nodeCount);
  cell.resize(elemCount);
  //face.resize(6 * elemCount);
  fin.close();
}

}/*model*/ }/*math*/ 
