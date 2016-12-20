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

#include "vtkIO.h"
#include <fstream>
#include <sstream>

namespace math{
namespace model{

void VtkIO::AddCellData( Vect<double>& data, std::string desc)
{
  cellData.push_back(&data);
  cellDataDesc.push_back(desc);
}

void VtkIO::AddPointData( Vect<double>& data, std::string desc)
{
  pointData.push_back(&data);
  pointDataDesc.push_back(desc);
}

void VtkIO::Clear()
{
  cellData.clear();
  pointData.clear();
  cellDataDesc.clear();
  pointDataDesc.clear();
}

void VtkIO::Write(std::string fname, std::string title)
{
  std::ofstream fout;
  fout.precision(10);
  fout.open( fname.c_str(), std::ios::out);
  fout << "# vtk DataFile Version 1.0" << std::endl;
  fout << title << std::endl;
  fout << "ASCII" << std::endl;
  fout << std::endl;
  
  fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
  fout << "POINTS " << dom->point.size() << " double" << std::endl;
  for (unsigned long long i = 0; i < dom->point.size(); ++i)
  {
    fout << std::fixed << dom->point[i].X() << " ";
    fout << std::fixed << dom->point[i].Y() << " ";
    fout << std::fixed << dom->point[i].Z();
    fout << std::endl;
  }
  
  fout << "CELLS " << dom->cell.size() << " ";
  unsigned long long totSize = 0;
  for (unsigned long long i = 0; i < dom->cell.size(); ++i)
  {
    totSize += dom->cell[i].pointID.size() + 1;
  }
  fout << totSize << std::endl;
  for (unsigned long long i = 0; i < dom->cell.size(); ++i)
  {
    fout << dom->cell[i].pointID.size();
    for (int j = 0; j < dom->cell[i].pointID.size(); ++j)
    {
      fout << " " << dom->cell[i].pointID[j]; 
    }
    fout << std::endl;
  }
  
  fout << "CELL_TYPES " << dom->cell.size() << std::endl;
  for (unsigned long long i = 0; i < dom->cell.size(); ++i)
  {
    if (dom->cell[i].pointID.size() == 4) fout << 10;
    else if (dom->cell[i].pointID.size() == 8) fout << 12;
    else if (dom->cell[i].pointID.size() == 6) fout << 13;
    else if (dom->cell[i].pointID.size() == 5) fout << 14;
    else if (dom->cell[i].pointID.size() == 3) fout << 5;
    fout << std::endl;
  }
  
  if (cellData.size() > 0)
  {
    fout << "CELL_DATA " << dom->cell.size() << std::endl;
    for (int j = 0; j < cellData.size(); ++j)
    {
      fout << "SCALARS " << cellDataDesc[j] << " double" << std::endl;
      fout << "LOOKUP_TABLE default" << std::endl;
      Vect<double>& data = (*cellData[j]);
      for (int i = 0; i < data.Size(); ++i)
      {
        fout << std::fixed << data[i];
        fout << std::endl;
      }
    }
  }
  
  if (pointData.size() > 0)
  {
    fout << "POINT_DATA " << dom->point.size() << std::endl;
    for (int j = 0; j < pointData.size(); ++j)
    {
      fout << "SCALARS " << pointDataDesc[j] << " double" << std::endl;
      fout << "LOOKUP_TABLE default" << std::endl;
      Vect<double>& data = (*pointData[j]);
      for (int i = 0; i < data.Size(); ++i)
      {
        fout << std::fixed << data[i];
        fout << std::endl;
      }
    }
  }
  fout.close();
}

}/*model*/ }/*math*/ 
