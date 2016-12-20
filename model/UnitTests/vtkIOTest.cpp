#include "../src/vtkIO.h"
#include <iostream>
#include <fstream>
#include <string>

namespace math{
namespace model{

class VtkIOTest
{
public:
  VtkIO vtkIO;
  Domain dom;
  
  int Exec()
  {
    std::string testDir = "/home/bk/prog/eng/modelv3/UnitTests";
    int i = 0;
    dom.point.push_back(geom::Point(0.0, 0.0, 0.0) );
    dom.point.push_back(geom::Point(1.0, 0.0, 0.0) );
    dom.point.push_back(geom::Point(2.0, 0.0, 0.0) );
    dom.point.push_back(geom::Point(0.0, 1.0, 0.0) );
    dom.point.push_back(geom::Point(1.0, 1.0, 0.0) );
    dom.point.push_back(geom::Point(2.0, 1.0, 0.0) );
    dom.point.push_back(geom::Point(0.0, 2.0, 0.0) );
    dom.point.push_back(geom::Point(1.0, 2.0, 0.0) );
    
    dom.cell.push_back(Cell());
    dom.cell.back().pointID.push_back(0);
    dom.cell.back().pointID.push_back(1);
    dom.cell.back().pointID.push_back(3);
    
    dom.cell.push_back(Cell());
    dom.cell.back().pointID.push_back(1);
    dom.cell.back().pointID.push_back(4);
    dom.cell.back().pointID.push_back(3);
    
    dom.cell.push_back(Cell());
    dom.cell.back().pointID.push_back(1);
    dom.cell.back().pointID.push_back(2);
    dom.cell.back().pointID.push_back(4);
    
    dom.cell.push_back(Cell());
    dom.cell.back().pointID.push_back(2);
    dom.cell.back().pointID.push_back(5);
    dom.cell.back().pointID.push_back(4);
    
    dom.cell.push_back(Cell());
    dom.cell.back().pointID.push_back(3);
    dom.cell.back().pointID.push_back(4);
    dom.cell.back().pointID.push_back(6);
    
    dom.cell.push_back(Cell());
    dom.cell.back().pointID.push_back(4);
    dom.cell.back().pointID.push_back(7);
    dom.cell.back().pointID.push_back(6);
    
    dom.ts.Resize(6);
    dom.ts[0] = 0.0;
    dom.ts[1] = 1.0;
    dom.ts[2] = 2.0;
    dom.ts[3] = 3.0;
    dom.ts[4] = 4.0;
    dom.ts[5] = 5.0;
    
    vtkIO.dom = &dom;
    vtkIO.AddCellData(dom.ts, "temperature");
    vtkIO.Write(testDir + "/outputs/test1.vtk");
    
    return i;
  }
  
  int FileMatch( std::string outFileName, std::string expFileName )
  {
    std::ifstream out, expected;
    std::string lineOut, lineExpected;
    out.open( outFileName );
    if ( out.fail() )
    {
      std::cout << "Failed open" << outFileName << std::endl;
      return 1;
    }
    expected.open( expFileName );
    if ( expected.fail() )
    {
      std::cout << "Failed open" << expFileName << std::endl;
      return 1;
    }
    while ( std::getline( expected, lineExpected ) )
    {
      std::getline( out, lineOut );
      if ( lineOut != lineExpected )
      {
        std::cout << "Failed " << outFileName << std::endl;
        return 1;
      }
    }
    if ( std::getline( out, lineOut ) )
    {
      std::cout << "Failed truncate " << outFileName << std::endl;
      return 1;
    }
    return 0;
  }

};

}/*model*/ }/*math*/ 

int main()
{
  math::model::VtkIOTest test;
  return test.Exec();
}
