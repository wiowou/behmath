#include "../src/thermal.h"
#include "../src/gmshIO.h"
#include "../src/vtkIO.h"
#include "../src/bc/heatFlux.h"
#include "../src/bc/HTC.h"
#include "../src/bc/temperature.h"
#include <iostream>
#include <fstream>
#include <string>
#include "eng/materialData.h"
#include "eng/isoSolid.h"


namespace math{
namespace model{

class ThermalTest
{
public:
  
  int Exec()
  {
    int i = 0;
    
    std::string testDir = "/home/bk/prog/eng/modelv3/UnitTests";
    GmshIO gmshIO;
    gmshIO.ReadFVM(testDir + "/inputs/sensor.msh");
    
    eng::MaterialData unobtaniumData;
    i += unobtaniumData.LoadJSON( testDir + "/inputs/unobtanium.txt" );
    eng::IsoSolid unobtanium;
    unobtanium.Set(unobtaniumData);
    
    Thermal modelTh;
    modelTh.Swap(gmshIO);
    
    std::map<unsigned long long, eng::Material*> matMap;
    matMap.insert(std::make_pair(28, &unobtanium) );
    modelTh.AssignMaterial(matMap);
    
    bc::HeatFlux<0> radFlux;
    radFlux.pID = 27; //tying this bc to the physical id of the faces in the model
    radFlux.hflux.val = 10000.0;
    modelTh.boundaryCondition.push_back(&radFlux);
    
    bc::HTC<0> htcCooling;
    htcCooling.pID = 25;
    htcCooling.htc.val = 100.0;
    htcCooling.Tamb.val = 450.0;
    modelTh.boundaryCondition.push_back(&htcCooling);
    
    bc::Temperature<0> tbase;
    tbase.pID = 26;
    tbase.Tfix.val = 350.0;
    modelTh.boundaryCondition.push_back(&tbase);
    
    modelTh.SolveFVM();
    VtkIO vtkIO;
    vtkIO.dom = &modelTh;
    vtkIO.AddCellData(modelTh.ts, "temperature");
    vtkIO.Write(testDir + "/outputs/sensorTemps.vtk");
    
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
  math::model::ThermalTest test;
  return test.Exec();
}
