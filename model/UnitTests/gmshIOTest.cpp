#include "../src/gmshIO.h"
#include <iostream>
#include <fstream>
#include <string>


std::ofstream fout, fout2, fout3, fout4;
std::ifstream in[3], expected[3];
    
namespace math{
namespace model{

class GmshIOTest
{
public:
  math::model::GmshIO gmshIO;
  std::string testDir;
  
  int Exec()
  {
    testDir = "/home/bk/prog/eng/modelv3/UnitTests";

    int i = 0;
    
    gmshIO.ReadFVM( testDir + "/inputs/box.msh");
    
    i += !( gmshIO.cell.capacity() == 24 );
    i += !( gmshIO.point.capacity() == 14 );
    //i += !( gmshIO.face.capacity() == 24 * 6 );
    
    //WriteGMFaces();
    //i += FileMatch();
    WriteNodes();
    WriteCells();
    WriteFaces();
    
    for (int j = 0; j < 3; ++j) i += FileMatch(j);
    return i;
  }
  
  /*
  void WriteGMFaces()
  {
    for ( const auto& gf : gmshIO.gmFace )
    {
      for ( const auto& id : gf.pointID )
      {
        fout << id << " ";
      }
      fout << "nbor: ";
      if (gf.cell[0] != nullptr ) fout << (gf.cell[0] - &gmshIO.cell[0]);
      else fout << -1;
      fout << " ";
      if (gf.cell[1] != nullptr ) fout << (gf.cell[1] - &gmshIO.cell[0]);
      else fout << -1;
      fout << std::endl;
    }
    fout.close();
  }
  */
  void WriteNodes()
  {
    fout2.precision(12);
    for ( auto& p : gmshIO.point )
    {
      for ( int i = 0; i < 3; ++i )
      {
        fout2 << std::fixed << p[i] << " ";
      }
      fout2 << std::endl;
    }
    fout2.close();
  }

  void WriteCells()
  {
    for ( auto& c : gmshIO.cell )
    {
      for ( int i = 0; i < c.nbor.size(); ++i )
      {
        if (c.nbor[i].cell == nullptr) fout3 << -1 << " ";
        else fout3 << c.nbor[i].cell - &gmshIO.cell[0] << " ";
      }
      fout3 << "faces: ";
      for ( int i = 0; i < c.nbor.size(); ++i )
      {
        if (c.nbor[i].face == nullptr) fout3 << -1 << " ";
        else fout3 << c.nbor[i].face - &gmshIO.face[0] << " ";
      }
      fout3 << "pID: " << c.pID << " ";
      fout3 << "ctr: ";
      for ( int i = 0; i < 3; ++i )
      {
        fout3 << c.ctr[i] << " ";
      }
      fout3 << std::endl;
    }
    fout3.close();
  }
  
  void WriteFaces()
  {
    for ( auto& f : gmshIO.face )
    {
      for ( int i = 0; i < 4; ++i )
      {
        if (f.point[i] == nullptr) fout4 << -1 << " ";
        else fout4 << f.point[i] - &gmshIO.point[0] << " ";
      }
      fout4 << "pID: " << f.pID << " ";
      fout4 << "area: " << f.area << " ";
      fout4 << "norm: ";
      for ( int i = 0; i < 3; ++i )
      {
        fout4 << f.norm[i] << " ";
      }
      fout4 << std::endl;
    }
    fout4 << std::endl;
    fout4 << "BCfaces" << std::endl;
    for ( auto& f : gmshIO.faceBC )
    {
      fout4 << "face: " << f - &gmshIO.face[0] << std::endl; 
    }
    fout4.close();
  }
  
  int FileMatch(int i)
  {
    std::string lineOut, lineExpected;
    while ( std::getline( expected[i], lineExpected ) )
    {
      std::getline( in[i], lineOut );
      if ( lineOut != lineExpected )
      {
        std::cout << "Failed " << std::endl;
        return 1;
      }
    }
    if ( std::getline( in[i], lineOut ) )
    {
      std::cout << "Failed truncate " << std::endl;
      return 1;
    }
    return 0;
  }

};

}/*model*/ }/*math*/ 

int main()
{
  fout.open( "outputs/gmfaces.txt", std::ios::out );
  fout2.open( "outputs/nodes.txt", std::ios::out );
  fout3.open( "outputs/cells.txt", std::ios::out );
  fout4.open( "outputs/faces.txt", std::ios::out );
  
  in[0].open("outputs/nodes.txt", std::ios::in);
  expected[0].open("expectedOutputs/nodes.txt", std::ios::in);
  in[1].open("outputs/cells.txt", std::ios::in);
  expected[1].open("expectedOutputs/cells.txt", std::ios::in);
  in[2].open("outputs/faces.txt", std::ios::in);
  expected[2].open("expectedOutputs/faces.txt", std::ios::in);
  math::model::GmshIOTest test;
  return test.Exec();
}
