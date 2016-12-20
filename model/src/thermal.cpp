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

#include "thermal.h"
#include "matrix/seidel.h"
#include "matrix/elimination.h"
#include "matrix/mult.h"
#include "matrix/add.h"

namespace math{
namespace model{

void Thermal::Resize()
{
  unsigned long long nc = cell.size();
  ts.Resize(nc);
  Q.Resize(nc);
  Q_a.Resize(nc);
  Q_b.Resize(nc);
  A_ts.Resize(nc, nc);
}

void Thermal::Initialize()
{
  unsigned long long nc = cell.size();
  A_ts.Resize(nc, nc);
  for (unsigned long long i = 0; i < nc; ++i)
  {
    Q_a[i] = 0.0;
    Q_b[i] = 0.0;
    Q[i] = 0.0;
    ts[i] = 1.0;
  }
}

void Thermal::SolveFVM( double endTime)
{
  if (endTime < 0.0) SolveFVMSteadyState();
  else SolveFVMTransient(endTime);
}

void Thermal::SolveFVMSteadyState()
{
  Resize();
  Initialize();
  CreateMatrix();
  ApplyBC();
  solver::Seidel<SparseVect> gauss;
  gauss.A(A_ts);
  gauss.B(Q);
  gauss.X(ts);
  gauss.Tolerance(0.00001);
  gauss.MaxIter(100000);
  gauss.Solve();
  unsigned long long niter = gauss.NIter();
  int dum = 5;
}

void Thermal::SolveFVMTransient( double endTime)
{
}

void Thermal::ApplyBC()
{
  vops::Mult<double> mult;
  vops::Add<double> add;
  
  for (int j = 0; j < boundaryCondition.size(); ++j)
  {
    boundaryCondition[j]->dom = this;      
    for (unsigned long long i = 0; i < faceBC.size(); ++i)
    {
      boundaryCondition[j]->Apply(faceBC[i]);
    }
  }
  
  mult(ts, Q_a, Q);
  add(Q, Q_b, Q);
}

void Thermal::CreateMatrix()
{
  for (unsigned long long m = 0; m < cell.size(); ++m)
  {
    double k_m = cell[m].mat->cond(ts[m]);
    double coeff_m = A_ts(m, m);
    for (int j = 0; j < cell[m].nbor.size(); ++j)
    {
      if (cell[m].nbor[j].cell == nullptr ) continue;
      unsigned long long n = Index(cell[m].nbor[j].cell);
      double k_n = cell[n].mat->cond(ts[n]);
      double dtot = geom::Dist(cell[m].ctr, cell[n].ctr);
      double d1 = geom::Dist(cell[m].ctr, cell[m].nbor[j].face->ctr);
      double area = cell[m].nbor[j].face->area;
      double kfact1 = d1 / dtot;
      double kfact2 = 1.0 - kfact1;
      double kavg = kfact1 * k_m + kfact2 * k_n;
      double coeff_n = A_ts(m, n);
      coeff_n += kavg * area / dtot;
      A_ts.Set(m, n, coeff_n);
      coeff_m -= kavg * area / dtot;
    }
    A_ts.Set(m, m, coeff_m);
  }
}

}/*model*/ }/*math*/ 
