//<license>

//</license>

#include "table.h"

namespace math{

Table::Table( double* d )
{
  m_d = d;
}

void Table::Data( double* d )
{
  m_d = d;
}

Uint Table::Size()
{
  Uint res = 1;
  Uint sum = 1;
  for ( Uint i = 1; i <= static_cast<Uint>(m_d[0]); ++i )
  {
    res *= static_cast<Uint>(m_d[i]);
    sum += static_cast<Uint>(m_d[i]);
  }
  return res + sum + static_cast<Uint>(m_d[0]);
}

double Table::operator()( double d1 )
{
  //m_d[1] = number of independent values 1
  //m_d[2] = start of independent values
  //m_d[2+ static_cast<Uint>(m_d[1])] = start of dependent values
  //m_d[ 1 + 2 * static_cast<Uint>(m_d[1]) ] = end of dependent values
  return Interp1( &m_d[1], &m_d[2], &m_d[2+ static_cast<Uint>(m_d[1])], d1 );
}

double Table::Interp1(  double size[], double valI1[], double valD[], double d1 )
{
  Sint idx = UpperBoundIdx( valI1, valI1 + static_cast<Uint>(size[0]), d1 );
  if ( idx <= 0 )
  {
    return valD[-idx];
  }
  double yr = ( d1 - valI1[ idx - 1 ] ) / ( valI1[ idx ] - valI1[ idx - 1 ] );
  return yr * ( valD[idx] - valD[idx - 1] ) + valD[idx - 1];
}

double Table::operator()( double d1, double d2 )
{
  //m_d[1] = number of independent values 1
  //m_d[2] = number of independent values 2
  //m_d[3] = start of independent values 1
  //m_d[3 + static_cast<Uint>(m_d[1]) ] = start of independent values 2
  //m_d[3 + static_cast<Uint>(m_d[1]) + static_cast<Uint>(m_d[2]) ] = start of dependent values
  //m_d[2 + static_cast<Uint>(m_d[1]) + static_cast<Uint>(m_d[2]) + static_cast<Uint>(m_d[1]) * static_cast<Uint>(m_d[2]) ] = end of dependent values
  return Interp2( &m_d[1], &m_d[3], &m_d[3 + static_cast<Uint>(m_d[1]) ], &m_d[3 + static_cast<Uint>(m_d[1]) + static_cast<Uint>(m_d[2]) ], d1, d2 );
}

double Table::Interp2( double size[], double valI1[], double valI2[], double valD[], double d1, double d2 )
{
  Sint idx[2];
  
  idx[0] = UpperBoundIdx( valI1, valI1 + static_cast<Uint>(size[0]), d1 );
  if ( idx[0] <= 0 )
  {
    return Interp1( &size[1], valI2, &valD[-idx[0] * static_cast<Uint>(size[1]) ], d2 );
  }
  
  double yr[2];
  
  yr[0] = ( d1 - valI1[ idx[0] - 1 ] ) / ( valI1[ idx[0] ] - valI1[ idx[0] - 1 ] );
  idx[1] = UpperBoundIdx( valI2, valI2 + static_cast<Uint>(size[1]), d2 );
  if ( idx[1] <= 0 )
  {
    return yr[0] * ( valD[ idx[0] * static_cast<Uint>(size[1]) - idx[1] ] - valD[ ( idx[0] - 1 ) * static_cast<Uint>(size[1]) - idx[1] ] ) + valD[ ( idx[0] - 1 ) * static_cast<Uint>(size[1]) - idx[1] ];
  }
  yr[1] = ( d2 - valI2[ idx[1] - 1 ] ) / ( valI2[ idx[1] ] - valI2[ idx[1] - 1 ] );
  
  double r[2];
  r[0] = yr[1] * ( valD[ ( idx[0] - 1 ) * static_cast<Uint>(size[1]) + idx[1] ] - valD[ ( idx[0] - 1 ) * static_cast<Uint>(size[1]) + idx[1] - 1 ] ) + valD[ ( idx[0] - 1 ) * static_cast<Uint>(size[1]) + idx[1] - 1 ];
  r[1] = yr[1] * ( valD[ idx[0] * static_cast<Uint>(size[1]) + idx[1] ] - valD[ idx[0] * static_cast<Uint>(size[1]) + idx[1] - 1 ] ) + valD[ idx[0] * static_cast<Uint>(size[1]) + idx[1] - 1 ];
  return yr[0] * ( r[1] - r[0] ) + r[0];
}

double Table::operator()( double d1, double d2, double d3 )
{
  //m_d[1] = number of independent values 1
  //m_d[2] = number of independent values 2
  //m_d[3] = number of independent values 3
  //m_d[4] = start of independent values 1
  //m_d[4 + static_cast<Uint>(m_d[1]) ] = start of independent values 2
  //m_d[4 + static_cast<Uint>(m_d[1]) + static_cast<Uint>(m_d[2]) ] = start of independent values 3
  //m_d[4 + static_cast<Uint>(m_d[1]) + static_cast<Uint>(m_d[2]) + static_cast<Uint>(m_d[3]) ] = start of dependent values
  //m_d[3 + static_cast<Uint>(m_d[1]) + static_cast<Uint>(m_d[2]) + static_cast<Uint>(m_d[3]) + static_cast<Uint>(m_d[1]) * static_cast<Uint>(m_d[2]) * static_cast<Uint>(m_d[3]) ] = end of dependent values
  return Interp3( &m_d[1], &m_d[4], &m_d[4 + static_cast<Uint>(m_d[1]) ], &m_d[4 + static_cast<Uint>(m_d[1]) + static_cast<Uint>(m_d[2]) ], &m_d[4 + static_cast<Uint>(m_d[1]) + static_cast<Uint>(m_d[2]) + static_cast<Uint>(m_d[3]) ], d1, d2, d3 );
}

double Table::Interp3( double size[], double valI1[], double valI2[], double valI3[], double valD[], double d1, double d2, double d3 )
{
  Sint idx[3];
  
  idx[0] = UpperBoundIdx( valI1, valI1 + static_cast<Uint>(size[0]), d1 );
  if ( idx[0] <= 0 )
  {
    return Interp2( &size[1], valI2, valI3, &valD[-idx[0] * static_cast<Uint>(size[1]) * static_cast<Uint>(size[2]) ], d2, d3 );
  }
  
  double yr[3];
  
  yr[0] = ( d1 - valI1[ idx[0] - 1 ] ) / ( valI1[ idx[0] ] - valI1[ idx[0] - 1 ] );
  idx[1] = UpperBoundIdx( valI2, valI2 + static_cast<Uint>(size[1]), d2 );
  double r[8];
  if ( idx[1] <= 0 )
  {
    r[0] = Interp1( &size[2], valI3, &valD[ (idx[0] - 1) * static_cast<Uint>(size[1]) * static_cast<Uint>(size[2]) - idx[1] * static_cast<Uint>(size[2]) ], d3 );
    r[1] = Interp1( &size[2], valI3, &valD[ idx[0] * static_cast<Uint>(size[1]) * static_cast<Uint>(size[2]) - idx[1] * static_cast<Uint>(size[2]) ], d3 );
    return yr[0] * ( r[1] - r[0] ) + r[0];
  }
  yr[1] = ( d2 - valI2[ idx[1] - 1 ] ) / ( valI2[ idx[1] ] - valI2[ idx[1] - 1 ] );
  idx[2] = UpperBoundIdx( valI3, valI3 + static_cast<Uint>(size[2]), d3 );
  if ( idx[2] <= 0 )
  {
    r[0] = valD[ (idx[0] - 1) * static_cast<Uint>(size[1]) * static_cast<Uint>(size[2]) + (idx[1] - 1) * static_cast<Uint>(size[2]) - (idx[2] ) ];
    r[1] = valD[ (idx[0] - 1) * static_cast<Uint>(size[1]) * static_cast<Uint>(size[2]) + (idx[1] ) * static_cast<Uint>(size[2]) - (idx[2] ) ];
    r[2] = valD[ (idx[0] ) * static_cast<Uint>(size[1]) * static_cast<Uint>(size[2]) + (idx[1] - 1) * static_cast<Uint>(size[2]) - (idx[2] ) ];
    r[3] = valD[ (idx[0] ) * static_cast<Uint>(size[1]) * static_cast<Uint>(size[2]) + (idx[1] ) * static_cast<Uint>(size[2]) - (idx[2] ) ];
    
    r[0] = yr[1] * ( r[1] - r[0] ) + r[0];
    r[2] = yr[1] * ( r[3] - r[2] ) + r[2];
    
    return yr[0] * ( r[2] - r[0] ) + r[0];
  }
  yr[2] = ( d3 - valI3[ idx[2] - 1 ] ) / ( valI3[ idx[2] ] - valI3[ idx[2] - 1 ] );
  r[0] = valD[ (idx[0] - 1) * static_cast<Uint>(size[1]) * static_cast<Uint>(size[2]) + (idx[1] - 1) * static_cast<Uint>(size[2]) + (idx[2] - 1) ];
  r[1] = valD[ (idx[0] - 1) * static_cast<Uint>(size[1]) * static_cast<Uint>(size[2]) + (idx[1] - 1) * static_cast<Uint>(size[2]) + (idx[2] ) ];
  r[2] = valD[ (idx[0] - 1) * static_cast<Uint>(size[1]) * static_cast<Uint>(size[2]) + (idx[1] ) * static_cast<Uint>(size[2]) + (idx[2] - 1) ];
  r[3] = valD[ (idx[0] - 1) * static_cast<Uint>(size[1]) * static_cast<Uint>(size[2]) + (idx[1] ) * static_cast<Uint>(size[2]) + (idx[2] ) ];
  r[4] = valD[ (idx[0] ) * static_cast<Uint>(size[1]) * static_cast<Uint>(size[2]) + (idx[1] - 1) * static_cast<Uint>(size[2]) + (idx[2] - 1) ];
  r[5] = valD[ (idx[0] ) * static_cast<Uint>(size[1]) * static_cast<Uint>(size[2]) + (idx[1] - 1) * static_cast<Uint>(size[2]) + (idx[2] ) ];
  r[6] = valD[ (idx[0] ) * static_cast<Uint>(size[1]) * static_cast<Uint>(size[2]) + (idx[1] ) * static_cast<Uint>(size[2]) + (idx[2] - 1) ];
  r[7] = valD[ (idx[0] ) * static_cast<Uint>(size[1]) * static_cast<Uint>(size[2]) + (idx[1] ) * static_cast<Uint>(size[2]) + (idx[2] ) ];
  
  for ( int i = 0; i < 7; i += 2)
  {
    r[i] = yr[2] * ( r[i+1] - r[i] ) + r[i];
  }
  for ( int i = 0; i < 7; i += 4)
  {
    r[i] = yr[1] * ( r[i+2] - r[i] ) + r[i];
  }
  return yr[2] * ( r[4] - r[0] ) + r[0];
}

Sint Table::UpperBoundIdx( double* beg, double* end, const double val )
{
  if ( val < *beg )
  {
    return 0;
  }
  --end;
  if ( val >= *end )
  {
    return -( end - beg );
  }
  double* zeroPos = beg;
  while ( ( end - beg ) > 2 )
  {
    double* mid = &beg[ ( end - beg ) >> 1 ]; // >> 1 is division by 2
    if ( *mid > val )
    {
      end = mid;
    }
    else 
    {
      beg = mid;
    }
  }
  Uint dist = end - beg;
  for ( Uint i = 0; i <= dist; ++i )
  {
    if ( *beg >= val )
    {
      if ( (*beg - val) / ( *beg - *(beg - 1) ) > 0.999 )
      {
        return -(beg - zeroPos - 1);
      }
      else if ( (*beg - val) / ( *beg - *(beg - 1) ) < 0.0001 )
      {
        return -(beg - zeroPos);
      }
      return (beg - zeroPos);
    }
    ++beg;    
  }
}

}/*math*/ 
