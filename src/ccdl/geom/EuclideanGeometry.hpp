#ifndef _ccdl_EuclideanGeometry_hpp_
#define _ccdl_EuclideanGeometry_hpp_

#include <cmath>

namespace ccdl
{

  template <class DblIn,class DblOut>
  inline void CrossProduct( DblIn const * a, DblIn const * b, DblOut * axb );

  template <class DblIn,class DblOut>
  inline void AccumulateCrossProduct( DblIn const * a, DblIn const * b, DblOut * axb );

  template <class DblIn,class DblOut>
  inline void OrthogonalVector( DblIn const * a, DblIn const * b, DblIn const * c, DblOut * n );

  double CosineAngle( double const * ca, 
		      double const * cb, 
		      double const * cc );

  double CosineAngle( double const * ca, 
		      double const * cb, 
		      double const * cc,
		      double * ga, 
		      double * gb, 
		      double * gc );

  double Bond( double const * ca, 
	       double const * cb );

  double Bond( double const * ca, 
	       double const * cb, 
	       double * ga, 
	       double * gb );

  double Angle( double const * ca, 
		double const * cb, 
		double const * cc );

  double Angle( double const * ca, 
		double const * cb, 
		double const * cc,
		double * ga, 
		double * gb, 
		double * gc );

  double DihedralAngle( double const * ca, 
			double const * cb, 
			double const * cc, 
			double const * cd, 
			bool & linear );
  
  double DihedralAngle( double const * ca, 
			double const * cb, 
			double const * cc, 
			double const * cd, 
			bool & linear,
			double * ga, 
			double * gb, 
			double * gc, 
			double * gd );

}


template <class DblIn,class DblOut>
inline void ccdl::CrossProduct
( DblIn const * a, 
  DblIn const * b,
  DblOut * axb )
{
  // http://en.wikipedia.org/wiki/Cross_product#Coordinate_notation
  axb[0] = a[1]*b[2] - a[2]*b[1];
  axb[1] = b[0]*a[2] - a[0]*b[2];
  axb[2] = a[0]*b[1] - a[1]*b[0];
}

template <class DblIn,class DblOut>
inline void ccdl::AccumulateCrossProduct
( DblIn const * a, 
  DblIn const * b,
  DblOut * axb )
{
  // http://en.wikipedia.org/wiki/Cross_product#Coordinate_notation
  axb[0] += a[1]*b[2] - a[2]*b[1];
  axb[1] += b[0]*a[2] - a[0]*b[2];
  axb[2] += a[0]*b[1] - a[1]*b[0];
}


template <class DblIn,class DblOut>
inline void ccdl::OrthogonalVector
( DblIn const * a, 
  DblIn const * b, 
  DblIn const * c, 
  DblOut * n )
{
  // http://en.wikipedia.org/wiki/Dihedral_angle#Methods_of_computation

  DblIn ba[] = { b[0]-a[0],b[1]-a[1],b[2]-a[2] };
  DblIn ca[] = { c[0]-a[0],c[1]-a[1],c[2]-a[2] };
  ccdl::CrossProduct(ba,ca,n);
}


#endif
