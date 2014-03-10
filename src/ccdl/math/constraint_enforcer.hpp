#ifndef _ccdl_constraint_enforcer_hpp_
#define _ccdl_constraint_enforcer_hpp_

namespace ccdl
{

  class constraint_enforcer
  {
  public:
    constraint_enforcer( ccdl::ge const conMat, ccdl::v const conVals );
    ccdl::v & enforce_constraints( ccdl::v & x );

    int nfast, nslow;
    std::vector<double> constraint_matrix;
    std::vector<double> enforcement_matrix;
    std::vector<double> constraint_values;
  };

}

inline ccdl::constraint_enforcer::constraint_enforcer
( ccdl::ge const conMat, ccdl::v const conVals )
  : nfast(conMat.nfast()),
    nslow(conMat.nslow()),
    constraint_matrix(conMat.begin(),conMat.end()),
    enforcement_matrix(nfast*nslow),
    constraint_values(conVals.begin(),conVals.end())
{
  assert( nslow < nfast );
  assert( conVals.size() == conMat.nslow() );

  // scratch space
  std::vector<double> scratch(nslow*nslow + ccdl::sy::query_svd_inverse(nslow));

  // wrappers
  ccdl::ge A(nfast,nslow,constraint_matrix.data());
  ccdl::ge A_dot_AtAinv(nfast,nslow,enforcement_matrix.data());
  ccdl::sy AtA(nslow,scratch.data());

  // A . (At.A)^{-1}
  AtA.dot(A.t(),A).svd_inverse( 1.e-10, scratch.size()-nslow*nslow, AtA.end() );
  A_dot_AtAinv.dot( A, AtA );
}



inline ccdl::v & ccdl::constraint_enforcer::enforce_constraints
( ccdl::v & x )
{
  assert( x.size() == nfast );

  // scratch space
  // notice that scratch is initialized to constraint values
  std::vector<double> scratch(constraint_values);

  // wrappers
  ccdl::ge A(nfast,nslow,constraint_matrix.data());
  ccdl::ge A_dot_AtAinv(nfast,nslow,enforcement_matrix.data());
  ccdl::v c_minus_Ax(nslow,scratch.data());

  // c - At.x
  c_minus_Ax.dot(-1.,A.t(),x,1.);

  // x := x + A . (At.A)^{-1} . ( c - At.x )
  return x.dot( 1., A_dot_AtAinv, c_minus_Ax, 1. );
}



#endif
