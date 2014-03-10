


  // Assuming A is othonormal...
  // The projected part of x on the constraints
  // (A.At).x
  // The part of x after removing the constraints
  // ( 1 - A.At ).x
  // If we add the constraints back in
  // ( 1 - A.At ).x + A.b
  // x + A.(b-At.x)
  //
  // If A is not orthonormal,
  // The projected part of x on the constraints
  // A.(At.A)^{-1}.At.x
  // The part of x after removing the constraints
  // ( 1 - A.(At.A)^{-1}.At ).x
  // If we add the constraints back in
  // ( 1 - A.(At.A)^{-1}.At ).x + A.(At.A)^{-1}.b
  // x + A.(At.A)^{-1}.( b - At.x )

// P = A.(At.A)^{-1} is called the "constraint_projector" of A
// so that
// x += P.(b-At.x)
// is an updated x that now enforces the constraints
// It can be computed by
// (P=A).constraint_projector();
//
// The update of x can be computed from
// x.enforce_constraints( A, P, b );


inline int ccdl::ge::query_constraint_projector( int const M, int const Nconstraints )
{
  return Nconstraints*Nconstraints + std::max( M*Nconstraints , ccdl::ge::query_svd_inverse(Nconstraints,Nconstraints) );
}

inline ccdl::ge & ccdl::ge::constraint_projector( int const nscr, double * scr )
{
  assert( nscr >= ccdl::ge::query_constraint_projector( nfast, nslow ) );
  ccdl::sy AtAinv( nslow, scr );
  AtAinv.ge().dot(this->t(),*this).svd_inverse( 1.e-10, nscr - std::distance(scr,AtAinv.end()), AtAinv.end() );
  dot( ccdl::ge(nfast,nslow,AtAinv.end()) = *this , AtAinv );
  return *this;
}


inline ccdl::ge & ccdl::ge::constraint_projector()
{
  std::vector<double> scr( ccdl::ge::query_constraint_projector( nfast, nslow ) );
  return constraint_projector( scr.size(), scr.data() );
}


inline int ccdl::v::query_enforce_constraints( int const nconstraints ) { return nconstraints; }



inline ccdl::v & ccdl::v::enforce_constraints( ccdl::ge & A, ccdl::ge const & constraint_projector_of_A, ccdl::v const & constraint_values, int const nscr, double * scr )
{
  assert( A.nfast == constraint_projector_of_A.nfast );
  assert( A.nslow == constraint_projector_of_A.nslow );
  assert( nfast == A.nfast );
  assert( A.nslow == constraint_values.nfast );
  assert( nscr >= A.nslow );
  ccdl::v c_minus_Ax(A.nslow,scr);
  (c_minus_Ax=constraint_values).dot(-1.,A.t(),*this,1.);
  dot( 1., constraint_projector_of_A,c_minus_Ax, 1. );
  return *this;
}


inline ccdl::v & ccdl::v::enforce_constraints( ccdl::ge & A, ccdl::ge const & constraint_projector_of_A, ccdl::v const & constraint_values )
{
  std::vector<double> scr( ccdl::v::query_enforce_constraints( A.nslow ) );
  return enforce_constraints(A,constraint_projector_of_A,constraint_values,scr.size(),scr.data());
}


inline ccdl::v & ccdl::v::enforce_constraints( ccdl::ge & A, ccdl::v const & constraint_values )
{
  int n_enforce = ccdl::v::query_enforce_constraints( A.nslow );
  int n_project = A.nfast*A.nslow + ccdl::ge::query_constraint_projector( A.nfast, A.nslow );
  std::vector<double> scr( std::max(n_enforce,n_project) );
  ccdl::ge P(A.nfast,A.nslow,scr.data());
  int nscr = scr.size() - A.nfast*A.nslow;
  (P=A).constraint_projector( nscr, P.end() );
  return enforce_constraints(A,P,constraint_values,nscr,P.end());
}


