template < class Matrix, class PreconditionerInverse >
bool 
CG( int const N, 
    Matrix & A, 
    ccdl::v & x, ccdl::v const & b,
    PreconditionerInverse const & Minv, 
    int & max_iter, 
    double & tol )
{
  double resid;
  
  std::vector<double> SCR( 4*N );
  ccdl::v p( N, SCR.data() );
  ccdl::v z( N, p.end() );
  ccdl::v q( N, z.end() );
  ccdl::v r( N, q.end() );
  
  /*
  std::vector<double> sp(N),sz(N),sq(N),sr(N);
  ccdl::v p( sp );
  ccdl::v z( sz );
  ccdl::v q( sq );
  ccdl::v r( sr );
  */

  double alpha, beta, rho, rho_1 = 1.;

  double normb = b.nrm2();
  r.axpy( -1., q.dot(A,x), b );

  if ( normb == 0.0 ) 
    normb = 1.;
  
  if ( (resid = r.nrm2() / normb) <= tol ) 
    {
      tol = resid;
      max_iter = 0;
      return true;
    }

  for ( int i = 1; i <= max_iter; i++ ) 
    {

      z.dot(Minv,r);
      rho = r.dot(z);
      if ( i == 1 )
	{
	  p = z;
	}
      else 
	{
	  beta = rho / rho_1;
	  //( p *= beta ) += z;
	  p *= beta;
	  p += z;
	}
      q.dot(A,p);
      alpha = rho / p.dot(q);
      x.axpy( alpha,p);
      r.axpy(-alpha,q);

      std::cout << std::scientific << r.nrm2() / normb << "\n";
      if ( (resid = r.nrm2() / normb ) <= tol) 
	{
	  tol = resid;
	  max_iter = i;
	  return true;     
	}

      rho_1 = rho;
    }
  tol = resid;
  return false;
}


