#ifndef _CCDL_GAUSSIANQUADRATUREMOD_HPP_
#define _CCDL_GAUSSIANQUADRATUREMOD_HPP_


namespace YGXX
{
  
  /** \brief Gauss-Laguerre quadrature points and weights
   *
   * Calculates the quadrature points and weights for Gauss-Laguerre
   *     integration with the parameter alpha of the Laguerre polynomials equal
   *     to alf.  The output arrays x (the coordinates of the quadrature 
   *     points) and w (the associated weights for each point) must be 
   *     allocated and of equal size.  
   *     The size of the output arrays determines the number of quadrature 
   *     points that will be calculated.
   *     See Numerical Recipes a detailed discussion of Gauss-Laguerre 
   *     integration.
   *
   * The weight function is \f$w(x)=x^{\alpha}\mbox{e}^{-x}\f$
   *
   * @param[in] rAlf = alpha exponent
   * @param[in] n = quadrature rule
   * @param[out] rVecX = vector of quadrature points
   * @param[out] rVecW = vector of quadrature weights
   *
   * \note The quadrature weights will need to have the weight function removed
   */
  void GaussLaguerreRule
  ( double const rAlf,
    int const n,
    double * rVecX, 
    double * rVecW );



  /** \brief Gauss-Legendre quadrature points and weights
   *
   *     Calculates the quadrature points and weights for Gauss-Legendre
   *     integration over the interval x1..x2.
   *     See Numerical Recipes a detailed discussion of Gauss-Legendre 
   *     integration.
   *
   * The weight function is \f$w(x)=1\f$
   *
   * @param[in] x1 = lower range
   * @param[in] x2 = upper range
   * @param[in] n = quadrature rule
   * @param[out] rVecX = vector of quadrature points
   * @param[out] rVecW = vector of quadrature weights
   *
   */
  void GaussLegendreRule
  ( double const x1,
    double const x2,
    int const n,
    double * rVecX, 
    double * rVecW );




  /** \brief Gauss-Hermite quadrature points and weights
   *
   *     Calculates the quadrature points and weights for Gauss-Hermite
   *     integration. 
   *     The order of the quadrature points returned in x and w are in
   *     DESENDING order. 
   *     The size of the output arrays determines the number of quadrature 
   *     points that will be calculated.
   *     See Numerical Recipes a detailed discussion of Gauss-Hermite 
   *     integration.
   *
   * The weight function is \f$w(x)=\mbox{e}^{-x^2}\f$
   *
   * @param[in] n = quadrature rule
   * @param[out] rVecX = vector of quadrature points
   * @param[out] rVecW = vector of quadrature weights
   *
   * \note The quadrature weights will need to have the weight function removed
   */
  void GaussHermiteRule
  ( int const n,
    double * rVecX, 
    double * rVecW );




  /** \brief Gauss-Jacobi quadrature points and weights
   *
   *     Calculates the quadrature points and weights for Gauss-Jacobi
   *     integration with the parameters alpha and beta of the Jacobi
   *     polynomials equal to alf and bet respectively. 
   *     See Numerical Recipes a detailed discussion of Gauss-Jacobi 
   *     integration.
   *
   * The weight function is \f$w(x)=(1-x)^{\alpha}(1+x)^{\beta}\f$
   *
   * @param[in] alf = alpha weight parameter
   * @param[in] bet = beta weight parameter
   * @param[in] n = quadrature rule
   * @param[out] rVecX = vector of quadrature points
   * @param[out] rVecW = vector of quadrature weights
   *
   * \note The quadrature weights will need to have the weight function removed
   */
  void GaussJacobiRule
  ( double const alf,
    double const bet,
    int const n,
    double * rVecX, 
    double * rVecW );


  /** \brief Lebedev angular quadrature rule
   *
   * Returns the Lebedev quadrature points and weights for a particular
   * integer rule.
   *
   * The valid rules are:
   *    6,   14,   26,   38,   50, 74
   *   86,  110,  146,  170,  194, 230
   *  266,  302,  350,  434,  590,
   *  770,  974, 1202, 1454, 1730,
   * 2030, 2354, 2702, 3074, 3470,
   * 3890, 4334, 4802, 5294, 5810
   *
   *
   * However, the following rules are NOT recommended: 74, 230, and 266
   *
   *
   * The following rules have special transfer relationships with regard
   * to spherical harmonics: 14, 38, 86, 110, 230, 302, 350, 590.
   *
   *
   *    Users of the Lebedev angular quadrature rules should cite:
   *
   *    V.I. Lebedev, and D.N. Laikov, A quadrature formula for the 
   *       sphere of the 131st algebraic order of accuracy,
   *       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   *
   *    LebedevGrid code based partially on code from 
   *    http://server.ccl.net/cca/software/SOURCES/FORTRAN/Lebedev-Laikov-Grids/Lebedev-Laikov.F
   *
   *   Lebedev Quadrature references:
   *
   *   (1) V.I. Lebedev
   *         Values of the nodes and weights of ninth to seventeenth 
   *         order Gauss-Markov quadrature formulae invariant under the
   *         octahedron group with inversion
   *         Computational Mathematics and Mathematical Physics, Vol. 15,
   *         1975, pp. 44-51.
   *
   *   (2) V.I. Lebedev
   *         Quadratures on a sphere
   *         Computational Mathematics and Mathematical Physics, Vol. 16,
   *         1976, pp. 10-24. 
   *
   *   (3) V.I. Lebedev
   *         Spherical quadrature formulas exact to orders 25-29
   *         Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107.
   *
   *   (4) V.I. Lebedev, and A.L. Skorokhodov
   *         Quadrature formulas of orders 41, 47, and 53 for the sphere
   *         Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   *
   *   (5) V.I. Lebedev
   *       A quadrature formula for the sphere of 59th algebraic
   *         order of accuracy
   *         Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   *
   *   (6) V.I. Lebedev, and D.N. Laikov
   *         A quadrature formula for the sphere of the 131st
   *         algebraic order of accuracy
   *         Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   *
   *   Other useful references:
   *
   *   (1) A. D. Becke, J. Chem. Phys. 88, 2547, (1988).
   *
   *   (2) B. Delley 
   *         High Order Integration Schemes on the Unit Sphere,
   *         Journal of Computational Chemistry 1996, vol. 17, pp. 1152-1155.
   *
   *   (3) A. H. Stroud, Approximate Calculations of Multiple Integrals, 
   *        Prentice-Hall, Englewood Cliffs (1971).
   *
   * SPECIAL NOTE
   *
   *   Lebedev-Laikov grids of order n=2m+1, where m=1,2,...,15, and additionally
   *   grids of order n=6m+5, where m=5,6,...,21 are implemented in this code.
   *   Three of the Lebedev angular quadrature rules (74, 230 and 266 pts) have points with 
   *   negative weights. While this is not necessisarily a bad thing, it causes COSMO some
   *   problems. I asked Dr. Laikov about whether there were other formulations of these three
   *   rules, and he responded:
   *
   *        As you know we have solved the equations for quadrature formulas
   *        of orders L=6m+5, m=0,...,21 which have very regular point distributions,
   *        but for orders not equal to 6m+5 the point distributions are quite
   *        unpredictable and negative weights are found in some cases - it is hard to solve for 
   *        high orders not equal to 6m+5. Moreover, there are multiple solutions in some cases, 
   *        either with all weights positive or having some negative weights. I have included in 
   *        my compilation only those solutions which have the smallest number of negative weights...
   *        ... I would recommend you to use only grids for L=6m+5. They all have beautiful point 
   *        distributions and all weights positive.
   *
   *
   *  @param[in] iRule = must correspond to a valid integer rule
   *  @param[out] rVecPt = array of points on a unit sphere
   *  @param[out] rVecWt = array of quadrature weights
   *
   */
  void LebedevRule
  ( unsigned int const & iRule,
    std::tr1::array<double,3> * rVecPt,
    double * rVecWt );


  double LebedevCosmoGaussianZetaScaleFactor( int const iRule );

}


#endif
