#ifndef _ccdl_NumericalBasis_hpp_
#define _ccdl_NumericalBasis_hpp_

namespace ccdl
{


  /** \overload */
   void HermiteBasis
   ( double const zetaExp,
     int const polyOrder,
     int const npt,
     double const * xVec,
     double *  basisVec );
  


  /** \overload */
   void HermiteBasis
   ( double const zetaExp,
     int const polyOrder,
     int const npt,
     double const *  xVec,
     double *  basisVec,
     double *  basisDerivVec );



  /** \overload */
   void HermiteBasis
   ( double const zetaExp,
     int const polyOrder,
     int const npt,
     double const *  xVec,
     double *  basisVec,
     double *  basisDerivVec,
     double *  basisSecDerivVec);
  
   




  /** \brief L2-normalized, 1-dimensional, Gaussian-decaying Hermite polynomial basis
   *
   *
   * This basis would typically be used to solve the 1-dimensional Schodinger equation
   * via the Rayleigh-Ritz variational method.
   *
   * This basis \f$B_n(x;\zeta)\f$ has the form
   * \f[
   *   B_n(x;\zeta) = N H_n(\zeta x) \mbox{e}^{-(\zeta x)^2/2}
   * \f]
   * where
   * \f[
   *   N = \sqrt{\frac{z}{\sqrt{\pi} n! 2^n}}
   * \f]
   * is a normalization constant such that
   *
   * \f[
   *   \int_{-\infty}^{\infty} B_m(x;\zeta) B_n(x;\zeta) dx =
   *       \delta_{m,n}
   * \f]
   *
   *  @param[in] zetaExp = \f$ \zeta \f$
   *  @param[in] polyOrder = order, i.e., n in the above equations
   *  @param[in] npt = number of x-values
   *  @param[in] xVec = an array of x-values
   *  @param[out] basisVec = the values of \f$ B_n(x;\zeta) \f$ at each x
   *  @param[out] basisDerivVec = the values of \f$ d/dx B_n(x;\zeta) \f$ at each x
   *  @param[out] basisSecDerivVec = the values of \f$ d^2/dx^2 B_n(x;\zeta) \f$ at each x
   *  @param[out] basisThirdDerivVec = the values of \f$ d^3/dx^3 B_n(x;\zeta) \f$ at each x
   *
   * \note For a perfect simple harmonic oscillator, this function returns the
   * exact wavefunction when \f$\zeta = \sqrt{ \frac{\sqrt{\mu k}}{\hbar} }\f$, where
   *  \f$\mu\f$ is the reduced mass (atomic units), k is the force constant (atomic units), i.e.,
   * \f$k = d^2/dx^2 V(x) \f$, and \f$\hbar\f$ is the reduced Planck's constant (which is
   * unity in atomic units).  Recall that a simple harmonic oscillator has the potential
   * \f$V(x) = kx^2/2\f$
   *
   */
   void HermiteBasis
   ( double const zetaExp,
     int const polyOrder,
     int const npt,
     double const *  xVec,
     double *  basisVec,
     double *  basisDerivVec,
     double *  basisSecDerivVec,
     double *  basisThirdDerivVec);
  
























  /** \overload */
   void ClmLaguerreBasis
   ( double const zetaExp,
     int const N,
     int const l,
     int const npt,
     double const *  xVec,
     double *  basisVec );
  


  /** \overload */
   void ClmLaguerreBasis
   ( double const zetaExp,
     int const N,
     int const l,
     int const npt,
     double const *  xVec,
     double *  basisVec,
     double *  basisDerivVec );
  


  /** \overload */
   void ClmLaguerreBasis
   ( double const zetaExp,
     int const N,
     int const l,
     int const npt,
     double const *  xVec,
     double *  basisVec,
     double *  basisDerivVec,
     double *  basisSecDerivVec);
  
   







  /** \brief Exponentially decaying Laguerre polynomial radial basis
   *
   *
   * This basis \f$B_n(r;n,l,\zeta)\f$ has the form
   * \f[
   *   B(r;n,l,\zeta) = N L_n^{(2l+2)}(2\zeta r) \mbox{e}^{-\zeta r}
   * \f]
   * where
   * \f[
   *   N = (2\zeta)^l \sqrt{ \frac{ (2\zeta)^3 n! }{ (n+2l+2)! } } \sqrt{ \frac{ 2l+1 }{ 4\pi } }
   * \f]
   * is a normalization constant such that
   *
   * \f[
   *   \int B(r;n,l,\zeta) C_{l,m}(\mathbf{r}) B(r;n',l',\zeta) C_{l',m'}(\mathbf{r}) d^3r =
   *       \delta_{n,n'} \delta_{l,l'} \delta_{m,m'}
   * \f]
   *
   * and
   *
   * \f[
   *   \int_{0}^{\infty} B(r;n,l,\zeta) B(r;n',l',\zeta) r^2 dr =
   *       \delta_{n,n'} \delta_{l,l'} \frac{ 2l+1 }{ 4\pi }
   * \f]
   *
   *  @param[in] zetaExp = \f$ \zeta \f$
   *  @param[in] N = the order of polynomial, n in the above notation.
   *  @param[in] l = principle angular momentum
   *  @param[in] npt = number of x-values
   *  @param[in] xVec = an array of x-values
   *  @param[out] basisVec = the values of \f$ B(r;n,l,\zeta) \f$ at each x
   *  @param[out] basisDerivVec = the values of \f$ d/dx B \f$ at each x
   *  @param[out] basisSecDerivVec = the values of \f$ d^2/dx^2 B \f$ at each x
   *  @param[out] basisThirdDerivVec = the values of \f$ d^3/dx^3 B \f$ at each x
   *
   */
   void ClmLaguerreBasis
   ( double const zetaExp,
     int const N,
     int const l,
     int const npt,
     double const *  xVec,
     double *  basisVec,
     double *  basisDerivVec,
     double *  basisSecDerivVec,
     double *  basisThirdDerivVec);
  



  /** \overload */
  void YlmLaguerreBasis
  ( double const zetaExp,
    int const N,
    int const l,
    int const npt,
    double const *  xVec,
    double *  basisVec );
  
  
  
  /** \overload */
  void YlmLaguerreBasis
  ( double const zetaExp,
    int const N,
    int const l,
    int const npt,
    double const *  xVec,
    double *  basisVec,
    double *  basisDerivVec );
  
  
  /** \overload */
  void YlmLaguerreBasis
  ( double const zetaExp,
    int const N,
    int const l,
    int const npt,
    double const *  xVec,
    double *  basisVec,
    double *  basisDerivVec,
    double *  basisSecDerivVec);
  
   


  /** \brief Exponentially decaying Laguerre polynomial radial basis
   *
   *
   * This basis \f$B_n(r;n,l,\zeta)\f$ has the form
   * \f[
   *   B(r;n,l,\zeta) = N L_n^{(2l+2)}(2\zeta r) \mbox{e}^{-\zeta r} r^l
   * \f]
   * where
   * \f[
   *   N = (2\zeta)^l \sqrt{ \frac{ (2\zeta)^3 n! }{ (n+2l+2)! } }
   * \f]
   * is a normalization constant such that
   *
   * \f[
   *   \int B(r;n,l,\zeta) Y_{l,m}(\Omega) B(r;n',l',\zeta) Y_{l',m'}(\Omega) d^3r =
   *       \delta_{n,n'} \delta_{l,l'} \delta_{m,m'}
   * \f]
   *
   * and
   *
   * \f[
   *   \int_{0}^{\infty} B(r;n,l,\zeta) B(r;n',l',\zeta) r^2 dr =
   *       \delta_{n,n'} \delta_{l,l'} 
   * \f]
   *
   *  @param[in] zetaExp = \f$ \zeta \f$
   *  @param[in] N = the order of polynomial, n in the above notation.
   *  @param[in] l = principle angular momentum
   *  @param[in] xVec = an array of x-values
   *  @param[out] basisVec = the values of \f$ B(r;n,l,\zeta) \f$ at each x
   *  @param[out] basisDerivVec = the values of \f$ d/dx B \f$ at each x
   *  @param[out] basisSecDerivVec = the values of \f$ d^2/dx^2 B \f$ at each x
   *  @param[out] basisThirdDerivVec = the values of \f$ d^3/dx^3 B \f$ at each x
   *
   */
  void YlmLaguerreBasis
  ( double const zetaExp,
    int const N,
    int const l,
    int const npt,
    double const * xVec,
    double *  basisVec,
    double *  basisDerivVec,
    double *  basisSecDerivVec,
    double *  basisThirdDerivVec );
  







  /** \overload */
  void ClmGaussianBasis
  ( double const zetaExp,
    int const l,
    int const npt,
    double const *  xVec,
    double *  basisVec );
  


  /** \overload */
  void ClmGaussianBasis
  ( double const zetaExp,
    int const l,
    int const npt,
    double const * xVec,
    double *  basisVec,
    double *  basisDerivVec );
  
  
  
  /** \overload */
  void ClmGaussianBasis
  ( double const zetaExp,
    int const l,
    int const npt,
    double const * xVec,
    double *  basisVec,
    double *  basisDerivVec,
    double *  basisSecDerivVec);
  
   


  /** \brief Gaussian radial basis
   *
   *
   * This basis \f$B_n(r;l,\zeta)\f$ has the form
   * \f[
   *   B(r;l,\zeta) = N \mbox{e}^{-\zeta r^2}
   * \f]
   * where
   * \f[
   *   N = \sqrt{ 8 \frac{ (4\zeta)^l }{ (2l+1)!! } \sqrt{ \frac{ (2\zeta)^l }{ \pi } } } \sqrt{ \frac{ 2l+1 }{ 4\pi } }
   * \f]
   * is a normalization constant such that
   *
   * \f[
   *   \int B(r;l,\zeta) C_{l,m}(\mathbf{r}) B(r;l',\zeta) C_{l',m'}(\mathbf{r}) d^3r =
   *       \delta_{l,l'} \delta_{m,m'}
   * \f]
   *
   *
   *  @param[in] zetaExp = \f$ \zeta \f$
   *  @param[in] l = principle angular momentum
   *  @param[in] npt = number of x-values
   *  @param[in] xVec = an array of x-values
   *  @param[out] basisVec = the values of \f$ B(r;n,l,\zeta) \f$ at each x
   *  @param[out] basisDerivVec = the values of \f$ d/dx B \f$ at each x
   *  @param[out] basisSecDerivVec = the values of \f$ d^2/dx^2 B \f$ at each x
   *  @param[out] basisThirdDerivVec = the values of \f$ d^3/dx^3 B \f$ at each x
   *
   */
   void ClmGaussianBasis
   ( double const zetaExp,
     int const l,
     int const npt,
     double const * xVec,
     double *  basisVec,
     double *  basisDerivVec,
     double *  basisSecDerivVec,
     double *  basisThirdDerivVec);
  




























  /** \overload */
   void YlmGaussianBasis
   ( double const zetaExp,
     int const l,
     int const npt,
     double const *  xVec,
     double *  basisVec );
  

  /** \overload */
   void YlmGaussianBasis
   ( double const zetaExp,
     int const l,
     int const npt,
     double const *  xVec,
     double *  basisVec,
     double *  basisDerivVec );



  /** \overload */
   void YlmGaussianBasis
   ( double const zetaExp,
     int const l,
     int const npt,
     double const *  xVec,
     double *  basisVec,
     double *  basisDerivVec,
     double *  basisSecDerivVec);

   



  /** \brief Gaussian radial basis
   *
   *
   * This basis \f$B_n(r;l,\zeta)\f$ has the form
   * \f[
   *   B(r;l,\zeta) = N \mbox{e}^{-\zeta r^2} r^l
   * \f]
   * where
   * \f[
   *   N = \sqrt{ 8 \frac{ (4\zeta)^l }{ (2l+1)!! } \sqrt{ \frac{ (2\zeta)^l }{ \pi } } }
   * \f]
   * is a normalization constant such that
   *
   * \f[
   *   \int B(r;l,\zeta) Y_{l,m}(\Omega) B(r;l',\zeta) Y_{l',m'}(\Omega) d^3r =
   *       \delta_{l,l'} \delta_{m,m'}
   * \f]
   *
   *
   *  @param[in] zetaExp = \f$ \zeta \f$
   *  @param[in] l = principle angular momentum
   *  @param[in] npt = number of x-values
   *  @param[in] xVec = an array of x-values
   *  @param[out] basisVec = the values of \f$ B(r;n,l,\zeta) \f$ at each x
   *  @param[out] basisDerivVec = the values of \f$ d/dx B \f$ at each x
   *  @param[out] basisSecDerivVec = the values of \f$ d^2/dx^2 B \f$ at each x
   *  @param[out] basisThirdDerivVec = the values of \f$ d^3/dx^3 B \f$ at each x
   *
   */
   void YlmGaussianBasis
   ( double const zetaExp,
     int const l,
     int const npt,
     double const * xVec,
     double *  basisVec,
     double *  basisDerivVec,
     double *  basisSecDerivVec,
     double *  basisThirdDerivVec);
  
}

#endif

