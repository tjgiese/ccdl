#ifndef _GaussianBasisDatabase_h_
#define _GaussianBasisDatabase_h_

#include <cstddef>

namespace GaussianBasis
{
   void augccpVDZ( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void augccpVQZ( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void augccpVTZ( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void ccpV5Z( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void ccpVDZ( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void ccpVQZ( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void ccpVTZ( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void MG3S( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void Pople321G( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void Pople321GSP( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void Pople6311G( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void Pople6311Gss( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void Pople6311ppGss( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void Pople631G( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void Pople631Gs( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void Pople631Gss( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void Pople631ppGss( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void STO3G( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
   void STO6G( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT );
}

#endif
