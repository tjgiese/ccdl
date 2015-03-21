#ifndef _GiesePointMultipole_H_
#define _GiesePointMultipole_H_

#include "GieseMultipole.hpp"

namespace ccdl
{
void SolidHarm_Ilm_S( double const * crd, double const r2, double * Ilm );
void SolidHarm_Ilm_P( double const * crd, double const r2, double * Ilm );
void SolidHarm_Ilm_D( double const * crd, double const r2, double * Ilm );
void SolidHarm_Ilm_F( double const * crd, double const r2, double * Ilm );
void SolidHarm_Ilm_G( double const * crd, double const r2, double * Ilm );
void SolidHarm_Ilm_H( double const * crd, double const r2, double * Ilm );
void SolidHarm_Ilm_I( double const * crd, double const r2, double * Ilm );
void SolidHarm_Ilm_J( double const * crd, double const r2, double * Ilm );
void SolidHarm_Ilm_K( double const * crd, double const r2, double * Ilm );
void SolidHarm_Ilm_L( double const * crd, double const r2, double * Ilm );
void SolidHarm_Ilm_M( double const * crd, double const r2, double * Ilm );
void SolidHarm_Ilm_N( double const * crd, double const r2, double * Ilm );
void SolidHarm_dIlm_S( double const * Ilm, double * dIlm );
void SolidHarm_dIlm_P( double const * Ilm, double * dIlm );
void SolidHarm_dIlm_D( double const * Ilm, double * dIlm );
void SolidHarm_dIlm_F( double const * Ilm, double * dIlm );
void SolidHarm_dIlm_G( double const * Ilm, double * dIlm );
void SolidHarm_dIlm_H( double const * Ilm, double * dIlm );
void SolidHarm_dIlm_I( double const * Ilm, double * dIlm );
void SolidHarm_dIlm_J( double const * Ilm, double * dIlm );
void SolidHarm_dIlm_K( double const * Ilm, double * dIlm );
void SolidHarm_dIlm_L( double const * Ilm, double * dIlm );
void SolidHarm_dIlm_M( double const * Ilm, double * dIlm );
void IlmInteraction_SS( double const * crd, double const r2, double * T );
void IlmInteraction_PS( double const * crd, double const r2, double * T );
void IlmInteraction_PP( double const * crd, double const r2, double * T );
void IlmInteraction_DS( double const * crd, double const r2, double * T );
void IlmInteraction_DP( double const * crd, double const r2, double * T );
void IlmInteraction_DD( double const * crd, double const r2, double * T );
void IlmInteraction_FS( double const * crd, double const r2, double * T );
void IlmInteraction_FP( double const * crd, double const r2, double * T );
void IlmInteraction_FD( double const * crd, double const r2, double * T );
void IlmInteraction_FF( double const * crd, double const r2, double * T );
void IlmInteraction_GS( double const * crd, double const r2, double * T );
void IlmInteraction_GP( double const * crd, double const r2, double * T );
void IlmInteraction_GD( double const * crd, double const r2, double * T );
void IlmInteraction_GF( double const * crd, double const r2, double * T );
void IlmInteraction_GG( double const * crd, double const r2, double * T );
void IlmInteraction_HS( double const * crd, double const r2, double * T );
void IlmInteraction_HP( double const * crd, double const r2, double * T );
void IlmInteraction_HD( double const * crd, double const r2, double * T );
void IlmInteraction_HF( double const * crd, double const r2, double * T );
void IlmInteraction_HG( double const * crd, double const r2, double * T );
void IlmInteraction_HH( double const * crd, double const r2, double * T );
void PtExpInt_SS( double const * crd, double const r2, double * T );
void PtExpInt_PS( double const * crd, double const r2, double * T );
void PtExpInt_PP( double const * crd, double const r2, double * T );
void PtExpInt_DS( double const * crd, double const r2, double * T );
void PtExpInt_DP( double const * crd, double const r2, double * T );
void PtExpInt_DD( double const * crd, double const r2, double * T );
void PtExpInt_FS( double const * crd, double const r2, double * T );
void PtExpInt_FP( double const * crd, double const r2, double * T );
void PtExpInt_FD( double const * crd, double const r2, double * T );
void PtExpInt_FF( double const * crd, double const r2, double * T );
void PtExpInt_GS( double const * crd, double const r2, double * T );
void PtExpInt_GP( double const * crd, double const r2, double * T );
void PtExpInt_GD( double const * crd, double const r2, double * T );
void PtExpInt_GF( double const * crd, double const r2, double * T );
void PtExpInt_GG( double const * crd, double const r2, double * T );
void PtExpInt_HS( double const * crd, double const r2, double * T );
void PtExpInt_HP( double const * crd, double const r2, double * T );
void PtExpInt_HD( double const * crd, double const r2, double * T );
void PtExpInt_HF( double const * crd, double const r2, double * T );
void PtExpInt_HG( double const * crd, double const r2, double * T );
void PtExpInt_HH( double const * crd, double const r2, double * T );
void PtExpGrd_SS( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_PS( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_PP( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_DS( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_DP( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_DD( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_FS( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_FP( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_FD( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_FF( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_GS( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_GP( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_GD( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_GF( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_GG( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_HS( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_HP( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_HD( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_HF( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_HG( double const * crd, double const r2, double * T, double * dT );
void PtExpGrd_HH( double const * crd, double const r2, double * T, double * dT );
void PtExpPot_SS( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_PS( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_PP( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_DS( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_DP( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_DD( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_FS( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_FP( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_FD( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_FF( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_GS( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_GP( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_GD( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_GF( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_GG( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_HS( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_HP( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_HD( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_HF( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_HG( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPot_HH( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb );
void PtExpPotGrd_SS( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_PS( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_PP( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_DS( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_DP( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_DD( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_FS( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_FP( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_FD( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_FF( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_GS( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_GP( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_GD( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_GF( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_GG( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_HS( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_HP( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_HD( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_HF( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_HG( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void PtExpPotGrd_HH( double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * g );
void SolidHarm_Ilm_Switch( int const lmax, double const * crd, double const r2, double * X );
void SolidHarm_dIlm_Switch( int const lmax, double const * Y, double * dY );
void IlmInteraction_Switch( int const la, int const lb, double const * crd, double const r2, double * X );
void PtExpInt_Switch( int const la, int const lb, double const * crd, double const r2, double * X );
void PtExpGrd_Switch( int const la, int const lb, double const * crd, double const r2, double * X, double * dX );
void PtExpPot_Switch( int const la, int const lb, double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb  );
void PtExpPotGrd_Switch( int const la, int const lb,  double const * crd, double const r2, double const * qa, double const * qb, double * pa, double * pb, double * grd  );

}


#include "GiesePointMultipoleInline.hpp"

#endif

