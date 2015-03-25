#ifndef _GieseMultipoleGen_hpp_
#define _GieseMultipoleGen_hpp_

namespace ccdl
{


  void AuxExpInt( int const la, int const lb, 
		  double const * crd,
		  double const * O,
		  double * T );

  void AuxExpGrd( int const la, int const lb, 
		  double const * crd,
		  double const * O,
		  double * T,
		  double * dT );

  void AuxExpPot( int const la, int const lb, 
		  double const * crd,
		  double const * O,
		  double const * qa,
		  double const * qb,
		  double * pa,
		  double * pb );

  void AuxExpPotGrd( int const la, int const lb, 
		     double const * crd,
		     double const * O,
		     double const * qa,
		     double const * qb,
		     double * pa,
		     double * pb,
		     double * g );


  void PrimGauAuxVec_Coulomb( int const order, 
			      double const zab, 
			      double const r2, 
			      double * O );

  void PrimGauAuxVec_Overlap( int const order, 
			      double const zab, 
			      double const r2, 
			      double * O );

  void PrimGauAuxVec_Ewald( int const order, 
			    double const sqrt_za, 
			    double const r2, 
			    double * O );


  void PrimGauExpInt_Coulomb( int const la, int const lb, 
			      double const zab,
			      double const * crd,
			      double const r2,
			      double * T );

  void PrimGauExpInt_Overlap( int const la, int const lb, 
			      double const zab,
			      double const * crd,
			      double const r2,
			      double * T );

  void PrimGauExpInt_Ewald( int const la, int const lb, 
			    double const sqrt_za,
			    double const * crd,
			    double const r2,
			    double * T );


  void PrimGauExpGrd_Coulomb( int const la, int const lb, 
			      double const zab,
			      double const * crd,
			      double const r2,
			      double * T,
			      double * dT );

  void PrimGauExpGrd_Overlap( int const la, int const lb, 
			      double const zab,
			      double const * crd,
			      double const r2,
			      double * T,
			      double * dT );

  void PrimGauExpGrd_Ewald( int const la, int const lb, 
			    double const zab,
			    double const * crd,
			    double const r2,
			    double * T,
			    double * dT );


  void PrimGauExpPot_Coulomb( int const la, int const lb, 
			      double const zab,
			      double const * crd,
			      double const r2,
			      double const * qa, double const * qb,
			      double * pa, double * pb );

  void PrimGauExpPot_Overlap( int const la, int const lb, 
			      double const zab,
			      double const * crd,
			      double const r2,
			      double const * qa, double const * qb,
			      double * pa, double * pb );

  void PrimGauExpPot_Ewald( int const la, int const lb, 
			    double const zab,
			    double const * crd,
			    double const r2,
			    double const * qa, double const * qb,
			    double * pa, double * pb );


  void PrimGauExpPotGrd_Coulomb( int const la, int const lb, 
				 double const zab,
				 double const * crd,
				 double const r2,
				 double const * qa, double const * qb,
				 double * pa, double * pb, double * g );

  void PrimGauExpPotGrd_Overlap( int const la, int const lb, 
				 double const zab,
				 double const * crd,
				 double const r2,
				 double const * qa, double const * qb,
				 double * pa, double * pb, double * g );

  void PrimGauExpPotGrd_Ewald( int const la, int const lb, 
			       double const zab,
			       double const * crd,
			       double const r2,
			       double const * qa, double const * qb,
			       double * pa, double * pb, double * g );
}

#endif
