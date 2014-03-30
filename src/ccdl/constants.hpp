#ifndef _CCDL_CONSTANTSMOD_HPP_
#define _CCDL_CONSTANTSMOD_HPP_


namespace ccdl
{
  /** \var PI
   * Fundamental unitless constant
   */
  double const PI          = 3.141592653589793238462643383279502884197;

  
  /** \var RAD_PER_DEG
   * Fundamental unitless constant  PI/180
   *
   * Converts degree to radients
   */
  double const RAD_PER_DEG    = PI/180.0;

  
  /** \var PIO2
   * Fundamental unitless constant
   */
  double const PIO2        = PI * 0.5;


  /** \var TWO_PI
   * Fundamental unitless constant
   */
  double const TWO_PI      = 2.0 * PI;


  /** \var FOUR_PI
   * Fundamental unitless constant
   */
  double const FOUR_PI     = 4.0 * PI;


  /** \var SQRT_PI
   * Fundamental unitless constant
   */
  double const SQRT_PI     = 1.772453850905516027298167483341;

  /** \var TWO_OVER_SQRT_PI
   * Fundamental unitless constant
   */
  double const TWO_OVER_SQRT_PI = 2./SQRT_PI;


  /** \var SQRT2
   * Fundamental unitless constant
   */
  double const SQRT2       = 1.41421356237309504880168872420969807856967;


  /** \var SQRT3
   * Fundamental unitless constant
   */
  double const SQRT3       = 1.73205080756887729352744634150587236694281;


  /** \var EULER
   * Fundamental unitless constant
   */
  double const EULER       = 0.5772156649015328606065120900824024310422;


  /** \var DEG_PER_RAD
   * Fundamental unitless constant
   *
   * Converts a number from radians to degrees
   */
  double const DEG_PER_RAD     = 180.0/PI;


  /** \var AVOGADRO_CONSTANT 
   * Unitless 
   */
  double const AVOGADRO_CONSTANT       = 6.02214179e+23;



  /** \var SPEED_OF_LIGHT_SI
   * SI
   */
  double const SPEED_OF_LIGHT_SI          = 299792458.0;



  /** \var BOLTZMANN_CONSTANT_SI
   * SI
   */
  double const BOLTZMANN_CONSTANT_SI      = 1.3806504e-23;



  /** \var PLANCK_CONSTANT_SI
   * SI
   */
  double const PLANCK_CONSTANT_SI         = 6.62606896e-34;



  /** \var HBAR_SI
   * SI
   */
  double const HBAR_SI                    = PLANCK_CONSTANT_SI/(2*PI);



  /** \var ELECTRON_CHARGE_SI
   * SI
   */
  double const ELECTRON_CHARGE_SI         = 1.602176487e-19;



  /** \var ELECTRON_REST_MASS_SI
   * SI
   */
  double const ELECTRON_REST_MASS_SI      = 9.10938215e-31;



  /** \var PROTON_REST_MASS_SI
   * SI
   */
  double const PROTON_REST_MASS_SI        = 1.67262171e-27;



  /** \var NEUTRON_REST_MASS_SI
   * SI
   */
  double const NEUTRON_REST_MASS_SI       = 1.67492729e-27;



  /** \var PERMITTIVITY_FREE_SPACE_SI
   * SI
   */
  double const PERMITTIVITY_FREE_SPACE_SI = 8.8541878176e-12;




  /** \var FARADAY_CONSTANT_SI
   * SI; derived
   */
  double const FARADAY_CONSTANT_SI        = AVOGADRO_CONSTANT*ELECTRON_CHARGE_SI;



  /** \var FINE_STRUCTURE_CONSTANT_SI
   * SI; derived
   */
  double const FINE_STRUCTURE_CONSTANT_SI = ELECTRON_CHARGE_SI*ELECTRON_CHARGE_SI / (2 * PERMITTIVITY_FREE_SPACE_SI * PLANCK_CONSTANT_SI * SPEED_OF_LIGHT_SI);



  /** \var RYDBERG_CONSTANT_SI
   * SI; derived
   */
  double const RYDBERG_CONSTANT_SI        =  FINE_STRUCTURE_CONSTANT_SI*FINE_STRUCTURE_CONSTANT_SI * ELECTRON_REST_MASS_SI * SPEED_OF_LIGHT_SI / (2 * PLANCK_CONSTANT_SI);
  
  
  
  
  //  ! REPLACED WITH CODATA FROM NIST
  //  !  double const FINE_STRUCTURE_CONSTANT = 7.2973525376E-3
  //  !  double const RYDBERG_CONSTANT = 10973731.568527 ! 1/m
  //  !
  //  ! Atomic Units in terms of SI
  //  !
  //  ! Independent quantities
  //  !

  /** \var MASS_SI_PER_AU
   * SI/AU (independent quantity)
   */
  double const MASS_SI_PER_AU             = ELECTRON_REST_MASS_SI;



  /** \var CHARGE_SI_PER_AU
   * SI/AU (independent quantity)
   */
  double const CHARGE_SI_PER_AU           = ELECTRON_CHARGE_SI;



  /** \var ANGULAR_MOMENTUM_SI_PER_AU
   * SI/AU (independent quantity)
   */
  double const ANGULAR_MOMENTUM_SI_PER_AU = HBAR_SI;



  /** \var PERMITTIVITY_OF_FREE_SPACE_SI_PER_AU
   * SI/AU (independent quantity)
   */
  double const PERMITTIVITY_OF_FREE_SPACE_SI_PER_AU     = 4*PI*PERMITTIVITY_FREE_SPACE_SI;


  //  !
  //  ! Dependent quantities
  //  !


  /** \var DISTANCE_SI_PER_AU
   * SI/AU (derived)
   */
  double const DISTANCE_SI_PER_AU = (PERMITTIVITY_OF_FREE_SPACE_SI_PER_AU*ANGULAR_MOMENTUM_SI_PER_AU*ANGULAR_MOMENTUM_SI_PER_AU)/(MASS_SI_PER_AU*CHARGE_SI_PER_AU*CHARGE_SI_PER_AU);



  /** \var ENERGY_SI_PER_AU
   * SI/AU (derived)
   */
  double const ENERGY_SI_PER_AU = (CHARGE_SI_PER_AU*CHARGE_SI_PER_AU)/(PERMITTIVITY_OF_FREE_SPACE_SI_PER_AU*DISTANCE_SI_PER_AU);



  /** \var VELOCTY_SI_PER_AU
   * SI/AU (derived)
   */
  double const VELOCTY_SI_PER_AU = (CHARGE_SI_PER_AU*CHARGE_SI_PER_AU)/(PERMITTIVITY_OF_FREE_SPACE_SI_PER_AU*ANGULAR_MOMENTUM_SI_PER_AU);



  /** \var TIME_SI_PER_AU
   * SI/AU (derived)
   */
  double const TIME_SI_PER_AU = (DISTANCE_SI_PER_AU/VELOCTY_SI_PER_AU);



  /** \var ELECTROSTATIC_POTENTIAL_SI_PER_AU
   * SI/AU (derived)
   */
  double const ELECTROSTATIC_POTENTIAL_SI_PER_AU = (CHARGE_SI_PER_AU)/(PERMITTIVITY_OF_FREE_SPACE_SI_PER_AU*DISTANCE_SI_PER_AU);



  /** \var FORCE_SI_PER_AU
   * SI/AU (derived)
   */
  double const FORCE_SI_PER_AU = (CHARGE_SI_PER_AU*CHARGE_SI_PER_AU)/(PERMITTIVITY_OF_FREE_SPACE_SI_PER_AU*DISTANCE_SI_PER_AU*DISTANCE_SI_PER_AU);



  /** \var MAGNETIC_DIPOLE_SI_PER_AU
   * SI/AU (derived)
   */
  double const MAGNETIC_DIPOLE_SI_PER_AU = (CHARGE_SI_PER_AU*ANGULAR_MOMENTUM_SI_PER_AU)/MASS_SI_PER_AU;



  /** \var PRESSURE_SI_PER_AU
   * SI/AU (derived)
   */
  double const PRESSURE_SI_PER_AU = FORCE_SI_PER_AU / (DISTANCE_SI_PER_AU*DISTANCE_SI_PER_AU);




  //  !
  //  ! Fundamental physical constants (AU)
  //  !

  /** \var SPEED_OF_LIGHT_AU
   * AU
   */
  double const SPEED_OF_LIGHT_AU     = SPEED_OF_LIGHT / VELOCTY_SI_PER_AU;



  /** \var BOLTZMANN_CONSTANT_AU
   * AU
   */
  double const BOLTZMANN_CONSTANT_AU = BOLTZMANN_CONSTANT_SI/ENERGY_SI_PER_AU;



  /** \var PLANCK_CONSTANT_AU
   * AU
   */
  double const PLANCK_CONSTANT_AU = PLANCK_CONSTANT_SI/(TIME_SI_PER_AU*ENERGY_SI_PER_AU);



  /** \var PROTON_REST_MASS_AU
   * AU
   */
  double const PROTON_REST_MASS_AU = PROTON_REST_MASS_SI / MASS_SI_PER_AU;



  /** \var NEUTRON_REST_MASS_AU
   * AU
   */
  double const NEUTRON_REST_MASS_AU = NEUTRON_REST_MASS_SI / MASS_SI_PER_AU;



  /** \var RYDBERG_CONSTANT_AU
   * AU
   */
  double const RYDBERG_CONSTANT_AU = RYDBERG_CONSTANT_SI * DISTANCE_SI_PER_AU;




  //  !
  //  ! Unit conversions (in terms of AU)
  //  !


  /** \var AU_PER_PASCAL
   * AU/PASCAL (Conversion factor)
   */
  double const AU_PER_PASCAL = 1.0 / PRESSURE_SI_PER_AU;



  /** \var AU_PER_ATMOSPHERE
   * AU/ATMOSPHERE  (Conversion factor)
   */
  double const AU_PER_ATMOSPHERE =  AU_PER_PASCAL * 101325.0;



  /** \var AU_PER_ELECTRON_VOLT
   * AU/ELECTRON_VOLT  (Conversion factor)
   */
  double const AU_PER_ELECTRON_VOLT = CHARGE_SI_PER_AU/ENERGY_SI_PER_AU;



  /** \var AU_PER_JOULE_PER_MOL
   * AU/JOULE_PER_MOL  (Conversion factor)
   */
  double const AU_PER_JOULE_PER_MOL = 1.0/(ENERGY_SI_PER_AU*AVOGADRO_CONSTANT);



  /** \var AU_PER_KJOULE_PER_MOL
   * AU/KJOULE_PER_MOL  (Conversion factor)
   */
  double const AU_PER_KJOULE_PER_MOL = 1000.0 * AU_PER_JOULE_PER_MOL;



  /** \var AU_PER_KCAL_PER_MOL
   * AU/KCAL_PER_MOL  (Conversion factor)
   */
  double const AU_PER_KCAL_PER_MOL = 4.184 * AU_PER_KJOULE_PER_MOL;



  /** \var AU_PER_KCAL_PER_MOL_ANG
   * AU/KCAL_PER_ANG_MOL  (Conversion factor)
   */
  double const AU_PER_KCAL_PER_MOL_ANG = 4.184e13 / (FORCE_SI_PER_AU * AVOGADRO_CONSTANT);



  /** \var AU_PER_INVERSE_CM
   * AU/INVERSE_CM  (Conversion factor)
   */
  double const AU_PER_INVERSE_CM = (100.*SPEED_OF_LIGHT_SI*PLANCK_CONSTANT_SI)/ENERGY_SI_PER_AU;



  /** \var AU_PER_ANGSTROM
   * AU/ANGSTROM  (Conversion factor)
   */
  double const AU_PER_ANGSTROM = 1.0e-10/DISTANCE_SI_PER_AU;



  /** \var AU_PER_DEBYE
   * AU/DEBYE  (Conversion factor)
   */
  double const AU_PER_DEBYE = 3.335641e-30/(CHARGE_SI_PER_AU*DISTANCE_SI_PER_AU);



  /** \var AU_PER_GRAMS_PER_MOL
   * AU/GRAMS_PER_MOL  (Conversion factor)
   */
  double const AU_PER_GRAMS_PER_MOL = 1.0 / (1000.0 * AVOGADRO_CONSTANT  * MASS_SI_PER_AU);



  /** \var AU_PER_CM3_PER_MOL
   * AU/CM3_PER_MOL  (Conversion factor)
   */
  double const AU_PER_CM3_PER_MOL =  AVOGADRO_CONSTANT * (100.0 * DISTANCE_SI_PER_AU)*(100.0 * DISTANCE_SI_PER_AU)*(100.0 * DISTANCE_SI_PER_AU);



  /** \var AU_PER_ATOMIC_MASS_UNIT
   * AU/ATOMIC_MASS_UNIT  (Conversion factor)
   */
  double const AU_PER_ATOMIC_MASS_UNIT =  1.0/(AVOGADRO_CONSTANT * 1000.0 * MASS_SI_PER_AU); //  ! au/amu


  //  !
  //  ! ENERGY UNIT OF RYDBERG in AU
  //  !

  /** \var AU_PER_RYDBERG
   * AU/RYDBERG  (Conversion factor)
   */
  double const AU_PER_RYDBERG = RYDBERG_CONSTANT_SI * AU_PER_INVERSE_CM / 100.;



  /** \var AU_PER_DYNE_PER_CM
   * AU/DYNE_PER_CM
   */
  double const AU_PER_DYNE_PER_CM =  0.001 * DISTANCE_SI_PER_AU/FORCE_SI_PER_AU;

}

#endif


