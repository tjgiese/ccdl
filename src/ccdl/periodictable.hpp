/**
 * \file PeriodicTable.hpp
 * \file PeriodicTable.hpp
 *
 * \brief Contains atom data tables
 * 
 * \author giese
 *
 * \version 0.1
 *
 * \date 2010/10/15
 * 
 * 
 */


#ifndef _ccdl_PERIODICTABLEMOD_H
#define _ccdl_PERIODICTABLEMOD_H


#include "constants.hpp"


/** \private */
namespace ccdl
{

  /** \private */
  namespace PeriodicTable
  {

    /** \private NUM_ATOMS */
    unsigned int const NUM_ATOMS = 55;



    /** \private ElementSymbolData */
    char const * const ElementSymbolData[] = {"XX","H ","He","Li","Be","B ","C ","N ","O ","F ","Ne","Na","Mg","Al","Si","P ","S ","Cl","Ar","K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe"};



    /** \private ElementNameData */
    char const * const ElementNameData[] = {"Dummy          ","Hydrogen       ","Helium         ","Lithium        ","Beryllium      ","Boron          ","Carbon         ","Nitrogen       ","Oxygen         ","Fluorine       ","Neon           ","Sodium         ","Magnesium      ","Aluminum       ","Silicon        ","Phosphorus     ","Sulfur         ","Chlorine       ","Argon          ","Potassium      ","Calcium        ","Scandium       ","Titanium       ","Vanadium       ","Chromium       ","Manganese      ","Iron           ","Cobalt         ","Nickel         ","Copper         ","Zinc           ","Gallium        ","Germanium      ","Arsenic        ","Selenium       ","Bromine        ","Krypton        ","Rubidium       ","Strontium      ","Yttrium        ","Zirconium      ","Niobium        ","Molybdenum     ","Technetium     ","Ruthenium      ","Rhodium        ","Palladium      ","Silver         ","Cadmium        ","Indium         ","Tin            ","Antimony       ","Tellurium      ","Iodine         ","Xenon          "};



    /** \private AtomicWeightData */
    double const AtomicWeightData[] = {   0.0, 1.0079000,   4.0026000,   6.9400000,   9.0121800,  10.8100000,  12.0110000,  14.0067000,  15.9994000,  18.9984000,  20.1790000,  22.9897700,  24.3050000,  26.9815400,  28.0855000,  30.9737600,  32.0600000,  35.4530000,  39.9480000,  39.0983000,  40.0800000,  44.9559000,  47.9000000,  50.9415000,  51.9960000,  54.9380000,  55.8470000,  58.9332000,  58.7100000,  63.5460000,  65.3800000,  69.7350000,  72.5900000,  74.9216000,  78.9600000,  79.9040000,  83.8000000,  85.4678000,  87.6200000,  88.9059000,  91.2200000,  92.9064000,  95.9400000,  98.9062000, 101.0700000, 102.9055000, 106.4000000, 107.8680000, 112.4100000, 114.8200000, 118.6900000, 121.7500000, 127.6000000, 126.9045000, 131.3000000};



    /** \private StandardIsotopeMassData */
    double const StandardIsotopeMassData[] = {   0.0,  1.0078250,   4.0026000,   7.0160000,   9.0121800,  11.0093100,  12.0000000,  14.0030700,  15.9949100,  18.9984000,  19.9924400,  22.9898000,  23.9850400,  26.9815300,  27.9769300,  30.9937600,  31.9720700,  34.9688500,  39.9627200,  38.9637100,  39.9625900,  44.9559200,  47.9000000,  50.9440000,  51.9405000,  54.9381000,  55.9349000,  58.9332000,  57.9353000,  62.9298000,  63.9291000,  68.9257000,  73.9219000,  74.9216000,  79.9165000,  78.9183000,  83.8000000,  84.9117000,  87.9056000,  88.9054000,  89.9043000,  92.9060000,  97.9055000,  98.9062000, 101.9037000, 102.9048000, 105.9032000, 106.9041000, 113.9036000, 114.9041000, 117.9018000, 120.9038000, 129.9067000, 126.9004000, 131.9042000};



    /** \private SecondIsotopeMassData */
    double const SecondIsotopeMassData[] = {   0.0,  2.0141020,   3.0160290,   6.0151220,   0.0000000,  10.0129370,  13.0033550,  15.0001090,  17.9991600,   0.0000000,  21.9913860,   0.0000000,  25.9825930,   0.0000000,  28.9764950,   0.0000000,  33.9678670,  36.9659030,  35.9675460,  40.9618260,  43.9554810,   0.0000000,  45.9526300,  49.9471630,  52.9406540,   0.0000000,  53.9396150,   0.0000000,  59.9307910,  64.9277940,  65.9260370,  70.9247050,  71.9220760,   0.0000000,  77.9173100,  80.9162910,  85.9106100,  86.9091840,  85.9092620,   0.0000000,  93.9063160,   0.0000000,  95.9046790,  97.9072160, 103.9054300,   0.0000000, 107.9038940, 108.9047560, 111.9027570, 112.9040610, 117.9016060, 122.9042160, 127.9044610,   0.0000000, 128.9047800};



    /** \private CovalentRadiusData */
    double const CovalentRadiusData[] = {   0.0,  0.3700000,   0.3200000,   1.3400000,   0.9000000,   0.8200000,   0.7700000,   0.7500000,   0.7300000,   0.7100000,   0.6900000,   1.5400000,   1.3000000,   1.1800000,   1.1100000,   1.0600000,   1.0200000,   0.9900000,   0.9700000,   1.9600000,   1.7400000,   1.4400000,   1.3600000,   1.2500000,   1.2700000,   1.3900000,   1.2500000,   1.2600000,   1.2100000,   1.3800000,   1.3100000,   1.2600000,   1.2200000,   1.1900000,   1.1600000,   1.1400000,   1.1000000,   2.1100000,   1.9200000,   1.6200000,   1.4800000,   1.3700000,   1.4500000,   1.5600000,   1.2600000,   1.3500000,   1.3100000,   1.5300000,   1.4800000,   1.4400000,   1.4100000,   1.3800000,   1.3500000,   1.3300000,   1.3000000};



    /** \private vdWRadiusData */
    double const vdWRadiusData[] = {   0.0,  1.2000000,   1.4000000,   1.8200000,   0.0000000,   0.0000000,   1.7000000,   1.5500000,   1.5200000,   1.4700000,   1.5400000,   2.2700000,   1.7300000,   0.0000000,   2.1000000,   1.8000000,   1.8000000,   1.7500000,   1.8800000,   2.7500000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,   1.6300000,   1.4000000,   1.3900000,   1.8700000,   0.0000000,   1.8500000,   1.9000000,   1.8500000,   2.0200000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,   1.6300000,   1.7200000,   1.5800000,   1.9300000,   2.1700000,   0.0000000,   2.0600000,   1.9800000,   2.1600000};



    /** \private ElectronegativityData */
    double const ElectronegativityData[] = {   0.0, 0.2638740,   0.4479860,   0.1106210,   0.1800810,   0.1576630,   0.2304300,   0.2682840,   0.2771040,   0.3825800,   0.3905280,   0.1047410,   0.1378170,   0.1187060,   0.1753030,   0.2065420,   0.2285920,   0.3050350,   0.2819480,   0.0889380,   0.0808530,   0.1227490,   0.1267920,   0.1323040,   0.1367140,   0.1367140,   0.1492100,   0.1580300,   0.1617050,   0.1646450,   0.1635430,   0.1176040,   0.1690550,   0.1947810,   0.2164650,   0.2789420,   0.2496070,   0.0859980,   0.0735020,   0.1172360,   0.1337740,   0.1470050,   0.1433300,   0.0000000,   0.1653800,   0.1580300,   0.1635430,   0.1631750,   0.1591330,   0.1139290,   0.1580300,   0.1782430,   0.2017640,   0.2484380,   0.2152490};



    /** \private HardnessData */
    double const HardnessData[] = {   0.0, 0.2363100,   0.4556300,   0.0878350,   0.1653800,   0.1473720,   0.1837560,   0.2657110,   0.2234470,   0.2576260,   0.4019940,   0.0845280,   0.1433300,   0.1018010,   0.1242190,   0.1793460,   0.1521500,   0.1719960,   0.2972370,   0.0705620,   0.1470050,   0.1176040,   0.1238520,   0.1139290,   0.1124590,   0.1367140,   0.1400220,   0.1323040,   0.1194410,   0.1194410,   0.1815510,   0.1065780,   0.1249540,   0.1653800,   0.1422270,   0.1550900,   0.2648950,   0.0679900,   0.1359790,   0.1172360,   0.1179710,   0.1102540,   0.1139290,   0.0000000,   0.1102540,   0.1161330,   0.1429620,   0.1153990,   0.1712610,   0.1029030,   0.1120910,   0.1396550,   0.1293640,   0.1356120,   0.2305380};



    /** \private IonizationPotentialData */
    double const IonizationPotentialData[] = {   0.0, 0.4998100,   0.9036170,   0.1981710,   0.3388570,   0.3049900,   0.4139050,   0.5342100,   0.5005330,   0.6403810,   0.7926480,   0.1888760,   0.2810290,   0.2200000,   0.2996190,   0.3854480,   0.3808000,   0.4766480,   0.5792760,   0.1595430,   0.2246860,   0.2411810,   0.2509710,   0.2479620,   0.2487240,   0.2732570,   0.2904760,   0.2896760,   0.2808000,   0.2840000,   0.3452950,   0.2204950,   0.2902860,   0.3607620,   0.3584760,   0.4342480,   0.5145900,   0.1535240,   0.2093330,   0.2285710,   0.2438480,   0.2484190,   0.2606860,   0.2674290,   0.2705520,   0.2741710,   0.3064380,   0.2784760,   0.3305900,   0.2126860,   0.2699430,   0.3177140,   0.3311620,   0.3841520,   0.4458670};



    /** \private ElectronAffinityData */
    double const ElectronAffinityData[] = {   0.0, 0.0277330,  -0.0076440,   0.0227050,   0.0000000,   0.0101710,   0.0586290,   0.0026670,   0.0537140,   0.1249520,  -0.0114670,   0.0201140,   0.0000000,   0.0161910,   0.0508950,   0.0274290,   0.0761900,   0.1329520,  -0.0152880,   0.0184380,   0.0000000,   0.0068950,   0.0028950,   0.0192760,   0.0244950,   0.0000000,   0.0059810,   0.0242670,   0.0426670,   0.0451050,   0.0000000,   0.0110100,   0.0453330,   0.0297140,   0.0742860,   0.1236570,  -0.0152880,   0.0178670,   0.0000000,   0.0112760,   0.0156570,   0.0328000,   0.0273900,   0.0201900,   0.0385900,   0.0417900,   0.0204570,   0.0478480,   0.0000000,   0.0110100,   0.0408760,   0.0393140,   0.0724570,   0.1124570,  -0.0152880};



    /** \private BraggSlaterRadiusData */
    double const BraggSlaterRadiusData[] = {   0.0, 0.2500000,   0.8500000,   1.4500000,   1.0500000,   0.8500000,   0.7000000,   0.6500000,   0.6000000,   0.5000000,   1.1500000,   1.8000000,   1.5000000,   1.2500000,   1.1000000,   1.0000000,   1.0000000,   1.0000000,   1.6000000,   2.2000000,   1.8000000,   1.6000000,   1.4000000,   1.3500000,   1.4000000,   1.4000000,   1.4000000,   1.3500000,   1.3500000,   1.3500000,   1.3500000,   1.3000000,   1.2500000,   1.1500000,   1.1500000,   1.1500000,   1.7500000,   2.3500000,   2.0000000,   1.8000000,   1.5500000,   1.4500000,   1.4500000,   1.3500000,   1.3000000,   1.3500000,   1.4000000,   1.6000000,   1.5500000,   1.5500000,   1.4500000,   1.4500000,   1.4000000,   1.4000000,   2.0000000};


  }
  




  /** \brief Naturally occuring isotopic average atomic weight
   *
   * @param[in] Z = atomic number
   *
   * \note Result is in atomic mass units (g/mol).  For atomic units of mass, convert with ccdl::ATOMIC_MASS_UNIT
   */
  double AtomicWeight( unsigned int const Z );




  /** \fn  double const StandardIsotopeMass( unsigned int const Z )
   * \brief Most naturally abundant isotope atomic weight
   *
   * @param[in] Z = atomic number
   *
   * \note Result is in atomic mass units (g/mol).  For atomic units of mass, convert with ccdl::ATOMIC_MASS_UNIT
   */
  double StandardIsotopeMass( unsigned int const Z );




  /** \brief Second-most naturally abundant isotope atomic weight
   *
   * @param[in] Z = atomic number
   *
   * \note Result is in atomic mass units (g/mol).  For atomic units of mass, convert with ccdl::ATOMIC_MASS_UNIT
   */
  double SecondIsotopeMass( unsigned int const Z );





  /** \brief Covalent radius
   *
   * @param[in] Z = atomic number
   *
   * \note Result is in atomic units (Bohr)
   * \note Values taken from www.webelements.com
   * \note R.T. Sanderson in Chemical Periodicity, Reinhold, New York, USA, 1962. 
   * \note L.E. Sutton (ed.) in Table of interatomic distances and configuration in molecules and ions, Supplement 1956-1959, Special publication No. 18, Chemical Society, London, UK, 1965. 
   * \note J.E. Huheey, E.A. Keiter, and R.L. Keiter, Inorganic Chemistry: Principles of structure and reactivity, 4th edition, HarperCollins, New York, 1993. 
   * \note W.W. Porterfield in Inorganic chemistry, a unified approach, Addison Wesley Publishing Co., Reading Massachusetts, USA, 1984. 
   * \note A.M. James and M.P. Lord in Macmillan's Chemical and Physical Data, Macmillan, London, UK, 1992. 
   */
  double CovalentRadius( unsigned int const Z );



  /** \brief van der Waals radius
   *
   * @param[in] Z = atomic number
   *
   * \note Result is in atomic units (Bohr)
   * \note Values taken from www.webelements.com
   * \note A. Bondi, J. Phys. Chem., 1964, 68, 441.
   */
  double vdWRadius( unsigned int const Z );





  /** \brief Bragg-Slater radius
   *
   * @param[in] Z = atomic number
   *
   * \note Result is in atomic units (Bohr)
   * \note J. C. Slater, J. Chem. Phys., 41, 3199 (1964)
   * \note Slater did not list rare gas values, so these are taken to be: He = (H+Li)/2; Ne = (F+Na)/2; Ar = (Cl+K)/2; Kr = (Br+Rb)/2; Xe = (I+Cs)/2
   */
  double BraggSlaterRadius( unsigned int const Z );



  /** \brief Electronegatvity (Pearson)
   *
   * @param[in] Z = atomic number
   *
   * \note Result is in atomic units (Hartree)
   * \note Noble gases: Klaus S. Lackner and George Zweig, Phys. Rev. D, 28:1671-1691, 1983.
   * \note All others: Pearson, Inorg. Chem. 27:734-740, 1988.
   * \note Data for Tc are currently missing.
   */
  double Electronegativity( unsigned int const Z );



  /** \brief Hardness (Pearson)
   *
   * @param[in] Z = atomic number
   *
   * \note Result is in atomic units (Hartree)
   * \note Noble gases: Klaus S. Lackner and George Zweig, Phys. Rev. D, 28:1671-1691, 1983.
   * \note All others: Pearson, Inorg. Chem. 27:734-740, 1988.
   * \note Data for Tc are currently missing.
   */
  double Hardness( unsigned int const Z );



  /** \brief Ionization potential
   *
   * @param[in] Z = atomic number
   *
   * \note Result is in atomic units (Hartree)
   * \note Taken from www.webelements.com
   * \note J.E. Huheey, E.A. Keiter, and R.L. Keiter, Inorganic Chemistry: Principles of structure and reactivity, 4th edition, HarperCollins, New York, 1993, ISBN 0-06-042995-X. 
   * \note D.R. Lide, (ed.) in Chemical Rubber Company handbook of chemistry and physics, CRC Press, Boca Raton, Florida, USA, 75th edition, 1994.
   * \note J.A. Dean (ed) in Lange's handbook of chemistry, McGraw-Hill, New York, USA, 14th edition, 1992. 
   */
  double IonizationPotential( unsigned int const Z );




  /** \brief Electron affinity
   *
   * @param[in] Z = atomic number
   *
   * \note Result is in atomic units (Hartree)
   * \note Noble gases: Klaus S. Lackner and George Zweig, Phys. Rev. D, 28:1671-1691, 1983.  
   * \note All others: www.webelements.com
   * \note J.E. Huheey, E.A. Keiter, and R.L. Keiter in Inorganic Chemistry: Principles of structure and reactivity, 4th edition, HarperCollins, New York, USA, 1993, ISBN 0-06-042995-X. 
   * \note A.M. James and M.P. Lord in Macmillan's Chemical and Physical Data, Macmillan, London, UK, 1992.
   */
  double ElectronAffinity( unsigned int const Z );



}




inline double ccdl::AtomicWeight( unsigned int const Z )
{
  return ccdl::PeriodicTable::AtomicWeightData[Z%100];
}

inline double ccdl::StandardIsotopeMass( unsigned int const Z )
{
  return ccdl::PeriodicTable::StandardIsotopeMassData[Z%100];
}

inline double ccdl::SecondIsotopeMass( unsigned int const Z )
{
  return ccdl::PeriodicTable::SecondIsotopeMassData[Z%100];
}


inline double ccdl::CovalentRadius( unsigned int const Z )
{
  return ccdl::PeriodicTable::CovalentRadiusData[Z%100] * ccdl::ANGSTROM;
}

inline double ccdl::vdWRadius( unsigned int const Z )
{
  return ccdl::PeriodicTable::vdWRadiusData[Z%100] * ccdl::ANGSTROM;
}

inline double ccdl::BraggSlaterRadius( unsigned int const Z )
{
  return ccdl::PeriodicTable::BraggSlaterRadiusData[Z%100] * ccdl::ANGSTROM;
}

inline double ccdl::Electronegativity( unsigned int const Z )
{
  return ccdl::PeriodicTable::ElectronegativityData[Z%100];
}

inline double ccdl::Hardness( unsigned int const Z )
{
  return ccdl::PeriodicTable::HardnessData[Z%100];
}

inline double ccdl::IonizationPotential( unsigned int const Z )
{
  return ccdl::PeriodicTable::IonizationPotentialData[Z%100];
}

inline double ccdl::ElectronAffinity( unsigned int const Z )
{
  return ccdl::PeriodicTable::ElectronAffinityData[Z%100];
}




#endif
