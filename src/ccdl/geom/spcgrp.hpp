#ifndef _spcgrp_hpp_
#define _spcgrp_hpp_

#include <vector>
#include <string>
#include <tr1/memory>

namespace spcgrp
{


   class symops
   {
      public:
         symops( char const * name, char const * lattice, int number, int nops, double const * begin, double const * end );
         virtual ~symops() {};
         virtual void SetLatVec( double a, double b, double c, double alp, double bet, double gam );
         double const * GetLatVec() const;
         void UpdateCell();
         std::vector<double> GenAtomPositions( double const * crd, bool is_fractional = false ) const;
         std::vector<double> GenForceSymmetries( double const * crd, double const * frc, bool is_fractional = false ) const;

      public:
         std::string name,lattice;
         int number,nops,lattice_code;
         std::vector<double> data;
         std::vector<double> latvec;
         std::vector<double> ucell,recip;
         std::vector<double> invdata,frcdata;
    };

}




namespace spcgrp
{
   std::tr1::shared_ptr< spcgrp::symops > make_symops( int space_group_number );
}




//////////////////////////////////////////////////////////////////////////////



namespace spcgrp
{


class triclinic : public symops
{
   public:
      triclinic( char const * name, int number, int nops, double const * begin, double const * end );
      virtual ~triclinic() {}
      virtual void SetLatVec( double a, double b, double c, double alp, double bet, double gam );
};


class monoclinic : public symops
{
   public:
      monoclinic( char const * name, int number, int nops, double const * begin, double const * end );
      virtual ~monoclinic() {}
      virtual void SetLatVec( double a, double b, double c, double alp, double bet, double gam );
};


class orthorhombic : public symops
{
   public:
      orthorhombic( char const * name, int number, int nops, double const * begin, double const * end );
      virtual ~orthorhombic() {}
      virtual void SetLatVec( double a, double b, double c, double alp, double bet, double gam );
};


class tetragonal : public symops
{
   public:
      tetragonal( char const * name, int number, int nops, double const * begin, double const * end );
      virtual ~tetragonal() {}
      virtual void SetLatVec( double a, double b, double c, double alp, double bet, double gam );
};


class trigonal : public symops
{
   public:
      trigonal( char const * name, int number, int nops, double const * begin, double const * end );
      virtual ~trigonal() {}
      virtual void SetLatVec( double a, double b, double c, double alp, double bet, double gam );
};


class hexagonal : public symops
{
   public:
      hexagonal( char const * name, int number, int nops, double const * begin, double const * end );
      virtual ~hexagonal() {}
      virtual void SetLatVec( double a, double b, double c, double alp, double bet, double gam );
};


class cubic : public symops
{
   public:
      cubic( char const * name, int number, int nops, double const * begin, double const * end );
      virtual ~cubic() {}
      virtual void SetLatVec( double a, double b, double c, double alp, double bet, double gam );
};


 class spcgrp_1 : public spcgrp::triclinic  { public: spcgrp_1(); virtual ~spcgrp_1() {} };
 class spcgrp_2 : public spcgrp::triclinic  { public: spcgrp_2(); virtual ~spcgrp_2() {} };
 class spcgrp_3 : public spcgrp::monoclinic  { public: spcgrp_3(); virtual ~spcgrp_3() {} };
 class spcgrp_4 : public spcgrp::monoclinic  { public: spcgrp_4(); virtual ~spcgrp_4() {} };
 class spcgrp_5 : public spcgrp::monoclinic  { public: spcgrp_5(); virtual ~spcgrp_5() {} };
 class spcgrp_6 : public spcgrp::monoclinic  { public: spcgrp_6(); virtual ~spcgrp_6() {} };
 class spcgrp_7 : public spcgrp::monoclinic  { public: spcgrp_7(); virtual ~spcgrp_7() {} };
 class spcgrp_8 : public spcgrp::monoclinic  { public: spcgrp_8(); virtual ~spcgrp_8() {} };
 class spcgrp_9 : public spcgrp::monoclinic  { public: spcgrp_9(); virtual ~spcgrp_9() {} };
 class spcgrp_10 : public spcgrp::monoclinic  { public: spcgrp_10(); virtual ~spcgrp_10() {} };
 class spcgrp_11 : public spcgrp::monoclinic  { public: spcgrp_11(); virtual ~spcgrp_11() {} };
 class spcgrp_12 : public spcgrp::monoclinic  { public: spcgrp_12(); virtual ~spcgrp_12() {} };
 class spcgrp_13 : public spcgrp::monoclinic  { public: spcgrp_13(); virtual ~spcgrp_13() {} };
 class spcgrp_14 : public spcgrp::monoclinic  { public: spcgrp_14(); virtual ~spcgrp_14() {} };
 class spcgrp_15 : public spcgrp::monoclinic  { public: spcgrp_15(); virtual ~spcgrp_15() {} };
 class spcgrp_16 : public spcgrp::orthorhombic  { public: spcgrp_16(); virtual ~spcgrp_16() {} };
 class spcgrp_17 : public spcgrp::orthorhombic  { public: spcgrp_17(); virtual ~spcgrp_17() {} };
 class spcgrp_18 : public spcgrp::orthorhombic  { public: spcgrp_18(); virtual ~spcgrp_18() {} };
 class spcgrp_19 : public spcgrp::orthorhombic  { public: spcgrp_19(); virtual ~spcgrp_19() {} };
 class spcgrp_20 : public spcgrp::orthorhombic  { public: spcgrp_20(); virtual ~spcgrp_20() {} };
 class spcgrp_21 : public spcgrp::orthorhombic  { public: spcgrp_21(); virtual ~spcgrp_21() {} };
 class spcgrp_22 : public spcgrp::orthorhombic  { public: spcgrp_22(); virtual ~spcgrp_22() {} };
 class spcgrp_23 : public spcgrp::orthorhombic  { public: spcgrp_23(); virtual ~spcgrp_23() {} };
 class spcgrp_24 : public spcgrp::orthorhombic  { public: spcgrp_24(); virtual ~spcgrp_24() {} };
 class spcgrp_25 : public spcgrp::orthorhombic  { public: spcgrp_25(); virtual ~spcgrp_25() {} };
 class spcgrp_26 : public spcgrp::orthorhombic  { public: spcgrp_26(); virtual ~spcgrp_26() {} };
 class spcgrp_27 : public spcgrp::orthorhombic  { public: spcgrp_27(); virtual ~spcgrp_27() {} };
 class spcgrp_28 : public spcgrp::orthorhombic  { public: spcgrp_28(); virtual ~spcgrp_28() {} };
 class spcgrp_29 : public spcgrp::orthorhombic  { public: spcgrp_29(); virtual ~spcgrp_29() {} };
 class spcgrp_30 : public spcgrp::orthorhombic  { public: spcgrp_30(); virtual ~spcgrp_30() {} };
 class spcgrp_31 : public spcgrp::orthorhombic  { public: spcgrp_31(); virtual ~spcgrp_31() {} };
 class spcgrp_32 : public spcgrp::orthorhombic  { public: spcgrp_32(); virtual ~spcgrp_32() {} };
 class spcgrp_33 : public spcgrp::orthorhombic  { public: spcgrp_33(); virtual ~spcgrp_33() {} };
 class spcgrp_34 : public spcgrp::orthorhombic  { public: spcgrp_34(); virtual ~spcgrp_34() {} };
 class spcgrp_35 : public spcgrp::orthorhombic  { public: spcgrp_35(); virtual ~spcgrp_35() {} };
 class spcgrp_36 : public spcgrp::orthorhombic  { public: spcgrp_36(); virtual ~spcgrp_36() {} };
 class spcgrp_37 : public spcgrp::orthorhombic  { public: spcgrp_37(); virtual ~spcgrp_37() {} };
 class spcgrp_38 : public spcgrp::orthorhombic  { public: spcgrp_38(); virtual ~spcgrp_38() {} };
 class spcgrp_39 : public spcgrp::orthorhombic  { public: spcgrp_39(); virtual ~spcgrp_39() {} };
 class spcgrp_40 : public spcgrp::orthorhombic  { public: spcgrp_40(); virtual ~spcgrp_40() {} };
 class spcgrp_41 : public spcgrp::orthorhombic  { public: spcgrp_41(); virtual ~spcgrp_41() {} };
 class spcgrp_42 : public spcgrp::orthorhombic  { public: spcgrp_42(); virtual ~spcgrp_42() {} };
 class spcgrp_43 : public spcgrp::orthorhombic  { public: spcgrp_43(); virtual ~spcgrp_43() {} };
 class spcgrp_44 : public spcgrp::orthorhombic  { public: spcgrp_44(); virtual ~spcgrp_44() {} };
 class spcgrp_45 : public spcgrp::orthorhombic  { public: spcgrp_45(); virtual ~spcgrp_45() {} };
 class spcgrp_46 : public spcgrp::orthorhombic  { public: spcgrp_46(); virtual ~spcgrp_46() {} };
 class spcgrp_47 : public spcgrp::orthorhombic  { public: spcgrp_47(); virtual ~spcgrp_47() {} };
 class spcgrp_48 : public spcgrp::orthorhombic  { public: spcgrp_48(); virtual ~spcgrp_48() {} };
 class spcgrp_49 : public spcgrp::orthorhombic  { public: spcgrp_49(); virtual ~spcgrp_49() {} };
 class spcgrp_50 : public spcgrp::orthorhombic  { public: spcgrp_50(); virtual ~spcgrp_50() {} };
 class spcgrp_51 : public spcgrp::orthorhombic  { public: spcgrp_51(); virtual ~spcgrp_51() {} };
 class spcgrp_52 : public spcgrp::orthorhombic  { public: spcgrp_52(); virtual ~spcgrp_52() {} };
 class spcgrp_53 : public spcgrp::orthorhombic  { public: spcgrp_53(); virtual ~spcgrp_53() {} };
 class spcgrp_54 : public spcgrp::orthorhombic  { public: spcgrp_54(); virtual ~spcgrp_54() {} };
 class spcgrp_55 : public spcgrp::orthorhombic  { public: spcgrp_55(); virtual ~spcgrp_55() {} };
 class spcgrp_56 : public spcgrp::orthorhombic  { public: spcgrp_56(); virtual ~spcgrp_56() {} };
 class spcgrp_57 : public spcgrp::orthorhombic  { public: spcgrp_57(); virtual ~spcgrp_57() {} };
 class spcgrp_58 : public spcgrp::orthorhombic  { public: spcgrp_58(); virtual ~spcgrp_58() {} };
 class spcgrp_59 : public spcgrp::orthorhombic  { public: spcgrp_59(); virtual ~spcgrp_59() {} };
 class spcgrp_60 : public spcgrp::orthorhombic  { public: spcgrp_60(); virtual ~spcgrp_60() {} };
 class spcgrp_61 : public spcgrp::orthorhombic  { public: spcgrp_61(); virtual ~spcgrp_61() {} };
 class spcgrp_62 : public spcgrp::orthorhombic  { public: spcgrp_62(); virtual ~spcgrp_62() {} };
 class spcgrp_63 : public spcgrp::orthorhombic  { public: spcgrp_63(); virtual ~spcgrp_63() {} };
 class spcgrp_64 : public spcgrp::orthorhombic  { public: spcgrp_64(); virtual ~spcgrp_64() {} };
 class spcgrp_65 : public spcgrp::orthorhombic  { public: spcgrp_65(); virtual ~spcgrp_65() {} };
 class spcgrp_66 : public spcgrp::orthorhombic  { public: spcgrp_66(); virtual ~spcgrp_66() {} };
 class spcgrp_67 : public spcgrp::orthorhombic  { public: spcgrp_67(); virtual ~spcgrp_67() {} };
 class spcgrp_68 : public spcgrp::orthorhombic  { public: spcgrp_68(); virtual ~spcgrp_68() {} };
 class spcgrp_69 : public spcgrp::orthorhombic  { public: spcgrp_69(); virtual ~spcgrp_69() {} };
 class spcgrp_70 : public spcgrp::orthorhombic  { public: spcgrp_70(); virtual ~spcgrp_70() {} };
 class spcgrp_71 : public spcgrp::orthorhombic  { public: spcgrp_71(); virtual ~spcgrp_71() {} };
 class spcgrp_72 : public spcgrp::orthorhombic  { public: spcgrp_72(); virtual ~spcgrp_72() {} };
 class spcgrp_73 : public spcgrp::orthorhombic  { public: spcgrp_73(); virtual ~spcgrp_73() {} };
 class spcgrp_74 : public spcgrp::orthorhombic  { public: spcgrp_74(); virtual ~spcgrp_74() {} };
 class spcgrp_75 : public spcgrp::tetragonal  { public: spcgrp_75(); virtual ~spcgrp_75() {} };
 class spcgrp_76 : public spcgrp::tetragonal  { public: spcgrp_76(); virtual ~spcgrp_76() {} };
 class spcgrp_77 : public spcgrp::tetragonal  { public: spcgrp_77(); virtual ~spcgrp_77() {} };
 class spcgrp_78 : public spcgrp::tetragonal  { public: spcgrp_78(); virtual ~spcgrp_78() {} };
 class spcgrp_79 : public spcgrp::tetragonal  { public: spcgrp_79(); virtual ~spcgrp_79() {} };
 class spcgrp_80 : public spcgrp::tetragonal  { public: spcgrp_80(); virtual ~spcgrp_80() {} };
 class spcgrp_81 : public spcgrp::tetragonal  { public: spcgrp_81(); virtual ~spcgrp_81() {} };
 class spcgrp_82 : public spcgrp::tetragonal  { public: spcgrp_82(); virtual ~spcgrp_82() {} };
 class spcgrp_83 : public spcgrp::tetragonal  { public: spcgrp_83(); virtual ~spcgrp_83() {} };
 class spcgrp_84 : public spcgrp::tetragonal  { public: spcgrp_84(); virtual ~spcgrp_84() {} };
 class spcgrp_85 : public spcgrp::tetragonal  { public: spcgrp_85(); virtual ~spcgrp_85() {} };
 class spcgrp_86 : public spcgrp::tetragonal  { public: spcgrp_86(); virtual ~spcgrp_86() {} };
 class spcgrp_87 : public spcgrp::tetragonal  { public: spcgrp_87(); virtual ~spcgrp_87() {} };
 class spcgrp_88 : public spcgrp::tetragonal  { public: spcgrp_88(); virtual ~spcgrp_88() {} };
 class spcgrp_89 : public spcgrp::tetragonal  { public: spcgrp_89(); virtual ~spcgrp_89() {} };
 class spcgrp_90 : public spcgrp::tetragonal  { public: spcgrp_90(); virtual ~spcgrp_90() {} };
 class spcgrp_91 : public spcgrp::tetragonal  { public: spcgrp_91(); virtual ~spcgrp_91() {} };
 class spcgrp_92 : public spcgrp::tetragonal  { public: spcgrp_92(); virtual ~spcgrp_92() {} };
 class spcgrp_93 : public spcgrp::tetragonal  { public: spcgrp_93(); virtual ~spcgrp_93() {} };
 class spcgrp_94 : public spcgrp::tetragonal  { public: spcgrp_94(); virtual ~spcgrp_94() {} };
 class spcgrp_95 : public spcgrp::tetragonal  { public: spcgrp_95(); virtual ~spcgrp_95() {} };
 class spcgrp_96 : public spcgrp::tetragonal  { public: spcgrp_96(); virtual ~spcgrp_96() {} };
 class spcgrp_97 : public spcgrp::tetragonal  { public: spcgrp_97(); virtual ~spcgrp_97() {} };
 class spcgrp_98 : public spcgrp::tetragonal  { public: spcgrp_98(); virtual ~spcgrp_98() {} };
 class spcgrp_99 : public spcgrp::tetragonal  { public: spcgrp_99(); virtual ~spcgrp_99() {} };
 class spcgrp_100 : public spcgrp::tetragonal  { public: spcgrp_100(); virtual ~spcgrp_100() {} };
 class spcgrp_101 : public spcgrp::tetragonal  { public: spcgrp_101(); virtual ~spcgrp_101() {} };
 class spcgrp_102 : public spcgrp::tetragonal  { public: spcgrp_102(); virtual ~spcgrp_102() {} };
 class spcgrp_103 : public spcgrp::tetragonal  { public: spcgrp_103(); virtual ~spcgrp_103() {} };
 class spcgrp_104 : public spcgrp::tetragonal  { public: spcgrp_104(); virtual ~spcgrp_104() {} };
 class spcgrp_105 : public spcgrp::tetragonal  { public: spcgrp_105(); virtual ~spcgrp_105() {} };
 class spcgrp_106 : public spcgrp::tetragonal  { public: spcgrp_106(); virtual ~spcgrp_106() {} };
 class spcgrp_107 : public spcgrp::tetragonal  { public: spcgrp_107(); virtual ~spcgrp_107() {} };
 class spcgrp_108 : public spcgrp::tetragonal  { public: spcgrp_108(); virtual ~spcgrp_108() {} };
 class spcgrp_109 : public spcgrp::tetragonal  { public: spcgrp_109(); virtual ~spcgrp_109() {} };
 class spcgrp_110 : public spcgrp::tetragonal  { public: spcgrp_110(); virtual ~spcgrp_110() {} };
 class spcgrp_111 : public spcgrp::tetragonal  { public: spcgrp_111(); virtual ~spcgrp_111() {} };
 class spcgrp_112 : public spcgrp::tetragonal  { public: spcgrp_112(); virtual ~spcgrp_112() {} };
 class spcgrp_113 : public spcgrp::tetragonal  { public: spcgrp_113(); virtual ~spcgrp_113() {} };
 class spcgrp_114 : public spcgrp::tetragonal  { public: spcgrp_114(); virtual ~spcgrp_114() {} };
 class spcgrp_115 : public spcgrp::tetragonal  { public: spcgrp_115(); virtual ~spcgrp_115() {} };
 class spcgrp_116 : public spcgrp::tetragonal  { public: spcgrp_116(); virtual ~spcgrp_116() {} };
 class spcgrp_117 : public spcgrp::tetragonal  { public: spcgrp_117(); virtual ~spcgrp_117() {} };
 class spcgrp_118 : public spcgrp::tetragonal  { public: spcgrp_118(); virtual ~spcgrp_118() {} };
 class spcgrp_119 : public spcgrp::tetragonal  { public: spcgrp_119(); virtual ~spcgrp_119() {} };
 class spcgrp_120 : public spcgrp::tetragonal  { public: spcgrp_120(); virtual ~spcgrp_120() {} };
 class spcgrp_121 : public spcgrp::tetragonal  { public: spcgrp_121(); virtual ~spcgrp_121() {} };
 class spcgrp_122 : public spcgrp::tetragonal  { public: spcgrp_122(); virtual ~spcgrp_122() {} };
 class spcgrp_123 : public spcgrp::tetragonal  { public: spcgrp_123(); virtual ~spcgrp_123() {} };
 class spcgrp_124 : public spcgrp::tetragonal  { public: spcgrp_124(); virtual ~spcgrp_124() {} };
 class spcgrp_125 : public spcgrp::tetragonal  { public: spcgrp_125(); virtual ~spcgrp_125() {} };
 class spcgrp_126 : public spcgrp::tetragonal  { public: spcgrp_126(); virtual ~spcgrp_126() {} };
 class spcgrp_127 : public spcgrp::tetragonal  { public: spcgrp_127(); virtual ~spcgrp_127() {} };
 class spcgrp_128 : public spcgrp::tetragonal  { public: spcgrp_128(); virtual ~spcgrp_128() {} };
 class spcgrp_129 : public spcgrp::tetragonal  { public: spcgrp_129(); virtual ~spcgrp_129() {} };
 class spcgrp_130 : public spcgrp::tetragonal  { public: spcgrp_130(); virtual ~spcgrp_130() {} };
 class spcgrp_131 : public spcgrp::tetragonal  { public: spcgrp_131(); virtual ~spcgrp_131() {} };
 class spcgrp_132 : public spcgrp::tetragonal  { public: spcgrp_132(); virtual ~spcgrp_132() {} };
 class spcgrp_133 : public spcgrp::tetragonal  { public: spcgrp_133(); virtual ~spcgrp_133() {} };
 class spcgrp_134 : public spcgrp::tetragonal  { public: spcgrp_134(); virtual ~spcgrp_134() {} };
 class spcgrp_135 : public spcgrp::tetragonal  { public: spcgrp_135(); virtual ~spcgrp_135() {} };
 class spcgrp_136 : public spcgrp::tetragonal  { public: spcgrp_136(); virtual ~spcgrp_136() {} };
 class spcgrp_137 : public spcgrp::tetragonal  { public: spcgrp_137(); virtual ~spcgrp_137() {} };
 class spcgrp_138 : public spcgrp::tetragonal  { public: spcgrp_138(); virtual ~spcgrp_138() {} };
 class spcgrp_139 : public spcgrp::tetragonal  { public: spcgrp_139(); virtual ~spcgrp_139() {} };
 class spcgrp_140 : public spcgrp::tetragonal  { public: spcgrp_140(); virtual ~spcgrp_140() {} };
 class spcgrp_141 : public spcgrp::tetragonal  { public: spcgrp_141(); virtual ~spcgrp_141() {} };
 class spcgrp_142 : public spcgrp::tetragonal  { public: spcgrp_142(); virtual ~spcgrp_142() {} };
 class spcgrp_143 : public spcgrp::trigonal  { public: spcgrp_143(); virtual ~spcgrp_143() {} };
 class spcgrp_144 : public spcgrp::trigonal  { public: spcgrp_144(); virtual ~spcgrp_144() {} };
 class spcgrp_145 : public spcgrp::trigonal  { public: spcgrp_145(); virtual ~spcgrp_145() {} };
 class spcgrp_146 : public spcgrp::trigonal  { public: spcgrp_146(); virtual ~spcgrp_146() {} };
 class spcgrp_147 : public spcgrp::trigonal  { public: spcgrp_147(); virtual ~spcgrp_147() {} };
 class spcgrp_148 : public spcgrp::trigonal  { public: spcgrp_148(); virtual ~spcgrp_148() {} };
 class spcgrp_149 : public spcgrp::trigonal  { public: spcgrp_149(); virtual ~spcgrp_149() {} };
 class spcgrp_150 : public spcgrp::trigonal  { public: spcgrp_150(); virtual ~spcgrp_150() {} };
 class spcgrp_151 : public spcgrp::trigonal  { public: spcgrp_151(); virtual ~spcgrp_151() {} };
 class spcgrp_152 : public spcgrp::trigonal  { public: spcgrp_152(); virtual ~spcgrp_152() {} };
 class spcgrp_153 : public spcgrp::trigonal  { public: spcgrp_153(); virtual ~spcgrp_153() {} };
 class spcgrp_154 : public spcgrp::trigonal  { public: spcgrp_154(); virtual ~spcgrp_154() {} };
 class spcgrp_155 : public spcgrp::trigonal  { public: spcgrp_155(); virtual ~spcgrp_155() {} };
 class spcgrp_156 : public spcgrp::trigonal  { public: spcgrp_156(); virtual ~spcgrp_156() {} };
 class spcgrp_157 : public spcgrp::trigonal  { public: spcgrp_157(); virtual ~spcgrp_157() {} };
 class spcgrp_158 : public spcgrp::trigonal  { public: spcgrp_158(); virtual ~spcgrp_158() {} };
 class spcgrp_159 : public spcgrp::trigonal  { public: spcgrp_159(); virtual ~spcgrp_159() {} };
 class spcgrp_160 : public spcgrp::trigonal  { public: spcgrp_160(); virtual ~spcgrp_160() {} };
 class spcgrp_161 : public spcgrp::trigonal  { public: spcgrp_161(); virtual ~spcgrp_161() {} };
 class spcgrp_162 : public spcgrp::trigonal  { public: spcgrp_162(); virtual ~spcgrp_162() {} };
 class spcgrp_163 : public spcgrp::trigonal  { public: spcgrp_163(); virtual ~spcgrp_163() {} };
 class spcgrp_164 : public spcgrp::trigonal  { public: spcgrp_164(); virtual ~spcgrp_164() {} };
 class spcgrp_165 : public spcgrp::trigonal  { public: spcgrp_165(); virtual ~spcgrp_165() {} };
 class spcgrp_166 : public spcgrp::trigonal  { public: spcgrp_166(); virtual ~spcgrp_166() {} };
 class spcgrp_167 : public spcgrp::trigonal  { public: spcgrp_167(); virtual ~spcgrp_167() {} };
 class spcgrp_168 : public spcgrp::hexagonal  { public: spcgrp_168(); virtual ~spcgrp_168() {} };
 class spcgrp_169 : public spcgrp::hexagonal  { public: spcgrp_169(); virtual ~spcgrp_169() {} };
 class spcgrp_170 : public spcgrp::hexagonal  { public: spcgrp_170(); virtual ~spcgrp_170() {} };
 class spcgrp_171 : public spcgrp::hexagonal  { public: spcgrp_171(); virtual ~spcgrp_171() {} };
 class spcgrp_172 : public spcgrp::hexagonal  { public: spcgrp_172(); virtual ~spcgrp_172() {} };
 class spcgrp_173 : public spcgrp::hexagonal  { public: spcgrp_173(); virtual ~spcgrp_173() {} };
 class spcgrp_174 : public spcgrp::hexagonal  { public: spcgrp_174(); virtual ~spcgrp_174() {} };
 class spcgrp_175 : public spcgrp::hexagonal  { public: spcgrp_175(); virtual ~spcgrp_175() {} };
 class spcgrp_176 : public spcgrp::hexagonal  { public: spcgrp_176(); virtual ~spcgrp_176() {} };
 class spcgrp_177 : public spcgrp::hexagonal  { public: spcgrp_177(); virtual ~spcgrp_177() {} };
 class spcgrp_178 : public spcgrp::hexagonal  { public: spcgrp_178(); virtual ~spcgrp_178() {} };
 class spcgrp_179 : public spcgrp::hexagonal  { public: spcgrp_179(); virtual ~spcgrp_179() {} };
 class spcgrp_180 : public spcgrp::hexagonal  { public: spcgrp_180(); virtual ~spcgrp_180() {} };
 class spcgrp_181 : public spcgrp::hexagonal  { public: spcgrp_181(); virtual ~spcgrp_181() {} };
 class spcgrp_182 : public spcgrp::hexagonal  { public: spcgrp_182(); virtual ~spcgrp_182() {} };
 class spcgrp_183 : public spcgrp::hexagonal  { public: spcgrp_183(); virtual ~spcgrp_183() {} };
 class spcgrp_184 : public spcgrp::hexagonal  { public: spcgrp_184(); virtual ~spcgrp_184() {} };
 class spcgrp_185 : public spcgrp::hexagonal  { public: spcgrp_185(); virtual ~spcgrp_185() {} };
 class spcgrp_186 : public spcgrp::hexagonal  { public: spcgrp_186(); virtual ~spcgrp_186() {} };
 class spcgrp_187 : public spcgrp::hexagonal  { public: spcgrp_187(); virtual ~spcgrp_187() {} };
 class spcgrp_188 : public spcgrp::hexagonal  { public: spcgrp_188(); virtual ~spcgrp_188() {} };
 class spcgrp_189 : public spcgrp::hexagonal  { public: spcgrp_189(); virtual ~spcgrp_189() {} };
 class spcgrp_190 : public spcgrp::hexagonal  { public: spcgrp_190(); virtual ~spcgrp_190() {} };
 class spcgrp_191 : public spcgrp::hexagonal  { public: spcgrp_191(); virtual ~spcgrp_191() {} };
 class spcgrp_192 : public spcgrp::hexagonal  { public: spcgrp_192(); virtual ~spcgrp_192() {} };
 class spcgrp_193 : public spcgrp::hexagonal  { public: spcgrp_193(); virtual ~spcgrp_193() {} };
 class spcgrp_194 : public spcgrp::hexagonal  { public: spcgrp_194(); virtual ~spcgrp_194() {} };
 class spcgrp_195 : public spcgrp::cubic  { public: spcgrp_195(); virtual ~spcgrp_195() {} };
 class spcgrp_196 : public spcgrp::cubic  { public: spcgrp_196(); virtual ~spcgrp_196() {} };
 class spcgrp_197 : public spcgrp::cubic  { public: spcgrp_197(); virtual ~spcgrp_197() {} };
 class spcgrp_198 : public spcgrp::cubic  { public: spcgrp_198(); virtual ~spcgrp_198() {} };
 class spcgrp_199 : public spcgrp::cubic  { public: spcgrp_199(); virtual ~spcgrp_199() {} };
 class spcgrp_200 : public spcgrp::cubic  { public: spcgrp_200(); virtual ~spcgrp_200() {} };
 class spcgrp_201 : public spcgrp::cubic  { public: spcgrp_201(); virtual ~spcgrp_201() {} };
 class spcgrp_202 : public spcgrp::cubic  { public: spcgrp_202(); virtual ~spcgrp_202() {} };
 class spcgrp_203 : public spcgrp::cubic  { public: spcgrp_203(); virtual ~spcgrp_203() {} };
 class spcgrp_204 : public spcgrp::cubic  { public: spcgrp_204(); virtual ~spcgrp_204() {} };
 class spcgrp_205 : public spcgrp::cubic  { public: spcgrp_205(); virtual ~spcgrp_205() {} };
 class spcgrp_206 : public spcgrp::cubic  { public: spcgrp_206(); virtual ~spcgrp_206() {} };
 class spcgrp_207 : public spcgrp::cubic  { public: spcgrp_207(); virtual ~spcgrp_207() {} };
 class spcgrp_208 : public spcgrp::cubic  { public: spcgrp_208(); virtual ~spcgrp_208() {} };
 class spcgrp_209 : public spcgrp::cubic  { public: spcgrp_209(); virtual ~spcgrp_209() {} };
 class spcgrp_210 : public spcgrp::cubic  { public: spcgrp_210(); virtual ~spcgrp_210() {} };
 class spcgrp_211 : public spcgrp::cubic  { public: spcgrp_211(); virtual ~spcgrp_211() {} };
 class spcgrp_212 : public spcgrp::cubic  { public: spcgrp_212(); virtual ~spcgrp_212() {} };
 class spcgrp_213 : public spcgrp::cubic  { public: spcgrp_213(); virtual ~spcgrp_213() {} };
 class spcgrp_214 : public spcgrp::cubic  { public: spcgrp_214(); virtual ~spcgrp_214() {} };
 class spcgrp_215 : public spcgrp::cubic  { public: spcgrp_215(); virtual ~spcgrp_215() {} };
 class spcgrp_216 : public spcgrp::cubic  { public: spcgrp_216(); virtual ~spcgrp_216() {} };
 class spcgrp_217 : public spcgrp::cubic  { public: spcgrp_217(); virtual ~spcgrp_217() {} };
 class spcgrp_218 : public spcgrp::cubic  { public: spcgrp_218(); virtual ~spcgrp_218() {} };
 class spcgrp_219 : public spcgrp::cubic  { public: spcgrp_219(); virtual ~spcgrp_219() {} };
 class spcgrp_220 : public spcgrp::cubic  { public: spcgrp_220(); virtual ~spcgrp_220() {} };
 class spcgrp_221 : public spcgrp::cubic  { public: spcgrp_221(); virtual ~spcgrp_221() {} };
 class spcgrp_222 : public spcgrp::cubic  { public: spcgrp_222(); virtual ~spcgrp_222() {} };
 class spcgrp_223 : public spcgrp::cubic  { public: spcgrp_223(); virtual ~spcgrp_223() {} };
 class spcgrp_224 : public spcgrp::cubic  { public: spcgrp_224(); virtual ~spcgrp_224() {} };
 class spcgrp_225 : public spcgrp::cubic  { public: spcgrp_225(); virtual ~spcgrp_225() {} };
 class spcgrp_226 : public spcgrp::cubic  { public: spcgrp_226(); virtual ~spcgrp_226() {} };
 class spcgrp_227 : public spcgrp::cubic  { public: spcgrp_227(); virtual ~spcgrp_227() {} };
 class spcgrp_228 : public spcgrp::cubic  { public: spcgrp_228(); virtual ~spcgrp_228() {} };
 class spcgrp_229 : public spcgrp::cubic  { public: spcgrp_229(); virtual ~spcgrp_229() {} };
 class spcgrp_230 : public spcgrp::cubic  { public: spcgrp_230(); virtual ~spcgrp_230() {} };

}

#endif
