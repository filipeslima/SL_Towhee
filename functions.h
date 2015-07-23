c     ******************************************************************
c     * functions.h                                                    *
c     * contains all of the type declarations for the functions in     *
c     * towhee                                                         *
c     *                                                                *
c     * originally written 11-12-2002 by M.G. Martin                   *
c     * last modified 08-15-2013 by M.G. Martin                        *
c     ******************************************************************
#ifdef FUNCTION_ARCCOS
      double precision twh_arccos
#undef FUNCTION_ARCCOS
#endif
c
#ifdef FUNCTION_CHECK_COMMENT
      logical twh_check_comment
#undef FUNCTION_CHECK_COMMENT
#endif
c
#ifdef FUNCTION_CHECK_LABEL
      logical twh_check_label
#undef FUNCTION_CHECK_LABEL
#endif
c
#ifdef FUNCTION_CMP_EQ
#ifdef SAFE_COMPARE
      logical twh_cmp_eq
#else
#define twh_cmp_eq(A,B) ((A).eq.(B))
#endif
#undef FUNCTION_CMP_EQ
#endif
c
#ifdef FUNCTION_CMP_GT
#ifdef SAFE_COMPARE
      logical twh_cmp_gt
#else
#define twh_cmp_gt(A,B) ((A).gt.(B))
#endif
#undef FUNCTION_CMP_GT
#endif
c
#ifdef FUNCTION_CMP_LT
#ifdef SAFE_COMPARE
      logical twh_cmp_lt
#else
#define twh_cmp_lt(A,B) ((A).lt.(B))
#endif
#undef FUNCTION_CMP_LT
#endif
c
#ifdef FUNCTION_DERF
      double precision twh_derf 
#undef FUNCTION_DERF
#endif
c
#ifdef FUNCTION_DERFC
      double precision twh_derfc
#undef FUNCTION_DERFC
#endif
c
#ifdef FUNCTION_DISTANCE
      double precision twh_distance
#undef FUNCTION_DISTANCE
#endif
c
#ifdef FUNCTION_DOTPRODUCT
      double precision twh_dotproduct
#undef FUNCTION_DOTPRODUCT
#endif
c
#ifdef FUNCTION_EWALD_CORRECT
      double precision twh_ewald_correct
#undef FUNCTION_EWALD_CORRECT
#endif
c
#ifdef FUNCTION_EWALD_SELF
      double precision twh_ewald_self
#undef FUNCTION_EWALD_SELF
#endif
c
#ifdef FUNCTION_EXPON
      double precision twh_expon 
#undef FUNCTION_EXPON
#endif
c
#ifdef FUNCTION_EXTRACTDENS
      double precision twh_extractdens 
#undef FUNCTION_EXTRACTDENS
#endif
c
#ifdef FUNCTION_EXTRACTEMBED
      double precision twh_extractembed
#undef FUNCTION_EXTRACTEMBED
#endif
c
#ifdef FUNCTION_EXTRACTPAIR
      double precision twh_extractpair 
#undef FUNCTION_EXTRACTPAIR
#endif
c
#ifdef FUNCTION_FEBIAS
      double precision twh_febias 
#undef FUNCTION_FEBIAS
#endif
c
#ifdef FUNCTION_GAUSSPROB
      double precision twh_gaussprob 
#undef FUNCTION_GAUSSPROB
#endif
c
#ifdef FUNCTION_GAUSSIAN
      double precision twh_gaussian 
#undef FUNCTION_GAUSSIAN
#endif
c
#ifdef FUNCTION_GETATOMNUM
      integer twh_getatomnum
#undef FUNCTION_GETATOMNUM
#endif
c
#ifdef FUNCTION_GETNBTYPE
      integer twh_getnbtype
#undef FUNCTION_GETNBTYPE
#endif
c
#ifdef FUNCTION_GETSIGN
      double precision twh_getsign 
#undef FUNCTION_GETSIGN
#endif
c
#ifdef FUNCTION_GET_AACOEFF
      double precision twh_get_aacoeff
#undef FUNCTION_GET_AACOEFF
#endif
c
#ifdef FUNCTION_GET_AAFF
      character*(FFNAMELEN) twh_get_aaff
#undef FUNCTION_GET_AAFF
#endif
c
#ifdef FUNCTION_GET_AANAMES
      character*(FFNAMELEN) twh_get_aanames
#undef FUNCTION_GET_AANAMES
#endif
c
#ifdef FUNCTION_GET_AASTYLE
      integer twh_get_aastyle
#undef FUNCTION_GET_AASTYLE
#endif
c
#ifdef FUNCTION_GET_ANGLESTYLE
      integer twh_get_anglestyle
#undef FUNCTION_GET_ANGLESTYLE
#endif
c
#ifdef FUNCTION_GET_ATOMLIST_CORD
      integer twh_get_atomlist_cord
#undef FUNCTION_GET_ATOMLIST_CORD
#endif
c
#ifdef FUNCTION_GET_BENCOEFF
      double precision twh_get_bencoeff
#undef FUNCTION_GET_BENCOEFF
#endif
c
#ifdef FUNCTION_GET_BENDFF
      character*(FFNAMELEN) twh_get_bendff
#undef FUNCTION_GET_BENDFF
#endif
c
#ifdef FUNCTION_GET_BENDNAMES
      character*(FFNAMELEN) twh_get_bendnames
#undef FUNCTION_GET_BENDNAMES
#endif
c
#ifdef FUNCTION_GET_BENPREFACT
      double precision twh_get_benprefact
#undef FUNCTION_GET_BENPREFACT
#endif
c
#ifdef FUNCTION_GET_BONDPATT
      character*(5) twh_get_bondpatt
#undef FUNCTION_GET_BONDPATT
#endif
c
#ifdef FUNCTION_GET_CANAME
      character*(FFNAMELEN) twh_get_caname
#undef FUNCTION_GET_CANAME
#endif
c
#ifdef FUNCTION_GET_CBNAME
      character*(FFNAMELEN) twh_get_cbname
#undef FUNCTION_GET_CBNAME
#endif
c
#ifdef FUNCTION_GET_CLASSICAL_MIXRULE
      character*(30) twh_get_classical_mixrule
#undef FUNCTION_GET_CLASSICAL_MIXRULE
#endif
c
#ifdef FUNCTION_GET_CLASSICAL_POTENTIAL
      character*(30) twh_get_classical_potential
#undef FUNCTION_GET_CLASSICAL_POTENTIAL
#endif
c
#ifdef FUNCTION_GET_CTNAME
      character*(FFNAMELEN) twh_get_ctname
#undef FUNCTION_GET_CTNAME
#endif
c
#ifdef FUNCTION_GET_DERIVATIVE
      double precision scp_get_derivative
#undef FUNCTION_GET_DERIVATIVE
#endif
c
#ifdef FUNCTION_GET_ELEMENTNAME
      character*2 twh_get_elementname
#undef FUNCTION_GET_ELEMENTNAME
#endif
c
#ifdef FUNCTION_GET_FOREIGN_LAMBDA
      double precision scp_get_foreign_lambda
#undef FUNCTION_GET_FOREIGN_LAMBDA
#endif
c
#ifdef FUNCTION_GET_GROWCOUNT
      integer twh_get_growcount
#undef FUNCTION_GET_GROWCOUNT
#endif
c
#ifdef FUNCTION_GET_GROWVALIDCOUNT
      integer twh_get_growvalidcount
#undef FUNCTION_GET_GROWVALIDCOUNT
#endif
c
#ifdef FUNCTION_GET_GROWVALIDLIST
      integer twh_get_growvalidlist
#undef FUNCTION_GET_GROWVALIDLIST
#endif
c
#ifdef FUNCTION_GET_HBONDCOEFF
      double precision twh_get_hbondcoeff
#undef FUNCTION_GET_HBONDCOEFF
#endif
c
#ifdef FUNCTION_GET_HBONDNAMES
      character*(FFNAMELEN) twh_get_hbondnames
#undef FUNCTION_GET_HBONDNAMES
#endif
c
#ifdef FUNCTION_GET_IBTBEN1
      integer twh_get_ibtben1
#undef FUNCTION_GET_IBTBEN1
#endif
c
#ifdef FUNCTION_GET_IBTBEN2
      integer twh_get_ibtben2
#undef FUNCTION_GET_IBTBEN2
#endif
c
#ifdef FUNCTION_GET_IJAA0
      integer twh_get_ijaa0
#undef FUNCTION_GET_IJAA0
#endif
c
#ifdef FUNCTION_GET_IJAA1
      integer twh_get_ijaa1
#undef FUNCTION_GET_IJAA1
#endif
c
#ifdef FUNCTION_GET_IJAA2
      integer twh_get_ijaa2
#undef FUNCTION_GET_IJAA2
#endif
c
#ifdef FUNCTION_GET_IJBEN2
      integer twh_get_ijben2
#undef FUNCTION_GET_IJBEN2
#endif
c
#ifdef FUNCTION_GET_IJBEN3
      integer twh_get_ijben3
#undef FUNCTION_GET_IJBEN3
#endif
c
#ifdef FUNCTION_GET_IJBNBN1
      integer twh_get_ijbnbn1
#undef FUNCTION_GET_IJBNBN1
#endif
c
#ifdef FUNCTION_GET_IJBNBN2
      integer twh_get_ijbnbn2
#undef FUNCTION_GET_IJBNBN2
#endif
c
#ifdef FUNCTION_GET_IJIMPROP2
      integer twh_get_ijimprop2
#undef FUNCTION_GET_IJIMPROP2
#endif
c
#ifdef FUNCTION_GET_IJIMPROP3
      integer twh_get_ijimprop3
#undef FUNCTION_GET_IJIMPROP3
#endif
c
#ifdef FUNCTION_GET_IJIMPROP4
      integer twh_get_ijimprop4
#undef FUNCTION_GET_IJIMPROP4
#endif
c
#ifdef FUNCTION_GET_IJOF5
      integer twh_get_ijof5
#undef FUNCTION_GET_IJOF5
#endif
c
#ifdef FUNCTION_GET_IJTOR2
      integer twh_get_ijtor2
#undef FUNCTION_GET_IJTOR2
#endif
c
#ifdef FUNCTION_GET_IJTOR3
      integer twh_get_ijtor3
#undef FUNCTION_GET_IJTOR3
#endif
c
#ifdef FUNCTION_GET_IJTOR4
      integer twh_get_ijtor4
#undef FUNCTION_GET_IJTOR4
#endif
c
#ifdef FUNCTION_GET_IJVIB
      integer twh_get_ijvib
#undef FUNCTION_GET_IJVIB
#endif
c
#ifdef FUNCTION_GET_IMPCOEFF
      double precision twh_get_impcoeff
#undef FUNCTION_GET_IMPCOEFF
#endif
c
#ifdef FUNCTION_GET_IMPFF
      character*(FFNAMELEN) twh_get_impff
#undef FUNCTION_GET_IMPFF
#endif
c
#ifdef FUNCTION_GET_IMPFORM
      integer twh_get_impform
#undef FUNCTION_GET_IMPFORM
#endif
c
#ifdef FUNCTION_GET_IMPNAMES
      character*(FFNAMELEN) twh_get_impnames
#undef FUNCTION_GET_IMPNAMES
#endif
c
#ifdef FUNCTION_GET_IMPSTYLE
      integer twh_get_impstyle
#undef FUNCTION_GET_IMPSTYLE
#endif
c
#ifdef FUNCTION_GET_INAA
      integer twh_get_inaa
#undef FUNCTION_GET_INAA
#endif
c
#ifdef FUNCTION_GET_INBEN
      integer twh_get_inben
#undef FUNCTION_GET_INBEN
#endif
c
#ifdef FUNCTION_GET_INBNBN
      integer twh_get_inbnbn
#undef FUNCTION_GET_INBNBN
#endif
c
#ifdef FUNCTION_GET_INIMPROP
      integer twh_get_inimprop
#undef FUNCTION_GET_INIMPROP
#endif
c
#ifdef FUNCTION_GET_INOF
      integer twh_get_inof
#undef FUNCTION_GET_INOF
#endif
c
#ifdef FUNCTION_GET_INTOR
      integer twh_get_intor
#undef FUNCTION_GET_INTOR
#endif
c
#ifdef FUNCTION_GET_INVIB
      integer twh_get_invib
#undef FUNCTION_GET_INVIB
#endif
c
#ifdef FUNCTION_GET_ITAA
      integer twh_get_itaa
#undef FUNCTION_GET_ITAA
#endif
c
#ifdef FUNCTION_GET_ITBEN
      integer twh_get_itben
#undef FUNCTION_GET_ITBEN
#endif
c
#ifdef FUNCTION_GET_ITBNBN
      integer twh_get_itbnbn
#undef FUNCTION_GET_ITBNBN
#endif
c
#ifdef FUNCTION_GET_ITIMPROP
      integer twh_get_itimprop
#undef FUNCTION_GET_ITIMPROP
#endif
c
#ifdef FUNCTION_GET_ITOF
      integer twh_get_itof
#undef FUNCTION_GET_ITOF
#endif
c
#ifdef FUNCTION_GET_ITSCALE
      double precision twh_get_itscale
#undef FUNCTION_GET_ITSCALE
#endif
c
#ifdef FUNCTION_GET_ITTOR
      integer twh_get_ittor
#undef FUNCTION_GET_ITTOR
#endif
c
#ifdef FUNCTION_GET_ITVIB
      integer twh_get_itvib
#undef FUNCTION_GET_ITVIB
#endif
c
#ifdef FUNCTION_GET_LAAHERE
      logical twh_get_laahere
#undef FUNCTION_GET_LAAHERE
#endif
c
#ifdef FUNCTION_GET_LBENDHERE
      logical twh_get_lbendhere
#undef FUNCTION_GET_LBENDHERE
#endif
c
#ifdef FUNCTION_GET_LBONANG
      logical twh_get_lbonang
#undef FUNCTION_GET_LBONANG
#endif
c
#ifdef FUNCTION_GET_LBONBON
      logical twh_get_lbonbon
#undef FUNCTION_GET_LBONBON
#endif
c
#ifdef FUNCTION_GET_LEXIST
      logical twh_get_lexist
#undef FUNCTION_GET_LEXIST
#endif
c
#ifdef FUNCTION_GET_LHERE
      logical twh_get_lhere
#undef FUNCTION_GET_LHERE
#endif
c
#ifdef FUNCTION_GET_LIMPHERE
      logical twh_get_limphere
#undef FUNCTION_GET_LIMPHERE
#endif
c
#ifdef FUNCTION_GET_LOFHERE
      logical twh_get_lofhere
#undef FUNCTION_GET_LOFHERE
#endif
c
#ifdef FUNCTION_GET_LOFTOR
      logical twh_get_loftor
#undef FUNCTION_GET_LOFTOR
#endif
c
#ifdef FUNCTION_GET_LTORHERE
      logical twh_get_ltorhere
#undef FUNCTION_GET_LTORHERE
#endif
c
#ifdef FUNCTION_GET_LVIBHERE
      logical twh_get_lvibhere
#undef FUNCTION_GET_LVIBHERE
#endif
c
#ifdef FUNCTION_GET_MASS
      double precision twh_get_mass
#undef FUNCTION_GET_MASS
#endif
c
#ifdef FUNCTION_GET_NAASAME
      integer twh_get_naasame
#undef FUNCTION_GET_NAASAME
#endif
c
#ifdef FUNCTION_GET_NATIVE_LAMBDA
      double precision scp_get_native_lambda
#undef FUNCTION_GET_NATIVE_LAMBDA
#endif
c
#ifdef FUNCTION_GET_NBCOEFF
      double precision twh_get_nbcoeff
#undef FUNCTION_GET_NBCOEFF
#endif
c
#ifdef FUNCTION_GET_NBOXI
      integer twh_get_nboxi
#undef FUNCTION_GET_NBOXI
#endif
c
#ifdef FUNCTION_GET_NCOEFF
      integer twh_get_ncoeff
#undef FUNCTION_GET_NCOEFF
#endif
c
#ifdef FUNCTION_GET_NBFF
      character*(FFNAMELEN) twh_get_nbff
#undef FUNCTION_GET_NBFF
#endif
c
#ifdef FUNCTION_GET_NBNAME
      character*(FFNAMELEN) twh_get_nbname
#undef FUNCTION_GET_NBNAME
#endif
c
#ifdef FUNCTION_GET_NBSAME
      integer twh_get_nbsame
#undef FUNCTION_GET_NBSAME
#endif
c
#ifdef FUNCTION_GET_NCMT
      integer twh_get_ncmt
#undef FUNCTION_GET_NCMT
#endif
c
#ifdef FUNCTION_GET_NIMPSAME
      integer twh_get_nimpsame
#undef FUNCTION_GET_NIMPSAME
#endif
c
#ifdef FUNCTION_GET_NTSAME
      integer twh_get_ntsame
#undef FUNCTION_GET_NTSAME
#endif
c
#ifdef FUNCTION_GET_NTORLOOP
      integer twh_get_ntorloop
#undef FUNCTION_GET_NTORLOOP
#endif
c
#ifdef FUNCTION_GET_NTYPE
      integer twh_get_ntype
#undef FUNCTION_GET_NTYPE
#endif
c
#ifdef FUNCTION_GET_NVSAME
      integer twh_get_nvsame
#undef FUNCTION_GET_NVSAME
#endif
c
#ifdef FUNCTION_GET_OFCOEFF
      double precision twh_get_ofcoeff
#undef FUNCTION_GET_OFCOEFF
#endif
c
#ifdef FUNCTION_GET_OFFF
      character*(FFNAMELEN) twh_get_offf
#undef FUNCTION_GET_OFFF
#endif
c
#ifdef FUNCTION_GET_OFNAMES
      character*(FFNAMELEN) twh_get_ofnames
#undef FUNCTION_GET_OFNAMES
#endif
c
#ifdef FUNCTION_GET_OFSTYLE
      integer twh_get_ofstyle
#undef FUNCTION_GET_OFSTYLE
#endif
c
#ifdef FUNCTION_GET_PM1BOXCOMSWITCH
      double precision twh_get_pm1boxcomswitch
#undef FUNCTION_GET_PM1BOXCOMSWITCH
#endif
c
#ifdef FUNCTION_GET_PMVLPR
      double precision twh_get_pmvlpr
#undef FUNCTION_GET_PMVLPR
#endif
c
#ifdef FUNCTION_GET_PMVOL
      double precision twh_get_pmvol
#undef FUNCTION_GET_PMVOL
#endif
c
#ifdef FUNCTION_GET_QBASEVALUE
      double precision twh_get_qbasevalue
#undef FUNCTION_GET_QBASEVALUE
#endif
c
#ifdef FUNCTION_GET_QBIFF
      character*(FFNAMELEN) twh_get_qbiff
#undef FUNCTION_GET_QBIFF
#endif
c
#ifdef FUNCTION_GET_QBINAMES
      character*(FFNAMELEN) twh_get_qbinames
#undef FUNCTION_GET_QBINAMES
#endif
c
#ifdef FUNCTION_GET_QBIVALUE
      double precision twh_get_qbivalue
#undef FUNCTION_GET_QBIVALUE
#endif
c
#ifdef FUNCTION_GET_QQATOM
      double precision twh_get_qqatom
#undef FUNCTION_GET_QQATOM
#endif
c
#ifdef FUNCTION_GET_RMVOL
      double precision twh_get_rmvol
#undef FUNCTION_GET_RMVOL
#endif
c
#ifdef FUNCTION_GET_SCALING_STYLE
      integer scp_get_scaling_style
#undef FUNCTION_GET_SCALING_STYLE
#endif
c
#ifdef FUNCTION_GET_SCALING_STYLE_STRING
      character*30 scp_get_scaling_style_string
#undef FUNCTION_GET_SCALING_STYLE_STRING
c
#endif
c
#ifdef FUNCTION_GET_STRING_LENGTH
      integer twh_get_string_length
#undef FUNCTION_GET_STRING_LENGTH
#endif
c
#ifdef FUNCTION_GET_TAVOL
      double precision twh_get_tavol
#undef FUNCTION_GET_TAVOL
#endif
c
#ifdef FUNCTION_GET_TORCOEFF
      double precision twh_get_torcoeff
#undef FUNCTION_GET_TORCOEFF
#endif
c
#ifdef FUNCTION_GET_TORFF
      character*(FFNAMELEN) twh_get_torff
#undef FUNCTION_GET_TORFF
#endif
c
#ifdef FUNCTION_GET_TORNAMES
      character*(FFNAMELEN) twh_get_tornames
#undef FUNCTION_GET_TORNAMES
#endif
c
#ifdef FUNCTION_GET_TORSTRING
      character*(FFNAMELEN) twh_get_torstring
#undef FUNCTION_GET_TORSTRING
#endif
c
#ifdef FUNCTION_GET_TORSTYLE
      integer twh_get_torstyle
#undef FUNCTION_GET_TORSTYLE
#endif
c
#ifdef FUNCTION_GET_VIBFF
      character*10 twh_get_vibff
#undef FUNCTION_GET_VIBFF
#endif
c
#ifdef FUNCTION_GET_VIBCOEFF
      double precision twh_get_vibcoeff
#undef FUNCTION_GET_VIBCOEFF
#endif
c
#ifdef FUNCTION_GET_VIBNAMES
      character*(FFNAMELEN) twh_get_vibnames
#undef FUNCTION_GET_VIBNAMES
#endif
c
#ifdef FUNCTION_GET_VIBRANG
      double precision twh_get_vibrang
#undef FUNCTION_GET_VIBRANG
#endif
c
#ifdef FUNCTION_IN_ATOMLIST
       logical twh_in_atomlist
#undef FUNCTION_IN_ATOMLIST
#endif
c
#ifdef FUNCTION_INTEGRATEDGAUSSPROB
      double precision twh_integratedgaussprob 
#undef FUNCTION_INTEGRATEDGAUSSPROB
#endif
c
#ifdef FUNCTION_INVERSELAWOFCOSINE
      double precision twh_inverselawofcosine
#undef FUNCTION_INVERSELAWOFCOSINE
#endif
c
#ifdef FUNCTION_LAWOFCOSINE
      double precision twh_lawofcosine
#undef FUNCTION_LAWOFCOSINE
#endif
c
#ifdef FUNCTION_LEN_TRIM
      integer twh_len_trim
#undef FUNCTION_LEN_TRIM
#endif
c
#ifdef FUNCTION_LINCLUDE
      logical twh_linclude
#undef FUNCTION_LINCLUDE
#endif
c
#ifdef FUNCTION_MAXBOXLENGTH
      double precision twh_maxboxlength 
#undef FUNCTION_MAXBOXLENGTH
#endif
c
#ifdef FUNCTION_MINBOXLENGTH
      double precision twh_minboxlength
#undef FUNCTION_MINBOXLENGTH
#endif
c
#ifdef FUNCTION_OLDGETATOMNUM
      integer twh_oldgetatomnum
#undef FUNCTION_OLDGETATOMNUM
#endif
c
#ifdef FUNCTION_ONEFIVETYPE
      integer twh_onefivetype
#undef FUNCTION_ONEFIVETYPE
#endif
c
#ifdef FUNCTION_PEEK_LABEL
      character*50 twh_peek_label
#undef FUNCTION_PEEK_LABEL
#endif
c
#ifdef FUNCTION_RANDOM
      double precision twh_random 
#undef FUNCTION_RANDOM
#endif
c
#ifdef FUNCTION_READ_DIR_STRING
      character*MAXDIRLENGTH twh_read_dir_string
#undef FUNCTION_READ_DIR_STRING
#endif
c
#ifdef FUNCTION_READ_FLOAT
      double precision twh_read_float
#undef FUNCTION_READ_FLOAT
#endif
c
#ifdef FUNCTION_READ_INTEGER
      integer twh_read_integer
#undef FUNCTION_READ_INTEGER
#endif
c
#ifdef FUNCTION_READ_LABELED_FLOAT
      double precision twh_read_labeled_float
#undef FUNCTION_READ_LABELED_FLOAT
#endif
c
#ifdef FUNCTION_READ_LABELED_INTEGER
      integer twh_read_labeled_integer
#undef FUNCTION_READ_LABELED_INTEGER
#endif
c
#ifdef FUNCTION_READ_LABELED_LOGICAL
      logical twh_read_labeled_logical
#undef FUNCTION_READ_LABELED_LOGICAL
#endif
c
#ifdef FUNCTION_READ_LOGICAL
      logical twh_read_logical
#undef FUNCTION_READ_LOGICAL
#endif
c
#ifdef FUNCTION_READ_STRING
      character*(50) twh_read_string
#undef FUNCTION_READ_STRING
#endif
c
#ifdef FUNCTION_SCALE_ATOMS
      logical scp_scale_atoms
#undef FUNCTION_SCALE_ATOMS
#endif
c
#ifdef FUNCTION_LIMITED_DOUBLE
      double precision twh_limited_double
#undef FUNCTION_LIMITED_DOUBLE
#endif
c
#ifdef FUNCTION_SAFE_DOUBLE
      double precision twh_safe_double
#undef FUNCTION_SAFE_DOUBLE
#endif
c
#ifdef FUNCTION_VANGANG
      double precision twh_vangang 
#undef FUNCTION_VANGANG
#endif
c
#ifdef FUNCTION_VANGLE
      double precision twh_vangle 
#undef FUNCTION_VANGLE
#endif
c
#ifdef FUNCTION_VBOND
      double precision twh_vbond 
#undef FUNCTION_VBOND
#endif
c
#ifdef FUNCTION_VBONBON
      double precision twh_vbonbon 
#undef FUNCTION_VBONBON
#endif
c
#ifdef FUNCTION_VCOULOMB
      double precision twh_vcoulomb 
#undef FUNCTION_VCOULOMB
#endif
c
#ifdef FUNCTION_VEEFONE
      double precision twh_veefone 
#undef FUNCTION_VEEFONE
#endif
c
#ifdef FUNCTION_VEMBED
      double precision twh_vembed
#undef FUNCTION_VEMBED
#endif
c
#ifdef FUNCTION_VFIELD
      double precision twh_vfield
#undef FUNCTION_VFIELD
#endif
c
#ifdef FUNCTION_VIMPROPER
      double precision twh_vimproper 
#undef FUNCTION_VIMPROPER
#endif
c
#ifdef FUNCTION_VONEFIVE
      double precision twh_vonefive
#undef FUNCTION_VONEFIVE
#endif
c
#ifdef FUNCTION_VSASA
      double precision twh_vsasa 
#undef FUNCTION_VSASA
#endif
c
#ifdef FUNCTION_VTHREEBODY
      double precision twh_vthreebody
#undef FUNCTION_VTHREEBODY
#endif
c
#ifdef FUNCTION_VTORSION
      double precision twh_vtorsion
#undef FUNCTION_VTORSION
#endif
c
#ifdef FUNCTION_VTWOBODY
      double precision twh_vtwobody 
#undef FUNCTION_VTWOBODY
#endif
c
#ifdef FUNCTION_WCOULOMB
      double precision twh_wcoulomb
#undef FUNCTION_WCOULOMB
#endif
c
#ifdef FUNCTION_WTWOBODY
      double precision twh_wtwobody
#undef FUNCTION_WTWOBODY
#endif
