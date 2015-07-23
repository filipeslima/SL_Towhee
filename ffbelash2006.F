#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2013-2013 Marcus G. Martin                       *
c     * see the file license.gpl for the full license information      *
c     *                                                                *
c     * This program is free software; you can redistribute it and/or  *
c     * modify it under the terms of the GNU General Public License    *
c     * as published by the Free Software Foundation; either version 2 *
c     * of the License, or (at your option) any later version.         *
c     *                                                                *
c     * This program is distributed in the hope that it will be useful,*
c     * but WITHOUT ANY WARRANTY; without even the implied warranty of *
c     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
c     * GNU General Public License for more details.                   *
c     *                                                                *
c     * You should have received a copy of the GNU General Public      *
c     * License along with this program; if not, write to the Free     *
c     * Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,*
c     * MA  02111-1307, USA.                                           *
c     ******************************************************************
      subroutine twh_ffbelash2006(lfailure,ntypes)
c     ******************************************************************
c     * sets up parameters for the belash2006 forcefield               *
c     *                                                                *
c     * originally written 07-15-2013 by M.G. Martin                   *
c     * last modified 07-15-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes
c     --- local variables
      integer ii,type,isame
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) ffname
c     --- double precision scalars
      double precision masslocal,dvalue
c     --- double precision arrays
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal

      ffname = 'Belash2006'
      call twh_classical_potential(GLB_SET,'Embedded Atom Method')
      call twh_classical_mixrule(GLB_SET,'Explicit')

c     sets up the mercury eam potential from
c     D.K. Belashchenko, "Application of the Embedded Atom Model to 
c     Liquid Metals: Liquid Mercury", High Temperature, 44, 5, (2006)
c     pp. 675-686

c     --- mercury
      type = 1
      masslocal = 200.59d0
      elelocal =  'Hg'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Hg'
      cbnlocal = 'Hg'
      canlocal = 'Hg'
      ctnlocal = 'Hg'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- Pair interactions with self
      ii = twh_getnbtype(type,type)
c     --- pair potential is tabulated
      call twh_eam_pair_style(GLB_SET,ii,'table')
c     --- tabulated
      isame = 0
c     --- r_1 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.400d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.483514d0)*CON_EVTOK)
c     --- r_2 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.450d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.351702d0)*CON_EVTOK)
c     --- r_3 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.500d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.247896d0)*CON_EVTOK)
c     --- r_4 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.550d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.126681d0)*CON_EVTOK)
c     --- r_5 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.600d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.127174d0)*CON_EVTOK)
c     --- r_6 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.650d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.663966d-1)*CON_EVTOK)
c     --- r_7 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.700d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.475004d-1)*CON_EVTOK)
c     --- r_8 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.750d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.161462d-1)*CON_EVTOK)
c     --- r_9 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.800d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.828018d-3)*CON_EVTOK)
c     --- r_10 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.850d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.178802d-1)*CON_EVTOK)
c     --- r_11 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.900d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.300109d-1)*CON_EVTOK)
c     --- r_12 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.950d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.404086d-1)*CON_EVTOK)
c     --- r_13 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.000d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.484684d-1)*CON_EVTOK)
c     --- r_14 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.050d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.547106d-1)*CON_EVTOK)
c     --- r_15 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.100d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.587955d-1)*CON_EVTOK)
c     --- r_16 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.150d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.613250d-1)*CON_EVTOK)
c     --- r_17 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.200d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.616199d-1)*CON_EVTOK)
c     --- r_18 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.250d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.606173d-1)*CON_EVTOK)
c     --- r_19 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.300d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.588065d-1)*CON_EVTOK)
c     --- r_20 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.350d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.564351d-1)*CON_EVTOK)
c     --- r_21 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.400d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.560472d-1)*CON_EVTOK)
c     --- r_22 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.450d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.565545d-1)*CON_EVTOK)
c     --- r_23 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.500d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.578759d-1)*CON_EVTOK)
c     --- r_24 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.550d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.580731d-1)*CON_EVTOK)
c     --- r_25 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.600d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.574882d-1)*CON_EVTOK)
c     --- r_26 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.650d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.563925d-1)*CON_EVTOK)
c     --- r_27 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.700d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.555761d-1)*CON_EVTOK)
c     --- r_28 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.750d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.542019d-1)*CON_EVTOK)
c     --- r_29 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.800d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.528850d-1)*CON_EVTOK)
c     --- r_30 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.850d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.508982d-1)*CON_EVTOK)
c     --- r_31 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.900d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.488353d-1)*CON_EVTOK)
c     --- r_32 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.950d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.461914d-1)*CON_EVTOK)
c     --- r_33 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.000d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.443331d-1)*CON_EVTOK)
c     --- r_34 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.050d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.418874d-1)*CON_EVTOK)
c     --- r_35 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.100d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.395831d-1)*CON_EVTOK)
c     --- r_36 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.150d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.368057d-1)*CON_EVTOK)
c     --- r_37 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.200d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.341666d-1)*CON_EVTOK)
c     --- r_38 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.250d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.308596d-1)*CON_EVTOK)
c     --- r_39 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.300d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.280620d-1)*CON_EVTOK)
c     --- r_40 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.350d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.257942d-1)*CON_EVTOK)
c     --- r_41 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.400d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.227871d-1)*CON_EVTOK)
c     --- r_42 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.450d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.203724d-1)*CON_EVTOK)
c     --- r_43 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.500d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.175871d-1)*CON_EVTOK)
c     --- r_44 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.550d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.153318d-1)*CON_EVTOK)
c     --- r_45 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.600d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.131229d-1)*CON_EVTOK)
c     --- r_46 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.650d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.115198d-1)*CON_EVTOK)
c     --- r_47 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.700d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.978748d-2)*CON_EVTOK)
c     --- r_48 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.750d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.818997d-2)*CON_EVTOK)
c     --- r_49 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.800d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.650784d-2)*CON_EVTOK)
c     --- r_50 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.850d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.489328d-2)*CON_EVTOK)
c     --- r_51 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.900d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.345595d-2)*CON_EVTOK)
c     --- r_52 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.950d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.187813d-2)*CON_EVTOK)
c     --- r_53 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.000d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.425142d-3)*CON_EVTOK)
c     --- r_54 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.050d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.110518d-2)*CON_EVTOK)
c     --- r_55 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.100d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.221891d-2)*CON_EVTOK)
c     --- r_56 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.150d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.332688d-2)*CON_EVTOK)
c     --- r_57 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.200d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.371251d-2)*CON_EVTOK)
c     --- r_58 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.250d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.429119d-2)*CON_EVTOK)
c     --- r_59 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.300d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.450934d-2)*CON_EVTOK)
c     --- r_60 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.350d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.454176d-2)*CON_EVTOK)
c     --- r_61 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.400d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.435775d-2)*CON_EVTOK)
c     --- r_62 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.450d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.436811d-2)*CON_EVTOK)
c     --- r_63 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.500d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.416698d-2)*CON_EVTOK)
c     --- r_64 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.550d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.370785d-2)*CON_EVTOK)
c     --- r_65 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.600d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.349032d-2)*CON_EVTOK)
c     --- r_66 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.650d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.338523d-2)*CON_EVTOK)
c     --- r_67 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.700d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.297338d-2)*CON_EVTOK)
c     --- r_68 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.750d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.268336d-2)*CON_EVTOK)
c     --- r_69 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.800d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.244778d-2)*CON_EVTOK)
c     --- r_70 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.850d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.212314d-2)*CON_EVTOK)
c     --- r_71 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.900d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.166614d-2)*CON_EVTOK)
c     --- r_72 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.950d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.962741d-3)*CON_EVTOK)
c     --- r_73 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.000d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.416439d-3)*CON_EVTOK)
c     --- r_74 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.050d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.348204d-3)*CON_EVTOK)
c     --- r_75 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.100d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.150444d-2)*CON_EVTOK)
c     --- r_76 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.150d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.232437d-2)*CON_EVTOK)
c     --- r_77 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.200d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.331693d-2)*CON_EVTOK)
c     --- r_78 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.250d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.430991d-2)*CON_EVTOK)
c     --- r_79 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.300d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.515263d-2)*CON_EVTOK)
c     --- r_80 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.350d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.569816d-2)*CON_EVTOK)
c     --- r_81 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.400d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.634562d-2)*CON_EVTOK)
c     --- r_82 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.450d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.681474d-2)*CON_EVTOK)
c     --- r_83 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.500d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.703123d-2)*CON_EVTOK)
c     --- r_84 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.550d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.720794d-2)*CON_EVTOK)
c     --- r_85 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.600d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.757106d-2)*CON_EVTOK)
c     --- r_86 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.650d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.789255d-2)*CON_EVTOK)
c     --- r_87 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.700d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.800038d-2)*CON_EVTOK)
c     --- r_88 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.750d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.829357d-2)*CON_EVTOK)
c     --- r_89 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.800d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.856128d-2)*CON_EVTOK)
c     --- r_90 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.850d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.872009d-2)*CON_EVTOK)
c     --- r_91 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.900d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.885832d-2)*CON_EVTOK)
c     --- r_92 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.950d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.894230d-2)*CON_EVTOK)
c     --- r_93 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.000d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.909542d-2)*CON_EVTOK)
c     --- r_94 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.050d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.895489d-2)*CON_EVTOK)
c     --- r_95 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.100d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.890581d-2)*CON_EVTOK)
c     --- r_96 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.150d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.884023d-2)*CON_EVTOK)
c     --- r_97 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.200d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.831443d-2)*CON_EVTOK)
c     --- r_98 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.250d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.796555d-2)*CON_EVTOK)
c     --- r_99 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.300d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.747324d-2)*CON_EVTOK)
c     --- r_100 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.350d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.695746d-2)*CON_EVTOK)
c     --- r_101 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.400d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.667990d-2)*CON_EVTOK)
c     --- r_102 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.450d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.604279d-2)*CON_EVTOK)
c     --- r_103 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.500d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.560247d-2)*CON_EVTOK)
c     --- r_104 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.550d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.516298d-2)*CON_EVTOK)
c     --- r_105 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.600d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.484546d-2)*CON_EVTOK)
c     --- r_106 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.650d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.457130d-2)*CON_EVTOK)
c     --- r_107 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.700d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.410413d-2)*CON_EVTOK)
c     --- r_108 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.750d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.386006d-2)*CON_EVTOK)
c     --- r_109 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.800d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.348579d-2)*CON_EVTOK)
c     --- r_110 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.850d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.291792d-2)*CON_EVTOK)
c     --- r_111 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.900d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.248405d-2)*CON_EVTOK)
c     --- r_112 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,7.950d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.201067d-2)*CON_EVTOK)
c     --- r_113 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,8.000d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.160872d-2)*CON_EVTOK)
c     --- r_114 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,8.050d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.119191d-2)*CON_EVTOK)
c     --- r_115 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,8.100d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.534409d-3)*CON_EVTOK)
c     --- r_116 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,8.150d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.159762d-3)*CON_EVTOK)
c     --- r_117 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,8.200d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.183094d-2)*CON_EVTOK)
c     --- r_118 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,8.250d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.183108d-2)*CON_EVTOK)
c     --- r_119 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,8.300d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.183121d-2)*CON_EVTOK)
c     --- r_120 and pair
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,8.350d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.0d0)*CON_EVTOK)
c     --- set table_npair
      call twh_table_npair(GLB_SET,ii,isame)

c     --- density with self
      call twh_eam_dens_style(GLB_SET,type,type,'exponential')
      isame = 0
c     --- density prefactor A (p1 in the paper)
      isame = isame + 1
      dvalue = (4.79860d0)
      call twh_eam_dens(GLB_SET,2,type,type,isame,dvalue)
c     --- density exponential factor B (-p2 in the paper)
      isame = isame + 1
      dvalue = -(1.30950d0)
      call twh_eam_dens(GLB_SET,2,type,type,isame,dvalue)
c     --- maximum density cutoff
      dvalue = 9.01d0
      call twh_eam_dens(GLB_SET,1,type,type,isame,dvalue)
c     --- set eam_ndens
      call twh_eam_ndens(GLB_SET,type,type,isame)

c     --- Embedding function 
      call twh_eam_embed_style(GLB_SET,type,'Belashchenko')
c     --- a1 + a2*(rho - rho0)^2 + a3*(rho - rho0)^3 for rho >= 0.8*rho0
c     --- alpha sqrt(rho) + beta*rho for rho < 0.8*rho0
      isame = 0
c     --- a1
      isame = isame + 1
      dvalue = (-0.1628d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,2,type,isame,dvalue)
c     --- a2
      isame = isame + 1
      dvalue = (0.7075)*CON_EVTOK
      call twh_eam_embed(GLB_SET,2,type,isame,dvalue)
c     --- a3
      isame = isame + 1
      dvalue = (1.3630d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,2,type,isame,dvalue)
c     --- alpha and beta
      isame = isame + 1
      dvalue = (-0.111472d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,1,type,isame,dvalue)
c     --- beta
      dvalue = (-0.057125d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,2,type,isame,dvalue)
c     --- rhomin (0.8*rho0) and rho0
      isame = isame + 1
c     --- default rho0 is 1
      dvalue = 1.0d0
      call twh_eam_embed(GLB_SET,2,type,isame,dvalue)
c     --- rhomin = 0.8*rho0
      dvalue = 0.8d0*dvalue
      call twh_eam_embed(GLB_SET,1,type,isame,dvalue)
c     --- set nembed
      call twh_eam_nembed(GLB_SET,type,isame)


      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end
