#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2011 Marcus G. Martin                       *
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
      subroutine twh_ffmend2003(lfailure,ntypes)
c     ******************************************************************
c     * sets up parameters for the Mend2003 forcefield                 *
c     *                                                                *
c     * originally written 03-21-2005 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes
c     --- local variables
      integer ii,type,isame
      double precision zzz,dvalue
      character*(FFNAMELEN) ffname
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal

      ffname = 'Mend2003'
      call twh_classical_potential(GLB_SET,'Embedded Atom Method')
      call twh_classical_mixrule(GLB_SET,'Explicit')

c     sets up iron potentials 2 and 4 from the appendix of           
c     M.I. Mendelev, S. Han, D.J. Srolovitz, G.J. Ackland, D.Y. Sun,
c     M. Asta, Philosophcal Magazine, 83, 3977-3994 (2003).

c     --- initialize type
      type = 0

c     --- iron potental 2
      type = type + 1
      masslocal = 55.847d0
      elelocal =  'Fe'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Fe-p2'
      cbnlocal = 'Fe'
      canlocal = 'Fe'
      ctnlocal = 'Fe'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
c     --- pair potential is Ackland
      call twh_eam_pair_style(GLB_SET,ii,'Ackland 3-part')
c     --- r_1
      nbclocal(1) = 1.00d0
c     --- r_2
      nbclocal(2) = 2.00d0
c     --- temporary variable Z (atomic number)
      zzz = 26.0d0
c     --- Z^2 q^2
      nbclocal(3) = (zzz**2)*(CON_EETOKCAL*CON_KCALTOK)
c     --- r_s
      nbclocal(4) = (0.88534d0)*CON_BOHRRADIUS
     &     /(dsqrt(2.0d0)*(zzz**(1.0d0/3.0d0)))
c     --- B_0
      nbclocal(5) = 6.4265260576348d0
c     --- B_1
      nbclocal(6) = 1.7900488524286d0
c     --- B_2
      nbclocal(7) = -4.5108316729807d0
c     --- B_3
      nbclocal(8) = 1.0866199373306d0
c     --- prefactor I added to convert this exponential from ev to K
      nbclocal(9) = CON_EVTOK
      call twh_set_nbcoeff(ii,nbclocal) 
c     --- now the tabulated part of this 3-part potential
      isame = 0
c     --- r_2 and a_2
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.2d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-24.028204854115d0)*CON_EVTOK)
c     --- r_3 and a_3
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(11.300691696477d0)*CON_EVTOK)
c     --- r_4 and a_4
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.4d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(5.3144495820462d0)*CON_EVTOK)
c     --- r_5 and a_5
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.5d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-4.6659532856049d0)*CON_EVTOK)
c     --- r_6 and a_6
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.6d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(5.9637758529194d0)*CON_EVTOK)
c     --- r_7 and a_7
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.7d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-1.7710262006061d0)*CON_EVTOK)
c     --- r_8 and a_8
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.8d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.85913830768731d0)*CON_EVTOK)
c     --- r_9 and a_9
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.0d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-2.1845362968261d0)*CON_EVTOK)
c     --- r_10 and a_10
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(2.6424377007466d0)*CON_EVTOK)
c     --- r_11 and a_11
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.7d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-1.0358345370208d0)*CON_EVTOK)
c     --- r_12 and a_12
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.2d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.33548264951582d0)*CON_EVTOK)
c     --- r_13 and a_13
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.7d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.046448582149334d0)*CON_EVTOK)
c     --- r_14 and a_14
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.0070294963048689d0)*CON_EVTOK)
c     --- set table_npair
      call twh_table_npair(GLB_SET,ii,isame)
c     --- density - store self density
      call twh_eam_dens_style(GLB_SET,type,type,'Ackland cubic sum')
      isame = 0
c     --- density r_1 and a_1
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,type,type,isame,2.4d0)
      call twh_eam_dens(GLB_SET,2,type,type,isame,11.686859407970d0)
c     --- density r_2 and a_2
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,type,type,isame,3.2d0)
      call twh_eam_dens(GLB_SET,2,type,type,isame,-0.014710740098830d0)
c     --- density r_3 and a_3
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,type,type,isame,4.2d0)
      call twh_eam_dens(GLB_SET,2,type,type,isame,0.47193527075943d0)
c     --- set eam_ndens
      call twh_eam_ndens(GLB_SET,type,type,isame)
c     --- Embedding function
      call twh_eam_embed_style(GLB_SET,type,'power 0.5 and 2')
      isame = 0
c     --- square root term
      isame = isame + 1
      dvalue = (-1.0d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,2,type,isame,dvalue)
c     --- power 2 term
      isame = isame + 1
      dvalue = (-0.00035387096579929d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,2,type,isame,dvalue)
c     --- set nembed
      call twh_eam_nembed(GLB_SET,type,isame)

c     --- Fake bit as a placeholder for the cross between these as 
c     --- you really should not use them together
      ii = twh_getnbtype(type,type+1)
c     --- no pair potential
      call twh_eam_pair_style(GLB_SET,ii,'none')
      call twh_table_npair(GLB_SET,ii,0)
c     --- cross densities also don't exist
      call twh_eam_dens_style(GLB_SET,type+1,type,'none')
      call twh_eam_ndens(GLB_SET,type+1,type,0)
      call twh_eam_dens_style(GLB_SET,type,type+1,'none')
      call twh_eam_ndens(GLB_SET,type,type+1,0)
c     --- end of fake bit

c     --- iron potental 4
      type = type + 1
      masslocal = 55.847d0
      elelocal =  'Fe'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Fe-p4'
      cbnlocal = 'Fe'
      canlocal = 'Fe'
      ctnlocal = 'Fe'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
c     --- pair potential is Ackland
      call twh_eam_pair_style(GLB_SET,ii,'Ackland 3-part')
c     --- r_1
      nbclocal(1) = 0.90d0
c     --- r_2
      nbclocal(2) = 1.95d0
c     --- temporary variable Z (atomic number)
      zzz = 26.0d0
c     --- Z^2 q^2
      nbclocal(3) = (zzz**2)*(CON_EETOKCAL*CON_KCALTOK)
c     --- r_s
      nbclocal(4) = (0.88534d0)*CON_BOHRRADIUS
     &     /(dsqrt(2.0d0)*(zzz**(1.0d0/3.0d0)))
c     --- B_0
      nbclocal(5) = 14.996917289290d0
c     --- B_1
      nbclocal(6) = -20.533174190155d0
c     --- B_2
      nbclocal(7) = 14.002591780752d0
c     --- B_3
      nbclocal(8) = -3.6473736591143d0
c     --- prefactor I added to convert this exponential from ev to K
      nbclocal(9) = CON_EVTOK
      call twh_set_nbcoeff(ii,nbclocal) 
c     --- now the tabulated part of this 3-part potential
      isame = 0
c     --- r_1 and a_1
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.1d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(195.92322853994d0)*CON_EVTOK)
c     --- r_2 and a_2
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.2d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(17.516698453315d0)*CON_EVTOK)
c     --- r_3 and a_3
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(1.4926525164290d0)*CON_EVTOK)
c     --- r_4 and a_4
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.4d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(6.4129476125197d0)*CON_EVTOK)
c     --- r_5 and a_5
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.5d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-6.8157461860553d0)*CON_EVTOK)
c     --- r_6 and a_6
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.6d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(9.6582581963600d0)*CON_EVTOK)
c     --- r_7 and a_7
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.7d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-5.3419002764419d0)*CON_EVTOK)
c     --- r_8 and a_8
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.8d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(1.7996558048346d0)*CON_EVTOK)
c     --- r_9 and a_9
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.0d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-1.4788966636288d0)*CON_EVTOK)
c     --- r_10 and a_10
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(1.8530435283665d0)*CON_EVTOK)
c     --- r_11 and a_11
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.7d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.64164344859316d0)*CON_EVTOK)
c     --- r_12 and a_12
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.2d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.24463630025168d0)*CON_EVTOK)
c     --- r_13 and a_13
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.7d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.057721650527383d0)*CON_EVTOK)
c     --- r_14 and a_14
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.023358616514826d0)*CON_EVTOK)
c     --- r_15 and a_15
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,6.0d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.0097064921265079d0)*CON_EVTOK)
c     --- set table_npair
      call twh_table_npair(GLB_SET,ii,isame)
c     --- density - stored in arrays depending only on type
      call twh_eam_dens_style(GLB_SET,type,type,'Ackland cubic sum')
      isame = 0
c     --- density r_1 and a_1
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,type,type,isame,2.4d0)
      call twh_eam_dens(GLB_SET,2,type,type,isame,11.686859407970d0)
c     --- density r_2 and a_2
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,type,type,isame,3.2d0)
      call twh_eam_dens(GLB_SET,2,type,type,isame,-0.014710740098830d0)
c     --- density r_3 and a_3
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,type,type,isame,4.2d0)
      call twh_eam_dens(GLB_SET,2,type,type,isame,0.47193527075943d0)
c     --- set eam_ndens
      call twh_eam_ndens(GLB_SET,type,type,isame)
c     --- Embedding function
      call twh_eam_embed_style(GLB_SET,type,'power 0.5 and 2')
      isame = 0
c     --- square root term
      isame = isame + 1
      dvalue = (-1.0d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,2,type,isame,dvalue)
c     --- square term
      isame = isame + 1
      dvalue = (-0.00034906178363530d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,2,type,isame,dvalue)
c     --- set nembed
      call twh_eam_nembed(GLB_SET,type,isame)

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end
