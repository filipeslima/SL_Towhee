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
      subroutine twh_ffackl2004(lfailure,ntypes)
c     ******************************************************************
c     * sets up parameters for the Ackl2004 forcefield                 *
c     *                                                                *
c     * originally written 03-26-2005 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes
c     --- local variables
      integer ii,type,isame,jtype
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) ffname
c     --- double precision scalars
      double precision zzz,masslocal,dvalue
c     --- double precision arrays
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal

      ffname = 'Ackl2004'
      call twh_classical_potential(GLB_SET,'Embedded Atom Method')
      call twh_classical_mixrule(GLB_SET,'Explicit')

c     sets up iron and phosphorus potentials from
c     G.J. Ackland, M.I. Mendelev, D.J. Srolovitz, S. Han, A.V. Barashev
c     J. Phys. Cond. Matt. 16 (2004) S2629-S2642

c     --- iron
      type = 1
      masslocal = 55.847d0
      elelocal =  'Fe'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Fe'
      cbnlocal = 'Fe'
      canlocal = 'Fe'
      ctnlocal = 'Fe'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- Embedding function 
      call twh_eam_embed_style(GLB_SET,type,'power 0.5, 2, and 4')
      isame = 0
c     --- square root factor
      isame = isame + 1
      dvalue = (-1.0d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,2,type,isame,dvalue)
c     --- square factor
      isame = isame + 1
      dvalue = (-6.7314115586063d-4)*CON_EVTOK
      call twh_eam_embed(GLB_SET,2,type,isame,dvalue)
c     --- 4th power factor
      isame = isame + 1
      dvalue = (7.6514905604792d-8)*CON_EVTOK
      call twh_eam_embed(GLB_SET,2,type,isame,dvalue)
c     --- set nembed
      call twh_eam_nembed(GLB_SET,type,isame)
c     --- density with self
      jtype = 1
      call twh_eam_dens_style(GLB_SET,jtype,type,'Ackland cubic sum')
      isame = 0
c     --- density r_1 and a_1
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,jtype,type,isame,2.4d0)
      call twh_eam_dens(GLB_SET,2,jtype,type,isame,11.686859407970d0)
c     --- density r_2 and a_2
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,jtype,type,isame,3.2d0)
      call twh_eam_dens(GLB_SET,2,jtype,type,isame,-0.014710740098830d0)
c     --- density r_3 and a_3
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,jtype,type,isame,4.2d0)
      call twh_eam_dens(GLB_SET,2,jtype,type,isame,0.47193527075943d0)
c     --- set eam_ndens
      call twh_eam_ndens(GLB_SET,jtype,type,isame)
c     --- density with phosphorous
      jtype = 2
      call twh_eam_dens_style(GLB_SET,jtype,type,'Ackland cubic sum')
      isame = 0
c     --- density r_1 and a_1
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,jtype,type,isame,2.4d0)
      dvalue = 11.686859407970d0*((21.0d0/24.0d0)**2)
      call twh_eam_dens(GLB_SET,2,jtype,type,isame,dvalue)
c     --- density r_2 and a_2
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,jtype,type,isame,3.2d0)
      dvalue = -0.014710740098830d0*((21.0d0/24.0d0)**2)
      call twh_eam_dens(GLB_SET,2,jtype,type,isame,dvalue)
c     --- density r_3 and a_3
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,jtype,type,isame,4.2d0)
      dvalue = 0.47193527075943d0*((21.0d0/24.0d0)**2)
      call twh_eam_dens(GLB_SET,2,jtype,type,isame,dvalue)
c     --- set eam_ndens
      call twh_eam_ndens(GLB_SET,jtype,type,isame)
c     --- Pair interactions with self
      ii = twh_getnbtype(type,type)
c     --- pair potential is Ackland 3-part
      call twh_eam_pair_style(GLB_SET,ii,'Ackland 3-part')
c     --- r_1
      nbclocal(1) = 1.00d0
c     --- r_2
      nbclocal(2) = 2.05d0
c     --- temporary variable Z (atomic number)
      zzz = 26.0d0
c     --- Z^2 q^2
      nbclocal(3) = (9.7342365892908E+03)*CON_EVTOK
c     --- r_s
      nbclocal(4) = (0.88534d0)*CON_BOHRRADIUS
     &     /(dsqrt(2.0d0)*(zzz**(1.0d0/3.0d0)))
c     --- B_0
      nbclocal(5) = 7.4122709384068d0
c     --- B_1
      nbclocal(6) = -0.64180690713367d0
c     --- B_2
      nbclocal(7) = -2.6043547961722d0
c     --- B_3
      nbclocal(8) = 0.6262539393123d0
c     --- prefactor I added to convert this exponential from ev to K
      nbclocal(9) = CON_EVTOK
      call twh_set_nbcoeff(ii,nbclocal) 
c     --- now the tabulated part of this 3-part potential
      isame = 0
c     --- r_2 and a_2
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.2d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-27.444805994228d0)*CON_EVTOK)
c     --- r_3 and a_3
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(15.738054058489d0)*CON_EVTOK)
c     --- r_4 and a_4
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.4d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(2.2077118733936d0)*CON_EVTOK)
c     --- r_5 and a_5
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.5d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-2.4989799053251d0)*CON_EVTOK)
c     --- r_6 and a_6
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.6d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(4.2099676494795d0)*CON_EVTOK)
c     --- r_7 and a_7
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.7d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.77361294129713d0)*CON_EVTOK)
c     --- r_8 and a_8
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,2.8d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.80656414937789d0)*CON_EVTOK)
c     --- r_9 and a_9
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.0d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-2.3194358924605d0)*CON_EVTOK)
c     --- r_10 and a_10
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(2.6577406128280d0)*CON_EVTOK)
c     --- r_11 and a_11
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,3.7d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-1.0260416933564d0)*CON_EVTOK)
c     --- r_12 and a_12
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.2d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(0.35018615891957d0)*CON_EVTOK)
c     --- r_13 and a_13
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,4.7d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.058531821042271d0)*CON_EVTOK)
c     --- r_14 and a_14
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-0.0030458824556234d0)*CON_EVTOK)
c     --- set table_npair
      call twh_table_npair(GLB_SET,ii,isame)

c     --- cross terms with phosphorous
      jtype = 2
      ii = twh_getnbtype(type,jtype)
c     --- pair potential is Ackland Power
      call twh_eam_pair_style(GLB_SET,ii,'Ackland Power')
c     --- r_1
      nbclocal(1) = 1.0d0
c     --- r_2
      nbclocal(2) = 2.0d0
c     --- prefactor
      nbclocal(3) = (5.6159057245908d03)*CON_EVTOK
c     --- exp-1 pre
      nbclocal(4) = 0.1818d0
c     --- exp-1 coeff
      nbclocal(5) = -2.6329090970098d01
c     --- exp-2 pre
      nbclocal(6) = 0.5099d0
c     --- exp-2 coeff
      nbclocal(7) = -7.7530945066009d0
c     --- exp-3 pre
      nbclocal(8) = 0.2802d0
c     --- exp-3 coeff
      nbclocal(9) = -3.3149971099538d0
c     --- exp-4 pre
      nbclocal(10) = 0.02817d0
c     --- exp-4 coeff
      nbclocal(11) = -1.6587327311161d0
c     --- B_0
      nbclocal(12) = 1.0761854424880d1
c     --- B_1
      nbclocal(13) = -1.0004045788895d1
c     --- B_2
      nbclocal(14) = 4.9854254472397d0
c     --- B_3
      nbclocal(15) = -1.2599788569372d0
c     --- prefactor I added to convert this exponential from ev to K
      nbclocal(16) = CON_EVTOK
      call twh_set_nbcoeff(ii,nbclocal) 
c     --- now the power series
      isame = 0
c     --- 4th power
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-3.3136605743629d0)*CON_EVTOK)
c     --- 5th power
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(1.2625238193602d1)*CON_EVTOK)
c     --- 6th power
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-2.0361693308072d1)*CON_EVTOK)
c     --- 7th power
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(1.7629292543942d1)*CON_EVTOK)
c     --- 8th power
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-8.8120728047659d0)*CON_EVTOK)
c     --- 9th power
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(2.5494288609989d0)*CON_EVTOK)
c     --- 10th power
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-3.9698390783403d-1)*CON_EVTOK)
c     --- 11th power
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(2.5779015833433d-2)*CON_EVTOK)
c     --- set table_npair
      call twh_table_npair(GLB_SET,ii,isame)

c     --- phosphorous
      type = 2
      masslocal = 30.974d0
      elelocal =  ' P'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'P'
      cbnlocal = 'P'
      canlocal = 'P'
      ctnlocal = 'P'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- Embedding function
      call twh_eam_embed_style(GLB_SET,type,'power 0.5 and 2')
      isame = 0
c     --- square root factor
      isame = isame + 1
      dvalue = (-1.0d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,2,type,isame,dvalue)
c     --- square factor
      isame = isame + 1
      dvalue = (0.001950274540243d0)*CON_EVTOK
      call twh_eam_embed(GLB_SET,2,type,isame,dvalue)
c     --- set nembed
      call twh_eam_nembed(GLB_SET,type,isame)
c     --- density with iron
      jtype = 1
      call twh_eam_dens_style(GLB_SET,jtype,type,'Ackland cubic sum')
      isame = 0
c     --- density r_1 and a_1
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,jtype,type,isame,2.4d0)
      dvalue = 11.686859407970d0*((21.0d0/24.0d0)**2)
      call twh_eam_dens(GLB_SET,2,jtype,type,isame,dvalue)
c     --- density r_2 and a_2
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,jtype,type,isame,3.2d0)
      dvalue = -0.014710740098830d0*((21.0d0/24.0d0)**2)
      call twh_eam_dens(GLB_SET,2,jtype,type,isame,dvalue)
c     --- density r_3 and a_3
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,jtype,type,isame,4.2d0)
      dvalue = 0.47193527075943d0*((21.0d0/24.0d0)**2)
      call twh_eam_dens(GLB_SET,2,jtype,type,isame,dvalue)
c     --- set eam_ndens
      call twh_eam_ndens(GLB_SET,jtype,type,isame)
c     --- density with self (P-P)
      jtype = 2
      call twh_eam_dens_style(GLB_SET,jtype,type,'Ackland cubic sum')
      isame = 0
c     --- density r_1 and a_1
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,jtype,type,isame,2.4d0)
      dvalue = 11.686859407970d0*((21.0d0/24.0d0)**4)
      call twh_eam_dens(GLB_SET,2,jtype,type,isame,dvalue)
c     --- density r_2 and a_2
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,jtype,type,isame,3.2d0)
      dvalue = -0.014710740098830d0*((21.0d0/24.0d0)**4)
      call twh_eam_dens(GLB_SET,2,jtype,type,isame,dvalue)
c     --- density r_3 and a_3
      isame = isame + 1
      call twh_eam_dens(GLB_SET,1,jtype,type,isame,4.2d0)
      dvalue = 0.47193527075943d0*((21.0d0/24.0d0)**4)
      call twh_eam_dens(GLB_SET,2,jtype,type,isame,dvalue)
c     --- set eam_ndens
      call twh_eam_ndens(GLB_SET,jtype,type,isame)
c     --- cross terms with self (P-P)
      jtype = 2
      ii = twh_getnbtype(type,jtype)
c     --- pair potential is Ackland Power
      call twh_eam_pair_style(GLB_SET,ii,'Ackland Power')
c     --- r_1
      nbclocal(1) = 0.9d0
c     --- r_2
      nbclocal(2) = 2.5d0
c     --- prefactor
      nbclocal(3) = (3.2399456103409d3)*CON_EVTOK
c     --- exp-1 pre
      nbclocal(4) = 0.1818d0
c     --- exp-1 coeff
      nbclocal(5) = -2.3822786399080d1
c     --- exp-2 pre
      nbclocal(6) = 0.5099d0
c     --- exp-2 coeff
      nbclocal(7) = -7.0150661324540d0
c     --- exp-3 pre
      nbclocal(8) = 0.2802d0
c     --- exp-3 coeff
      nbclocal(9) = -2.9994377000591d0
c     --- exp-4 pre
      nbclocal(10) = 0.02817d0
c     --- exp-4 coeff
      nbclocal(11) = -1.5008355431420d0
c     --- B_0
      nbclocal(12) = 9.9382842499617d0
c     --- B_1
      nbclocal(13) = -8.5637164272526d0
c     --- B_2
      nbclocal(14) = +3.4519627285990d0
c     --- B_3
      nbclocal(15) = -6.1453831350215d-01
c     --- prefactor I added to convert this exponential from ev to K
      nbclocal(16) = CON_EVTOK
      call twh_set_nbcoeff(ii,nbclocal) 
c     --- now the power series
      isame = 0
c     --- 4th power
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(-7.8293794709143d-02)*CON_EVTOK)
c     --- 5th power
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,5.3d0)
      call twh_table_pair(GLB_SET,2,ii,isame
     &     ,(+3.7557214911646d-02)*CON_EVTOK)
c     --- set table_npair
      call twh_table_npair(GLB_SET,ii,isame)

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end
