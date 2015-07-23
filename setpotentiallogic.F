#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2012 Marcus G. Martin                       *
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
      subroutine twh_setpotentiallogic(lfailure,ldomixrule,ldolshift
     &,ldoltailc,ldormin,ldorcut,ldorcutin,ldointerpolatestyle
     &,ldorpd,ldoscaling)
c     ******************************************************************
c     * sets up all of the input and output logical related to the     *
c     * classical_potential                                            *
c     *                                                                *
c     * split out of readclassical on 04-13-2005 by M.G. Martin        *
c     * last modified 08-31-2012 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical ldomixrule,ldolshift,ldoltailc,ldormin,ldorcut
      logical ldorcutin,ldointerpolatestyle,ldorpd,lfailure
      logical ldoscaling
c     --- local variables
      integer cpnumber,maxchain

c     --- get the classical potential number
      call twh_cp_number(GLB_GET,cpnumber)

c     --- initialize potential dependent reads
      ldomixrule = .false.
      ldolshift = .false.
      ldoltailc = .false.
      ldormin = .false.
      ldorcut = .false.
      ldorcutin = .false.
      ldointerpolatestyle = .false.
      ldorpd = .false.
      ldoscaling = .false.
      
      if ( cpnumber .eq. CP_LENNARD_JONES .or.
     &     cpnumber .eq. CP_9_6 .or.
     &     cpnumber .eq. CP_12_6_PLUS_SOLVATION .or.
     &     cpnumber .eq. CP_12_6_PLUS_12_10_H_BOND .or.
     &     cpnumber .eq. CP_12_9_6 .or.
     &     cpnumber .eq. CP_GORDON_N_6 .or.
     &     cpnumber .eq. CP_UFF_12_6 .or.
     &     cpnumber .eq. CP_DRUKKER_H_BOND ) then
         ldomixrule = .true.
         ldolshift = .true.
         ldoltailc = .true.
         ldormin = .true.
         ldorcut = .true.
         ldorcutin = .true.
      elseif ( cpnumber .eq. CP_SCALED_LENNARD_JONES .or.
     &        cpnumber .eq. CP_WEEKS_CHANDLER_ANDERSON ) then
         ldomixrule = .true.
         ldolshift = .true.
         ldormin = .true.
         ldorcut = .true.
         ldorcutin = .true.
         ldoscaling = .true.
      elseif ( cpnumber .eq. CP_EXPONENTIAL_6 ) then
         ldomixrule = .true.
         ldolshift = .true.
         ldoltailc = .true.
         ldorcut = .true.
         ldorcutin = .true.
      elseif ( cpnumber .eq. CP_HARD_SPHERE .or.
     &        cpnumber .eq. CP_SQUARE_WELL ) then
         ldomixrule = .true.
         ldorpd = .true.
      elseif ( cpnumber .eq. CP_REPULSIVE_SPHERE .or.
     &        cpnumber .eq. CP_REPULSIVE_WELL .or.
     &        cpnumber .eq. CP_MULTIWELL .or.
     &        cpnumber .eq. CP_REPULSIVE_MULTIWELL ) then
c     --- repulsive sphere or repulsive well, or multiwell
         ldomixrule = .true.
      elseif ( cpnumber .eq. CP_EXPONENTIAL_12_6 ) then
c        --- exponential 6 plus 12-6
         ldomixrule = .true.
         ldolshift = .true.
         ldoltailc = .true.
         ldormin = .true.
         ldorcut = .true.
         ldorcutin = .true.
      elseif ( cpnumber .eq. CP_STILLINGER_WEBER 
     &        .or. cpnumber .eq. CP_SW_PAIR_ONLY ) then
c        --- stillinger-weber potential
         ldomixrule = .true.
      elseif ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD .or.
     &        cpnumber .eq. CP_EAM_PAIR_ONLY ) then
c        --- EAM potential
c        --- set mixrule to explicit
         call twh_classical_mixrule(GLB_SET,'Explicit')
         ldointerpolatestyle = .true.
         ldorcut = .true.
c        --- allocate the arrays for the EAM potential
c        --- maxchain is the number of chains + 1
         call twh_nchain(GLB_GET,maxchain)
         maxchain = maxchain + 1
         call twh_allocate_eam(maxchain)
      elseif ( cpnumber .eq. CP_TABULATED_PAIR ) then
c        --- tabulated pair potential
c        --- set mixrule to explicit
         call twh_classical_mixrule(GLB_SET,'Explicit')
         ldointerpolatestyle = .true.
      elseif ( cpnumber .eq. CP_HARD_2580_MULTISTEP ) then
c        --- inherently short ranged potential
         ldomixrule = .true.
         ldorpd = .true.
      elseif ( cpnumber .eq. CP_REPULSIVE_2580_MULTISTEP ) then
c        --- inherently short ranged potential
         ldomixrule = .true.
      elseif ( cpnumber .eq. CP_BUFFERED_14_7 ) then
c        --- buffered 14-7
         ldomixrule = .true.
         ldoltailc = .true.
         ldormin = .true.
         ldorcut = .true.
         ldorcutin = .true.
      elseif ( cpnumber .eq. CP_DOUBLE_EXPONENTIAL ) then
c        --- double exponential
         ldomixrule = .true.
         ldoltailc = .true.
         ldolshift = .true.
         ldormin = .true.
         ldorcut = .true.
         ldorcutin = .true.
      else
         write(6,*) 'SETPOTENTIALLOGIC: Invalid cpnumber',cpnumber
         lfailure = .true.
         return
      endif

      return
      end
