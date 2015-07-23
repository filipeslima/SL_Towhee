#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2000-2014 Marcus G. Martin                       *
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
      function twh_vbond(lfailure,type,length,imolty,iunit,junit
     &     ,lovrlap )
c     ******************************************************************
c     * Computes the bond energy for a given type and length           *
c     * type: the forcefield number of the bond vibration              *
c     * length: the distance between the two atoms                     *
c     * imolty: the molecule type - used for nonbonded vib style       *
c     * iunit: the first unit number - used for nonbonded vib style    *
c     * junit: the second unit number - used for nonboned vib style    *
c     * lovrlap: returned as true if there is a hard overlap in the    *
c     *          bond interaction - this occurs with nonbonded vib     *
c     *          style, or if a fixed bond length is > 1% from its mean*
c     *                                                                *
c     * originally written 1-24-2000 M.G. Martin                       *
c     * last modified 07-17-2014 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_GT
#define FUNCTION_CMP_LT
#define FUNCTION_EXPON
#define FUNCTION_GET_VIBCOEFF
#define FUNCTION_VBOND
#include "functions.h"
c     --- variables passed to/from function
      logical lfailure,lovrlap
      integer type,imolty,iunit,junit
      double precision length
c     --- local variables
      integer style
c     --- double precision scalars
      double precision dist,distsq,max_bond_length,onethird

c     --- retrieve constants
      call twh_vibstyle(GLB_GET,type,style)
      call twh_max_bond_length(GLB_GET,max_bond_length)

c     --- for all bond lengths, check against maximum bond length
      if ( length .gt. max_bond_length ) then
         twh_vbond=0.0d0
         lovrlap = .true.
         return
      endif

c     --- for valid bond lengths, check the style to determine the
c     --- bond energy
      if ( style .eq. BOND_STYLE_FIXED ) then
c        --- fixed bond length
c        --- zero bond stretch energy unless outside tolerance
         dist = dabs( (length-twh_get_vibcoeff(type,0))
     &        /twh_get_vibcoeff(type,0)) 
         if ( twh_cmp_gt(dist, 0.01d0) ) then
c           --- this fixed bond length differs too much from its 
c           --- target value
            lovrlap = .true.
         endif
         twh_vbond = 0.0d0
      elseif (style .eq. BOND_STYLE_HARMONIC ) then
c        --- standard harmonic
         twh_vbond = twh_get_vibcoeff(type,1) 
     &        * (length - twh_get_vibcoeff(type,0))**2
      elseif ( style .eq. BOND_STYLE_GROMOS ) then
c        --- GROMOS bond stretch potential
         twh_vbond = twh_get_vibcoeff(type,1)*(length*length 
     &        - twh_get_vibcoeff(type,0)*twh_get_vibcoeff(type,0) )**2

      elseif ( style .eq. BOND_STYLE_NONLINEAR ) then
c        --- nonlinear (same as LAMMPS 2001 implementation)
c        --- LAMMPS references Rector, van Swol, Henderson
c        --- Molecular Physics 82 1009 (1994)
c        --- U = twh_get_vibcoeff(1) * (r-r0)^2 
c        --- / ( twh_get_vibcoeff(2) - (r-r0)^2 )
         dist = (length - twh_get_vibcoeff(type,0))**2
         twh_vbond = twh_get_vibcoeff(type,1)*dist 
     &        / (twh_get_vibcoeff(type,2)-dist)

      elseif ( style .eq. BOND_STYLE_MMTWO ) then
c        --- MM2 style
c        ---- this bond potential is not well formed and has a 
c        --- sensible minimum when the length = vibcoeff(0), but then
c        --- has a nonsense maximum when the length = vibcoeff + 1/3
c        --- and then decreases in energy towards -infinity afterwards
c        --- In order to make this sensible enforce a +infinity energy
c        --- beyond that distance
         dist = length - twh_get_vibcoeff(type,0)
         call twh_constant_onethird(GLB_GET,onethird)
         if ( dist .gt. onethird ) then
c           --- infinite energy
            twh_vbond=0.0d0
            lovrlap = .true.
         else
            twh_vbond = twh_get_vibcoeff(type,1)*dist*dist
     &           *(1.0d0 - 2.0d0*dist)
         endif
      elseif ( style .eq. BOND_STYLE_COMPASS ) then
c        --- Compass quartic potential
c        --- also used by MMFF94
         dist = length - twh_get_vibcoeff(type,0)
         distsq = dist*dist
         twh_vbond = twh_get_vibcoeff(type,1)*(distsq)
     &        + twh_get_vibcoeff(type,2)*(dist*distsq)
     &        + twh_get_vibcoeff(type,3)*(distsq*distsq)
      elseif ( style .eq. BOND_STYLE_NONBOND ) then
c        --- nonbonded interactions used as a bonding potential
c        --- need to initialize bond energy and set logical
         call vnonbondinbond(lfailure,lovrlap,imolty,iunit,junit
     &        ,length,twh_vbond)
         if ( lfailure .or. lovrlap ) return
      elseif ( style .eq. BOND_STYLE_NONE ) then
c         --- no interaction
         twh_vbond = 0.0d0
      elseif ( style .eq. BOND_STYLE_MORSE ) then
c         --- Morse potential
         dist = twh_get_vibcoeff(type,2)
     &        *(length - twh_get_vibcoeff(type,0))
         twh_vbond = twh_get_vibcoeff(type,1)
     &        *((twh_expon(dist) - 1.0d0)**2)
      elseif ( style .eq. BOND_STYLE_SQUARE_WELL ) then
c         --- infinite square well
         if ( twh_cmp_lt(length, twh_get_vibcoeff(type,1)) .or.
     &        twh_cmp_gt(length, twh_get_vibcoeff(type,2)) ) then
c            --- outside of the square well so overlap
            twh_vbond=0.0d0
            lovrlap = .true.
         else
c           ---- inside the square well so assign energy
            twh_vbond = twh_get_vibcoeff(type,3)
         endif
      elseif ( style .eq. BOND_STYLE_HARMONIC_NONBOND ) then
c        --- harmonic plus nonbond
         call vnonbondinbond(lfailure,lovrlap,imolty,iunit,junit
     &        ,length,twh_vbond)
         if ( lfailure .or. lovrlap ) return
c        --- standard harmonic
         twh_vbond = twh_get_vibcoeff(type,1) 
     &        * (length - twh_get_vibcoeff(type,0))**2
      elseif ( style .eq. BOND_STYLE_FENE ) then
c        --- Frank Heilmann, BASF AG
c        --- FENE bond potential for CG models
c        --- similar to to P. Virnau's thesis:
c        --- V_fene = - 33.75 * eps * ln( 1 - ( r/(1.5*sig) )**2 )
         if ( length .le. twh_get_vibcoeff(type,1) ) then
c           --- nonbond portion
            call vnonbondinbond(lfailure,lovrlap,imolty,iunit,junit
     &           ,length,twh_vbond)
            if ( lfailure .or. lovrlap ) return
c           --- distance based FENE term
c           --- vibcoeff(3) = -0.5 * vibcoeff(2) * vibcoeff(1)^2
            twh_vbond = twh_vbond + twh_get_vibcoeff(type,3)
     &           * dlog(1.0d0 -(length/(twh_get_vibcoeff(type,1)))**2)
         else
c           --- infinite energy for overextended bond
            twh_vbond = 0.0d0
            lovrlap = .true.
         endif

      else
         write(6,*) 'vibration style ',style,' is unknown.',type
         twh_vbond=0.0d0
         lfailure = .true.
         return
      endif

c      write(6,*) 'iunit,junit,length,vbond',iunit,junit
c     &     ,length,twh_vbond/CON_KCALTOK

      return
      end

      subroutine vnonbondinbond(lfailure,lovrlap,imolty,iunit,junit
     &     ,length,vnonbond)
c     ******************************************************************
c     * computes nonbond terms in the bond energy                      *
c     *                                                                *
c     * split from vbond 07-18-2014 by M.G. Martin                     *
c     * last modified 07-18-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_QQATOM
#define FUNCTION_VTWOBODY
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure,lovrlap
      integer imolty,iunit,junit
      double precision length,vnonbond
c     --- local variables
      integer ntii,ntjj
      double precision qqfact,distsq,dv
c     --- double precision arrays
      double precision vf
      dimension vf(MAX_FOREIGN_LAMBDA)

c     --- compute the nonbonded energy between these atoms
      call twh_ntype(GLB_GET,imolty,iunit,ntii)
      call twh_ntype(GLB_GET,imolty,junit,ntjj)
      distsq = length*length
c     --- nonbond pair potential
      vnonbond = twh_vtwobody(lfailure,distsq,ntii,ntjj,0,lovrlap,vf,dv)
      if ( lfailure .or. lovrlap ) return
c     --- get the nonbonded coulombic energy (only the real space 
c     --- part)
      call twh_constant_qqfact(GLB_GET,qqfact)
      vnonbond = vnonbond + qqfact*twh_get_qqatom(imolty,iunit)
     &     *twh_get_qqatom(imolty,junit)/length

      return
      end
