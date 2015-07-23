#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2014 Marcus G. Martin                       *
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
      subroutine twh_writebond(index)
c     ******************************************************************
c     * writes all of the bond parameters used in this simulation to   *
c     * the file index                                                 *
c     * index: file index for output                                   *
c     *                                                                *
c     * originally written 02-22-2002 by M.G. Martin                   *
c     * last modified 07-18-2014 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LVIBHERE
#define FUNCTION_GET_VIBCOEFF
#include "functions.h"

c     --- variables passed to/from the subroutine
      integer index

c     --- local variables
      logical lfound
      integer type,vstyle

      write(index,'(a10)') 'Bond Types'

      lfound = .false.
      do type = 1,TVIBMAX
         if ( twh_get_lvibhere(type) )  then
            lfound = .true.
            call twh_vibstyle(GLB_GET,type,vstyle)
            if ( vstyle .eq. BOND_STYLE_FIXED ) then
               write(index,10) type,twh_get_vibcoeff(type,0)
 10            format('Type: ',i6,' Style: Fixed  Length: ',f6.4)
            elseif ( vstyle .eq. BOND_STYLE_HARMONIC ) then
               write(index,20) type,twh_get_vibcoeff(type,0)
     &              ,twh_get_vibcoeff(type,1)
 20            format('Type: ',i6,1x,'Style: Standard Harmonic',1x
     &              ,'Length: ',f6.4,1x,'Constant: ',f10.1)
            elseif ( vstyle .eq. BOND_STYLE_GROMOS ) then
               write(index,30) type,twh_get_vibcoeff(type,0)
     &              ,twh_get_vibcoeff(type,1)
 30            format('Type: ',i6,1x,'Style: Gromos  Length: ',f6.4
     &              ,1x,'Constant: ',f10.1)
            elseif ( vstyle .eq. BOND_STYLE_NONLINEAR ) then
c              --- nonlinear
               write(index,40) type,twh_get_vibcoeff(type,0)
     &              ,twh_get_vibcoeff(type,1)
     &              ,twh_get_vibcoeff(type,2)
 40            format('Type: ',i6,1x,'Style: Nonlinear Length: ',f6.4
     &              ,1x,'k1: ',f10.3,'k2: ',f10.3)
            elseif ( vstyle .eq. BOND_STYLE_MMTWO ) then
c              --- MM2
               write(index,50) type,twh_get_vibcoeff(type,0)
     &              ,twh_get_vibcoeff(type,1)
 50            format('Type: ',i6,1x,'Style: MM2 Triplet Length: ',f6.4
     &              ,1x,'Constant: ',f10.3)
            elseif ( vstyle .eq. BOND_STYLE_COMPASS ) then
               write(index,60) type,twh_get_vibcoeff(type,0)
     &              ,twh_get_vibcoeff(type,1)
     &              ,twh_get_vibcoeff(type,2),twh_get_vibcoeff(type,3)
 60            format('Type: ',i6,1x,'Style: Quartic   Length: ',f6.4,/
     &              ,10x,'Constants:',3(1x,f10.1))
            elseif ( vstyle .eq. BOND_STYLE_NONBOND ) then
               write(index,70) type
 70            format('Type: ',i6,1x,'Style: Nonbonded Interactions')
            elseif ( vstyle .eq. BOND_STYLE_NONE ) then
               write(index,80) type
 80            format('Type: ',i6,1x,'Style: No Interactions')
            elseif ( vstyle .eq. BOND_STYLE_MORSE ) then
               write(index,90) type,twh_get_vibcoeff(type,0)
     &              ,twh_get_vibcoeff(type,1)
     &              ,twh_get_vibcoeff(type,2)
 90            format('Type: ',i6,1x
     &              ,'Style: Morse  u(r) = k1 * {exp[k2(r-d)] - 1}^2',/
     &              ,10x,'d:',1x,f6.4,1x,'k1:',1x,f10.1,1x,'k2:'
     &              ,1x,f10.1)
            elseif ( vstyle .eq. BOND_STYLE_SQUARE_WELL ) then
               write(index,100) type,twh_get_vibcoeff(type,1)
     &              ,twh_get_vibcoeff(type,2)
     &              ,twh_get_vibcoeff(type,3)
 100           format('Type: ',i6,1x,'Style: Infinte Square Well',/,10x
     &              ,'Minimum: ',f10.3,1x,'Maximum:',1x,f10.3,1x
     &              ,'Energy:',f10.5)
            elseif ( vstyle .eq. BOND_STYLE_HARMONIC_NONBOND ) then
               write(index,110) type,twh_get_vibcoeff(type,0)
     &              ,twh_get_vibcoeff(type,1)
 110            format('Type: ',i6,1x
     &              ,'Style: Standard Harmonic plus Nonbond',1x
     &              ,'Length: ',f6.4,1x,'Constant: ',f10.1)
            elseif ( vstyle .eq. BOND_STYLE_FENE ) then
c              --- Frank Heilmann, BASF AG
c              --- put out FENE-related stuff
               write(index,120) type,twh_get_vibcoeff(type,1)
     &              ,twh_get_vibcoeff(type,2)
 120            format('Type: ',i6,1x
     &              ,'Style: FENE for CG models',/,10x
     &              ,'R0: ',f6.4,2x,'k: ',f10.3,f10.3)
            else
               write(index,*) 'WRITEBOND: unknown vibstyle',vstyle
               stop
            endif
         endif
      enddo
      if ( .not. lfound ) write(index,*) '  No Bond Types'

      return
      end
