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
      subroutine twh_writeangle(index)
c     ******************************************************************
c     * writes all of the angle parameters used in this simulation to  *
c     * the file index                                                 *
c     * index: file index for output                                   *
c     *                                                                *
c     * originally written 02-22-2002 by M.G. Martin                   *
c     * last modified 07-21-2014 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_ANGLESTYLE
#define FUNCTION_GET_BENCOEFF
#define FUNCTION_GET_LBENDHERE
#define FUNCTION_GET_LBONANG
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer index

c     --- local variables
      logical lfound
      integer type,style
      character*(FFNAMELEN) namea,nameb
      character*15 order
      double precision radtodeg

c     --- retrieve constant
      call twh_constant_radtodeg(GLB_GET,radtodeg)

      write(index,'(a11)') 'Angle Types'
      lfound = .false.
      do type = 1,TBENMAX
         if ( twh_get_lbendhere(type) )  then
            lfound = .true.
            style = twh_get_anglestyle(type)
            if ( style .eq. ANG_STYLE_RIGID ) then
               write(index,10) type,radtodeg*twh_get_bencoeff(type,0)
 10            format('Type: ',i3,1x,'Style: Fixed Angle  Angle: '
     &              ,f10.3)
            elseif ( style .eq. ANG_STYLE_HARMONIC ) then
               write(index,11) type,radtodeg*twh_get_bencoeff(type,0)
     &              ,twh_get_bencoeff(type,1)
 11            format('Type: ',i3,1x,'Style: Standard Harmonic'
     &              ,'  Angle: ',f10.3,1x,'Constant: ',f10.1)
            elseif ( style .eq. ANG_STYLE_DREIDING ) then
               write(index,21) type,radtodeg*twh_get_bencoeff(type,0)
     &              ,twh_get_bencoeff(type,1)
 21            format('Type: ',i3,1x,'Style: 1 - Cosine'
     &              ,'  Angle: ',f10.3,1x,'Constant: ',f10.1)
            elseif ( style .eq. ANG_STYLE_HARM_COS ) then
               write(index,31) type,radtodeg*twh_get_bencoeff(type,0)
     &              ,twh_get_bencoeff(type,1)
 31            format('Type: ',i3,1x,'Style: Harmonic Cosine Angle: '
     &              ,f10.3,1x,'Constant: ',f10.1)
            elseif ( style .eq. ANG_STYLE_COMPASS_AUTO ) then
               call twh_bendnames(GLB_GET,type,1,1,namea)
               call twh_bendnames(GLB_GET,type,1,3,nameb)
               write(index,41) type,radtodeg*twh_get_bencoeff(type,0)
     &              ,twh_get_bencoeff(type,1),twh_get_bencoeff(type,2)
     &              ,twh_get_bencoeff(type,3),namea
     &              ,twh_get_bencoeff(type,4),nameb
     &              ,twh_get_bencoeff(type,5),twh_get_bencoeff(type,6)
 41            format('Type: ',i3,1x,'Style: Compass  Angle: ',f10.3
     &              ,/,10x,'Constants: ',3f10.1,/,11x,a,f10.3,2x,a,f10.3
     &              ,/,10x,'Bond-Bond Constant ',f10.3)
            elseif ( style .eq. ANG_STYLE_CHARMM_UB ) then
               write(index,51) type,radtodeg*twh_get_bencoeff(type,0)
     &              ,twh_get_bencoeff(type,1),twh_get_bencoeff(type,2)
     &              ,twh_get_bencoeff(type,3)
 51            format('Type: ',i3,1x,'Style: Charmm   Angle: ',f10.3
     &              ,1x,'Constant: ',f10.1,/,10x,'UB Length: ',f6.4
     &              ,1x,'Constant: ',f10.1)
            elseif ( style .eq. ANG_STYLE_NONBOND ) then
               write(index,61) type
 61            format('Type: ',i3,1x,'Style: Nonbonded Interaction')
            elseif ( style .eq. ANG_STYLE_NB_HARM ) then
               write(index,71) type,radtodeg*twh_get_bencoeff(type,0)
     &              ,twh_get_bencoeff(type,1)
 71            format('Type: ',i3,1x,'Style: Harmonic plus Nonbonded '
     &              ,'Interaction'
     &              ,/,10x,'Angle: ',f10.3,'Constant: ',f10.1)
            elseif ( style .eq. ANG_STYLE_COMPASS_EX ) then
               write(index,81) type,radtodeg*twh_get_bencoeff(type,0)
     &              ,twh_get_bencoeff(type,1),twh_get_bencoeff(type,2)
     &              ,twh_get_bencoeff(type,3),twh_get_bencoeff(type,4)
     &              ,twh_get_bencoeff(type,5),twh_get_bencoeff(type,6)
     &              ,twh_get_bencoeff(type,7),twh_get_bencoeff(type,8)
     &              ,twh_get_bencoeff(type,9),twh_get_bencoeff(type,10)
 81            format('Type: ',i3,1x,'Style: Compass  Angle: ',f10.3
     &              ,/,10x,'Constants: ',3f10.1
     &              ,/,10x,'Bond-Angle Constants '
     &              ,'k1 :',f10.3,' d1:',f10.3,' k2 :',f10.3
     &              ,' d2:',f10.3,/
     &              ,10x,'Bond-Bond Constants k: ',f10.3,' d1:',f10.3
     &              ,' d2:',f10.3)
            elseif ( style .eq. ANG_STYLE_FOURIER ) then
               write(index,91)'Type:',type
     &              ,'Style: UFF One Term Fourier Expansion'
     &              ,'K (1 + m Cos(n theta))'
     &              ,'Angle:'
     &              ,radtodeg*twh_get_bencoeff(type,0)
     &              ,'n:',twh_get_bencoeff(type,1)
     &              ,'m:',twh_get_bencoeff(type,2)
 91            format(a5,1x,i4,1x,a37,1x,a22,/
     &              ,11x,a6,1x,f10.3,1x,a2,1x,f5.2,1x,a2,1x,f4.1)
               call writeanglek(type)
            elseif ( style .eq. ANG_STYLE_TTFOURIER ) then
               write(index,101) 'Type:',type
     &              ,'Style: UFF Three Term Fourier Expansion','Angle:'
     &              ,radtodeg*twh_get_bencoeff(type,0)
     &              ,'C0:',twh_get_bencoeff(type,1)
     &              ,'C1:',twh_get_bencoeff(type,2)
     &              ,'C2:',twh_get_bencoeff(type,3)
 101           format(a5,1x,i4,1x,a39,/,11x,a6,1x,f7.3
     &              ,3(1x,a3,1x,f10.3))
               call writeanglek(type)
            elseif ( style .eq. ANG_STYLE_NONE ) then
               write(index,111) type
 111            format('Type: ',i3,1x,'Style: No Interaction')
            elseif ( style .eq. ANG_STYLE_SIXTH_AUTO ) then
               call twh_bendnames(GLB_GET,type,1,1,namea)
               call twh_bendnames(GLB_GET,type,1,3,nameb)
               write(index,121) type,radtodeg*twh_get_bencoeff(type,0)
     &              ,twh_get_bencoeff(type,1),twh_get_bencoeff(type,2)
     &              ,twh_get_lbonang(type)
     &              ,namea,nameb
     &              ,twh_get_bencoeff(type,3)
 121           format('Type: ',i3,1x,'Style: MM2  Angle: ',f10.3,/
     &              ,10x,'Constants: ',2f10.3,/
     &              ,10x,'Angle-Bond Logical',1x,l1,1x
     &              ,'Atom A:',1x,a10,1x,'Atom B:',1x,a10,1x,/
     &              ,10x,'Angle-Bond Constant:',f10.3)
            elseif ( style .eq. ANG_STYLE_SQUARE_WELL ) then
               write(index,131) type,radtodeg*twh_get_bencoeff(type,0)
     &              ,twh_get_bencoeff(type,1),twh_get_bencoeff(type,2)
     &              ,twh_get_bencoeff(type,3)
 131           format('Type:',1x,i3,1x,'Style: Infinite Square Well'
     &              ,1x,'Angle:',1x,f10.3,/
     &              ,10x,'Minimum Distance:',1x,f6.3
     &              ,1x,'Maximum Distance:',1x,f6.3,1x,'Well Depth:'
     &              ,1x,f10.5)
            elseif ( style .eq. ANG_STYLE_MULTI_RIGID ) then
               write(index,141) type,radtodeg*twh_get_bencoeff(type,0)
     &              ,radtodeg*twh_get_bencoeff(type,1)
     &              ,twh_get_bencoeff(type,2)
 141           format('Type:',1x,i3,1x,'Style: Multiple Rigid',/
     &              ,10x,'Angle(1):',1x,f10.3,/
     &              ,10x,'Angle(2):',1x,f10.3,/
     &              ,10x,'Tolerance:',1x,e10.5)
            elseif ( style .eq. ANG_STYLE_MMFF ) then
               write(index,151) type,radtodeg*twh_get_bencoeff(type,0)
     &              ,twh_get_bencoeff(type,1),twh_get_bencoeff(type,2)
 151           format('Type:',1x,i4,1x,'Style: MMFF Cubic',/
     &              ,10x,'Angle:',1x,f10.3,1x,/
     &              ,10x,'Constants:',2(1x,f10.3))
               write(index,152) twh_get_lbonang(type)
 152           format(10x,'Bond-Angle Logical',1x,l1)
               if ( twh_get_lbonang(type) ) then
                  write(index,153) twh_get_bencoeff(type,3)
     &                 ,twh_get_bencoeff(type,4)
 153              format(10x,'Bond-Angle Constants:',2(1x,f10.2))
               endif
               call twh_angleorder(GLB_GET,type,order)
               write(index,154) order
 154           format(10x,'Angle Order:',1x,a15)
            elseif ( style .eq. ANG_STYLE_HARM_COS_OT ) then
               write(index,161) type,radtodeg*twh_get_bencoeff(type,0)
     &              ,twh_get_bencoeff(type,1)
 161           format('Type: ',i3,1x
     &              ,'Style: Harmonic Cosine Angle plus 1-3 nonbond: '
     &              ,f10.3,1x,'Constant: ',f10.1)
            else
               write(index,*) 'WRITEANGLE: angle type ',type,
     &              'has an invalid anglestyle',style
               stop
            endif
         endif
      enddo

      if ( .not. lfound ) write(index,*) '  No Angle Types'

      return
      end
    
      subroutine writeanglek(type)
c     *****************************************************************
c     * writes out the energetic parameters for the angles            *
c     *                                                               *
c     * last modified 10-13-2008 by M.G. Martin                       *
c     *****************************************************************
      implicit none
#define FUNCTION_GET_BENPREFACT
#define FUNCTION_GET_IJBEN2
#define FUNCTION_GET_IJBEN3
#define FUNCTION_GET_INBEN
#define FUNCTION_GET_ITBEN
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer type
c     --- local variables 
      integer imolty,iunit,ibend,nmolty,nunit

c     --- the prefactor for this depend upon the actual atoms
c     --- run through all of the molecules and print out the
c     --- constants for each instance
      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            do ibend = 1,twh_get_inben(imolty,iunit)
               if ( twh_get_itben(imolty,iunit,ibend) .eq. type ) then
c                 --- output the atoms and prefactor
                  write(6,102) 'Molecule:',imolty,'Atoms:',
     &                 iunit,twh_get_ijben2(imolty,iunit,ibend)
     &                 ,twh_get_ijben3(imolty,iunit,ibend)
     &                 ,'K:',twh_get_benprefact(imolty,iunit,ibend)
 102              format(11x,a9,1x,i3,1x,a6,3(1x,i4),1x,a2,1x
     &                 ,f10.3)
               endif
            enddo
         enddo
      enddo
      return
      end
