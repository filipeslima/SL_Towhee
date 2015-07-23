#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2011 Marcus G. Martin                       *
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
      subroutine twh_writetorsion(index)
c     ******************************************************************
c     * writes all of the torsion parameters used in this simulation to*
c     * the file index                                                 *
c     * index: file index for output                                   *
c     *                                                                *
c     * originally written 02-22-2002 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LOFTOR
#define FUNCTION_GET_LTORHERE
#define FUNCTION_GET_NTORLOOP
#define FUNCTION_GET_TORCOEFF
#define FUNCTION_GET_TORFF
#define FUNCTION_GET_TORSTRING
#define FUNCTION_GET_TORSTYLE
#include "functions.h"

c     --- variables passed to/from the subroutine
      integer index

c     --- local variables
      logical lfound
      integer type,i,itemp,tstyle,tloop
      double precision radtodeg

c     --- retrieve constants
      call twh_constant_radtodeg(GLB_GET,radtodeg)

      write(index,'(a13)') 'Torsion Types'
      lfound = .false.
      do type = 1,TTORMAX
         if ( twh_get_ltorhere(type) ) then
            lfound = .true.
            tstyle = twh_get_torstyle(type)
            if ( tstyle .eq. 1 ) then
               write(index,10) type,twh_get_torcoeff(type,1)
     &              ,twh_get_torcoeff(type,0)
 10            format('Type: ',i3,1x,'Style: Harmonic  Angle:'
     &              ,f6.4,1x,'k ',f10.1)
            elseif ( tstyle .eq. 2 ) then
               write(index,20) type,(twh_get_torcoeff(type,i),i=1,3)
 20            format('Type: ',i3
     &              ,1x,'Style: OPLS Cosine Series',/
     &              ,10x,'k1:',f10.1,1x,'k2:',f10.1,1x,'k3:'
     &              ,f10.1)
            elseif ( tstyle .eq. 3 ) then
               write(index,30) type,twh_get_ntorloop(type)
 30            format('Type: ',i3,1x
     &              ,'Style: Amber/Charmm/Gromos Cosine Series'
     &              ,' n: ',i3)
               do i = 1,twh_get_ntorloop(type)
                  itemp = 3*(i-1)
                  write(index,33) i,twh_get_torcoeff(type,itemp+1)
     &                 ,twh_get_torcoeff(type,itemp+2)
     &                 ,twh_get_torcoeff(type,itemp+3)
               enddo
 33            format(10x,'Index',i3,3x,' k:',f10.1,' m:',f3.1
     &              ,' d:',f10.7)
            elseif ( tstyle .eq. 4 ) then
               write(index,40) type,twh_get_ntorloop(type)
 40            format('Type: ',i3,1x
     &              ,'Style: Amber/Charmm/Gromos Cosine Series'
     &              ,' n: ',i3)
               do i = 1,twh_get_ntorloop(type)
                  itemp = 3*(i-1)
                  write(index,33) i,twh_get_torcoeff(type,itemp+1)
     &                 ,twh_get_torcoeff(type,itemp+2)
     &                 ,twh_get_torcoeff(type,itemp+3)
               enddo
               itemp = 3*twh_get_ntorloop(type)
               write(index,41) twh_get_torcoeff(type,itemp+2)*radtodeg
     &              ,twh_get_torcoeff(type,itemp+1)
 41            format('With additional Harmonic   Angle:'
     &              ,f8.4,1x,'k ',f10.1)
            elseif ( tstyle .eq. 5 ) then
               write(index,62) type,(twh_get_torstring(type,i),i=1,4)
     &              ,(twh_get_torcoeff(type,i),i=0,19)
            elseif ( tstyle .eq. 6 ) then
               write(index,63) type,(twh_get_torcoeff(type,i),i=0,2)
               
            elseif ( tstyle .eq. 7 ) then
               write(index,70) type,(twh_get_torcoeff(type,i),i=0,1)
 70            format('Type: ',i3,1x,'Style: TraPPE simple cosine '
     &              ,/,10x,'e0 ',f10.2,'e1 ',f10.2)
            elseif ( tstyle .eq. 8 ) then
               write(index,73) type,twh_get_torff(type)
               
            elseif ( tstyle .eq. 9 ) then
c               --- class 2 with explicit terms
               write(index,91) type,(twh_get_torcoeff(type,i),i=0,31)
               
            elseif ( tstyle .eq. 10 ) then
               write(index,100)type,twh_get_ntorloop(type)
 100           format('Type: ',i3,1x,'Style: Cosine Power Series '
     &              ,'Sum x=1,Ntor a(x) Cos(phi)^x'
     &              ,/,10x,'Ntor:',1x,i3)
               do i=0,twh_get_ntorloop(type)
                  write(index,101) i,twh_get_torcoeff(type,i)
               enddo
 101           format (10x,'a(',i1,'):',1x,g12.6)
            elseif ( tstyle .eq. 11 ) then
               write(index,50) type,(twh_get_torcoeff(type,i),i=0,3)
            elseif ( tstyle .eq. 12 ) then
               write(index,'(a5,1x,i3,1x,a28,1x,a2,1x,i3)') 
     &              'Type:',type
     &              ,'Style: Sum2003 Cosine Series','n:'
     &              ,twh_get_ntorloop(type)
               write(index,'(a34,1x,f6.2)') 
     &              'Torsion shift angle [torcoeff(0)]:'
     &              ,twh_get_torcoeff(type,0)*radtodeg
               do i = 0,twh_get_ntorloop(type)
                  itemp = i+1
                  write(index,'(a2,i3,1x,a9,f15.5)')
     &                 'n:',i,'Constant:',twh_get_torcoeff(type,itemp)
               enddo
            elseif ( tstyle .eq. 13 ) then
               write(index,'(a5,1x,i3,1x,a27)') 'Type:',type
     &              ,'Style: 2 term Cosine Series'
               do i = 1,2
                  write(index,'(a9,i1,a2,f15.5)') 
     &                 'twh_get_torcoeff(',i,'):'
     &                 ,twh_get_torcoeff(type,i)
               enddo
            elseif ( tstyle .eq. 14 ) then
               write(index,140) 'Type:',type
     &              ,'Style: UFF K[1 - Cn Cos(n phi)]'
     &              ,'K:',twh_get_torcoeff(type,1),'Cn:'
     &              ,twh_get_torcoeff(type,2)
     &              ,'n:',twh_get_torcoeff(type,3)
 140           format(a5,1x,i4,1x,a31,/,11x,a2,1x,f10.3,1x,a3,f10.3,1x
     &              ,a2,1x,f5.2)
            elseif ( tstyle .eq. 15 ) then
               write(index,150) 'Type:',type
     &              ,'Style: K[1 - Cos(n (phi - phi_0))]'
     &              ,'K:',twh_get_torcoeff(type,1),'n:'
     &              ,twh_get_torcoeff(type,2)
     &              ,'phi_0:',twh_get_torcoeff(type,3)
 150           format(a5,1x,i4,1x,a34,/,11x,a2,1x,f10.3,1x,a2,f10.3,1x
     &              ,a6,1x,f10.3)
            elseif ( tstyle .eq. 16 ) then
               write(index,160) 'Type:',type
     &              ,'Style: K[1 - Cos(2 phi)]','K:'
     &              ,twh_get_torcoeff(type,1)
 160           format(a5,1x,i4,1x,a24,/,10x,a2,1x,f10.3)
            elseif ( tstyle .eq. 17 ) then
               write(index,170) type,'c1:',twh_get_torcoeff(type,1)
     &              ,'c2:',twh_get_torcoeff(type,2),'f:'
     &              ,twh_get_torcoeff(type,3)
 170           format('Type:',i4,1x
     &              ,'Style: c1[1 + Cos(phi + f)] + c2[1 - Cos^2(phi)]'
     &              ,/,10x,2(a3,1x,f10.3,1x),a2,1x
     &              ,f10.5)
            elseif ( tstyle .eq. 18 ) then
               write(index,180) type,twh_get_torcoeff(type,1)
     &              ,twh_get_torcoeff(type,2)
     &              ,twh_get_torcoeff(type,3)
 180           format('Type:',i4,1x,'Style: Square Well',/,10x
     &              ,'Low Wall:',1x,f10.3,1x,'High Wall:',1x,f10.3,1x
     &              ,'Well Depth:',f10.3)
            elseif ( tstyle .eq. 19 ) then
               write(index,190) type,twh_get_ntorloop(type)
 190           format('Type: ',i3,1x
     &              ,'Style: Amber Cosine Series / Total Torsions'
     &              ,' n: ',i3)
               do i = 1,twh_get_ntorloop(type)
                  itemp = 3*(i-1)
                  write(index,33) i,twh_get_torcoeff(type,itemp+1)
     &                 ,twh_get_torcoeff(type,itemp+2)
     &                 ,twh_get_torcoeff(type,itemp+3)
               enddo
            elseif ( tstyle .eq. 20 ) then
               write(index,200) type,(twh_get_torcoeff(type,i),i=1,4)
 200           format('Type: ',i3
     &              ,1x,'Style: OPLS 4 Parameter Cosine Series',/
     &              ,10x,'k1:',f10.1,1x,'k2:',f10.1,1x,'k3:'
     &              ,f10.1,1x,'k4:',f10.1)
            elseif ( tstyle .eq. 21 ) then
               tloop = twh_get_ntorloop(type)
               write(index,210) type,tloop
 210           format('Type:',1x,i3,1x,'Style: Multiple Rigid Torsions'
     &              ,/,10x,'Number of allowed Torsions:',1x,i3)
               do i = 1,tloop
                  write(index,211) i,twh_get_torcoeff(type,i)
               enddo
 211           format(10x,'Torsion(',i1,'):',1x,f10.5)
            elseif ( tstyle .eq. 22 ) then
               write(index,220) type,(twh_get_torcoeff(type,i),i=0,5)
 220           format('Type:',1x,i3,1x,'Style: Siepmann Fluoroalkane'
     &              ,/,'a0 + a1(1-cos(phi)) + a2(1-Cos(3 phi))',1x
     &              ,'+ a3(1-Cos(phi))^5 + b0 exp[-b1 phi^2]'
     &              ,/,10x,'a0:',f10.1,1x,'a1:',f10.1,1x,'a2:',f10.1
     &              ,1x,'a3:',f10.1
     &              ,/,10x,'b0:',f10.1,1x,'b1:',f10.1)
            elseif ( tstyle .eq. 23 ) then
               write(index,230) type,(twh_get_torcoeff(type,i),i=0,1)
 230           format('Type:',1x,i3,1x,'Style: k1 (1 - Cos(phi - k0))'
     &              ,/,'k0:',1x,g12.6
     &              ,/,'k1:',1x,g12.6)
            else
               write(6,*) 'nondefined torsion style of '
     &              ,tstyle
     &              ,'for torsion number ',type
               stop
            endif
            if ( twh_get_loftor(type) ) then
               write(6,11) twh_get_torcoeff(type,-1)
            endif
 11         format(10x,'with 1-4 vdw and scaled (',f5.3,') 1-4 coulomb')
         endif
      enddo
      if ( .not. lfound ) write(index,*) '  No Torsion Types'

c     --- format statements
      
 50   format('Type: ',i3
     &     ,1x,'Style: Old UA OPLS Cosine Series',/
     &     ,10x,'k0:',f10.1,1x,'k1:',f10.1,1x,'k2:'
     &     ,f10.1,1x,'k3:',f10.1)
 62   format('Type: ',i3,1x,'Style: Compass with Class 2 cross terms'
     &     ,/,10x,4(a,1x)
     &     ,/,10x,'Normal Torsion Constants',3f10.2
     &     ,/,10x,'Bond(ij) - Torsion Constants',3f10.2
     &     ,/,10x,'Bond(jk) - Torsion Constants',3f10.2
     &     ,/,10x,'Bond(kl) - Torsion Constants',3f10.2
     &     ,/,10x,'Angle(ijk) - Torsion Constants',3f10.2
     &     ,/,10x,'Angle(jkl) - Torsion Constants',3f10.2
     &     ,/,10x,'Angle(ijk) - Angle(jkl) - Torsion Constant',f10.2
     &     ,/,10x,'Bond(ij) - Bond(kl) Constant',f10.2
     &     )
 63   format('Type: ',i3,1x,'Style: Compass without cross terms'
     &     ,/,10x,'Normal Torsion Constants',3f10.2)
 73   format('Type: ',i3,1x,'Style: ',a10,' Nonbonded terms only')

 91   format('Type: ',i3,1x,'Style: Class 2 with explicit cross terms'
     &     ,/,10x,'Normal Torsion Constants',6f10.2
     &     ,/,10x,'Bond(jk) - Torsion Constants',4f10.2
     &     ,/,10x,'Bond(ij) - Torsion Constants',4f10.2
     &     ,/,10x,'Bond(kl) - Torsion Constants',4f10.2
     &     ,/,10x,'Angle(ijk) - Torsion Constants',4f10.2
     &     ,/,10x,'Angle(jkl) - Torsion Constants',4f10.2
     &     ,/,10x,'Angle(ijk) - Angle(jkl) - Torsion Constants',3f10.2
     &     ,/,10x,'Bond(ij) - Bond(kl) Constants',3f10.2
     &     )

      return
      end
    
