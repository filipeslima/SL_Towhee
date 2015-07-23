#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2013 Marcus G. Martin                       *
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
      subroutine twh_writenonbond(ifile,lerror)
c     ******************************************************************
c     * writes all of the information about the nonbonded interactions *
c     * to the file index                                              *
c     * ifile: file index for output                                   *
c     * lerror: logical used if there is a problem                     *
c     *                                                                *
c     * originally split out of readtowhee 08-16-2004 by M.G. Martin   *
c     * last modified 07-17-2013 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_HBONDCOEFF
#define FUNCTION_GET_LHERE
#define FUNCTION_GET_NATIVE_LAMBDA
#define FUNCTION_GET_NBCOEFF
#define FUNCTION_GET_NBNAME
#define FUNCTION_GET_SCALING_STYLE
#define FUNCTION_LIMITED_DOUBLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lerror
      integer ifile
c     --- local variables
c     --- logical scalars
      logical ltesti,ltestj,lshift,ltailc
c     --- character string scalars
      character*20 eamdstyle,eampstyle,eamestyle
c     --- integer scalars
      integer itype,jtype,ij,itemp,icoeff,hbtype,cpnumber,iswmax
      integer potentialstyle,tabnpair
c     --- double precsion scalars
      double precision cvalue,dvalue,sfactor,tpairone,tpairtwo
      double precision rcutin

c     --- retrieve constants
      call twh_potentialstyle(GLB_GET,potentialstyle)
c     --- nothing to output for an external potentialstyle
      if ( potentialstyle .eq. POT_EXTERNAL ) return

c     --- get the classical potential number
      call twh_cp_number(GLB_GET,cpnumber)

c     --- output information about the nonbonded potential parameters
      write(ifile,'(a21)') 'Nonbonded Force Field'
      if ( cpnumber .eq. CP_9_6 ) then
         write(ifile,'(a27)') 'Lennard-Jones 9-6 potential'
      elseif ( cpnumber .eq. CP_12_6_PLUS_SOLVATION ) then
         write(ifile,'(a,1x,a)') 'Lennard-Jones 12-6 potential for EEF1'
     &        ,'implicit solvent'
      elseif ( cpnumber .eq. CP_12_6_PLUS_12_10_H_BOND ) then
         write(ifile,'(a)') '12-6 LJ plus 12-10 Hydrogen bond potential'
      elseif ( cpnumber .eq. CP_12_9_6 ) then
         write(ifile,'(a)') '12-9-6 potential'
      elseif ( cpnumber .eq. CP_BUFFERED_14_7 ) then
         write(ifile,'(a)') 'Buffered 14-7 potential'
      elseif ( cpnumber .eq. CP_DOUBLE_EXPONENTIAL ) then
         write(ifile,'(a)') 'Double Exponential'
      elseif ( cpnumber .eq. CP_EAM_PAIR_ONLY ) then
         write(ifile,'(a,1x,a)') 'Embedded-Atom pair term without'
     &        ,'the embeddeding term'
      elseif ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD ) then
         write(ifile,'(a)') 'Embedded-Atom potential'
      elseif ( cpnumber .eq. CP_EXPONENTIAL_6 ) then
         write(ifile,'(a)') 'Exponential-6 potential'
      elseif ( cpnumber .eq. CP_EXPONENTIAL_12_6 ) then
         write(ifile,'(a)') 'Exponential + LJ 12 + LJ 6 potential'
      elseif ( cpnumber .eq. CP_GORDON_N_6 ) then
         write(6,'(a)') 'Modified n-6 potential of Peter Gordon'
      elseif ( cpnumber .eq. CP_HARD_2580_MULTISTEP ) then
         write(ifile,'(a43)') 
     &        'Hard Core plus 1.2, 1.5, 1.8 and 2.0  Wells'
      elseif ( cpnumber .eq. CP_HARD_SPHERE ) then
         write(ifile,'(a)') 'Hard Sphere potential'
      elseif ( cpnumber .eq. CP_LENNARD_JONES ) then
         write(ifile,'(a)') 'Lennard-Jones 12-6 potential'
      elseif ( cpnumber .eq. CP_SCALED_LENNARD_JONES ) then
         if (scp_get_scaling_style() .eq. SCP_SC ) then
            write(ifile,'(a)') 'Soft Core Lennard Jones potential'
         elseif (scp_get_scaling_style() .eq. SCP_NONE ) then
            write(ifile,'(a)') 'Lennard-Jones 12-6 potential'
         else 
            write(ifile,*) 'WRITENONBOND: undefined rescaling style'
            lerror = .true.
            return
         endif
      elseif ( cpnumber .eq. CP_MULTIWELL ) then
         write(ifile,'(a)') 'Multiple square well potential'
      elseif ( cpnumber .eq. CP_REPULSIVE_2580_MULTISTEP ) then
         write(ifile,'(a48)') 
     &        'Repulsive Core plus 1.2, 1.5, 1.8 and 2.0  Wells'
      elseif ( cpnumber .eq. CP_REPULSIVE_SPHERE ) then
         write(ifile,'(a)') 
     &        'Repulsive Sphere potential, used for equilibration'
      elseif ( cpnumber .eq. CP_REPULSIVE_MULTIWELL ) then
         write(ifile,'(a)') 'Repulsive Multiwell'
      elseif ( cpnumber .eq. CP_REPULSIVE_WELL ) then
         write(ifile,'(a)') 'Repulsive well potential'
      elseif ( cpnumber .eq. CP_SQUARE_WELL ) then
         write(ifile,'(a)') 'Square well potential'
      elseif ( cpnumber .eq. CP_STILLINGER_WEBER ) then
         write(ifile,'(a)') 'Stillinger-Weber potential'
      elseif ( cpnumber .eq. CP_SW_PAIR_ONLY ) then
         write(ifile,'(a)') 'Stillinger-Weber using only the pair term'
      elseif ( cpnumber .eq. CP_TABULATED_PAIR ) then
         write(ifile,'(a)') 'Tabulated pair potential'
      elseif ( cpnumber .eq. CP_UFF_12_6 ) then
         write(6,'(a)') 'UFF 12-6 potential'
      elseif ( cpnumber .eq. CP_WEEKS_CHANDLER_ANDERSON ) then
         if (scp_get_scaling_style() .eq. SCP_WCA ) then
            write(ifile,'(a)') 'Weeks-Chandler-Anderson potential'
         elseif (scp_get_scaling_style() .eq. SCP_NONE ) then
            write(ifile,'(a)') 'Lennard-Jones 12-6 potential'
         else 
            write(ifile,*) 'WRITENONBOND: undefined rescaling style'
            lerror = .true.
            return
         endif
      else
         write(ifile,*) 'WRITENONBOND: undefined potential type'
         write(6,*) 'cpnumber:',cpnumber
         lerror = .true.
         return
      endif

      call twh_ltailc(GLB_GET,ltailc)
      if ( ltailc ) write(ifile,'(4x,a)') 'with tail corrections'
      call twh_lshift(GLB_GET,lshift)
      if ( lshift ) write(ifile,'(4x,a)') 
     &     'with a shift to zero at the cut-off'


      if ( cpnumber .eq. CP_9_6 ) then
         write(ifile,'(a)')
     &        'u(r) = epsilon*[2*(Ro/r)^9 - 3*(Ro/r)^6] - S'
         write(ifile,'(a4,1x,a7,3x,1x,a4,1x,a7,3x,3(1x,a7))') 
     &        'Num.','Atom(i)','Num.','Atom(j)','Ro','epsilon','S'
         do itype = 1, NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  write(ifile,'(i4,1x,A10,1x,i4,1x,A10,3(1x,f7.3))') 
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,twh_get_nbcoeff(1,ij)**(1.0d0/3.0d0)
     &                 ,twh_get_nbcoeff(2,ij),twh_get_nbcoeff(0,ij)
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_12_9_6 ) then
         write(ifile,'(a,1x,a)')
     &        'u(r) = nbcoeff(1)/r^12 + nbcoeff(2)/r^9'
     &        ,'+ nbcoeff(3)/r^6 - shift'
         write(ifile,'(a4,1x,a7,3x,1x,a4,1x,a7,3x,4(1x,a11))') 
     &        'Num.','Atom(i)','Num.','Atom(j)'
     &        ,'nbcoeff(1)','nbcoeff(2)','nbcoeff(3)','shift'
         do itype = 1,NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  write(ifile,'(i4,1x,a10,1x,i4,1x,a10,4(1x,e11.5))') 
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,twh_get_nbcoeff(1,ij),twh_get_nbcoeff(2,ij)
     &                 ,twh_get_nbcoeff(3,ij),twh_get_nbcoeff(0,ij)
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_BUFFERED_14_7 ) then
         write(ifile,'(a,a)') 
     &        'u(r) = epsilon*[1.07*Rstar/(r + 0.07*Rstar)]^7'
     &        ,' * [(1.12*Rstar^7 / (r^7 + 0.12*Rstar^7)) -2 ]'
         write(ifile,'(a4,1x,a7,4x,a4,1x,a7,3x,2(1x,a10))') 
     &        'Num.','Atom(i)','Num.','Atom(j)','Rstar','epsilon'
         do itype = 1,NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  write(ifile,'(i4,1x,a10,1x,i4,1x,a10,2(1x,f10.4))') 
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,twh_get_nbcoeff(5,ij)
     &                 ,twh_get_nbcoeff(6,ij)
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_DOUBLE_EXPONENTIAL ) then
         write(ifile,'(a)') 'u(r) = a exp(-b r) - c exp(-d r) - shift'
         write(ifile,'(a4,1x,a7,4x,a4,1x,a7,3x,5(3x,a5,5x))') 
     &        'Num.','Atom(i)','Num.','Atom(j)','a','b','c','d','shift'
         do itype = 1,NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  write(ifile,'(i4,1x,a10,1x,i4,1x,a10,5(1x,g12.5))') 
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,(twh_get_nbcoeff(icoeff,ij),icoeff=1,4)
     &                 ,twh_get_nbcoeff(0,ij)
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_EAM_PAIR_ONLY ) then
         write(ifile,'(a)') 
     &        'U = eam_pair(r)'
      elseif ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD ) then
         write(ifile,'(a)') 
     &        'U = eam_pair(r) + eam_embed[ Sum[eam_dens(r)] ]'
         write(ifile,'(a4,1x,a7,4x,a4,1x,a7,4x,a17)') 
     &        'Num.','Atom(i)','Num.','Atom(j)','eam_dens(r) style'
         do itype = 1,NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  call twh_eam_dens_style(GLB_GET,itype,jtype,eamdstyle)
                  write(ifile,'(i4,1x,a10,1x,i4,1x,a10,1x,a20)') 
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,eamdstyle
               endif
            enddo
         enddo
         write(ifile,'(a4,1x,a7,4x,a4,1x,a7,4x,a17)') 
     &        'Num.','Atom(i)','Num.','Atom(j)','eam_pair(r) style'
         do itype = 1,NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  call twh_eam_pair_style(GLB_GET,ij,eampstyle)
                  write(ifile,'(i4,1x,a10,1x,i4,1x,a10,1x,a20)') 
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,eampstyle
               endif
            enddo
         enddo
         write(ifile,'(a4,1x,a7,20x,a)') 
     &        'Num.','Atom(i)','eam_embed(rho) style'
         do itype = 1,NNTYPE
            ltesti = twh_get_lhere(itype)
            if ( ltesti ) then
               call twh_eam_embed_style(GLB_GET,itype,eamestyle)
               write(ifile,'(i4,1x,a10,17x,a20)') 
     &              itype,twh_get_nbname(itype),eamestyle
            endif
         enddo
      elseif ( cpnumber .eq. CP_EXPONENTIAL_12_6 ) then
         write(ifile,'(a)') 'u(r) = A/r^6 + B exp(C r) + D/r^12 - S '
         write(ifile,'(a4,1x,a7,3x,1x,a4,1x,a7,3x,2(1x,a12)
     &        ,1x,a8,1x,a12,1x,a8)') 
     &        'Num.','Atom(i)','Num.','Atom(j)','A','B','C','D','S'
         do itype = 1,NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  write(ifile,'(i4,1x,A10,1x,i4,1x,A10,2(1x,f12.1)
     &                 ,1x,f8.5,1x,f12.1,1x,f8.5)')
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,twh_get_nbcoeff(1,ij),twh_get_nbcoeff(3,ij)
     &                 ,twh_get_nbcoeff(4,ij)
     &                 ,twh_get_nbcoeff(2,ij),twh_get_nbcoeff(0,ij)
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_EXPONENTIAL_6 ) then
         write(ifile,'(a)') 'u(r) = infinity for r < expmin'
         write(ifile,'(a)')
     &        'u(r) = A/r^6 + B exp(C r) - S for r > expmin'
         write(ifile,'(a4,1x,a7,3x,1x,a4,1x,a7,4x,4(1x,a11),1x,a7)') 
     &        'Num.','Atom(i)','Num.','Atom(j)','A','B'
     &        ,'C','S','expmin'
         do itype = 1,NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  write(ifile,
     &                 '(i4,1x,a10,1x,i4,1x,a10,1x
     &                 ,4(1x,e11.4),1x,f7.3)')
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,twh_get_nbcoeff(1,ij),twh_get_nbcoeff(2,ij)
     &                 ,twh_get_nbcoeff(3,ij),twh_get_nbcoeff(0,ij)
     &                 ,dsqrt(twh_get_nbcoeff(5,ij))
c                 --- make sure rcutin is larger than expmin
                  call twh_rcutin(GLB_GET,rcutin)
                  if ( rcutin .lt. dsqrt(twh_get_nbcoeff(5,ij)) ) then
                     write(ifile,*) 'WRITENONBOND: rcutin < expmin '
                     write(ifile,*) 'this is not allowed'
                     lerror = .true.
                     return
                  endif
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_GORDON_N_6 ) then
         write(ifile,'(a,1x,a)')
     &        'u(r) = 4*epsilon*c(n)[{sigma/(r-a(n)}^n'
     &        ,'- {sigma/(r-a(n)}^n] - shift'
         write(ifile,'(a4,1x,a7,3x,1x,a4,1x,a7,3x,4(1x,a7))') 
     &        'Num.','Atom(i)','Num.','Atom(j)','sigma','epsilon'
     &        ,'n','shift'
         do itype = 1,NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  write(ifile,'(i4,1x,a10,1x,i4,1x,a10,4(1x,f7.3))') 
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,twh_get_nbcoeff(1,ij),twh_get_nbcoeff(2,ij)
     &                 ,twh_get_nbcoeff(3,ij),twh_get_nbcoeff(0,ij)
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_HARD_2580_MULTISTEP .or.
     &        cpnumber .eq. CP_REPULSIVE_2580_MULTISTEP ) then
         write(ifile,'(a39,1x,a26)') 
     &        'Energy computed for the first well that'
     &        ,'satisfies these conditions'
         if ( cpnumber .eq. CP_HARD_2580_MULTISTEP ) then
            write(ifile,'(a29)') 'u(r) = infinity for r < sigma'
         else
            write(ifile,'(a41)') 
     &           'u(r) = 1d10*(sigma^2 - r^2) for r < sigma'
         endif            
         write(ifile,'(a35)') 'u(r) = nbcoeff(2) for r < 1.2 sigma'
         write(ifile,'(a35)') 'u(r) = nbcoeff(3) for r < 1.5 sigma'
         write(ifile,'(a35)') 'u(r) = nbcoeff(4) for r < 1.8 sigma'
         write(ifile,'(a35)') 'u(r) = nbcoeff(5) for r < 2.0 sigma'

         write(ifile,'(a4,1x,a7,3x,1x,a4,1x,a7,3x,1x,a10,4(1x,a15))') 
     &        'Num.','Atom(i)','Num.','Atom(j)','sigma'
     &        ,'nbcoeff(2)','nbcoeff(3)','nbcoeff(4)','nbcoeff(5)'
         do itype = 1, NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
c                 --- output the hard core
                  write(ifile,'(2(i4,1x,a10,1x),f10.5,4(1x,f15.5))')
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,twh_get_nbcoeff(1,ij),twh_get_nbcoeff(2,ij)
     &                 ,twh_get_nbcoeff(3,ij)
     &                 ,twh_get_nbcoeff(4,ij),twh_get_nbcoeff(5,ij)
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_HARD_SPHERE ) then
         write(ifile,'(a)') 'u(r) = infinity for r < sigma'
         write(ifile,'(a)') 'u(r) = 0 for r > sigma'
         write(ifile,'(a4,1x,a7,3x,1x,a4,1x,a7,3x,1x,a7)') 
     &        'Num.','Atom(i)','Num.','Atom(j)','sigma'
         do itype = 1, NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  write(ifile,'(i4,1x,a10,1x,i4,1x,a10,1x,f7.3)') 
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,twh_get_nbcoeff(1,ij)
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_LENNARD_JONES .or. 
     &     cpnumber .eq. CP_12_6_PLUS_SOLVATION ) then
         write(ifile,'(a,a)') 'u(r) = 4*epsilon[(sigma/r)^12'
     &        ,' - (sigma/r)^6] - shift'
         write(ifile,'(a4,1x,a7,4x,a4,1x,a7,3x,5(1x,a10))') 
     &        'Num.','Atom(i)','Num.','Atom(j)','sigma','epsilon'
     &        ,'shift','1-4sig','1-4eps'
         do itype = 1,NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  write(ifile,'(i4,1x,a10,1x,i4,1x,a10,5(1x,f10.4))') 
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,(twh_get_nbcoeff(13,ij))**(1.0/6.0)
     &                 , twh_get_nbcoeff(14,ij)/4.0d0
     &                 ,twh_limited_double(twh_get_nbcoeff(0,ij),1.0d-6)
     &                 , (twh_get_nbcoeff(15,ij))**(1.0/6.0)
     &                 ,twh_get_nbcoeff(16,ij)/4.0
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_SCALED_LENNARD_JONES ) then
         write(ifile,'(a,1x,a)')
     &        'u(r) = 4*epsilon[1/(alpha(1-lambda)^2*+(r/sigma)^6)^2'
     &        ,'- 1/(alpha(1-lambda)^2+(r/sigma)^6)] - shift'
         write(ifile,'(2(1x,a,1x,g12.6))') 'alpha ='
     &        ,scp_get_native_lambda(ALPHA_LJ),'lambda ='
     &        ,scp_get_native_lambda(LAMBDA_LJ)
         write(ifile,'(a4,1x,a7,3x,1x,a4,1x,a7,3x,5(1x,a7))') 
     &        'Num.','Atom(i)','Num.','Atom(j)','sigma','epsilon'
     &        ,'shift','1-4sig','1-4eps'
         do itype = 1,NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  write(ifile,'(i4,1x,a10,1x,i4,1x,a10,5(1x,f7.3))') 
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,(twh_get_nbcoeff(13,ij))**(1.0/6.0)
     &                 , twh_get_nbcoeff(14,ij)/4.0d0
     &                 ,twh_get_nbcoeff(0,ij)
     &                 , (twh_get_nbcoeff(15,ij))**(1.0/6.0)
     &                 ,twh_get_nbcoeff(16,ij)/4.0
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_12_6_PLUS_12_10_H_BOND ) then
         write(ifile,'(a,1x,a)')
     &        'u(r) = 4*epsilon[(sigma/r)^12 - (sigma/r)^6]'
     &        ,'+ C/(r^12) - D/(r^10) - shift'
         write(ifile,'(a4,1x,a7,3x,1x,a4,1x,a7,3x,5(1x,a7),4(1x,a9))') 
     &        'Num.','Atom(i)','Num.','Atom(j)','sigma','epsilon'
     &        ,'shift','1-4 sig','1-4 eps','C','D','1-4 C','1-4 D'
         do itype = 1,NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  call twh_nb_to_hbtype(GLB_GET,ij,hbtype)
                  if ( hbtype .eq. 0 ) then
                     cvalue = 0.0d0
                     dvalue = 0.0d0
                     sfactor = 0.0d0
                  else
                     cvalue = twh_get_hbondcoeff(hbtype,1)
                     dvalue = twh_get_hbondcoeff(hbtype,2)
                     sfactor = twh_get_hbondcoeff(hbtype,3)
                  endif
                  write(ifile,'(i4,1x,a10,1x,i4,1x,a10,5(1x,f7.3)
     &                 ,4(1x,e9.3))') 
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,(twh_get_nbcoeff(13,ij))**(1.0/6.0)
     &                 , twh_get_nbcoeff(14,ij)/4.0d0
     &                 ,twh_get_nbcoeff(0,ij)
     &                 , (twh_get_nbcoeff(15,ij))**(1.0/6.0)
     &                 ,twh_get_nbcoeff(16,ij)/4.0,cvalue,dvalue
     &                 ,sfactor*cvalue,sfactor*dvalue
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_MULTIWELL ) then
         write(ifile,'(a)') 'u(r) = infinity for r < table_pair(1,1)'
         write(ifile,'(a)') 'u(r) = table_pair(2,n)'
         write(ifile,'(3x,a)') 
     &        'for table_pair(1,n-1) < r < table_pair(1,n)'
         write(ifile,'(a4,1x,a7,3x,1x,a4,1x,a7,3x,1x,a5,2(1x,a15))') 
     &        'Num.','Atom(i)','Num.','Atom(j)','npair'
     &        ,'table_pair(1,n)','table_pair(2,n)'
         do itype = 1, NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
c                 --- output the hard core
                  icoeff = 1
                  call twh_table_pair(GLB_GET,1,ij,icoeff,tpairone)
                  write(ifile,'(i4,1x,a10,1x,i4,1x,a10,1x,i5
     &                 ,1x,f15.3,1x,7x,a8)')
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,icoeff,tpairone,'infinity'
c                 --- output all of the wells
                  call twh_table_npair(GLB_GET,ij,tabnpair)
                  do icoeff = 2,tabnpair
                     call twh_table_pair(GLB_GET,1,ij,icoeff,tpairone)
                     call twh_table_pair(GLB_GET,2,ij,icoeff,tpairtwo)
                     write(ifile,'(32x,i5,2(1x,f15.3))')
     &                    icoeff,tpairone,tpairtwo
                  enddo
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_REPULSIVE_MULTIWELL ) then
         write(ifile,'(a)') 'u(r) = 1d5 * [table_pair(2,1)-r]'
     &        ,' for r < table_pair(1,1)'
         write(ifile,'(a)') 'u(r) = table_pair(2,n)'
         write(ifile,'(3x,a)') 
     &        'for table_pair(1,n-1) < r < table_pair(1,n)'
         write(ifile,'(a4,1x,a7,3x,1x,a4,1x,a7,3x,1x,a5,2(1x,a15))') 
     &        'Num.','Atom(i)','Num.','Atom(j)','npair'
     &        ,'table_pair(1,n)','table_pair(2,n)'
         do itype = 1, NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
c                 --- output the hard core
                  icoeff = 1
                  call twh_table_pair(GLB_GET,1,ij,icoeff,tpairone)
                  write(ifile,'(i4,1x,a10,1x,i4,1x,a10,1x,i5
     &                 ,1x,f15.3,1x,7x,a9)')
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,icoeff,tpairone,'repulsion'
c                 --- output all of the wells
                  call twh_table_npair(GLB_GET,ij,tabnpair)
                  do icoeff = 2,tabnpair
                     call twh_table_pair(GLB_GET,1,ij,icoeff,tpairone)
                     call twh_table_pair(GLB_GET,2,ij,icoeff,tpairtwo)
                     write(ifile,'(32x,i5,2(1x,f15.3))')
     &                    icoeff,tpairone,tpairtwo
                  enddo
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_REPULSIVE_SPHERE ) then
         write(ifile,'(a)')
     &        'u(r) = 1d5 + 1d5*(sigma^2 - r^2) for r < sigma'
         write(ifile,'(a)') 'u(r) = 0 for r > sigma'
         write(ifile,'(a4,1x,a7,3x,1x,a4,1x,a7,3x,1x,a7)') 
     &        'Num.','Atom(i)','Num.','Atom(j)','sigma'
         do itype = 1, NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  write(ifile,'(i4,1x,a10,1x,i4,1x,a10,1x,f7.3)') 
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,twh_get_nbcoeff(1,ij)
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_REPULSIVE_WELL ) then
         write(ifile,'(a)') 
     &        'u(r) = 1d5 + 1d5*(nb(1)^2 - r^2) for r < nbcoeff(1)'
         write(ifile,'(a)') 
     &        'u(r) = -nbcoeff(3) for nbcoeff(1) < r < nbcoeff(2)'
         write(ifile,'(a4,1x,a7,3x,1x,a4,1x,a7,3x,3(1x,a10))') 
     &        'Num.','Atom(i)','Num.','Atom(j)','nbcoeff(1)'
     &        ,'nbcoeff(2)','nbcoeff(3)'
         do itype = 1, NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  write(ifile,'(i4,1x,a10,1x,i4,1x,a10,3(1x,f10.3))')
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,dsqrt(twh_get_nbcoeff(1,ij))
     &                 ,dsqrt(twh_get_nbcoeff(2,ij))
     &                 ,twh_get_nbcoeff(3,ij)
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_SQUARE_WELL ) then
         write(ifile,'(a)') 'u(r) = infinity for r < nbcoeff(1)'
         write(ifile,'(a)') 
     &        'u(r) = -nbcoeff(3) for nbcoeff(1) < r < nbcoeff(2)'
         write(ifile,'(a4,1x,a7,3x,1x,a4,1x,a7,3x,3(1x,a10))') 
     &        'Num.','Atom(i)','Num.','Atom(j)','nbcoeff(1)'
     &        ,'nbcoeff(2)','nbcoeff(3)'
         do itype = 1, NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  write(ifile,'(i4,1x,a10,1x,i4,1x,a10,3(1x,f10.3))')
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,dsqrt(twh_get_nbcoeff(1,ij))
     &                 ,dsqrt(twh_get_nbcoeff(2,ij))
     &                 ,twh_get_nbcoeff(3,ij)
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_STILLINGER_WEBER .or.
     &        cpnumber .eq. CP_SW_PAIR_ONLY ) then
         if ( cpnumber .eq. CP_STILLINGER_WEBER ) then
            write(ifile,'(a,1x,a)') 
     &           'U = nbcoeff(1)*[nbcoeff(2)*Sum u2(rij)'
     &           ,'+ nbcoeff(7)*Sum u3(rij,rjk)]'
         else
c           --- just the twobody part
            write(ifile,'(a)') 'U = nbcoeff(1)*[nbcoeff(2)*Sum u2(rij)]'
         endif
         write(ifile,'(a,a,1x,a,1x,a)')
     &        'u2(r_ij) = [nbcoeff(3)*(r_ij/nbcoeff(4)^'
     &        ,'-nbcoeff(5) - 1]'
     &        ,'* Exp{1/[(r_ij/nbcoeff(4) - nbcoeff(6))]}'
     &        ,'* Heaviside(nbcoeff(6) - [r_ij/nbcoeff(4)])'
         if ( cpnumber .eq. CP_STILLINGER_WEBER ) then
            write(ifile,'(a,1x,a,1x,a,1x,a,/,3x,a,1x,a)') 
     &           'u3(r_ij,r_jk) = exp[ nbcoeff(8)/{r_ij/nbcoeff(4)'
     &           ,'- nbcoeff(6)}'
     &           ,'+ nbcoeff(8)/{r_jk/nbcoeff(4) - nbcoeff(6)}]'
     &           ,'* (cos(theta_ijk)-nbcoeff(9))^2'
     &           ,'* Heaviside(nbcoeff(6) - r_ij/nbcoeff(4))'
     &           ,'* Heaviside(nbcoeff(6) - r_jk/nbcoeff(4))'
         endif
         write(ifile,*)

         if ( cpnumber .eq. CP_STILLINGER_WEBER ) then
            write(ifile,'(a,1x,a)')
     &           '-i- -j- nbcoeff(1) nb(2) nb(3) nb(4) nb(5) nb(6)'
     &           ,'nb(7) nb(8) nb(9)'
            iswmax = 9
         else
            write(ifile,'(a)')
     &           '-i- -j- nbcoeff(1) nb(2) nb(3) nb(4) nb(5) nb(6)'
            iswmax = 6
         endif
         do itype = 1,NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  write(ifile,130) itype,jtype
     &                 ,(twh_get_nbcoeff(itemp,ij),itemp=1,iswmax)
               endif
            enddo
         enddo
 130     format(i3,1x,i3,1x,e10.4,8(1x,f5.2)) 
      elseif ( cpnumber .eq. CP_TABULATED_PAIR ) then
         write(ifile,'(a)') 'U = eam_pair(r) '

      elseif ( cpnumber .eq. CP_UFF_12_6 ) then
         write(ifile,'(a,1x,a)')
     &        'u(r) = 4*epsilon[(sigma/r)^12 - (sigma/r)^6]'
     &        ,'- shift'
         write(ifile,'(a4,1x,a7,3x,1x,a4,1x,a7,3x,6(1x,a7))') 
     &        'Num.','Atom(i)','Num.','Atom(j)','sigma','epsilon'
     &        ,'shift','1-4sig','1-4eps','Z'
         do itype = 1,NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  write(ifile,'(i4,1x,a10,1x,i4,1x,a10,6(1x,f7.3))') 
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,(twh_get_nbcoeff(13,ij))**(1.0/6.0)
     &                 , twh_get_nbcoeff(14,ij)/4.0d0
     &                 ,twh_get_nbcoeff(0,ij)
     &                 , (twh_get_nbcoeff(15,ij))**(1.0/6.0)
     &                 ,twh_get_nbcoeff(16,ij)
     &                 /4.0,twh_get_nbcoeff(10,ij)
               endif
            enddo
         enddo
         write(ifile,'(a4,7x,8(a8,1x))') 
     &        'Atom','r','theta','x','D','zeta','Z','chi','V'
         do itype = 1,NNTYPE
            ltesti = twh_get_lhere(itype)
            if ( ltesti ) then
               ij = twh_getnbtype(itype,itype)
               write(ifile,'(a10,1x,8(f8.4,1x))') twh_get_nbname(itype)
     &              ,(twh_get_nbcoeff(itemp,ij),itemp=5,12)
            endif
         enddo
      elseif ( cpnumber .eq. CP_WEEKS_CHANDLER_ANDERSON ) then
         write(ifile,'(a,1x,a)')
     &        'u(r) = 4*epsilon[(sigma/r)^12 - (sigma/r)^6]'
     &        ,'+ epsilon - S {for rij > 2**(1/6)sigma, else 0}'
         write(ifile,'(a4,1x,a7,4x,a4,1x,a7,3x,5(1x,a10))') 
     &        'Num.','Atom(i)','Num.','Atom(j)','sigma','epsilon'
     &        ,'shift','1-4sig','1-4eps'
         do itype = 1,NNTYPE
            ltesti = twh_get_lhere(itype)
            do jtype = itype,NNTYPE
               ltestj = twh_get_lhere(jtype)
               if ( ltesti .and. ltestj ) then
                  ij = twh_getnbtype(itype,jtype)
                  write(ifile,'(i4,1x,a10,1x,i4,1x,a10,5(1x,f10.4))') 
     &                 itype,twh_get_nbname(itype)
     &                 ,jtype,twh_get_nbname(jtype)
     &                 ,(twh_get_nbcoeff(13,ij))**(1.0/6.0)
     &                 , twh_get_nbcoeff(14,ij)/4.0d0
     &                 ,twh_get_nbcoeff(0,ij)
     &                 , (twh_get_nbcoeff(15,ij))**(1.0/6.0)
     &                 ,twh_get_nbcoeff(16,ij)/4.0
               endif
            enddo
         enddo
      else
         write(ifile,*) 'READTOWHEE: unknown classical potential of '
     &        ,cpnumber
         lerror = .true.
         return
      endif

      return
      end
    
