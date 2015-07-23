#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2013 Marcus G. Martin                       *
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
      subroutine twh_bondorder(lfailure,ntii,ntjj,rij,bij,bji,lzero)
c     ******************************************************************
c     * computes the bondorder between two atoms in a multi-body force *
c     * field                                                          *
c     *                                                                *
c     * originally written 06-26-2002 by  M.G. Martin                  *
c     * last update 07-21-2013 by M.G. Martin                          *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_NBCOEFF
#define FUNCTION_EXTRACTDENS
#define FUNCTION_EXPON
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure,lzero
      integer ntii,ntjj
      double precision rij,bij,bji
c     --- local variables
c     --- logical scalars
      logical lcontinue
c     --- integer scalars
      integer ntij,looptype_a,looptype_b,loop,icoeff,potnum,eamdensnum
c     --- double precision scalars
      double precision bloop,x,eamdensone,eamdenstwo,rcut,maxdist

c     --- retrieve constants
      call twh_rcut(GLB_GET,rcut)

c     --- initialize lzero
      lzero = .false.
      call twh_cp_number(GLB_GET,potnum)
      if ( potnum .eq. CP_STILLINGER_WEBER ) then
         ntij = twh_getnbtype(ntii,ntjj)
c        --- stillinger-weber uses a heaviside to compute bond
c        --- order so it is either zero or one
         if ( rij/twh_get_nbcoeff(4,ntij) 
     &        .lt. twh_get_nbcoeff(6,ntij)) then
c           --- positive value in a heaviside results in 1.0
            bij = 1.0d0
            lzero = .false.
         else
c           --- negative value in a heaviside results in 0.0
            bij = 0.0d0
            lzero = .true.
         endif
      elseif ( potnum .eq. CP_EMBEDDED_ATOM_METHOD ) then
c        --- get the two densities
         if ( rij .ge. rcut ) then
            bij = 0.0d0
            bji = 0.0d0
         else
            do loop = 1,2
               if ( loop .eq. 1 ) then
                  looptype_a = ntii
                  looptype_b = ntjj
               else
                  looptype_a = ntjj
                  looptype_b = ntii
               endif
               call twh_eam_dens_num(GLB_GET,looptype_a,looptype_b
     &              ,eamdensnum)
               if ( eamdensnum .eq. EAM_DENS_TABLE ) then
                  bloop = twh_extractdens(lfailure,INTERP_NORMAL
     &                 ,looptype_a,looptype_b
     &                 ,rij)
                  if ( lfailure ) return
               elseif ( eamdensnum .eq. EAM_DENS_EXPON ) then
c                 --- check maximum cutoff
                  call twh_eam_dens(GLB_GET,1,looptype_a,looptype_b,2
     &                 ,maxdist)
                  if ( rij .le. maxdist ) then
c                    --- compute density
                     call twh_eam_dens(GLB_GET,2,looptype_a,looptype_b,1
     &                    ,eamdensone)
                     call twh_eam_dens(GLB_GET,2,looptype_a,looptype_b,2
     &                    ,eamdenstwo)
                     bloop = eamdensone * twh_expon(eamdenstwo*rij)
                  else
                     bloop = 0.0d0
                  endif
               elseif ( eamdensnum .eq. EAM_DENS_ACKLAND ) then
c                 --- cubic sum with heaviside (using a while loop to 
c                 --- recreate the heaviside and save some expense)
                  bloop = 0.0d0
                  call twh_eam_ndens(GLB_GET,looptype_a,looptype_b
     &                 ,icoeff)
                  if ( icoeff .gt. 0 ) then
                     lcontinue = .true.
                  else
                     lcontinue = .false.
                  endif
                  do while ( lcontinue )
                     call twh_eam_dens(GLB_GET,1,looptype_a,looptype_b
     &                    ,icoeff,eamdensone)
                     x = eamdensone - rij
                     if ( x .gt. 0.0d0 ) then
c                       --- compute density contribution
                        call twh_eam_dens(GLB_GET,2,looptype_a
     &                       ,looptype_b,icoeff,eamdenstwo)
                        bloop = bloop + eamdenstwo*x*x*x
c                       --- decrement for the next round
                        icoeff = icoeff - 1
                        if ( icoeff .eq. 0 ) lcontinue = .false.
                     else
c                       --- exceeded Heavyside functions, this loop
c                       --- is finished
                        lcontinue = .false.
                     endif
                  enddo
               elseif ( eamdensnum .eq. EAM_DENS_NONE ) then
c                  --- no contribution to the density
                  bloop = 0.0d0
               else
                  write(6,*) 'BONDORDER: unknown eam_dens_num'
                  write(6,*) 'looptype(a,b)',looptype_a,looptype_b
     &                 ,'eam_dens_num:',eamdensnum
                  lfailure = .true.
                  return
               endif
               
               if ( loop .eq. 1 ) then
                  bij = bloop
               else 
                  bji = bloop
               endif
            enddo
         endif
      else
         write(6,*) 'BONDORDER: unknown classical_potential ',potnum
         lfailure = .true.
         return
      endif

      return
      end
