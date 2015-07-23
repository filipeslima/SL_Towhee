#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2003-2013 Marcus G. Martin                       *
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
      subroutine twh_stresstensor(lfailure,stress,ibox)
c     ******************************************************************
c     * computes the stress tensor using the molecular virial          *
c     * From Allen and Tildesley equation 2.54                         *
c     * PV = NkT + <W>                                                 *
c     * equation 2.61                                                  *
c     * W = -1/3 Sum{i}Sum{j>i} w(r_ij)                                *
c     * and equation 2.63                                              *
c     * w(r_ij) = Sum{a}Sum{b} w_ab(r_ij)/r_ab^2 (r_ab*r_ij)           *
c     * where a and b are the individual atoms in each molecule and    *
c     * i and j are the centers of mass of the molecules               *
c     *                                                                *
c     * The stress tensor is symetric so only the upper half is stored *
c     * Total of 6 elements, plus a 0 element that is the trace        *
c     * 0: trace of the diagonal elements (ie the isotropic pressure)  *
c     * 1: x-x stress                                                  *
c     * 2: y-y stress                                                  *
c     * 3: z-z stress                                                  *
c     * 4: x-y stress                                                  *
c     * 5: x-z stress                                                  *
c     * 6: y-z stress                                                  *
c     *                                                                *
c     * rewritten from pressure 12-17-2003 by M.G. Martin              *
c     * last modified 07-26-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_GT
#define FUNCTION_GET_NBOXI
#define FUNCTION_WCOULOMB
#define FUNCTION_WTWOBODY
#include "functions.h"

#define DEBUG_STRESSTENSOR 0
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ibox
      double precision stress
      dimension stress(0:7)

c     --- local variables
      logical lovrlap,leval_c,leval_v,lperiodic,scp_scale_twobody
      logical luse_ci,luse_cj,luse_vi,luse_vj,ltailc,lcoulomb
      logical lewald
c     --- integer scalars
      integer ichain,jchain,indx,ensemble
      integer ii,jj,ntii,ntjj,imolty,jmolty,nchain,nmolty
      integer numiunit,numjunit,nchbox,ncmti,ncmtj
c     --- double precision scalars
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,rx_ab,ry_ab,rz_ab
      double precision rijsq,rhosq,rij,w_ab,xcmi,ycmi,zcmi
      double precision rcm,rcmsq,rcmi,volsq,temperature,corp
      double precision rx_ij,ry_ij,rz_ij,xcmj,ycmj,zcmj,rcmj,rcelect
      double precision boxvolume
c     --- double precision arrays
      double precision repress
      dimension repress(6)

c     --- retrieve constants
      call twh_rcelect(GLB_GET,ibox,rcelect)
      call twh_lcoulomb(GLB_GET,lcoulomb)
      call twh_lewald(GLB_GET,lewald)
      call twh_boxvolume(GLB_GET,ibox,boxvolume)

c     --- initialize stress tensor
      do indx = 0,6
         stress(indx) = 0.0d0
      enddo
c     --- stress tensor meaningless in box 0 if grand canonical
      call twh_ensemble(GLB_GET,ensemble)
      if ( ensemble .eq. ENS_muVT .and. ibox .eq. 0 ) return

c     --- retrieve constants
      call twh_temperature(GLB_GET,temperature)
      call twh_scp_scale_twobody(GLB_GET,scp_scale_twobody)

c     --- get the total number of molecules
      call twh_nchain(GLB_GET,nchain)

c     --- make sure there are no periodic molecules in this box as the
c     --- molecular virial does not make sense for molecules that do
c     --- not have a center-of-mass
      do ichain = 1,nchain
         if ( twh_get_nboxi(ichain) .eq. ibox ) then
            call twh_lperiodic(GLB_GET,ichain,lperiodic)
            if ( lperiodic ) then
c              --- cannot compute the virial pressure in this box
c              --- so leave it as zero and return
               return
            endif
         endif
      enddo

c     --- loop over all molecules
      do ichain = 1, nchain - 1
c        --- check if i is in relevant box 
         if ( twh_get_nboxi(ichain) .eq. ibox ) then
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,numiunit)
            call twh_com(GLB_GET,CRD_REAL,ichain,xcmi,ycmi,zcmi)
            call twh_rcmu(GLB_GET,ichain,rcmi)
c           --- loop over all chains j with j>i 
            do jchain = ichain + 1, nchain
c              --- check for simulation box
               if ( twh_get_nboxi(jchain) .eq. ibox ) then
                  call twh_moltyp(GLB_GET,jchain,jmolty)
                  call twh_nunit(GLB_GET,jmolty,numjunit)
c                 --- need to compute the center of mass separations
                  call twh_com(GLB_GET,CRD_REAL,jchain,xcmj,ycmj,zcmj)
                  call twh_mimage (ibox,xcmi,ycmi,zcmi,xcmj,ycmj,zcmj
     &                 ,rx_ij,ry_ij,rz_ij )
                  rijsq = rx_ij*rx_ij + ry_ij*ry_ij + rz_ij*rz_ij
c                 --- set rcm to the nonbond cutoff
                  call twh_rcut(GLB_GET,rcm)
                  if ( lcoulomb ) then
c                    --- compare with the electrostatic cutoff
                     if ( rcelect .gt. rcm ) rcm = rcelect
                  endif
c                 --- add in the maximum distances from the COM
                  call twh_rcmu(GLB_GET,jchain,rcmj)
                  rcm = rcm + rcmi + rcmj
c                 --- square the distance
                  rcmsq = rcm*rcm
c                  if ( rijsq .le. rcmsq ) then
                  if (twh_cmp_gt(rcmsq, rijsq)) then
c                    --- loop over all beads of ichain  
                     do ii = 1, numiunit
                        call twh_ntype(GLB_GET,imolty,ii,ntii)
                        call twh_coordinates(GLB_GET,CRD_REAL,ichain,ii
     &                       ,rxui,ryui,rzui)
                        call twh_luse_c(GLB_GET,imolty,ii,luse_ci)
                        call twh_luse_v(GLB_GET,imolty,ii,luse_vi)
c                       --- loop over all beads of jchain  
                        do jj = 1, numjunit
                           call twh_ntype(GLB_GET,jmolty,jj,ntjj)
c                          --- check if nonzero interaction sites
                           call twh_luse_v(GLB_GET,jmolty,jj,luse_vj)
                           leval_v = luse_vi .and. luse_vj
                           if ( lcoulomb ) then
                              call twh_luse_c(GLB_GET,jmolty,jj,luse_cj)
                              leval_c = luse_ci .and. luse_cj
                           else
                              leval_c = .false.
                           endif
                           call twh_coordinates(GLB_GET,CRD_REAL
     &                          ,jchain,jj,rxuj,ryuj,rzuj)
c                          --- minimum image the pair separations
                           call twh_mimage ( ibox,rxui,ryui,rzui
     &                          ,rxuj,ryuj,rzuj,rx_ab,ry_ab,rz_ab )
                           rijsq = rx_ab*rx_ab + ry_ab*ry_ab 
     &                          + rz_ab*rz_ab
c                          --- initialize the pair virial function
                           w_ab = 0.0d0
                           if ( leval_c ) then
c                             --- compute virial from coulombic terms
                              rij = dsqrt(rijsq)
                              w_ab = w_ab + twh_wcoulomb(ibox,rij,
     &                             imolty,ii,jmolty,jj)
                           endif
                           if ( leval_v ) then
c                             --- compute pair virial from vdw terms
                              w_ab = w_ab + twh_wtwobody(lfailure,rijsq
     &                             ,ntii,ntjj,lovrlap)
                              if ( lfailure ) return
                              if ( lovrlap ) then
                                 write(6,*) 'STRESSTENSOR: overlap'
     &                                ,' in twobody'
c                                --- overlap is OK if we're scaling
                                 if (.not. scp_scale_twobody) then
                                    write(6,*) 'rijsq,ntii,ntjj:'
     &                                   ,rijsq,ntii,ntjj
                                    write(6,*) 'imolty,ii:',imolty,ii
                                    write(6,*) 'jmolty,jj:',jmolty,jj
                                    write(6,*) 'ichain,jchain:'
     &                                   ,ichain,jchain
                                    write(6,*) 'luse_v(imolty,ii):'
     &                                   ,luse_vi
                                    write(6,*) 'luse_v(jmolty,jj):'
     &                                   ,luse_vj
                                    write(6,*) 'luse_c(imolty,ii):'
     &                                   ,luse_ci
                                    write(6,*) 'luse_c(jmolty,jj):'
     &                                   ,luse_cj
                                    lfailure = .true.
                                    return
                                 endif
                              endif
                           endif
                           stress(1) = stress(1) + w_ab * rx_ab * rx_ij
                           stress(2) = stress(2) + w_ab * ry_ab * ry_ij
                           stress(3) = stress(3) + w_ab * rz_ab * rz_ij
c                          --- use the average of the two possible 
c                          --- combinations of the different indices
                           stress(4) = stress(4) + w_ab * 
     &                          0.5d0*(rx_ab * ry_ij + ry_ab * rx_ij)
                           stress(5) = stress(5) + w_ab * 
     &                          0.5d0*(rx_ab * rz_ij + rz_ab * rx_ij)
                           stress(6) = stress(6) + w_ab * 
     &                          0.5d0*(ry_ab * rz_ij + rz_ab * ry_ij)
                        enddo
                     enddo
                  endif
               endif
            enddo
         endif
      enddo

      if ( lewald ) then
c        --- Compute the reciprocal space contribution
c        --- by using the thermodynamic definition
         call twh_recippress(ibox,repress)
         do indx = 1,6
            stress(indx) = stress(indx) - repress(indx)
         enddo
      endif

c     --- get isotropic pressure from the trace
      stress(0) = (stress(1) + stress(2) + stress(3))/3.0d0

#if DEBUG_STRESSTENSOR
      c = CON_PRESSCONST / boxvolume
      write(6,*) "STRESSTENSOR: pressure components"
      call twh_nchbox(GLB_GET,ibox,nchbox)
      write(6,*) "Ideal:      ", c*nchbox*temperature
      write(6,*) "isotropic:  ", -c*stress(0)
      write(6,*) "x-x stress: ", -c*stress(1)
      write(6,*) "y-y stress: ", -c*stress(2)
      write(6,*) "z-z stress: ", -c*stress(3)
      write(6,*) "x-y stress: ", -c*stress(4)
      write(6,*) "x-z stress: ", -c*stress(5)
      write(6,*) "y-z stress: ", -c*stress(6)
#endif
c
c     --- turn <W> into pressure
      call twh_nchbox(GLB_GET,ibox,nchbox)
      do indx = 0,6
         stress(indx) = CON_PRESSCONST*(( dble(nchbox)*temperature) 
     &        - ( stress(indx) )) / boxvolume
      enddo

c     --- isotropic pressure tail correction for pair potential
      call twh_nmolty(GLB_GET,nmolty)
      stress(7) = 0.0d0
      call twh_ltailc(GLB_GET,ltailc)
      if ( ltailc ) then
c        -- add tail corrections for the VDW energy
         volsq = ( boxvolume )**2
         do imolty=1, nmolty
            call twh_ncmt(GLB_GET,ibox,imolty,ncmti)
            do jmolty=1, nmolty
               call twh_ncmt(GLB_GET,ibox,jmolty,ncmtj)
               rhosq = dble(ncmti)*dble(ncmtj) / volsq
               call twh_corp(GLB_GET,imolty,jmolty,corp)
               stress(7) = stress(7) + rhosq*corp
            enddo
         enddo
         stress(7) = CON_PRESSCONST*stress(7)
      endif

c     --- only adding into the isotropic pressure right now as not 
c     --- sure how this would work otherwise
      stress(0) = stress(0) + stress(7)

c      write(6,*) 'end stresstensor'

      return
      end
