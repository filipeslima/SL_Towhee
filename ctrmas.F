#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2013 Marcus G. Martin                       *
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
      subroutine twh_ctrmas (lfailure, iall, ibox, ichain, mtype)
c     ******************************************************************
c     * finds the center of mass of a chain and returns it to the      *
c     * periodic box if it has left the box.                           *
c     * iall 0 means compute ctrmas for all molecules in this box      *
c     *   CRD_REAL compute ctrmas for one molecule using real coords   *
c     *   CRD_TEMP compute ctrmas for one molecule using temp coords   *
c     *   CRD_SINGLE compute ctrmas for one molecule using single crds *
c     * ibox is the box of the chain                                   *
c     * ichain is the chain number.                                    *
c     * originally written in 1996 by M.G. Martin                      *
c     * last modified 07-15-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_GET_MASS
#define FUNCTION_GET_NBOXI
#define FUNCTION_MAXBOXLENGTH
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer iall,ibox,ichain,mtype

c     --- local variables
c     --- logical scalars
      logical lrescale,lwarn,lnewrcmu,ladjust,lfound,lperiodic
      logical lxx,lyy,lzz,lreset,lbondgraph
c     --- integer scalars
      integer iunit,imolty,nunit,stt,edd,flagon,itype,iflag,nchain
      integer testchain,cenunit,itogo,itest,junit,ivib,nvib,endunit
      integer testunit,testring,istorage
c     --- double precision parameters
      double precision tolerance
      parameter (tolerance = 1.0d-10)
c     --- double precision scalars
      double precision dx,dy,dz,nxcm,nycm,nzcm,xold,yold,zold
      double precision dmaxsq,rxuij,ryuij,rzuij,rijsq
      double precision rxui,ryui,rzui,rcmu
      double precision xtemp,ytemp,ztemp,gyration
      double precision boxcenx,boxceny,boxcenz,arba,arbb,arbc
      double precision shortest,altdist,dist,masstotal

c      write(6,*) 'start CTRMAS',ibox

      flagon = 0
      if ( iall .eq. 0 ) then
         stt = 1
         call twh_nchain(GLB_GET,nchain)
         edd = nchain
         iflag = CRD_REAL
      elseif ( iall .eq. CRD_REAL .or. iall .eq. CRD_TEMP 
     &        .or. iall .eq. CRD_SINGLE ) then
         stt = ichain
         edd = ichain
         iflag = iall
      else
         write(6,*) 'CTRMAS: invalid iall'
         write(6,*) 'iall',iall
         lfailure = .true.
         return
      endif
c     --- store in either real or temp
      if ( iflag .eq. CRD_REAL ) then
         istorage = CRD_REAL
      else
         istorage = CRD_TEMP
      endif

      if ( mtype .eq. CTR_INITIAL ) then 
c        --- initial readtowhee move, warn on rescale, recompute rcmu
         lwarn = .true.
         lnewrcmu = .true.
      else
c        --- move that can displace COM and changes rcmu
         lwarn = .false.
         lnewrcmu = .true.
      endif

c     --- determine whether we should alter the coordinates to put them
c     --- back into the central box.  We do not want to do this if we
c     --- are in the midst of a move, only at the end of an accepted
c     --- move.
      if ( mtype .eq. CTR_ENERGY .or. mtype .eq. CTR_INITCONF ) then
         lrescale = .false.
      else
         lrescale = .true.
      endif

      do testchain = stt, edd
         call twh_moltyp(GLB_GET,testchain,imolty)
         call twh_nunit(GLB_GET,imolty,nunit)
         call twh_masstotal(GLB_GET,imolty,masstotal)
         call twh_lbondgraph(GLB_GET,imolty,lbondgraph)
c         ----- Check if testchain is in the correct box
         if ( twh_get_nboxi(testchain) .eq. ibox) then
            if ( mtype .eq. CTR_INITIAL ) then
c              --- make sure all of the atoms are in the same reference
c              --- frame.  This could get messed up for molecules
c              --- generated outside of Towhee and therefore could
c              --- occur at the start of a simulation
c              --- find the center of the box
               arba = 0.5d0
               arbb = 0.5d0
               arbc = 0.5d0
               call twh_arbtocart(ibox,arba,arbb,arbc
     &              ,boxcenx,boxceny,boxcenz)
c              --- find the atom that is the closest to the center of 
c              --- the box, initialize shortest
               shortest = 1.0d10
               do iunit = 1,nunit
c                 --- get the coordinates for this atom
                  call twh_coordinates(GLB_GET,iflag,testchain,iunit
     &                 ,xtemp,ytemp,ztemp)
c                 --- compute the distance to the center of the box
                  call twh_mimage ( ibox,xtemp,ytemp,ztemp
     &                 ,boxcenx,boxceny,boxcenz ,rxuij,ryuij,rzuij )
                  dist = twh_distance(rxuij,ryuij,rzuij)
                  if ( iunit .eq. 1 .or. dist .lt. shortest ) then
c                    --- default to this atom to start the search
                     cenunit = iunit
                     shortest = dist
                  endif
               enddo
c              --- initialize lexist
               call twh_lexist(GLB_INIT,nunit,.false.)
c              --- set glist to all of the atoms in the molecule except
c              --- for cenunit
               itogo = 0
               do iunit = 1,nunit
                  if ( iunit .ne. cenunit ) then
                     itogo = itogo + 1
                     call twh_glist(GLB_SET,itogo,iunit)
                  endif
               enddo
c              --- fetch the center atom coordinates
               call twh_coordinates(GLB_GET,iflag,testchain,cenunit
     &              ,xtemp,ytemp,ztemp)
c              --- this atom exists in the proper reference frame
               call twh_lexist(GLB_SET,cenunit,.true.)
c              --- repeatedly loop through all of the atoms until they
c              --- are all put into the same reference frame as an atom
c              --- that exists.
               itest = 0
               lreset = .false.
               do while ( itogo .gt. 0 ) 
                  itest = itest + 1
                  call twh_glist(GLB_GET,itest,testunit)
                  if ( lbondgraph ) then
c                    --- for a bondgraph molecule, see if the testunit
c                    --- has any neighbors that exist
                     call twh_invib(GLB_GET,imolty,testunit,nvib)
                     lfound = .false.
                     ivib = 0
                     do while ( .not. lfound .and. ivib .lt. nvib )
                        ivib = ivib + 1
                        call twh_ijvib(GLB_GET,imolty,testunit,ivib
     &                       ,junit)
                        call twh_lexist(GLB_GET,junit,lfound)
                     enddo
                  else
c                    --- for a non bondgraph molecule always place in
c                    --- the same reference as cenunit
                     junit = cenunit
                     lfound = .true.
                  endif

                  if ( lfound ) then
c                    --- junit is a neighbor that exists, put this atom
c                    --- into the same reference frame
                     call twh_coordinates(GLB_GET,iflag,testchain
     &                    ,testunit,xtemp,ytemp,ztemp)
                     call twh_coordinates(GLB_GET,iflag,testchain
     &                    ,junit,rxui,ryui,rzui)
c                    --- compute the non minimum imaged distance
                     xold = xtemp - rxui
                     yold = ytemp - ryui
                     zold = ztemp - rzui
                     altdist = twh_distance(xold,yold,zold)
c                    --- compute the minimum imaged distance
                     call twh_mimage(ibox,xtemp,ytemp,ztemp
     &                    ,rxui,ryui,rzui,rxuij,ryuij,rzuij)
                     shortest = twh_distance(rxuij,ryuij,rzuij)
c                    --- compare the two distances
                     if ( (altdist - shortest) .gt. tolerance ) then
c                        write(6,*) 'CTRMAS: adjusting initial'
c     &                       ,' coordinates as molecule spans the'
c     &                       ,' periodic boundaries in an odd manner'
c                        write(6,*) 'Chain number: ',testchain
c                        write(6,*) 'Unit number: ',testunit
c                        write(6,*) 'pre-coords:',xtemp,ytemp,ztemp
c                       --- set back into the junit frame
                        xtemp = rxui + rxuij
                        ytemp = ryui + ryuij
                        ztemp = rzui + rzuij
                        call twh_coordinates(GLB_SET,iflag,testchain
     &                       ,testunit,xtemp,ytemp,ztemp)
c                        write(6,*) 'post-coords:',xtemp,ytemp,ztemp
                     endif
c                    --- this atom now exists in the proper frame
                     call twh_lexist(GLB_SET,testunit,.true.)
                     lreset = .true.
                  endif
                  if ( nvib .eq. 0 ) then
c                    --- this atom is not bonded to any others
c                    --- remove it from the list
                     lreset = .true.
                  endif
                  if ( lreset ) then
c                    --- remove this atom from the glist and replace
c                    --- with the end of the list
                     call twh_glist(GLB_GET,itogo,endunit)
                     call twh_glist(GLB_SET,itest,endunit)
                     itogo = itogo - 1
c                    --- reset itest to 0 to start searching again
                     itest = 0
                     lreset = .false.
                  endif
               enddo
            endif

c           --- it is possible that this molecule is bonded through
c           --- the periodic boundaries in such a way that it is 
c           --- impossible to put all of the atoms in the same
c           --- reference frame, and therefore its center-of-mass is
c           --- not well defined
            lperiodic = .false.
            iunit = 0
            do while ( .not. lperiodic .and. iunit .lt. nunit )
               iunit = iunit + 1
c              --- see if any of the atoms are still bonded through
c              --- a periodic boundary.  If they are then the molecule
c              --- is currently bonded in such a manner that it has
c              --- an infinite extent
c              --- only need to test ring atoms
               call twh_inring(GLB_GET,imolty,iunit,testring)
               if ( testring .ne. 0 ) then
c                 --- fetch the coordinates for iunit
                  call twh_coordinates(GLB_GET,iflag,testchain,iunit
     &                 ,rxui,ryui,rzui)
                  call twh_invib(GLB_GET,imolty,iunit,nvib)
                  do ivib = 1,nvib
                     call twh_ijvib(GLB_GET,imolty,iunit,ivib,junit)
                     call twh_coordinates(GLB_GET,iflag,testchain,junit
     &                    ,xtemp,ytemp,ztemp)
c                    --- compute the direct bond length
                     xold = xtemp - rxui
                     yold = ytemp - ryui
                     zold = ztemp - rzui
                     altdist = twh_distance(xold,yold,zold)
c                    --- compute the minimum imaged distance
                     call twh_mimage(ibox,xtemp,ytemp,ztemp
     &                    ,rxui,ryui,rzui,rxuij,ryuij,rzuij)
                     shortest = twh_distance(rxuij,ryuij,rzuij)
c                    --- compare the two distances
                     if ( (altdist - shortest) .gt. tolerance ) then
c                       --- this molecule is periodically bonded
                        lperiodic = .true.
                     endif
                  enddo
               endif
            enddo

c           ----- Determine new center of mass for chain testchain
            nxcm = 0.0d0
            nycm = 0.0d0
            nzcm = 0.0d0
            if ( lperiodic ) then
c              --- set the "COM" to the middle of the box (this is used
c              --- for some of the energy calcaultions to save time)
               arba = 0.50d0
               arbb = 0.50d0
               arbc = 0.50d0
               call twh_arbtocart(ibox,arba,arbb,arbc,nxcm,nycm,nzcm)
            else
               do iunit = 1, nunit
                  call twh_ntype(GLB_GET,imolty,iunit,itype)
                  call twh_coordinates(GLB_GET,iflag,testchain,iunit
     &                 ,xtemp,ytemp,ztemp)
                  if ( masstotal .gt. 0.0d0 ) then
c                    --- mass weighting
                     nxcm = nxcm + xtemp*twh_get_mass(itype)
                     nycm = nycm + ytemp*twh_get_mass(itype)
                     nzcm = nzcm + ztemp*twh_get_mass(itype)
                  else
c                    --- determine the positional center without regard
c                    --- to mass
                     nxcm = nxcm + xtemp
                     nycm = nycm + ytemp
                     nzcm = nzcm + ztemp
                  endif
               enddo
               if ( masstotal .gt. 0.0d0 ) then
c                 --- divide out the total mass
                  nxcm = nxcm / masstotal
                  nycm = nycm / masstotal
                  nzcm = nzcm / masstotal
               endif
            endif

            if ( lrescale ) then
c              --- initialize ladjust to false
               ladjust = .false.
               if ( lperiodic ) then
c                 --- place each atom individually into the central box
                  do iunit = 1,nunit
                     call twh_coordinates(GLB_GET,iflag,testchain
     &                    ,iunit,xtemp,ytemp,ztemp)
                     rxui = xtemp
                     ryui = ytemp
                     rzui = ztemp
                     call twh_putcartinbox(ibox,rxui,ryui,rzui)
                     dx = rxui - xtemp
                     dy = ryui - ytemp
                     dz = rzui - ztemp
                     call twh_com_safecompare(dx,dy,dz)
                     if ( dx*dx .gt. tolerance ) then
                        lxx = .true.
                     else
                        lxx = .false.
                     endif
                     if ( dy*dy .gt. tolerance ) then
                        lyy = .true.
                     else
                        lyy = .false.
                     endif
                     if ( dz*dz .gt. tolerance ) then
                        lzz = .true.
                     else
                        lzz = .false.
                     endif

                     if ( lxx .or. lyy .or. lzz ) then
                        if ( lwarn ) then
c                           write(6,*) 'CTRMAS: testchain,iunit:'
c     &                          ,testchain,iunit
c                           write(6,*) 'Initial coords:',xtemp,ytemp
c     &                          ,ztemp
                        endif
                        ladjust = .true.
                        xtemp = rxui
                        ytemp = ryui
                        ztemp = rzui
c                       --- set the new coordinates
                        call twh_coordinates(GLB_SET,iflag,testchain
     &                       ,iunit,xtemp,ytemp,ztemp)
                        if ( lwarn ) then
c                           write(6,*) 'New coords:',xtemp,ytemp,ztemp
                        endif
                     endif
                  enddo
               else
c                 --- move the entire molecule based upon the COM
                  xold = nxcm
                  yold = nycm
                  zold = nzcm
c                 --- find the coordinates inside the central box
                  call twh_putcartinbox(ibox,nxcm,nycm,nzcm)
c                 --- find the difference of the inside and old coords
                  dx = nxcm - xold
                  dy = nycm - yold
                  dz = nzcm - zold
                  call twh_com_safecompare(dx,dy,dz)
               
                  if ( dx*dx + dy*dy + dz*dz .gt. tolerance ) then
                     ladjust = .true.
                     do iunit = 1,nunit
                        call twh_coordinates(GLB_GET,iflag,testchain
     &                       ,iunit,xtemp,ytemp,ztemp)
                        xtemp = xtemp + dx
                        ytemp = ytemp + dy
                        ztemp = ztemp + dz
                        call twh_coordinates(GLB_SET,CRD_REAL,testchain
     &                       ,iunit,xtemp,ytemp,ztemp)
                     enddo
                  endif
               endif

               if( ladjust .and. lwarn ) then
                  if ( mtype .eq. CTR_SWAP ) then
                     write(6,*) 'swap move'
                  endif
                  if ( mtype .eq. CTR_TRANATOM ) then
                     write(6,*) 'translate ATOM move'
                  endif
                  if ( mtype .eq. CTR_INITIAL ) then
                     write(6,*) 'initial readtowhee'
                  endif
                  if ( mtype .eq. CTR_CBREGROW ) then
                     write(6,*) 'cbregrow move'
                  endif
                  if ( mtype .eq. CTR_ENERGY ) then
                     write(6,*) 'engmolec call'
                  endif
                  if ( mtype .eq. CTR_INITCONF ) then
                     write(6,*) 'initconf'
                  endif
                  if ( mtype .eq. CTR_CONROT ) then
                     write(6,*) 'conrot move'
                  endif
                  if ( mtype .eq. CTR_PIVOT ) then
                     write(6,*) 'pivot move'
                  endif
                  if ( .not. lperiodic ) then
                     write(6,'(1x,a,3(1x,i12))') 'ibox,testchain,nunit'
     &                    ,ibox,testchain,nunit
                     write(6,'(1x,a,3(1x,g12.6))') 'xold,yold,zold'
     &                    ,xold,yold,zold
                     write(6,'(1x,a,3(1x,g12.6))') 'dx,dy,dz',dx,dy,dz
                     if ( mtype .ne. CTR_INITIAL ) then
                        lfailure = .true.
                        return
                     endif
                  endif
               endif
            endif
c           --- assign the new center of mass
            call twh_com(GLB_SET,istorage,testchain,nxcm,nycm,nzcm)
c           --- assign the new maximum distance from any atom to the
c           --- center of mass
            if ( lnewrcmu ) then
               dmaxsq = 0.0d0 
c              --- initialize radius of gyration
               gyration = 0.0d0
               if ( lperiodic ) then
c                 --- these values are not well defined for a periodic
c                 --- molecule
c                 --- default to a large value of rcmu for the energy
c                 --- routines
                  rcmu = 10.0d0*twh_maxboxlength(ibox)
c                 --- leave radius of gyration as zero
               else
                  do iunit=1,nunit
c                    --- minimum image the ctrmas pair separations
                     call twh_coordinates(GLB_GET,iflag,testchain,iunit
     &                    ,rxui,ryui,rzui)
                     call twh_mimage ( ibox,rxui,ryui,rzui
     &                    ,nxcm,nycm,nzcm,rxuij,ryuij,rzuij )
                     rijsq = rxuij*rxuij + ryuij*ryuij  + rzuij*rzuij
c                    --- test if this is the largest distance
                     if ( rijsq .gt. dmaxsq ) dmaxsq = rijsq
c                    --- update radius of gyration
                     call twh_ntype(GLB_GET,imolty,iunit,itype)
                     gyration = gyration + rijsq*twh_get_mass(itype)
                  enddo                  
                  rcmu = dsqrt(dmaxsq)+ 1.0d-5
c                 --- divide out the mass for the Rg
                  if ( masstotal .gt. 0.0d0 ) then
                     gyration = dsqrt(gyration / masstotal)
                  endif
               endif
c              --- store the new periodic flag
               call twh_lperiodic(GLB_SET,testchain,lperiodic)
c              --- store the new maximum distance
               call twh_rcmu(GLB_SET,testchain,rcmu)
c              --- store the new radius of gyration
               call twh_gyration(GLB_SET,testchain,gyration)
            endif
         else
            if ( iall .ne. 0 ) write(6,*)'prob with box in ctrmas'
         endif
      enddo

c      write(6,*) 'end CTRMAS',ibox
      
      return
      end

      subroutine twh_com_safecompare(dx,dy,dz)
c     ******************************************************************
c     * checks to see if the distance changes are significant          *
c     *                                                                *
c     * originally split from ctrmas 12-15-2007 by M.G. Martin         *
c     * last modified 12-15-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_EQ
#include "functions.h"
c     --- variables passed to/from subroutine
      double precision dx,dy,dz

#ifdef SAFE_COMPARE
c     --- if d* are not numerically different from zero then
c     --- set them to zero
      if (twh_cmp_eq(dx, 0.0d0)) dx = 0.0d0
      if (twh_cmp_eq(dy, 0.0d0)) dy = 0.0d0
      if (twh_cmp_eq(dz, 0.0d0)) dz = 0.0d0
#endif

      return
      end
