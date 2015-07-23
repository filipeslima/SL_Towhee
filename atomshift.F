#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2008 Marcus G. Martin                       *
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
      subroutine twh_atomshift(lfailure,movetype)
c     ******************************************************************
c     * wrapper around the atomshift move that helps set up some data  *
c     * structures                                                     *
c     *                                                                *
c     * originally written 06-25-2009 by M.G. Martin                   *
c     * last modified 06-25-2009 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from subroutine
      logical lfailure
      integer movetype
c     --- local variables
      integer nchain
      call twh_nchain(GLB_GET,nchain)

      call twh_atomshift_main(lfailure,movetype,nchain)

      return
      end

      subroutine twh_atomshift_main(lfailure,movetype,nchain)
c     ******************************************************************
c     * movetype=1: finds a new position for an entire plane of atoms. *
c     * This shift is done according to the arbitrary coordinate       *
c     * systems. a new center of the plane is chosen at random from the*
c     * entire box length.                                             *
c     * movetype=2: finds a new position for an entire row of atoms    *
c     *                                                                *
c     * originally written as atomshift 01-10-2002 by M.G. Martin      *
c     * renamed atomshift_main 06-25-2009 by M.G. Martin               *
c     * last modified 04-21-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_EXPON
#define FUNCTION_GET_NBOXI
#define FUNCTION_RANDOM
#include "functions.h"
c     --- debug toggle
#define DEBUG_ATOMSHIFT 0
c     --- variables passed to/from the subroutine
      logical lfailure
      integer movetype,nchain
c     --- local variables
c     --- logical scalars
      logical lovrlap,lfound,laccept,lovrold,lewald
c     --- logical arrays
      logical lfnd
      dimension lfnd(2)
c     --- integer scalars
      integer arbindex,ibox,ichain,icond
      integer imolty,iunit,idim,imove,tempdim
      integer ntomove,nunit
c     --- integer arrays
      integer ncond
      dimension ncond(2)
      integer tardim
      dimension tardim(2)
      integer movelist
      dimension movelist(nchain)
c     --- double precision scalars
      double precision blength,rm,upper,lower,rcenter,range
      double precision xtemp,ytemp,ztemp,xdisp,ydisp,zdisp,burn
      double precision vnewtotal,vnewinter,vnewintra,vnewext,vnewelect
      double precision vnewquantum,voldquantum
      double precision voldtotal,voldinter,voldintra,voldext,voldelect
      double precision vduma,vdumb,vdumc,vdumd,dele,ptest,vdumz
      double precision vnewsolv,voldsolv,xxx,yyy,zzz
      double precision newdvl,newdvc,olddvl,olddvc
      double precision xcmi,ycmi,zcmi,beta,dv,planewidth,rowwidth
      double precision hm1,hm2,hm3
c     --- double precision arrays
      double precision bounds
      dimension bounds(2,2,2)
      double precision arbcord,displace
      dimension arbcord(3),displace(3)
      double precision oldvf,newvf
      dimension oldvf(MAX_FOREIGN_LAMBDA),newvf(MAX_FOREIGN_LAMBDA)

#if DEBUG_ATOMSHIFT
      write(6,*) 'start ATOMSHIFT'
#endif
c     --- retrieve constants
      call twh_lewald(GLB_GET,lewald)

c     --- select a box
      rm = twh_random()
      lfound = .false.
      ibox = 0
      do while ( .not. lfound ) 
         ibox = ibox + 1
         if ( movetype .eq. 1 ) then
            call twh_pmplanebox(GLB_GET,ibox,ptest)
         elseif ( movetype .eq. 2 ) then
            call twh_pmrowbox(GLB_GET,ibox,ptest)
         else
            write(6,*) 'ATOMSHIFT: invalid movetype'
            lfailure = .true.
            return
         endif
         if ( rm .lt. ptest ) lfound = .true.
      enddo

c     --- select a dimension for translation
c     --- in the case of plane shift this is the dimension that is 
c     --- perpendicular to the plane shift
c     --- in the case of row shift this is the dimension that is 
c     --- parallel to the row
      arbindex = int(3.0d0*twh_random()) + 1

c     --- increment number of attempted moves
      call twh_bnshift(GLB_INCR,movetype,ibox,arbindex,1)

      if ( movetype .eq. 1 ) then
c        --- compute the displacements of the atoms of the plane
c        --- this is in arbitrary units so this displacement is 
c        --- equivalent to selecting a new position for the center of 
c        --- the plane uniformly from the availble choices
         if ( arbindex .eq. 1 ) then
            displace(1) = 0.0d0
            displace(2) = twh_random()
            displace(3) = twh_random()
         elseif ( arbindex .eq. 2 ) then
            displace(1) = twh_random()
            displace(2) = 0.0d0
            displace(3) = twh_random()
         else
            displace(1) = twh_random()
            displace(2) = twh_random()
            displace(3) = 0.0d0
         endif

c        --- select a position at random
         rcenter = twh_random()
c        --- get the boxlength of the target dimension
         call twh_hmatrix(GLB_GET,ibox,arbindex,1,hm1)
         call twh_hmatrix(GLB_GET,ibox,arbindex,2,hm2)
         call twh_hmatrix(GLB_GET,ibox,arbindex,3,hm3)
         blength = twh_distance(hm1,hm2,hm3)
c        --- compute the range
         call twh_planewidth(GLB_GET,planewidth)
         range = 0.5d0*planewidth/blength
c        --- compute upper and lower
         upper = rcenter + range
         lower = rcenter - range
         if ( upper .gt. 1.0d0 .and. lower .lt. 0.0d0 ) then
c           --- somehow we are spanning the entire simulation box
c           --- this really is not a good idea as the translation has
c           --- no purpose if we move every atom
            write(6,*) 'ibox,planewidth,blength',ibox,planewidth,blength
            write(6,*) 'you have choosen a planewidth so large that it '
            write(6,*) 'covers every atom in the box'
            write(6,*) 'such a move serves no purpose'
            lfailure = .true.
            return
         elseif ( upper .gt. 1.0d0 ) then
c           --- need to check two different conditions
            ncond(1) = 2
            bounds(1,1,1) = lower
            bounds(1,1,2) = 1.0d0
            bounds(1,2,1) = 0.0d0
            bounds(1,2,2) = upper - 1.0d0
         elseif ( lower .lt. 0.0d0 ) then
c           --- need to check two different conditions
            ncond(1) = 2
            bounds(1,1,1) = 0.0d0
            bounds(1,1,2) = upper
            bounds(1,2,1) = lower + 1.0d0
            bounds(1,2,2) = 1.0d0
         else
c           --- just one condition to check
            ncond(1) = 1
            bounds(1,1,1) = lower
            bounds(1,1,2) = upper
         endif
      elseif ( movetype .eq. 2 ) then
c        --- compute the displacements of the atoms of the row
c        --- this is in arbitrary units so this displacement is 
c        --- equivalent to selecting a new position for the center of 
c        --- the row uniformly from the availble choices
         do idim = 1,3
            displace(idim) = 0.0d0
         enddo
c        --- displace this index
         displace(arbindex) = twh_random()

         if ( arbindex .eq. 1 ) then
            tardim(1) = 2
            tardim(2) = 3
         elseif ( arbindex .eq. 2 ) then
            tardim(1) = 1
            tardim(2) = 3
         else
            tardim(1) = 1
            tardim(2) = 2
         endif

         call twh_rowwidth(GLB_GET,rowwidth)
         do idim = 1,2
c           --- select a position at random
            rcenter = twh_random()
c           --- store the target dimension
            tempdim = tardim(idim)
c           --- get the boxlength of the target dimension
            call twh_hmatrix(GLB_GET,ibox,tempdim,1,hm1)
            call twh_hmatrix(GLB_GET,ibox,tempdim,2,hm2)
            call twh_hmatrix(GLB_GET,ibox,tempdim,3,hm3)
            blength = twh_distance(hm1,hm2,hm3)
c           --- compute the range
            range = 0.5d0*rowwidth/blength
c           --- compute upper and lower
            upper = rcenter + range
            lower = rcenter - range
            if ( upper .gt. 1.0d0 .and. lower .lt. 0.0d0 ) then
c              --- somehow we are spanning the entire simulation box
c              --- this really is not a good idea as the translation has
c              --- no purpose if we move every atom
               write(6,*) 'ibox,rowwidth,blength',ibox,rowwidth,blength
               write(6,*) 'you have choosen a rowwidth so large that it'
               write(6,*) 'covers every atom in the box'
               write(6,*) 'such a move serves no purpose'
               lfailure = .true.
               return
            elseif ( upper .gt. 1.0d0 ) then
c              --- need to check two different conditions
               ncond(idim) = 2
               bounds(idim,1,1) = lower
               bounds(idim,1,2) = 1.0d0
               bounds(idim,2,1) = 0.0d0
               bounds(idim,2,2) = upper - 1.0d0
            elseif ( lower .lt. 0.0d0 ) then
c              --- need to check two different conditions
               ncond(idim) = 2
               bounds(idim,1,1) = 0.0d0
               bounds(idim,1,2) = upper
               bounds(idim,2,1) = lower + 1.0d0
               bounds(idim,2,2) = 1.0d0
            else
c              --- just one condition to check
               ncond(idim) = 1
               bounds(idim,1,1) = lower
               bounds(idim,1,2) = upper
            endif
         enddo
      endif

c     --- make a list of the molecules to move based on the COM position
      ntomove = 0
      do ichain = 1,nchain
c        --- initialize lmoveme
         call twh_lmoveme(GLB_SET,ichain,.false.)
         if ( twh_get_nboxi(ichain) .eq. ibox ) then
c           --- transform the COM into arbitrary units
            call twh_com(GLB_GET,CRD_REAL,ichain,xcmi,ycmi,zcmi)
            call twh_carttoarb(ibox,xcmi,ycmi,zcmi,arbcord(1),arbcord(2)
     &           ,arbcord(3))

c           --- check the bounds
            if ( movetype .eq. 1 ) then 
               lfound = .false.
               do icond = 1,ncond(1)
                  if ( bounds(1,icond,1) .le. arbcord(arbindex) 
     &                 .and. bounds(1,icond,2) .ge. arbcord(arbindex) )
     &                 lfound = .true.
               enddo
            elseif ( movetype .eq. 2 ) then
               lfound = .true.
               do idim = 1,2
                  lfnd(idim) = .false.
                  do icond = 1,ncond(idim)
                     if ( bounds(idim,icond,1) 
     &                    .le. arbcord(tardim(idim))
     &                    .and. arbcord(tardim(idim)) 
     &                    .le. bounds(idim,icond,2)) lfnd(idim) = .true.
                  enddo
                  if ( .not. lfnd(idim) ) lfound = .false.
               enddo
            endif

            if ( lfound ) then
c              --- add this molecule to the list
               ntomove = ntomove + 1
               movelist(ntomove) = ichain
               call twh_lmoveme(GLB_SET,ichain,.true.)
c              --- store the old COM in cart coords
               call twh_com(GLB_STORE,CRD_R_TO_T,ichain,xcmi,ycmi,zcmi)
c              --- apply the displacement to the COM
               do idim = 1,3
                  arbcord(idim) = arbcord(idim) + displace(idim)
               enddo
c              --- put the new COM back into the central box
               call twh_putarbinbox(arbcord(1),arbcord(2),arbcord(3))
c              --- get the new COM in cartesian coords
               call twh_arbtocart(ibox,arbcord(1),arbcord(2),arbcord(3)
     &              ,xtemp,ytemp,ztemp)
c              --- compute the displacement of the cart COM
               call twh_com(GLB_GET,CRD_TEMP,ichain,xcmi,ycmi,zcmi)
               xdisp = xtemp - xcmi
               ydisp = ytemp - ycmi
               zdisp = ztemp - zcmi
c              --- set the new COM
               call twh_com(GLB_SET,CRD_REAL,ichain,xtemp,ytemp,ztemp)
c              --- compute molecule type
               call twh_moltyp(GLB_GET,ichain,imolty)
               call twh_nunit(GLB_GET,imolty,nunit)
               do iunit = 1,nunit
c                 --- store the old cartesian coords for each unit
                  call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &                 ,xtemp,ytemp,ztemp)
                  call twh_coordinates(GLB_SET,CRD_TEMP,ichain,iunit
     &                 ,xtemp,ytemp,ztemp)
c                 --- apply the displacement to get the trial coords
                  xtemp = xtemp + xdisp
                  ytemp = ytemp + ydisp
                  ztemp = ztemp + zdisp
c                 --- set up the new trial coordinates
                  call twh_coordinates(GLB_SET,CRD_REAL,ichain,iunit
     &                 ,xtemp,ytemp,ztemp)
               enddo
            endif
         endif
      enddo

      if ( ntomove .eq. 0 ) then
c        --- no atoms to move, nothing has happened, just return
         return
      endif

c     --- store the old values for the ewald sum
      if ( lewald ) call twh_recip(ibox,3)

c     --- compute the new energy 
      lovrlap = .false.
      call twh_engtotal(lfailure,lovrlap,vnewtotal,vnewinter,vdumz,vduma
     &     ,vnewintra,vdumb,vdumc,vdumd,vnewext,vnewelect,vnewsolv
     &     ,vnewquantum,newvf,newdvl,newdvc,ibox,.true.,.true.,.false.
     &     ,.false.)
      if ( lfailure ) return
c     --- reset all of the coordinates to their old values and then
c     --- set *cordtmp to the new values
      do imove = 1,ntomove
         ichain = movelist(imove)
         call twh_moltyp(GLB_GET,ichain,imolty)
c        --- exchange the two storage arrays
         call twh_com(GLB_EXCHANGE,CRD_R_TO_T,ichain,xcmi,ycmi,zcmi)
c        --- now reset all of the atoms
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            call twh_coordinates(GLB_GET,CRD_TEMP,ichain,iunit
     &           ,xtemp,ytemp,ztemp)
c           --- now reset the *cordtmp variables to the new coords
            call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &           ,xxx,yyy,zzz)
            call twh_coordinates(GLB_SET,CRD_TEMP,ichain,iunit
     &           ,xxx,yyy,zzz)
c           --- set the r*u variables to the old coordinates
            call twh_coordinates(GLB_SET,CRD_REAL,ichain,iunit
     &           ,xtemp,ytemp,ztemp)
         enddo
      enddo

      if ( lovrlap ) then
c        --- we had an overlap in the new conformation, move rejected
c        --- set the ewald stuff back to the old values and exit
         if ( lewald ) call twh_recip(ibox,4)
         return
      else
c        --- we are about to compute the old energy, this will recompute
c        --- all of the ewald terms so now we want to store the new
c        --- ewald terms in the same way we just stored the old above
         if ( lewald ) call twh_recip(ibox,3)
      endif

c     --- we need to compute the old energy to determine acceptance
c     --- compute the old energy
      call twh_engtotal(lfailure,lovrold,voldtotal,voldinter,vdumz,vduma
     &     ,voldintra,vdumb,vdumc,vdumd,voldext,voldelect,voldsolv
     &     ,voldquantum,oldvf,olddvl,olddvc,ibox,.true.,.true.,.false.
     &     ,.false.)
      if ( lfailure ) return
      if ( lovrold ) then
         write(6,*) 'ATOMSHIFT: disaster overlap in old conformation'
         lfailure = .true.
         return
      endif

c     --- determine acceptance
      laccept = .false.
c     --- compute the energy change, use the old total energy without
c     --- the bonded parts
      dele = vnewtotal - voldtotal

      if ( dele .lt. 0.0d0 ) then
         laccept = .true.
c        --- burn a random number
         burn = twh_random()
      else
         call twh_beta(GLB_GET,beta)
         if ( twh_random() .lt. twh_expon(-beta*dele) ) then
            laccept = .true.
         endif
      endif

#if DEBUG_ATOMSHIFT
      write(6,*) 'dele:',dele
      write(6,*) 'laccept:',laccept
#endif
      if ( laccept ) then
c        --- move accepted, update counter and energies
         call twh_bashift(GLB_INCR,movetype,ibox,arbindex,1)
         dv = vnewtotal - voldtotal
         call twh_boxvclassic(GLB_INCR,ibox,dv)
         dv = vnewinter - voldinter
         call twh_boxvinter(GLB_INCR,ibox,dv)
         dv = vnewintra - voldintra
         call twh_boxvintra(GLB_INCR,ibox,dv)
         dv = vnewext - voldext
         call twh_boxvfield(GLB_INCR,ibox,dv)
         dv = vnewelect - voldelect
         call twh_boxvelect(GLB_INCR,ibox,dv)
         dv = vnewsolv - voldsolv
         call twh_boxvsolvation(GLB_INCR,ibox,dv)

c        In principle, the new values, as they're calculated by
c        engtotal, should have correct new energies and derivs.  
         call scp_set_derivative(ibox, LAMBDA_C, newdvc)
         call scp_set_derivative(ibox, LAMBDA_LJ, newdvl)

         call scp_set_foreign_energies(ibox, newvf)

c        --- reset all of the coordinates to the new values 
c        --- they were back to the old values for the energy 
c        --- computation
c        --- so even though it looks like we are setting these to the
c        --- old values, we are really setting them to the new
         if ( lewald ) call twh_recip(ibox,4)
         do imove = 1,ntomove
            ichain = movelist(imove)
c           --- reset COM to the temp values
            call twh_com(GLB_STORE,CRD_T_TO_R,ichain,xcmi,ycmi,zcmi)
c           ---- reset atoms
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
c               --- get old coordinates
               call twh_coordinates(GLB_GET,CRD_TEMP,ichain,iunit,xtemp,
     &              ytemp,ztemp)
c              --- restore old coordinates
               call twh_coordinates(GLB_SET,CRD_REAL,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
            enddo
         enddo
      endif

#if DEBUG_ATOMSHIFT
      write(6,*) 'end ATOMSHIFT'
#endif
      return
      end
