#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2009 Marcus G. Martin                       *
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
      subroutine twh_composite(lfailure)
c     ******************************************************************
c     * Perform MC move which is a composite of translation and        *
c     * rotation                                                       *
c     * the maximum displacement is controlled by rmcomp               *
c     * number of attempts stored in acncomp                           *
c     * number of successful trial moves is stored in acscomp          *
c     *                                                                *
c     * based closely on the code trancom.F and rotate.F               *
c     * originally written 03-23-2005 by MAW                           *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EXPON
#define FUNCTION_RANDOM
#define FUNCTION_MINBOXLENGTH
#include "functions.h"
c     --- debugging variable
#define DEBUG_COMPOSITE 0
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
c     --- logical scalars
      logical laccept,lgetcom,lintra,lbonded,lpartial,lperiodic
      logical scp_eval_derivative,lewald
c     --- logical arrays
      logical loverlap
      dimension loverlap(2)
c     --- integer scalars
      integer ichain,ibox,totunits,imolty,icbu,nchainmove,iunit,nmolty
      integer ncmt,ipoint,nmtype,potentialstyle,ensemble
c     --- integer arrays
      integer chaininfo
      dimension chaininfo(1,3)
      integer iboxa
      dimension iboxa(1)
c     --- double precision scalars
      double precision rx,ry,rz,dchain,ddx,ddy,ddz,rchain,burn,beta
      double precision xtemp,ytemp,ztemp
      double precision dgammax,dgammay,dgammaz
      double precision dvclassical_trial,dvintra_trial,dvinter_trial
      double precision dvext_trial,dvelect_trial,dvsolv_trial
      double precision dvtor_trial,dvbend_trial,dvvib_trial
      double precision scaledist,dvquantum_trial,dv_trial
      double precision ddvdl_lj,ddvdl_c,rmcomrot
      double precision rxorig,ryorig,rzorig,rcmu,pmtest,rmcomtra
c     --- double precision arrays
      double precision dvf
      dimension dvf(MAX_FOREIGN_LAMBDA)

#if DEBUG_COMPOSITE
      write(6,*) 'start COMPOSITE'
#endif
c     --- retrieve constants
      call twh_lewald(GLB_GET,lewald)
c     --- get the potentialstyle
      call twh_potentialstyle(GLB_GET,potentialstyle)
c     --- select a molecule type at random
      rchain  = twh_random()
      call twh_nmolty(GLB_GET,nmolty)
      do icbu = 1,nmolty
         call twh_pmcomt(GLB_GET,icbu,pmtest)
         if ( rchain .lt. pmtest ) then
            imolty = icbu
            rchain = 2.0d0
         endif
      enddo

c     --- for grand canonical ensemble only want to perform move on
c     --- molecules actually in the 'real' system box (box = 1)
      call twh_ensemble(GLB_GET,ensemble)
      if (ensemble .eq. ENS_muVT) then
         ibox = 1         
c        --- see if there are any molecules of this type in active box
         call twh_ncmt(GLB_GET,ibox,imolty,ncmt)
         if (ncmt .eq. 0 ) return
c        --- select a chain of the proper type from this box
         dchain = dble(ncmt)
         ipoint = int(dchain*twh_random()+1)
         call twh_chainlist(GLB_GET,ipoint,ibox,imolty,ichain)
      else
         call twh_nmolectyp(GLB_GET,imolty,nmtype)
         dchain = dble(nmtype)
         ipoint = int( dchain*twh_random() + 1 )
c        --- see note about parall in readtowhee.F
c        --- below, ichain is index of randomly selected chain of type
c        --- imolty.
         call twh_parall(GLB_GET,imolty,ipoint,ichain)
      endif
      call twh_nboxi(GLB_GET,ichain,ibox)
      if ( ibox .eq. 0 ) then
         write(6,*) 'COMPOSITE: problem selecting proper box'
         lfailure = .true.
         return
      endif
c     --- increment move attempt counter
      call twh_acncomp(GLB_INCR,imolty,ibox,1.0d0)
c     --- check to make sure this is not a periodically bonded molecule
      call twh_lperiodic(GLB_GET,ichain,lperiodic)
      if ( lperiodic ) then
         write(6,*) 'COMPOSITE: attempted move on a periodically bonded'
     &        ,' molecule'
         write(6,*) 'Molecule type:',imolty
         write(6,*) 'Chain number:',ichain
         return
      endif

c     --- store number of units of ichain in totunits 
      call twh_nunit(GLB_GET,imolty,totunits)

c     --- store the current atomic coordinates as the trial coords
      do iunit = 1,totunits
         call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &        ,xtemp,ytemp,ztemp)
         call twh_coordinates(GLB_SET,CRD_TEMP,ichain,iunit
     &        ,xtemp,ytemp,ztemp)
      enddo

c     --- perform a rotation
      call twh_rmcomrot(GLB_GET,imolty,ibox,rmcomrot)
c     --- choose a random angular displacement over the 3 axes
      dgammax = ( 2.0d0*twh_random() - 1.0d0 )* rmcomrot
      dgammay = ( 2.0d0*twh_random() - 1.0d0 )* rmcomrot
      dgammaz = ( 2.0d0*twh_random() - 1.0d0 )* rmcomrot

c     --- Use the center of mass for rotation
      call twh_com(GLB_GET,CRD_REAL,ichain,rxorig,ryorig,rzorig)
c     --- Determine the rotation coordinates
      call twh_rotationmatrix(CRD_TEMP,ichain,rxorig,ryorig,rzorig
     &     ,dgammax,dgammay,dgammaz)

c     --- perform a translation
c     --- translate on the unit sphere (rx,ry,rz reset here)
      call twh_sphere(rx,ry,rz)
c     --- compute a random scale distance between 0 and maximum disp
      call twh_rmcomtra(GLB_GET,imolty,ibox,rmcomtra)
      scaledist = twh_random()*rmcomtra
c     --- scale by the maximum displacement
      rx = rx*scaledist
      ry = ry*scaledist
      rz = rz*scaledist

c     --- find out if we need to displace the atoms to keep them in 
c     --- the central box
      ddx = rxorig + rx
      ddy = ryorig + ry
      ddz = rzorig + rz
      call twh_putcartinbox(ibox,ddx,ddy,ddz)
c     --- if the coordinates were shifted then compute the shift
      ddx = ddx - rxorig
      ddy = ddy - ryorig
      ddz = ddz - rzorig

c     --- update the move coordinates, which already have been modified
c     --- by rotation
      do iunit = 1,totunits
         call twh_coordinates(GLB_GET,CRD_TEMP,ichain,iunit
     &        ,xtemp,ytemp,ztemp)
         xtemp = xtemp + ddx
         ytemp = ytemp + ddy
         ztemp = ztemp + ddz
         call twh_coordinates(GLB_SET,CRD_TEMP,ichain,iunit
     &        ,xtemp,ytemp,ztemp)
      enddo

      call twh_rcmu(GLB_GET,ichain,rcmu)
      if ( rcmu .gt. twh_minboxlength(ibox)*0.25d0 ) then
c        --- molecule is large enough that it might be interacting with
c        --- iself through the periodic images (occurs for some very 
c        --- short range potentials like hard spheres).   Need to 
c        --- compute the intramolecular and bonded interactions as these
c        --- could change in this case
         lintra = .true.
         lbonded = .true.
      else
c        --- molecule does not interact with itself through the periodic
c        --- images.  don't bother computing intramolecular or bonded
         lintra = .false.
         lbonded = .false.
      endif

#if DEBUG_COMPOSITE
      write(6,*) 'COMPOSITE: proposed rotation: ', 
     &      dgammax, dgammay, dgammaz
      write(6,*) 'COMPOSITE: proposed translation: ', ddx, ddy, ddz
#endif

      nchainmove = 1
      chaininfo(1,CRD_OLD) = CRD_REAL
      chaininfo(1,CRD_NEW) = CRD_TEMP
      chaininfo(1,3) = ichain
      lbonded = .false.
      lgetcom = .true.
      lpartial = .false.

      call twh_energy_change(lfailure,ibox,nchainmove,chaininfo,loverlap
     &     ,lpartial,lgetcom,lintra,lbonded,dvclassical_trial
     &     ,dvintra_trial,dvinter_trial,dvext_trial,dvelect_trial
     &     ,dvsolv_trial,dvtor_trial,dvbend_trial
     &     ,dvvib_trial,dvquantum_trial,dvf,ddvdl_lj,ddvdl_c)
      if ( lfailure ) return

c     --- check new hard overlap, if overlap move rejected
      if ( loverlap(CRD_NEW) ) then
#if DEBUG_COMPOSITE
         write(6,*) "COMPOSITE: move rejected due to overlap"
#endif
         return
      endif
c     --- check old hard overlap, if overlap then big problems
      if (loverlap(CRD_OLD)) then
         write(6,*) 'COMPOSITE: disaster loverlap in old conf'
         lfailure = .true.
         return
      endif

      if ( potentialstyle .eq. POT_INTERNAL ) then
         dv_trial = dvclassical_trial
      elseif ( potentialstyle .eq. POT_EXTERNAL ) then
         dv_trial = dvquantum_trial
      else
         write(6,*) 'COMPOSITE: unknown potentialstyle'
         write(6,*) potentialstyle
         lfailure = .true.
         return
      endif
c     --- check for acceptance
      laccept = .false.

      if ( dv_trial .le. 0.0d0 ) then
c        --- accept move
         laccept = .true.
c        --- burn a random number
         burn = twh_random()
      else
         call twh_beta(GLB_GET,beta)
         if ( twh_random() .lt. twh_expon(-beta*dv_trial) ) then
            laccept = .true.
         endif
      endif

#if DEBUG_COMPOSITE
      write(6,*)  'COMPOSITE: dv_trial ',dv_trial
#endif
      if ( .not. laccept ) return
#if DEBUG_COMPOSITE
      write(6,*) 'COMPOSITE: accepted chain',ichain
#endif
c     --- update acceptance accumulators
      call twh_acscomp(GLB_INCR,imolty,ibox,1.0d0)
c     --- update energies
      if ( potentialstyle .eq. POT_INTERNAL ) then
         call twh_boxvclassic(GLB_INCR,ibox,dvclassical_trial)
         call twh_boxvinter(GLB_INCR,ibox,dvinter_trial)
         call twh_boxvintra(GLB_INCR,ibox,dvintra_trial)
         call twh_boxvfield(GLB_INCR,ibox,dvext_trial)
         call twh_boxvelect(GLB_INCR,ibox,dvelect_trial)
         call twh_boxvsolvation(GLB_INCR,ibox,dvsolv_trial)
         call twh_boxvtorsion(GLB_INCR,ibox,dvtor_trial)
         call twh_boxvangle(GLB_INCR,ibox,dvbend_trial)
         call twh_boxvbond(GLB_INCR,ibox,dvvib_trial)
         call twh_scp_eval_derivative(GLB_GET,scp_eval_derivative)
         if (scp_eval_derivative) then
            call scp_increment_derivative(ibox, LAMBDA_C, ddvdl_c)
            call scp_increment_derivative(ibox, LAMBDA_LJ, ddvdl_lj)
         endif
         call scp_increment_foreign_energies(ibox, dvf)

      elseif ( potentialstyle .eq. POT_EXTERNAL ) then
         call twh_boxvquantum(GLB_INCR,ibox,dvquantum_trial)
      else
         write(6,*) 'COMPOSITE: unknown potentialstyle'
         write(6,*) potentialstyle
         lfailure = .true.
         return
      endif
c     --- store new coordinates
      do iunit = 1,totunits
         call twh_coordinates(GLB_GET,CRD_TEMP,ichain,iunit
     &        ,xtemp,ytemp,ztemp)
         call twh_coordinates(GLB_SET,CRD_REAL,ichain,iunit
     &        ,xtemp,ytemp,ztemp)
      enddo

      if ( lewald ) then
c        --- update reciprocal-space sum
         call twh_recip(ibox,2)
      endif

c     --- update chain center of mass
c     --- no need to recompute as we simply shifted the center of mass
      call twh_com(GLB_INCR,CRD_REAL,ichain,ddx,ddy,ddz)
c     --- update embedding density
      iboxa(1) = ibox
      call twh_reseteamrho(1,iboxa)

#if DEBUG_COMPOSITE
      write(6,*) 'end COMPOSITE'
#endif
      return
      end
