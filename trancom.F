#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2011 Marcus G. Martin                       *
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
      subroutine twh_trancom(lfailure)
c     ******************************************************************
c     * makes a translational movement in a random direction           *
c     * the maximum displacement is controlled by rmtrac               *
c     * number of attempts stored in bntractemp                        *
c     * number of successful trial moves is stored in bstractemp       *
c     *                                                                *
c     * originally written prior to 1995                               *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EXPON
#define FUNCTION_RANDOM
#define FUNCTION_MINBOXLENGTH
#include "functions.h"
c     --- debugging variable
#define DEBUG_TRANCOM 0
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
c     --- local logical scalars
      logical laccept,lgetcom,lintra,lbonded,lpartial
      logical scp_eval_derivative,lewald
c     --- local logical arrays
      logical loverlap
      dimension loverlap(2)
c     --- local integer scalars
      integer ichain,ibox,totunit,iunit,imolty,icbu,nchainmove,nmolty
      integer ncmt,ipoint,nmtype,potentialstyle,ensemble
c     --- local integer arrays
      integer chaininfo
      dimension chaininfo(1,3)
c     --- local double precision scalars
      double precision rx,ry,rz,dchain,ddx,ddy,ddz,rchain,burn,beta
      double precision xtemp,ytemp,ztemp,xcmi,ycmi,zcmi
      double precision dvclassical_trial,dvintra_trial,dvinter_trial
      double precision dvext_trial,dvelect_trial,dvsolv_trial
      double precision dvtor_trial,dvbend_trial,dvvib_trial
      double precision scaledist,dvquantum_trial,dv_trial
      double precision ddvdl_lj,ddvdl_c,rcmu,rmtrac,pmtest
c     --- local double precision arrays
      double precision dvf
      dimension dvf(MAX_FOREIGN_LAMBDA)
      integer iboxa
      dimension iboxa(1)

c     PLEASE NOTE: changes/bugfixes made here might also need to be made
c     to the composite move (composite.F)

#if DEBUG_TRANCOM
      write(6,*) 'start TRANCOM'
#endif
c     --- retrieve constants
      call twh_potentialstyle(GLB_GET,potentialstyle)
      call twh_lewald(GLB_GET,lewald)

c     --- select a molecule type at random
      rchain  = twh_random()
      call twh_nmolty(GLB_GET,nmolty)
      do icbu = 1,nmolty
         call twh_pmtcmt(GLB_GET,icbu,pmtest)
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
         call twh_parall(GLB_GET,imolty,ipoint,ichain)
      endif
      call twh_nboxi(GLB_GET,ichain,ibox)
      if ( ibox .eq. 0 ) then
         write(6,*) 'TRANCOM: problem selecting proper box'
         stop
      endif

#if DEBUG_TRANCOM
      write(6,*) 'TRANCOM box:',ibox
#endif

c     --- store number of units of i in totunit 
      call twh_nunit(GLB_GET,imolty,totunit)

c     --- translate on the unit sphere
      call twh_sphere(rx,ry,rz)
c     --- compute a random scale distance between 0 and rmtrac
      call twh_rmtrac(GLB_GET,imolty,ibox,rmtrac)
      scaledist = twh_random()*rmtrac
c     --- scale by the maximum displacement
      rx = rx*scaledist
      ry = ry*scaledist
      rz = rz*scaledist

c     --- increment move attempt counter
      call twh_bntrac(GLB_INCR,imolty,ibox,1.0d0)

c     --- find out if we need to displace the atoms to keep them in 
c     --- the central box
      call twh_com(GLB_GET,CRD_REAL,ichain,xcmi,ycmi,zcmi)
      ddx = xcmi + rx
      ddy = ycmi + ry
      ddz = zcmi + rz
      call twh_putcartinbox(ibox,ddx,ddy,ddz)
c     --- set the new COM to the temp setting
      call twh_com(GLB_SET,CRD_TEMP,ichain,ddx,ddy,ddz)
c     --- if the coordinates were shifted then compute the shift
      call twh_com(GLB_GET,CRD_REAL,ichain,xcmi,ycmi,zcmi)
      ddx = ddx - xcmi
      ddy = ddy - ycmi
      ddz = ddz - zcmi

      do iunit = 1,totunit
         call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit,rx,ry,rz)
         rx = rx + ddx
         ry = ry + ddy
         rz = rz + ddz
         call twh_coordinates(GLB_SET,CRD_TEMP,ichain,iunit,rx,ry,rz)
      enddo

#if DEBUG_TRANCOM
      write(6,*) 'TRANCOM: ichain:',ichain
      write(6,*) 'TRANCOM: totunit:',totunit
      do iunit = 1,totunit
         write(6,*) 'TRANCOM: iunit:',iunit
         write(6,*) 'TRANCOM: dd*:',ddx,ddy,ddz
         call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit,rx,ry,rz)
         write(6,*) 'TRANCOM: old coords:',rx,ry,rz
         call twh_coordinates(GLB_GET,CRD_TEMP,ichain,iunit,rx,ry,rz)
         write(6,*) 'TRANCOM: new coords:',rx,ry,rz
      enddo
#endif

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

#if DEBUG_TRANCOM
      write(6,*) 'TRANCOM: proposed translation: ', ddx, ddy, ddz
#endif
      nchainmove = 1
      chaininfo(1,CRD_OLD) = CRD_REAL
      chaininfo(1,CRD_NEW) = CRD_TEMP
      chaininfo(1,3) = ichain
      lpartial = .false.
      lgetcom = .false.
      call twh_energy_change(lfailure,ibox,nchainmove,chaininfo,loverlap
     &     ,lpartial,lgetcom,lintra,lbonded,dvclassical_trial
     &     ,dvintra_trial,dvinter_trial,dvext_trial,dvelect_trial
     &     ,dvsolv_trial,dvtor_trial,dvbend_trial
     &     ,dvvib_trial,dvquantum_trial,dvf,ddvdl_lj,ddvdl_c)
      if ( lfailure ) return
#if DEBUG_TRANCOM
      write(6,*) 'TRANCOM: dvinter_trial:',dvinter_trial
#endif
c     --- check new hard overlap, if overlap move rejected
      if ( loverlap(2) ) return
c     --- check old hard overlap, if overlap then big problems
      if (loverlap(1)) stop 'TRANCOM: disaster loverlap in old conf'

      if ( potentialstyle .eq. POT_INTERNAL ) then
         dv_trial = dvclassical_trial
      elseif ( potentialstyle .eq. POT_EXTERNAL ) then
         dv_trial = dvquantum_trial
      else
         write(6,*) 'TRANCOM: unknown potentialstyle'
         write(6,*) potentialstyle
         stop
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

#if DEBUG_TRANCOM
      write(6,*) 'TRANCOM: dv_trial ',dv_trial
      write(6,*) 'TRANCOM: laccept:',laccept
#endif
      if ( .not. laccept ) return
#if DEBUG_TRANCOM
      write(6,*) 'TRANCOM: accepted chain',ichain
#endif
c     --- update acceptance accumulators
      call twh_batrac(GLB_INCR,imolty,ibox,1.0d0)
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
         write(6,*) 'TRANCOM: unknown potentialstyle'
         write(6,*) potentialstyle
         stop
      endif
c     --- store new coordinates
      do iunit = 1,totunit
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
      iboxa(1)=ibox
      call twh_reseteamrho(1,iboxa)

#if DEBUG_TRANCOM
      write(6,*) 'end TRANCOM'
#endif
      return
      end
