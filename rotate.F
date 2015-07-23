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
      subroutine twh_rotate(lfailure)
c     ******************************************************************
c     * makes a rotational movement around "x" space-fixed axis.       *
c     * the maximum displacement is controlled by rmrot(xyz) and the   *
c     * number of successful rotation is given by bsrot(xyz).          *
c     * rotate chooses one of the three space-fixed axes at random     *
c     * and rotates the molecule around this axis.                     *
c     *                                                                *
c     * originally written prior to 1995                               *
c     * last modifed in 08-09-2011 by M.G. Martin                      *
c     ******************************************************************
      implicit none
#define FUNCTION_EXPON
#define FUNCTION_MINBOXLENGTH
#define FUNCTION_RANDOM
#include "functions.h"

#define DEBUG_ROTATE 0
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
c     --- local logical scalars
      logical laccept,lgetcom,lintra,lbonded,lpartial,lperiodic
      logical scp_eval_derivative,lewald
c     --- local logical arrays
      logical loverlap
      dimension loverlap(2)
c     --- local integer scalars
      integer ichain,ibox,totunits,imolty,icbu,iunit,nmolty,ncmt,ipoint
      integer nmtype,ensemble
c     --- local integer arrays
      integer nchainmove,chaininfo
      dimension chaininfo(1,3)
      integer iboxa
      dimension iboxa(1)
c     --- local double precision scalars
      double precision dchain,rchain,burn,rxtemp,rytemp,rztemp
      double precision dgammax,dgammay,dgammaz,rxorig,ryorig,rzorig
      double precision dvtor,dvbend,dvvib,dvquantum
      double precision dvtot,dvintra,dvinter,dvext,dvsolv,dvelect
      double precision ddvdl_lj,ddvdl_c,beta,rcmu,pmtest,rmrot
c     --- local double precision arrays
      double precision dvf
      dimension dvf(MAX_FOREIGN_LAMBDA)

c     PLEASE NOTE: changes/bugfixes made here might also need to be made
c     to the composite move (composite.F)

#if DEBUG_ROTATE
      write(6,*) 'start ROTATE'
#endif
c     --- retrieve constants
      call twh_lewald(GLB_GET,lewald)

c     --- select a chain type at random 
      rchain  = twh_random()
      call twh_nmolty(GLB_GET,nmolty)
      do icbu = 1,nmolty
         call twh_pmromt(GLB_GET,icbu,pmtest)
         if ( rchain .lt. pmtest ) then
            imolty = icbu
            rchain = 2.0d0
         endif
      enddo
#if DEBUG_ROTATE
      write(6,*) 'ROTATE imolty:',imolty
#endif
      
      call twh_ensemble(GLB_GET,ensemble)
      if (ensemble .eq. ENS_muVT) then
c        --- for grand canonical ensemble only want to perform move on
c        --- molecules actually in the 'real' system box (box = 1)
         ibox = 1         
c        --- check to see if there are any molecules of this type in
c        --- the active box
         call twh_ncmt(GLB_GET,ibox,imolty,ncmt)
         if (ncmt .eq. 0) return
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
c     --- increment move attemp counter
      call twh_bnrot(GLB_INCR,imolty,ibox,1.0d0)
#if DEBUG_ROTATE
      write(6,*) 'ROTATE ichain,ibox:',ichain,ibox
#endif
c     --- if this is a periodic molecule then it does not have a 
c     --- valid center of mass to rotate about
      call twh_lperiodic(GLB_GET,ichain,lperiodic)
      if ( lperiodic ) then
         write(6,*) 'ROTATE: attemped rotation move on a periodicly'
     &        ,' bonded molecule'
         write(6,*) 'Molecule type:',imolty
         write(6,*) 'Chain number:',ichain
         return
      endif
c     --- store number of units of ichain in totunits
      call twh_nunit(GLB_GET,imolty,totunits)

c     --- store current positions in old-new array 
      do iunit = 1,totunits
c        --- fetch the current coordinates
         call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &        ,rxtemp,rytemp,rztemp)
c        --- store them as the trial coordinates
         call twh_coordinates(GLB_SET,CRD_TEMP,ichain,iunit
     &        ,rxtemp,rytemp,rztemp)
# if DEBUG_ROTATE
         write(6,*) 'old coordinates iunit:',iunit,rxtemp,rytemp,rztemp
#endif
      enddo

c     --- initialize displacements
      dgammax = 0.0d0
      dgammay = 0.0d0
      dgammaz = 0.0d0

c     --- choose a random angular displacement over the 3 axes
      call twh_rmrot(GLB_GET,imolty,ibox,rmrot)
      dgammax = ( 2.0d0*twh_random() - 1.0d0 )* rmrot
      dgammay = ( 2.0d0*twh_random() - 1.0d0 )* rmrot
      dgammaz = ( 2.0d0*twh_random() - 1.0d0 )* rmrot

c     --- Determine the rotation coordinates
c     --- Use the center of mass for rotation
      call twh_com(GLB_GET,CRD_REAL,ichain,rxorig,ryorig,rzorig)
c     --- new com is the same as the old
      call twh_com(GLB_SET,CRD_TEMP,ichain,rxorig,ryorig,rzorig)
c     --- set the new coordinates
      call twh_rotationmatrix(CRD_TEMP,ichain,rxorig,ryorig,rzorig
     &     ,dgammax,dgammay,dgammaz)

c     --- get the energy difference
      call twh_rcmu(GLB_GET,ichain,rcmu)
      if ( rcmu .gt. twh_minboxlength(ibox)*0.25d0 ) then
c        --- molecule is large enough that it might be interacting with
c        --- iself through the periodic images (occurs for some very 
c        --- short range potentials like hard spheres).   Need to 
c        --- compute bonded interactions as these could change in this 
c        --- case
         lbonded = .true.
      else
c        --- molecule does not interact with itself through the periodic
c        --- images.  don't bother computing  bonded
         lbonded = .false.
      endif

#if DEBUG_ROTATE
      write(6,*) 'ROTATE: proposed rotation: ', 
     &      dgammax, dgammay, dgammaz
#endif
c     --- get the energy difference
      nchainmove = 1
      chaininfo(1,CRD_OLD) = CRD_REAL
      chaininfo(1,CRD_NEW) = CRD_TEMP
      chaininfo(1,3) = ichain
      lpartial = .false.
      lgetcom = .false.
      lintra = .true.
      call twh_energy_change(lfailure,ibox,nchainmove,chaininfo,loverlap
     &     ,lpartial,lgetcom,lintra,lbonded,dvtot,dvintra
     &     ,dvinter,dvext,dvelect,dvsolv,dvtor,dvbend,dvvib
     &     ,dvquantum,dvf,ddvdl_lj,ddvdl_c)
      if ( lfailure ) return
#if DEBUG_ROTATE
      write(6,*) 'ROTATE: energy_change dvelect:',dvelect
#endif
c     --- new conformation hard overlap check
      if ( loverlap(CRD_NEW) ) return
c     --- old conformation hard overlap check
      if ( loverlap(CRD_OLD) ) 
     &      stop 'disaster- overlap for old conf in ROTATE'

c     --- check for acceptance
      laccept = .false.
      if ( dvtot .lt. 0.0d0 ) then
c        --- move accepted
         laccept = .true.
c        --- burn a random number
         burn = twh_random()
      else
         call twh_beta(GLB_GET,beta)
         if ( twh_random() .lt. twh_expon(-beta*dvtot) ) then
            laccept = .true.
         endif
      endif

#if DEBUG_ROTATE
      write(6,*) 'ROTATE: dvtot',dvtot
#endif
      if ( laccept ) then
#if DEBUG_ROTATE
          write(6,*) 'ROTATE: accepted, chain',ichain
#endif
c        --- increment the successful counter array
         call twh_barot(GLB_INCR,imolty,ibox,1.0d0)
c        --- the move was accepted so update the arrays
         call twh_boxvclassic(GLB_INCR,ibox,dvtot)
         call twh_boxvinter(GLB_INCR,ibox,dvinter)
         if ( lintra ) then
            call twh_boxvintra(GLB_INCR,ibox,dvintra)
         endif
         call twh_boxvfield(GLB_INCR,ibox,dvext)
         call twh_boxvelect(GLB_INCR,ibox,dvelect)
         call twh_boxvsolvation(GLB_INCR,ibox,dvsolv)
         if ( lbonded ) then
            call twh_boxvtorsion(GLB_INCR,ibox,dvtor)
            call twh_boxvangle(GLB_INCR,ibox,dvbend)
            call twh_boxvbond(GLB_INCR,ibox,dvvib)
         endif
         call twh_scp_eval_derivative(GLB_GET,scp_eval_derivative)
         if (scp_eval_derivative) then
            call scp_increment_derivative(ibox, LAMBDA_C, ddvdl_c)
            call scp_increment_derivative(ibox, LAMBDA_LJ, ddvdl_lj)
         endif
         call scp_increment_foreign_energies(ibox, dvf)

         do iunit = 1,totunits
            call twh_coordinates(GLB_GET,CRD_TEMP,ichain,iunit
     &           ,rxtemp,rytemp,rztemp)
            call twh_coordinates(GLB_SET,CRD_REAL,ichain,iunit
     &           ,rxtemp,rytemp,rztemp)
#if DEBUG_ROTATE
            write(6,*) 'new coordinates iunit:',iunit,rxtemp,rytemp
     &           ,rztemp
#endif
         enddo

c        --- update the reciprical space sum
         if ( lewald ) call twh_recip(ibox,2)
c        --- update embedding density
         iboxa(1)=ibox
         call twh_reseteamrho(1,iboxa)
c        --- no need to update center of mass as the COM did not move
c        --- the distance from the atoms to the COM did not move
      endif

#if DEBUG_ROTATE
      write(6,*) 'end ROTATE'
#endif

      return
      end
