#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2012 Marcus G. Martin                       *
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
      subroutine twh_tranatom(lfailure)
c     ******************************************************************
c     * makes a translational movement randomly on a unit sphere       *
c     * the maximum displacement is controlled by rmtraa               *
c     * then number of attemped trial moves is stored in bntraatemp    *
c     * then number of successful trial moves is stored in bstraatemp  *
c     *                                                                *
c     * originally written 01-26-2001 by M.G. Martin                   *
c     * last modified 01-27-2012 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EXPON
#define FUNCTION_RANDOM
#include "functions.h"
#define DEBUG_TRANATOM 0
c     --- varibles passed to/from the subroutine
      logical lfailure
c     --- local variables
c     --- logical scalars
      logical laccept,scp_eval_derivative,lewald
c     --- logical arrays
      logical loverlap
      dimension loverlap(2)
c     --- integer scalars
      integer ichain,ibox,imolty,icbu,ioldnew,atom,units,natommove
      integer nmolty,ncmt,ipoint,nmtype,potentialstyle,ensemble
c     --- integer arrays
      integer atomlist
      dimension atomlist(1,3)
      integer globalcord
      dimension globalcord(2)
      integer iboxa
      dimension iboxa(1)
c     --- double precision scalars
      double precision dchain,rchain,burn,dvtot,dvintra,dvinter,dvext
      double precision dvelect,dvvib,dvbend,dvtors,dvsolv,dvquantum
      double precision dvclassical
      double precision xtemp,ytemp,ztemp,rxdisp,rydisp,rzdisp
      double precision scaledist,beta
      double precision ddvdl_lj,ddvdl_c,pmtest,rmtraa
c     --- double precision arrays
      double precision dvf
      dimension dvf(MAX_FOREIGN_LAMBDA)

c      write(6,*) 'start TRANATOM'
c     --- retrieve constants
      call twh_potentialstyle(GLB_GET,potentialstyle)
      call twh_lewald(GLB_GET,lewald)

c     --- select a chain at random
      rchain  = twh_random()
      call twh_nmolty(GLB_GET,nmolty)
      do icbu = 1,nmolty
         call twh_pmtamt(GLB_GET,icbu,pmtest)
         if ( rchain .lt. pmtest ) then
            imolty = icbu
            rchain = 2.0d0
         endif
      enddo

#if DEBUG_TRANATOM
      write(6,*) 'tranatom imolty:',imolty
#endif

      call twh_ensemble(GLB_GET,ensemble)
      if (ensemble .eq. ENS_muVT) then
c        --- grand canonical ensemble box 1 is the actual system
         ibox = 1
c        --- see if there are any molecules of this type in active box
         call twh_ncmt(GLB_GET,ibox,imolty,ncmt)
         if (ncmt .eq. 0) return
c        --- get a chain of type imolty from box ibox
         ipoint = idint(dble(ncmt)*twh_random())+1
         call twh_chainlist(GLB_GET,ipoint,1,imolty,ichain)
      else
         call twh_nmolectyp(GLB_GET,imolty,nmtype)
         dchain = dble(nmtype)
         ipoint = idint( dchain*twh_random() + 1 )
         call twh_parall(GLB_GET,imolty,ipoint,ichain)
      endif
      call twh_nboxi(GLB_GET,ichain,ibox)

c     --- store number of units in the molecule
      call twh_nunit(GLB_GET,imolty,units)
c     --- choose an atom at random from the molecule
      atom = idint(twh_random()*units) + 1

#if DEBUG_TRANATOM
      write(6,*) 'TRANATOM atom:',atom
#endif

c     --- move atom in a random direction
      call twh_sphere(rxdisp,rydisp,rzdisp)

c     --- compute a random displacement distance between 0 and rmtraa
      call twh_rmtraa(GLB_GET,imolty,ibox,rmtraa)
      scaledist = twh_random()*rmtraa

c     --- rescale the unit vector
      rxdisp = rxdisp*scaledist
      rydisp = rydisp*scaledist
      rzdisp = rzdisp*scaledist

c     --- move the atom by the displacement
      call twh_coordinates(GLB_GET,CRD_REAL,ichain,atom
     &     ,xtemp,ytemp,ztemp)
      xtemp = xtemp + rxdisp
      ytemp = ytemp + rydisp
      ztemp = ztemp + rzdisp
      call twh_coordinates(GLB_SET,CRD_TEMP,ichain,atom
     &     ,xtemp,ytemp,ztemp)

c     --- increment move attemp counter
      call twh_bntraa(GLB_INCR,imolty,ibox,1.0d0)

c     --- construct the atomlist
      natommove = 1
      atomlist(1,CRD_OLD) = CRD_REAL
      atomlist(1,CRD_NEW) = CRD_TEMP
      atomlist(1,3) = atom
c     --- use the real coordinates for the rest of the atoms for both
c     --- old and new
      globalcord(CRD_OLD) = CRD_REAL
      globalcord(CRD_NEW) = CRD_REAL
c     --- initialize loverlap
      loverlap(CRD_NEW) = .false.
      loverlap(CRD_OLD) = .false.

      call twh_ec_atomlist_nonbond(lfailure,ichain,ibox,natommove
     &     ,atomlist,globalcord,loverlap,dvquantum,dvclassical,dvinter
     &     ,dvintra,dvelect,dvsolv,dvext,dvf,ddvdl_lj,ddvdl_c)
      if ( lfailure ) return
      do ioldnew = CRD_NEW,CRD_OLD,(CRD_OLD-CRD_NEW)
         if ( loverlap(ioldnew) ) then
            if ( ioldnew .eq. CRD_NEW ) then
c              --- move rejected
               return
            else
               write(6,*) 'TRANTOM: Overlap in old conformation'
               lfailure = .true.
               return
            endif
         endif
      enddo

      if ( potentialstyle .eq. POT_INTERNAL ) then
c        --- initialize intramolecular bonded energies
         dvvib = 0.0d0
         dvbend = 0.0d0
         dvtors = 0.0d0

c        --- vibration energy
         call twh_ec_atomlist_bond(lfailure,ichain,ibox,natommove
     &        ,atomlist,globalcord,loverlap,dvvib)
         if ( lfailure ) return
         if ( loverlap(CRD_NEW) ) then
c           --- move rejected
            return
         elseif ( loverlap(CRD_OLD) ) then
            write(6,*) 'TRANATOM: Overlap in old conformation for bond'
            lfailure = .true.
            return
         endif

c        --- bond-bond energy
         call twh_ec_atomlist_bond_bond(ichain,ibox,natommove
     &        ,atomlist,globalcord,dvvib)

c        --- regular angle energy
         call twh_ec_atomlist_angle(lfailure,ichain,ibox,natommove
     &        ,atomlist,globalcord,loverlap,dvbend)
         if ( lfailure ) return
         if ( loverlap(CRD_NEW) ) then
c           --- move rejected
            return
         elseif ( loverlap(CRD_OLD) ) then
            write(6,*) 'TRANATOM: Overlap in old conformation for angle'
            lfailure = .true.
            return
         endif

c        --- angle-angle energy
         call twh_ec_atomlist_angle_angle(ichain,ibox,natommove
     &        ,atomlist,globalcord,dvbend)

c        --- normal torsion energy
         call twh_ec_atomlist_torsion(lfailure,ichain,ibox,natommove
     &        ,atomlist,globalcord,loverlap,dvtors)
         if ( lfailure ) return
         if ( loverlap(CRD_NEW) ) then
c           --- move rejected
            return
         elseif ( loverlap(CRD_OLD) ) then
            write(6,*) 'TRANATOM: Overlap in old conformation for'
     &           ,' regular torsions'
            lfailure = .true.
            return
         endif

c        --- improper torsion energy
         call twh_ec_atomlist_improper(ichain,ibox,natommove
     &        ,atomlist,globalcord,dvtors)

         dvclassical = dvclassical + dvvib + dvbend + dvtors
      endif

c     --- set the total energy to use for acceptance
      if ( potentialstyle .eq. POT_INTERNAL ) then
         dvtot = dvclassical
      elseif ( potentialstyle .eq. POT_EXTERNAL ) then
         dvtot = dvquantum
      else
         write(6,*) 'TRANATOM: unknown potentialstyle'
         lfailure = .true.
         return
      endif
 
c     --- test for acceptance
      laccept = .false. 

      if ( dvtot .lt. 0.0d0 ) then
c        --- accept moves if the energy goes down
         laccept = .true.
c        --- burn a random number
         burn = twh_random()
      else
         call twh_beta(GLB_GET,beta)
         if ( twh_random() .lt. twh_expon(-beta*dvtot) ) then
            laccept = .true.
         endif
      endif
c     --- else the move is rejected as the probability is so low that
c     --- we can get underflow errors
      if ( laccept ) then
#if DEBUG_TRANATOM
         write(6,*) 'TRANATOM move accepted',ichain,atom
#endif
c        --- update acceptance accumulators
         call twh_batraa(GLB_INCR,imolty,ibox,1.0d0)
c        --- update energies
         call twh_boxvclassic(GLB_INCR,ibox,dvtot)
         if ( potentialstyle .eq. POT_INTERNAL ) then
            call twh_boxvinter(GLB_INCR,ibox,dvinter)
            call twh_boxvintra(GLB_INCR,ibox,dvintra)
            call twh_boxvfield(GLB_INCR,ibox,dvext)
            call twh_boxvelect(GLB_INCR,ibox,dvelect)
            call twh_boxvbond(GLB_INCR,ibox,dvvib)
            call twh_boxvangle(GLB_INCR,ibox,dvbend)
            call twh_boxvtorsion(GLB_INCR,ibox,dvtors)
            call twh_boxvsolvation(GLB_INCR,ibox,dvsolv)
            call twh_scp_eval_derivative(GLB_GET,scp_eval_derivative)
            if (scp_eval_derivative) then
               call scp_increment_derivative(ibox, LAMBDA_C, ddvdl_c)
               call scp_increment_derivative(ibox, LAMBDA_LJ, ddvdl_lj)
            endif
            call scp_increment_foreign_energies(ibox, dvf)
         endif

c        --- assign new positions for the single atom that moved
         call twh_coordinates(GLB_GET,CRD_TEMP,ichain,atom
     &        ,xtemp,ytemp,ztemp)
         call twh_coordinates(GLB_SET,CRD_REAL,ichain,atom
     &        ,xtemp,ytemp,ztemp)
         if ( lewald ) then
c           --- update reciprocal-space sum
            call twh_recip(ibox,2)
         endif

c        --- update chain center of mass
         call twh_ctrmas(lfailure,CRD_REAL,ibox,ichain,CTR_TRANATOM)
         if ( lfailure ) return
c        --- update embedding densities
         iboxa(1) = ibox
         call twh_reseteamrho(1,iboxa)

      endif

c      write(6,*) 'end TRANATOM'
      return
      end
