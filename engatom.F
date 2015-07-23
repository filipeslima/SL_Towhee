#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2009 Marcus G. Martin                       *
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
      subroutine twh_engatom(lfailure,testcord,globalcord,engastyle
     &     ,lovrlap,ichain,imolty,ibox,ichoi,iufrom,ntogrow,maxlen
     &     ,lnbterm,dimatoms,dimtrials,rxp,ryp,rzp,vctry,vtrext
     &     ,vtrintra,vtrinter,vtrelect,vtrewald,vqtry)
c     ******************************************************************
c     * calculates the potential energy and the boltzmann factor for   *
c     * atoms that are being inserted in a configurational-bias move   *
c     *                                                                *
c     * testcord: controls which set of coordinates are used for the   *
c     *            ichain                                              *
c     * globalcord: controls which set of coordinates are used for all *
c     *             other chains in the system                         *
c     * engastyle: 1: single atom inserted that is not bonded to any   *
c     *                other existing atoms                            *
c     *            2: one or more atoms inserted that are all bonded to*
c     *               iufrom with a maximum distance of maxlen         *
c     *            3: inserting all of the atoms in an entire molecule *
c     * lovrlap is true if this subroutine runs into a walk termination*
c     * ichain is the molecule number - used if this is the old conf   *
c     * imolty is the molecule type                                    *
c     * ibox is the box to use when computing intermolecular energies  *
c     * ichoi is the number of trial positions for each of the atoms   *
c     * iufrom the the atom we are growing the new atoms from          *
c     * ntogrow is the number of atoms we are growing this step        *
c     * glist is the list of the atoms we are growing this step        *
c     * maxlen is the maximum length of any of the atoms from iufrom   *
c     * lnbterm is a logical array that is set by the calling routine  *
c     * to indicate whether this trial has already been terminated.    *
c     * it also gets returns with the overlap information in case a    *
c     * trial results in an overlap in this routine                    *
c     * this trial should be computed                                  *
c     * dimatoms is the first dimension in the r*p arrays              *
c     * dimtrials is the second dimension in the r*p arrays            *
c     * rxp(dimatoms,dimtrials), ryp and rzp are the coordinates for   *
c     * each atom and trial                                            *
c     *                                                                *
c     * rewritten from old boltz and boltzu 09-23-1997 by M.G. Martin  *
c     * last modified 06-25-2009 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer dimatoms,dimtrials
      logical lfailure,lovrlap,lnbterm
      dimension lnbterm(dimtrials)
      integer ichain,imolty,ibox,ichoi,iufrom,ntogrow,engastyle
      integer testcord,globalcord
      double precision maxlen
      double precision rxp,ryp,rzp
      dimension rxp(dimatoms,dimtrials)
      dimension ryp(dimatoms,dimtrials)
      dimension rzp(dimatoms,dimtrials)
      double precision vctry
      dimension vctry(dimtrials)
      double precision vtrext
      dimension vtrext(dimtrials)
      double precision vtrintra
      dimension vtrintra(dimtrials)
      double precision vtrinter
      dimension vtrinter(dimtrials)
      double precision vtrelect
      dimension vtrelect(dimtrials)
      double precision vtrewald
      dimension vtrewald(dimtrials)
      double precision vqtry
      dimension vqtry(dimtrials)
c     --- local variables
      integer nchain,potentialstyle

c     --- retrieve constants
      call twh_potentialstyle(GLB_GET,potentialstyle)
      if ( potentialstyle .eq. POT_INTERNAL ) then
         call twh_nchain(GLB_GET,nchain)
         call twh_classical_engatom(lfailure,testcord,globalcord
     &        ,engastyle,lovrlap,nchain,ichain,imolty,ibox,ichoi,iufrom
     &        ,ntogrow,maxlen,lnbterm,dimatoms,dimtrials,rxp,ryp,rzp
     &        ,vctry,vtrext,vtrintra,vtrinter,vtrelect,vtrewald)
         if ( lfailure ) return
      elseif ( potentialstyle .eq. POT_EXTERNAL ) then
c        --- do external stuff
         call twh_external_engatom(lfailure,testcord,ichain,ibox,ichoi
     &        ,ntogrow,dimtrials,vqtry)
         if ( lfailure ) return
      else
         write(6,*) 'ENGATOM: unknown potentialstyle'
         write(6,*) potentialstyle
         lfailure = .true.
         return
      endif
      return
      end

      subroutine twh_external_engatom(lfailure,testcord,ichain,ibox
     &     ,ichoi,ntogrow,dimtrials,vqtry)
c     ******************************************************************
c     * computes the external energy for each trial in a               *
c     * configurational bias style of move                             *
c     * this does not appear to be currently functional so beware if   *
c     * you are using this routine                                     *
c     *                                                                *
c     * originally written 10-24-2005 by M.G. Martin                   *
c     * last modified 12-02-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer testcord,ichain,ibox,ichoi,ntogrow,dimtrials
      double precision vqtry
      dimension vqtry(dimtrials)
#define DEBUG_EXTERNAL_ENGATOM 0
c     --- local variables
      integer itrial
      double precision vquantum

      do itrial = 1,ichoi
         write(6,*) 'engatom ibox,ichain,testcord,ntogrow,itrial'
     &      ,ibox,ichain,testcord,ntogrow,itrial
c     --- not currently functional, need to rewrite call
c         call twh_quantumenergy(lfailure,ibox,ichain,testcord,ntogrow
c     &        ,itrial,vquantum)
         vquantum = 0.0d0
         if ( lfailure ) return
         vqtry(itrial) = vquantum
      enddo

      return
      end


      subroutine twh_classical_engatom(lfailure,testcord,globalcord
     &     ,engastyle,lovrlap,nchain,ichain,imolty,ibox,ichoi,iufrom
     &     ,ntogrow,maxlen,lnbterm,dimatoms,dimtrials,rxp,ryp,rzp,vctry
     &     ,vtrext,vtrintra,vtrinter,vtrelect,vtrewald)
c     ******************************************************************
c     * computes the classical energy on an atom by atom basis for the *
c     * configurational-bias class of moves                            *
c     *                                                                *
c     * split from engatom 10-24-2005 by M.G. Martin                   *
c     * last modified 06-25-2009 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LEXIST
#define FUNCTION_EWALD_CORRECT
#define FUNCTION_EWALD_SELF
#define FUNCTION_LINCLUDE
#define FUNCTION_ONEFIVETYPE
#define FUNCTION_VCOULOMB
#define FUNCTION_VFIELD
#define FUNCTION_VONEFIVE
#define FUNCTION_VTWOBODY
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer dimatoms,dimtrials
      logical lfailure,lovrlap,lovrtest
      logical lnbterm
      dimension lnbterm(dimtrials)
      integer nchain,ichain,imolty,ibox,ichoi,iufrom,ntogrow,engastyle
      integer testcord,globalcord
      double precision maxlen
      double precision rxp,ryp,rzp
      dimension rxp(dimatoms,dimtrials)
      dimension ryp(dimatoms,dimtrials)
      dimension rzp(dimatoms,dimtrials)
      double precision vctry
      dimension vctry(dimtrials)
      double precision vtrext
      dimension vtrext(dimtrials)
      double precision vtrintra
      dimension vtrintra(dimtrials)
      double precision vtrinter
      dimension vtrinter(dimtrials)
      double precision vtrelect
      dimension vtrelect(dimtrials)
      double precision vtrewald
      dimension vtrewald(dimtrials)
c     --- local variables
c     --- logical scalars
      logical lcompute,lusecom,lcoulomb,lewald,ldualcutoff
c     --- logical arrays
      logical lcmno
      dimension lcmno(nchain)
c     --- integer scalars
      integer count,iunit,cutstyle,jlist,ioftype,nmolty,nfield
      integer itrial,ntii,jchain,jj,ntjj
      integer cnt,iu,jmolty,ii,numjunits,ncmt
c     --- double precision scalars
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij
      double precision rcm,rcmsq,rijsq,rij,qqfact,xcmj,ycmj,zcmj,rcmj
      double precision vinter,vintra,vext,velect,vewald,dvtemp,rcelect
      double precision rcutin,rcutinsq
c     --- double precision arrays
      double precision cord
      dimension cord(3)
      double precision vf
      dimension vf(MAX_FOREIGN_LAMBDA)

c      write(6,*) 'start ENGATOM',ntogrow

c     --- retrieve constants
      call twh_constant_qqfact(GLB_GET,qqfact)
      call twh_nmolty(GLB_GET,nmolty)
      call twh_nunit(GLB_GET,imolty,iunit)
      call twh_lusecom(GLB_GET,lusecom)
      if ( ibox .ne. 0 ) call twh_rcelect(GLB_GET,ibox,rcelect)
      call twh_lcoulomb(GLB_GET,lcoulomb)
      call twh_lewald(GLB_GET,lewald)
      call twh_rcutin(GLB_GET,rcutin)
      call twh_rcutinsq(GLB_GET,rcutinsq)
      call twh_ldualcutoff(GLB_GET,ldualcutoff)

c     --- set up lcmno array for intermolecular interactions
      if ( engastyle .eq. 2 .and. ibox .ne. 0 ) then
c        --- check previous bead (iufrom) COM for each molecule in the 
c        --- box
         call twh_coordinates(GLB_GET,testcord,ichain,iufrom
     &        ,rxui,ryui,rzui)

         do jmolty = 1,nmolty
            call twh_ncmt(GLB_GET,ibox,jmolty,ncmt)
            do jlist = 1,ncmt
               call twh_chainlist(GLB_GET,jlist,ibox,jmolty,jchain)
               lcmno(jchain) = .false.
               if ( ichain .ne. jchain ) then
c                 --- minimum image the pseudo-ctrmas pair separation
                  call twh_com(GLB_GET,CRD_REAL,jchain,xcmj,ycmj,zcmj)
                  call twh_mimage( ibox,rxui,ryui,rzui,xcmj,ycmj,zcmj
     &                 ,rxuij,ryuij,rzuij )
                  rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                  call twh_rcmu(GLB_GET,jchain,rcmj)
                  if ( ldualcutoff ) then
                     rcm = rcutin + rcmj + maxlen
                  else
c                    --- when not using dual cut-off start with vdw
                     call twh_rcut(GLB_GET,rcm)
                     if ( lcoulomb ) then
c                       --- switch to rcelect if it is bigger than rcm
                        if ( rcelect .gt. rcm ) rcm = rcelect
                     endif
c                    --- addin the com radius and max bond length
                     rcm = rcm + rcmj + maxlen
                  endif
                  rcmsq = rcm*rcm
                  if (rijsq .gt. rcmsq ) lcmno(jchain) = .true.
               endif
            enddo
         enddo
      endif

c     --- begin energy computations
      do itrial = 1, ichoi
c        --- initialze the variables
         lovrtest = .false.
         vinter = 0.0d0
         vintra = 0.0d0
         vext = 0.0d0
         velect = 0.0d0
         vewald = 0.0d0 
         if ( lnbterm(itrial) ) then
c           --- this trial was already determined to be an overlap
c           --- in the subroutine that called us
            lovrtest = .true.
            goto 19
         endif
         do count = 1,ntogrow
            call twh_glist(GLB_GET,count,ii)
            
            if ( lewald .and. engastyle .ne. 3 .and. ibox .ne. 0 ) then
c           -- This part does not change for fixed charge moves, but is
c           -- used in the cbswap rosenbluth weight. - ewald self term
               vewald = vewald + twh_ewald_self(ibox,imolty,ii)
            endif
         enddo

c        --- no intramolecular interactions if this is the first bead
         if ( engastyle .eq. 1 .or. engastyle .eq. 3 ) goto 18

C        *****************************************
C        *** INTRACHAIN BEAD-BEAD INTERACTIONS***
C        *****************************************
c        --- cycle through molecule and check bead by bead
         do 17 iu = 1, iunit
c           --- see if iu exists in the new chain yet
            if (.not. twh_get_lexist(iu)) goto 17
            
c           --- loop over all the grown beads
            do count = 1,ntogrow
               call twh_glist(GLB_GET,count,ii)

c              --- see if iu has nonbonded intramolecular interaction 
c              --- with ii must compute the distances if linclude is 
c              -- true or lcoulomb
               if ( lcoulomb .or. twh_linclude(imolty,ii,iu) ) then
c                 --- assign bead type for ii,iu, and the cross term
                  call twh_ntype(GLB_GET,imolty,ii,ntii)
                  call twh_ntype(GLB_GET,imolty,iu,ntjj)
c                 --- determine distances
                  call twh_coordinates(GLB_GET,testcord,ichain,iu
     &                 ,rxui,ryui,rzui)
                  call twh_mimage ( ibox,rxui,ryui,rzui
     &                 ,rxp(count,itrial),ryp(count,itrial)
     &                 ,rzp(count,itrial),rxuij,ryuij,rzuij )
                  rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
               endif

               if ( twh_linclude(imolty,ii,iu) ) then
c                 --- non-coulombic nonbond
                  cutstyle = 0
                  ioftype = twh_onefivetype(imolty,ii,iu)
                  if ( ioftype .eq. 0 ) then
c                    --- regular nonbond
                     vintra = vintra + twh_vtwobody(lfailure,rijsq,ntii
     &                    ,ntjj,cutstyle,lovrtest,vf,dvtemp)
                     if ( lfailure ) return
                  else
c                    --- special one-five nonbond
                     vintra = vintra + twh_vonefive(lfailure,rijsq
     &                    ,ioftype,cutstyle,lovrtest)
                     if ( lfailure ) return
                  endif
                  if ( lovrtest ) goto 19

                  if ( lcoulomb ) then
c                    --- intramolecular electrostatics
                     rij = dsqrt(rijsq)
                     velect = velect + twh_vcoulomb(ibox,rij,imolty,ii
     &                    ,imolty,iu,.false.,vf,dvtemp,lovrtest)
                     if ( lovrtest ) goto 19
                  endif
               elseif ( lewald .and. ibox .ne. 0 ) then
c                 --- ewald sum correction term (not needed if linclude)
c                 --- note that this term is applied under the 
c                 --- assumption that these molecule have no non-bond 
c                 --- coulombic interactions in the main box.  Some of 
c                 --- these atoms actually do have a non-bond coulombic
c                 --- 1-4 interaction which is accounted for in 
c                 --- vtorsion.  The combination of vtorsion with this 
c                 --- correction term results in the proper answer
                  rij = dsqrt(rijsq)
                  vewald = vewald
     &                 + twh_ewald_correct(ibox,imolty,ii,iu,rij)
               endif
            enddo
 17      continue

         if ( lewald .and. ibox .ne. 0 ) then
c          --- ewald sum correction term for interactions of the 
c          --- growing beads with each other
            do cnt = 1,ntogrow-1
               call twh_glist(GLB_GET,cnt,iu)
               do count = cnt+1,ntogrow
                  call twh_glist(GLB_GET,count,ii)
c                 --- assign bead type for ii,iu
                  call twh_ntype(GLB_GET,imolty,ii,ntii)
                  call twh_ntype(GLB_GET,imolty,iu,ntjj)
c                 --- determine distances - use trial chain coordinates
                  call twh_mimage ( ibox,rxp(cnt,itrial),ryp(cnt,itrial)
     &                 ,rzp(cnt,itrial),rxp(count,itrial)
     &                 ,ryp(count,itrial),rzp(count,itrial)
     &                 ,rxuij,ryuij,rzuij )
                  rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij

c                 --- ewald sum correction term
                  rij = dsqrt(rijsq)
                  vewald = vewald
     &                 + twh_ewald_correct(ibox,imolty,ii,iu,rij)
               enddo
            enddo
         endif

c     --- if no intramolecular interactions then come here to 18
 18      continue

C     **************************************************************
C     *** CALCULATION OF INTERACTION ENERGY WITH EXTERNAL FIELD ***
C     **************************************************************
 
         call twh_nfield(GLB_GET,nfield)
         if ( nfield .ne. 0 ) then
            do count = 1,ntogrow
               call twh_glist(GLB_GET,count,ii)
               call twh_ntype(GLB_GET,imolty,ii,ntii)
               cord(1) = rxp(count,itrial)
               cord(2) = ryp(count,itrial)
               cord(3) = rzp(count,itrial)
               vext = vext + twh_vfield(lfailure,lovrtest,ibox,ntii
     &              ,imolty,ichain,ii,cord)
               if ( lfailure ) return
               if ( lovrtest ) goto 19
            enddo
         endif

c        *******************************
c        *** INTERCHAIN INTERACTIONS ***
c        *******************************
c        --- grand canonical ensemble: do need to do interchain energy 
c        --- calculations for box 0 as that is the ideal gas reservoir
         if ( ibox .eq. 0 ) goto 19

c        ---  loop over all chains except ichain
         do jmolty = 1,nmolty
            call twh_nunit(GLB_GET,jmolty,numjunits)
            call twh_ncmt(GLB_GET,ibox,jmolty,ncmt)
            do 98 jlist = 1,ncmt
               call twh_chainlist(GLB_GET,jlist,ibox,jmolty,jchain)
c              --- see if we compute this interation
c              --- avoid interaction with self
               if ( ichain .eq. jchain ) goto 98
               if ( engastyle .eq. 2 ) then
c                 --- check COM table calculated above
                  if (lcmno(jchain)) goto 98
               endif
c              --- loop over all beads of molecule i grown this step
               do count = 1,ntogrow
c                 --- assign bead type for ii
                  call twh_glist(GLB_GET,count,ii)
                  call twh_ntype(GLB_GET,imolty,ii,ntii)
c                 --- assign positions to r*ui
                  rxui = rxp(count,itrial)
                  ryui = ryp(count,itrial)
                  rzui = rzp(count,itrial)                     
                  if ( lusecom .and. engastyle .eq. 1 ) then
c                      --- check if ctrmas within rcmsq
c                      --- minimum image the ctrmas pair separations
                     call twh_com(GLB_GET,CRD_REAL,jchain
     &                    ,xcmj,ycmj,zcmj)
                     call twh_mimage ( ibox,rxui,ryui,rzui,xcmj
     &                    ,ycmj,zcmj,rxuij,ryuij,rzuij)
                     rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
c                    --- determine cutoff
                     call twh_rcmu(GLB_GET,jchain,rcmj)
                     if ( ldualcutoff ) then
c                       --- must be engastyle 1, so no previous bead
                        rcm = rcutin + rcmj
                     else
                        call twh_rcut(GLB_GET,rcm)
                        if ( lcoulomb ) then
c                          --- check the coulombic cutoff
                           if ( rcelect .gt. rcm ) rcm = rcelect
                        endif
c                       --- add in the com radius
                        rcm = rcm + rcmj
                     endif
                     rcmsq = rcm*rcm
c                    --- check if interaction distance is greater than
c                    --- cutoff
                     if ( rijsq .gt. rcmsq ) goto 98
                  endif

c                 --- loop over all beads jj of chain j
                  do jj = 1,numjunits
                     call twh_ntype(GLB_GET,jmolty,jj,ntjj)
                     call twh_coordinates(GLB_GET,globalcord,jchain,jj
     &                    ,rxuj,ryuj,rzuj)
c                    --- minimum image the pair separations
                     call twh_mimage ( ibox,rxui,ryui,rzui
     &                    ,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij )
                     rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
c                    --- compute vinter (ex. lennard-jones)
                     cutstyle = 1
                     vinter = vinter + twh_vtwobody(lfailure,rijsq,ntii
     &                    ,ntjj,cutstyle,lovrtest,vf,dvtemp)
                     if ( lfailure ) return
                     if ( lovrtest ) goto 19

c                    --- compute velect (coulomb and ewald)
                     if ( lcoulomb ) then
                        if ( ldualcutoff ) then
                           if ( rijsq .lt. rcutinsq ) then
c                             --- compute this here
                              lcompute = .true.
                           else
c                             --- will compute this later
                              lcompute = .false.
                           endif
                        else
c                          --- always compute if not dual cutoff
                           lcompute = .true.
                        endif
                        if ( lcompute ) then
                           rij = dsqrt(rijsq)
                           velect = velect + twh_vcoulomb(ibox,rij
     &                          ,imolty,ii,jmolty,jj,.false.,vf,dvtemp
     &                          ,lovrtest)
                           if ( lovrtest ) goto 19
                        endif
                     endif
                  enddo
               enddo
 98         continue
         enddo

C        *********************************************
C        ***  CALCULATION OF TOTAL POTENTIAL ENERGY***
C        *********************************************
 
 19      lnbterm(itrial) = lovrtest
         if ( .not. lnbterm(itrial) ) then
c           --- add in the appropriate factor for coulombic intractions
            velect = velect*qqfact
            vewald = vewald*qqfact
            vctry(itrial) = vinter + vintra + vext + velect + vewald
            vtrintra(itrial) = vintra
            vtrext(itrial)   = vext
            vtrinter(itrial) = vinter
            vtrelect(itrial) = velect
            vtrewald(itrial) = vewald
         endif

      enddo

      lovrlap = .true.
      do itrial = 1, ichoi
         if ( .not. lnbterm(itrial)) lovrlap = .false.
      enddo

c      write(6,*) 'end ENGATOM'

      return

      end
