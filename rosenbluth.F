#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2010 Marcus G. Martin                       *
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
      subroutine twh_rosenbluth(lfailure,lnew,lterm,testcord,globalcord
     &     ,ichain,imolty,isteps,ibox,igrow,dpscale,wrosen,vrostot
     &     ,vrosvib,vrosang,vrostor,vrosext,vrosintra,vrosinter
     &     ,vroselect,vrosewald)
c     ******************************************************************
c     * a shell around the main rosenbluth routine that fetches some   *
c     * of the trial values for use dimensioning arrays in the main    *
c     * routine                                                        *
c     *                                                                *
c     * written in current form 12-06-2008 by M.G. Martin              *
c     * last modified 11-16-2010 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to and from the subroutine
      logical lnew,lterm,lfailure
      integer ichain,imolty,isteps,ibox,igrow,testcord,globalcord
      double precision dpscale,wrosen
      double precision vrostot,vrosvib,vrosang,vrostor,vrosext
      double precision vrosintra,vrosinter,vroselect,vrosewald
c     --- local variables
      integer nch_nb,nch_pre_nb,nch_tor,nch_tor_connect
      integer nch_bend_a,nch_bend_b,max_nch_bend
      integer nch_vib

c     --- retrieve appropriate trial constants for this moltyp
      call twh_nch_nb(GLB_GET,imolty,nch_nb)
      call twh_nch_pre_nb(GLB_GET,imolty,nch_pre_nb)
      call twh_nch_tor(GLB_GET,imolty,nch_tor)
      call twh_nch_tor_connect(GLB_GET,imolty,nch_tor_connect)
      call twh_nch_bend_a(GLB_GET,imolty,nch_bend_a)
      call twh_nch_bend_b(GLB_GET,imolty,nch_bend_b)
      max_nch_bend = nch_bend_a
      if ( nch_bend_b .gt. max_nch_bend ) then
         max_nch_bend = nch_bend_b
      endif
      call twh_nch_vib(GLB_GET,imolty,nch_vib)
      call twh_rosenbluth_main (lfailure,lnew,lterm,testcord
     &     ,globalcord,ichain,imolty,isteps,ibox,igrow
     &     ,nch_nb,nch_pre_nb,nch_tor,nch_tor_connect
     &     ,nch_bend_a,nch_bend_b,max_nch_bend,nch_vib
     &     ,dpscale,wrosen
     &     ,vrostot,vrosvib,vrosang,vrostor,vrosext,vrosintra,vrosinter
     &     ,vroselect,vrosewald)
      
      return
      end

      subroutine twh_rosenbluth_main (lfailure,lnew,lterm,testcord
     &     ,globalcord,ichain,imolty,isteps,ibox,igrow
     &     ,nch_nb,nch_pre_nb,nch_tor,nch_tor_connect
     &     ,nch_bend_a,nch_bend_b,max_nch_bend,nch_vib
     &     ,dpscale,wrosen
     &     ,vrostot,vrosvib,vrosang,vrostor,vrosext,vrosintra,vrosinter
     &     ,vroselect,vrosewald)
c     ******************************************************************
c     * performs a configurational bias move for linear or branched    *
c     * molecules                                                      *
c     *                                                                *
c     * lnew: true if this is a growth generating all new coordinates  *
c     *       false if we should trace along the old coordinates       *
c     * lterm: returned as true if the growth terminates               *
c     * testcord: place the coordinates are stored for the growth chain*
c     * globalcord: coordinates to use for all other chains            *
c     * ichain: the number of the molecule to use for coordinates      *
c     * imolty: the molecule type                                      *
c     * isteps: the number of growth steps                             *
c     * ibox: the box the molecule is being grown into                 *
c     * igrow: the number of units in the molecule                     *
c     * dpscale: number of times the wieght was scaled by scalelog     *
c     *                                                                *
c     * originally rewritten from rosnewbr and rosoldbr in 1998        *
c     *            by M.G. Martin                                      *
c     * renamed rosenbluth_main 12-06-2008 by M.G. Martin              *
c     * last modified 11-06-2012 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_RANDOM
#include "functions.h"
#define DEBUG_ROSENBLUTH 0
c     --- variables passed to and from the subroutine
      logical lnew,lterm,lfailure
      integer ichain,imolty,isteps,ibox,igrow,testcord,globalcord
      integer nch_nb,nch_pre_nb,nch_bend_a,nch_bend_b,max_nch_bend
      integer nch_tor,nch_tor_connect,nch_vib
      double precision dpscale,wrosen
      double precision vrostot,vrosvib,vrosang,vrostor,vrosext
      double precision vrosintra,vrosinter,vroselect,vrosewald
c     --- local variables
c     --- logical scalars
      logical lovrlap,lfound,lsoftovr,ldihedterm,luseoldone,lovrangle
      logical ldihed_direct,ldihed_indirect,ltemp,ldelta_dist
c     --- logical arrays
      logical lovrprenb
      dimension lovrprenb(nch_pre_nb)
      logical lnbterm
      dimension lnbterm(nch_nb)
c     --- integer scalars
      integer engastyle,iouter,iuprev,iufrom,ntogrow,count,torkey
      integer iu,iw,ip,ichtor,angstart,iwalk,istoreflag,potentialstyle
      integer grownbon,nxangle,nxbias,cbmc_form_code
c     --- double precision scalars
      double precision alpha,gamma,xub,yub,zub,xxx,yyy,zzz,length,dpmin
      double precision wmin,gcba_vvib,gcba_vbend,gcba_vtor
      double precision maxlen,rbf,bs,vvibtr,w_nb
      double precision vcbvib,vcbbend,vcbtor,bfcbbias
      double precision w_cbbond,gcba_weight,w_cbtor,phishift
      double precision dpscale_bond,dpscale_angle,dpscale_dihed
      double precision dpscale_nb,scalelolog
      double precision rxfrom,ryfrom,rzfrom,rxprev,ryprev,rzprev
      double precision xtemp,ytemp,ztemp
c     --- double precision arrays
      double precision  dpscale_adjust,dpscale_tor
      dimension dpscale_adjust(nch_nb),dpscale_tor(nch_nb)
      double precision dpscale_prenb
      dimension dpscale_prenb(nch_pre_nb)
      double precision xtorcord,ytorcord,ztorcord
      dimension xtorcord(NNBOND),ytorcord(NNBOND),ztorcord(NNBOND)
      double precision bfadjust,vmin
      dimension bfadjust(nch_nb)
      double precision bfactor
      dimension bfactor(nch_nb)
      double precision w_tor,boltzbias
      dimension w_tor(nch_nb),boltzbias(nch_nb)
      double precision vxvibtr,vxbendtr,vtgtr
      dimension vxvibtr(nch_nb),vxbendtr(nch_nb)
      dimension vtgtr(nch_nb)
      double precision w_tout,vvib_tout,vbend_tout,vtor_tout,bias_tout
      dimension w_tout(nch_pre_nb),vvib_tout(nch_pre_nb)
      dimension vbend_tout(nch_pre_nb),vtor_tout(nch_pre_nb)
      dimension bias_tout(nch_pre_nb)
      double precision rx_tout,ry_tout,rz_tout
      dimension rx_tout(nch_pre_nb,NNBOND)
      dimension ry_tout(nch_pre_nb,NNBOND)
      dimension rz_tout(nch_pre_nb,NNBOND)
      double precision rxp,ryp,rzp
      dimension rxp(NNBOND,nch_nb)
      dimension ryp(NNBOND,nch_nb)
      dimension rzp(NNBOND,nch_nb)
      double precision vctry
      dimension vctry(nch_nb)
      double precision vtrext
      dimension vtrext(nch_nb)
      double precision vtrintra
      dimension vtrintra(nch_nb)
      double precision vtrinter
      dimension vtrinter(nch_nb)
      double precision vtrelect
      dimension vtrelect(nch_nb)
      double precision vtrewald
      dimension vtrewald(nch_nb)
      double precision vqtry
      dimension vqtry(nch_nb)

#if DEBUG_ROSENBLUTH
      write(6,*) 'start ROSENBLUTH',lnew,testcord
#endif
c     --- retrieve constants
      call twh_potentialstyle(GLB_GET,potentialstyle)
      call twh_scalelolog(GLB_GET,scalelolog)
      call twh_cbmc_form_code(GLB_GET,cbmc_form_code)

      alpha = 0.0d0
      gamma = 0.0d0
      lterm = .false.
      ldihedterm = .false.
c     --- rosenbluth always uses engastyle 2
      engastyle = 2
c     --- return without doing anything if this is the ideal gas 
c     --- resivour
      if ( ibox .eq. 0 ) return

c    *******************************************
c    * Rosenbluth weight of trial conformation*
c    *******************************************
c     --- initialize dpscale
      dpscale = 0.0d0
c     --- set lexist to lexsched
      do iu = 1,igrow
         call twh_lexsched(GLB_GET,iu,ltemp)
         call twh_lexist(GLB_SET,iu,ltemp)
      enddo

c     ******************************************************************
c     * perform the biased selection for the entire walk               *
c     ******************************************************************
      do iw = 1, isteps
         call twh_growfrom(GLB_GET,iw,iufrom)
         call twh_grownum(GLB_GET,iw,ntogrow)
         call twh_growprev(GLB_GET,iw,iuprev)
c        --- set up the cone if we have a previous bead
         if ( iuprev .ne. 0 ) then
c           --- if we have a previous bead that is already in place
c           --- then we can set up the iuprev, iufrom cone just once
c           --- and use it throughout.  If not, then we will have to 
c           --- set up this cone every time as we are choosing a place
c           --- for the iuprev vector
c           --- set up the cone based on iuprev and iufrom
            call twh_coordinates(GLB_GET,testcord,ichain,iufrom
     &           ,rxfrom,ryfrom,rzfrom)
            call twh_coordinates(GLB_GET,testcord,ichain,iuprev
     &           ,rxprev,ryprev,rzprev)
            call twh_mimage(ibox,rxfrom,ryfrom,rzfrom,rxprev,ryprev
     &           ,rzprev,xxx,yyy,zzz)
c           --- compute vector length
            length = twh_distance(xxx,yyy,zzz)
c           --- turn vectors into unit vector
            xxx = xxx / length 
            yyy = yyy / length
            zzz = zzz / length
            call twh_cone(lfailure,1,xxx,yyy,zzz,alpha,gamma
     &           ,xub,yub,zub)
            if ( lfailure ) return
         endif

c        --- initialize maxlen
         maxlen = 0.0d0

         if ( cbmc_form_code .eq. CBMC_FORM_MS1999_MT2004 ) then
c           --- bond and angle are decoupled from the rest
c           --- only need to select them once and then use in the 
c           --- dihedral selections down below

            if ( .not. lnew ) then
c              --- want to use the old conformation for the first 
c              --- selection
               luseoldone = .true.
            else
c              --- all selections generated from scratch
               luseoldone = .false.
            endif
c           --- perform the biased selection to get lengths
            call twh_getcbbond(lfailure,luseoldone,testcord,iw,ichain
     &           ,imolty,iuprev,nch_vib,vvibtr,maxlen,w_cbbond
     &           ,dpscale_bond,ibox)
            if ( lfailure ) return
#if DEBUG_ROSENBLUTH
            write(6,*) 'w_cbbond:',w_cbbond
#endif
c           --- perform the biased selection to get bond angles
            call twh_getcbangle(lfailure,luseoldone,lovrangle,testcord
     &           ,iw,ichain
     &           ,imolty,angstart,iuprev,nch_bend_a,nch_bend_b
     &           ,max_nch_bend
     &           ,gcba_vvib,gcba_vbend,gcba_vtor,gcba_weight
     &           ,dpscale_angle,ibox)
            if ( lfailure ) return
            if ( lovrangle ) then
               lterm = .true.
               if ( .not. lnew ) then
                  write(6,*) 'ROSENBLUTH: angle selction problem'
                  write(6,*) 'hard overlap in old angle selection'
                  lfailure = .true.
                  return
               endif
               return
            endif
#if DEBUG_ROSENBLUTH
            write(6,*) 'gcba_weight:',gcba_weight
#endif
c           --- setup all of the information needed to do the dihedral
c           --- selection
            call twh_setcbdihed(lfailure,ldihed_direct,ldihed_indirect
     &           ,ldelta_dist,ntogrow,iw,ichain,imolty,iuprev,iufrom
     &           ,nch_tor,nch_tor_connect,ichtor,testcord,ibox,torkey
     &           ,angstart,phishift)
            if ( lfailure ) return
         elseif ( cbmc_form_code .eq. CBMC_FORM_MF2006 ) then
c           --- bond and angle and torsion_in are decoupled from each
c           --- other, but all of them are coupled in with torsion_out
c           --- we will make all of these calls down below as we 
c           --- cannot take advantage of only doing them once
         else
            write(6,*) 'ROSENBLUTH: unknown cbmc_form_code:'
     &           ,cbmc_form_code
            lfailure = .true.
            return
         endif
         do ip = 1,nch_nb
            if ( cbmc_form_code .eq. CBMC_FORM_MS1999_MT2004 ) then
c              --- just have to call a decoupled torsion selection 
c              --- for each of the nonbonded choices (coupling only
c              --- between nonbond and torsion).  There is no outer
c              --- torsion for this coupled-decoupled form

               if ( .not. lnew .and. ip .eq. 1 ) then
                  luseoldone = .true.
               else
                  luseoldone = .false.
               endif

c              --- perform the biased selection of the torsional angle
               call twh_getcbdihed(lfailure,ldihed_direct
     &              ,ldihed_indirect,luseoldone,ldihedterm,ldelta_dist
     &              ,testcord,iw,ichain,ichtor,iufrom,ibox,angstart
     &              ,ntogrow,imolty,dpscale_dihed,phishift,w_cbtor
     &              ,vcbvib,vcbbend,vcbtor,bfcbbias,xtorcord,ytorcord
     &              ,ztorcord,torkey)
               if ( lfailure ) return
#if DEBUG_ROSENBLUTH
      write(6,*) 'w_cbtor',w_cbtor
#endif      
               if ( ldihedterm ) then
c                 --- unsucessful attempt
                  w_tor(ip) = 0.0d0
                  dpscale_tor(ip) = 0.0d0
                  lnbterm(ip) = .true.
               else
                  lnbterm(ip) = .false.
c                 --- store the coordinates
                  do count = 1,ntogrow
                     rxp(count,ip) = xtorcord(count)
                     ryp(count,ip) = ytorcord(count)
                     rzp(count,ip) = ztorcord(count)
                  enddo
c                 --- store the total energies for each type.  This 
c                 --- includes energies computed up above in getcbbond
c                 --- and getcbangle which are the same for each nb
c                 --- choice
                  vxvibtr(ip) = vvibtr + gcba_vvib + vcbvib
                  vxbendtr(ip) = gcba_vbend + vcbbend
                  vtgtr(ip) = gcba_vtor + vcbtor
c                 --- save bias bolztman factor of this trial position
                  boltzbias(ip) = bfcbbias
c                 --- save boltzmann sum for torsion
                  w_tor(ip) = w_cbtor
                  dpscale_tor(ip) = dpscale_dihed
               endif
            elseif ( cbmc_form_code .eq. CBMC_FORM_MF2006 ) then
c              --- bond, angle, and inner torsion are decoupled from
c              --- each other, but all are coupled into outer torsion
c              --- outer torsion is also coupled to nonbond
               lfound = .false.
               do iouter = 1,nch_pre_nb
                  if ( .not. lnew .and. ip .eq. 1 .and. iouter .eq. 1
     &                 ) then
c                    --- select the old conformation
                     luseoldone = .true.
                  else
C                    --- all trial sites are newly generated
                     luseoldone = .false.
                  endif

c                 --- perform the biased selection to get lengths
                  call twh_getcbbond(lfailure,luseoldone,testcord,iw
     &                 ,ichain,imolty,iuprev,nch_vib,vvibtr,maxlen
     &                 ,w_cbbond,dpscale_bond,ibox)
                  if ( lfailure ) return
c                 --- perform the biased selection to get bond angles
                  call twh_getcbangle(lfailure,luseoldone,lovrangle
     &                 ,testcord,iw
     &                 ,ichain,imolty,angstart,iuprev,nch_bend_a
     &                 ,nch_bend_b,max_nch_bend,gcba_vvib,gcba_vbend
     &                 ,gcba_vtor,gcba_weight,dpscale_angle,ibox)
                  if ( lfailure ) return
                  if ( lovrangle ) then
c                    --- angle selection is failed so this whole trial
c                    --- is an overlap, kickout past dihedral selection
                     goto 100
                  endif
#if DEBUG_ROSENBLUTH
            write(6,*) 'gcba_weight:',gcba_weight
#endif
c                 --- setup all of the information needed to do the 
c                 --- dihedral selection
                  call twh_setcbdihed(lfailure,ldihed_direct
     &                 ,ldihed_indirect,ldelta_dist,ntogrow,iw,ichain
     &                 ,imolty,iuprev,iufrom,nch_tor,nch_tor_connect
     &                 ,ichtor,testcord,ibox,torkey,angstart
     &                 ,phishift)
                  if ( lfailure ) return
c                 --- perform the biased selection to get inner dihed
                  call twh_getcbdihed(lfailure,ldihed_direct
     &                 ,ldihed_indirect,luseoldone,ldihedterm
     &                 ,ldelta_dist,testcord,iw,ichain,ichtor,iufrom
     &                 ,ibox,angstart,ntogrow,imolty,dpscale_dihed
     &                 ,phishift,w_cbtor,vcbvib,vcbbend,vcbtor,bfcbbias
     &                 ,xtorcord,ytorcord,ztorcord,torkey)
                  if ( lfailure ) return
c                 --- check for walk termination
 100              if ( lovrangle .or. ldihedterm ) then
                     lovrprenb(iouter) = .true.
                  else
                     lovrprenb(iouter) = .false.
c                    --- store the weight products
                     w_tout(iouter) = w_cbbond*gcba_weight*w_cbtor
                     dpscale_prenb(iouter) = dpscale_bond 
     &                    + dpscale_angle + dpscale_dihed
c                    --- adjust dpscale_tor if needed
                     if ( .not. lfound ) then
c                       --- set dpscale_tor to dpscale_prenb
                        dpscale_tor(ip) = dpscale_prenb(iouter)
                     else
c                       --- check if dpscale_prenb is larger than 
c                       --- dpscale_tor
                        if ( dpscale_tor(ip) .gt. dpscale_prenb(iouter) 
     &                       ) then
                           dpscale_tor(ip) = dpscale_prenb(iouter)
                        endif
                     endif
c                    --- set logical so we know we found one trial that
c                    --- is not a hard overlap
                     lfound = .true.
                  endif
c                 --- store the coordinates
                  do count = 1,ntogrow
                     rx_tout(iouter,count) = xtorcord(count)
                     ry_tout(iouter,count) = ytorcord(count)
                     rz_tout(iouter,count) = ztorcord(count)
                  enddo
c                 --- store the extra bias factors
                  bias_tout(iouter) = bfcbbias
c                 --- store the energies
                  vvib_tout(iouter) = vvibtr + gcba_vvib + vcbvib
                  vbend_tout(iouter) = gcba_vbend + vcbbend
                  vtor_tout(iouter) = vcbtor + gcba_vtor
               enddo
c              --- select the outer torsion in a biased fashion
c              --- initialze weight for this selection
               w_tor(ip) = 0.0d0
               if ( .not. lfound ) then
c                 --- all trials were hard overlaps
                  lnbterm(ip) = .true.
               else
                  lnbterm(ip) = .false.
                  do iouter = 1,nch_pre_nb
c                    --- incorporate the adjusted dpscale_tor into the
c                    --- weights
                     if ( lovrprenb(iouter) ) then
                        w_tout(iouter) = 0.0d0
                     else
                        w_tout(iouter) = w_tout(iouter)
     &                       *scalelolog**(dpscale_prenb(iouter)
     &                       - dpscale_tor(ip))
                     endif
c                    --- add to w_tor
                     w_tor(ip) = w_tor(ip) + w_tout(iouter)
                  enddo
                  if ( .not. lnew .and. ip .eq. 1 ) then
c                    --- choose the first selection
                     iouter = 1
                  else
c                    --- randomly select one of the trials
                     rbf = w_tor(ip)*twh_random()
                     bs = 0.0d0
                     iouter = 0
                     do while ( iouter .le. nch_pre_nb
     &                    .and. rbf .gt. bs )
                        iouter = iouter + 1
                        bs = bs + w_tout(iouter)
                     enddo
c                    --- safety check
                     if ( iouter .gt. nch_pre_nb ) then
                        write(6,*) 'ROSENBLUTH: problem in pre-nonbond'
                        lfailure = .true.
                        return
                     endif
                  endif
c                 --- store the energies for the nonbonded selections
                  vxvibtr(ip) = vvib_tout(iouter)
                  vxbendtr(ip) = vbend_tout(iouter)
                  vtgtr(ip) = vtor_tout(iouter)
c                 --- save bias bolztman factor of this trial position
                  boltzbias(ip) = bias_tout(iouter)
c                 --- store the coordinates
                  do count = 1,ntogrow
                     rxp(count,ip) = rx_tout(iouter,count)
                     ryp(count,ip) = ry_tout(iouter,count)
                     rzp(count,ip) = rz_tout(iouter,count)
                  enddo
c                 --- divide out the number of outer trials
                  w_tor(ip) = w_tor(ip) / dble(nch_pre_nb)
               endif
            else
               write(6,*) 'ROSENBLUTH: unknown cbmc_form_code:'
     &              ,cbmc_form_code
               lfailure = .true.
               return
            endif
         enddo

c        --- now that we have the trial sites, need to compute 
c        --- non-bond energy
         call twh_engatom(lfailure,testcord,globalcord,engastyle,lovrlap
     &        ,ichain,imolty,ibox,nch_nb,iufrom,ntogrow,maxlen,lnbterm
     &        ,NNBOND,nch_nb,rxp,ryp,rzp,vctry,vtrext,vtrintra
     &        ,vtrinter,vtrelect,vtrewald,vqtry)
         if ( lfailure ) return
         if ( lovrlap ) then
#if DEBUG_ROSENBLUTH
      write(6,*) 'engatom termination'
#endif
            lterm = .true.
            if ( .not. lnew ) then
               write(6,*) 'major problem in rosenbluth'
               write(6,*) 'hard overlap in intermolecular'
               lfailure = .true.
               return
            endif
            return
         endif
c        --- compute bfadjust using getweight
         lfound = .false.
         do ip = 1, nch_nb
            if ( lnbterm(ip) ) then
c              --- hard overlap, zero bfadjust
               bfadjust(ip) = 0.0d0
            else
               vmin = vctry(ip)
               lsoftovr = .false.
               call twh_getweight(lfailure,0,vmin,wmin,dpmin,lsoftovr)
               if ( lfailure ) return
               if ( lsoftovr ) then
c                 --- the energies were far too high in getweight
c                 --- treat this as an overlap
                  lnbterm(ip) = .true.
                  bfadjust(ip) = 0.0d0
               else
                  bfadjust(ip) = wmin
                  dpscale_adjust(ip) = dpmin
               endif
            endif
c           --- set lfound to true if we don't have a hard or soft 
c           --- overlap
            if ( .not. lnbterm(ip) ) lfound = .true.
         enddo
         if ( .not. lfound ) then
            if ( .not. lnew ) then
               write(6,*) 'ROSENBLUTH: nonbond overlap for old'
               lfailure = .true.
               return
            endif
            lterm = .true.
            return
         endif
c        --- find the lowest number of log scaling (avoid the lovr)
         lfound = .false.
         dpscale_nb = 0.0d0
         do ip = 1,nch_nb
            if ( .not. lnbterm(ip) ) then
               if ( lfound ) then
c                 --- we have already assigned a minumum, do comparison
                  if (dpscale_nb .gt. 
     &                 dpscale_tor(ip)+dpscale_adjust(ip) ) then
                     dpscale_nb = dpscale_tor(ip)+dpscale_adjust(ip)
                  endif
               else
                  dpscale_nb = dpscale_tor(ip)+dpscale_adjust(ip)
                  lfound = .true.
               endif
            endif
         enddo
c        --- compute the nonbonded rosenbluth weight
         w_nb = 0.0d0
         lterm = .true.
         do ip = 1,nch_nb
            if ( lnbterm(ip) ) then
               bfactor(ip) = 0.0d0
            else
               lterm = .false.
               bfactor(ip) = w_tor(ip)*bfadjust(ip)
     &              *scalelolog**(dpscale_tor(ip)+dpscale_adjust(ip)
     &              -dpscale_nb)
            endif
            w_nb = w_nb + bfactor(ip)
#if DEBUG_ROSENBLUTH
            write(6,*) 'ip,bfactor(ip):',ip,bfactor(ip)
#endif
         enddo
c        --- possible that all were hard overlaps
         if ( lterm ) then
#if DEBUG_ROSENBLUTH
      write(6,*) 'nonbond termination'
#endif
            if ( .not. lnew ) write(6,*) 'ROSENBLUTH: overlap in old'
            return
         endif
         if ( w_nb .lt. 1d-300 ) then
            write(6,*) 'ROSENBLUTH: '
     &           ,'Danger of underflow in W_nonbond',w_nb
            write(6,*) 'lnew',lnew
            do count = 1,ntogrow
               call twh_grownbon(GLB_GET,iw,count,grownbon)
               write(6,*) 'iw,count,grownbon',iw,count,grownbon
               call twh_nxangle(GLB_GET,count,nxangle)
               write(6,*) 'count,nxangle',count,nxangle
               call twh_nxbias(GLB_GET,count,nxbias)
               write(6,*) 'count,nxbias',count,nxbias
            enddo
            do ip = 1,nch_nb
               write(6,*) 'ip,bfadjust,w_tor,bfactor',ip,bfadjust(ip)
     &              ,w_tor(ip),bfactor(ip)
               write(6,*) 'lnbterm,vctry',lnbterm(ip),vctry(ip)
               write(6,*) 'vtrintra,vtrelect',vtrintra(ip),vtrelect(ip)
            enddo
         endif
c        --- perform the walk according to the availibility of the 
c        --- choices and calculate the correct weight for the trial walk
         if ( lnew ) then
c           --- select one position at random according to boltzman sum
            rbf = w_nb*twh_random()
            bs = 0.0d0 
            do ip = 1, nch_nb
               bs = bs + bfactor(ip)
               if ( rbf .lt. bs ) then
c                 --- select ip position ---
                  iwalk = ip
                  if ( lnbterm(ip) ) then
                     write(6,*) 'problem in Rosenbluth'
                     lfailure = .true.
                     return
                  endif
                  goto 180
               endif
            enddo
c           --- if we got to this point the w_nb sum must be 0
            lterm = .true.
            write(6,*) 'ROSENBLUTH: final w_nb all zero'
            write(6,*) 'lnew,w_nb,rbf',lnew,w_nb,rbf
            lfailure = .true.
            return
         else
            iwalk = 1
         endif

 180     continue

c        --- divide out the biasing boltzmann factor
         w_nb = w_nb / boltzbias(iwalk)
c        --- divide out the number of nonbond trials
         w_nb = w_nb / dble(nch_nb)

#if DEBUG_ROSENBLUTH
         write(6,*) 'Rosenbluth w_nb:',w_nb
#endif

         if ( cbmc_form_code .eq. CBMC_FORM_MS1999_MT2004 ) then
c           --- update rosenbluth weight - include angle and
c           --- bond weights
            wrosen = wrosen * w_nb * gcba_weight * w_cbbond 
c           --- add in the scaling
            dpscale = dpscale + dpscale_nb + dpscale_bond 
     &           + dpscale_angle
         elseif ( cbmc_form_code .eq. CBMC_FORM_MF2006 ) then
c           --- update rosenbluth weight - only nonbonded wieght
            wrosen = wrosen * w_nb
            dpscale = dpscale + dpscale_nb
         else
            write(6,*) 'ROSENBLUTH: unknown cbmc_form_code:'
     &           ,cbmc_form_code
            lfailure = .true.
            return
         endif

c        --- rescale the weight
         call twh_getweight(lfailure,1,vmin,wrosen,dpmin,lterm)
         if ( lfailure ) return
         dpscale = dpscale + dpmin
         if ( lterm ) then
            if ( lnew ) then
               write(6,*) 'Rosenbluth: safety return'
               return
            else
               write(6,*) '###old weight too low'
               write(6,*) 'w_nb,wmin',w_nb,wmin
            endif
         endif
#if DEBUG_ROSENBLUTH
      write(6,*) 'iwalk,vctry',iwalk,vctry(iwalk)
      write(6,*) 'vtgtr',vtgtr(iwalk)
      write(6,*) 'vxvibtr',vxvibtr(iwalk)
      write(6,*) 'vxbendtr',vxbendtr(iwalk)
      write(6,*) 'vtrext',vtrext(iwalk)
      write(6,*) 'vtrintra',vtrintra(iwalk)
      write(6,*) 'vtrinter',vtrinter(iwalk)
      write(6,*) 'vtrelect',vtrelect(iwalk)
      write(6,*) 'vtrewald',vtrewald(iwalk)
#endif
c        --- update new trial energies
         if ( lnew ) then
            istoreflag = CRD_NEW
         else
            istoreflag = CRD_OLD
         endif
         if ( potentialstyle .eq. POT_INTERNAL ) then
            vrostot = vrostot + vctry(iwalk)  
     &           +vtgtr(iwalk) + vxvibtr(iwalk) + vxbendtr(iwalk)
            vrosvib = vrosvib + vxvibtr(iwalk)
            vrosang = vrosang + vxbendtr(iwalk)
            vrostor = vrostor + vtgtr(iwalk)
            vrosext = vrosext + vtrext(iwalk)
            vrosintra = vrosintra + vtrintra(iwalk)
            vrosinter = vrosinter + vtrinter(iwalk)
            vroselect = vroselect + vtrelect(iwalk)
            vrosewald = vrosewald + vtrewald(iwalk)
         elseif ( potentialstyle .eq. POT_EXTERNAL ) then
            vrostot = vrostot + vqtry(iwalk)
         else
            write(6,*) 'ROSENBLUTH: unknown potentialstyle'
            write(6,*) potentialstyle
            lfailure = .true.
            return
         endif

         do count = 1,ntogrow
            call twh_growlist(GLB_GET,iw,count,iu)
c           --- set lexist(iu) to true so ewald sum computed properly
            call twh_lexist(GLB_SET,iu,.true.)
            if ( lnew ) then
c              --- assign new positions to r*new
               xtemp = rxp(count,iwalk)
               ytemp = ryp(count,iwalk)
               ztemp = rzp(count,iwalk)
               call twh_coordinates(GLB_SET,testcord,ichain,iu
     &              ,xtemp,ytemp,ztemp)
#if DEBUG_ROSENBLUTH
      write(6,*) 'new iu,iwalk,coords',iu,iwalk,xtemp,ytemp,ztemp
#endif
            endif
         enddo
      enddo
#if DEBUG_ROSENBLUTH
      write(6,*) 'end ROSENBLUTH',lnew
#endif
      return
      end
