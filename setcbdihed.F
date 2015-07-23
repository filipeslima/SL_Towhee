#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2011 Marcus G. Martin                       *
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
      subroutine twh_setcbdihed(lfailure,ldihed_direct,ldihed_indirect
     &     ,ldelta_dist,ntogrow,iw,ichain,imolty,iuprev,iufrom
     &     ,nch_tor,nch_tor_connect,ichtor
     &     ,testcord,ibox,torkey,angstart,phishift)
c     ******************************************************************
c     * sets up the logic arrays needed to do the dihedral selection   *
c     * for configurational-bias                                       *
c     *                                                                *
c     * originally split from rosenbluth 08-07-2002 by M.G. Martin     *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_BONDPATT
#define FUNCTION_GET_IJBEN2
#define FUNCTION_GET_IJBEN3
#define FUNCTION_GET_IJIMPROP2
#define FUNCTION_GET_IJIMPROP3
#define FUNCTION_GET_IJIMPROP4
#define FUNCTION_GET_IJTOR2
#define FUNCTION_GET_IJTOR3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_INBEN
#define FUNCTION_GET_INIMPROP
#define FUNCTION_GET_INTOR
#define FUNCTION_GET_LEXIST
#define FUNCTION_GET_NTYPE
#define FUNCTION_VTORSION
#include "functions.h"
c     --- debug toggle
#define DEBUG_SETCBDIHED 0
c     --- variables passed to/from the subroutine
      logical lfailure,ldihed_direct,ldihed_indirect,ldelta_dist
      integer ntogrow,iw,imolty,iuprev,iufrom,nch_tor,nch_tor_connect
      integer ichtor
      integer testcord,ichain,ibox,torkey,angstart
      double precision phishift
c     --- local variables
c     --- logical scalars
      logical lstore,lbias,lcompute,lovrlap
c     --- character string scalars
      character*5 hyfrom,hyprev
c     --- integer scalars
      integer count,iu,ibend,jut1,jut2,jut3,jut4,ixtra,ibias,it,imp
      integer itemp,iutwo,iuthree,iufour,ibond,iubond,iaa,ibead
      integer itest,ttype,iang,ivalue,rimptotal,total,growprev,iflip
      integer grownbon,xtracount,nxangle,nxbias,naabead,nimpbead
      integer cbmc_dihedgen_code,tortype,torstyle
c     --- integer arrays
      integer torbead
      dimension torbead(4)
c     --- double precision scalars
      double precision phinew,thetanew,dum,ux,uy,uz,length,xxx,yyy,zzz
      double precision rxfrom,ryfrom,rzfrom,xpos,ypos,zpos
c     --- double precision arrays
      double precision xttmp,yttmp,zttmp
      dimension xttmp(4),yttmp(4),zttmp(4)

#if DEBUG_SETCBDIHED
      write(6,*) 'starting setcbdihed:',lfailure,ldihed_direct
     &,ldihed_indirect,ntogrow,iw,ichain,imolty,iuprev,iufrom,ichtor
     &,testcord,ibox,phishift,torkey,angstart
#endif

c     --- initialize torkey and phishift
      torkey = 0
      phishift = 0.0d0
c     --- initialize logicals
      ldihed_direct = .false.
      ldihed_indirect = .false.
      ldelta_dist = .false.
c     --- pretabulate all of the information that is independent 
c     --- of the torsional selection process below
      do count = 1,ntogrow
         call twh_growlist(GLB_GET,iw,count,iu)
c        --- initialize nxangle and nxbias
         call twh_nxangle(GLB_SET,count,0)
         call twh_nxbias(GLB_SET,count,0)
c        --- set up growbond and xtra angle information
c        --- where iu is one of the end beads
         do ibend = 1,twh_get_inben(imolty,iu)
c           --- account for angles with the forms
c           --- iu - growbond - xx
c           --- iu - iufrom - growxtra 
            lstore = .false.
            lbias = .false.
            jut2 = twh_get_ijben2(imolty,iu,ibend)
            jut3 = twh_get_ijben3(imolty,iu,ibend)
            
            if ( twh_get_lexist(jut3) .and. jut3 .ne. iuprev 
     &           .and. jut3 .ne. iufrom ) then
               if ( jut2 .eq. iufrom ) then
c                --- see if jut3 is one of the xtra beads
                  call twh_growxtra(GLB_GET,iw,1,xtracount)
                  do ixtra = 2,xtracount+1
                     call twh_growxtra(GLB_GET,iw,ixtra,itest)
                     if ( jut3 .eq. itest ) lstore = .true.
                  enddo
               elseif ( twh_get_lexist(jut2) ) then
c                --- jut2 must be one of the growbond atoms
                  lstore = .true.
               else
c                --- jut2 does not yet exist, but we want to use 
c                --- this angle to do some biasing
                  lbias = .true.
               endif
            endif
            if ( lstore ) then
               call twh_nxangle(GLB_INCR,count,1)
               call twh_nxangle(GLB_GET,count,nxangle)
               iflip =  -count
               call twh_xangle(GLB_SET,1,count,nxangle,iflip)
               call twh_xangle(GLB_SET,2,count,nxangle,jut2)
               call twh_xangle(GLB_SET,3,count,nxangle,jut3)
               call twh_txangle(GLB_SET,count,nxangle,ibend)
            elseif ( lbias ) then
c              --- add this to the list of biasing sites
               call twh_nxbias(GLB_INCR,count,1)
               call twh_nxbias(GLB_GET,count,ibias)
               call twh_xbias(GLB_SET,1,count,ibias,1)
               call twh_xbias(GLB_SET,2,count,ibias,ibend)
            endif
         enddo

c        --- set up growbond and xtra angle information
c        --- where iu is the middle bead
         do ibend = 1,twh_get_inben(imolty,iufrom)
c           --- account for angles with the form
c           --- iufrom - iu - growbond
            jut1 = twh_get_ijben3(imolty,iufrom,ibend)
            if ( iu .eq. twh_get_ijben2(imolty,iufrom,ibend) .and.
     &           twh_get_lexist(jut1) ) then
c              --- add this to the data structure
               call twh_nxangle(GLB_INCR,count,1)
               call twh_nxangle(GLB_GET,count,nxangle)
               call twh_xangle(GLB_SET,1,count,nxangle,iufrom)
               iflip = -count
               call twh_xangle(GLB_SET,2,count,nxangle,iflip)
               call twh_xangle(GLB_SET,3,count,nxangle,jut1)
               call twh_txangle(GLB_SET,count,nxangle,ibend)
            endif
         enddo

c        --- set up the biasing information for torsions which
c        --- do not yet exist, but have both endpoints in place
         do it = 1, twh_get_intor(imolty,iu)
            jut2 = twh_get_ijtor2(imolty,iu,it)
            jut3 = twh_get_ijtor3(imolty,iu,it)
            jut4 = twh_get_ijtor4(imolty,iu,it)
            if ( (.not. twh_get_lexist(jut2))
     &           .and. (.not. twh_get_lexist(jut3))
     &           .and. twh_get_lexist(jut4) ) then
               call twh_nxbias(GLB_INCR,count,1)
               call twh_nxbias(GLB_GET,count,ibias)
               call twh_xbias(GLB_SET,1,count,ibias,2)
               call twh_xbias(GLB_SET,2,count,ibias,it)
            endif
         enddo
      enddo

c     --- set up the improper torsion information
c     --- add the impropers centered on iufrom into the data struc
      call twh_nimpbead(GLB_SET,0)
      call twh_cbmcrimptype(GLB_GET,0,rimptotal)
      do imp = 1,rimptotal
c        --- increment the number of improper torsions
         call twh_nimpbead(GLB_INCR,1)
         call twh_nimpbead(GLB_GET,nimpbead)
c        --- store the information about this improper
         call twh_cbmcrimptype(GLB_GET,imp,itemp)
c        --- sign convention is opposite from geometry's
         call twh_impbead(GLB_SET,nimpbead,1,iufrom)
         call twh_cbmcimproper(GLB_GET,2,itemp,ivalue)
         iflip = -ivalue
         call twh_impbead(GLB_SET,nimpbead,2,iflip)
         call twh_cbmcimproper(GLB_GET,3,itemp,ivalue)
         iflip = -ivalue
         call twh_impbead(GLB_SET,nimpbead,3,iflip)
         call twh_cbmcimproper(GLB_GET,4,itemp,ivalue)
         iflip = -ivalue
         call twh_impbead(GLB_SET,nimpbead,4,iflip)
         call twh_impropernum(GLB_SET,nimpbead,itemp)
      enddo

c     --- see if any of the growing beads fills in the center of 
c     --- an improper torsion 
      do count = 1,ntogrow
         call twh_growlist(GLB_GET,iw,count,iu)
         do imp = 1,twh_get_inimprop(imolty,iu)
            iutwo = twh_get_ijimprop2(imolty,iu,imp)
            iuthree = twh_get_ijimprop3(imolty,iu,imp)
            iufour = twh_get_ijimprop4(imolty,iu,imp)
            if ( twh_get_lexist(iutwo) .and. twh_get_lexist(iuthree) 
     &           .and. twh_get_lexist(iufour) ) then
c              --- increment the number of improper torsions
               call twh_nimpbead(GLB_INCR,1)
               call twh_nimpbead(GLB_GET,nimpbead)
c              --- store the information about this improper
               iflip = -count
               call twh_impbead(GLB_SET,nimpbead,1,iflip)
               call twh_impbead(GLB_SET,nimpbead,2,iutwo)
               call twh_impbead(GLB_SET,nimpbead,3,iuthree)
               call twh_impbead(GLB_SET,nimpbead,4,iufour)
               call twh_impropernum(GLB_SET,nimpbead,imp)
            endif
         enddo
      enddo

c     --- see if the growbond atoms are centers of imp torsions
      do count = 1,ntogrow
         call twh_growlist(GLB_GET,iw,count,iu)
         call twh_grownbon(GLB_GET,iw,count,grownbon)
         do ibond = 1,grownbon
            call twh_growbond(GLB_GET,iw,count,ibond,iubond)
            do imp = 1,twh_get_inimprop(imolty,iubond)
               iutwo = twh_get_ijimprop2(imolty,iubond,imp) 
               iuthree = twh_get_ijimprop3(imolty,iubond,imp)
               iufour = twh_get_ijimprop4(imolty,iubond,imp)
               if ( iutwo .eq. iu .and. twh_get_lexist(iuthree) 
     &              .and. twh_get_lexist(iufour) ) then
c                 --- increment the number of improper torsions
                  call twh_nimpbead(GLB_INCR,1)
                  call twh_nimpbead(GLB_GET,nimpbead)
c                 --- store the information about this improper
                  call twh_impbead(GLB_SET,nimpbead,1,iubond)
                  iflip = -count
                  call twh_impbead(GLB_SET,nimpbead,2,iflip)
                  call twh_impbead(GLB_SET,nimpbead,3,iuthree)
                  call twh_impbead(GLB_SET,nimpbead,4,iufour)
                  call twh_impropernum(GLB_SET,nimpbead,imp)
               elseif ( twh_get_lexist(iutwo) .and. iuthree .eq. iu 
     &                 .and. twh_get_lexist(iufour) ) then
c                 --- increment the number of improper torsions
                  call twh_nimpbead(GLB_INCR,1)
                  call twh_nimpbead(GLB_GET,nimpbead)
c                 --- store the information about this improper
                  call twh_impbead(GLB_SET,nimpbead,1,iubond)
                  call twh_impbead(GLB_SET,nimpbead,2,iutwo)
                  iflip = -count
                  call twh_impbead(GLB_SET,nimpbead,3,iflip)
                  call twh_impbead(GLB_SET,nimpbead,4,iufour)
                  call twh_impropernum(GLB_SET,nimpbead,imp)
               elseif ( twh_get_lexist(iutwo)
     &                 .and. twh_get_lexist(iuthree)
     &                 .and. iufour .eq .iu ) then
c                 --- increment the number of improper torsions
                  call twh_nimpbead(GLB_INCR,1)
                  call twh_nimpbead(GLB_GET,nimpbead)
c                 --- store the information about this improper
                  call twh_impbead(GLB_SET,nimpbead,1,iubond)
                  call twh_impbead(GLB_SET,nimpbead,2,iutwo)
                  call twh_impbead(GLB_SET,nimpbead,3,iuthree)
                  iflip = -count
                  call twh_impbead(GLB_SET,nimpbead,4,iflip)
                  call twh_impropernum(GLB_SET,nimpbead,imp)
               endif
            enddo
         enddo
      enddo

c     --- set up the angle-angle information
c     --- add the angle-angles centered on iufrom into the data struc
      call twh_naabead(GLB_SET,0)
      do count = 1,ntogrow
         call twh_cbmcaatot(GLB_GET,count,3,total)
         do iaa = 1,total
c           --- increment the number of angle-angle terms
            call twh_naabead(GLB_INCR,1)
            call twh_naabead(GLB_GET,naabead)
c           --- store the information about this angle-angle
            call twh_cbmcaatype(GLB_GET,count,3,iaa,itemp)
            call twh_aabead(GLB_SET,naabead,1,iufrom)
            call twh_cbmcaalocl(GLB_GET,0,itemp,ivalue)
            call twh_aabead(GLB_SET,naabead,2,ivalue)
            call twh_cbmcaalocl(GLB_GET,1,itemp,ivalue)
            call twh_aabead(GLB_SET,naabead,3,ivalue)
            call twh_cbmcaalocl(GLB_GET,2,itemp,ivalue)
            call twh_aabead(GLB_SET,naabead,4,ivalue)
            call twh_aaind(GLB_SET,naabead,itemp)
         enddo
      enddo
c     --- reprocess the angle-angle terms above so that they follow
c     --- the convention that units which already exist are positive
c     --- integers with the integer corresponding to the unit number
c     --- atoms grown during this step are negative number with the
c     --- absolute value of the number coresponding to the count 
c     --- number
      call twh_naabead(GLB_GET,naabead)
      do iaa = 1,naabead
         do ibead = 2,4
            call twh_aabead(GLB_GET,iaa,ibead,itemp)
            if ( itemp .eq. -1 ) then
c              --- this is the iuprev
               call twh_aabead(GLB_SET,iaa,ibead,iuprev)
            elseif ( itemp .lt. -1 ) then
c              --- one of the xtra atoms
               iflip = -itemp
               call twh_growxtra(GLB_GET,iw,iflip,xtracount)
               call twh_aabead(GLB_SET,iaa,ibead,xtracount)
            else
c              --- one of the count atoms, change sign
               iflip = -itemp
               call twh_aabead(GLB_SET,iaa,ibead,iflip)
            endif
         enddo
      enddo

c     --- see if there are any interactions that need to be computed
c     --- this can come either from torsions along iufrom-iuprev, 
c     --- , or bonds, torsions, and angles that involve the xtra or
c     --- growbond atoms.
      call twh_growprev(GLB_GET,iw,growprev)
      if ( growprev .ne. 0 ) then
c        --- see if there is the possiblity for existing torsions
         do iang = 1,twh_get_inben(imolty,iufrom)
            if ( twh_get_ijben2(imolty,iufrom,iang) .eq. iuprev 
     &           .and. 
     &           twh_get_lexist(twh_get_ijben3(imolty,iufrom,iang)) 
     &           ) ldihed_direct = .true.
         enddo
      endif
      do count = 1,ntogrow
         call twh_nxangle(GLB_GET,count,nxangle)
         if ( nxangle .ne. 0 ) ldihed_indirect = .true.
         call twh_nxbias(GLB_GET,count,nxbias)
         if ( nxbias .ne. 0 ) ldihed_indirect = .true.
      enddo
      call twh_nimpbead(GLB_GET,nimpbead)
      if ( nimpbead .ne. 0 ) ldihed_indirect = .true.
      call twh_naabead(GLB_GET,naabead)
      if ( naabead .ne. 0 ) ldihed_indirect = .true.

c     --- determine how many torsional selections to attempt
c     --- default is the regular torsion trial count
         ichtor = nch_tor
      call twh_growxtra(GLB_GET,iw,1,xtracount)
      if ( xtracount .ne. 0 ) then
c        --- connecting back to existing atoms, use nch_tor_connect
         ichtor = nch_tor_connect
      endif
      do count = 1,ntogrow
         call twh_grownbon(GLB_GET,iw,count,grownbon)
         call twh_nxbias(GLB_GET,count,nxbias)
         if ( grownbon .ne. 0 .or. nxbias .ne. 0) then
            ichtor = nch_tor_connect
         endif
      enddo

c     --- if there are no interaction we only need to run through
c     --- one time
      if ( (.not. ldihed_direct) .and. (.not. ldihed_indirect) ) then
         ichtor = 1
      endif
      call twh_cbmc_dihedgen_code(GLB_GET,cbmc_dihedgen_code)
      if ( ldihed_direct ) then
c        --- determine function to generate trial positions.  
c        --- Need to know the phishift that brings the
c        --- x1-iuprev-iufrom-count(1) angle to be zero
c        --- this is not needed for a zero growprev as in that case
c        --- we just use uniform
         if ( cbmc_dihedgen_code .eq. CBMC_GLOBAL_GAUSSIAN ) then
c           --- store the hybridizations of iufrom and iuprev
            hyfrom = twh_get_bondpatt(twh_get_ntype(imolty,iufrom))
            call twh_hyfrom(GLB_SET,hyfrom)
            hyprev = twh_get_bondpatt(twh_get_ntype(imolty,iuprev))
            call twh_hyprev(GLB_SET,hyprev)
         endif
c        --- determine the dihedral angle between the first 
c        --- atom that exists bonded to iuprev which is not 
c        --- iufrom, and the angstart atom.  This will be used
c        --- to bias the selection of the torsional angles
         count = angstart
         call twh_growlist(GLB_GET,iw,count,iu)
         call twh_cbphi(GLB_GET,count,phinew)
         call twh_bendang(GLB_GET,count,thetanew)
         call twh_cone(lfailure,2,xxx,yyy,zzz,thetanew,phinew,ux,uy,uz)
         if ( lfailure ) return
c        --- store the cartesian postions
         call twh_bondlen(GLB_GET,count,length)
         call twh_coordinates(GLB_GET,testcord,ichain,iufrom,rxfrom
     &        ,ryfrom,rzfrom)
c        --- store these locally as *pos
         xpos = rxfrom + ux*length
         ypos = ryfrom + uy*length
         zpos = rzfrom + uz*length
c        --- initialize lcompute
         lcompute = .false.
c        --- check all of the torsions starting at iu
         it = 0
         do while ( it .lt. twh_get_intor(imolty,iu) 
     &        .and. .not. lcompute ) 
            it = it + 1
c           --- assign the unit types for the torsion
            torbead(1) = iu
            torbead(2) = twh_get_ijtor2(imolty,iu,it)
            torbead(3) = twh_get_ijtor3(imolty,iu,it)
            torbead(4) = twh_get_ijtor4(imolty,iu,it)
c           --- determine if this torsion should be computed
            if ( torbead(2) .eq. iufrom .and.
     &           torbead(3) .eq. iuprev .and.
     &           twh_get_lexist(torbead(2)) .and.
     &           twh_get_lexist(torbead(3)) .and.
     &           twh_get_lexist(torbead(4)) ) then
c              --- compute this torsion
               lcompute =.true.
c              --- temporarily switch back to count
               torbead(1) = -count
            endif
            if ( lcompute ) then
c              --- assign torkey
               torkey = it
               if ( torkey .gt. MAXTOR ) then
                  write(6,*) 'SETCBDIHED: invalid torkey'
                  write(6,*) 'torkey',torkey
                  write(6,*) 'torbead:',torbead
                  write(6,*) 'imolty,iu,intor',imolty,iu
     &                 ,twh_get_intor(imolty,iu)
               endif
c              --- assign the positions
               do ibead = 1,4
                  itest = torbead(ibead)
                  if ( itest .lt. 0 ) then
c                    --- use the count coordinates
                     iflip = -itest
                     xttmp(ibead) = xpos
                     yttmp(ibead) = ypos
                     zttmp(ibead) = zpos
c                    --- reset torbead
                     call twh_growlist(GLB_GET,iw,iflip,torbead(ibead))
                  else
c                    --- use old/new unit coordinates
                     call twh_coordinates(GLB_GET,testcord,ichain,itest
     &                    ,xttmp(ibead),yttmp(ibead),zttmp(ibead))
                  endif
               enddo
c              --- store unit numbers
               ttype = torbead(1)
c              --- compute the shift angle
               phishift = twh_vtorsion(lfailure,ibox,0,imolty
     &              ,ttype,it,xttmp,yttmp,zttmp,dum,lovrlap)
               if ( lfailure ) return
c              --- check to see if this is a multiple rigid angle
               call twh_ittor(GLB_GET,imolty,iu,it,tortype)
               call twh_torstyle(GLB_GET,tortype,torstyle)
               if ( torstyle .eq. 21 ) then
c                 --- multiple rigid, override other options
                  ldelta_dist = .true.
c                 --- reset the number of torsion trials to the number
c                 --- of values in the delta distribution
                  call twh_ntorloop(GLB_GET,tortype,ichtor)
               endif
            endif
         enddo
      endif

      return
      end
