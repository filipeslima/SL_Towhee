#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2012 Marcus G. Martin                       *
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
      subroutine twh_getcbdihed(lfailure,ldihed_direct,ldihed_indirect
     &     ,luseoldone,lterm,ldelta_dist,testcord,iw,ichain,ichtor
     &     ,iufrom,ibox,angstart,ntogrow,imolty,dpcbtor,phishift
     &     ,w_cbtor,vcbvib,vcbbend,vcbtor,bfcbbias
     &     ,xtorcord,ytorcord,ztorcord,torkey)
c     ******************************************************************
c     * performs a biased dihedral selection                           *
c     *                                                                *
c     * written by M.G. Martin 08-08-2002 from rosenbluth              *
c     * last modified 08-30-2012 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_EXPON
#define FUNCTION_FEBIAS
#define FUNCTION_GAUSSIAN
#define FUNCTION_GAUSSPROB
#define FUNCTION_GET_IJBEN3
#define FUNCTION_GET_IJTOR2
#define FUNCTION_GET_IJTOR3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_IJVIB
#define FUNCTION_GET_INTOR
#define FUNCTION_GET_INVIB
#define FUNCTION_GET_LEXIST
#define FUNCTION_RANDOM
#define FUNCTION_VANGANG
#define FUNCTION_VANGLE
#define FUNCTION_VBOND
#define FUNCTION_VIMPROPER
#define FUNCTION_VTORSION
#include "functions.h"
#define DEBUG_GETCBDIHED 0
c     --- variables passed to/from the subroutine
      logical lfailure,ldihed_direct,ldihed_indirect,lterm,luseoldone
      logical ldelta_dist
      integer iw,ichtor,iufrom,ibox,testcord,angstart,ntogrow,imolty
      integer ichain,torkey
      double precision w_cbtor,vcbvib,vcbbend,vcbtor,bfcbbias
      double precision dpcbtor,phishift
      double precision xtorcord,ytorcord,ztorcord
      dimension xtorcord(NNBOND),ytorcord(NNBOND),ztorcord(NNBOND)
c     --- local variables
c     --- logical scalars
      logical lovrlap,lbadstereo,lmin,lfirstold,lcompute,lfound,lgetdist
c     --- logical arrays
      logical lovrtor
      dimension lovrtor(ichtor)
c     --- integer scalars
      integer itor,count,iu,ixtra,jut2,ivib,itype,ibend,ibead,itest,iaa
      integer ttype,it,iend,iimprop,ibias,iunit,biastype,index,iuother
      integer imin,growprev,iflip,growtest,grownbon,nxangle,nxbias
      integer naabead,nimpbead,style_code,iwinner,npeak,delta_old,ipeak
c     --- integer arrays
      integer torbead
      dimension torbead(4)
      integer iutmp
      dimension iutmp(3)
c     --- double precision scalars
      double precision vmin,pdf_true,pdf_bias,x,y,z,length,phidisp
      double precision xub,yub,zub,dum,lolimit,hilimit,mean
      double precision phinew,thetanew,rxuij,ryuij,rzuij,anglec,dumb
      double precision dumc,vdha,vimprop,bs,wmin,std_dev,onepi
      double precision prob_g_region,prob_u_region,xpos,ypos,zpos
      double precision rxfrom,ryfrom,rzfrom,rxjut,ryjut,rzjut
      double precision rxgrow,rygrow,rzgrow,twopi,beta,rancheck,cbphi
      double precision ideal_frac,delta_diff
c     --- double precision arrays
      double precision bf_tor,vtors,ran_tor,vxbend,vxvib,vxaa
      dimension bf_tor(ichtor),vtors(ichtor)
      dimension vxbend(ichtor),vxvib(ichtor),vxaa(ichtor)
      double precision bf_bias
      dimension bf_bias(ichtor)
      double precision xatmp,yatmp,zatmp
      dimension xatmp(3),yatmp(3),zatmp(3)
      double precision xttmp,yttmp,zttmp
      dimension xttmp(4),yttmp(4),zttmp(4)
      double precision xcord,ycord,zcord
      dimension xcord(3),ycord(3),zcord(3)
      double precision dumangle
      dimension dumangle(2)
      double precision delta_dist
      dimension delta_dist(0:CBMCMAXTOR)
      double precision xcoord
      dimension xcoord(ichtor,NNBOND)
      double precision ycoord
      dimension ycoord(ichtor,NNBOND)
      double precision zcoord
      dimension zcoord(ichtor,NNBOND)

#if DEBUG_GETCBDIHED
      write(6,*) 'start GETCBDIHED'
#endif

c     --- retrieve constants
      call twh_constant_pi(GLB_GET,onepi)
      call twh_constant_twopi(GLB_GET,twopi)
      call twh_beta(GLB_GET,beta)

c     --- initialize lmin to show we haven't found a min yet
      lmin = .false.
      lterm = .false.
c     --- set a default vmin of zero
      vmin = 0.0d0
c     --- set delta_old to null
      delta_old = 0

c     --- set pdf_true to the uniform distribution
      pdf_true = 1.0d0/twopi
      do itor = 1,ichtor
c        --- initialize overlap and style
         lovrlap = .false.
         style_code = DIST_NULL
         if ( luseoldone .and. itor .eq. 1) then
            lfirstold = .true.
c           --- old conformation - set phidisp to 0.0d0
            phidisp = 0.0d0
         else
            lfirstold = .false.
         endif
c        --- initialize all energies and bias variables
         vxvib(itor) = 0.0d0
         vxbend(itor) = 0.0d0
         vxaa(itor) = 0.0d0
         vtors(itor) = 0.0d0
         bf_bias(itor) = 1.0d0
c        --- we have already set up the iufrom-iuprev cone except
c        --- for the case where growprev is 0
         call twh_growprev(GLB_GET,iw,growprev)
         if ( growprev .eq. 0 ) then
c           --- set up the cone based on iuprev (which is an atom
c           --- we grew this step).  Select a position for iuprev
c           --- from the unit sphere
            count = 1
            call twh_bondlen(GLB_GET,count,length)
            call twh_coordinates(GLB_GET,testcord,ichain,iufrom
     &           ,rxfrom,ryfrom,rzfrom)
            if ( lfirstold ) then
c              --- use old molecule position
               call twh_growlist(GLB_GET,iw,count,iu)
c              --- compute unit vector to be used in cone
               call twh_coordinates(GLB_GET,testcord,ichain,iu
     &              ,rxgrow,rygrow,rzgrow)
               call twh_mimage(ibox,rxgrow,rygrow,rzgrow
     &              ,rxfrom,ryfrom,rzfrom,x,y,z)
               x = x/length 
               y = y/length 
               z = z/length 
            else
c              --- choose randomly on the unit sphere
               call twh_sphere(x,y,z)
            endif
c           --- store the coordinates in xcoord
            xcoord(itor,count) = rxfrom + x*length
            ycoord(itor,count) = ryfrom + y*length
            zcoord(itor,count) = rzfrom + z*length
c           ---set up the cone 
            xub = -x 
            yub = -y 
            zub = -z 
            call twh_cone(lfailure,1,xub,yub,zub,dum,dum,dum,dum,dum )
            if ( lfailure ) return
c           --- choose phidisp uniformly
            style_code = DIST_UNIFORM
         else
            if ( lfirstold ) then
c              --- add in the phishift to get the angle we 
c              --- would have chosen from the biasing procedure
               phidisp = phishift
c              --- make sure this is still on (-Pi,Pi)
               if ( phidisp .gt. onepi ) then
                  phidisp = phidisp - twopi
               elseif ( phidisp .lt. -onepi ) then
                  phidisp = phidisp + twopi
               endif
            endif

c           --- determine which type of biasing to use to select
c           --- the phidisp
            call twh_growlist(GLB_GET,iw,1,iu)
            call twh_dihedral_distribution(lfailure,ldihed_direct
     &           ,ldelta_dist,style_code,npeak,delta_dist,lolimit
     &           ,hilimit,mean,std_dev,phidisp,prob_g_region,lfirstold
     &           ,imolty,iu,torkey)
            if ( lfailure ) return
         endif
#if DEBUG_GETCBDIHED
         write(6,*) 'style_code: ',style_code
#endif
         if ( style_code .eq. DIST_UNIFORM ) then
            if ( .not. lfirstold ) then
c              --- choose a displacement uniformly
               phidisp = twopi*twh_random()
            endif
            pdf_bias = 1.0d0/twopi
         elseif ( style_code .eq. DIST_GAUSSIAN ) then
c            --- gaussian distribution
            if ( .not. lfirstold ) then
               phidisp = twh_gaussian(lfailure,mean,std_dev
     &              ,lolimit,hilimit)
               if ( lfailure ) return
            endif
c           --- safety check
            if ( phidisp .lt. lolimit .or. phidisp .gt. hilimit ) then
               write(6,*) 'GETCBDIHED: phidisp out of bounds'
               lfailure = .true.
               return
            endif
            pdf_bias = twh_gaussprob(phidisp,mean,std_dev,lolimit
     &           ,hilimit)*prob_g_region
         elseif ( style_code .eq. DIST_UNIFORM_GAUSSIAN ) then
            call twh_dihedral_ideal_fraction(GLB_GET,ideal_frac)
            if ( .not. lfirstold ) then
               rancheck = twh_random()
               if ( rancheck .lt. ideal_frac ) then
c                 --- choose a displacement uniformly
                  phidisp = twopi*twh_random() - onepi
c                 --- find the appropriate bounds for biasing
                  lgetdist = .true.
                  call twh_dihedral_distribution(lfailure,ldihed_direct
     &                 ,ldelta_dist,style_code,npeak,delta_dist,lolimit
     &                 ,hilimit,mean,std_dev,phidisp,prob_u_region
     &                 ,lgetdist,imolty,iu,torkey)
                  if ( lfailure ) return
               else
c                 --- choose from a gaussian distribution
                  phidisp = twh_gaussian(lfailure,mean,std_dev
     &                 ,lolimit,hilimit)
                  if ( lfailure ) return
#if DEBUG_GETCBDIHED
                  write(6,*) 'mean,std_dev,lolimit,hilimit:'
     &                 ,mean,std_dev,lolimit,hilimit
#endif
c                 --- safety check
                  if ( phidisp .lt. lolimit .or. 
     &                 phidisp .gt. hilimit ) then
                     write(6,*) 'GETCBDIHED: phidisp out of bounds'
                     lfailure = .true.
                     return
                  endif
               endif
            endif
            pdf_bias = ideal_frac/twopi + (1.0d0-ideal_frac)
     &           *twh_gaussprob(phidisp,mean,std_dev,lolimit,hilimit)
     &           *prob_g_region
         elseif ( style_code .eq. DIST_DELTA ) then
            if ( lfirstold ) then
c              --- swap the order so that the actual old conformation
c              --- is the first entry
               ipeak = 0
               lfound = .false.
               do while ( ipeak .lt. npeak .and. .not. lfound )
                  ipeak = ipeak + 1
                  delta_diff = dabs( delta_dist(ipeak) - phidisp)
c                 --- found if they are equal, or exactly 2Pi apart
                  if ( (delta_diff .lt. delta_dist(0)) .or.
     &                 (dabs(delta_diff-twopi) .lt. delta_dist(0)) 
     &                 ) then
                     lfound = .true.
                     delta_old = ipeak
                  endif
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'GETCBDIHED: problem finding old'
                  lfailure = .true.
                  write(6,*) 'phidisp:',phidisp*360.0d0/twopi
                  do ipeak = 1,npeak
                     write(6,*) 'ipeak,delta_dist:',ipeak
     &                    ,delta_dist(ipeak)*360.0d0/twopi
                  enddo
                  return
               endif
            else
c              --- assign the delta distribution value for this trial
c              --- there are the same number of value as there are 
c              --- torsion trials
               if ( itor .eq. delta_old ) then
c                 --- need to swap values with the first entry on this
c                 --- one as we sampled this value for the old entry
                  phidisp = delta_dist(1)
               else
                  phidisp = delta_dist(itor)
               endif
            endif
            pdf_bias = 1.0d0/dble(npeak)
         else
            write(6,*) 'GETCBDIHED: unknown style'
            lfailure = .true.
            return
         endif

c        --- store the ratio of the pdf
         bf_tor(itor) = pdf_true/pdf_bias
#if DEBUG_GETCBDIHED
         write(6,*) 'phidisp,phishift:',phidisp,phishift
         write(6,*) 'pdf_true,pdf_bias:',pdf_true,pdf_bias
#endif
c        --- initialize lovrtor(itor)
         lovrtor(itor) = .false.
c        --- apply the phishift
         phidisp = phidisp - phishift

c        --- store the new positions
         call twh_coordinates(GLB_GET,testcord,ichain,iufrom
     &        ,rxfrom,ryfrom,rzfrom)
         do count = angstart,ntogrow
            call twh_cbphi(GLB_GET,count,cbphi)
            phinew = cbphi + phidisp
            call twh_bendang(GLB_GET,count,thetanew)
            call twh_cone(lfailure,2,dum,dum,dum,thetanew,phinew,x,y,z)
            if ( lfailure ) return
c           --- store the cartesian postions
            call twh_bondlen(GLB_GET,count,length)
            xcoord(itor,count) = rxfrom + x*length
            ycoord(itor,count) = ryfrom + y*length
            zcoord(itor,count) = rzfrom + z*length
         enddo

         if ( ldihed_direct .or. ldihed_indirect ) then
c           --- compute the energies
c           --- compute the bond vibration energies for growbond
            do count = 1,ntogrow
               call twh_growlist(GLB_GET,iw,count,iu)
               call twh_grownbon(GLB_GET,iw,count,grownbon)
c              --- fetch the coords
               xpos = xcoord(itor,count)
               ypos = ycoord(itor,count)
               zpos = zcoord(itor,count)
               do ixtra = 1,grownbon
                  call twh_growbond(GLB_GET,iw,count,ixtra,jut2)
                  call twh_coordinates(GLB_GET,testcord,ichain,jut2
     &                 ,rxjut,ryjut,rzjut)
                  call twh_mimage(ibox,xpos,ypos,zpos,rxjut,ryjut,rzjut
     &                 ,rxuij,ryuij,rzuij)
                  length = twh_distance(rxuij,ryuij,rzuij)
c                 --- find the vibration type
                  do ivib = 1,twh_get_invib(imolty,iu)
                     if ( twh_get_ijvib(imolty,iu,ivib) .eq. jut2 ) then
                        call twh_itvib(GLB_GET,imolty,iu,ivib,itype)
                        goto 66
                     endif
                  enddo
c                 --- we can only get here if there is a problem
                  write(6,*) 'GETCBDIHED: type not found'
                  lfailure = .true.
                  return
 66               vxvib(itor) = vxvib(itor) + twh_vbond(lfailure,itype
     &                 ,length,imolty,iu,jut2,lovrlap)
                  if ( lfailure ) return
                  if ( lovrlap ) goto 90
               enddo
            enddo

c           --- compute the xtra and growbond angle energies
            do count = 1,ntogrow
               call twh_nxangle(GLB_GET,count,nxangle)
               do ibend = 1,nxangle
                  do ibead = 1,3
                     call twh_xangle(GLB_GET,ibead,count,ibend,itest)
                     if ( itest .lt. 0 ) then
                        iflip = -itest
                        xatmp(ibead) = xcoord(itor,iflip)
                        yatmp(ibead) = ycoord(itor,iflip)
                        zatmp(ibead) = zcoord(itor,iflip)
                        call twh_growlist(GLB_GET,iw,iflip,iutmp(ibead))
                     else
                        call twh_coordinates(GLB_GET,testcord,ichain
     &                       ,itest,xatmp(ibead),yatmp(ibead)
     &                       ,zatmp(ibead))
                        iutmp(ibead) = itest
                     endif
                  enddo
                  call twh_txangle(GLB_GET,count,ibend,itype)
                  vxbend(itor) = vxbend(itor) + twh_vangle(lfailure,
     &                 lovrlap,ibox,imolty,iutmp(1),itype,1
     &                 ,anglec,dumb,dumc,xatmp,yatmp,zatmp)
                  if ( lfailure ) return
                  if ( lovrlap ) goto 90
               enddo
            enddo
c           --- compute the xtra angle-angle energies
            call twh_naabead(GLB_GET,naabead)
            do iaa = 1,naabead
c              --- extract positional information
               do ibead = 1,4
                  call twh_aabead(GLB_GET,iaa,ibead,itest)
                  if ( itest .lt. 0 ) then
                     iflip = -itest
                     xttmp(ibead) = xcoord(itor,iflip)
                     yttmp(ibead) = ycoord(itor,iflip)
                     zttmp(ibead) = zcoord(itor,iflip)
                     if ( ibead .eq. 1 ) then
                        call twh_growlist(GLB_GET,iw,iflip,ttype)
                     endif
                  else
                     call twh_coordinates(GLB_GET,testcord,ichain,itest
     &                    ,xttmp(ibead),yttmp(ibead),zttmp(ibead))
                     if ( ibead .eq. 1 ) ttype = itest
                  endif
               enddo
c              --- assign angle-angle index
               call twh_aaind(GLB_GET,iaa,it)
               vxaa(itor) = vxaa(itor) 
     &              + twh_vangang(AA_COORD,ibox,imolty,ttype,it,dumangle
     &              ,xttmp,yttmp,zttmp)
            enddo

c           --- set energies of trial position to zero ---
            vdha = 0.0d0
            vimprop = 0.0d0
c           --- compute torsion energy for trial conformation
            do count = 1,ntogrow
               call twh_grownbon(GLB_GET,iw,count,grownbon)
               do iend = -1,grownbon
                  if ( iend .eq. -1 ) then
c                    --- compute the torsions which have iu as the
c                    --- end and don't include any other growlist
c                    --- atoms
                     call twh_growlist(GLB_GET,iw,count,iu)
                  elseif ( iend .eq. 0 ) then
c                    --- compute the torsions which have iufrom as
c                    --- the end bead, growlist as the second, and
c                    --- growbond as the third
                     iu = iufrom
                  else
c                    --- compute the torsions which have growbond
c                    --- as the first, count as the second, iufrom
c                    --- as the third and either growlist, iuprev
c                    --- or ixtra as the 4th
                     call twh_growbond(GLB_GET,iw,count,iend,iu)
                  endif

c                 --- check all of the torsions starting at iu
                  do it = 1, twh_get_intor(imolty,iu)
c                    --- initialize lcompute
                     lcompute = .false.
c                    --- assign the unit types for the torsion
                     torbead(1) = iu
                     torbead(2) = twh_get_ijtor2(imolty,iu,it)
                     torbead(3) = twh_get_ijtor3(imolty,iu,it)
                     torbead(4) = twh_get_ijtor4(imolty,iu,it)
c                    --- determine if this torsion should be computed
                     if ( iend .eq. -1 ) then
                        if ( twh_get_lexist(torbead(2)) 
     &                       .and. twh_get_lexist(torbead(3))
     &                       .and. twh_get_lexist(torbead(4)) ) then
c                          --- compute this torsion
                           lcompute =.true.
c                          --- temporarily switch back to count
                           torbead(1) = -count
                        endif
                     elseif ( iend .eq. 0 ) then
                        call twh_growlist(GLB_GET,iw,count,growtest)
                        if ( growtest .eq. torbead(2)
     &                       .and. twh_get_lexist(torbead(3)) 
     &                       .and. twh_get_lexist(torbead(4)) ) then
c                          --- compute this torsion
                           lcompute =.true.
c                          --- the second bead is count
                           torbead(2) = -count
                        endif
                     else
                        call twh_growlist(GLB_GET,iw,count,growtest)
                        if ( growtest .eq. torbead(2)
     &                       .and. twh_get_lexist(torbead(3)) ) then
c                          --- see if torbead(4) either exists, or
c                          --- is one of the other growlist atoms
                           if ( twh_get_lexist(torbead(4)) ) then
c                             --- compute this torsion
                              lcompute = .true.
c                             --- torbead 2 is count
                              torbead(2) = -count
                           else
c                             --- search through the growlist
                              do itest = 1,ntogrow
                                 call twh_growlist(GLB_GET,iw,itest
     &                                ,growtest)
                                 if ( growtest .eq. torbead(4) ) then
c                                   --- compute this torsion
                                    lcompute = .true.
c                                   --- torbead 2,4 are growlist
                                    torbead(2) = -count
                                    torbead(4) = -itest
                                 endif
                              enddo
                           endif
                        endif   
                     endif
                     
                     if ( lcompute ) then
c                       --- assign the positions
                        do ibead = 1,4
                           itest = torbead(ibead)
                           if ( itest .lt. 0 ) then
c                             --- use the count coordinates
                              iflip = -itest
                              xttmp(ibead) = xcoord(itor,iflip)
                              yttmp(ibead) = ycoord(itor,iflip)
                              zttmp(ibead) = zcoord(itor,iflip)
c                             --- reset torbead
                              call twh_growlist(GLB_GET,iw,iflip
     &                             ,torbead(ibead))
                           else
c                             --- use old/new unit coordinates
                              call twh_coordinates(GLB_GET,testcord
     &                             ,ichain,itest,xttmp(ibead)
     &                             ,yttmp(ibead),zttmp(ibead))
                           endif
                        enddo
c                       --- store unit numbers
                        ttype = torbead(1)
c                       --- add torsion energy to vdha
                        vdha = vdha + twh_vtorsion(lfailure,ibox,1
     &                       ,imolty,ttype,it,xttmp,yttmp,zttmp,dum
     &                       ,lovrlap)
                        if ( lfailure ) return
                        if ( lovrlap ) goto 90
                     endif
                  enddo
               enddo
            enddo

c           --- compute the improper torsion energy and add it to vdha
            call twh_nimpbead(GLB_GET,nimpbead)
            do iimprop = 1,nimpbead
c              --- extract positional information
               do ibead = 1,4
                  call twh_impbead(GLB_GET,iimprop,ibead,itest)
                  if ( itest .lt. 0 ) then
                     iflip = -itest
                     xttmp(ibead) = xcoord(itor,iflip)
                     yttmp(ibead) = ycoord(itor,iflip)
                     zttmp(ibead) = zcoord(itor,iflip)
                     if ( ibead .eq. 1 ) then
                        call twh_growlist(GLB_GET,iw,iflip,ttype)
                     endif
                  else
                     call twh_coordinates(GLB_GET,testcord,ichain,itest
     &                    ,xttmp(ibead),yttmp(ibead),zttmp(ibead))
                     if ( ibead .eq. 1 ) ttype = itest
                  endif
               enddo
c              --- assign improper torsion type
               call twh_impropernum(GLB_GET,iimprop,it)
               vimprop = vimprop + twh_vimproper(ibox,imolty
     &              ,ttype,it,xttmp,yttmp,zttmp,dum,lbadstereo)
            enddo
            
c           --- compute any additional biasing energies which are 
c           --- used to help with fecbmc regrowth
            do count = 1,ntogrow
               call twh_nxbias(GLB_GET,count,nxbias)
               do ibias = 1,nxbias
                  call twh_growlist(GLB_GET,iw,count,iunit)
                  call twh_xbias(GLB_GET,1,count,ibias,biastype)
                  call twh_xbias(GLB_GET,2,count,ibias,index)
                  if ( biastype .eq. 1 ) then
c                    --- angle biasing
                     iuother = twh_get_ijben3(imolty,iunit,index)
                  elseif ( biastype .eq. 2 ) then
c                    --- torsional biasing
                     iuother = twh_get_ijtor4(imolty,iunit,index)
                  else
                     write(6,*) 'GETCBDIHED: unknown biastype',biastype
                     lfailure = .true.
                     return
                  endif
c                 --- store the cordinates in xyzcord 
c                 --- iunit is in slot 1
                  xcord(1) = xcoord(itor,count)
                  ycord(1) = ycoord(itor,count)
                  zcord(1) = zcoord(itor,count)
c                 --- iuother in slot 2
                  call twh_coordinates(GLB_GET,testcord,ichain,iuother
     &                 ,xcord(2),ycord(2),zcord(2))
c                 --- iufrom in slot 3
                  call twh_coordinates(GLB_GET,testcord,ichain,iufrom
     &                 ,xcord(3),ycord(3),zcord(3))
c                 --- call bias function
                  bf_bias(itor) = bf_bias(itor)*
     &                 twh_febias(lfailure,ibox,biastype,imolty,iunit
     &                 ,index,iufrom,xcord,ycord,zcord)
                  if ( lfailure ) return
               enddo
            enddo

c           --- store vtors for this trial
            vtors(itor) = vdha + vimprop

c           --- compute boltzmann factor and add it to w_cbtor
 90         if ( lovrlap ) then
               lovrtor(itor) = .true.
            else
               if ( lmin ) then
c                 --- see if this energy is lower than current
                  if ( vtors(itor)+vxbend(itor)+vxvib(itor)
     &                 +vxaa(itor) .lt. vmin ) then
                     vmin  = vtors(itor)+vxbend(itor)+vxvib(itor)
     &                    + vxaa(itor)
                     imin = itor
                  endif
               else
c                 --- this is our first nonoverlap energy set vmin
                  lmin = .true.
                  vmin = vtors(itor)+vxbend(itor)+vxvib(itor)+vxaa(itor)
                  imin = itor
               endif
            endif
         endif
      enddo

c     --- Begin loop to select the torsional angle
c     --- initialize w_cbtor
      w_cbtor = 0.0d0
      lfound = .false.
      do itor = 1,ichtor
         if ( lovrtor(itor) ) then
            bf_tor(itor) = 0.0d0
         else
            lfound = .true.
            bf_tor(itor) = bf_tor(itor)*bf_bias(itor)
     &           *twh_expon( -beta*
     &           (vtors(itor)+vxbend(itor)+vxaa(itor)+vxvib(itor)
     &           -vmin))
         endif
         w_cbtor = w_cbtor + bf_tor(itor)
#if DEBUG_GETCBDIHED
         write(6,*) 'itor,bf_tor,w_cbtor:',itor,bf_tor(itor),w_cbtor
#endif
      enddo

#if DEBUG_GETCBDIHED
      write(6,*) 'luseoldone:',luseoldone,' lfound:',lfound
     &,' lterm:',lterm
      write(6,*) 'w_cbtor:',w_cbtor
#endif
      if ( luseoldone ) then
c        --- make sure there is not a hard overlap
         if ( .not. lfound ) then
c           --- this should never happen for the old conformation
            write(6,*) 'GETCBDIHED: none found for old'
            lfailure = .true.
            return
         endif
c        --- select the actual old conformation
         iwinner = 1
      else
c        --- new growth
         if ( .not. lfound ) then
c           --- all sites had hard overlap
            lterm = .true.
         elseif ( ichtor .eq. 1 ) then
c           --- no need for a random number as there is only 1 option
            iwinner = 1
         else
c           --- choose one of the trial sites in a biased fashion
            ran_tor = twh_random()*w_cbtor
            bs = 0.0d0
            iwinner = 0
            do while ( ran_tor .ge. bs .and. iwinner .le. ichtor ) 
               iwinner = iwinner + 1
               bs = bs + bf_tor(iwinner)
            enddo
c           --- safety check
            if ( iwinner .gt. ichtor ) then
c              --- none of the sites had a non-zero boltzmann weight
               lterm = .true.
            endif
         endif

         if ( lterm ) then
            if ( lmin ) then
               write(6,*) 'GETCBDIHED: all zero boltzmann'
               write(6,*) 'lmin,vmin',lmin,vmin
               write(6,*) 'expon(-beta vmin)',twh_expon(-beta*vmin)
               write(6,*) 'w_cbtor:',w_cbtor
               write(6,*) 'ldihed_direct,ldihed_indirect',ldihed_direct
     &              ,ldihed_indirect
               write(6,*) 'should not happen when lmin is true'
               write(6,*) 'lfound ',lfound
               write(6,*) 'bs:',bs
               write(6,*) 'ran_tor:',ran_tor
               write(6,*) 'ichtor:',ichtor
               do itor = 1,ichtor
                  write(6,*) 'itor,lovrtor,bf_tor',itor,lovrtor(itor)
     &                 ,bf_tor(itor)
                  write(6,*) 'bf_bias(itor)',bf_bias(itor)
                  write(6,*) 'expon(stuff)',
     &                 twh_expon( -beta*
     &                 (vtors(itor)+vxbend(itor)+vxaa(itor)+vxvib(itor)
     &                 -vmin))
                  write(6,*) 'vtors:',vtors(itor)
                  write(6,*) 'vxbend:',vxbend(itor)
                  write(6,*) 'vxaa:',vxaa(itor)
                  write(6,*) 'vxvib(itor):',vxvib(itor)
                  write(6,*) 'vmin:',vmin
               enddo
               lfailure = .true.
               return
            endif
            return
         endif
      endif

c     --- done with phidisp selection 
c     --- accepted coordinates, save them in *torcord
      do count = 1,ntogrow
         xtorcord(count) = xcoord(iwinner,count)
         ytorcord(count) = ycoord(iwinner,count)
         ztorcord(count) = zcoord(iwinner,count)
      enddo

c     --- divide bsum by ichtor
      w_cbtor = w_cbtor / dble(ichtor)
c     --- rescale the boltzman weight to account for vmin
      call twh_getweight(lfailure,0,vmin,wmin,dpcbtor,lterm)
      if ( lfailure ) return
      if ( lterm ) then
c        --- extremely high energy compared to kT - growth failed
         w_cbtor = 0.0d0
         if ( luseoldone ) then
c           --- really should not have termination for old
            write(6,*) 'GETCBDIHED: old termination'
            write(6,*) 'vmin,wmin,dpcbtor,lterm',vmin,wmin,dpcbtor
     &           ,lterm
         endif
      else
c        --- rescasle the rosenbluth weight
         w_cbtor = w_cbtor*wmin
c        --- save energes of this trial position
         vcbvib = vxvib(iwinner)
         vcbbend = vxbend(iwinner)+vxaa(iwinner)
         vcbtor = vtors(iwinner)
c        --- save bias bolztman factor of this trial position
         bfcbbias = bf_bias(iwinner)
      endif

#if DEBUG_GETCBDIHED
      write(6,*) 'finish GETCBDIHED:'
     &     ,'ldihed_direct: ',ldihed_direct
     &     ,'lidhed_indirect: ',ldihed_indirect
     &     ,'luseoldone: ',luseoldone
     &     ,'lterm: ',lterm
     &     ,'testcord: ',testcord
     &     ,'iw: ',iw
     &     ,'ichain: ',ichain
     &     ,'ichtor: ',ichtor
     &     ,'iufrom: ',iufrom
     &     ,'ibox: ',ibox
     &     ,'angstart: ',angstart
     &     ,'ntogrow: ',ntogrow
     &     ,'imolty: ',imolty
     &     ,'dpcbtor: ',dpcbtor
     &     ,'phishift: ',phishift
     &     ,'w_cbtor: ',w_cbtor
     &     ,'vcbvib: ',vcbvib
     &     ,'vcbbend: ',vcbbend
     &     ,'vcbtor: ',vcbtor
     &     ,'bfcbbias: ',bfcbbias
     &     ,'xtorcord: ',(xtorcord(count),count=1,ntogrow)
     &     ,'ytorcord: ',(ytorcord(count),count=1,ntogrow)
     &     ,'ztorcord: ',(ztorcord(count),count=1,ntogrow)
     &     ,'torkey: ',torkey
#endif
      return
      end
