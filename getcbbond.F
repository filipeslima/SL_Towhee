#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2014 Marcus G. Martin                       *
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
      subroutine twh_getcbbond(lfailure,luseoldone,testcord,iw,ichain
     &     ,imolty,iuprev,nch_vib,vvibtr,maxlen,w_cbbond,dpcbbond,ibox)
c     ******************************************************************
c     * determines the bond lengths in a biased fashion                *
c     *                                                                *
c     * luseoldone: T if first selection uses the old coordinate       *
c     *             F if all selections generated from scratch         *
c     * testcord: flag for fetching existing coordinates               *
c     * iw: number for the walk - used with the variables that are set *
c     *     in schedule.F                                              *
c     * ichain: chain number                                           *
c     * imolty: molecule type number                                   *
c     * iuprev: unit number for that atom that is the previous bead in *
c     *         the cbmc move                                          *
c     * vvibtr: bond (including bond-bond) energy for the selected     *
c     *         coordinates (returned from this subroutine)            *
c     * maxlen: maximum bond length of the accepted atom coords.  This *
c     *         NOTE: maxlen is not initialized in this subroutine as  *
c     *         it keeps track of the maxlen for any of the bonding    *
c     *         selections of a single walk                            *
c     * w_cbbond: rosenbluth weight (normalized) for the bond selection*
c     *           (returned by this subroutine)                        *
c     * dpcbbond: scaling factor for number of time the weight         *
c     *          had to be scaled (returned by this subroutine         *
c     * ibox: the box in which the growth is performed                 *
c     *                                                                *
c     * rewritten from geometry 07-26-2002 by M.G. Martin              *
c     * last modified 07-16-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_EXPON
#define FUNCTION_GAUSSIAN
#define FUNCTION_GAUSSPROB
#define FUNCTION_GET_IJBNBN1
#define FUNCTION_GET_IJBNBN2
#define FUNCTION_GET_IJVIB
#define FUNCTION_GET_INBNBN
#define FUNCTION_GET_INVIB
#define FUNCTION_GET_ITBNBN
#define FUNCTION_GET_VIBCOEFF
#define FUNCTION_MINBOXLENGTH
#define FUNCTION_RANDOM
#define FUNCTION_VBOND
#define FUNCTION_VBONBON
#include "functions.h"
#define DEBUG_GETCBBOND 0
c     --- variables passed to/from the subroutine
      logical lfailure,luseoldone
      integer iw,imolty,iuprev,testcord,ibox,ichain,nch_vib
      double precision vvibtr,maxlen,dpcbbond
      double precision w_cbbond
c     --- local variables
c     --- logical scalars
      logical lovrlap,lmin,lterm
c     --- logical arrays
      logical lovrvib
      dimension lovrvib(nch_vib)
c     --- integer scalars
      integer count,ntogrow,iugrow,iufrom,iv,juvib,jtvib,ib
      integer ivib,iwinner,iut,ixtra,inum,ic,it,nchoice
      integer vibnumber,growprev,xtracount,dist_code
c     --- integer arrays
      integer bblocl
      dimension bblocl(2,MAXBEND)
      integer bbtot
      dimension bbtot(NNBOND)
      integer bbtype
      dimension bbtype(NNBOND,MAXBEND)
c     --- double precision scalars
      double precision length,dpmin,beta,distance,bsum_try,rbf,bs
      double precision rxfrom,ryfrom,rzfrom,rxprev,ryprev,rzprev
      double precision rxtry,rytry,rztry,rxgrow,rygrow,rzgrow
      double precision wmin,vmin,rxuij,ryuij,rzuij,minboxdim
      double precision pdf_true,pdf_bias
      double precision xvecprev,yvecprev,zvecprev
      double precision mean,std_dev,lolimit,hilimit,onethird
      double precision arbnorm,truenorm
c     --- double precision arrays
      double precision vvib
      dimension vvib(nch_vib)
      double precision blfact,vib_trial
      dimension blfact(nch_vib),vib_trial(nch_vib)
      double precision bblen
      dimension bblen(3)

c      write(6,*) 'START GETCBBOND'
c     --- retrieve constants
      call twh_constant_onethird(GLB_GET,onethird)
      call twh_minboxdim(GLB_GET,minboxdim)

c     --- initialize the local lterm which is used in getweight
      lterm = .false.
c     --- initialize trial energies
      vvibtr = 0.0d0
c     --- initialize w_cbbond and dpcbbond
      w_cbbond = 1.0d0
      dpcbbond = 0.0d0
c     --- assign grownum and growfrom to local variables
      call twh_grownum(GLB_GET,iw,ntogrow)
      call twh_growfrom(GLB_GET,iw,iufrom)
c     ************************************************************
c     ** Begin Bond length biased selection                     **
c     ************************************************************

c     --- store r*ui positions for unit iufrom
      call twh_coordinates(GLB_GET,testcord,ichain,iufrom,rxfrom
     &     ,ryfrom,rzfrom)

      call twh_growprev(GLB_GET,iw,growprev)
      if ( growprev .ne. 0 ) then
c        --- compute the vector from iufrom to iuprev
         call twh_coordinates(GLB_GET,testcord,ichain,iuprev
     &        ,rxprev,ryprev,rzprev)
         call twh_mimage ( ibox,rxprev,ryprev,rzprev
     &        ,rxfrom,ryfrom,rzfrom,xvecprev,yvecprev,zvecprev)
c        --- store the vecprev for use in getcbangle
         call twh_vecprev(GLB_SET,xvecprev,yvecprev,zvecprev)
         distance = twh_distance(xvecprev,yvecprev,zvecprev)
         call twh_distprev(GLB_SET,1,distance)
      endif

c     --- initialize bond-bond list
      do count = 1,ntogrow
         bbtot(count) = 0
      enddo

c     --- check for IJ-JK bond-bond and set up list
      do ic = 1,twh_get_inbnbn(imolty,iufrom)
         do ib = 1,2
            bblocl(ib,ic) = 0
         enddo
c        --- determine which (if any) of the slots is growprev
         if ( growprev .ne. 0 ) then 
            if (twh_get_ijbnbn1(imolty,iufrom,ic) .eq. growprev 
     &           ) then
               bblocl(1,ic) = -1
            elseif (twh_get_ijbnbn2(imolty,iufrom,ic).eq. growprev
     &              ) then
               bblocl(2,ic) = -1 
            endif
         endif
c        --- account for any growxtra atoms in the bond-bond
         call twh_growxtra(GLB_GET,iw,1,xtracount)
         do count = 2,xtracount+1
            if (twh_get_ijbnbn1(imolty,iufrom,ic) .eq. xtracount ) then
               bblocl(1,ic) = -count
            elseif (twh_get_ijbnbn2(imolty,iufrom,ic) 
     &              .eq. xtracount ) then
               bblocl(2,ic) = -count 
            endif
         enddo
c        --- set up the logic for bond-bond
         do count = 1,ntogrow
            call twh_growlist(GLB_GET,iw,count,iugrow)
            if (twh_get_ijbnbn1(imolty,iufrom,ic) .eq. iugrow) then
               bblocl(1,ic) = count 
            elseif (twh_get_ijbnbn2(imolty,iufrom,ic) .eq. iugrow) then 
               bblocl(2,ic) = count 
            endif
            if ( (bblocl(1,ic)*bblocl(2,ic)) .ne. 0 ) then
c              --- we now have enough atoms to check the bond-bond
               bbtot(count) = bbtot(count) + 1
               bbtype(count,bbtot(count)) = ic
c              --- get out of loop so we don't count this again
               goto 14
            endif
         enddo
         write(6,*) 'GETCBBOND: bond-bond logic messed up'
         lfailure = .true.
         return
 14      continue
      enddo

c     --- set up xtra distances for use in bond-bond and in rosenbluth
      call twh_growxtra(GLB_GET,iw,1,xtracount)
      if ( xtracount .ne. 0 ) then
         do ixtra = 2,xtracount+1
            call twh_growxtra(GLB_GET,iw,ixtra,iut)
            call twh_coordinates(GLB_GET,testcord,ichain,iut
     &           ,rxtry,rytry,rztry)
            call twh_mimage( ibox,rxtry,rytry,rztry
     &           ,rxfrom,ryfrom,rzfrom,rxuij,ryuij,rzuij)
            distance = twh_distance(rxuij,ryuij,rzuij)
            call twh_distprev(GLB_SET,ixtra,distance)
         enddo
      endif

c     --- determine individual bond lengths
      do count = 1,ntogrow
c        --- initialize lmin which tells us if we have minimum yet
         lmin = .false.
         call twh_growlist(GLB_GET,iw,count,iugrow)
         call twh_glist(GLB_SET,count,iugrow)
            
c        --- determine the vibration (bond) type as jtvib
         do iv = 1, twh_get_invib(imolty,iugrow)
            juvib = twh_get_ijvib(imolty,iugrow,iv)
            if ( juvib .eq. iufrom ) then
               call twh_itvib(GLB_GET,imolty,iugrow,iv,jtvib)
               vibnumber = iv
            endif
         enddo

c        --- set default number of choices
         nchoice = nch_vib
         call twh_cbbonddist(GLB_GET,imolty,iugrow,vibnumber,dist_code)
         if ( dist_code .eq. DIST_DELTA ) then
c           --- one rigid bond length
            nchoice = 1
         elseif ( dist_code .eq. DIST_GAUSSIAN ) then
c           --- set the limits of the gaussian bond generation
            lolimit = 0.0d0
            call twh_max_bond_length(GLB_GET,hilimit)
            call twh_cbmcbondinfo(GLB_GET,imolty,iugrow,vibnumber
     &           ,3,mean)
            call twh_cbmcbondinfo(GLB_GET,imolty,iugrow,vibnumber
     &           ,4,std_dev)
         elseif ( dist_code .eq. DIST_R_SQ ) then
            call twh_cbmcbondinfo(GLB_GET,imolty,iugrow,vibnumber,3
     &           ,lolimit)
            call twh_cbmcbondinfo(GLB_GET,imolty,iugrow,vibnumber,4
     &           ,hilimit)
         else
c           --- unknown dist_code
            write(6,*) 'GETCBBOND: unknown dist_code'
            write(6,*) dist_code
            write(6,*) 'imolty,iugrow,iufrom,jtvib:'
     &           ,imolty,iugrow,iufrom,jtvib
            lfailure = .true.
            return
         endif
c        --- retrive the arbitrary and true normalization factors
         call twh_cbmcbondinfo(GLB_GET,imolty,iugrow,vibnumber,1
     &        ,truenorm)
         call twh_cbmcbondinfo(GLB_GET,imolty,iugrow,vibnumber,2
     &        ,arbnorm)

         do ivib = 1,nchoice
c           --- initialize lovrlap
            lovrlap = .false.
            if ( ivib .eq. 1 .and. luseoldone ) then
c              --- first ivib is old bond length
               call twh_coordinates(GLB_GET,testcord,ichain,iugrow
     &              ,rxgrow,rygrow,rzgrow)
               call twh_mimage(ibox,rxgrow,rygrow,rzgrow
     &              ,rxfrom,ryfrom,rzfrom,rxuij,ryuij,rzuij)
               length = twh_distance(rxuij,ryuij,rzuij)
               if ( length .gt. 0.5d0*twh_minboxlength(ibox) ) then
                  write(6,*) 'GETCBBOND: old bond length '
                  write(6,*) 'exceeds 0.5*minboxlength'
                  write(6,*) 'ichain:',ichain
                  write(6,*) 'rxgrow,rygrow,rzgrow:'
     &                 ,rxgrow,rygrow,rzgrow
                  write(6,*) 'rxfrom,ryfrom,rzfrom:'
     &                 ,rxfrom,ryfrom,rzfrom
                  write(6,*) 'rxuij,ryuij,rzuij:',rxuij,ryuij,rzuij
                  write(6,*) 'imolty,iugrow,jtvib',imolty,iugrow,jtvib
                  write(6,*) 'length ',length
                  lfailure = .true.
                  return
               endif
            else
               if ( dist_code .eq. DIST_DELTA ) then
c                 --- fixed bond lengths.
                  call twh_cbmcbondinfo(GLB_GET,imolty,iugrow
     &                 ,vibnumber,3,length)
               elseif ( dist_code .eq. DIST_GAUSSIAN ) then
c                 --- sample from a gaussian distribution
                  length = twh_gaussian(lfailure,mean,std_dev
     &                 ,lolimit,hilimit)
                  if ( lfailure ) return
               elseif ( dist_code .eq. DIST_R_SQ ) then
c                 --- sample according to a bounded r^2 distribution
                  length = ((lolimit+twh_random()*hilimit))**onethird
               else
                  write(6,*) 'GETCBBOND: unknown dist_code'
                  write(6,*) dist_code
                  lfailure = .true.
                  return
               endif
            endif

c           --- set pdf_bias
            if ( dist_code .eq. DIST_GAUSSIAN ) then
c              --- sampling from a gaussian distribution
               pdf_bias = twh_gaussprob(length,mean,std_dev,lolimit
     &              ,hilimit)
            elseif ( dist_code .eq. DIST_DELTA ) then
c              --- rigid bond, constant
               pdf_bias = arbnorm
            elseif ( dist_code .eq. DIST_R_SQ ) then
c              --- sampling from ideal r^2 distribution
               pdf_bias = length*length*arbnorm
            else
               write(6,*) 'GETCBBOND: unknown dist_code'
               write(6,*) dist_code
               lfailure = .true.
               return
            endif
c           --- store the trial length
            vib_trial(ivib) = length
c           --- store the ratio of the trial pdf factors
            if ( pdf_bias .le. 0.0d0 ) then
c              --- this is an extremely unlikely bond, print a warning
               write(6,*) 'GETCBBOND: pdf_bias probability too low'
               write(6,*) 'dist_code:',dist_code
               write(6,*) 'luseoldone',luseoldone
               write(6,*) 'ivib,pdf_bias',ivib,pdf_bias
               write(6,*) 'length,mean,std_dev',length,mean,std_dev
               write(6,*) 'equilibrium length: '
     &              ,twh_get_vibcoeff(jtvib,0)
               write(6,*) 'imolty,iugrow,vibnumber',imolty,iugrow
     &              ,vibnumber
               write(6,*) 'ijvib',twh_get_ijvib(imolty,iugrow,vibnumber)
               write(6,*) 'unless this is near the start of a run'
     &              ,' you have a big problem'
               write(6,*) 'resetting pdf_bias to 1.0d20 in the hope of '
     &              ,'accepting a move away from this bond length'
               pdf_bias = 1.0d20
            endif
            if ( dist_code .eq. DIST_DELTA ) then
c              --- finite number of non-infinite energy lengths
               pdf_true = truenorm
            else
c              --- continuous distribution is proportional to r squared
               pdf_true = length*length*truenorm
            endif
            blfact(ivib) = pdf_true/pdf_bias
c           --- calculate the bond vibration energy
            vvib(ivib) = twh_vbond(lfailure,jtvib,length,imolty,iufrom
     &           ,iugrow,lovrlap)
            if ( lfailure ) return
c           --- calculate the IJ-JK bond-bond energy
            do ib = 1,bbtot(count)
               ic = bbtype(count,ib)
               do inum = 1,2
                  if ( bblocl(inum,ic) .eq. count ) then
                     bblen(inum) = length
                  elseif ( bblocl(inum,ic) .lt. 0 ) then
                     iut = -bblocl(inum,ic)
                     call twh_distprev(GLB_GET,iut,distance)
                     bblen(inum) = distance
                  else
                     call twh_bondlen(GLB_GET,bblocl(inum,ic),distance)
                     bblen(inum) = distance
                  endif
               enddo
               it = twh_get_itbnbn(imolty,iufrom,ic)
               vvib(ivib) = vvib(ivib) 
     &              + twh_vbonbon(imolty,iufrom,ic,bblen)
            enddo

            if ( lovrlap ) then
               lovrvib(ivib) = .true.
            else
               lovrvib(ivib) = .false.
               if ( lmin ) then
c                 --- see if this energy is lower than our current 
c                 --- minimum
                  if ( vvib(ivib) .lt. vmin ) vmin = vvib(ivib)
               else
c                 --- this is the first energy that is not an overlap
c                 --- set vmin to this energy and change the flag
                  lmin = .true.
                  vmin = vvib(ivib)
               endif
            endif
         enddo

         if ( .not. lmin ) then
c           --- all vibrations were overlaps, this is odd
            write(6,*) 'all attempted vibrations were overlaps'
            lfailure = .true.
            return
         endif

c        --- initialize bsum_try
         bsum_try = 0.0d0
c        --- compute the rosenbluth weight after subtracting the 
c        --- minimum energy
         do ivib = 1,nchoice
#if DEBUG_GETCBBOND
      write(6,*) 'ivib,ntogrow,iufrom,iugrow,vib_trial,vvib',ivib
     &,ntogrow,iufrom,iugrow
     &,vib_trial(ivib),vvib(ivib)
#endif
            if ( lovrvib(ivib) ) then
               blfact(ivib) = 0.0d0
            else
               call twh_beta(GLB_GET,beta)
               blfact(ivib)= blfact(ivib)
     &              *twh_expon(-beta*(vvib(ivib)-vmin))
            endif
            bsum_try = bsum_try + blfact(ivib)
         enddo
         
         if ( .not. luseoldone ) then
c           --- select one of the bond lengths in a biased fashion
            if ( nchoice .eq. 1 ) then
               if ( .not. lovrvib(1) ) then
c                 --- no need for random number, accept our only option
                  iwinner = 1
                  goto 5
               endif
            else
               rbf = twh_random()*bsum_try
               bs = 0.0d0
               do ivib = 1,nchoice
                  bs = bs + blfact(ivib)
                  if ( rbf .lt. bs ) then
                     iwinner = ivib
                     goto 5
                  endif
               enddo
            endif
c           --- if we got here then all of the bond lengths must have
c           --- been hard overlaps.  This is unusual so print a warning
            write(6,*) 'No good bond lengths in in bond selection ',bs
            lfailure = .true.
            return
 5          continue
         else
c           --- select the old bond length
            iwinner = 1
         endif

c        --- account for vmin
         call twh_getweight(lfailure,0,vmin,wmin,dpmin,lterm)
         if ( lfailure ) return
         if ( lterm ) then
c           --- this really should never happen as the bond length
c           --- energies should be pretty well behaved. 
            write(6,*) 'GETCBBOND: getweight returned lterm in bond'
            lfailure = .true.
            return
         endif

c        --- propagate the rosenbluth weight and dpcbbond
         w_cbbond = w_cbbond * wmin * bsum_try / dble(nchoice)
         dpcbbond = dpcbbond + dpmin

c        --- assign bondlength and add up vibrational energy
         distance = vib_trial(iwinner)
         call twh_bondlen(GLB_SET,count,distance)
         vvibtr = vvibtr + vvib(iwinner)
c        --- adjust maximum bond length of those being grown
         if ( distance .gt. maxlen ) maxlen = distance

      enddo

c     ***********************************************************
c     ** End Bond length biased selection                      **
c     ***********************************************************

c      write(6,*) 'FINISH GETCBBOND'
      return
      end
