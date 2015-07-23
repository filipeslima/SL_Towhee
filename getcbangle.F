#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2014 Marcus G. Martin                       *
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
      subroutine twh_getcbangle(lfailure,luseoldone,lovrangle,testcord
     &     ,iw,ichain,imolty,angstart,iuprev,nch_bend_a,nch_bend_b
     &     ,max_nch_bend,gcba_vvib,gcba_vbend,gcba_vtor,gcba_weight
     &     ,dpcbangle,ibox)
c     ******************************************************************
c     * determines the new geometry of the angles in a cbmc growth     *
c     *                                                                *
c     * luseoldone: T use the old coordinates for the first trial in   *
c     *               each of the selections                           *
c     *             F generate new coordaintes for each trial          *
c     * testcord: flag for fetching existing coordinates               *
c     * iw: walk index used for the cbmc variables set in schedule.F   *
c     * ichain: chain index of the target molecule                     *
c     * imolty: molecule type number                                   *
c     * angstart: the number of the grow atoms that forms the first    *
c     *           angle with the prev bead.  Could be 1 or 2 depending *
c     *           on whether the prev bead already existed, or if it is*
c     *           set from one of the grown beads                      *
c     * iuprev: the unit number of the prev bead                       *
c     * gcba_vvib: the vibration energy of the accepted trial          *
c     * gbca_vbend: the bond bending energy of the accepted trial      *
c     * gcba_vtor: the improper and regular torsion energy of the      *
c     *            accepted trial                                      *
c     * gcba_weight: the rosenbluth weight (normalized) from this angle*
c     *            selection                                           *
c     * dpcbangle: the scale factor for the angle weight               *
c     * ibox: the box in which the growth is performed                 *
c     *                                                                *
c     * written by M.G. Martin 07-10-1998 from geomnew and geomold     *
c     * last modified 07-24-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_ARCCOS
#define FUNCTION_DISTANCE
#define FUNCTION_DOTPRODUCT
#define FUNCTION_EXPON
#define FUNCTION_GAUSSIAN
#define FUNCTION_GAUSSPROB
#define FUNCTION_GET_BONDPATT
#define FUNCTION_GET_IJAA0
#define FUNCTION_GET_IJAA1
#define FUNCTION_GET_IJAA2
#define FUNCTION_GET_IJBEN2
#define FUNCTION_GET_IJBEN3
#define FUNCTION_GET_IJIMPROP2
#define FUNCTION_GET_IJIMPROP3
#define FUNCTION_GET_IJIMPROP4
#define FUNCTION_GET_IMPFORM
#define FUNCTION_GET_INAA
#define FUNCTION_GET_INBEN
#define FUNCTION_GET_INIMPROP
#define FUNCTION_GET_NTYPE
#define FUNCTION_GET_ITBEN
#define FUNCTION_GET_ITIMPROP
#define FUNCTION_LAWOFCOSINE
#define FUNCTION_INVERSELAWOFCOSINE
#define FUNCTION_RANDOM
#define FUNCTION_VANGANG
#define FUNCTION_VANGLE
#define FUNCTION_VBOND
#define FUNCTION_VIMPROPER
#define FUNCTION_VTORSION
#include "functions.h"
c     --- debog toggle
#define DEBUG_GETCBANGLE 0
c     --- variables passed to/from the subroutine
      logical luseoldone,lfailure,lovrangle
      integer iw,ichain,imolty,angstart,iuprev,testcord,ibox
      integer nch_bend_a,nch_bend_b,max_nch_bend
      double precision gcba_weight,gcba_vvib,gcba_vbend,gcba_vtor
      double precision dpcbangle
c     --- local variables
c     --- logical scalars
      logical lingeom,lovrlap,lbadstereo,lrigid,lrandelta,lthreemem
c     --- logical arrays
      logical lovrben,lmin,lterm,lfound,lfirstold
      dimension lovrben(max_nch_bend)
c     --- character string scalars
      character*5 hybrid
c     --- integer scalars
      integer count,ntogrow,iugrow,iufrom,iu2back,rimptotal,total
      integer iczero,icone,ictwo,iulast,type,aaa,iuxtra,iangle,nangle
      integer start_ang,index,delta_key
      integer iwinner,ixtra,ic,itemp,itor,step,other,imin,angle_style
      integer iuone,ipeak,npeaks,ivalue,xtracount
      integer testzero,testone,testtwo,testthree,testfour,growprev
      integer cbmc_bendgen_code,style_code,nchoice,iutest
      integer ivib,numvib,vibtype,vibstyle,nbend,ibend
      integer fgp_index,gpf_index,ifind,itrial,itest,numtor
      integer threecode,fourcode,ntmtorsion
      integer cordstyle,angleadist
c     --- integer arrays
      integer imptot
      dimension imptot(NNBOND)
      integer tmpimptype
      dimension tmpimptype(NNBOND,MAXIMPROP)
      integer tmpimpform
      dimension tmpimpform(NNBOND,MAXIMPROP)
c     --- double precision scalars
      double precision thetaone,thetatwo,theta,phione,phitwo,std_dev
      double precision distg,vphi,rdecide,dpmin,ang_gamma
      double precision bsum_try,rbf,bs,onepi,ang_alpha,ang_beta
      double precision rxfrom,ryfrom,rzfrom,rxprev,ryprev,rzprev
      double precision rxgrow,rygrow,rzgrow,beta,angletemp
      double precision xvecgrow,yvecgrow,zvecgrow,distgrow,anglec
      double precision xub,yub,zub,dum,ux,uy,uz,alpha,gamma,dumb
      double precision icr,icphi,ictheta,wmin,vmin,twopi,distance
      double precision xvecprev,yvecprev,zvecprev,bendang,cbphi
      double precision bend_a_ideal_fraction,bend_b_ideal_fraction
      double precision delta_shift,distlo,disthi,lobound,hibound
      double precision testlen,vtotal,distother,distaaa
      double precision fgp_anglec,gpf_anglec,vtest,phic
c     --- double precison arrays
      double precision ang_trial
      dimension ang_trial(max_nch_bend)
      double precision bfactor
      dimension bfactor(max_nch_bend)
      double precision vang
      dimension vang(max_nch_bend)
      double precision vvib
      dimension vvib(max_nch_bend)
      double precision xcord,ycord,zcord
      dimension xcord(4),ycord(4),zcord(4)
      double precision ddum3
      dimension ddum3(3)
      double precision ddum4
      dimension ddum4(4)
      double precision vtg
      dimension vtg(max_nch_bend),vphi(max_nch_bend)
      double precision cosang,crossang
      dimension crossang(0:NNBOND,0:NNBOND)
      dimension cosang(2)
      double precision pdf_true,pdf_bias,angle,pdf_ratio
      double precision mean,lolimit,hilimit
      double precision delta_angle
      dimension delta_angle(CBMC_MAXBENDPEAKS)
      double precision delta_anglec
      dimension delta_anglec(CBMC_MAXBENDPEAKS)
      double precision angsq
      dimension angsq(CBMC_MAXBENDPEAKS)

c      write(6,*) 'START GETCBANGLE'

c     --- retrieve constants
      call twh_constant_pi(GLB_GET,onepi)
      call twh_constant_twopi(GLB_GET,twopi)
      call twh_beta(GLB_GET,beta)
      call twh_cbmc_bendgen_code(GLB_GET,cbmc_bendgen_code)
      call twh_bend_a_ideal_fraction(GLB_GET,bend_a_ideal_fraction)
      call twh_bend_b_ideal_fraction(GLB_GET,bend_b_ideal_fraction)
c     --- assign grownum and growfrom to local variables
      call twh_grownum(GLB_GET,iw,ntogrow)
      call twh_growfrom(GLB_GET,iw,iufrom)

c     --- initialize trial energies
      gcba_vvib = 0.0d0
      gcba_vbend = 0.0d0
      gcba_vtor = 0.0d0
c     --- initialize the local lterm which is used in getweight
      lterm = .false.
c     --- initialize overlap for the full angle selection process
      lovrangle = .false.
      hibound = 0.0d0
      lobound = 0.0d0
      rdecide = 0.0d0

      hybrid = twh_get_bondpatt(twh_get_ntype(imolty,iufrom))

c     --- initialize gcba_weight and dpcbangle
      gcba_weight = 1.0d0
      dpcbangle = 0.0d0

c     ***********************************************************
c     ** Begin Bond angle biased selection                     **
c     ***********************************************************
c     --- store r*ui positions for unit iufrom
      call twh_coordinates(GLB_GET,testcord,ichain,iufrom
     &     ,rxfrom,ryfrom,rzfrom)

c     --- fetch the vecprev
      call twh_vecprev(GLB_GET,xvecprev,yvecprev,zvecprev)
c     --- determine iuprev
      call twh_growprev(GLB_GET,iw,growprev)
      if ( growprev .eq. 0 ) then
c        --- need to choose one bead to grow on unit sphere
         call twh_growlist(GLB_GET,iw,1,iuprev)
         angstart = 2

         if ( luseoldone ) then
c           --- compute the vector from iufrom to iuprev
            call twh_coordinates(GLB_GET,testcord,ichain,iuprev
     &           ,rxprev,ryprev,rzprev)
            call twh_mimage( ibox,rxprev,ryprev,rzprev
     &           ,rxfrom,ryfrom,rzfrom,xvecprev,yvecprev,zvecprev)
            distance = twh_distance(xvecprev,yvecprev,zvecprev)
            call twh_distprev(GLB_SET,1,distance)
         else
c           --- set distprev to the first bondlength
            call twh_bondlen(GLB_GET,1,distance)
            call twh_distprev(GLB_SET,1,distance)
         endif
      else
c        --- already did all this stuff back up in getcbbond
         angstart = 1
      endif

c     --- initialize angle-angle list
      do count = 1,ntogrow
         do step = 1,3 
            call twh_cbmcaatot(GLB_SET,count,step,0)
         enddo
      enddo

c     --- check for angle-angle and set up list
      do ic = 1,twh_get_inaa(imolty,iufrom)
         do ibend = 0,2
            call twh_cbmcaalocl(GLB_SET,ibend,ic,0)
         enddo
c        --- unless ijaa0 is iuprev or we have an xtra 
c        --- we will do this on step 2
         step = 2
c        --- determine which (if any) of the slots is iuprev
         if ( twh_get_ijaa0(imolty,iufrom,ic) .eq. iuprev ) then
            call twh_cbmcaalocl(GLB_SET,0,ic,-1)
            step = 1
         elseif ( twh_get_ijaa1(imolty,iufrom,ic) .eq. iuprev ) then
            call twh_cbmcaalocl(GLB_SET,1,ic,-1)
         elseif ( twh_get_ijaa2(imolty,iufrom,ic) .eq. iuprev ) then
            call twh_cbmcaalocl(GLB_SET,2,ic,-1)
         endif

c        --- determine which (if any) of the slots is filled by xtra
         call twh_growxtra(GLB_GET,iw,1,xtracount)
         do ixtra = 2,xtracount+1
            call twh_growxtra(GLB_GET,iw,ixtra,iuxtra)
            if ( twh_get_ijaa0(imolty,iufrom,ic) .eq. iuxtra ) then
               call twh_cbmcaalocl(GLB_SET,0,ic,-ixtra)
               step = 3
            elseif ( twh_get_ijaa1(imolty,iufrom,ic) .eq. iuxtra ) then
               call twh_cbmcaalocl(GLB_SET,1,ic,-ixtra)
               step = 3
            elseif ( twh_get_ijaa2(imolty,iufrom,ic) .eq. iuxtra ) then
               call twh_cbmcaalocl(GLB_SET,2,ic,-ixtra)
               step = 3
            endif
         enddo

         do count = angstart,ntogrow
            call twh_growlist(GLB_GET,iw,count,iugrow)
            if ( twh_get_ijaa0(imolty,iufrom,ic) .eq. iugrow ) then
               call twh_cbmcaalocl(GLB_SET,0,ic,count)
            endif
            if ( twh_get_ijaa1(imolty,iufrom,ic) .eq. iugrow ) then
               call twh_cbmcaalocl(GLB_SET,1,ic,count)
            endif
            if ( twh_get_ijaa2(imolty,iufrom,ic) .eq. iugrow ) then
               call twh_cbmcaalocl(GLB_SET,2,ic,count)
            endif

            call twh_cbmcaalocl(GLB_GET,0,ic,testzero)
            call twh_cbmcaalocl(GLB_GET,1,ic,testone)
            call twh_cbmcaalocl(GLB_GET,2,ic,testtwo)
            if ( (testzero*testone*testtwo) .ne. 0 ) then
c              --- we now have enough atoms to check the angle-angle
               call twh_cbmcaatot(GLB_INCR,count,step,1)
               call twh_cbmcaatot(GLB_GET,count,step,ivalue)
               call twh_cbmcaatype(GLB_SET,count,step,ivalue,ic)
c              --- get out of loop so we don't count this again
               goto 12
            endif
         enddo
         write(6,*) 'GETCBANGLE: angle-angle logic messed up'
         lfailure = .true.
         return
 12      continue
      enddo

c     --- determine the iugrow-iufrom-iuprev angles 
      do count = angstart,ntogrow
c        --- initialize lmin which tells us if we have a minimum yet
         lmin = .false.
         call twh_growlist(GLB_GET,iw,count,iugrow)
c        --- default number of trial selections is nch_bend_a (this can
c        --- be changed later depending on the generation distribution)
         nchoice = nch_bend_a
c        --- set up distances for angle calculations
         call twh_bondlen(GLB_GET,count,distgrow)
         call twh_distprev(GLB_GET,1,distance)
c        ---- check to see if this is also a bond (3-mem ring)
         lthreemem = .false.
         call twh_invib(GLB_GET,imolty,iugrow,numvib)
         do ivib = 1,numvib
            call twh_ijvib(GLB_GET,imolty,iugrow,ivib,iutest)
            if ( iutest .eq. iuprev ) then
c              --- also forms a 3-mem ring
               lthreemem = .true.
c              --- get the type of this extra iuprev-iugrow bond
               call twh_itvib(GLB_GET,imolty,iugrow,ivib,vibtype)
c              --- get the vibration style
               call twh_vibstyle(GLB_GET,vibtype,vibstyle)
            endif
         enddo
         if ( lthreemem ) then
c           --- store the angle information for the 2 other angles
c           --- completed by the 3-membered ring
c           --- find the from-grow-prev
            lfound = .false.
            ifind = 0
            call twh_inben(GLB_GET,imolty,iufrom,nbend)
            do while ( .not. lfound .and. (ifind .lt. nbend) )
               ifind = ifind + 1
               call twh_ijben2(GLB_GET,imolty,iufrom,ifind,iulast)
               if ( iulast .eq. iugrow ) then
                  call twh_ijben3(GLB_GET,imolty,iufrom,ifind,iu2back)
                  if ( iu2back .eq. iuprev ) then
                     lfound = .true.
                     fgp_index = ifind
                  endif
               endif
            enddo
            if ( .not. lfound ) then
               write(6,*) 'GETCBANGLE: problem finding f-g-p'
               lfailure = .true.
               return
            endif
c           --- find the grow-prev-from
            lfound = .false.
            ifind = 0
            call twh_inben(GLB_GET,imolty,iugrow,nbend)
            do while ( .not. lfound .and. (ifind .lt. nbend) )
               ifind = ifind + 1
               call twh_ijben2(GLB_GET,imolty,iugrow,ifind,iulast)
               if ( iulast .eq. iuprev ) then
                  call twh_ijben3(GLB_GET,imolty,iugrow,ifind,iu2back)
                  if ( iu2back .eq. iufrom ) then
                     lfound = .true.
                     gpf_index = ifind
                  endif
               endif
            enddo
            if ( .not. lfound ) then
               write(6,*) 'GETCBANGLE: problem finding g-p-f'
               lfailure = .true.
               return
            endif
         endif
c        --- find the angle type for iugrow-iufrom-iuprev
         lfound = .false.
         ibend = 0
         style_code = DIST_NULL
         do while ( .not. lfound 
     &        .and. ibend .lt. twh_get_inben(imolty,iugrow) )
            ibend = ibend + 1
            iulast = twh_get_ijben2(imolty,iugrow,ibend)
            if ( iulast .eq. iufrom ) then
               iu2back = twh_get_ijben3(imolty,iugrow,ibend)
               if ( iu2back .eq. iuprev ) then
                  lfound = .true.
                  index = ibend
                  type = twh_get_itben(imolty,iugrow,ibend)
               endif
            endif
         enddo
         if ( .not. lfound ) then
            write(6,*) 'GETCBANGLE: angle not found in getcbangle'
            write(6,*) 'count,iugrow,iufrom,iuprev',count,iugrow
     &           ,iufrom,iuprev
            lfailure = .true.
            return
         endif

         call twh_anglestyle(GLB_GET,type,angle_style)
         call twh_cbmcbenddist(GLB_GET,imolty,iugrow,index,index
     &        ,angleadist)
         if ( angleadist .eq. DIST_TRIANGLE ) then
c           --- rigid connects the two end of this angle to make a 
c           --- triangle.  Use a delta distribution based on the rigid
c           --- bond length to compute the allowed angle
            style_code = DIST_DELTA
            call twh_cbmcbendcount(GLB_GET,imolty,iugrow,index,index
     &           ,npeaks)
            do ipeak = 1,npeaks
c              --- combine the allowed bond lengths with existing
c              --- lengths to compute the allowed angle
               call twh_vibcoeff(GLB_GET,vibtype,0,testlen)
               call twh_cbmcbendinfo(GLB_GET,imolty,iugrow,index,index
     &              ,ipeak,1,testlen)
               anglec = twh_inverselawofcosine(testlen,distance
     &              ,distgrow)
               angle = twh_arccos(anglec)
               delta_angle(ipeak) = angle
               delta_anglec(ipeak) = anglec
            enddo
         elseif ( angleadist .eq. DIST_DELTA ) then
c           --- rigid angle, so we flip the usual process of
c           --- sampling based on the angle and weighting on
c           --- the Boltzman weight to sampling based upon
c           --- a discrete distribution consisting of every 
c           --- angle that has a non-zero Boltzman weight
c           --- and then weighting based upon the angle
c           --- distribution plus the Boltzman weight
            style_code = DIST_DELTA
            call twh_cbmcbendcount(GLB_GET,imolty,iugrow,index,index
     &           ,npeaks)
            do ipeak = 1,npeaks
               call twh_cbmcbendinfo(GLB_GET,imolty,iugrow,index
     &              ,index,ipeak,1,angle)
               anglec = dcos(angle)
               delta_angle(ipeak) = angle
               delta_anglec(ipeak) = anglec
            enddo
c           --- no reason to sample more than the finite number of
c           --- possible angles
            if ( nchoice .ge. npeaks ) then
               nchoice = npeaks
               lrandelta = .false.
            endif
         elseif ( angleadist .eq. DIST_BOUNDED_SINE ) then
c           --- narrow range of allowed angles depending on
c           --- the bond lengths and allowed distances
            style_code = DIST_BOUNDED_SINE
            call twh_cbmcbendcount(GLB_GET,imolty,iugrow,index,index
     &           ,npeaks)
            if ( npeaks .ne. 1 ) then
               write(6,*) 'GETCBANGLE: too many peaks for bounded_sine'
               write(6,*) 'imolty,iugrow,index,npeaks:',imolty,iugrow
     &              ,index,npeaks
               lfailure = .true.
               return
            endif
            call twh_cbmcbendinfo(GLB_GET,imolty,iugrow,index,index
     &           ,npeaks,1,distlo)
            call twh_cbmcbendinfo(GLB_GET,imolty,iugrow,index,index
     &           ,npeaks,2,disthi)
            lobound = twh_inverselawofcosine(distance,
     &           distgrow,distlo)
            hibound = twh_inverselawofcosine(distance,
     &           distgrow,disthi)
         elseif ( angleadist .eq. DIST_SINE ) then
            style_code = DIST_SINE
            npeaks = 0
         elseif ( angleadist .eq. DIST_GAUSSIAN ) then
            if ( cbmc_bendgen_code
     &           .eq. CBMC_IDEAL_PLUS_AUTOFIT_GAUSSIAN ) then
               style_code = DIST_SINE_GAUSSIAN
            else
               style_code = DIST_GAUSSIAN
            endif
            call twh_cbmcbendcount(GLB_GET,imolty,iugrow,index
     &           ,index,npeaks)
            if ( npeaks .ne. 1 ) then
               write(6,*) 'GETCBANGLE: angle A unsure how to handle'
     &              ,' multiple npeaks for Gaussian'
               write(6,*) 'npeaks:',npeaks
               lfailure = .true.
               return
            endif
            call twh_cbmcbendinfo(GLB_GET,imolty,iugrow,index
     &           ,index,npeaks,1,mean)
            call twh_cbmcbendinfo(GLB_GET,imolty,iugrow,index
     &           ,index,npeaks,2,std_dev)
            call twh_cbmcbendinfo(GLB_GET,imolty,iugrow,index
     &           ,index,npeaks,3,lolimit)
            call twh_cbmcbendinfo(GLB_GET,imolty,iugrow,index
     &           ,index,npeaks,4,hilimit)
         endif

c        --- compute trial angles and energies
         do itrial = 1,nchoice
c           --- initialize lovrlap
            lovrlap = .false.
            if ( luseoldone .and. itrial .eq. 1) then
c              --- first itrial is the old conformation
c              --- compute vector from iufrom to iugrow
               call twh_coordinates(GLB_GET,testcord,ichain,iugrow
     &              ,rxgrow,rygrow,rzgrow)
               call twh_mimage(ibox,rxgrow,rygrow,rzgrow
     &              ,rxfrom,ryfrom,rzfrom,xvecgrow,yvecgrow,zvecgrow)
c              --- dot product divided by lengths gives cos(angle)
               anglec = twh_dotproduct(xvecprev,yvecprev,zvecprev
     &              ,xvecgrow,yvecgrow,zvecgrow,dum,dumb)
               angle = twh_arccos(anglec)
            else
               if ( style_code .eq. DIST_DELTA ) then
c                 --- delta functions, choose from a finite number of 
c                 --- peaks.
                  if ( lrandelta ) then
c                    --- fewer trials than finite allowed angles, select
c                    --- one at random
                     ipeak = idint(twh_random()*dble(npeaks)) + 1
                  else
c                    --- No need for randomness as we can sample
c                    --- the full distribution
                     ipeak = itrial
                  endif
                  angle = delta_angle(ipeak)
                  anglec = delta_anglec(ipeak)
               else
                  if ( style_code .eq. DIST_SINE_GAUSSIAN ) then
                     rdecide = twh_random()
                  endif

                  if ( style_code .eq. DIST_SINE .or.
     &                 (style_code .eq. DIST_SINE_GAUSSIAN 
     &                 .and. rdecide .le. bend_a_ideal_fraction )
     &                 ) then
c                    --- choose an angle with probability sin(theta)
                     anglec = 2.0*twh_random() - 1.0d0
                     angle = twh_arccos(anglec)
                  elseif ( style_code .eq. DIST_GAUSSIAN .or.
     &                    ( style_code .eq. DIST_SINE_GAUSSIAN 
     &                    .and. rdecide .gt. bend_a_ideal_fraction )
     &                    ) then
c                    --- choose an angle based on gaussian distribution
                     angle = twh_gaussian(lfailure,mean,std_dev
     &                    ,lolimit,hilimit)
                     if ( lfailure ) return
                     anglec = dcos(angle)
                  elseif ( style_code .eq. DIST_BOUNDED_SINE ) then
                     anglec = (hibound-lobound)*twh_random() + lobound
                     angle = twh_arccos(anglec)
                  else
                     write(6,*) 'GETCBANGLE: bend A unknown style_code'
                     write(6,*) style_code
                     lfailure = .true.
                     return
                  endif
               endif
            endif

            if ( style_code .eq. DIST_DELTA ) then
c              --- unlike the others, we are sampling from the complete
c              --- list of angles that have a non-zero Boltzman weight
c              --- (usually just one rigid angle).  Our sampling
c              --- distribution is discrete with an equal probability
c              --- of chosing any of the non-infinite energy angles.
c              --- if there ever was more than one peak this would have
c              --- to be modified to have the ratio of the angle
c              --- probability factors
               pdf_ratio = 1.0d0
            elseif ( style_code .eq. DIST_SINE ) then
               pdf_bias = 0.5d0*dsin(angle)
            elseif ( style_code .eq. DIST_GAUSSIAN ) then
               pdf_bias = twh_gaussprob(angle,mean,std_dev,lolimit
     &              ,hilimit)/dble(npeaks)
            elseif ( style_code .eq. DIST_SINE_GAUSSIAN ) then
               pdf_bias = bend_a_ideal_fraction*0.5d0*dsin(angle)
     &              + (1.0d0-bend_a_ideal_fraction)
     &              *twh_gaussprob(angle,mean,std_dev,lolimit,hilimit)
     &              /dble(npeaks)
            elseif ( style_code .eq. DIST_BOUNDED_SINE ) then
               pdf_bias = dsin(angle)/(lobound-hibound)
            else
               write(6,*) 'GETCBANGLE: unknown cbmc_bendgen_code'
               write(6,*) cbmc_bendgen_code
               lfailure = .true.
               return
            endif
c           --- safety check pdf_bias
            if ( style_code .ne. DIST_DELTA 
     &           .and. pdf_bias .le. 0.0d0 ) then
               write(6,*) 'GETCBANGLE: invalid bend A pdf_bias'
               write(6,*) 'lfirstold:',lfirstold
               write(6,*) 'pdf_bias:',pdf_bias
               if ( cbmc_bendgen_code .eq. CBMC_IDEAL ) then
                  write(6,*) 'angle,dsin(angle)',angle,dsin(angle)
               else
                  write(6,*) 'angle,mean,std_dev',angle,mean,std_dev
               endif
               lfailure = .true.
               return
            endif
c           --- initialize vang
            vang(itrial) = 0.0d0
c           --- assign ang_trial
            ang_trial(itrial) = anglec
            if ( lthreemem ) then
c              --- calcule the prev-grow bond energy using the law
c              --- of cosine
               testlen = twh_lawofcosine(distgrow,distance,anglec)
               vvib(itrial) = twh_vbond(lfailure,vibtype,testlen
     &              ,imolty,iugrow,iuprev,lovrlap)
               if ( lfailure ) return
c              --- connecting this bond also computes 2 additional
c              --- angles for the rest of the triangle
c              --- get the angle for the from-grow-prev
               fgp_anglec = twh_inverselawofcosine(distgrow,testlen
     &              ,distance)
               vang(itrial) = vang(itrial)
     &              + twh_vangle(lfailure,lovrlap,ibox,imolty,iufrom
     &              ,fgp_index,0,fgp_anglec,distgrow,testlen
     &              ,ddum3,ddum3,ddum3)
               if ( lfailure ) return
c              --- get the angle for the grow-prev-from
               gpf_anglec = twh_inverselawofcosine(testlen,distance
     &              ,distgrow)
               vang(itrial) = vang(itrial)
     &              + twh_vangle(lfailure,lovrlap,ibox,imolty,iugrow
     &              ,gpf_index,0,gpf_anglec,testlen,distance
     &              ,ddum3,ddum3,ddum3)
               if ( lfailure ) return
            else
c              --- no bond energy
               vvib(itrial) = 0.0d0
            endif
c           --- calculate the bond angle energy
            vang(itrial) = vang(itrial)
     &           + twh_vangle(lfailure,lovrlap,ibox,imolty,iugrow
     &           ,index,0,anglec,distgrow,distance,ddum3,ddum3,ddum3)
            if ( lfailure ) return
c           --- no need to compute any more interactions if we 
c           -- have an infinite energy overlap: otherwise compute angang
            if ( .not. lovrlap ) then
c              --- calculate angle-angle cross terms
               call twh_cbmcaatot(GLB_GET,count,1,total)
               do itor = 1, total
                  call twh_cbmcaatype(GLB_GET,count,1,itor,ic)
c                 --- on this step iuzero must be iuprev so we use
c                 --- bendang
                  call twh_cbmcaalocl(GLB_GET,1,ic,testone)
                  if ( testone .eq. count ) then
c                    --- normal order
                     cosang(1) = anglec
                     call twh_cbmcaalocl(GLB_GET,2,ic,other)
                     call twh_bendang(GLB_GET,other,bendang)
                     cosang(2) = dcos(bendang)
                  else
c                    --- reverse order
                     cosang(2) = anglec
                     call twh_cbmcaalocl(GLB_GET,1,ic,other)
                     call twh_bendang(GLB_GET,other,bendang)
                     cosang(1) = dcos(bendang)
                  endif
                  vang(itrial) = vang(itrial)
     &                 + twh_vangang(AA_ANGLE,ibox,imolty
     &                 ,iufrom,ic,cosang,ddum4,ddum4,ddum4)
               enddo
            endif

            if ( lovrlap ) then
               lovrben(itrial) = .true.
            else
               lovrben(itrial) = .false.
               vtotal = vvib(itrial) + vang(itrial)
               if ( lmin ) then
c                 --- see if this is the new minimum energy
                  if ( vtotal .lt. vmin ) vmin = vtotal
               else
c                 --- this is our first non-overlap energy, set vmin
                  vmin = vtotal
                  lmin = .true.
               endif
            endif
c           --- compute the ideal probability density
            if ( style_code .eq. DIST_DELTA ) then
c              --- the ratio is used directly in a discrte dist.
               bfactor(itrial) = pdf_ratio
            else
c              --- sampling from a continuous distribution, use ideal
               pdf_true = 0.5d0*dsin(angle)
               if ( pdf_true .lt. SMALLEST ) then
c                 --- bad news, reset to a very small value instead
                  pdf_true = SMALLEST
               endif
c              --- store the ratio of the true probability density to 
c              --- the arbitrary one that was used to generate the trial
               bfactor(itrial) = pdf_true/pdf_bias
            endif
         enddo

c        --- initialize bsum_try and lovrlap
         lovrangle = .true.
         bsum_try = 0.0d0
         do itrial = 1,nchoice
            if ( lovrben(itrial) ) then
               bfactor(itrial) = 0.0d0
            else
               vtotal = vvib(itrial) + vang(itrial)
               bfactor(itrial) 
     &              = bfactor(itrial)*twh_expon(-beta*(vtotal-vmin))
               lovrangle = .false.
            endif
            bsum_try = bsum_try + bfactor(itrial)
         enddo

         if ( lovrangle ) return

c        --- initialize the winner variable
         iwinner = 0
         if ( luseoldone ) then
c           --- select the old conformation
            iwinner = 1
         elseif ( nchoice .eq. 1 ) then
c           --- only one choice, no need for a random number
c           --- accept unless it is an overlap
            if ( .not. lovrben(1) ) iwinner = 1
         else
c           --- select one of the trial sites via bias
            rbf = twh_random()*bsum_try
            bs = 0.0d0
            do while ( bs .lt. rbf .and. iwinner .lt. nchoice )
               iwinner = iwinner + 1
               bs = bs + bfactor(iwinner)
            enddo
         endif

c        --- account for vmin
         call twh_getweight(lfailure,0,vmin,wmin,dpmin,lterm)
         if ( lfailure ) return
         if ( lterm ) then
c           --- this should not happen
            write(6,*) 'GETCBANGLE: getweight returned lterm in bend'
            lfailure = .true.
            return
         endif

c        --- propagate the rosenbluth weight and dpcbangle
         gcba_weight = gcba_weight*wmin*bsum_try/dble(nchoice)
         dpcbangle = dpcbangle + dpmin

c        --- store accepted angle for use later
         bendang = twh_arccos(ang_trial(iwinner))
         call twh_bendang(GLB_SET,count,bendang)

c        --- store angle as a cross with iuprev for angle-angle term
         crossang(count,0) = ang_trial(iwinner)
         crossang(0,count) = ang_trial(iwinner)

c        --- increment energy
         gcba_vvib = gcba_vvib + vvib(iwinner)
         gcba_vbend = gcba_vbend + vang(iwinner)
#if DEBUG_GETCBANGLE
      write(6,*) 'winning vvib, vang:',vvib(iwinner),vang(iwinner)
#endif

      enddo

c     --- start of the "B" selection
c     --- initialize torsion list
      do count = 1,ntogrow
         imptot(count) = 0
      enddo
      rimptotal = 0
      call twh_cbmcrimptype(GLB_SET,0,rimptotal)

c     --- check for improper torsions and set up list
      do ic = 1,twh_get_inimprop(imolty,iufrom)
c        --- initialize lingeom
         lingeom = .true.
c        --- initialize cbmcimproper
         do ibend = 2,4
            call twh_cbmcimproper(GLB_SET,ibend,ic,0)
         enddo
c        --- see if any of the atoms needed for this torsion are iuprev
         if ( twh_get_ijimprop2(imolty,iufrom,ic) .eq. iuprev ) then
            ivalue = -iuprev
            call twh_cbmcimproper(GLB_SET,2,ic,ivalue)
         endif
         if ( twh_get_ijimprop3(imolty,iufrom,ic) .eq. iuprev ) then
            ivalue = -iuprev
            call twh_cbmcimproper(GLB_SET,3,ic,ivalue)
         endif
         if ( twh_get_ijimprop4(imolty,iufrom,ic) .eq. iuprev ) then
            ivalue = -iuprev 
            call twh_cbmcimproper(GLB_SET,4,ic,ivalue)
         endif

c        --- see if any of the atoms needed are in growxtra
         call twh_growxtra(GLB_GET,iw,1,xtracount)
         do ixtra = 2,xtracount+1
            call twh_growxtra(GLB_GET,iw,ixtra,iuxtra)
            if ( twh_get_ijimprop2(imolty,iufrom,ic) .eq. iuxtra ) then
               ivalue = -iuxtra
               call twh_cbmcimproper(GLB_SET,2,ic,ivalue)
               lingeom = .false.
            endif
            if ( twh_get_ijimprop3(imolty,iufrom,ic) .eq. iuxtra ) then
               ivalue = -iuxtra
               call twh_cbmcimproper(GLB_SET,3,ic,ivalue)
               lingeom = .false.
            endif
            if ( twh_get_ijimprop4(imolty,iufrom,ic) .eq. iuxtra ) then
               ivalue = -iuxtra
               call twh_cbmcimproper(GLB_SET,4,ic,ivalue)
               lingeom = .false.
            endif
         enddo

c        --- see if any of the atoms needed are being grown
         do count = angstart,ntogrow
            call twh_growlist(GLB_GET,iw,count,iugrow)
            if ( twh_get_ijimprop2(imolty,iufrom,ic) .eq. iugrow ) then
               call twh_cbmcimproper(GLB_SET,2,ic,count)
            endif
            if ( twh_get_ijimprop3(imolty,iufrom,ic) .eq. iugrow ) then
               call twh_cbmcimproper(GLB_SET,3,ic,count)
            endif
            if ( twh_get_ijimprop4(imolty,iufrom,ic) .eq. iugrow ) then
               call twh_cbmcimproper(GLB_SET,4,ic,count)
            endif
c           --- get the product of the cbmcimproper values
            call twh_cbmcimproper(GLB_GET,2,ic,testtwo)
            call twh_cbmcimproper(GLB_GET,3,ic,testthree)
            call twh_cbmcimproper(GLB_GET,4,ic,testfour)
            if ( ( testtwo .ne. 0 ) .and. ( testthree .ne. 0 )
     &           .and. ( testfour .ne. 0) ) then
c              --- we now have enough atoms to check the improper charmm
               if ( lingeom ) then
c                 --- all of the atoms exist in this subroutine
                  imptot(count) = imptot(count) + 1
                  tmpimptype(count,imptot(count)) = ic
                  tmpimpform(count,imptot(count)) 
     &                 = twh_get_impform(
     &                 twh_get_itimprop(imolty,iufrom,ic))
               else
c                 --- we need to pass this back to rosenbluth
                  rimptotal = rimptotal + 1
                  call twh_cbmcrimptype(GLB_SET,0,rimptotal)
                  call twh_cbmcrimptype(GLB_SET,rimptotal,ic)
               endif
c              --- get out of loop so we don't count this again
               goto 22
            endif
         enddo
         write(6,*) 'improper torsion logic messed up in getcbangle'
         lfailure = .true.
         return
 22      continue
      enddo

c     --- skip angstart in the loop below
      start_ang = angstart+1
c     --- initialize cbphi(angstart) to 0.0
      call twh_cbphi(GLB_SET,angstart,0.0d0)
      if ( luseoldone ) then
c        --- need to determine angstart
         if ( angstart .le. ntogrow) then
c           --- set up the cone using iuprev
            call twh_distprev(GLB_GET,1,distance)
            xub = -xvecprev/distance
            yub = -yvecprev/distance
            zub = -zvecprev/distance
            call twh_cone(lfailure,1,xub,yub,zub,dum,dum,dum,dum,dum)
            if ( lfailure ) return

            call twh_growlist(GLB_GET,iw,angstart,iugrow)
c           --- compute vector from iufrom to iugrow
            call twh_coordinates(GLB_GET,testcord,ichain,iugrow
     &           ,rxgrow,rygrow,rzgrow)
            call twh_mimage(ibox,rxgrow,rygrow,rzgrow
     &           ,rxfrom,ryfrom,rzfrom,xvecgrow,yvecgrow,zvecgrow)
            call twh_bondlen(GLB_GET,angstart,distgrow)
c           --- turn this into a unit vector
            ux = xvecgrow/distgrow
            uy = yvecgrow/distgrow
            uz = zvecgrow/distgrow
            call twh_bendang(GLB_GET,angstart,alpha)
c           --- compute gamma (phi) for these unit vectors
            call twh_cone(lfailure,3,dum,dum,dum,alpha,gamma,ux,uy,uz)
            if ( lfailure ) return
            call twh_cbphi(GLB_SET,angstart,gamma)
         endif
      endif

c     --- determine the angles of the grown beads relative to anglestart
      do count = start_ang,ntogrow
c        --- initialize lmin which tells us we havn't found the min yet
         lmin = .false.
         call twh_growlist(GLB_GET,iw,count,iugrow)
c        --- check to see if there are any bonds (3-mem ring), if they
c        --- they are rigid, and if there are any rigid angle types
c        --- between the- current value of count and any
c        --- atom already grown during this step
         lrigid = .false.
         lthreemem = .false.
         ntmtorsion = 0
         do aaa = angstart,count-1
            call twh_growlist(GLB_GET,iw,aaa,iutest)
c           --- see if there is a bond between iutest and iugrow
            call twh_invib(GLB_GET,imolty,iutest,numvib)
            do ivib = 1,numvib
               call twh_ijvib(GLB_GET,imolty,iutest,ivib,ifind)
               if ( ifind .eq. iugrow ) then
c                 --- bond between growing beads found
                  lthreemem = .true.
c                 --- check the bond type
                  call twh_itvib(GLB_GET,imolty,iutest,ivib,vibtype)
                  call twh_vibstyle(GLB_GET,vibtype,vibstyle)
                  if ( vibstyle .eq. BOND_STYLE_FIXED ) then
                     lrigid = .true.
c                    --- get the rigid bond length
                     call twh_vibcoeff(GLB_GET,vibtype,0,testlen)
                     call twh_bondlen(GLB_GET,angstart,distother)
c                    --- convert to angle
                     anglec = twh_inverselawofcosine(distgrow
     &                    ,distother,testlen)
                     ang_gamma = twh_arccos(anglec)
                  endif
               endif
            enddo
            if ( .not. lrigid ) then
               call twh_growlist(GLB_GET,iw,aaa,iutest)
c              --- check iutest - iufrom - iugrow to see if it is rigid
               call twh_inben(GLB_GET,imolty,iugrow,nangle)
               lfound = .false.
               iangle = 0
               do while ( .not. lfound .and. (iangle .lt. nangle))
                  iangle = iangle + 1
                  call twh_ijben2(GLB_GET,imolty,iugrow,iangle,testtwo)
                  if ( testtwo .eq. iufrom ) then
c                    --- check the third atom
                     call twh_ijben3(GLB_GET,imolty,iugrow,iangle
     &                    ,testthree)
                     if ( testthree .eq. iutest ) then
c                       --- match found
                        lfound = .true.
c                       --- check the style of the angle
                        call twh_itben(GLB_GET,imolty,iugrow,iangle
     &                       ,type)
                        call twh_anglestyle(GLB_GET,type,angle_style)
                        if ( angle_style .eq. ANG_STYLE_RIGID ) then
                           lrigid = .true.
c                          --- all rigid angles, will not change 
                           delta_key = aaa
c                          --- set the delta angle to the equilibrium
                           call twh_bencoeff(GLB_GET,type,0,ang_gamma)
                        endif
                     endif
                  endif
               enddo
            endif
         enddo
c        --- set the default number of choices to nch_bend_b
         nchoice = nch_bend_b
         if ( lrigid ) then
c           --- setup the delta distribution
c           --- fetch ang_alpha
            call twh_bendang(GLB_GET,delta_key,ang_alpha)
c           --- fetch ang_beta
            call twh_bendang(GLB_GET,count,ang_beta)
c           --- already grabbed ang_gamma so compute the delta dist
            call twh_angle_b_rigid(ang_alpha,ang_beta,ang_gamma
     &           ,npeaks,delta_angle)
c           --- set the number of choices to npeaks
            if ( npeaks .le. nchoice ) then
c              --- sample each possible allowed angle once
               nchoice = npeaks
               lrandelta = .false.
            else
c              --- randomly choose from the list of possible angles
               lrandelta = .true.
            endif
            if ( nchoice .eq. 0 ) then
               write(6,*) 'GETCBANGLE: failure to generate valid angle'
               write(6,*) 'for step B using the DIST_DELTA'
               lfailure = .true.
               return
            endif
         endif
c        --- compute trial energies and weights
         do itrial = 1,nchoice
            lovrlap = .false.
            if ( luseoldone .and. itrial .eq. 1) then
c              --- compute vector from iufrom to iugrow
               call twh_coordinates(GLB_GET,testcord,ichain,iugrow
     &              ,rxgrow,rygrow,rzgrow)
               call twh_mimage(ibox,rxgrow,rygrow,rzgrow
     &              ,rxfrom,ryfrom,rzfrom,xvecgrow,yvecgrow,zvecgrow)
               call twh_bondlen(GLB_GET,count,distgrow)
c              --- turn this into a unit vector
               ux = xvecgrow/distgrow
               uy = yvecgrow/distgrow
               uz = zvecgrow/distgrow
               call twh_bendang(GLB_GET,count,alpha)
c              --- compute gamma (phi) for these unit vectors
               call twh_cone(lfailure,3,dum,dum,dum,alpha,gamma
     &              ,ux,uy,uz)
               if ( lfailure ) return
               phitwo = gamma
c              --- set the flag so we know this is the first old one
               lfirstold = .true.
            else
               lfirstold = .false.
            endif

            if ( lfirstold ) then
c              --- need to "unshift" the first angle to get it in the
c              --- proper relative terms for the biasing
               call twh_cbphi(GLB_GET,angstart,cbphi)
               phitwo = phitwo - cbphi
c              --- make sure this angle is in (-Pi,Pi)
               if ( phitwo .lt. -onepi ) then
                  phitwo = phitwo + twopi
               elseif ( phitwo .gt. onepi ) then
                  phitwo = phitwo - twopi
               endif
            endif

            style_code = DIST_NULL
c           --- determine a value of phitwo
            if ( lrigid ) then
c              --- delta function
               style_code = DIST_DELTA
c              --- get the shift from cbphi based upon the delta_key
               call twh_cbphi(GLB_GET,delta_key,delta_shift)
               call twh_cbphi(GLB_GET,angstart,cbphi) 
               delta_shift = delta_shift - cbphi
               if ( lfirstold ) then
c                 --- if lfirstold then reset the order of the 
c                 --- delta_angle to make sure the real old one is the
c                 --- first entry
                  imin = 0
                  do ipeak = 1,npeaks
                     angsq(ipeak) = ((delta_angle(itrial) + delta_shift)
     &                    - phitwo)**2
                     if ( ipeak .eq. 1 ) then
                        imin = ipeak
                     elseif ( angsq(ipeak) .lt. angsq(imin) ) then
                        imin = ipeak
                     endif
                  enddo
c                 --- switch 1 with imin
                  angletemp = delta_angle(1)
                  delta_angle(1) = delta_angle(imin)
                  delta_angle(imin) = angletemp
               endif
            elseif ( cbmc_bendgen_code .eq. CBMC_IDEAL  ) then
c              --- uniform distribution
               style_code = DIST_UNIFORM
            elseif ( cbmc_bendgen_code .eq. CBMC_GLOBAL_GAUSSIAN ) then
c              --- global standard deviation
               call twh_sdevbenb(GLB_GET,std_dev)
c              --- we are computing this value relative to angstart
c              --- value.  Depends on the number of atoms that are 
c              --- bonded to the central atom.
               if ( hybrid .eq. 'sp3' ) then
c                 --- sp3 like bonding structure.  This is 
c                 --- threefold symmetric and there is not much 
c                 --- reason to sample very hard in the vicinity of 
c                 --- the first angle.  Do this using two different 
c                 --- gaussian peaks that are centered at 2Pi/3 and 
c                 --- 4Pi/3.  Choose randomly between the two centers
c                 --- and do not take options that cross into the 
c                 --- other region
                  style_code = DIST_GAUSSIAN
                  npeaks = 2
                  if ( lfirstold ) then
c                    --- rdecide is determined based on the phitwo
c                    --- computed up above
                     if ( phitwo .lt. 0.0d0 ) then
c                       --- mean of -2Pi/3
                        rdecide = 0.25
                     else
c                       --- mean of 2Pi/3
                        rdecide = 0.75
                     endif
                  else
c                    --- choose randomly
                     rdecide = twh_random()
                  endif
                  if ( rdecide.lt. 0.5d0 ) then
c                    --- centered at -2Pi/3, and between 0 and -Pi
                     lolimit = -onepi
                     hilimit = 0.0d0
                     mean = -twopi/3.0
                  else
c                    --- centered at 2Pi/3, and between 0 and Pi
                     lolimit = 0.0d0
                     hilimit = onepi
                     mean = twopi/3.0
                  endif
               elseif ( hybrid .eq. 'sp2' .or. hybrid .eq. 'arom') then
c                 --- sp2 like bonding structure.  This is 
c                 --- twofold symmetric and there is not much 
c                 --- reason to sample very hard in the vicinity of 
c                 --- the first angle.  Do this using a single
c                 --- gaussian peak centered at Pi, between 0 and 2Pi
                  style_code = DIST_GAUSSIAN
                  npeaks = 1
                  lolimit = 0.0d0
                  hilimit = twopi
                  mean = onepi
                  if ( lfirstold .and. phitwo .lt. 0.0 ) then
c                    --- switch the mean to negative pi as cone returns
c                    --- angles between -Pi and Pi
                     mean = mean - twopi
                     lolimit = lolimit - twopi
                     hilimit = hilimit - twopi
                  endif
               else
c                 --- have not figured out what to do with this 
c                 --- particular bonding pattern so just use uniform
                  style_code = DIST_UNIFORM
                  npeaks = 0
               endif
            elseif ( cbmc_bendgen_code .eq. CBMC_AUTOFIT_GAUSSIAN .or.
     &              cbmc_bendgen_code 
     &              .eq. CBMC_IDEAL_PLUS_AUTOFIT_GAUSSIAN ) then
               if ( cbmc_bendgen_code .eq. CBMC_AUTOFIT_GAUSSIAN ) then
                  style_code = DIST_GAUSSIAN
               else
                  style_code = DIST_UNIFORM_GAUSSIAN
               endif
c              --- get number of sections based on the bend numbers for
c              --- iuprev - iufrom - growlist(iw,angstart) and
c              --- iuprev - iufrom - iugrow
               icone = 0
               ictwo = 0
               call twh_growlist(GLB_GET,iw,angstart,iuone)
               do ibend = 1,twh_get_inben(imolty,iuprev)
                  if ( twh_get_ijben2(imolty,iuprev,ibend) .eq. iufrom
     &                 .and.
     &                 twh_get_ijben3(imolty,iuprev,ibend) .eq. iuone
     &                 ) then
                     icone = ibend
                  endif
                  if ( twh_get_ijben2(imolty,iuprev,ibend) .eq. iufrom
     &                 .and.
     &                 twh_get_ijben3(imolty,iuprev,ibend) .eq. iugrow
     &                 ) then
                     ictwo = ibend
                  endif
               enddo
               if ( icone .eq. 0 .or. ictwo .eq. 0 ) then
                  write(6,*) 'GETCBANGLE: angle B types not found'
                  lfailure = .true.
                  return
               endif

               call twh_cbmcbendcount(GLB_GET,imolty,iuprev,icone
     &              ,ictwo,npeaks)
               if ( lfirstold ) then
                  call twh_angle_b_distribution(style_code,lolimit
     &                 ,hilimit,mean,std_dev,phitwo,imolty,iuprev
     &                 ,icone,ictwo,npeaks,lfailure)
                  if ( lfailure ) return
               else
c                 --- select a peak at random
                  ipeak = int(twh_random()*dble(npeaks)) + 1
                  call twh_cbmcbendinfo(GLB_GET,imolty,iuprev,icone
     &                 ,ictwo,ipeak,1,mean)
                  call twh_cbmcbendinfo(GLB_GET,imolty,iuprev,icone
     &                 ,ictwo,ipeak,2,std_dev)
                  call twh_cbmcbendinfo(GLB_GET,imolty,iuprev,icone
     &                 ,ictwo,ipeak,3,lolimit)
                  call twh_cbmcbendinfo(GLB_GET,imolty,iuprev,icone
     &                 ,ictwo,ipeak,4,hilimit)
               endif
            else
               write(6,*) 'GETCBANGLE: unknown cbmc_bendgen_code'
               write(6,*) cbmc_bendgen_code
               lfailure = .true.
               return
            endif

            if ( .not. lfirstold ) then
               if ( style_code .eq. DIST_UNIFORM_GAUSSIAN ) then
                  rdecide = twh_random()
               endif
c              --- select an angle
               if ( style_code .eq. DIST_UNIFORM ) then
c                 --- uniform distribution
                  call twh_cbphi(GLB_GET,angstart,cbphi)
                  phitwo = twh_random()*twopi - cbphi
               elseif ( style_code .eq. DIST_DELTA ) then
c                 --- sampling of a multipeak delta distribution
                  if ( lrandelta ) then
c                    --- select from one of the options
                     ipeak = idint(twh_random()*dble(npeaks)) + 1
                  else
c                    --- sampling the full list of angles in order
                     ipeak = itrial
                  endif
                  phitwo = delta_angle(ipeak) + delta_shift
               elseif ( style_code .eq. DIST_UNIFORM_GAUSSIAN 
     &              .and. rdecide .le. bend_b_ideal_fraction  ) then
c                 --- uniform distribution, but need gaussian bounds
c                 --- to properly account for the bias
                  call twh_cbphi(GLB_GET,angstart,cbphi)
                  phitwo = twh_random()*twopi - cbphi
                  call twh_angle_b_distribution(style_code,lolimit
     &                 ,hilimit,mean,std_dev,phitwo,imolty,iuprev
     &                 ,icone,ictwo,npeaks,lfailure)
                  if ( lfailure ) return
               elseif ( style_code .eq. DIST_GAUSSIAN .or.
     &                 ( style_code .eq. DIST_UNIFORM_GAUSSIAN
     &                 .and. rdecide .gt. bend_b_ideal_fraction )
     &                 ) then
c                 --- gaussian distribution, mean and stdev set above
                  phitwo = twh_gaussian(lfailure,mean,std_dev
     &                 ,lolimit,hilimit)
                  if ( lfailure ) return
               else
                  write(6,*) 'GETCBANGLE: unknown bend B style'
                  write(6,*) 'style_code:',style_code
                  lfailure = .true.
                  return
               endif
            endif

c           --- assign bias
            if ( style_code .eq. DIST_UNIFORM ) then
               pdf_bias = 1.0d0/twopi
            elseif ( style_code .eq. DIST_GAUSSIAN ) then
               pdf_bias = twh_gaussprob(phitwo,mean,std_dev,lolimit
     &              ,hilimit)/dble(npeaks)
            elseif ( style_code .eq. DIST_UNIFORM_GAUSSIAN ) then
c              --- linear combination of uniform and gaussian
               pdf_bias = bend_b_ideal_fraction/twopi 
     &              + (1.0d0-bend_b_ideal_fraction)
     &              *twh_gaussprob(phitwo,mean,std_dev,lolimit,hilimit)
     &              /dble(npeaks)
            elseif ( style_code .eq. DIST_DELTA ) then
c              --- discrete sampling of every non-infinite energy angle
c              --- every peak is equally likely, use ratio
               pdf_ratio = 1.0d0
            else
               write(6,*) 'GETCBANGLE: unknown bend B style'
               write(6,*) 'style_code:',style_code
               lfailure = .true.
               return
            endif

c           --- continuous true probablity density
            pdf_true = 1.0d0/twopi

c           --- now shift this value relative to the angstart
            call twh_cbphi(GLB_GET,angstart,cbphi)
            phitwo = phitwo + cbphi

c           --- safety check pdf_bias
            if ( style_code .ne. DIST_DELTA 
     &           .and. pdf_bias .le. 0.0d0 ) then
               write(6,*) 'GETCBANGLE: invalid bend B pdf_bias'
               write(6,*) 'lfirstold',lfirstold
               write(6,*) 'pdf_bias:',pdf_bias
               write(6,*) 'phitwo,mean,std_dev',phitwo,mean,std_dev
               write(6,*) 'lolimit,hilimit',lolimit,hilimit
               lfailure = .true.
               return
            endif

c           --- initialize energies
            vvib(itrial) = 0.0d0
            vphi(itrial) = 0.0d0
            vtg(itrial) = 0.0d0

            do aaa = angstart,count-1
               call twh_growlist(GLB_GET,iw,aaa,iuone)
               call twh_cbphi(GLB_GET,aaa,phione)
               call twh_bendang(GLB_GET,aaa,thetaone)
               call twh_bendang(GLB_GET,count,thetatwo)
               call twh_coneangle(thetaone,phione,thetatwo,phitwo
     &              ,anglec)
               angle = twh_arccos(anglec)
c              --- store cross angle for angle-angle computation
               crossang(count,aaa) = anglec
               crossang(aaa,count) = anglec

c              --- fetch distances for 1-3 and bond-angle
               call twh_bondlen(GLB_GET,aaa,distaaa)
               call twh_bondlen(GLB_GET,count,distg)
c              --- find the angle type for iugrow-iufrom-iuone
               do ibend = 1, twh_get_inben(imolty,iugrow)
                  iulast = twh_get_ijben2(imolty,iugrow,ibend)
                  if ( iulast .eq. iufrom ) then
                     iu2back = twh_get_ijben3(imolty,iugrow,ibend)
                     if ( iu2back .eq. iuone ) then
c                       --- calculate the bond angle energy
                        vphi(itrial) = vphi(itrial) 
     &                       + twh_vangle(lfailure,lovrlap,ibox,imolty
     &                       ,iugrow,ibend,0,anglec,distg,distaaa
     &                       ,ddum3,ddum3,ddum3)
                        if ( lfailure ) return
                        if ( lovrlap ) goto 13
                     endif
                  endif
               enddo
               if ( lthreemem ) then
c                 --- 3-mem rings are present, see if iugrow-iuone is
c                 --- a valid bond
                  call twh_invib(GLB_GET,imolty,iugrow,numvib)
                  do ivib = 1,numvib
                     call twh_ijvib(GLB_GET,imolty,iugrow,ivib,itest)
                     if ( itest .eq. iuone ) then
c                       --- bond between growing atoms found
                        testlen = twh_lawofcosine(distg,distaaa,anglec)
                        call twh_itvib(GLB_GET,imolty,iugrow,ivib
     &                       ,vibtype)
                        vvib(itrial) = twh_vbond(lfailure,vibtype
     &                       ,testlen,imolty,iugrow,iuone,lovrlap)
                        if ( lfailure ) return
                        if ( lovrlap ) goto 13
                     endif
                  enddo
               endif
            enddo
c           --- calculate angle-angle cross terms (step = 2)
            call twh_cbmcaatot(GLB_GET,count,2,total)
            do itor = 1, total
               call twh_cbmcaatype(GLB_GET,count,2,itor,ic)
c              --- on this step aalocl(0,ic) cannot be iuprev
               call twh_cbmcaalocl(GLB_GET,0,ic,iczero)
c              --- determine icone
               call twh_cbmcaalocl(GLB_GET,1,ic,testone)
               if ( testone .eq. -1 ) then
c                 --- this is iuprev
                  icone = 0
               else
                  call twh_cbmcaalocl(GLB_GET,1,ic,icone)
               endif
c              --- determine ictwo
               call twh_cbmcaalocl(GLB_GET,2,ic,testtwo)
               if ( testtwo .eq. -1 ) then
c                 --- this is iuprev
                  ictwo = 0
               else
                  call twh_cbmcaalocl(GLB_GET,2,ic,ictwo)
               endif
               cosang(1) = crossang(iczero,icone)
               cosang(2) = crossang(iczero,ictwo)
               vphi(itrial) = vphi(itrial) 
     &              + twh_vangang(AA_ANGLE,ibox,imolty,iufrom,ic,cosang
     &              ,ddum4,ddum4,ddum4)
            enddo
c           --- calculate energy for regular torsions that involve a
c           --- 3-membered ring that we just formed
            if ( lthreemem ) then
               call twh_intor(GLB_GET,imolty,iuprev,numtor)
               do itor = 1,numtor
                  call twh_ijtor2(GLB_GET,imolty,iuprev,itor,testtwo)
                  if ( testtwo .eq. iufrom ) then
c                    --- check that 3 and 4 are growing beads
                     call twh_ijtor3(GLB_GET,imolty,iuprev,itor
     &                    ,testthree)
                     call twh_ijtor4(GLB_GET,imolty,iuprev,itor
     &                    ,testfour)
                     threecode = 0
                     fourcode = 0
                     if ( testthree .eq. iugrow ) then
                        threecode = count
                     elseif ( testfour .eq. iugrow ) then
                        fourcode = count
                     endif
                     if ( threecode + fourcode .ne. 0 ) then
                        do aaa = angstart,count-1
                           call twh_growlist(GLB_GET,iw,aaa,iuone)
                           if ( testthree .eq. iuone ) then
                              threecode = aaa
                           elseif ( testfour .eq. iuone ) then
                              fourcode = aaa
                           endif
                        enddo
                     endif
                     if ( threecode*fourcode .ne. 0 ) then
c                       --- torsion found, compute the torsion from the
c                       --- angles
c                       --- alpha = prev - from - three
                        call twh_bendang(GLB_GET,threecode,ang_alpha)
c                       --- beta = four - from - three
                        ang_beta = crossang(threecode,fourcode)
c                       --- gamma = prev - from - four
                        call twh_bendang(GLB_GET,fourcode,ang_gamma)
                        call twh_angle_b_rigid_cosine(ang_alpha,ang_beta
     &                       ,ang_gamma,phic)
                        cordstyle = -1
                        vtest = twh_vtorsion(lfailure,ibox,cordstyle
     &                       ,imolty,iuprev,itor,xcord,ycord,zcord
     &                       ,phic,lovrlap)
                        if ( lfailure ) return
                        if ( lovrlap ) goto 13
                        vtg(itrial) = vtg(itrial) + vtest
                     endif
                  endif
               enddo
            endif
c           --- calculate energy for improper torsions
            do itor = 1, imptot(count)
               ic = tmpimptype(count,itor)
c              --- assign all of the vectors relative to the iufrom
c              --- as the origin
               xcord(1) = 0.0d0
               ycord(1) = 0.0d0
               zcord(1) = 0.0d0
               do ibend = 2,4
                  call twh_cbmcimproper(GLB_GET,ibend,ic,itemp)
                  if ( itemp .lt. 0 ) then
                     call twh_distprev(GLB_GET,1,icr)
                     icphi = 0.0d0
                     ictheta = 0.0d0
                  else
                     if ( itemp .eq. count ) then
                        icphi = phitwo
                     else
                        call twh_cbphi(GLB_GET,itemp,icphi)
                     endif   
                     call twh_bondlen(GLB_GET,itemp,icr)
                     call twh_bendang(GLB_GET,itemp,ictheta)
                  endif
c                 --- transform from spherical to cartesian coords
                  xcord(ibend) = icr*dsin(ictheta)*dcos(icphi)
                  ycord(ibend) = icr*dsin(ictheta)*dsin(icphi)
                  zcord(ibend) = icr*dcos(ictheta)
               enddo
               vtg(itrial) = vtg(itrial) + twh_vimproper(ibox,imolty
     &              ,iufrom,ic,xcord,ycord,zcord,theta,lbadstereo)
            enddo

c           --- store the boltzmann factors and phi
 13         if ( lovrlap ) then
               lovrben(itrial) = .true.
            else
               lovrben(itrial) = .false.
               vtotal = vvib(itrial) + vphi(itrial) + vtg(itrial)
               if ( lmin ) then
c                 --- see if this energy is lower than the min
                  if ( vtotal .lt. vmin ) vmin = vtotal
               else
c                 --- set vmin
                  vmin = vtotal
                  lmin = .true.
               endif
            endif
c           --- store ang_trial
            ang_trial(itrial) = phitwo
c           --- store the ratio of the pdf factors in bfactor.  We will
c           --- update this with the energy once we have determined the
c           --- minimum energy
            if ( style_code .eq. DIST_DELTA ) then
c              --- using ratio directly
               bfactor(itrial) = pdf_ratio
            else
c              --- include the bias terms
               bfactor(itrial) = pdf_true/pdf_bias
            endif
         enddo
c        --- initialize bsum_try
         bsum_try = 0.0d0
         lovrangle = .true.
         do itrial = 1,nchoice
c           --- compute bsum_try and the adjusted bfactors
            if ( lovrben(itrial) ) then
               bfactor(itrial) = 0.0d0
            else
c              --- non-overlap found
               lovrangle = .false.
c              --- compute boltzmann factor
               vtotal = vvib(itrial) + vphi(itrial) + vtg(itrial) - vmin
               bfactor(itrial) = bfactor(itrial)*twh_expon(-beta*vtotal)
            endif
            bsum_try = bsum_try + bfactor(itrial)
         enddo
c        --- if all overlaps then selection failed
         if ( lovrangle ) return

         if ( .not. luseoldone ) then
c           --- select a value of phitwo in a biased fashion
            rbf = twh_random()*bsum_try
            bs = 0.0d0
            iwinner = 0
            do while ( bs .lt. rbf .and. iwinner .lt. nchoice )
               iwinner = iwinner + 1
               bs = bs + bfactor(iwinner)
            enddo
            if ( bs .le. rbf ) then
c              --- somehow all of the choices had infinite energy
               write(6,*) 'GETCBANGLE: all choices have infinite energy'
               write(6,*) 'for the angle B selection'
               write(6,*) 'molecule type:',imolty
               write(6,*) 'iufrom,iuprev,ntogrow:',iufrom,iuprev,ntogrow
               lfailure = .true.
               return
            endif
         else
c           --- select the OLD value of phitwo
            iwinner = 1
         endif

c        --- account for vmin
#if DEBUG_GETCBANGLE
      write(6,*) 'bend b vmin',vmin
#endif
         call twh_getweight(lfailure,0,vmin,wmin,dpmin,lterm)
         if ( lfailure ) return
         if ( lterm ) then
c           --- this should not happen
            write(6,*) 'GETCBANGLE: getweight returned lterm in bend b'
            lfailure = .true.
            return
         endif
c        --- propagate angle weight and dpcbangle
#if DEBUG_GETCBANGLE
      write(6,*) 'wmin,bsum_try,B nchoice)',wmin,bsum_try,nchoice
#endif
         gcba_weight = gcba_weight*wmin*(bsum_try/dble(nchoice))
         dpcbangle = dpcbangle + dpmin

c        --- store the angle for thetatwo
         cbphi = ang_trial(iwinner)
         call twh_cbphi(GLB_SET,count,cbphi)
c        --- increment energies
         gcba_vvib = gcba_vvib + vvib(iwinner)
         gcba_vbend = gcba_vbend + vphi(iwinner)
         gcba_vtor = gcba_vtor + vtg(iwinner)
c        --- recompute and store the cross angles
         do aaa = angstart,count-1
            call twh_cbphi(GLB_GET,count,phitwo)
            call twh_growlist(GLB_GET,iw,aaa,iuone)
            call twh_cbphi(GLB_GET,aaa,phione)
            call twh_bendang(GLB_GET,aaa,thetaone)
            call twh_bendang(GLB_GET,count,thetatwo)
            call twh_coneangle(thetaone,phione,thetatwo,phitwo
     &           ,anglec)
c           --- store cross angle for angle-angle computation
            crossang(count,aaa) = anglec
            crossang(aaa,count) = anglec
         enddo

      enddo

c     **********************************************************
c     ** End Bond angle biased selection                     **
c     **********************************************************

c      write(6,*) 'FINISH GETCBANGLE'
      return
      end
