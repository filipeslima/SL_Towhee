#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2008 Marcus G. Martin                       *
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
      subroutine twh_dihedral_distribution(lfailure,ldihed_direct
     &     ,ldelta_dist,style_code,npeak,delta_dist,lolimit,hilimit
     &     ,mean,std_dev,phidisp,prob_subregion,lgetdist,imolty,iunit
     &     ,torkey)
c     ******************************************************************
c     * performs a biased dihedral selection                           *
c     *                                                                *
c     * rewritten by M.G. Martin 08-09-2005 from getcbdihed            *
c     * last modified 10-20-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- debug toggle
#define DEBUG_DIHEDRAL_DISTRIBUTION 0
c     --- function declarations
#define FUNCTION_RANDOM
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure,ldihed_direct,lgetdist,ldelta_dist
      integer style_code,imolty,iunit,torkey,npeak
      double precision lolimit,hilimit,mean,std_dev,phidisp
      double precision prob_subregion
      double precision delta_dist
      dimension delta_dist(0:CBMCMAXTOR)
c     --- local variables
c     --- logical scalars
      logical lfound
c     --- character string scalars
      character*5 hyfrom,hyprev
c     --- integer scalars
      integer ipeak,ngausspeaks,cbmc_dihedgen_code,tortype,torstyle
c     --- double precision scalars
      double precision sixthpi,rdecide,onepi,dvalue,rtest,rselect
      double precision twopi,onethird

#if DEBUG_DIHEDRAL_DISTRIBUTION
      write(6,*) 'starting dihedral_distribution'
#endif

c     --- retrieve constants
      call twh_constant_pi(GLB_GET,onepi)
      call twh_constant_twopi(GLB_GET,twopi)
      call twh_constant_onethird(GLB_GET,onethird)
      call twh_cbmc_dihedgen_code(GLB_GET,cbmc_dihedgen_code)

c     --- determine which type of biasing to use to select the phidisp
      if ( ldelta_dist ) then
         call twh_ittor(GLB_GET,imolty,iunit,torkey,tortype)
         call twh_torstyle(GLB_GET,tortype,torstyle)
         if ( torstyle .eq. 21 ) then
c           --- multiple rigid dihedrals, use a delta distribution for
c           --- any value of cbmc_dihedgen_code
            style_code = DIST_DELTA
            call twh_ntorloop(GLB_GET,tortype,npeak)
            do ipeak = 0,npeak
               call twh_torcoeff(GLB_GET,tortype,ipeak
     &              ,delta_dist(ipeak))
            enddo
         else
            write(6,*) 'DIHEDRAL_DISTRIBUTION: unsure how to handle a'
     &           ,' delta distribution for torstyle:',torstyle
            lfailure = .true.
            return
         endif
      elseif ( .not. ldihed_direct .or. 
     &        cbmc_dihedgen_code .eq. CBMC_IDEAL ) then
c        --- uniform distribution
         style_code = DIST_UNIFORM
         prob_subregion = 1.0d0
      elseif ( cbmc_dihedgen_code .eq. CBMC_GLOBAL_GAUSSIAN ) then
c        --- choose a displacement angle from some 
c        --- gaussians, or uniform if we do not know
c        --- what gaussians to use
         call twh_sdevtor(GLB_GET,std_dev)
         call twh_hyfrom(GLB_GET,hyfrom)
         call twh_hyprev(GLB_GET,hyprev)
         if ( hyfrom .eq. 'sp3' .and. hyprev .eq. 'sp3' ) then
c           --- use a gaussian distribution with 3 peaks
            style_code = DIST_GAUSSIAN
            prob_subregion = 1.0d0/3.0d0
c           --- choose which third we are sampling
            if ( lgetdist ) then
c              --- select the proper third based on phidisp
               if ( phidisp .lt. -twopi*onethird ) then
                  rdecide = 0.1d0
               elseif ( phidisp .lt. 0.0d0 ) then
                  rdecide = 0.5d0
               elseif ( phidisp .lt. twopi*onethird ) then
                  rdecide = 0.9d0
               else
                  rdecide = 0.1d0
               endif
            else
c              --- choose a third randomly
               rdecide = twh_random()
            endif
            if ( rdecide .lt. onethird ) then
c              --- a mean of Pi or -Pi (same thing)
c              --- we will choose the positive version as it
c              --- will just wrap around to the negative
               lolimit = twopi*onethird
               hilimit = 2.0*lolimit
               mean = onepi
               if ( lgetdist ) then
c                 --- make sure our mean matches up with 
c                 --- the sign of phidisp
                  if ( phidisp .lt. 0.0d0 ) then
                     mean = mean - twopi
                     lolimit = lolimit - twopi
                     hilimit = hilimit - twopi
                  endif
               endif
            elseif ( rdecide .lt. 2.0*onethird ) then
               mean = -onepi*onethird
               lolimit = -twopi*onethird
               hilimit = 0.0d0
            else
               mean = onepi*onethird
               lolimit = 0.0d0
               hilimit = twopi*onethird
            endif
         elseif ( hyfrom .eq. 'sp2' .and. hyprev .eq. 'sp2' ) then
c           --- use a gaussian distribution with 4 peaks
c           --- normally just the planar ones will be 
c           --- relevent, but if these two sp2 style 
c           --- atoms are connected by a single bond then
c           --- the perpendicular arrangement might be favored.
            style_code = DIST_GAUSSIAN
            prob_subregion = 0.25d0
c           --- choose which quadrant we are sampling
            if ( lgetdist ) then
c              --- select the proper quad based on phidisp
               if ( phidisp .lt. -0.75d0*onepi ) then
                  rdecide = 0.001d0
               elseif ( phidisp .lt. -0.25*onepi ) then
                  rdecide = 0.251d0
               elseif ( phidisp .lt. 0.25*onepi ) then
                  rdecide = 0.501d0
               elseif ( phidisp .lt. 0.75*onepi ) then
                  rdecide = 0.751d0
               else
                  rdecide = 0.001d0
               endif
            else
c              --- choose a quad randomly
               rdecide = twh_random()
            endif
            if ( rdecide .lt. 0.25d0 ) then
c              --- a mean of Pi or -Pi (same thing)
c              --- we will choose the positive version as it
c              --- will just wrap around to the negative
               lolimit = 0.75d0*onepi
               hilimit = 1.25d0*onepi
               mean = onepi
               if ( lgetdist ) then
c                 --- make sure our mean matches up with 
c                 --- the sign of phidisp
                  if ( phidisp .lt. 0.0d0 ) then
                     mean = mean - twopi
                     lolimit = lolimit - twopi
                     hilimit = hilimit - twopi
                  endif
               endif
            elseif ( rdecide .lt. 0.50d0 ) then
c              --- centered at -Pi/2
               mean = -0.5d0*onepi
               lolimit = -0.75d0*onepi
               hilimit = -0.25d0*onepi
            elseif ( rdecide .lt. 0.75d0 ) then
c              --- centered at 0
               mean = 0.0d0
               lolimit = -0.25d0*onepi
               hilimit = 0.25d0*onepi
            else
c              --- centered at Pi/2
               mean = 0.5d0*onepi
               lolimit = 0.25d0*onepi
               hilimit = 0.75d0*onepi
            endif
         elseif ( hyfrom .eq. 'arom' .and. hyprev .eq. 'arom' ) then
c           --- multiple aromatic atoms, either 0 or Pi
            style_code = DIST_GAUSSIAN
            prob_subregion = 0.5d0
c           --- choose which quadrant we are sampling
            if ( lgetdist ) then
c              --- select the proper quad based on phidisp
               if ( phidisp .lt. -0.5d0*onepi 
     &              .or. phidisp .ge. 0.5d0*onepi ) then
                  rdecide = 0.001d0
               else
                  rdecide = 0.501d0
               endif
            else
c              --- choose a section randomly
               rdecide = twh_random()
            endif
            if ( rdecide .lt. 0.5d0 ) then
c              --- a mean of Pi or -Pi (same thing)
c              --- we will choose the positive version as it
c              --- will just wrap around to the negative
               lolimit = 0.5d0*onepi
               hilimit = 1.5d0*onepi
               mean = onepi
               if ( lgetdist ) then
c                 --- make sure our mean matches up with 
c                 --- the sign of phidisp
                  if ( phidisp .lt. 0.0d0 ) then
                     mean = mean - twopi
                     lolimit = lolimit - twopi
                     hilimit = hilimit - twopi
                  endif
               endif
            else
c              --- centered at 0
               mean = 0.0d0
               lolimit = -0.5d0*onepi
               hilimit = 0.5d0*onepi
            endif
         elseif ( (hyfrom .eq. 'sp3' .and. hyprev .eq. 'sp2') .or. 
     &           (hyprev .eq. 'sp3' .and. hyfrom .eq. 'sp2')) then
c           --- use a gaussian distribution with 6 peaks
            style_code = DIST_GAUSSIAN
            prob_subregion = 1.0d0/6.0d0
            sixthpi = onepi/6.0d0
c           --- choose which sixth we are sampling
            if ( lgetdist ) then
c              --- select the proper quad based on phidisp
               if ( phidisp .lt. -5.0*sixthpi ) then
c                 --- first sixth -5/6 Pi to 5/6 Pi
                  rdecide = 0.5d0
               elseif ( phidisp .lt. -0.5*onepi ) then
c                 --- second sixth -1/2 Pi to -5/6 Pi
                  rdecide = 1.5d0
               elseif ( phidisp .lt. -sixthpi ) then
c                 --- third sixth -1/6 Pi to -1/2 Pi
                  rdecide = 2.5d0
               elseif ( phidisp .lt. sixthpi ) then
c                 --- fourth sixth 1/6 Pi to -1/6 Pi
                  rdecide = 3.5d0
               elseif ( phidisp .lt. 0.5d0*onepi ) then
c                 --- fifth sixth 1/2 Pi to 1/6 Pi
                  rdecide = 4.5d0
               elseif ( phidisp .lt. 5.0d0*sixthpi ) then
c                 --- sixth sixth 5/6 Pi to 1/2 Pi
                  rdecide = 5.5d0
               else
c                 --- back to the first sixth
                  rdecide = 0.5d0
               endif
            else
c              --- choose a sixth randomly
               rdecide = 6.0d0*twh_random()
            endif
            if ( rdecide .lt. 1.0d0 ) then
c              --- a mean of Pi or -Pi (same thing)
c              --- we will choose the positive version as it
c              --- will just wrap around to the negative
               lolimit = 5.0d0*sixthpi
               mean = onepi
               hilimit = 7.0d0*sixthpi
               if ( lgetdist ) then
c                 --- make sure our mean matches up with 
c                 --- the sign of phidisp
                  if ( phidisp .lt. 0.0d0 ) then
                     mean = mean - twopi
                     lolimit = lolimit - twopi
                     hilimit = hilimit - twopi
                  endif
               endif
            elseif ( rdecide .lt. 2.0d0 ) then
c              --- centered at -4/6 Pi
               lolimit = -5.0d0*sixthpi
               mean = -4.0d0*sixthpi
               hilimit = -3.0d0*sixthpi
            elseif ( rdecide .lt. 3.0d0 ) then
c              --- centered at -2/6 Pi
               lolimit = -3.0d0*sixthpi
               mean = -2.0d0*sixthpi
               hilimit = -sixthpi
            elseif ( rdecide .lt. 4.0d0 ) then
c              --- centered at 0 Pi
               lolimit = -sixthpi
               mean = 0.0d0
               hilimit = sixthpi
            elseif ( rdecide .lt. 5.0d0 ) then
c              --- centered at 2/6 Pi
               lolimit = sixthpi
               mean = 2.0d0*sixthpi
               hilimit = 3.0d0*sixthpi
            else
c              --- centered at 4/6 Pi
               lolimit = 3.0d0*sixthpi
               mean = 4.0d0*sixthpi
               hilimit = 5.0d0*sixthpi
            endif
         else
c           --- have not come up with a strategy for this
c           --- pattern yet so just use uniform
c                  write(6,*) 'hyfrom ',hyfrom,' hyprev',hyprev
            style_code = DIST_UNIFORM
            prob_subregion = 1.0d0
         endif
      elseif ( cbmc_dihedgen_code .eq. CBMC_AUTOFIT_GAUSSIAN .or.
     &        cbmc_dihedgen_code .eq. CBMC_IDEAL_PLUS_AUTOFIT_GAUSSIAN
     &        ) then
         call twh_cbmctorcount(GLB_GET,imolty,iunit,torkey,ngausspeaks)
c        --- check for uniform distribution kickout flag
         if ( torkey .eq. 0 ) then
c           --- no dihedrals, but some other angles and such
            style_code = DIST_UNIFORM
            prob_subregion = 1.0d0
         elseif ( ngausspeaks .eq. 0 ) then
            style_code = DIST_UNIFORM
            prob_subregion = 1.0d0
         else
            if ( cbmc_dihedgen_code .eq. CBMC_AUTOFIT_GAUSSIAN ) then
               style_code = DIST_GAUSSIAN
            else
               style_code = DIST_UNIFORM_GAUSSIAN
            endif
            if ( lgetdist ) then
c              --- find the location of the phidisp
               ipeak = 0
               lfound = .false.
               do while ( .not. lfound .and. ipeak .lt. ngausspeaks )
                  ipeak = ipeak + 1
                  call twh_cbmctorinfo(GLB_GET,imolty,iunit,torkey
     &                 ,ipeak,1,mean)
                  call twh_cbmctorinfo(GLB_GET,imolty,iunit,torkey
     &                 ,ipeak,3,lolimit)
                  call twh_cbmctorinfo(GLB_GET,imolty,iunit,torkey
     &                 ,ipeak,4,hilimit)
                  if ( lolimit .le. phidisp .and. 
     &                 phidisp .le. hilimit ) then
                     lfound = .true.
                  elseif ( lolimit+twopi .le. phidisp .and.
     &                    phidisp .le. hilimit+twopi ) then
                     lfound = .true.
                     mean = mean + twopi
                     lolimit = lolimit + twopi
                     hilimit = hilimit + twopi
                  endif
               enddo
c              --- assign the rest of the parameters
               call twh_cbmctorinfo(GLB_GET,imolty,iunit,torkey
     &              ,ipeak,2,std_dev)
               call twh_dihedral_peak_weights(GLB_GET,imolty,iunit
     &              ,torkey,ipeak,prob_subregion)
               if ( .not. lfound ) then
                  write(6,*) 'GETCBDIHED: match not found autofit'
                  lfailure = .true.
                  return
               endif
            else
c              --- select a peak at random
               rtest = twh_random()
               rselect = 0.0d0
               ipeak = 0
               do while ( rtest .gt. rselect )
                  ipeak = ipeak + 1
                  call twh_dihedral_peak_weights(GLB_GET,imolty,iunit
     &                 ,torkey,ipeak,dvalue)
                  rselect = rselect + dvalue
               enddo
c              --- we have now seleted a peak, so assign values
               call twh_cbmctorinfo(GLB_GET,imolty,iunit,torkey
     &              ,ipeak,1,mean)
               call twh_cbmctorinfo(GLB_GET,imolty,iunit,torkey
     &              ,ipeak,2,std_dev)
               call twh_cbmctorinfo(GLB_GET,imolty,iunit,torkey
     &              ,ipeak,3,lolimit)
               call twh_cbmctorinfo(GLB_GET,imolty,iunit,torkey
     &              ,ipeak,4,hilimit)
               prob_subregion = dvalue
            endif
         endif
      else
         write(6,*) 'GETCBDIHED: unknown cbmc_dihedgen_code'
         write(6,*) cbmc_dihedgen_code
         lfailure = .true.
         return
      endif

#if DEBUG_DIHEDRAL_DISTRIBUTION
      write(6,*) 'finished dihedral_distribution'
      write(6,*) 'lfailure:',lfailure
      write(6,*) 'ldihed_direct:',ldihed_direct
      write(6,*) 'style_code:',style_code
      write(6,*) 'lolimit:',lolimit
      write(6,*) 'hilimit:',hilimit
      write(6,*) 'mean:',mean
      write(6,*) 'std_dev:',std_dev
      write(6,*) 'phidisp:',phidisp
      write(6,*) 'prob_subregion:',prob_subregion
      write(6,*) 'lgetdist:',lgetdist
      write(6,*) 'imolty:',imolty
      write(6,*) 'iunit:',iunit
      write(6,*) 'torkey:',torkey
#endif

      return
      end
