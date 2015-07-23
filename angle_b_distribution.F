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
      subroutine twh_angle_b_distribution(style_code,lolimit,hilimit
     &     ,mean,std_dev,phidisp,imolty,iuprev,icone,ictwo,ngausspeaks
     &     ,lfailure)
c     ******************************************************************
c     * sets the appropriate bounds for the angle_b distribution       *
c     *                                                                *
c     * originally written 09-20-2005 by M.G. Martin                   *
c     * last modified 10-15-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer style_code,imolty,iuprev,icone,ictwo,ngausspeaks
      double precision lolimit,hilimit,mean,std_dev,phidisp
      logical lfailure
c     --- local variables
      logical lfound
      integer ipeak
      double precision twopi,radtodeg

      lfound = .false.
c     --- retrieve constants
      call twh_constant_twopi(GLB_GET,twopi)
      call twh_constant_radtodeg(GLB_GET,radtodeg)

      if ( phidisp .gt. twopi ) then
         phidisp = phidisp - twopi
      elseif ( phidisp .lt. -twopi ) then
         phidisp = phidisp + twopi
      endif

      ipeak = 0
      if ( style_code .eq. DIST_GAUSSIAN
     &     .or. style_code .eq. DIST_UNIFORM_GAUSSIAN ) then
         do while ( .not. lfound .and. ipeak .lt. ngausspeaks )
            ipeak = ipeak + 1
            call twh_cbmcbendinfo(GLB_GET,imolty,iuprev,icone,ictwo
     &           ,ipeak,1,mean)
            call twh_cbmcbendinfo(GLB_GET,imolty,iuprev,icone,ictwo
     &           ,ipeak,2,std_dev)
            call twh_cbmcbendinfo(GLB_GET,imolty,iuprev,icone,ictwo
     &           ,ipeak,3,lolimit)
            call twh_cbmcbendinfo(GLB_GET,imolty,iuprev,icone,ictwo
     &           ,ipeak,4,hilimit)
            if ( lolimit .le. phidisp 
     &           .and. phidisp .lt. hilimit ) then
               lfound = .true.
            elseif ( lolimit+twopi .le. phidisp 
     &              .and. phidisp .lt. hilimit+twopi ) then
               lfound = .true.
               mean = mean + twopi
               lolimit = lolimit + twopi
               hilimit = hilimit + twopi
            elseif ( lolimit - twopi .le. phidisp 
     &              .and. phidisp .lt. hilimit-twopi ) then
               lfound = .true.
               mean = mean - twopi
               lolimit = lolimit - twopi
               hilimit = hilimit - twopi
            endif
         enddo

         if ( .not. lfound ) then
            write(6,*) 'ANGLE_B_DISTRIBUTION: match not found'
            write(6,*) 'style_code,imolty,iuprev,icone,ictwo'
     &           ,style_code,imolty,iuprev,icone,ictwo
            write(6,*) 'phidisp:',phidisp*radtodeg
            write(6,*) 'ngausspeaks:',ngausspeaks
            do ipeak = 1,ngausspeaks
               write(6,*) 'ipeak,lolimit,hilimit,mean',ipeak
     &              ,lolimit*radtodeg,hilimit*radtodeg,mean*radtodeg
            enddo
            lfailure = .true.
            return
         endif
      else
c        --- unknown style
         write(6,*) 'ANGLE_B_DISTRIBUTION: unknown style_code'
         write(6,*) style_code
         lfailure = .true.
         return
      endif

      return
      end

      subroutine twh_angle_b_rigid(ang_alpha,ang_beta,ang_gamma,nvalid
     &     ,validlist)
c     ******************************************************************
c     * computes the number of valid angle-b angles given the values   *
c     * of ang_one, ang_two and ang_three                              *
c     * designed for use with rigid bond angles in CBMC growths        *
c     * ang_alpha and ang_beta are the current-grow-prev angles        *
c     * ang_gamma is the regular angle current-grow-other current      *
c     * Thanks to http://whistleralley.com/polyhedra/derivations.htm   *
c     * for this derivation                                            *
c     *                                                                *
c     * originally written 03-02-2007 by M.G. Martin                   *
c     * last modified 03-17-2009 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_ARCCOS
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision ang_alpha,ang_beta,ang_gamma
      integer nvalid
      double precision validlist
      dimension validlist(CBMC_MAXBENDPEAKS)
c     --- local variables
      double precision angle_one,angle_two,twopi,cosdihedral

c     --- retrieve 2Pi
      call twh_constant_twopi(GLB_GET,twopi)

      call twh_angle_b_rigid_cosine(ang_alpha,ang_beta,ang_gamma
     &     ,cosdihedral)
      angle_one = twh_arccos(cosdihedral)
      
      angle_two = twopi - angle_one
      if ( angle_one .eq. angle_two ) then
         nvalid = 1
         validlist(1) = angle_one
      else
         nvalid = 2
         validlist(1) = angle_one
         validlist(2) = angle_two
      endif

      return
      end

      subroutine twh_angle_b_rigid_cosine(ang_alpha,ang_beta,ang_gamma
     &     ,cosdihedral)
c     ******************************************************************
c     * computes the cosine of the dihedral angle described by 3 bend  *
c     * angles.                                                        *
c     *                                                                *
c     * split from angle_b_rigid 01-30-2011 by M.G. Martin             *
c     * last modified 01-30-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision ang_alpha,ang_beta,ang_gamma,cosdihedral

      cosdihedral = (dcos(ang_gamma) -dcos(ang_alpha)*dcos(ang_beta))
     &     /(dsin(ang_alpha)*dsin(ang_beta))
      return
      end
