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
      subroutine twh_checkhmatrix (ibox,lwarn,lbad)
c     ******************************************************************
c     * checks the vectors of the hmatrix to make sure they are all    *
c     * still between 45 and 135 degrees apart from each other         *
c     * this needs to be true for the assumptions used in mimage to be *
c     * valid                                                          *
c     *                                                                *
c     * originally written 01-03-2002 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_DOTPRODUCT
#include "functions.h"

c     --- variables passed to/from the subroutine
      logical lwarn,lbad
      integer ibox

c     --- local variables
      double precision xone,xtwo,xthr,yone,ytwo,ythr,zone,ztwo,zthr
      double precision angle,anglec,dista,distb
      double precision loang,total,onepi,twopi,minboxdim

c     --- retrieve constants
      call twh_minboxdim(GLB_GET,minboxdim)
      call twh_constant_pi(GLB_GET,onepi)
      call twh_constant_twopi(GLB_GET,twopi)

c     --- set the low angle
c     --- note that the symmetry of the cosine means this test will 
c     --- pass if the angle is between loang and twopi-loang
      loang = 45.0d0*onepi/180.0d0 - 1.0d-10

c     --- set the values
      call twh_hmatrix(GLB_GET,ibox,1,1,xone)
      call twh_hmatrix(GLB_GET,ibox,1,2,yone)
      call twh_hmatrix(GLB_GET,ibox,1,3,zone)
      call twh_hmatrix(GLB_GET,ibox,2,1,xtwo)
      call twh_hmatrix(GLB_GET,ibox,2,2,ytwo)
      call twh_hmatrix(GLB_GET,ibox,2,3,ztwo)
      call twh_hmatrix(GLB_GET,ibox,3,1,xthr)
      call twh_hmatrix(GLB_GET,ibox,3,2,ythr)
      call twh_hmatrix(GLB_GET,ibox,3,3,zthr)

c     --- initialize lbad
      lbad = .false.
c     --- make sure the absolute values of the elements add up to more
c     --- than twice the non-bond cutoff in each cartesian direction
      total = dabs(xone)+dabs(xtwo)+dabs(xthr)
      if ( total .lt. minboxdim ) then
         lbad = .true.
         if ( lwarn ) then
            write(6,*) 'x-dimension absolute value sum ',total
            write(6,*) 'minboxdim ',minboxdim
            write(6,*) 'x-dimensions < minboxdim '
         endif
      endif
      total = dabs(yone)+dabs(ytwo)+dabs(ythr)
      if ( total .lt. minboxdim ) then
         lbad = .true.
         if ( lwarn ) then
            write(6,*) 'y-dimension absolute value sum ',total
            write(6,*) 'minboxdim ',minboxdim
            write(6,*) 'y-dimensions < minboxdim '
         endif
      endif
      total = dabs(zone)+dabs(ztwo)+dabs(zthr)
      if ( total .lt. minboxdim ) then
         lbad = .true.
         if ( lwarn ) then
            write(6,*) 'z-dimension absolute value sum ',total
            write(6,*) 'minboxdim ',minboxdim
            write(6,*) 'z-dimensions < minboxdim '
         endif
      endif
      if ( lbad ) return
c     --- make sure each of the vector lengths is twice minboxdim
      total = twh_distance(xone,yone,zone)
      if ( total .lt. minboxdim ) then
         lbad = .true.
         if ( lwarn ) then
            write(6,*) 'a-dimension vector length ',total
            write(6,*) 'minboxdim ',minboxdim
            write(6,*) 'a-dimensions < minboxdim '
         endif
         return
      endif
      total = twh_distance(xtwo,ytwo,ztwo)
      if ( total .lt. minboxdim ) then
         lbad = .true.
         if ( lwarn ) then
            write(6,*) 'b-dimension vector length ',total
            write(6,*) 'minboxdim ',minboxdim
            write(6,*) 'b-dimensions < minboxdim '
         endif
         return
      endif
      total = twh_distance(xthr,ythr,zthr)
      if ( total .lt. minboxdim ) then
         lbad = .true.
         if ( lwarn ) then
            write(6,*) 'c-dimension vector length ',total
            write(6,*) 'minboxdim ',minboxdim
            write(6,*) 'c-dimensions < minboxdim '
         endif
         return
      endif

c     --- check the angle between vectors one and two
      anglec = twh_dotproduct(xone,yone,zone,xtwo,ytwo,ztwo,dista,distb)
      angle = dacos(anglec)
      if ( angle .lt. loang ) then
         lbad = .true.
         if ( lwarn ) then
            write(6,*) 'hmatrix in box ',ibox
     &           ,' vectors 1 and 2 have a bad angle '
            write(6,*) '   the angle should be between ',loang,' and '
     &           ,twopi-loang
            write(6,*) '   actual value of angle is ',angle
         endif
         return
      endif

c     --- check the angle between vectors two and three
      anglec = twh_dotproduct(xtwo,ytwo,ztwo,xthr,ythr,zthr,dista,distb)
      angle = dacos(anglec)
      if ( angle .lt. loang ) then
         lbad = .true.
         if ( lwarn ) then
            write(6,*) 'hmatrix in box ',ibox
     &           ,' vectors 2 and 3 have a bad angle '
            write(6,*) '  the angle should be between ',loang,' and '
     &           ,twopi-loang
            write(6,*) '  actual value of angle is ',angle
         endif
         return
      endif

c     --- check the angle between vectors one and three
      anglec = twh_dotproduct(xone,yone,zone,xthr,ythr,zthr,dista,distb)
      angle = dacos(anglec)
      if ( angle .lt. loang ) then
         lbad = .true.
         if ( lwarn ) then
            write(6,*) 'hmatrix in box ',ibox
     &           ,' vectors 1 and 3 have a bad angle '
            write(6,*) '   the angle should be between ',loang,' and '
     &           ,twopi-loang
            write(6,*) '   actual value of angle is ',angle
         endif
         return
      endif

      return
      end
