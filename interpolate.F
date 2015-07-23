#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2013 Marcus G. Martin                       *
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
      function twh_extractpair(lfailure,interpflag,cross,xin,lovrlap)
c     ******************************************************************
c     * extrapolates the eam pair potential for a cross type and dist  *
c     * returns table_pair value computed with cubic spline            *
c     *                                                                *
c     * originally taken from Peter Schultz (Quest) on 11-18-2002      *
c     * last modified 07-15-2015 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EXTRACTPAIR
#include "functions.h"
c     --- variables passed to/from the function
      logical lfailure,lovrlap
      integer cross,interpflag
      double precision xin
c     --- local variables
      integer ntot,klo,khi,k
      character*20 interpolatestyle
      double precision xmin,xmax,dx,a,b
      double precision tpairlo1,tpairhi1,tpairlo2,tpairhi2,tpairlo3
      double precision tpairhi3,tpair

      call twh_table_npair(GLB_GET,cross,ntot)
      call twh_table_pair(GLB_GET,1,cross,1,xmin)
      call twh_table_pair(GLB_GET,1,cross,ntot,xmax)

c     --- out of range test
      if( xin .ge. xmax )then
         call twh_table_pair(GLB_GET,2,cross,ntot,twh_extractpair)
      elseif( xin .lt. xmin )then
         twh_extractpair=0.0d0
         lovrlap = .true.
      else
c        --- in range do binary search
         klo = 1
         khi = ntot
         do while ( khi-klo .gt. 1 )
            k = (khi+klo)/2
            call twh_table_pair(GLB_GET,1,cross,k,tpair)
            if( tpair .gt. xin )then
               khi = k
            else
               klo = k
            endif
         enddo
         call twh_interpolatestyle(GLB_GET,interpolatestyle)
         if ( interpolatestyle .eq. 'cubicspline' ) then
c           --- use a cubic spline
            call twh_table_pair(GLB_GET,1,cross,khi,tpairhi1)
            call twh_table_pair(GLB_GET,1,cross,klo,tpairlo1)
            dx = tpairhi1 - tpairlo1
            a = (tpairhi1 - xin)/dx
            b = 1.d0 - a
            call twh_table_pair(GLB_GET,2,cross,klo,tpairlo2)
            call twh_table_pair(GLB_GET,2,cross,khi,tpairhi2)
            call twh_table_pair(GLB_GET,3,cross,klo,tpairlo3)
            call twh_table_pair(GLB_GET,3,cross,khi,tpairhi3)

            call twh_interpolatecubic(lfailure,interpflag,a,b
     &           ,dx,tpairlo2,tpairhi2,tpairlo3,tpairhi3
     &           ,twh_extractpair)
            if ( lfailure ) return

         elseif ( interpolatestyle .eq. 'linear' ) then
c           --- linear interpolation
            call twh_table_pair(GLB_GET,1,cross,khi,tpairhi1)
            call twh_table_pair(GLB_GET,1,cross,klo,tpairlo1)
            call twh_table_pair(GLB_GET,2,cross,khi,tpairhi2)
            call twh_table_pair(GLB_GET,2,cross,klo,tpairlo2)

            call twh_interpolatelinear(lfailure,interpflag,xin
     &           ,tpairlo1,tpairhi1,tpairlo2,tpairhi2,twh_extractpair)
            if ( lfailure ) return

         else
            write(6,*) 'EXTRACTPAIR: bad interpolatestyle'
     &           ,interpolatestyle
            twh_extractpair=0.0d0
            lfailure = .true.
            return
         endif

      endif

      return
      end

      function twh_extractdens(lfailure,interpflag,itype,jtype,xin)
c     ******************************************************************
c     * extrapolates the eam density for an atom type and dist         *
c     * returns eam_den value computed with cubic spline               *
c     *                                                                *
c     * originally taken from Peter Schultz (Quest) on 11-18-2002      *
c     * last modified 07-15-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EXTRACTDENS
#include "functions.h"
c     --- variables passed to/from the function
      logical lfailure
      integer itype,jtype,interpflag
      double precision xin
c     --- local variables
      character*20 interpolatestyle
      integer ntot,klo,khi,k
      double precision xmin,xmax,dx,a,b,test
      double precision hione,loone,hitwo,lotwo,hithr,lothr

      call twh_eam_ndens(GLB_GET,itype,jtype,ntot)
      call twh_eam_dens(GLB_GET,1,itype,jtype,1,xmin)
      call twh_eam_dens(GLB_GET,1,itype,jtype,ntot,xmax)

c     --- out of range test
      if( xin .ge. xmax )then
         call twh_eam_dens(GLB_GET,2,itype,jtype,ntot,twh_extractdens)
      elseif( xin .lt. xmin )then
         call twh_eam_dens(GLB_GET,2,itype,jtype,1,twh_extractdens)
      else
c        --- in range do binary search
         klo = 1
         khi = ntot
         do while ( khi-klo .gt. 1 )
            k = (khi+klo)/2
            call twh_eam_dens(GLB_GET,1,itype,jtype,k,test)
            if ( test .gt. xin )then
               khi = k
            else
               klo = k
            endif
         enddo
         call twh_interpolatestyle(GLB_GET,interpolatestyle)
         if ( interpolatestyle .eq. 'cubicspline' ) then
c           --- use a cubic spline
            call twh_eam_dens(GLB_GET,1,itype,jtype,khi,hione)
            call twh_eam_dens(GLB_GET,1,itype,jtype,klo,loone)
            dx = hione - loone
            a = (hione - xin)/dx
            b = 1.d0 - a
            call twh_eam_dens(GLB_GET,2,itype,jtype,klo,lotwo)
            call twh_eam_dens(GLB_GET,2,itype,jtype,khi,hitwo)
            call twh_eam_dens(GLB_GET,3,itype,jtype,klo,lothr)
            call twh_eam_dens(GLB_GET,3,itype,jtype,khi,hithr)

            call twh_interpolatecubic(lfailure,interpflag,a,b
     &           ,dx,lotwo,hitwo,lothr,hithr,twh_extractdens)
            if ( lfailure ) return

         elseif ( interpolatestyle .eq. 'linear' ) then
c           --- linear interpolation
            call twh_eam_dens(GLB_GET,1,itype,jtype,khi,hione)
            call twh_eam_dens(GLB_GET,1,itype,jtype,klo,loone)
            call twh_eam_dens(GLB_GET,2,itype,jtype,khi,hitwo)
            call twh_eam_dens(GLB_GET,2,itype,jtype,klo,lotwo)

            call twh_interpolatelinear(lfailure,interpflag,xin
     &           ,loone,hione,lotwo,hitwo,twh_extractdens)
            if ( lfailure ) return

         else
            write(6,*) 'EXTRACTDENS: invalid interpolatestyle'
     &           ,interpolatestyle
            lfailure = .true.
            return
         endif

      endif
c
      return
      end

      function twh_extractembed(lfailure,interpflag,itype,xin)
c     ******************************************************************
c     * extrapolates the eam embeding energy for itype and dist        *
c     * returns eam_embed value computed with cubic spline             *
c     *                                                                *
c     * originally taken from Peter Schultz (Quest) on 11-18-2002      *
c     * last modified 07-15-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EXTRACTEMBED
#include "functions.h"
c     --- variables passed to/from the function
      logical lfailure
      integer itype,interpflag
      double precision xin
c     --- local variables
      character*20 interpolatestyle
      integer ntot,klo,khi,k
      double precision xmin,xmax,dx,a,b
      double precision embedhione,embedhitwo,embedloone,embedlotwo
      double precision embedhithr,embedlothr

      call twh_eam_nembed(GLB_GET,itype,ntot)
      call twh_eam_embed(GLB_GET,1,itype,1,xmin)
      call twh_eam_embed(GLB_GET,1,itype,ntot,xmax)

c     --- out of range test
      if( xin .ge. xmax )then
         call twh_eam_embed(GLB_GET,2,itype,ntot,twh_extractembed)
      elseif( xin .lt. xmin )then
         call twh_eam_embed(GLB_GET,2,itype,1,twh_extractembed)
      else
c        --- in range do binary search
         klo = 1
         khi = ntot
         do while ( khi-klo .gt. 1 )
            k = (khi+klo)/2
            call twh_eam_embed(GLB_GET,1,itype,k,embedhione)
            if( embedhione .gt. xin )then
               khi = k
            else
               klo = k
            endif
         enddo
         call twh_interpolatestyle(GLB_GET,interpolatestyle)
         if ( interpolatestyle .eq. 'cubicspline' ) then
c           --- interpolate using a cubic spline
            call twh_eam_embed(GLB_GET,1,itype,khi,embedhione)
            call twh_eam_embed(GLB_GET,1,itype,klo,embedloone)
            dx = embedhione - embedloone
            a = (embedhione - xin)/dx
            b = 1.d0 - a
            call twh_eam_embed(GLB_GET,2,itype,klo,embedlotwo)
            call twh_eam_embed(GLB_GET,2,itype,khi,embedhitwo)
            call twh_eam_embed(GLB_GET,3,itype,klo,embedlothr)
            call twh_eam_embed(GLB_GET,3,itype,khi,embedhithr)

            call twh_interpolatecubic(lfailure,interpflag,a,b
     &           ,dx,embedlotwo,embedhitwo,embedlothr,embedhithr
     &           ,twh_extractembed)
            if ( lfailure ) return

         elseif ( interpolatestyle .eq. 'linear' ) then
c           --- linear interpolation
            call twh_eam_embed(GLB_GET,1,itype,khi,embedhione)
            call twh_eam_embed(GLB_GET,1,itype,klo,embedloone)
            call twh_eam_embed(GLB_GET,2,itype,khi,embedhitwo)
            call twh_eam_embed(GLB_GET,2,itype,klo,embedlotwo)

            call twh_interpolatelinear(lfailure,interpflag,xin
     &           ,embedloone,embedhione,embedlotwo,embedhitwo
     &           ,twh_extractembed)
            if ( lfailure ) return
         else
            write(6,*) 'EXTRACTEMBED: invalid interpolatestyle'
     &           ,interpolatestyle
            lfailure = .true.
            return
         endif

      endif

      return
      end

      subroutine twh_interpolatecubic(lfailure,interpflag,a,b
     &     ,dx,lotwo,hitwo,lothree,hithree,dvalue)
c     ******************************************************************
c     * performs either a normal or derivative interpolation using a   *
c     * cubic spline.                                                  *
c     *                                                                *
c     * Split out of the individual eam routines 07-15-2013            *
c     * last modified 07-15-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer interpflag
      double precision a,b,dx,lotwo,hitwo,lothree,hithree,dvalue
c     --- local variables
      double precision c,d,dxsq

      dxsq = (dx**2)/6.d0
      if ( interpflag .eq. INTERP_NORMAL ) then
         c = (a**3 - a)*dxsq
         d = (b**3 - b)*dxsq
         dvalue = a*lotwo + b*hitwo + c*lothree + d*hithree
      elseif ( interpflag .eq. INTERP_DERIV ) then
         c = 3.0d0*(a**2 - 1.0d0)*dxsq
         d = 3.0d0*(b**2 - 1.0d0)*dxsq
         dvalue = (-lotwo + hitwo - c*lothree + d*hithree ) / dx
      else
         write(6,*) 'TWH_INTERPOLATECUBIC: unknown interpflag'
     &        ,interpflag
         lfailure = .true.
         return
      endif

      return
      end

      subroutine twh_interpolatelinear(lfailure,interpflag,xin
     &     ,loone,hione,lotwo,hitwo,dvalue)
c     ******************************************************************
c     * linear interpolation, either normal or derivatives             *
c     *                                                                *
c     * split out of eam routines 07-15-2013 by M.G. Martin            *
c     * last modified 07-15-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer interpflag
      double precision xin,loone,hione,lotwo,hitwo,dvalue
c     --- local variables
      double precision dx,dy,b

      dx = hione - loone
      if ( interpflag .eq. INTERP_NORMAL ) then
         b = (xin - loone)/dx
         dvalue = b*(hitwo - lotwo) + lotwo
      elseif ( interpflag .eq. INTERP_DERIV ) then
         dy = hitwo - lotwo
         dvalue = dy/dx
      else
         write(6,*) 'TWH_INTERPOLATELINEAR: unknown interpflag'
     &        ,interpflag
         lfailure = .true.
         return
      endif
      return
      end
