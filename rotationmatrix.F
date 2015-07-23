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
      subroutine twh_rotationmatrix(rmoveflag,ichain,rxorig,ryorig
     &     ,rzorig,dgammax,dgammay,dgammaz)
c     ******************************************************************
c     * rotates the coordinates of a molecule about the three axes     *
c     * which are parallel to the cartesian coordinate vectors         *
c     * Updates the global coordinates for the flag rmoveflag          *
c     * None of the passed arguments are modified                      *
c     * r*orig give the coordinates of the origin (often, COM) about   *
c     *  which the rotation takes place.                               *
c     * dgamma* are the euler angles about which rotation occurs       *
c     *                                                                *
c     * originally rewritten from rotate 05-08-2002 by M.G. Martin     *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- debug toggle
#define DEBUG_ROTATIONMATRIX 0
c     --- variables passed to/from the subroutine
      integer rmoveflag,ichain
      double precision rxorig,ryorig,rzorig,dgammax,dgammay,dgammaz
c     --- local variables
      integer iunit,totunit,imolty
      double precision cosdg,sindg,rx,ry,rz,rxnew,rynew,rznew

#if DEBUG_ROTATIONMATRIX
      write(6,*) 'starting rotationmatrix'
#endif

c     --- get molecule type
      call twh_moltyp(GLB_GET,ichain,imolty)
c     --- set number of units
      call twh_nunit(GLB_GET,imolty,totunit)

c     --- rotate around the x-axis
c     --- set up the rotation marix
      cosdg = dcos(dgammax)
      sindg = dsin(dgammax)
      do iunit = 1, totunit
c        --- fetch current coordinates
         call twh_coordinates(GLB_GET,rmoveflag,ichain,iunit
     &        ,rx,ry,rz)
         ry = ry - ryorig
         rz = rz - rzorig
         rxnew = rx
         rynew = cosdg * ry + sindg * rz + ryorig
         rznew = cosdg * rz - sindg * ry + rzorig
c        --- set new coordinates
         call twh_coordinates(GLB_SET,rmoveflag,ichain,iunit
     &        ,rxnew,rynew,rznew)
      enddo
C     --- rotate around the y-axis
c     --- set up the rotation marix
      cosdg = dcos(dgammay)
      sindg = dsin(dgammay)
      do iunit = 1,totunit
c        --- fetch current coordinates
         call twh_coordinates(GLB_GET,rmoveflag,ichain,iunit
     &        ,rx,ry,rz)
         rx = rx - rxorig
         rz = rz - rzorig
         rxnew = cosdg * rx - sindg * rz + rxorig
         rynew = ry
         rznew = cosdg * rz + sindg * rx + rzorig
c        --- set new coordinates
         call twh_coordinates(GLB_SET,rmoveflag,ichain,iunit
     &        ,rxnew,rynew,rznew)
      enddo
C     --- rotate around the z-axis
c     --- set up the rotation marix
      cosdg = dcos(dgammaz)
      sindg = dsin(dgammaz)
      do iunit = 1,totunit
c        --- fetch current coordinates
         call twh_coordinates(GLB_GET,rmoveflag,ichain,iunit
     &        ,rx,ry,rz)
         rx = rx - rxorig
         ry = ry - ryorig
         rxnew = cosdg * rx + sindg * ry + rxorig
         rynew = cosdg * ry - sindg * rx + ryorig
         rznew = rz
c        --- set new coordinates
         call twh_coordinates(GLB_SET,rmoveflag,ichain,iunit
     &        ,rxnew,rynew,rznew)
      enddo

#if DEBUG_ROTATIONMATRIX
      do iunit = 1,totunit
         call twh_coordinates(GLB_GET,rmoveflag,ichain,iunit,rx,ry,rz)
         write(6,*) 'trial iunit r*:',iunit,rx,ry,rz
      enddo
#endif

#if DEBUG_ROTATIONMATRIX
      write(6,*) 'finished rotationmatrix'
#endif

      return
      end
