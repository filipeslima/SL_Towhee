#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2008 Marcus G. Martin                       *
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
      function twh_vfield(lfailure,lovrlap,ibox,ntii,imolty,ichain
     &,iunit,coord)
c     ******************************************************************
c     * computes the energy from all of the external fields with the   *
c     * atom in box ibox                                               *
c     *                                                                *
c     * originally written 11-17-2004 by M.G. Martin                   *
c     * last update 12-06-2008 by M.G. Martin                          *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_NBCOEFF
#define FUNCTION_VFIELD
#include "functions.h"

c     --- variables passed to/from the function
      logical lfailure,lovrlap
      integer ibox,ntii,imolty,ichain,iunit,nbii
      double precision coord
      dimension coord(3)
c     --- local variables
c     --- logical scalars
      logical linfield,lmatch
c     --- character scalars
      character*11 fchar
c     --- integer scalars
      integer ifield,ixyz,nfield,fieldbox,fielddir,fieldtype
c     --- double precision scalars
      double precision dum,dist,dumxx,dumyy,dumzz,invdist,invdist3
      double precision invdist9,distx,disty,distz,xinit,yinit,zinit
      double precision szone,sztwo,szfour,szten,energy
      double precision dpfieldone,dpfieldtwo,dpfieldthree,dpfieldfour
      double precision fnbone,fnbtwo,fnbthr,fnbfou,fconstant
c     --- double precision arrays
      double precision localcoord
      dimension localcoord(3)

      twh_vfield = 0.0d0
      call twh_nfield(GLB_GET,nfield)
      do ifield = 1,nfield
c        --- check box (always fieldint(1) for all fields)
         call twh_fieldint(GLB_GET,ifield,1,fieldbox)
         if ( fieldbox .eq. ibox ) then
c           --- make sure the local coordinate is in the box
            localcoord(1) = coord(1)
            localcoord(2) = coord(2)
            localcoord(3) = coord(3)
            call twh_putcartinbox(ibox,localcoord(1),localcoord(2)
     &           ,localcoord(3))
            call twh_fieldtype(GLB_GET,ifield,fieldtype)
c           --- harmonic attractor field
            if ( fieldtype .eq. FLD_HARMONIC_ATTRACTOR ) then
c              --- harmonic attractor is applied to certain atoms 
c              --- in order to keep them fixed in position relative
c              --- to either initial positions or some pre-entered 
c              --- position
c              --- v = k_fha*(xyz - xyz_0)^2

c              --- see if there is a match for this molecule type and
c              --- atom type
               call twh_vfield_match(ifield,imolty,ntii,lmatch)
               if ( lmatch ) then
                  call twh_fieldchar(GLB_GET,ifield,1,fchar)
                  if ( fchar .eq. 'Global') then
                     call twh_fielddp(GLB_GET,ifield,2,xinit)
                     call twh_fielddp(GLB_GET,ifield,3,yinit)
                     call twh_fielddp(GLB_GET,ifield,4,zinit)
                  elseif ( fchar .eq. 'Initial') then
                     call twh_fieldcoordinates(GLB_GET,ichain,iunit
     &                    ,xinit,yinit,zinit)
                  else
                     write(6,*) 'VFIELD: unknown hafrefpos'
                     lfailure = .true.
                     return
                  endif
                  call twh_mimage(ibox,localcoord(1),localcoord(2)
     &                 ,localcoord(3)
     &                 ,xinit,yinit,zinit,distx,disty,distz)
                  dist = twh_distance(distx,disty,distz)
                  call twh_fielddp(GLB_GET,ifield,1,fconstant)
                  twh_vfield = twh_vfield + (fconstant*(dist*dist))
               endif
            elseif ( fieldtype .eq. FLD_HARD_BALL ) then
c              --- get field parameters.
               call twh_fielddp(GLB_GET,ifield,FLD_HRBPOSX,xinit)
               call twh_fielddp(GLB_GET,ifield,FLD_HRBPOSY,yinit)
               call twh_fielddp(GLB_GET,ifield,FLD_HRBPOSZ,zinit)
               call twh_fielddp(GLB_GET,ifield,FLD_HRBRAD,dpfieldone)
               call twh_fielddp(GLB_GET,ifield,FLD_HRBFAC,dpfieldtwo)
               call twh_fielddp(GLB_GET,ifield,FLD_HRBENERGY,energy)
c              --- hard ball is placed inside the simulation box.
c              --- calculate the distance
               call twh_mimage(ibox
     &           ,localcoord(1),localcoord(2),localcoord(3)
     &           ,xinit,yinit,zinit
     &           ,distx,disty,distz)
               dist = twh_distance(distx,disty,distz)
c              --- calculate atom radius 
               nbii = twh_getnbtype(ntii,ntii)
c              --- perform overlap test
               dum = dpfieldone
     &               + 0.5d0*twh_get_nbcoeff(1,nbii)
               dum = dum * dpfieldtwo
               if ( dabs(dist) .lt. dum ) then
c                 --- check the overlap style
                  call twh_fieldchar(GLB_GET,ifield,1,fchar)
                  if ( fchar .eq. 'infinite' ) then
                     lovrlap = .true.
c                    --- infinite energy, no sense computing other
c                    --- fields
                     return
                  elseif ( fchar .eq. 'finite' ) then
c                    --- apply the penalty energy
                     twh_vfield = twh_vfield + energy
                  endif
               endif
            else
c              --- determine plane of the field
               call twh_fieldchar(GLB_GET,ifield,1,fchar)
               if ( fchar .eq. 'x' ) then
                  ixyz = 1
               elseif ( fchar .eq. 'y' ) then
                  ixyz = 2
               elseif ( fchar .eq. 'z' ) then
                  ixyz = 3 
               else
                  write(6,*) 'VFIELD: unknown fieldchar',fchar
                  write(6,*) 'field ',ifield
                  lfailure = .true.
                  return
               endif
c              --- determine the type of the field
               if ( fieldtype .eq. FLD_HARD_WALL ) then
                  dum = 0.0d0
                  call twh_fielddp(GLB_GET,ifield,1,dpfieldone)
                  if ( ixyz .eq. 1 ) then
c                    --- xcoord field
                     call twh_mimage(ibox,localcoord(ixyz),dum,dum
     &                    ,dpfieldone,dum,dum,dist,dumyy,dumzz)
                  elseif ( ixyz .eq. 2 ) then
c                    --- ycoord field
                     call twh_mimage(ibox,dum,localcoord(ixyz),dum
     &                    ,dum,dpfieldone,dum,dumxx,dist,dumzz)
                  elseif ( ixyz .eq. 3 ) then
c                    --- zcoord field
                     call twh_mimage(ibox,dum,dum,localcoord(ixyz)
     &                    ,dum,dum,dpfieldone,dumxx,dumyy,dist)
                  endif
                  linfield = .false.
                  call twh_fieldchar(GLB_GET,ifield,2,fchar)
                  if ( fchar .eq. 'centers' ) then
                     call twh_fielddp(GLB_GET,ifield,2,dpfieldtwo)
                     if ( dabs(dist) .lt. dpfieldtwo ) then
                        linfield = .true.
                     endif
                  elseif ( fchar .eq. 'hard radii' ) then
                     nbii = twh_getnbtype(ntii,ntii)
                     call twh_fielddp(GLB_GET,ifield,2,dpfieldtwo)
                     if ( dabs(dist) .lt. dpfieldtwo
     &                    + 0.5d0*twh_get_nbcoeff(1,nbii) ) then
                        linfield = .true.
                     endif
                  else
                     write(6,*) 'VFIELD: unknown hrd_repulsion_style '
     &                    ,fchar
                     lfailure = .true.
                     return
                  endif

                  if ( linfield ) then
                     call twh_fieldchar(GLB_GET,ifield,3,fchar)
                     if ( fchar .eq. 'infinite' ) then
                        lovrlap = .true.
c                       --- infinite energy, no sense computing other
c                       --- fields
                        return
                     elseif ( fchar .eq. 'finite' ) then
c                       --- set to the designated finite energy
                        call twh_fielddp(GLB_GET,ifield,4,dpfieldfour)
                        twh_vfield = twh_vfield + dpfieldfour
                     else
                        write(6,*) 'VFIELD: unknown hrd_field_type'
                        write(6,*) fchar
                        lfailure = .true.
                        return
                     endif
                  endif
               elseif ( fieldtype .eq. FLD_LJ_9_3_WALL ) then
c                 --- interaction energy of an atom with a 9-3 LJ wall.
c                 v = (2 pi eps_wf rho_wall sig_wf^3)/3 
c                     * [ 2/15 (sig_wf/z)^9 - (sig_wf/z)^3 ]
c                 ---  computed here as
c                 v = ljfrep(ntii)/(z^9) + ljfatr(ntii)/(z^3)
c                 --- the prefactors were computed in readclassical
                  call twh_fieldint(GLB_GET,ifield,2,fielddir)
                  call twh_fielddp(GLB_GET,ifield,1,dpfieldone)
                  if ( fielddir  .eq. -1 ) then
c                    --- particle interacts only with left face of wall
                     dist = dpfieldone - localcoord(ixyz)
                  else
c                    --- particle interacts only with right face of wall
                     dist = localcoord(ixyz) - dpfieldone
                  endif

c                 --- negative distances do not interact with this wall
                  call twh_fielddp(GLB_GET,ifield,2,dpfieldtwo)
                  if ( dist .ge. 0.0d0 .and. dist .le. dpfieldtwo ) then
                     invdist = 1.0d0/dist
                     invdist3 = invdist*invdist*invdist
                     invdist9 = invdist3*invdist3*invdist3
                     call twh_fieldnbparam(GLB_GET,ifield,ntii,1,fnbone)
                     call twh_fieldnbparam(GLB_GET,ifield,ntii,2,fnbtwo)
                     call twh_fieldnbparam(GLB_GET,ifield,ntii,3,fnbthr)
                     twh_vfield = twh_vfield + (fnbone*invdist9)
     &                    + (fnbtwo*invdist3) + fnbthr
                  endif
               elseif ( fieldtype .eq. FLD_HOOPER_UMBRELLA ) then
c                 --- computes the interaction energy of an atom with a 
c                 --- 4th power "umbrella" potential
c                 v(d) = umba * [ (d - umbcenter)/umbcenter) ]^4
                  dist = localcoord(ixyz)
                  call twh_fielddp(GLB_GET,ifield,1,dpfieldone)
                  call twh_fielddp(GLB_GET,ifield,2,dpfieldtwo)
                  twh_vfield = twh_vfield + dpfieldtwo
     &                 * ((dist-dpfieldone)/dpfieldone)**4
               elseif ( fieldtype .eq. FLD_STEELE_WALL ) then
c                 --- interaction energy of an atom with a 10-4 Steele
c                 --- wall
c                 v = eps_w * [ 2/5 (sig_wf/z)^10 - (sig_wf/z)^4
c                     - sig_wf^4 / 3 delta (z + 0.61 delta)^3 ]
c                 ---  computed here as
c                 v = nb(1) * [ 2/5 (nb(2)/z)^10 - (nb(2)/z)^4
c                     - nb(3) / 3 fdp(3) (z + 0.61 fdp(3))^3 ]
c                 --- the prefactors were computed in readclassical
                  call twh_fieldint(GLB_GET,ifield,2,fielddir)
                  call twh_fielddp(GLB_GET,ifield,1,dpfieldone)
                  if ( fielddir .eq. -1 ) then
c                    --- particle interacts only with left face of wall
                     dist = dpfieldone - localcoord(ixyz)
                  else
c                    --- particle interacts only with right face of wall
                     dist = localcoord(ixyz) - dpfieldone
                  endif
c                 --- exactly zero distance is overlap
                  if ( dist .eq. 0.0d0 ) then
                     lovrlap = .true.
                     return
                  endif
c                 --- negative distances do not interact with this wall
                  call twh_fielddp(GLB_GET,ifield,2,dpfieldtwo)
                  call twh_fielddp(GLB_GET,ifield,3,dpfieldthree)
                  if ( dist .gt. 0.0d0 .and. 
     &                 dist .le. dpfieldtwo ) then
                     call twh_fieldnbparam(GLB_GET,ifield,ntii,2,fnbtwo)
                     szone = fnbtwo/dist
                     sztwo = szone*szone
                     szfour = sztwo*sztwo
                     szten = szfour*szfour*sztwo
                     call twh_fieldnbparam(GLB_GET,ifield,ntii,1,fnbone)
                     call twh_fieldnbparam(GLB_GET,ifield,ntii,3,fnbthr)
                     call twh_fieldnbparam(GLB_GET,ifield,ntii,4,fnbfou)
                     twh_vfield = twh_vfield + 
     &                    fnbone*(
     &                    (2.0d0/5.0d0)*szten - szfour - fnbthr/
     &                    ( 3.0d0*dpfieldthree
     &                    *(dist + 0.61d0*dpfieldthree)**3))
     &                    + fnbfou
                  
                  endif                           
               else
                  write(6,*) 'VFIELD: unknown fieldtype'
                  write(6,*) fieldtype
                  lfailure = .true.
                  return
               endif
            endif
         endif
      enddo

      return
      end

      subroutine twh_vfield_match(ifield,imolty,ntii,lmatch)
c     ******************************************************************
c     * checks the fieldint data structure to see if there is a match  *
c     * for both molecule type and atom type                           *
c     *                                                                *
c     * originally written 12-02-2007 by M.G. Martin                   *
c     * last modified 12-02-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer ifield,imolty,ntii
      logical lmatch
c     --- local variables
      integer ientry,ivalue,testmolty,testtype

c     --- get the total number of possible values
      call twh_fieldint(GLB_GET,ifield,2,ivalue)
      do ientry = 1,ivalue
         call twh_fieldnbint(GLB_GET,ifield,ientry,1,testmolty)
         if ( testmolty .eq. imolty ) then
            call twh_fieldnbint(GLB_GET,ifield,ientry,2,testtype)
            if ( testtype .eq. ntii ) then
               lmatch = .true.
               return
            endif
         endif
      enddo

c     --- if we got this far then there was no match
      lmatch = .false.
      return
      end
