#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2011 Marcus G. Martin                       *
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
      subroutine twh_buildhelix(lfailure,imolty,testcord)
c     ******************************************************************
c     * places certain key atoms onto a helix for molecule type imolty *
c     * stores the results in the coordinates according to testcord    *
c     *                                                                *
c     * originally written 08-03-2004 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_ELEMENTNAME
#define FUNCTION_GET_NBNAME
#include "functions.h"
c     --- debug toggle
#define DEBUG_BUILDHELIX 0
c     ---- variables passed to/from the subroutine
      logical lfailure
      integer imolty,testcord
c     --- local variables
c     --- logical scalars
      logical lfound 
c     --- character string scalars
      character*2 key_element
      character*4 key_pdbname,pdbname
      character*10 key_nbname,helix_keyname,helix_keytype
c     --- integer scalars
      integer count,iunit,ntii,ichain,nunit
c     --- double precision scalars
      double precision helix_radius,helix_angle,helix_conlen,helix_phase
c     --- double precision arrays
      double precision positions
      dimension positions(3)

#if DEBUG_BUILDHELIX
      write(6,*) 'start BUILDHELIX'
#endif

c     --- ichain does not matter, just set to zero
      ichain = 0
      count = 0
c     --- set the key
      call twh_helix_keyname(GLB_GET,imolty,helix_keyname)
      call twh_helix_keytype(GLB_GET,imolty,helix_keytype)
      if ( helix_keytype .eq. 'element' ) then
         key_element = helix_keyname(1:2)
      elseif ( helix_keytype .eq. 'pdbname' ) then
         key_pdbname = helix_keyname(1:4)
      elseif ( helix_keytype .eq. 'nbname' ) then
         key_nbname = helix_keyname
      else
         write(6,*) 'BUILDHELIX: unknown helix_keyname for type',imolty
         write(6,*) 'keyname =',helix_keyname,'='
         lfailure = .true.
         return
      endif

      call twh_nunit(GLB_GET,imolty,nunit)
      do iunit = 1,nunit
         call twh_ntype(GLB_GET,imolty,iunit,ntii)
         lfound = .false.
         if ( helix_keytype .eq. 'element' ) then
c           --- test against the element
            if ( key_element .eq. twh_get_elementname(ntii) ) then
               lfound = .true.
            endif
         elseif ( helix_keytype .eq. 'pdbname' ) then
c           --- test against the pdbname
            call twh_pdbname(GLB_GET,imolty,iunit,pdbname)
            if ( key_pdbname .eq. pdbname ) lfound = .true.
         elseif ( helix_keytype .eq. 'nbname' ) then
c           --- test against the atomname
            if ( key_nbname .eq. twh_get_nbname(ntii) ) lfound = .true.
         else
            write(6,*) 'BUILDHELIX: unknown helix_keytype',imolty
            write(6,*) 'keyname =',helix_keytype,'='
            lfailure = .true.
            return
         endif

         if ( lfound ) then
c           --- we have a match, increment count
            count = count + 1
            call twh_helix_radius(GLB_GET,imolty,helix_radius)
            call twh_helix_angle(GLB_GET,imolty,helix_angle)
            call twh_helix_conlen(GLB_GET,imolty,helix_conlen)
            call twh_helix_phase(GLB_GET,imolty,helix_phase)
            call twh_helixpos(lfailure,helix_radius,helix_angle
     &           ,helix_conlen,helix_phase,count,positions)
            if ( lfailure ) return
            call twh_coordinates(GLB_SET,testcord,ichain,iunit
     &           ,positions(1),positions(2),positions(3))
            call twh_cbmc_atomname_list(GLB_SET,imolty,count,iunit)
#if DEBUG_BUILDHELIX
            write(6,*) 'count,iunit,positions',count,iunit,positions
#endif
         endif
      enddo
      call twh_cbmc_atomname_tot(GLB_SET,imolty,count)
c     --- make sure we had at least 1 atom
      if ( count .eq. 0 ) then
         write(6,*) 'BUILDHELIX: no matches found'
         write(6,*) 'Molecule type:',imolty
         write(6,*) 'helix_keyname=',helix_keyname,'='
         lfailure = .true.
         return
      endif

#if DEBUG_BUILDHELIX
      write(6,*) 'finished BUILDHELIX'
#endif

      return
      end
