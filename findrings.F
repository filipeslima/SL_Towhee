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
      subroutine twh_findrings( imolty )
c     ******************************************************************
c     * this subroutine searches through the structure of the molecule *
c     * and determines which atoms are a part of a cyclic structure    *
c     * by using a depth first search algorithm                        *
c     * also determines the index of each molecule for use in an       *
c     * interior regrowth                                              *
c     *                                                                *
c     * completely rewritten 10-09-2002 by C.D. Lorenz                 *
c     * last modified 07-19-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit  none
#define FUNCTION_GET_INVIB
#define FUNCTION_GET_LEXIST
#include "functions.h"
c     --- debugging logical
#define DEBUG_FINDRINGS 0
c     --- variables passed to and from the subroutine
      integer imolty

c     --- local variables
c     --- logical scalars
      logical ldone,lringfound
c     --- logical arrays
      logical lring
      dimension lring(NUMAX)
c     --- integer scalars
      integer iunit,startunit,lowringnum,currentunit,nextunit,ivib,itest
      integer nunit,testring,aminonum
c     --- integer arrays
      integer vibpointer,prevunit,visitcount
      dimension vibpointer(NUMAX),prevunit(NUMAX),visitcount(NUMAX)

c     --- other stuff
      integer nrings,ibackbone
      character*4 name

      write(6,'(a,1x,i10)') 
     &     'Determining cyclic subunits for molecule type'
     &     ,imolty

c     --- initialize the number of distinct rings
      nrings = 0
      call twh_nunit(GLB_GET,imolty,nunit)

c     --- loop through all of the atoms in the molecule and for each
c     --- starting atom find all sequences of atoms that lead back to
c     --- the starting atom without retracing their steps.
c     --- similar to a self-avoiding path back to the starting atom
      startunit = 0
      do while ( startunit .lt. nunit )
         startunit = startunit + 1
#if DEBUG_FINDRINGS
         write(6,*) 'startunit:',startunit
#endif
c        --- initialize lexist
         call twh_lexist(GLB_INIT,nunit,.false.)
         do iunit = 1,nunit
c           --- initialize the lring logical
            lring(iunit) = .false.
c           --- initialize the vibpointer
            vibpointer(iunit) = 0
c           --- initialize the prevunit
            prevunit(iunit) = 0
c           --- initialize visitcount
            visitcount(iunit) = 0
         enddo
         if ( twh_get_invib(imolty,startunit) .gt. 1 ) then
            currentunit = startunit
            call twh_lexist(GLB_SET,currentunit,.true.)
            ldone = .false.
         else
c           --- cannot make a ring with an atom that only has 1 bond
            ldone = .true.
            currentunit = 0
         endif
         do while ( .not. ldone )
            vibpointer(currentunit) = vibpointer(currentunit) + 1
            if ( vibpointer(currentunit) .gt. 
     &           twh_get_invib(imolty,currentunit) ) then
c              --- we are finished exploring this atom, kick out
               if ( prevunit(currentunit) .eq. 0 ) then
c                 --- completely finished with this startunit
                  ldone = .true.
               else
c                 --- exit down to the previous atom
                  call twh_lexist(GLB_SET,currentunit,.false.)
                  currentunit = prevunit(currentunit)
               endif
            else
c              --- check to see if this completes a ring
               call twh_ijvib(GLB_GET,imolty,currentunit
     &              ,vibpointer(currentunit),itest)
               if ( itest .eq. prevunit(currentunit) ) then
c                 --- we do not want to explore the bond that leads 
c                 --- directly back to the previous unit
c                 --- nothing to do here and the algorithm will go to
c                 --- the next vibration
               elseif ( itest .eq. startunit .or. lring(itest) ) then
c                 --- the ring is complete, set all atoms in the list
c                 --- to true for ring detection
                  iunit = startunit
                  lring(iunit) = .true.
                  do while ( iunit .ne. currentunit )
                     ivib = vibpointer(iunit)
                     call twh_ijvib(GLB_GET,imolty,iunit,ivib,iunit)
                     lring(iunit) = .true.
                  enddo
               elseif ( twh_get_lexist(itest) ) then
c                 --- if itest is not the start unit and it already
c                 --- exists then we are looping back onto the chain
c                 --- while this clearly has found a cyclic subunit,
c                 --- we are only interested in cyclic subunits that 
c                 --- involve the startunit so we termininate this line
c                 --- of inquiry
               else
c                 --- ring not formed yet, continue search
                  call twh_ijvib(GLB_GET,imolty,currentunit
     &                 ,vibpointer(currentunit),nextunit)
                  prevunit(nextunit) = currentunit
                  currentunit = nextunit
                  call twh_lexist(GLB_SET,currentunit,.true.)
                  vibpointer(currentunit) = 0
                  visitcount(currentunit) = visitcount(currentunit)+1
               endif
            endif
            if ( visitcount(currentunit) .gt. NNBOND*NUMAX ) then
c              --- we must be stuck in some strange periodic ring
c              --- structure.  Kick out and hope we figure it out from
c              --- the next startunit
               ldone = .true.
            endif
         enddo
c        --- find the lowest existing ring number
         lowringnum = nrings + 1
         lringfound = .false.
         do iunit = 1,nunit
            if ( lring(iunit) ) then
               lringfound = .true.
               call twh_inring(GLB_GET,imolty,iunit,testring)
               if ( testring .ne. 0 ) then
                  if ( testring .lt. lowringnum ) then
                     lowringnum = testring
                  endif
               endif
            endif
         enddo
         if ( lringfound ) then
c           --- set all atoms in this list to the lowest ring number
            do iunit = 1,nunit
               if ( lring(iunit) ) then
                  call twh_inring(GLB_SET,imolty,iunit,lowringnum)
               endif
            enddo
c           --- update nrings if required
            if ( lowringnum .gt. nrings ) nrings = nrings + 1
c           --- see if there are any atoms that are not in a ring
            ldone = .true.
            do iunit = 1,nunit
               if ( .not. lring(iunit) ) ldone = .false.
            enddo
            if ( ldone ) startunit = nunit
         endif
      enddo

#if DEBUG_FINDRINGS
      do iunit = 1,nunit
         call twh_inring(GLB_GET,imolty,iunit,testring)
         if ( testring .ne. 0 ) then
            write(6,*) 'iunit,inring:',iunit,testring
         endif
      enddo
#endif

c     --- until I can make this more general I am setting up the list
c     --- of atoms to use in the torsion fixed-endpoint move by using
c     --- the pdbnames to identify the backbone.
      ibackbone = 0
      do iunit = 1,nunit
         call twh_lbackbone(GLB_SET,imolty,iunit,.false.)
         call twh_aminonum(GLB_GET,imolty,iunit,aminonum)
         if ( aminonum .ne. 0 ) then
c           --- we have a protein, set up trickery
            call twh_pdbname(GLB_GET,imolty,iunit,name)
            if ( name .eq. ' C' 
     &           .or. name .eq. ' CA'
     &           .or. name .eq. ' N' ) then
               call twh_lbackbone(GLB_SET,imolty,iunit,.true.)
               ibackbone = ibackbone + 1
               call twh_backlist(GLB_SET,imolty,ibackbone,iunit)
            endif
         endif
      enddo
c     --- save the number of backbone atoms
      call twh_backlist(GLB_SET,imolty,0,ibackbone)

      return
      end
