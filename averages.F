#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2012 Marcus G. Martin                       *
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
      subroutine twh_averages(lfailure,updatetype,keycode,inbox,inmolty
     &     ,value,count)
c     ******************************************************************
c     * this subroutine keeps track of all simulation averages         *
c     *                                                                *
c     * originally written 10-25-2004 by M.G. Martin                   *
c     * last modified 09-21-2012 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_LT
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer updatetype,keycode,inbox,inmolty,count
      double precision value

c     --- local saved variables
      integer bacode
      dimension bacode(MAXBAPROP)
      save bacode
      integer bamolty,babox
      dimension bamolty(MAXBAPROP)
      dimension babox(MAXBAPROP)
      save bamolty,babox
      integer nkey,nblock
      save nkey,nblock
      double precision batotal,bacount
      dimension batotal(MAXBAPROP)
      dimension bacount(MAXBAPROP)
      save batotal,bacount
      double precision oldcount,oldtotal
      dimension oldcount(MAXBAPROP)
      dimension oldtotal(MAXBAPROP)
      save oldcount,oldtotal

c     --- local variables
      logical lfound
      integer ikey,iblock
      double precision tempcount,temptotal,mean,msq,dvalue

c     --- find the proper key for certain updatetype
      ikey = 0
      if ( updatetype .eq. AVG_INCREMENT
     &     .or. updatetype .eq. AVG_EXTRACT_AVERAGE
     &     .or. updatetype .eq. AVG_EXTRACT_BLOCK_VAL
     &     .or. updatetype .eq. AVG_SET_BLOCK_VAL
     &     .or. updatetype .eq. AVG_EXTRACT_BLOCK_AVG
     &     .or. updatetype .eq. AVG_EXTRACT_BLOCK_STD
     &     .or. updatetype .eq. AVG_EXTRACT_COUNT
     &     ) then
c        --- find the proper key
         lfound = .false.
         do while ( ikey .lt. nkey .and. .not. lfound )
            ikey = ikey + 1
            if ( bacode(ikey) .eq. keycode
     &           .and. babox(ikey) .eq. inbox
     &           .and. bamolty(ikey) .eq. inmolty ) lfound = .true.
         enddo
         if ( .not. lfound ) then
c           --- problem, we are incrementing something that has not
c           --- been assigned
            write(6,*) 'AVERAGES: no match found'
            write(6,*) 'updatetype=',updatetype
            write(6,*) 'keycode=',keycode
            write(6,*) 'inbox',inbox
            write(6,*) 'inmolty',inmolty
            lfailure = .true.
            return
         endif
      endif

c     --- do the work, depending on the updatetype
      if ( updatetype .eq. AVG_INITIALIZE ) then
c        --- initialize everything to do with the averages
         nkey = 0
         nblock = 0
         do ikey = 1,MAXBAPROP
            bacode(ikey) = -100
            bamolty(ikey) = -1
            batotal(ikey) = 0.0d0
            bacount(ikey) = 0.0d0
            oldtotal(ikey) = 0.0d0
            oldcount(ikey) = 0.0d0
         enddo

      elseif ( updatetype .eq. AVG_ASSIGN ) then
c        --- make sure this is not already assigned
         do ikey = 1,nkey
            if ( keycode .eq. bacode(ikey) .and. inbox .eq. babox(ikey)
     &           .and. inmolty .eq. bamolty(ikey) ) then
c               --- error, this is already in use
               write(6,*) 'AVERAGES: error requested assignment already'
     &              ,' exists'
               write(6,*) keycode,inmolty,inbox
               write(6,*) bacode(ikey),bamolty(ikey),babox(ikey)
               lfailure = .true.
               return
            endif
         enddo
c        --- increment nkey
         nkey = nkey + 1
c        --- safety check nkey
         if ( nkey .gt. MAXBAPROP ) then
            write(6,*) 'AVERAGES: MAXBAPROP exceeded'
            lfailure = .true.
            return
         endif
c        --- assign properties
         bacode(nkey) = keycode
         babox(nkey) = inbox
         bamolty(nkey) = inmolty
      
      elseif ( updatetype .eq. AVG_INCREMENT ) then
c        --- add information into the accumulators
         batotal(ikey) = batotal(ikey) + value
         bacount(ikey) = bacount(ikey) + dble(count)

      elseif ( updatetype .eq. AVG_EXTRACT_AVERAGE ) then
c        --- return the average
         if ( bacount(ikey) .lt. 0.5d0 ) then
c           --- did not have any increments, return zero
            value = 0.0d0
         else
c           --- divide out the number of counts
            value = batotal(ikey)/bacount(ikey)
         endif

      elseif ( updatetype .eq. AVG_EXTRACT_COUNT ) then
c        --- return the count
         count = idnint(bacount(ikey))

      elseif ( updatetype .eq. AVG_END_OF_BLOCK ) then
c        --- increment nblock
         nblock = nblock + 1
c        --- store the block averages
         do ikey = 1,nkey
c           --- get the value for the current block
            tempcount = bacount(ikey) - oldcount(ikey)
            temptotal = batotal(ikey) - oldtotal(ikey)
            if ( tempcount .gt. 0.5d0 ) then
               dvalue = temptotal/tempcount
            else
               dvalue = 0.0d0
            endif
            call twh_blockvalue(GLB_SET,ikey,nblock,dvalue)
c           --- store the current values as the old for next block
            oldcount(ikey) = bacount(ikey)
            oldtotal(ikey) = batotal(ikey)
         enddo
c        --- return nblock
         count = nblock

      elseif ( updatetype .eq. AVG_SET_BLOCK_VAL ) then
c        --- sets the block value.  this is used in special cases 
c        --- where the quantity is difficult to assign on the fly
         if ( count .gt. nblock .or. count .lt. 1 ) then
            write(6,*) 'AVERAGES: invalid extract block request'
            write(6,*) 'nonexistent block ',count
            lfailure = .true.
            return
         endif
         call twh_blockvalue(GLB_SET,ikey,count,value)

      elseif ( updatetype .eq. AVG_EXTRACT_BLOCK_VAL ) then
c        --- return the data for the specified block and key
c        --- block number is stored in count
         if ( count .gt. nblock .or. count .lt. 1 ) then
            write(6,*) 'AVERAGES: invalid extract block request'
            write(6,*) 'nonexistent block ',count
            lfailure = .true.
            return
         endif
         call twh_blockvalue(GLB_GET,ikey,count,value)

      elseif ( updatetype .eq. AVG_EXTRACT_BLOCK_STD ) then
c        --- return the standard deviation from the block averages
         if ( nblock .eq. 0 ) then
c           --- no blocks, return zero
            value = 0.0d0
         else
            mean = 0.0d0
            msq = 0.0d0
            do iblock = 1,nblock
               call twh_blockvalue(GLB_GET,ikey,iblock,dvalue)
               mean = mean + dvalue
               msq = msq + dvalue**2
            enddo
c           --- divide out the number of blocks
            mean = mean / dble(nblock)
            msq = msq / dble(nblock)

            if (twh_cmp_lt(mean*mean,msq)) then
c              --- standard deviation is the square root of variance
               value = dsqrt( msq - mean*mean )
            else
c              --- values are essentially zero but with the noise
c              --- we have something unphysical.  set to zero
               value = 0.0d0
            endif

         endif

      elseif ( updatetype .eq. AVG_EXTRACT_BLOCK_AVG ) then
c        --- return the mean from the block averages
         if ( nblock .eq. 0 ) then
c           --- no blocks, return zero
            value = 0.0d0
         else
            mean = 0.0d0
            do iblock = 1,nblock
               call twh_blockvalue(GLB_GET,ikey,iblock,dvalue)
               mean = mean + dvalue
            enddo
c           --- divide out the number of blocks
            value = mean / dble(nblock)
         endif

      elseif ( updatetype .eq. AVG_EXTRACT_NBLOCK ) then
c        --- return the number of completed blocks
         count = nblock

      else
c        --- unknown updatetype
         write(6,*) 'AVERAGES: unknown updatetype'
         write(6,*) updatetype
         lfailure = .true.
         return
      endif

      return
      end

