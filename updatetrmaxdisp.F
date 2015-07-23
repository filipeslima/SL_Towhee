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
      subroutine twh_updatetrmaxdisp()
c     ******************************************************************
c     * updated the maximum translation and rotational displacements   *
c     *                                                                *
c     * split out of the pupdate subroutine 10-12-2005 by M.G. Martin  *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- local variables
      character*20 routput
      integer ibox,imolty,nmolty,numboxes
      double precision rarot,rtrot,ratran,rttran,datraa,dntraa
      double precision datrac,dntrac,darot,dnrot,onepi,rmrot
      double precision minvalue,tatraa,rmtraa,tatrac,rmtrac,tarot
      double precision minboxdim

c     --- get the runoutput style
      call twh_runoutput(GLB_GET,routput)

c     --- we need to perform some updates
c     --- adjust maximum translational displacement
      if ( routput .eq. 'full' .or. routput .eq. 'updates' ) then
         write(6,'(a55,a14)')
     &        'Updating maximum translational/rotational'
     &        ,' displacements'
      endif
c     --- retrieve constants
      call twh_constant_pi(GLB_GET,onepi)
      call twh_nmolty(GLB_GET,nmolty)
      call twh_numboxes(GLB_GET,numboxes)
      call twh_minboxdim(GLB_GET,minboxdim)

      do ibox = 1,numboxes
         do imolty = 1,nmolty
c           --- update rmtraa if any atom translation moves have
c           --- been attempted
            call twh_bntraa(GLB_GET,imolty,ibox,dntraa)
            call twh_batraa(GLB_GET,imolty,ibox,datraa)
            if ( dntraa .ne. 0.0d0 ) then
               call twh_tatraa(GLB_GET,tatraa)
c              --- compute possible new ratio
               ratran = datraa/dntraa
               call twh_rmtraa(GLB_GET,imolty,ibox,rmtraa)
               rttran = rmtraa* ratran / tatraa
               if ( rttran .gt. minboxdim ) then 
c                 --- maximum translational displacement
                  rmtraa = minboxdim
               elseif (rttran .lt. 1.0d-10 ) then  
c                 --- ratio must have been zero, so divide by 10
                  rmtraa = rmtraa/10.0d0
               else 
c                 --- accept new ratio
                  rmtraa = rttran
               endif
               minvalue = 0.001d0
               if ( rmtraa .lt. minvalue ) rmtraa = minvalue
               call twh_rmtraa(GLB_SET,imolty,ibox,rmtraa)
            endif

c           --- update rmtrac if any COM translation moves have been
c           --- attempted
            call twh_bntrac(GLB_GET,imolty,ibox,dntrac)
            call twh_batrac(GLB_GET,imolty,ibox,datrac)
            if ( dntrac .gt. 0.5d0 ) then
c              --- compute possible new ratio
               ratran = datrac/dntrac
               call twh_tatrac(GLB_GET,tatrac)
               call twh_rmtrac(GLB_GET,imolty,ibox,rmtrac)
               rttran = rmtrac* ratran / tatrac
               if ( rttran .gt. minboxdim ) then 
c                 --- maximum translational displacement
                  rmtrac = minboxdim
               elseif (rttran .lt. 1.0d-10 ) then  
c                 --- ratio must have been zero, so divide by 10
                  rmtrac = rmtrac/10.0d0
               else 
c                 --- accept new ratio
                  rmtrac = rttran
               endif
               minvalue = 0.01d0
               if ( rmtrac .lt. minvalue ) then
                  rmtrac = minvalue
               endif
               call twh_rmtrac(GLB_SET,imolty,ibox,rmtrac)
            endif
 
c           --- update rmrot if any rotation moves have been
c           --- attempted
            call twh_bnrot(GLB_GET,imolty,ibox,dnrot)
            call twh_barot(GLB_GET,imolty,ibox,darot)
            if ( dnrot .gt. 0.5d0 ) then
c              --- compute possible new ratio
               rarot = darot/dnrot
               call twh_tarot(GLB_GET,tarot)
               call twh_rmrot(GLB_GET,imolty,ibox,rmrot)
               rtrot = rmrot * rarot / tarot
               if ( rtrot .gt. onepi ) then 
c                 --- maximum translational displacement
                  rmrot = onepi
               elseif (rtrot .lt. 1.0d-10 ) then  
c                 --- ratio must have been zero, so divide by 10
                  rmrot = rmrot / 10.0d0
               else 
c                 --- accept new ratio
                  rmrot = rtrot
               endif
               call twh_rmrot(GLB_SET,imolty,ibox,rmrot)
            endif

c           --- write some ratio update information if we updated 
c           --- and of these values
            if ( ( dntraa + dntrac + dnrot .gt. 0.5d0) .and.
     &           ( routput .eq. 'full' 
     &           .or. routput .eq. 'updates' ) ) then
               write(6,'(1x,a4,i2,1x,a9,i3,2x,a8,3x,a8,3x,a16)') 
     &              'Box:',ibox,'Molecule:',imolty
     &              ,'Attempts','Accepted','New Displacement'
c              --- atom translation move info
               call twh_rmtraa(GLB_GET,imolty,ibox,rmtraa)
               if ( dntraa .gt. 0.5d0 ) then
                  write(6,'(1x,a14,5x,f10.0,1x,f10.0,1x,f10.6)')
     &                 'Translate Atom',dntraa,datraa,rmtraa
               endif
c              --- COM translation move info
               if ( dntrac .gt. 0.5d0 ) then
                  call twh_rmtrac(GLB_GET,imolty,ibox,rmtrac)
                  write(6,'(1x,a13,6x,f10.0,1x,f10.0,1x,f10.6)')
     &                 'Translate COM',dntrac,datrac,rmtrac
               endif
c              --- rotation move info
               if ( dnrot .gt. 0.5d0 ) then
                  call twh_rmrot(GLB_GET,imolty,ibox,rmrot)
                  write(6,'(1x,a14,5x,f10.0,1x,f10.0,1x,f10.6)')
     &                 'Rotate',dnrot,darot,rmrot
               endif
            endif

c           --- increment the total arrays
c           --- atom translation
            call twh_acstraa(GLB_INCR,imolty,ibox,datraa)
            call twh_acntraa(GLB_INCR,imolty,ibox,dntraa)
c           --- COM translation
            call twh_acstrac(GLB_INCR,imolty,ibox,datrac)
            call twh_acntrac(GLB_INCR,imolty,ibox,dntrac)
c           --- rotation
            call twh_acsrot(GLB_INCR,imolty,ibox,darot)
            call twh_acnrot(GLB_INCR,imolty,ibox,dnrot)
         enddo
      enddo
c     --- reinitialize the temporary counters
      call twh_batraa(GLB_INIT,nmolty,numboxes,0.0d0)
      call twh_bntraa(GLB_INIT,nmolty,numboxes,0.0d0)
      call twh_batrac(GLB_INIT,nmolty,numboxes,0.0d0)
      call twh_bntrac(GLB_INIT,nmolty,numboxes,0.0d0)
      call twh_barot(GLB_INIT,nmolty,numboxes,0.0d0)
      call twh_bnrot(GLB_INIT,nmolty,numboxes,0.0d0)

      return
      end
