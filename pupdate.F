#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2011 Marcus G. Martin                       *
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
      subroutine twh_pupdate(lfailure)
c     ******************************************************************
c     * this subroutine performs all of the end of block operations    *
c     *                                                                *
c     * orignally taken out of old monola in 1998 by M.G. Martin       *
c     * last modified 08-08-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
      character*20 routput
      character*30 keyname
c     --- integer scalars
      integer imolty,ibox,imove,nblock,iener,nmolty,numboxes,nener
      integer pressurestyle,keycode
c     --- double precision scalars
      double precision baver,invalue,beta

c     --- retrieve constants
      call twh_beta(GLB_GET,beta)
      call twh_numboxes(GLB_GET,numboxes)
      call twh_pressurestyle(GLB_GET,pressurestyle)

c     --- calculation of block averages
c     --- end this block
      call twh_averages(lfailure,AVG_END_OF_BLOCK,keycode,ibox,imolty
     &     ,invalue,nblock)

      if ( lfailure ) return
      call twh_runoutput(GLB_GET,routput)
      if ( routput .eq. 'full' .or. routput .eq. 'blocks' ) then
c        --- retrieve number of molecule types
         call twh_nmolty(GLB_GET,nmolty)
c        --- retrieve and output information about this block
         write(6,'(a29,1x,i5,1x)') 'Block Averages (BA) for block'
     &        ,nblock
         do ibox = 1,numboxes
c           --- volume
            keycode = AC_VOLUME
            imolty = 0
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL,keycode
     &           ,ibox,imolty,baver ,nblock)
            if ( lfailure ) return
            write(6,'(a7,i2,1x,a12,1x,e14.8)') 'BA Box:',ibox 
     &           ,'Volume [A^3]' ,baver
c           --- volume squared
            keycode = AC_VOLUME_SQ
            imolty = 0
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL,keycode
     &           ,ibox,imolty,baver ,nblock)
            if ( lfailure ) return
            write(6,'(a7,i2,1x,a9,1x,e14.8)') 'BA Box:',ibox 
     &           ,'V^2 [A^6]' ,baver
c           --- specific density
            keycode = AC_SPECIFIC_DENSITY
            imolty = 0
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL,keycode
     &           ,ibox,imolty,baver ,nblock)
            if ( lfailure ) return
            write(6,'(a7,i2,1x,a23,1x,e14.8)') 'BA Box:',ibox 
     &           ,'Specific density [g/ml]' ,baver
            if ( pressurestyle .eq. PRES_VIRIAL ) then
c              --- total virial pressure
               keycode = AC_PRESSURE
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL,keycode
     &              ,ibox,imolty,baver,nblock)
               if ( lfailure ) return
               call twh_pcomp(GLB_GET,0,keyname)
               write(6,'(a7,i2,1x,a6,1x,a13,1x,a5,1x,e14.8)') 
     &              'BA Box:',ibox,'Virial',keyname,'[kPa]',baver
            endif
c           --- energies
            call twh_nener(GLB_GET,nener)
            do iener = 1,nener
               keycode = iener
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL,keycode
     &              ,ibox,imolty,baver,nblock)
               if ( lfailure ) return
               call twh_vname(GLB_GET,iener,keyname)
               write(6,'(a7,i2,1x,a15,1x,e14.8)') 'BA Box:',ibox
     &              ,keyname,baver
            enddo

c           --- gibbs chemical potential
            keycode = AC_GIBBS_TOT_CHEMPOT
            do imolty = 1, nmolty
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL,keycode
     &              ,ibox,imolty,baver,nblock)
               if ( lfailure ) return
c              --- convert into chemical potential
               if ( baver .ne. 0.0d0 ) then
                  baver = ((-1.0d0)/beta)*dlog(baver)
               endif
               write(6,'(a7,i2,1x,a24,1x,i3,1x,e14.8)')
     &              'BA Box:',ibox,'u (Gibbs Total) [K] Type'
     &              ,imolty,baver
            enddo
c           --- number density
            keycode = AC_NUMBER_DENSITY
            do imolty = 1, nmolty
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL,keycode
     &              ,ibox,imolty,baver,nblock)
               if ( lfailure ) return
               write(6,'(a7,i2,1x,a26,1x,i3,1x,e14.8)')
     &              'BA Box:',ibox,'Number density [nm-3] Type'
     &              ,imolty,baver
            enddo
c           --- mol fraction
            keycode = AC_MOL_FRACTION
            do imolty = 1, nmolty
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL,keycode
     &              ,ibox,imolty,baver,nblock)
               if ( lfailure ) return
               write(6,'(a7,i2,1x,a17,1x,i3,1x,e14.8)')
     &              'BA Box:',ibox,'Mol Fraction Type',imolty
     &              ,baver
            enddo
            if ( pressurestyle .eq. PRES_VIRIAL ) then
c              --- virial stress tensor components
               imolty = 0
               do imove = 1,7
                  keycode = AC_PRESSURE - imove
                  call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL
     &                 ,keycode,ibox,imolty,baver,nblock)
                  if ( lfailure ) return
                  call twh_pcomp(GLB_GET,imove,keyname)
                  write(6
     &                 ,'(a7,i2,1x,a13,1x,a6,1x,a13,1x,a5,1x,e14.8)')
     &                 'BA Box:',ibox,'Stress Tensor','Virial'
     &                 ,keyname,'[kPa]',baver
               enddo
            endif
c           --- radius of gyration
            keycode = AC_RADIUS_OF_GYRATION
            do imolty = 1,nmolty
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL,keycode
     &              ,ibox,imolty,baver,nblock)
               if ( lfailure ) return
               write(6,'(a7,i2,1x,a18,1x,a5,i3,1x,f14.5)')
     &              'BA Box:',ibox,'Radius of Gyration','Type:'
     &              ,imolty,baver
            enddo
         enddo
      endif

      return
      end
