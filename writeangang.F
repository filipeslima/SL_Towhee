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
      subroutine twh_writeangang(index)
c     ******************************************************************
c     * writes all of the angle-angle parameters used in this          *
c     * simulation to the file index                                   *
c     * index: file index for output                                   *
c     *                                                                *
c     * originally written 02-22-2002 by M.G. Martin                   *
c     * last modified 08-09-2011 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_AACOEFF
#define FUNCTION_GET_AASTYLE
#define FUNCTION_GET_LAAHERE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer index
c     --- local variables
      logical lfound
      integer type
      double precision onepi

      call twh_constant_pi(GLB_GET,onepi)

      lfound = .false.
      do type = 1,TAAMAX
         if ( twh_get_laahere(type) )  then
            if ( .not. lfound ) then
c              --- just found the first angle-angle type
               write(index,'(a17)') 'Angle-Angle Types'
               lfound = .true.
            endif
c           --- output information about this interaction
            if ( twh_get_aastyle(type) .eq. 1 ) then
c              --- Compass angle-angle
               write(index,64) type,twh_get_aacoeff(type,0)
            elseif ( twh_get_aastyle(type) .eq. 2 ) then
c              --- Compass angle-angle with explicit angle terms
               write(index,65) type,twh_get_aacoeff(type,0)
     &              ,twh_get_aacoeff(type,1)*180.0d0/onepi
     &              ,twh_get_aacoeff(type,2)*180.0d0/onepi
            else
               write(6,*) 'WRITEANGANG: nondefined angle-angle type '
     &              ,twh_get_aastyle(type)
               stop
            endif
         endif
      enddo
 64   format('Type: ',i3,' Style: Compass  k:',2x,f10.1)
 65   format('Type: ',i3,' Style: Compass  k:',2x,f10.1
     &     ,1x,' Angle(1) ',f10.3,' Angle(2) ',f10.3)
      
      return
      end
    
