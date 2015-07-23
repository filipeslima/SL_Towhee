#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2011 Marcus G. Martin                       *
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
      function twh_minboxlength(ibox)
c     ******************************************************************
c     * finds the minimum dimension of a simulation box                *
c     *                                                                *
c     * originally written 12-14-2001 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_MINBOXLENGTH
#include "functions.h"

c     --- variables passed to/from the subroutine
      integer ibox

c     --- local variables
      integer jdim
      double precision boxlength,hbx

      if ( ibox .eq. 0 ) then
c        --- ideal gas box, infinite length, set to large value
         twh_minboxlength = 1.0d10
      else
         call twh_halfbox(GLB_GET,ibox,1,hbx)
         twh_minboxlength = 2.0d0 * hbx
         do jdim = 2,3
            call twh_halfbox(GLB_GET,ibox,jdim,hbx)
            boxlength = 2.0d0 * hbx
            if ( boxlength .lt. twh_minboxlength ) then
               twh_minboxlength = boxlength
            endif
         enddo
      endif

      return
      end
