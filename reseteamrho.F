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
      subroutine twh_reseteamrho(nbcount, boxes )
c     ******************************************************************
c     * resets all of the current eam_rho values to the temp_eam_rho   *
c     * values for the boxes specified                                 *
c     *                                                                *
c     * originally written 04-23-2004 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer nbcount,boxes
      dimension boxes(nbcount)
c     --- local variables
      integer icf,jchain,jlist,jmolty,cpnumber,nmolty,ncmt
      double precision density

c     --- only performed for EAM potential
      call twh_cp_number(GLB_GET,cpnumber)
      if ( cpnumber .ne. CP_EMBEDDED_ATOM_METHOD ) return

      call twh_nmolty(GLB_GET,nmolty)
      do icf = 1,nbcount
         do jmolty = 1,nmolty
            call twh_ncmt(GLB_GET,boxes(icf),jmolty,ncmt)
            do jlist = 1,ncmt
               call twh_chainlist(GLB_GET,jlist,boxes(icf),jmolty
     &              ,jchain)
               if ( boxes(icf) .eq. 0 ) then
c                 --- set density to zero
                  call twh_eam_rho(GLB_SET,CRD_REAL,jchain,0.0d0)
               else
c                 --- set density to temp_eam_rho
                  call twh_eam_rho(GLB_STORE,CRD_T_TO_R,jchain,density)
               endif
            enddo
         enddo
      enddo

      return
      end
