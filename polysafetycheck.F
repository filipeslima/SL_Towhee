#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2003-2011 Marcus G. Martin                       *
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
      subroutine twh_polysafetycheck(lfailure,iatom,total)
c     ******************************************************************
c     * makes sure we are not going to go over any NUMAX bound arrays  *
c     *                                                                *
c     * originally written 07-29-2003 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     ---- variables passed to/from the subroutine
      logical lfailure
      integer iatom,total

      if ( iatom + total .gt. NUMAX ) then
         write(6,*) 'POLYSAFETYCHECK: total atoms requested'
     &        ,' exceeds NUMAX'
         write(6,*) 'atoms requested:',iatom
         write(6,*) 'current atoms allocated:',total
         write(6,*) 'NUMAX:',NUMAX
         lfailure = .true.
         return
      endif
      
      return
      end
