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
      function twh_getnbtype ( ntii,ntjj )
c     ******************************************************************
c     * turns the integer codes for two single atom types into the     *
c     * integer code for the mixed type.  This is used to reference    *
c     * the nbcoeffs array                                             *
c     *                                                                *
c     * originally written 11-15-2002 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ntii,ntjj
c     --- local variables
      integer low,high,istep,classical_pot_num

      if ( ntii .lt. ntjj ) then
         low = ntii
         high = ntjj
      else
         low = ntjj
         high = ntii
      endif

      istep = low-1
      
      call twh_cp_number(GLB_GET,classical_pot_num)
      if (classical_pot_num .lt. 100) then
c        --- all non-tabular potentials are < 100
         istep = istep*NNTYPE - (istep*(istep-1))/2
      elseif ( classical_pot_num .eq. CP_EMBEDDED_ATOM_METHOD .or.
     &        classical_pot_num .eq. CP_EAM_PAIR_ONLY .or.
     &        classical_pot_num .eq. CP_TABULATED_PAIR .or.
     &        classical_pot_num .eq. CP_MULTIWELL .or.
     &        classical_pot_num .eq. CP_REPULSIVE_MULTIWELL) then
c        --- tabular potential, use different index
         istep = istep*MAXTABTYPE - (istep*(istep-1))/2
      else
         write(6,*) 'GETNBTYPE: unknown classical_pot_num'
         write(6,*) classical_pot_num
         stop
      endif
      
      twh_getnbtype = istep + high - low + 1

      return
      end
