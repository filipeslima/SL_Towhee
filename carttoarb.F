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
      subroutine twh_carttoarb (ibox,rxuij,ryuij,rzuij,arba,arbb,arbc)
c     ******************************************************************
c     * transforms from the cartesian coordinate system to the         *
c     * arbitrary coordinate system (on 0,1)                           *
c     *                                                                *
c     * originally written 12-14-2001 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer ibox
      double precision arba,arbb,arbc,rxuij,ryuij,rzuij
c     --- local variables
      double precision hinv1,hinv2,hinv3

c      write(6,*) 'c to a r*uij',rxuij,ryuij,rzuij
      call twh_hinverse(GLB_GET,ibox,1,1,hinv1)
      call twh_hinverse(GLB_GET,ibox,2,1,hinv2)
      call twh_hinverse(GLB_GET,ibox,3,1,hinv3)
      arba = rxuij*hinv1 + ryuij*hinv2 + rzuij*hinv3
      call twh_hinverse(GLB_GET,ibox,1,2,hinv1)
      call twh_hinverse(GLB_GET,ibox,2,2,hinv2)
      call twh_hinverse(GLB_GET,ibox,3,2,hinv3)
      arbb = rxuij*hinv1 + ryuij*hinv2 + rzuij*hinv3
      call twh_hinverse(GLB_GET,ibox,1,3,hinv1)
      call twh_hinverse(GLB_GET,ibox,2,3,hinv2)
      call twh_hinverse(GLB_GET,ibox,3,3,hinv3)
      arbc = rxuij*hinv1 + ryuij*hinv2 + rzuij*hinv3

c      write(6,*) 'c to a arb* ',arba,arbb,arbc

      return
      end
