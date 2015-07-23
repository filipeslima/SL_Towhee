#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2011 Marcus G. Martin                       *
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
      subroutine twh_ffmgmstereo(lfailure,nimprotypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the mgmstereo    *
c     * forcefield                                                     *
c     *                                                                *
c     * originally written 02-11-2005 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer nimprotypes
c     --- local variables
      integer type
      character*(FFNAMELEN) ffname
      double precision iclocal
      dimension iclocal(MINIMPCOEFF:MAXIMPCOEFF)
      integer nislocal 

      ffname = 'mgmstereo'
      call twh_classical_potential(GLB_SET,'wild')
      call twh_classical_mixrule(GLB_SET,'wild')
      type = 0

c     --- references
c     This is just something that I made up in order to enforce 
c     stereochemistry using improper torsions.  Configurational-bias
c     is extremely good at changing stereochemistry and this prevents
c     such a thing from happening.
c     M.G. Martin 02-11-2005

c     ******************************************************************
c     * improper torsion                                               *
c     ******************************************************************

      type = type + 1
      call twh_impform(GLB_SET,type,1)
      call twh_impstyle(GLB_SET,type,3)
      iclocal(0) =  50.0d0
      iclocal(1) =  100.0d0
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 0
      call twh_nimpsame(GLB_SET,type,nislocal)
c     --- this is not assigned in the conventional search manner
c     --- instead there is a hard-wired check in Towhee

      nimprotypes = type
      call twh_ffcheck(lfailure,4,type)
      if ( lfailure ) return

      return
      end
