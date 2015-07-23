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
      subroutine twh_fflast1993(lfailure,ntypes)
c     ******************************************************************
c     * sets up parameters for dinitrogen gas treates as a single bead *
c     * from a study by Lastoskie, Gubbins, and Quirke                 *
c     *                                                                *
c     * originally written 11-22-2004 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes
c     --- local variables
      integer ii,type
      character*(FFNAMELEN) ffname
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal

      ffname = 'Last1993'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'LB or Geometric')
      type = 0

c     ******************************************************************
c     * Parameter taken from                                           *
c     C. Lastoskie; K.E. Gubbins; N. Qurike; "Pore Size Heterogeneity
c     and the Carbon Slit Pore: A Density Functional Theory Model",
c     Langmuir, 9, 2693-2702 (1993)
c     ******************************************************************

c     --- dinitrogen gas (N2)
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.572d0
      nbclocal(2) = 93.98d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 28.014d0
      elelocal = ' N'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='N2'
      cbnlocal = 'N2'
      canlocal = 'N2'
      ctnlocal = 'N2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end
