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
      subroutine twh_resetcbmc(arraydim,wrosen,vrostot,vrosvib,vrosang
     &     ,vrostor,vrosext,vrosintra,vrosinter,vroselect,vrosewald)
c     ******************************************************************
c     * initializes all of the variables used during a cbmc move       *
c     *                                                                *
c     * originally written 05-08-2002 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer arraydim
      double precision wrosen
      dimension wrosen(arraydim)
      double precision vrostot
      dimension vrostot(arraydim)
      double precision vrosvib
      dimension vrosvib(arraydim)
      double precision vrosang
      dimension vrosang(arraydim)
      double precision vrostor
      dimension vrostor(arraydim)
      double precision vrosext
      dimension vrosext(arraydim)
      double precision vrosintra
      dimension vrosintra(arraydim)
      double precision vrosinter
      dimension vrosinter(arraydim)
      double precision vroselect
      dimension vroselect(arraydim)
      double precision vrosewald
      dimension vrosewald(arraydim)
c     --- local variables
      integer ioldnew

      do ioldnew = 1,arraydim
c        --- set the initial weight to unity
         wrosen(ioldnew) = 1.0d0
c        --- set the energies of the trials to zero
         vrostot(ioldnew) = 0.0d0
         vrosvib(ioldnew) = 0.0d0
         vrosang(ioldnew) = 0.0d0
         vrostor(ioldnew) = 0.0d0
         vrosext(ioldnew) = 0.0d0
         vrosintra(ioldnew) = 0.0d0
         vrosinter(ioldnew) = 0.0d0
         vroselect(ioldnew) = 0.0d0
         vrosewald(ioldnew) = 0.0d0
      enddo
      call twh_hyfrom(GLB_SET,'null')
      call twh_hyprev(GLB_SET,'null')

      return
      end
