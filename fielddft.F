#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2008 Marcus G. Martin                       *
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
      subroutine twh_fielddft(lfailure,iflag,activebox,nchainmove
     &     ,chaininfo,vdft)
c     ******************************************************************
c     * computes the solvation free energy for one simulation box      *
c     * this is done by calling the Tramonto density functional code   *
c     * that code is not yet generally available, so this won't work   *
c     * for anyone other than the code developers at this time         *
c     * if you have the tramonto code and want to use this then set    *
c     *                                                                *
c     * originally written 06-18-2002 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer iflag,activebox,nchainmove
      integer chaininfo
      dimension chaininfo(nchainmove,3)
      double precision vdft
#ifdef USETRAMONTO
      double precision tramonto,temperature
#endif
      vdft = 0.0d0

      call twh_writetramonto(lfailure,iflag,activebox,nchainmove
     &    ,chaininfo)
      if ( lfailure ) return
c     --- call tramonto, tramonto returns a unitless energy
c     --- need to multiply by the temperature to get Kelvin energies
#ifdef USETRAMONTO
      call twh_temperature(GLB_GET,temperature)
      vdft = tramonto()*temperature
#endif
      return
      end

