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
      function twh_vonefive(lfailure,rijsq,ioftype,cutstyle,lovrlap )
c     ******************************************************************
c     * Computes the nonbond two-body energy for a given one-five type *
c     * and distance                                                   *
c     * cutstyle has several different values which do different things*
c     * -1: setting up the shifted potential, only done at begining    *
c     *  0: use the normal vdw nonbonded cutoff of rcutsq              *
c     *  1: use the inner cuttoff, unless rcut < rcutin                *
c     *  3: compute 1-4 interaction using rcutin and special 1-4 params*
c     *                                                                *
c     * originally written as vnonbond 2-18-2000 by M.G. Martin        *
c     * last update 08-09-2011 by M.G. Martin                          *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_OFCOEFF
#define FUNCTION_GET_OFSTYLE
#define FUNCTION_VONEFIVE
#include "functions.h"
c     --- debug toggle
#define DEBUG_VONEFIVE 0

c     --- variables passed to/from the function
      logical lfailure,lovrlap
      integer ioftype,cutstyle
      double precision rijsq

c     --- local variables
      logical ldummy
      double precision rmaxsq,srsix,rtwelve,rsix,rminsq

#if DEBUG_VONEFIVE
      write(6,*) 'starting VONEFIVE'
#endif
c     --- retrieve constants
      call twh_rminsq(GLB_GET,rminsq)

c     --- initialize vonefive
      twh_vonefive = 0.0d0

c     --- determine the cutoff
      call twh_getnbinfo(lfailure,cutstyle,rijsq,rmaxsq,ldummy)
      if ( lfailure ) return

c     -- compute the nonbond energy according to the classical potential
      if ( twh_get_ofstyle(ioftype) .eq. 1 ) then
c        --- lennard jones 12-6
c        --- using pretabulated sigma^6 and 4*epsilon
         if ( rijsq .lt. rminsq ) then
            lovrlap = .true.
         elseif ( rijsq .lt. rmaxsq ) then
            srsix = twh_get_ofcoeff(ioftype,3) / (rijsq**3)
            twh_vonefive = srsix*(srsix-1.0d0)
     &           *twh_get_ofcoeff(ioftype,4)
     &           - twh_get_ofcoeff(ioftype,0)
         endif
      elseif ( twh_get_ofstyle(ioftype) .eq. 2 ) then
c        --- a/r^12
         if ( rijsq .lt. rminsq ) then
            lovrlap = .true.
         elseif ( rijsq .lt. rmaxsq ) then
            rsix = rijsq*rijsq*rijsq
            rtwelve = rsix*rsix
            twh_vonefive = twh_get_ofcoeff(ioftype,1)/rtwelve
         endif
      else
         write(6,*) 'VONEFIVE: non defined ofstyle'
         lfailure = .true.
         return
      endif

#if DEBUG_VONEFIVE
      write(6,*) 'twh_vonefive:',twh_vonefive
      write(6,*) 'finished vonefive'
#endif
      return
      end

      function twh_onefivetype(imolty,iunit,junit)
c     ******************************************************************
c     * returns the type number of the one-five interaction between    *
c     * atoms iunit and junit on molecule type imolty                  *
c     * if there is no one-five interaction then returns zero          *
c     *                                                                *
c     * originally written 12-15-2004 by M.G. Martin                   *
c     * last modified 04-22-2009 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJOF5
#define FUNCTION_GET_INOF
#define FUNCTION_ONEFIVETYPE
#include "functions.h"
c     --- variables passed to/from the function
      integer imolty,iunit,junit
c     --- local variables
      integer iof

c     --- initialize onefivetype
      twh_onefivetype = 0
c     --- see if there is a special one-five interaction between these
c     --- atoms on this molecule type
      do iof = 1,twh_get_inof(imolty,iunit)
         if ( twh_get_ijof5(imolty,iunit,iof) .eq. junit ) then
c           --- match found
            call twh_itof(GLB_GET,imolty,iunit,iof,twh_onefivetype)
            return
         endif
      enddo
      
      return
      end

