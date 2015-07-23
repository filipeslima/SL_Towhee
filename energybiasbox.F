#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2006-2012 Marcus G. Martin                       *
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
      subroutine twh_energybiasbox(lfailure,lgenerate,ibox,xxx,yyy,zzz
     &     ,bias)
c     ******************************************************************
c     * handles the non-uniform generation of an atom position using   *
c     * the energy-biasing map                                         *
c     *                                                                *
c     * split from swapmoves 11-02-2007 by M.G. Martin                 *
c     * last modified 09-18-2012 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_RANDOM
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure,lgenerate
      integer ibox
      double precision xxx,yyy,zzz,bias
c     --- local variables
      integer iaa,ibb,icc,cubea,cubeb,cubec,cubeindex
      double precision rndeb,cweight,asize,bsize,csize
      double precision totbias,aaa,bbb,ccc,dcubes

c     --- retrieve biasing values
      call twh_cubelet(GLB_GET,cubea,cubeb,cubec)
c     --- compute the cubelet sizes
      asize = 1.0d0/dble(cubea)
      bsize = 1.0d0/dble(cubeb)
      csize = 1.0d0/dble(cubec)
c     --- compute number of cubelets for biasing
      dcubes = dble(cubea*cubeb*cubec)

      if ( lgenerate ) then
c        --- select a position according to the energy bias
         totbias = 0.0d0
         rndeb = twh_random()
         do iaa = 1,cubea
            do ibb = 1,cubeb
               do icc = 1,cubec
c                 --- get the index
                  call twh_cubelet_index(iaa,ibb,icc,cubeindex)
                  call twh_cubeletweight(GLB_GET,cubeindex,cweight)
                  totbias = totbias + cweight
                  if ( rndeb .le. totbias ) then
c                    --- generate in arbitrary coordinates
                     aaa = asize*(dble(iaa-1) + twh_random())
                     bbb = bsize*(dble(ibb-1) + twh_random())
                     ccc = csize*(dble(icc-1) + twh_random())
c                    --- transform to cartesian coordinates
                     call twh_arbtocart(ibox,aaa,bbb,ccc,xxx,yyy,zzz)
c                    --- bias factor back to uniform (in this case 1/V)
c                    --- p_uniform = 1/Vbox
c                    --- p_gen = cweight* (1/Vcubelet)
c                    --- p_uniform/p_gen = 1 / (cweight * Vbox/Vcubelet)
c                    --- Vbox / Vcubelet = number of cubes
c                    --- bias = 1 / (cweight * cubes)
                     bias = 1.0d0/(cweight*dcubes)
                     return
                  endif
               enddo
            enddo
         enddo
      else
c        --- transform into arb coords
         call twh_carttoarb(ibox,xxx,yyy,zzz,aaa,bbb,ccc)
c        --- make sure this is in the central box
         call twh_putarbinbox(aaa,bbb,ccc)
c        --- evaluate the biasing correction for the old position
         iaa = idint(aaa/asize) + 1
         ibb = idint(bbb/bsize) + 1
         icc = idint(ccc/csize) + 1
         call twh_cubelet_index(iaa,ibb,icc,cubeindex)
         call twh_cubeletweight(GLB_GET,cubeindex,cweight)
c        --- normalize back to uniform
         bias = 1.0d0/(cweight*dcubes)
         return
      endif

c     --- if we got here then we did not find the correct biasing box
c     --- this is a fatal error
      write(6,*) 'ENERGYBOXBIAS: matching subcube not found'
      write(6,*) 'lgenerate:',lgenerate
      write(6,*) 'cubea,cubeb,cubec:',cubea,cubeb,cubec
      write(6,*) 'this is a fatal error'
      lfailure = .true.

      return
      end
