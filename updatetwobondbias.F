#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2008 Marcus G. Martin                       *
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
      subroutine twh_updatetwobondbias(nnn)
c     ******************************************************************
c     * updates the parameters used in the two bond fixed endpoint     *
c     * biasing                                                        *
c     *                                                                *
c     * originally written 10-13-2005 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJBEN3
#define FUNCTION_GET_INBEN
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer nnn
c     --- local variables
      character*50 twobbias_style
      integer ibox,imolty,ichain,iunit,ibend,kunit,iparam,nchain,nmolty
      integer nunit,two_bond_bias_update_frequency
      integer two_bond_bias_compute_frequency
      double precision rxui,ryui,rzui,rxuk,ryuk,rzuk,rxuik,ryuik,rzuik
      double precision rik,riksq,samples,newmean,newdsq,newsdev
      double precision variance,twoboldfrac,dvalue

      call twh_two_bond_bias_compute_frequency(GLB_GET
     &     ,two_bond_bias_compute_frequency)
      if ( two_bond_bias_compute_frequency .eq. 0 ) then
c        --- updates are never performed for this variable
         return
      endif

      if ( mod(nnn,two_bond_bias_compute_frequency) .eq. 0 ) then
c        --- run through all of the molecules that exist in a real box
c        --- and compute the average distances needed for the two_bond
c        --- biasing
         call twh_two_bond_fe_bias_style(GLB_GET,twobbias_style)
         call twh_nchain(GLB_GET,nchain)
         do ichain = 1,nchain
            call twh_nboxi(GLB_GET,ichain,ibox)
            if ( ibox .ne. 0 ) then
               call twh_moltyp(GLB_GET,ichain,imolty)
               call twh_nunit(GLB_GET,imolty,nunit)
               do iunit = 1,nunit
                  call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &                 ,rxui,ryui,rzui)
                  do ibend = 1,twh_get_inben(imolty,iunit)
                     if ( twobbias_style .eq.
     &                    'self adapting gaussian using 1-3 distance'
     &                    ) then
                        kunit = twh_get_ijben3(imolty,iunit,ibend)
                        call twh_coordinates(GLB_GET,CRD_REAL,ichain
     &                       ,kunit,rxuk,ryuk,rzuk)
c                       --- compute i-k distance
                        call twh_mimage(ibox,rxui,ryui,rzui,rxuk,ryuk
     &                       ,rzuk,rxuik,ryuik,rzuik)
                        riksq = rxuik*rxuik + ryuik*ryuik 
     &                       + rzuik*rzuik
                        rik = dsqrt(riksq)
c                       --- increment count
                        call twh_two_bond_bias_info(GLB_INCR,imolty
     &                       ,iunit,ibend,3,1.0d0)
c                       --- increment sum of distances
                        call twh_two_bond_bias_info(GLB_INCR,imolty
     &                       ,iunit,ibend,4,rik)
c                       --- increment sum of distances squared
                        call twh_two_bond_bias_info(GLB_INCR,imolty
     &                       ,iunit,ibend,5,riksq)
                     else
                        write(6,*) 'UPDATETWOBONDBIAS: unknown '
     &                       ,'two_bond_fixed_endpoint_bias_style'
                        write(6,*) twobbias_style
                        stop
                     endif
                  enddo
               enddo
            endif
         enddo
      endif

      call twh_two_bond_bias_update_frequency(GLB_GET
     &     ,two_bond_bias_update_frequency)
      if ( mod(nnn,two_bond_bias_update_frequency) .eq. 0 ) then
c        --- mix the running results for the biasing with the current
c        --- values and then reset the running results to zero
         call twh_nmolty(GLB_GET,nmolty)
         call twh_two_bond_fe_bias_style(GLB_GET,twobbias_style)
         call twh_two_bond_bias_old_fraction(GLB_GET,twoboldfrac)
         do imolty = 1,nmolty
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
               do ibend = 1,twh_get_inben(imolty,iunit)
                  call twh_two_bond_bias_info(GLB_GET,imolty,iunit
     &                 ,ibend,3,samples)
                  if ( samples.gt. 0.5d0 ) then
                     if ( twobbias_style .eq.
     &                    'self adapting gaussian using 1-3 distance'
     &                    ) then
c                       --- compute the new mean and standard deviation
                        call twh_two_bond_bias_info(GLB_GET,imolty
     &                       ,iunit,ibend,4,dvalue)
                        newmean = dvalue / samples
                        call twh_two_bond_bias_info(GLB_GET,imolty
     &                       ,iunit,ibend,5,dvalue)
                        newdsq = dvalue / samples
                        variance = newdsq - newmean*newmean
                        if ( variance .gt. 0.0d0 ) then
                           newsdev = dsqrt(variance)
                        else
                           newsdev = 0.0d0
                        endif
c                       --- mix the old and the new
                        call twh_two_bond_bias_info(GLB_GET,imolty,iunit
     &                       ,ibend,1,dvalue)
                        dvalue = twoboldfrac * dvalue
     &                       + ( 1.0d0 - twoboldfrac ) * newmean
                        call twh_two_bond_bias_info(GLB_SET,imolty,iunit
     &                       ,ibend,1,dvalue)

                        call twh_two_bond_bias_info(GLB_GET,imolty,iunit
     &                       ,ibend,2,dvalue)
                        dvalue = twoboldfrac * dvalue
     &                       + ( 1.0d0 - twoboldfrac ) * newsdev
                        call twh_two_bond_bias_info(GLB_SET,imolty,iunit
     &                       ,ibend,2,dvalue)
                     else
                        write(6,*) 'UPDATETWOBONDBIAS: unknown '
     &                       ,'two_bond_fixed_endpoint_bias_style'
                        write(6,*) twobbias_style
                        stop
                     endif
                  endif
c                 --- reinitialize accumulators
                  do iparam = 3,MAXTWOBOND
                     call twh_two_bond_bias_info(GLB_SET,imolty,iunit
     &                    ,ibend,iparam,0.0d0)
                  enddo
               enddo
            enddo
         enddo
      endif

      return
      end
