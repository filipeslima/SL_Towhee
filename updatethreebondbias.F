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
      subroutine twh_updatethreebondbias(nnn)
c     ******************************************************************
c     * updates the parameters used in the three bond fixed endpoint   *
c     * biasing                                                        *
c     *                                                                *
c     * originally written 10-14-2005 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_INTOR
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer nnn
c     --- local variables
      character*50 threebbias_style
c     --- integer scalars
      integer ibox,imolty,ichain,iunit,idihed,lunit,iparam,nchain
      integer nmolty,nunit,compute_freq,update_freq
c     --- double precision scalars
      double precision rxui,ryui,rzui,rxul,ryul,rzul,rxuil,ryuil,rzuil
      double precision ril,rilsq,samples,newmean,newdsq,newsdev
      double precision variance,dvalue,three_oldfrac

      call twh_three_bond_bias_compute_freq(GLB_GET,compute_freq)
      if ( compute_freq .eq. 0 ) then
c        --- no updates are ever performed for this variable
         return
      endif

      if ( mod(nnn,compute_freq) .eq. 0 ) then
c       --- run through all of the molecules that exist in a real box
c       --- and compute the average distances needed for the three_bond
c       --- biasing
         call twh_three_bond_fe_bias_style(GLB_GET,threebbias_style)
         call twh_nchain(GLB_GET,nchain)
         do ichain = 1,nchain
            call twh_nboxi(GLB_GET,ichain,ibox)
            if ( ibox .ne. 0 ) then
               call twh_moltyp(GLB_GET,ichain,imolty)
               call twh_nunit(GLB_GET,imolty,nunit)
               do iunit = 1,nunit
                  call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &                 ,rxui,ryui,rzui)
                  do idihed = 1,twh_get_intor(imolty,iunit)
                     if ( threebbias_style .eq.
     &                    'self adapting gaussian using 1-4 distance'
     &                    ) then
                        lunit = twh_get_ijtor4(imolty,iunit,idihed)
                        call twh_coordinates(GLB_GET,CRD_REAL,ichain
     &                       ,lunit,rxul,ryul,rzul)
c                       --- compute i-l distance
                        call twh_mimage(ibox,rxui,ryui,rzui,rxul,ryul
     &                       ,rzul,rxuil,ryuil,rzuil)
                        rilsq = rxuil*rxuil + ryuil*ryuil 
     &                       + rzuil*rzuil
                        ril = dsqrt(rilsq)
c                       --- increment count
                        call twh_three_bond_bias_info(GLB_INCR,imolty
     &                       ,iunit,idihed,3,1.0d0)
c                       --- increment sum of distances
                        call twh_three_bond_bias_info(GLB_INCR,imolty
     &                       ,iunit,idihed,4,ril)
c                       --- increment sum of distances squared
                        call twh_three_bond_bias_info(GLB_INCR,imolty
     &                       ,iunit,idihed,5,rilsq)
                     else
                        write(6,*) 'UPDATETHREEBONDBIAS: unknown '
     &                       ,'three_bond_fixed_endpoint_bias_style'
                        write(6,*) threebbias_style
                        stop
                     endif
                  enddo
               enddo
            endif
         enddo
      endif

      call twh_three_bond_bias_update_freq(GLB_GET,update_freq)
      if ( mod(nnn,update_freq) .eq. 0 ) then
c        --- mix the running results for the biasing with the current
c        --- values and then reset the running results to zero
         call twh_three_bond_fe_bias_style(GLB_GET,threebbias_style)
         call twh_three_bond_bias_old_fraction(GLB_GET,three_oldfrac)
         call twh_nmolty(GLB_GET,nmolty)
         do imolty = 1,nmolty
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
               do idihed = 1,twh_get_intor(imolty,iunit)
                  call twh_three_bond_bias_info(GLB_GET,imolty,iunit
     &                 ,idihed,3,samples)
                  if ( samples.gt. 0.5d0 ) then
                     if ( threebbias_style .eq.
     &                    'self adapting gaussian using 1-4 distance'
     &                    ) then
c                       --- compute the new mean and standard deviation
                        call twh_three_bond_bias_info(GLB_GET,imolty
     &                       ,iunit,idihed,4,dvalue)
                        newmean = dvalue / samples
                        call twh_three_bond_bias_info(GLB_GET,imolty
     &                       ,iunit,idihed,5,dvalue)
                        newdsq = dvalue / samples
                        variance = newdsq - newmean*newmean
                        if ( variance .gt. 0.0d0 ) then
                           newsdev = dsqrt(variance)
                        else
                           newsdev = 0.0d0
                        endif
c                       --- mix the old and the new
                        call twh_three_bond_bias_info(GLB_GET,imolty
     &                       ,iunit,idihed,1,dvalue)
                        dvalue = three_oldfrac * dvalue
     &                       + ( 1.0d0 - three_oldfrac ) * newmean
                        call twh_three_bond_bias_info(GLB_SET,imolty
     &                       ,iunit,idihed,1,dvalue)

                        call twh_three_bond_bias_info(GLB_GET,imolty
     &                       ,iunit,idihed,2,dvalue)
                        dvalue = three_oldfrac * dvalue
     &                       + ( 1.0d0 - three_oldfrac ) * newsdev
                        call twh_three_bond_bias_info(GLB_SET,imolty
     &                       ,iunit,idihed,2,dvalue)
                     else
                        write(6,*) 'UPDATETHREEBONDBIAS: unknown '
     &                       ,'three_bond_fixed_endpoint_bias_style'
                        write(6,*) threebbias_style
                        stop
                     endif
                  endif
c                 --- reinitialize accumulators
                  do iparam = 3,MAXTHREEBOND
                     call twh_three_bond_bias_info(GLB_SET,imolty,iunit
     &                    ,idihed,iparam,0.0d0)
                  enddo
               enddo
            enddo
         enddo
      endif

      return
      end
