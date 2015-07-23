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
      subroutine twh_accumulate_averages(lfailure,speden)
c     ******************************************************************
c     * accumulates the averages taken at the end of each move         *
c     *                                                                *
c     * split from mainloop 10-20-2005 by M.G. Martin                  *
c     * last modified 10-14-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_DERIVATIVE
#define FUNCTION_GET_NCMT
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision speden
      dimension speden(MAXBOX)
      logical lfailure
c     --- local variables
c     --- local scalars
      logical scp_eval_derivative
c     --- integer scalars
      integer ibox,iener,imolty,incount,ivap,iliq,nmolty,nchbox
      integer potentialstyle,numboxes,keycode
c     --- double precision scalars
      double precision value,uliquid,uvapor,totnumden,beta,masstotal
      double precision boxvolume
c     --- double precision arrays
      double precision energy
      dimension energy(2)

c     --- retrieve constants
      call twh_potentialstyle(GLB_GET,potentialstyle)
      call twh_nmolty(GLB_GET,nmolty)
      call twh_numboxes(GLB_GET,numboxes)

      do ibox = 1, numboxes
c        --- accumulate the single count averages
         incount = 1
c        --- accumulate molecule type independent values
         imolty = 0
c        --- accumulate energies and energies squared
c        --- total energy
         iener = 0
         if ( potentialstyle .eq. POT_INTERNAL ) then
            iener = iener + 1
            keycode = iener
            call twh_boxvclassic(GLB_GET,ibox,value)
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
c           --- intermolecular vdw
            iener = iener + 1
            keycode = iener
            call twh_boxvinter(GLB_GET,ibox,value)
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
c           --- bond bending
            iener = iener + 1
            keycode = iener
            call twh_boxvangle(GLB_GET,ibox,value)
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
c           --- torsion
            iener = iener + 1
            keycode = iener
            call twh_boxvtorsion(GLB_GET,ibox,value)
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
c           --- intramolecular
            iener = iener + 1
            keycode = iener
            call twh_boxvintra(GLB_GET,ibox,value)
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
c           --- external field
            iener = iener + 1
            keycode = iener
            call twh_boxvfield(GLB_GET,ibox,value)
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
c           --- vibrations
            iener = iener + 1
            keycode = iener
            call twh_boxvbond(GLB_GET,ibox,value)
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
c           --- coulombic
            iener = iener + 1
            keycode = iener
            call twh_boxvelect(GLB_GET,ibox,value)
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
c           --- vdw tail
            iener = iener + 1
            keycode = iener
            call twh_boxvtail(GLB_GET,ibox,value)
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
c           --- solvation
            iener = iener + 1
            keycode = iener
            call twh_boxvsolvation(GLB_GET,ibox,value)
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
c           --- accumulate dU/d[lambda] if scp_eval_derivative
            call twh_scp_eval_derivative(GLB_GET,scp_eval_derivative)
            if (scp_eval_derivative) then
c              --- dU/d[lambda_lj]
               keycode = AC_DU_DLAMBDA_LJ
               value = scp_get_derivative(ibox, LAMBDA_LJ)
               call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox
     &              ,imolty,value,incount)
               if ( lfailure ) return
c              --- dU/d[lambda_c]
               keycode = AC_DU_DLAMBDA_C
               value = scp_get_derivative(ibox, LAMBDA_C)
               call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox
     &              ,imolty,value,incount)
               if ( lfailure ) return
            endif
         endif
         if ( potentialstyle .eq. POT_EXTERNAL ) then
c           --- external
            iener = iener + 1
            keycode = iener
            call twh_boxvquantum(GLB_GET,ibox,value)
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
         endif
c        --- total energy squared
         keycode = AC_TOT_ENERGY_SQ
         if ( potentialstyle .eq. POT_INTERNAL ) then
            call twh_boxvclassic(GLB_GET,ibox,value)
            value = value**2
         elseif ( potentialstyle .eq. POT_EXTERNAL ) then
            call twh_boxvquantum(GLB_GET,ibox,value)
            value = value**2
         endif
         call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &        ,value,incount)
         if ( lfailure ) return
c        --- average volume (in cubic angstroms)
         keycode = AC_VOLUME
         call twh_boxvolume(GLB_GET,ibox,boxvolume)
         value =  boxvolume
         call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &        ,value,incount)
         if ( lfailure ) return
c        --- average volume squared (in cubic angstroms squared)
         keycode = AC_VOLUME_SQ
         value =  boxvolume*boxvolume
         call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &        ,value,incount)
         if ( lfailure ) return
c        --- accumulate molecule type dependent quantities
c        --- accumulate number of molecules
         do imolty = 1, nmolty
            value = dble(twh_get_ncmt(ibox,imolty))
            keycode = AC_NUM_MOLECULES
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
         enddo
c        --- accumulate number density (and compute spec den)
         speden(ibox) = 0.0d0
         totnumden = 0.0d0
         do imolty = 1, nmolty
            value = 1.0d3*dble(twh_get_ncmt(ibox,imolty))/boxvolume
            call twh_masstotal(GLB_GET,imolty,masstotal)
            speden(ibox) = speden(ibox) + masstotal*value/0.6022045d3
            totnumden = totnumden + value
            keycode = AC_NUMBER_DENSITY
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
         enddo
c        --- total number density (from the number densities)
         imolty = 0
         keycode = AC_NUMBER_DENSITY
         value = totnumden
         call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &        ,value,incount)
         if ( lfailure ) return
c        --- specific density (from the number densities)
         imolty = 0
         keycode = AC_SPECIFIC_DENSITY
         value = speden(ibox)
         call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &        ,value,incount)
         if ( lfailure ) return
c        --- accumulate mol fractions
         do imolty = 1, nmolty
            call twh_nchbox(GLB_GET,ibox,nchbox)
            if ( nchbox .ne. 0 ) then
               value = dble(twh_get_ncmt(ibox,imolty))/dble(nchbox)
            else
               value = 0.0d0
            endif
            keycode = AC_MOL_FRACTION
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
         enddo
      enddo
      if ( numboxes .eq. 2 .and. nmolty .eq. 1 ) then
c        --- heat of vaporization (pV/n=RT)
         keycode = AC_H_VAP_PV_NRT
         ibox = 0
         imolty = 1
         if ( twh_get_ncmt(1,imolty) * twh_get_ncmt(2,imolty)
     &        .ne. 0 ) then
c           --- determine the vapor and liquid boxes
            if ( speden(1) .gt. speden(2) ) then
               ivap = 2
               iliq = 1
            else
               ivap = 1
               iliq = 2
            endif
            if ( potentialstyle .eq. POT_INTERNAL ) then
               call twh_boxvclassic(GLB_GET,ivap,energy(ivap))
               call twh_boxvclassic(GLB_GET,iliq,energy(iliq))
            elseif ( potentialstyle .eq. POT_EXTERNAL ) then
               call twh_boxvquantum(GLB_GET,ivap,energy(ivap))
               call twh_boxvquantum(GLB_GET,iliq,energy(iliq))
            else
               write(6,*) 'ACCUMULATEAVERAGES: unknown potentialstyle'
               write(6,*) potentialstyle
               lfailure = .true.
               return
            endif
            uvapor = energy(ivap)*8.314d0/1000.0d0
            uliquid = energy(iliq)*8.314d0/1000.0d0
            call twh_beta(GLB_GET,beta)
            value = (uvapor/dble(twh_get_ncmt(ivap,imolty)))
     &           - (uliquid/dble(twh_get_ncmt(iliq,imolty)))
     &           + 8.314d0/(beta*1000.0d0)
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
         endif
      endif

      return
      end
