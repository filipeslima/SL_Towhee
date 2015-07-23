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
      subroutine twh_initaverages(lfailure)
c     ******************************************************************
c     * initializes all of the quantities tracked by the averages      *
c     * subroutine                                                     *
c     *                                                                *
c     * split from mainloop 05-13-2005 by M.G. Martin                  *
c     * last modified 08-02-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
c     --- logical scalars
      logical scp_eval_derivative
c     --- integer scalars
      integer ibox,iener,imolty,incount,istress,nmolty,nmtype,numboxes
      integer nener,pressurestyle,keycode
c     --- double precision scalars
      double precision value

c     --- retrieve constants
      call twh_nmolty(GLB_GET,nmolty)
      call twh_numboxes(GLB_GET,numboxes)
      call twh_pressurestyle(GLB_GET,pressurestyle)

c     --- assign all of the averages
      if ( pressurestyle .eq. PRES_VIRIAL ) then
c        --- assign pressure and stress tensor
         imolty = 0
         do istress = 0,7
            keycode = AC_PRESSURE - istress
            do ibox = 1,numboxes
               call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty
     &              ,value,incount)
               if ( lfailure ) return
            enddo
         enddo
      elseif ( pressurestyle .eq. PRES_RADIAL ) then
c        --- assign isotropic pressure
         imolty = 0
         keycode = AC_PRESSURE
         do ibox = 1,numboxes
            call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
         enddo
      endif
c     --- assign dU/dV for each box
      imolty = 0
      keycode = AC_DU_DV
      do ibox = 1,numboxes
         call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty,value
     &        ,incount)
         if ( lfailure ) return
      enddo
c     --- assign thermodynamic pressure for each box
      imolty = 0
      keycode = AC_THERMO_PRESSURE
      do ibox = 1,numboxes
         call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty,value
     &        ,incount)
         if ( lfailure ) return
      enddo
c     --- assign energies
      imolty = 0
      call twh_nener(GLB_GET,nener)
      do iener = 1,nener
         keycode = iener
         do ibox = 1,numboxes
            call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
         enddo
      enddo
c     --- assign specific density
      keycode = AC_SPECIFIC_DENSITY
      imolty = 0
      do ibox = 1,numboxes
         call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty,value
     &        ,incount)
         if ( lfailure ) return
      enddo
c     --- assign total energy squared
      keycode = AC_TOT_ENERGY_SQ
      imolty = 0
      do ibox = 1,numboxes
         call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty,value
     &        ,incount)
         if ( lfailure ) return
      enddo
c     --- assign volume
      keycode = AC_VOLUME
      imolty = 0
      do ibox = 1,numboxes
         call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty,value
     &        ,incount)
         if ( lfailure ) return
      enddo
c     --- assign volume squared
      keycode = AC_VOLUME_SQ
      imolty = 0
      do ibox = 1,numboxes
         call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty,value
     &        ,incount)
         if ( lfailure ) return
      enddo
c     --- assign enthalpy
      keycode = AC_ENTHALPY
      imolty = 0
      do ibox = 1,numboxes
         call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty,value
     &        ,incount)
         if ( lfailure ) return
      enddo
c     --- assign enthalpy squared
      keycode = AC_ENTHALPY_SQ
      imolty = 0
      do ibox = 1,numboxes
         call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty,value
     &        ,incount)
         if ( lfailure ) return
      enddo
c     --- assign pressure*volume
      keycode = AC_PV
      imolty = 0
      do ibox = 1,numboxes
         call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty,value
     &        ,incount)
         if ( lfailure ) return
      enddo
c     --- assign Gibbs chemical potential
      keycode = AC_GIBBS_TOT_CHEMPOT
      do ibox = 1,numboxes
         do imolty = 1,nmolty
            call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
         enddo
      enddo
c     --- assign npt insertion chemical potential
      keycode = AC_NPT_INSERTION_CHEMPOT
      do ibox = 1,numboxes
         do imolty = 1,nmolty
            call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
         enddo
      enddo
c     --- assign nvt insertion chemical potential
      keycode = AC_NVT_INSERTION_CHEMPOT
      do ibox = 1,numboxes
         do imolty = 1,nmolty
            call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
         enddo
      enddo
c     --- assign isolated molecule chemical potential
      ibox = MAXBOX
      do imolty = 1,nmolty
         call twh_nmolectyp(GLB_GET,imolty,nmtype)
         if ( nmtype .eq. 0 ) then
            keycode = AC_ISOLATION_CHEMPOT
            call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
         endif
      enddo
c     --- assign number densities
      keycode = AC_NUMBER_DENSITY
      do ibox = 1,numboxes
         do imolty = 0,nmolty
            call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
         enddo
      enddo
c     --- assign mol fractions
      keycode = AC_MOL_FRACTION
      do ibox = 1,numboxes
         do imolty = 1,nmolty
            call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
         enddo
      enddo
c     --- assign radius of gyration
      keycode = AC_RADIUS_OF_GYRATION
      do ibox = 1,numboxes
         do imolty = 1,nmolty
            call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
         enddo
      enddo
c     --- assign number of molecules
      keycode = AC_NUM_MOLECULES
      do ibox = 1,numboxes
         do imolty = 1,nmolty
            call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
         enddo
      enddo
      ibox = 0
      imolty = 1
c     --- assign heat of vaporization (Direct)
      keycode = AC_H_VAP_DIRECT
      call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty,value
     &     ,incount)
      if ( lfailure ) return
c     --- assign heat of vaporization (vapor p)
      keycode = AC_H_VAP_VAPOR_P
      call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty,value
     &     ,incount)
      if ( lfailure ) return
c     --- assign heat of vaporization (pV/n=RT)
      keycode = AC_H_VAP_PV_NRT
      call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty,value
     &     ,incount)
      if ( lfailure ) return

c     --- assign dU/d[lambda_*] for each box of interest
      call twh_scp_eval_derivative(GLB_GET,scp_eval_derivative)
      if (scp_eval_derivative) then
         imolty = 0
         keycode = AC_DU_DLAMBDA_LJ
         do ibox = 1,numboxes
            call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
         enddo

         keycode = AC_DU_DLAMBDA_C
         do ibox = 1,numboxes
            call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox,imolty
     &           ,value,incount)
            if ( lfailure ) return
         enddo
      endif

      return
      end
