#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2013 Marcus G. Martin                       *
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
      subroutine twh_writeaverages(lfailure)
c     ******************************************************************
c     * outputs all of the final running and block averages            *
c     *                                                                *
c     * split from mainloop on 05-12-2005 by M.G. Martin               *
c     * last modified 07-24-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
c     --- logicals
      logical lhenry_good,scp_eval_derivative
c     --- character strings
      character*30 filename,pressurestr,keyname
c     --- integer scalars
      integer ibox,imolty,incount,nblock,iblock,ivap,iliq,istress,iener
      integer jbox,outcount,pcount,tpcount,jmolty,zz,istart,iend,stemp
      integer nmolty,nunit,nmtype,numboxes,ensemble,nener,pressurestyle
      integer keycode
c     --- integer arrays
      integer success
      dimension success(NTMAX)
c     --- double precision scalars
      double precision value,ostwald,stdost,gconst,c_v,c_p,stdev
      double precision henry_convert,henry_coeff,chempotential
      double precision vdum,voltemp,temperature,beta,debroglie
c     --- double precision arrays
      double precision molvol,speden,acvol
      dimension molvol(MAXBOX),speden(MAXBOX),acvol(MAXBOX)
      double precision baver
      dimension baver(8+NTMAX)
      double precision tpressure
      dimension tpressure(MAXBOX)
      double precision idealp
      dimension idealp(MAXBOX)
      double precision dudv
      dimension dudv(MAXBOX)
      double precision gibbschem,nvtchem,nptchem,isochem,residchem
      dimension isochem(NTMAX)
      dimension gibbschem(MAXBOX,NTMAX),residchem(MAXBOX,NTMAX)
      dimension nvtchem(MAXBOX,NTMAX),nptchem(MAXBOX,NTMAX)
      double precision compression
      dimension compression(MAXBOX)
      double precision acenthalpy,acenthalpysq
      dimension acenthalpy(MAXBOX),acenthalpysq(MAXBOX)
      double precision numdenchem
      dimension numdenchem(MAXBOX,NTMAX)
      double precision ntotal,enthalpy,uterm,pvterm,hpermolec
      dimension ntotal(MAXBOX),hpermolec(MAXBOX)
      dimension uterm(MAXBOX)
      dimension enthalpy(MAXBOX),pvterm(MAXBOX)
      double precision freeenergy,entropy
      dimension freeenergy(MAXBOX),entropy(MAXBOX)
      double precision acrg
      dimension acrg(MAXBOX,NTMAX)
      double precision acv,acvboxsq
      dimension acv(MAXENER,MAXBOX)
      dimension acvboxsq(MAXBOX)
      double precision acstress
      dimension acstress(MAXBOX,0:7)
      double precision acdens,molfra
      dimension acdens(MAXBOX,NTMAX),molfra(MAXBOX,NTMAX)
      double precision acnbox
      dimension acnbox(MAXBOX,NTMAX)
      double precision totnumden
      dimension totnumden(MAXBOX)
      double precision dudl_lj, dudl_c
      dimension dudl_lj(MAXBOX), dudl_c(MAXBOX)
      double precision acvolsq
      dimension acvolsq(MAXBOX)

c     --- retrieve constants
      call twh_nmolty(GLB_GET,nmolty)
      call twh_temperature(GLB_GET,temperature)
      call twh_beta(GLB_GET,beta)
      call twh_scp_eval_derivative(GLB_GET,scp_eval_derivative)
      call twh_numboxes(GLB_GET,numboxes)
      call twh_ensemble(GLB_GET,ensemble)
      call twh_nener(GLB_GET,nener)
      call twh_pressurestyle(GLB_GET,pressurestyle)

c     --- set the henry's law conversion k_b * nm-3 to m-3 * Pa to MPa
      henry_convert = (1.38066d-23)*(1.0d27)*(1.0d-6)

c     --- check the vname(1) average to see if we ever did any as this
c     --- is saved every step and will tell us if we have done any moves
      keycode = 1
      ibox = 1
      imolty = 0
      call twh_averages(lfailure,AVG_EXTRACT_COUNT,keycode,ibox,imolty
     &     ,value,incount)
      if ( lfailure ) return
c     --- if incount is zero then there are no averages to report
      if ( incount .eq. 0 ) return

c     --- set pressurestr for output
      if ( pressurestyle .eq. PRES_VIRIAL) then
         pressurestr = 'Virial'
      elseif ( pressurestyle .eq. PRES_RADIAL) then
         pressurestr = 'Radial'
      elseif ( pressurestyle .eq. PRES_THERMO ) then
         pressurestr = 'Thermo'
      else
         stop 'Unknown pressure style in writeaverages.F'
      endif


c     --- on the fly pressure count
      if ( pressurestyle .eq. PRES_VIRIAL .or.
     &     pressurestyle .eq. PRES_RADIAL ) then
         keycode = AC_PRESSURE
         ibox = 1
         imolty = 0
         call twh_averages(lfailure,AVG_EXTRACT_COUNT,keycode,ibox
     &        ,imolty,value,pcount)
         if ( lfailure ) return
      else
         pcount = 0
      endif
c     --- set istart and iend for the stress/pressure depending upon
c     --- the on the fly pressure style
      if ( pressurestyle .eq. PRES_VIRIAL ) then
         istart = 0
         iend = 7
      elseif ( pressurestyle .eq. PRES_RADIAL ) then
         istart = 0
         iend = 0
      else
         istart = 0
         iend = 0
      endif
c     --- calculate and write out the running averages if we have 
c     --- performed any steps

c     --- do the final computations on the running averages
      do ibox = 1, numboxes
c        --- volume
         keycode = AC_VOLUME
         imolty = 0
         call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &        ,imolty,value,incount)
         if ( lfailure ) return
         acvol(ibox) = value
c        --- volume
         keycode = AC_VOLUME_SQ
         imolty = 0
         call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &        ,imolty,value,incount)
         if ( lfailure ) return
         acvolsq(ibox) = value
         if ( pcount .ne. 0 ) then
c           --- virial stress tensor
            do istress = istart,iend
               keycode = AC_PRESSURE - istress
               call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode
     &              ,ibox,imolty,value,incount)
               if ( lfailure ) return
               acstress(ibox,istress) = value
            enddo
c           --- enthalpy
            keycode = AC_ENTHALPY
            call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &           ,imolty,value,incount)
            if ( lfailure ) return
            acenthalpy(ibox) = value
            keycode = AC_ENTHALPY_SQ
            call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &           ,imolty,value,incount)
            if ( lfailure ) return
            acenthalpysq(ibox) = value
         endif
c        --- total energy squared
         keycode = AC_TOT_ENERGY_SQ
         call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &        ,imolty,value,incount)
         if ( lfailure ) return
         acvboxsq(ibox) = value
c        --- energies
         do iener = 1,nener
            keycode = iener
            call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &           ,imolty,value,incount)
            if ( lfailure ) return
            acv(iener,ibox) = value
         enddo

c        --- averages that are moltype dependent
         do imolty = 1, nmolty
c           --- number of molecules
            keycode = AC_NUM_MOLECULES
            call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &           ,imolty,value,incount)
            if ( lfailure ) return
            acnbox(ibox,imolty) = value
c           --- molfraction
            keycode = AC_MOL_FRACTION
            call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &           ,imolty,value,incount)
            if ( lfailure ) return
            molfra(ibox,imolty) = value
         enddo

         do imolty = 1, nmolty
c           --- number density
            keycode = AC_NUMBER_DENSITY
            call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &           ,imolty,value,incount)
            if ( lfailure ) return
            acdens(ibox,imolty) = value
         enddo
c        --- total number density
         keycode = AC_NUMBER_DENSITY
         imolty = 0
         call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &        ,imolty,value,incount)
         if ( lfailure ) return
         totnumden(ibox) = value
c        --- molar volume
c        --- sum over all types of molecules
         if ( totnumden(ibox) .gt. 1.0d-10 ) then
            molvol(ibox) = 602.2045d0 / totnumden(ibox)
         else
            molvol(ibox) = 0.0d0
         endif
c        --- specific density
         keycode = AC_SPECIFIC_DENSITY
         imolty = 0
         call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &        ,imolty,value,incount)
         if ( lfailure ) return
         speden(ibox) = value

c        --- radius of gyration
         keycode = AC_RADIUS_OF_GYRATION
         do imolty = 1,nmolty
            call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &           ,imolty,value,incount)
            if ( lfailure ) return
            acrg(ibox,imolty) = value
         enddo

c        --- dU/d[lambda] if we are calculating it
         if (scp_eval_derivative) then
            imolty = 0
            keycode = AC_DU_DLAMBDA_LJ
            call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &           ,imolty,value,incount)
            if ( lfailure ) return
            dudl_lj(ibox) = value

            keycode = AC_DU_DLAMBDA_C
            call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &           ,imolty,value,incount)
            if ( lfailure ) return
            dudl_c(ibox) = value
         endif
   

c        --- chemical potential
         do imolty = 1, nmolty
c           --- get the averages for the chemical potential
            keycode = AC_GIBBS_TOT_CHEMPOT
            call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &           ,imolty,value,incount)
            if ( lfailure ) return
            gibbschem(ibox,imolty) = value
            keycode = AC_NPT_INSERTION_CHEMPOT
            call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &           ,imolty,value,incount)
            if ( lfailure ) return
            nptchem(ibox,imolty) = value
            keycode = AC_NVT_INSERTION_CHEMPOT
            call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &           ,imolty,value,incount)
            if ( lfailure ) return
            nvtchem(ibox,imolty) = value
c           --- initialize insertion and numden chem potentials
            numdenchem(ibox,imolty) = 0.0d0
c           --- compute the number of successful growths into ibox 
c           --- from a box other than ibox as we cannot compute chem 
c           --- pot for ibox to ibox
            success(imolty) = 0
            do jbox = 0,numboxes
               if ( ibox .ne. jbox ) then
                  call twh_bscbswap(GLB_GET,imolty,ibox,jbox,stemp)
                  success(imolty) = success(imolty) + stemp
               endif
            enddo
            if ( success(imolty) .ne. 0 ) then
c              --- transform the gibbs chemical potential
               if ( gibbschem(ibox,imolty) .gt. 0.0d0 ) then
                  gibbschem(ibox,imolty) = -temperature* 
     &                 dlog(gibbschem(ibox,imolty))
               else
c                 --- infinitely high chemical potential, set to zero
                  gibbschem(ibox,imolty) = 0.0d0
               endif
c              --- transform the npt insertion chempot and 
c              --- add in the average volume term
c              --- from Eqn 2.69b from allen and tildesley
               if ( nptchem(ibox,imolty) .gt. 0.0d0 ) then
                  nptchem(ibox,imolty)  = -temperature* 
     &                 dlog(nptchem(ibox,imolty)/acvol(ibox))
               else
                  nptchem(ibox,imolty) = 0.0d0
               endif
c              --- perform the final rescaling of the nvt chemical
c              --- potential including the final average densities
c              --- store the insertion chemical potential
               if ( nvtchem(ibox,imolty) .gt. 0.0d0 ) then
                  nvtchem(ibox,imolty) = -temperature* 
     &                 dlog( nvtchem(ibox,imolty) )
               else
                  nvtchem(ibox,imolty) = 0.0d0
               endif
               if ( acnbox(ibox,imolty) .ne. 0.0d0 ) then
                  vdum = acvol(ibox)/acnbox(ibox,imolty)
c                 --- numden chemical potential
                  call twh_debroglie(GLB_GET,imolty,debroglie)
                  numdenchem(ibox,imolty) 
     &                 = -temperature
     &                 * dlog( vdum /(debroglie**3))
               else
c                 --- technically this would be an infinite chempot
c                 --- but set to zero instead
                  numdenchem(ibox,imolty) = 0.0d0
               endif
            endif
         enddo
      enddo

c     --- formatting to make this look good for most common cases
      if ( numboxes .eq. 1 ) then
         write(6,1215) (ibox,ibox=1,numboxes)
 1215    format('Averages',15x,'Units',1x,'Type',1x,1(6x,'Box ',i2))
      elseif ( numboxes .eq. 2 ) then
         write(6,1216) (ibox,ibox=1,numboxes)
 1216    format('Averages',15x,'Units',1x,'Type',1x,2(6x,'Box ',i2))
      else
         write(6,1217) (ibox,ibox=1,numboxes)
 1217    format('Averages',15x,'Units',1x,'Type'
     &        ,1x,MAXBOX(6x,'Box ',i2))
      endif

c     --- acvol is stored in cubic angstroms, output in cubic nm
      write(6,1206) (acvol(ibox)/1000.0d0,ibox=1,numboxes)
 1206 format(1x,'Volume',16x,' nm^3',6x,MAXBOX(1x,e11.5))
c     --- acvolsq is stored in hexic angstroms, output in hexic nm
      write(6,1207) (acvolsq(ibox)/(1000.0d0**2),ibox=1,numboxes)
 1207 format(1x,'Volume^2',14x,' nm^6',6x,MAXBOX(1x,e11.5))
         
      do imolty = 1, nmolty
         write(6,1201) imolty, (acnbox(ibox,imolty),ibox=1,numboxes)
      enddo
 1201 format(1x,'Molecule Number',7x,6x,i4,1x,MAXBOX(f12.3))
      
      write(6,1204) (molvol(ibox),ibox=1,numboxes)
 1204 format(1x,'Molar Volume',9x,'ml/mol',6x,MAXBOX(1x,e11.5))
      write(6,1205) (speden(ibox),ibox=1,numboxes)
 1205 format(1x,'Specific Density',6x,' g/ml',6x,MAXBOX(1x,f11.8))
      
      do imolty = 1, nmolty
         write(6,1203) imolty,(acdens(ibox,imolty),ibox=1,numboxes)
      enddo
 1203 format(1x,'Number Density',8x,' nm-3',1x,i4,1x,MAXBOX(f12.5))
      
      do imolty = 1, nmolty
         write(6,1212) imolty,(molfra(ibox,imolty),ibox=1,numboxes)
      enddo
 1212 format(1x,'Mole Fraction',9x,6x,i4,1x,MAXBOX(f12.7))
      
      do imolty = 1,nmolty
         write(6,1213) imolty,(acrg(ibox,imolty),ibox=1,numboxes)
      enddo
 1213 format(1x,'Radius of Gyration',4x,'    A',1x,i4,1x
     &     ,MAXBOX(f12.7))

      if (scp_eval_derivative) then
         write(6,1225) (dudl_lj(ibox),ibox=1,numboxes)
 1225    format(1x,'dU/d[lambda_lj]',5x,'[K]   ',6x,MAXBOX(1x,e11.5))
         write(6,1226) (dudl_c(ibox),ibox=1,numboxes)
 1226    format(1x,'dU/d[lambda_c] ',5x,'[K]   ',6x,MAXBOX(1x,e11.5))
      endif
      
c     --- output the pressure only if we ever measured this quantity
      if ( pcount .ne. 0 ) then
         do istress = istart,iend
            call twh_pcomp(GLB_GET,istress,keyname)
            write(6,1220) pressurestr,keyname
     &           ,(acstress(ibox,istress),ibox=1,numboxes)
         enddo
 1220    format(1x,a6,1x,a13,2x,'  kPa',6x,MAXBOX(1x,e11.5))
c        --- partial pressure
         do imolty = 1, nmolty
            write(6,1295) pressurestr,imolty,(molfra(ibox,imolty)
     &           *acstress(ibox,0),ibox=1,numboxes)
         enddo
 1295    format(1x,a6,1x,'p_i <x_i><p_v>',3x,'kPa'
     &        ,1x,i4,1x,MAXBOX(1x,e11.5))
      endif
c     --- compute the ideal total pressure 
      do ibox = 1,numboxes
c        --- set ideal total pressure
         idealp(ibox) = totnumden(ibox)*temperature*(1.0d27*8.314d0)
     &        /(6.02d23*1000.0d0)
      enddo
      write(6,1219) 'Ideal Pressure','kPa'
     &     ,(idealp(ibox),ibox=1,numboxes)
 1219 format(1x,a14,8x,2x,a3,6x,MAXBOX(1x,e11.5))
c     --- ideal partial pressure
      do imolty = 1, nmolty
         write(6,1294) imolty,(acdens(ibox,imolty)*temperature
     &        *(1.0d27*8.314d0)/(6.02d23*1000.0d0)
     &        ,ibox=1,numboxes)
      enddo
 1294 format(1x,'Ideal p_i <N/V>kT',7x,'kPa'
     &     ,1x,i4,1x,MAXBOX(1x,e11.5))

c     --- output the thermodynamic pressure only if we measured it
      keycode = AC_DU_DV
      tpcount = 0
      imolty = 0
      do ibox = 1,numboxes
         call twh_averages(lfailure,AVG_EXTRACT_COUNT,keycode,ibox
     &        ,imolty,value,incount)
         if ( lfailure ) return
         tpcount = tpcount + incount
      enddo
      if ( tpcount .ne. 0 ) then
c        --- extract and output the thermodynamic pressures
         do ibox = 1,numboxes
            call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &           ,imolty,value,incount)
            if ( lfailure ) return
c           --- convert dU/dV into kPa
            dudv(ibox) = 
     &           value*(8.314d0)*(1.0d30)/((6.02d23)*(1000d0))
c           --- virial pressure is the ideal bit minus du/dv
            tpressure(ibox) = idealp(ibox) - dudv(ibox)
         enddo
         write(6,1221) '-<dU/dV>','kPa'
     &        ,(-dudv(ibox),ibox=1,numboxes)
 1221    format(1x,a8,14x,2x,a3,6x,MAXBOX(1x,e11.5))
         write(6,1222) 'Thermodynamic Pressure','kPa'
     &        ,(tpressure(ibox),ibox=1,numboxes)
 1222    format(1x,a22,2x,a3,6x,MAXBOX(1x,e11.5))

c        --- thermodynamic partial pressure
         do imolty = 1, nmolty
            write(6,1251) imolty,(molfra(ibox,imolty)
     &           *tpressure(ibox),ibox=1,numboxes)
         enddo
 1251    format(1x,'Thermo p_i <x_1><p_t>',3x,'kPa'
     &        ,1x,i4,1x,MAXBOX(1x,e11.5))
      endif

c     --- energies (measured every step)
      do iener = 1,nener
         call twh_vname(GLB_GET,iener,keyname)
         write(6,1300) keyname,(acv(iener,ibox),ibox=1,numboxes)
      enddo
 1300 format(1x,a15,7x,'    K',6x,MAXBOX(e12.4))

c     --- chemical potentials
      do imolty = 1, nmolty
         if ( success(imolty) .ne. 0 ) then
            write(6,1210) 'u (Density)           ','K',imolty
     &           ,(numdenchem(ibox,imolty),ibox=1,numboxes)
            write(6,1210) 'u (NVT Insertion)     ','K',imolty
     &           ,(nvtchem(ibox,imolty),ibox=1,numboxes)
            write(6,1210) 'u (NpT Insertion)     ','K',imolty
     &           ,(nptchem(ibox,imolty),ibox=1,numboxes)
            write(6,1210) 'u (Den. + NVT Insert) ','K',imolty
     &           ,(nvtchem(ibox,imolty)+numdenchem(ibox,imolty)
     &           ,ibox=1,numboxes)
            write(6,1210) 'u (Den. + NpT Insert) ','K',imolty
     &           ,(nptchem(ibox,imolty)+numdenchem(ibox,imolty)
     &           ,ibox=1,numboxes)
            write(6,1210) 'u (Gibbs Total)       ','K',imolty
     &           ,(gibbschem(ibox,imolty),ibox=1,numboxes)
         endif
      enddo
 1210 format(1x,a22,4x,a1,1x,i4,1x,MAXBOX(f12.3)) 

c     --- free energy (from chemical potential)
      do ibox = 1,numboxes
         freeenergy(ibox) = 0.0d0
         do imolty = 1,nmolty
            freeenergy(ibox) = freeenergy(ibox) + 
     &           acnbox(ibox,imolty)*8.314d0*gibbschem(ibox,imolty)
     &           /1000.0d0
         enddo
      enddo
      write(6,1299) 'G: Sum{<u_i><N_i>}'
     &     ,(freeenergy(ibox),ibox=1,numboxes)
 1299 format(1x,a18,3x,'kJ/mol',6x,MAXBOX(e12.4))

c     --- potential energy and enthalpy
      do ibox = 1,numboxes
         ntotal(ibox) = 0.0d0
         do imolty = 1,nmolty
            ntotal(ibox) = ntotal(ibox) + acnbox(ibox,imolty)
         enddo
         uterm(ibox) = acv(1,ibox)*8.314d0/1000.0d0
         pvterm(ibox) = 6.02d23*acstress(ibox,0)*acvol(ibox)
     &        /1.0d30
         enthalpy(ibox) = (uterm(ibox) + pvterm(ibox))
         if ( ntotal(ibox) .gt. 0.0d0 ) then
            hpermolec(ibox) = enthalpy(ibox)/ntotal(ibox)
         else
c           --- set to zero
            hpermolec(ibox) = 0.0d0
         endif
      enddo

      write(6,1302) 'U',(uterm(ibox),ibox=1,numboxes)
 1302 format(1x,a1,20x,'kJ/mol',6x,MAXBOX(e12.4))

c     --- all of the following require pressure so only output
c     --- if we measured the pressure
      if ( pcount .ne. 0 ) then
         write(6,1304) 'pV: <p><V>',(pvterm(ibox),ibox=1,numboxes)
 1304    format(1x,a10,11x,'kJ/mol',6x,MAXBOX(e12.4))

c        --- enthalpy from the averages of u, p, and v
         write(6,1305) 'H: <U> + <p><V>'
     &        ,(enthalpy(ibox),ibox=1,numboxes)
 1305    format(1x,a15,6x,'kJ/mol',6x,MAXBOX(e12.4))

c        --- direct enthalpy averaged every time we measured the 
c        --- pressure
         write(6,1333) 'H: <U + pV>',(acenthalpy(ibox)
     &        ,ibox=1,numboxes)
 1333    format(1x,a11,10x,'kJ/mol',6x,MAXBOX(e12.4))
         
c        --- entropy (derived as difference from G and H)
         do ibox = 1,numboxes
            entropy(ibox) = beta*(enthalpy(ibox) - freeenergy(ibox))
         enddo
         write(6,1306) 'S: (<H> - <G>)/T'
     &        ,(entropy(ibox),ibox=1,numboxes)
 1306    format(1x,a16,3x,'kJ/K mol',6x,MAXBOX(e12.4))

c        --- Compression factor
         do ibox = 1,numboxes
            if ( ntotal(ibox) .gt. 0.0d0 ) then
               compression(ibox) = beta*pvterm(ibox)*1000.0d0
     &              /(8.314d0*ntotal(ibox))
            else
               compression(ibox) = 0.0d0
            endif
         enddo
         write(6,1320) 'Z: <p><V>/<N>RT'
     &        ,(compression(ibox),ibox=1,numboxes)
 1320    format(1x,a15,7x,11x,MAXBOX(f12.6))
      endif

      if ( numboxes .eq. 1 .and. ensemble .eq. ENS_NVT ) then
c        --- 1 box canonical ensemble constant volume heat capacity
         c_v = beta*beta*(acvboxsq(1) - acv(1,1)**2)*8.314d0
         write(6,1325) 'C_v: T^2<H^2>-<H>^2',c_v
      elseif ( numboxes .eq. 1 .and. ensemble .eq. ENS_NPT 
     &        .and. pcount .ne. 0 ) then
c        --- 1 box isobaric-isothermal ensemble constant pressure 
c        --- heat capacity
         c_p = beta*beta*1000.0d0
     *        *(acenthalpysq(1) - acenthalpy(1)**2)/8.314d0
         write(6,1325) 'C_p: T^2<H^2>-<H>^2',c_p
 1325    format(1x,a19,1x,'J/K mol',6x,e12.4)
      endif

c     --- enthalpy of vaporization (Direct)
      keycode = AC_H_VAP_DIRECT
      imolty = 1
      ibox = 0
      call twh_averages(lfailure,AVG_EXTRACT_COUNT,keycode,ibox,imolty
     &     ,value,outcount)
      if ( lfailure ) return
      if ( outcount .gt. 0 ) then
         call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &        ,imolty,value,outcount)
         if ( lfailure ) return
         write(6,'(1x,a14,7x,a6,6x,e12.4)') 'H_vap (Direct)'
     &        ,'kJ/mol',value
      endif
c     --- enthalpy of vaporization (vapor p)
      keycode = AC_H_VAP_VAPOR_P
      imolty = 1
      ibox = 0
      call twh_averages(lfailure,AVG_EXTRACT_COUNT,keycode,ibox,imolty
     &     ,value,outcount)
      if ( lfailure ) return
      if ( outcount .gt. 0 ) then
         call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &        ,imolty,value,outcount)
         if ( lfailure ) return
         write(6,'(1x,a15,6x,a6,6x,e12.4)') 'H_vap (vapor p)'
     &        ,'kJ/mol',value
      endif
c     --- enthalpy of vaporization (pV/n=RT)
      keycode = AC_H_VAP_PV_NRT
      imolty = 1
      ibox = 0
      call twh_averages(lfailure,AVG_EXTRACT_COUNT,keycode,ibox,imolty
     &     ,value,outcount)
      if ( lfailure ) return
      if ( outcount .gt. 0 ) then
         call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &        ,imolty,value,outcount)
         if ( lfailure ) return
         write(6,'(1x,a15,6x,a6,6x,e12.4)') 'H_vap (pV/n=RT)'
     &        ,'kJ/mol',value
      endif
c     --- heat of vaporization (thermodynamic pressure and final avg)
      if ( nmolty .eq. 1 .and. numboxes .eq. 2 ) then
         if (speden(1) .lt. speden(2)) then
            ivap = 1
            iliq = 2
         else
            ivap = 2
            iliq = 1
         endif
c        --- uterm was computed above
         value = (uterm(ivap)/acnbox(ivap,1))
     &        - (uterm(iliq)/acnbox(iliq,1))
     &        + (6.02d23/1.0d30)*tpressure(ivap)
     &        *( (acvol(ivap)/acnbox(ivap,1))
     &        - (acvol(iliq)/acnbox(iliq,1)) )
         write(6,'(1x,a16,5x,a6,6x,e12.4)') 'H_vap (thermo p)'
     &        ,'kJ/mol',value
      endif

c     --- chemical potential of isolated molecule
      do imolty = 1, nmolty
         call twh_nmolectyp(GLB_GET,imolty,nmtype)
         if ( nmtype .eq. 0 ) then
            ibox = MAXBOX
            jbox = 0
c           --- get the averages for isolated chemical potential
            keycode = AC_ISOLATION_CHEMPOT
            call twh_averages(lfailure,AVG_EXTRACT_AVERAGE,keycode,ibox
     &           ,imolty,value,incount)
            if ( lfailure ) return
            isochem(imolty) = value
c           --- compute the number of successful growths into ibox 
c           --- from a box other than ibox as we cannot compute chem 
c           --- pot for ibox to ibox
            success(imolty) = 0
            call twh_bscbswap(GLB_GET,imolty,ibox,jbox,stemp)
            success(imolty) = success(imolty) + stemp
            if ( success(imolty) .ne. 0 ) then
c              --- transform the chemical potential
               isochem(imolty) = -temperature
     &              *dlog(isochem(imolty))
                  
               write(6,'(1x,a13,9x,4x,a1,1x,i4,1x,f12.3)') 
     &              'u (Isolation)','K',imolty,isochem(imolty)
               do ibox = 1,numboxes
                  if ( ensemble .eq. ENS_NVT ) then
                     residchem(ibox,imolty) = nvtchem(ibox,imolty)
     &                    - isochem(imolty)
                  elseif ( ensemble .eq. ENS_NPT ) then
                     residchem(ibox,imolty) = nptchem(ibox,imolty)
     &                    - isochem(imolty)
                  endif
                  write(6,1210) 'u (Residual)          ','K',imolty
     &                 ,(residchem(jbox,imolty),jbox=1,numboxes)
               enddo
            endif
         endif
      enddo
c     --- Henry's Law
      do imolty = 1, nmolty
         ibox = 1
         call twh_nmolectyp(GLB_GET,imolty,nmtype)
         if ( nmtype .eq. 0 ) then
c           --- we want the residual chemical potential
            henry_coeff = totnumden(ibox)*temperature*henry_convert
     &           *dexp(residchem(ibox,imolty)*beta)
            write(6,'(1x,a20,1x,3x,a3,1x,i4,1x,e12.5)') 
     &           'Henry Law (Residual)','MPa',imolty,henry_coeff
         endif
      enddo
      
c     --- get the final number of blocks
      call twh_averages(lfailure,AVG_EXTRACT_NBLOCK,keycode,ibox,imolty
     &     ,value,nblock)
      if ( lfailure ) return

c     --- calculate statistical errors 
      if ( nblock .ge. 2 ) then
         write(6,*)
         write(6,1301) nblock
 1301    format('Block Averages (',i3,' blocks)',1x,'Units',1x
     &        ,'Type',1x,'Box',1x,'Average',6x,'Standard Deviation')
c        --- specific density
         keycode = AC_SPECIFIC_DENSITY
         imolty = 0
         do ibox = 1, numboxes
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &           ,ibox,imolty,value,nblock)
            if ( lfailure ) return
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &           ,ibox,imolty,stdev,nblock)
            if ( lfailure ) return
            write(6,1331) ibox,value,stdev
         enddo
 1331    format(' Specific Density',11x,' g/ml',1x,4x,1x,i3
     &        ,2(1x,e12.5))
         if ( pcount .ne. 0 ) then
c           --- virial stress tensor
            do ibox = 1, numboxes
c              --- total pressure
               do istress = istart,iend
                  keycode = AC_PRESSURE - istress
                  call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG
     &                 ,keycode,ibox,imolty,value,nblock)
                  if ( lfailure ) return
                  call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD
     &                 ,keycode,ibox,imolty,stdev,nblock)
                  if ( lfailure ) return
                  call twh_pcomp(GLB_GET,istress,keyname)
                  write(6,1341) pressurestr,keyname,ibox,value,stdev
               enddo
            enddo
 1341       format(1x,a6,1x,a13,7x,'  kPa',1x,4x,1x,i3,2(1x,e12.5))
         endif

         if ( tpcount .ne. 0 ) then
c           --- thermodynamic pressure is the product of the total
c           ---- number density and du/dv
c           --- transform the block averages of these quantities into 
c           --- a new block average
            imolty = 0
            do ibox = 1,numboxes
               do iblock = 1,nblock
c                 --- total number density
                  keycode = AC_NUMBER_DENSITY
                  call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL
     &                 ,keycode,ibox,imolty,baver(1),iblock)
                  if ( lfailure ) return
c                 --- du/dv
                  keycode = AC_DU_DV
                  call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL
     &                 ,keycode,ibox,imolty,baver(2),iblock)
                  if ( lfailure ) return
c                 --- compute the thermodynamic pressure for this block
                  value = baver(1)*temperature
     &                 *(1.0d27*8.314d0)/(6.02d23*1000.0d0)
     &                 -baver(2)*(8.314d0)*(1.0d30)/((6.02d23)*(1000d0))
c                 --- reset the block value
                  keycode = AC_THERMO_PRESSURE
                  call twh_averages(lfailure,AVG_SET_BLOCK_VAL,keycode
     &                 ,ibox,imolty,value,iblock)
                  if ( lfailure ) return
               enddo
            enddo
c           --- now extract the thermodynamic pressure block averages
            do ibox = 1, numboxes
               keycode = AC_THERMO_PRESSURE
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &              ,ibox,imolty,value,nblock)
               if ( lfailure ) return
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &              ,ibox,imolty,stdev,nblock)
               if ( lfailure ) return
               write(6,100) 'Thermodynamic Pressure',ibox,value,stdev
            enddo
 100        format(1x,a22,7x,'kPa',1x,4x,1x,i3,2(1x,e12.5))
         endif
c        --- energies
         do ibox = 1, numboxes
            do iener = 1,nener
               keycode = iener
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &              ,ibox,imolty,value,nblock)
               if ( lfailure ) return
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &              ,ibox,imolty,stdev,nblock)
               if ( lfailure ) return
               call twh_vname(GLB_GET,iener,keyname)
               write(6,1311) keyname,ibox,value,stdev
            enddo
         enddo
 1311    format(1x,a15,12x,'    K',1x,4x,1x,i3,2(1x,e12.5))

c        --- chemical potentials are stored in a different format 
c        --- then the actual final answer.  Transform the individual
c        --- blocks before computing the block averages and stdevs
         do imolty = 1,nmolty
            call twh_nmolectyp(GLB_GET,imolty,nmtype)
            if ( nmtype .eq. 0 ) then
               ibox = MAXBOX
               do iblock = 1,nblock
c                 --- isolation chemical potential
                  keycode = AC_ISOLATION_CHEMPOT
                  call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL
     &                 ,keycode,ibox,imolty,value,iblock)
                  if ( lfailure ) return
c                 --- transform average into chemical potential
                  if ( value .ne. 0.0d0 ) then
                     value = -temperature*dlog(value)
                  endif
c                 --- reset the block value
                  call twh_averages(lfailure,AVG_SET_BLOCK_VAL,keycode
     &                 ,ibox,imolty,value,iblock)
                  if ( lfailure ) return
               enddo
            endif
            do ibox = 1,numboxes
               do iblock = 1,nblock
c                 --- Gibbs total chemical potential
                  keycode = AC_GIBBS_TOT_CHEMPOT
                  call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL
     &                 ,keycode,ibox,imolty,value,iblock)
                  if ( lfailure ) return
c                 --- transform average into chemical potential
                  if ( value .ne. 0.0d0 ) then
                     value = -temperature*dlog(value)
                  endif
c                 --- reset the block value
                  call twh_averages(lfailure,AVG_SET_BLOCK_VAL,keycode
     &                 ,ibox,imolty,value,iblock)
                  if ( lfailure ) return
c                 --- NpT insertion chemical potential
                  keycode = AC_NPT_INSERTION_CHEMPOT
                  call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL
     &                 ,keycode,ibox,imolty,value,iblock)
                  if ( lfailure ) return
c                 --- also need the average volume for this block
                  keycode = AC_VOLUME
                  call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL
     &                 ,keycode,ibox,0,voltemp,iblock)
                  if ( lfailure ) return
c                 --- transform average into chemical potential
                  if ( voltemp .ne. 0.0d0 .and. value .ne. 0.0d0 ) then
                     value = -temperature*dlog(value/voltemp)
                  else
                     value = 0.0d0
                  endif
c                 --- reset the block value
                  keycode = AC_NPT_INSERTION_CHEMPOT
                  call twh_averages(lfailure,AVG_SET_BLOCK_VAL,keycode
     &                 ,ibox,imolty,value,iblock)
                  if ( lfailure ) return
c                 --- NVT insertion chemical potential
                  keycode = AC_NVT_INSERTION_CHEMPOT
                  call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL
     &                 ,keycode,ibox,imolty,value,iblock)
                  if ( lfailure ) return
c                 --- transform average into chemical potential
                  if ( value .ne. 0.0d0 ) then
                    value = -temperature*dlog(value)
                  endif
c                 --- reset the block value
                  call twh_averages(lfailure,AVG_SET_BLOCK_VAL,keycode
     &                 ,ibox,imolty,value,iblock)
                  if ( lfailure ) return
               enddo
            enddo
         enddo

c        --- gibbs total configurational-bias chemical potential
         keycode = AC_GIBBS_TOT_CHEMPOT
         do imolty = 1, nmolty
            do ibox = 1, numboxes
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &              ,ibox,imolty,value,nblock)
               if ( lfailure ) return
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &              ,ibox,imolty,stdev,nblock)
               if ( lfailure ) return
               if ( value .ne. 0.0d0 ) then
                  write(6,1351) imolty,ibox,value,stdev
               endif
            enddo
         enddo
 1351    format(' u (Gibbs Total)',12x,'    K',1x,i4,1x,i3,2(1x,f12.3))

c        --- NpT insertion configurational-bias chemical potential
         keycode = AC_NPT_INSERTION_CHEMPOT
         do imolty = 1, nmolty
            do ibox = 1, numboxes
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &              ,ibox,imolty,value,nblock)
               if ( lfailure ) return
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &              ,ibox,imolty,stdev,nblock)
               if ( lfailure ) return
               if ( value .ne. 0.0d0 ) then
                  write(6,1352) imolty,ibox,value,stdev
               endif
            enddo
         enddo
 1352    format(' u (NpT Insertion)',10x,'    K',1x,i4,1x,i3
     &        ,2(1x,f12.3))

c        --- NVT insertion configurational-bias chemical potential
         keycode = AC_NVT_INSERTION_CHEMPOT
         do imolty = 1, nmolty
            do ibox = 1, numboxes
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &              ,ibox,imolty,value,nblock)
               if ( lfailure ) return
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &              ,ibox,imolty,stdev,nblock)
               if ( lfailure ) return
               if ( value .ne. 0.0d0 ) then
                  write(6,1353) imolty,ibox,value,stdev
               endif
            enddo
         enddo
 1353    format(' u (NVT Insertion)',10x,'    K',1x,i4,1x,i3
     &        ,2(1x,f12.3))

c        --- number density
         keycode = AC_NUMBER_DENSITY
         do imolty = 1, nmolty
            do ibox = 1, numboxes
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &              ,ibox,imolty,value,nblock)
               if ( lfailure ) return
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &              ,ibox,imolty,stdev,nblock)
               if ( lfailure ) return
               write(6,1361) imolty,ibox,value,stdev
            enddo
         enddo
 1361    format(' Number Density',13x,' nm-3',1x,i4,1x,i3,2(1x,e12.5))

c        --- henry law inferred from the insertion chempot and total
c        --- number density
         if ( numboxes .eq. 1 .and. ensemble .ne. ENS_muVT ) then
            ibox = 1
c           --- compute henry law coeffs for each block from the blocks
c           --- of number density (computed above) and chemical pot.
            do imolty = 1,nmolty
               call twh_nmolectyp(GLB_GET,imolty,nmtype)
               if ( nmtype .eq. 0 ) then
                  lhenry_good = .true.
                  keycode = AC_HENRY_LAW
                  call twh_averages(lfailure,AVG_ASSIGN,keycode,ibox
     &                 ,imolty,value,incount)
                  if ( lfailure ) return
                  do iblock = 1,nblock
c                    --- get the total number density for this block
                     totnumden(ibox) = 0.0d0
                     do jmolty = 1,nmolty
                        keycode = AC_NUMBER_DENSITY
                        call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL
     &                       ,keycode,ibox,jmolty,value,iblock)
                        if ( lfailure ) return
                        totnumden(ibox) = totnumden(ibox) + value
                     enddo
c                    --- set the blocks from the block values
                     if ( ensemble .eq. ENS_NVT ) then
                        keycode = AC_NVT_INSERTION_CHEMPOT
                     elseif ( ensemble .eq. ENS_NPT ) then
                        keycode = AC_NPT_INSERTION_CHEMPOT
                     else
                        write(6,*) 'WRITEAVERAGES: unknown ensemble in'
     &                       ,' Henry law calc',ensemble
                        stop
                     endif
c                    --- get the insertion chemical potential
                     call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL
     &                    ,keycode,ibox,imolty,chempotential,iblock)
                     if ( lfailure ) return
c                    --- get the isolation chemical potential
                     keycode = AC_ISOLATION_CHEMPOT
                     call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL
     &                    ,keycode,MAXBOX,imolty,value,iblock)
                     if ( lfailure ) return
c                    --- subtract the Isolation chemical potential
                     chempotential = chempotential - value
                     if ( chempotential .ne. 0.0d0 ) then
c                       --- combine with the number density in the box
                        henry_coeff = totnumden(ibox)*temperature
     &                       *henry_convert*dexp(chempotential*beta)
c                       --- set the block value to this answer
                        keycode = AC_HENRY_LAW
                        call twh_averages(lfailure,AVG_SET_BLOCK_VAL
     &                       ,keycode,ibox,imolty,henry_coeff,iblock)
                        if ( lfailure ) return
                     else
c                       --- bad block
                        lhenry_good = .false.
                     endif
                  enddo
c                 --- now extract the average and standard deviation
c                 --- using the blocks we just assigned
                  if ( lhenry_good ) then
                     keycode = AC_HENRY_LAW
                     call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG
     &                    ,keycode,ibox,imolty,value,nblock)
                     if ( lfailure ) return
                     call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD
     &                    ,keycode,ibox,imolty,stdev,nblock)
                     if ( lfailure ) return
                     write(6,1355) imolty,ibox,value,stdev
                  endif
               endif
            enddo
         endif
 1355    format(' Henry Coeff (Residual)',5x,'  MPa',1x,i4,1x
     &        ,i3,2(1x,e12.5))

c        --- molfraction
         keycode = AC_MOL_FRACTION
         do imolty = 1, nmolty
            do ibox = 1, numboxes
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &              ,ibox,imolty,value,nblock)
               if ( lfailure ) return
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &              ,ibox,imolty,stdev,nblock)
               if ( lfailure ) return
               write(6,1371) imolty,ibox,value,stdev
            enddo
         enddo
 1371    format(' Mole Fraction',14x,5x,1x,i4,1x,i3,2(1x,f12.7))
c        --- molarity (from number density)
         keycode = AC_NUMBER_DENSITY
         do imolty = 1, nmolty
            do ibox = 1, numboxes
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &              ,ibox,imolty,value,nblock)
               if ( lfailure ) return
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &              ,ibox,imolty,stdev,nblock)
               if ( lfailure ) return
               write(6,1380) imolty,ibox,1d1*value/6.02d0
     &              ,1.0d1*stdev/6.02d0
            enddo
         enddo
 1380    format(' Molarity',19x,'    M',1x,i4,1x,i3,2(1x,e12.5))

c        --- radius of gyration
         keycode = AC_RADIUS_OF_GYRATION
         do imolty = 1, nmolty
            do ibox = 1, numboxes
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &              ,ibox,imolty,value,nblock)
               if ( lfailure ) return
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &              ,ibox,imolty,stdev,nblock)
               if ( lfailure ) return
               write(6,1385) imolty,ibox,value,stdev
            enddo
         enddo
 1385    format(' Radius of Gyration',9x,'    A',1x,i4,1x,i3
     &        ,2(1x,f12.5))

c        --- dU/d[lambda]
         if (scp_eval_derivative) then
            keycode = AC_DU_DLAMBDA_LJ
            imolty = 0
            do ibox = 1, numboxes
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &              ,ibox,imolty,value,nblock)
               if ( lfailure ) return
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &              ,ibox,imolty,stdev,nblock)
               if ( lfailure ) return
               write(6,1386) ibox,value,stdev
            enddo

            keycode = AC_DU_DLAMBDA_C
            do ibox = 1, numboxes
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &              ,ibox,imolty,value,nblock)
               if ( lfailure ) return
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &              ,ibox,imolty,stdev,nblock)
               if ( lfailure ) return
               write(6,1387) ibox,value,stdev
            enddo
         endif
 1386    format(' dU/d[lambda_lj] ',11x,'  [K]',1x,4x,1x,i3
     &        ,2(1x,e12.5))
 1387    format(' dU/d[lambda_lj] ',11x,'  [K]',1x,4x,1x,i3
     &        ,2(1x,e12.5))
         
c        --- heat of vaporization (Direct)
         keycode = AC_H_VAP_DIRECT
         ibox = 0
         imolty = 1
c        --- see if this was ever computed
         call twh_averages(lfailure,AVG_EXTRACT_COUNT,keycode,ibox
     &        ,imolty,value,outcount)
         if ( lfailure ) return
         if ( outcount .gt. 0 ) then
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &           ,ibox,imolty,value,nblock)
            if ( lfailure ) return
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &           ,ibox,imolty,stdev,nblock)
            if ( lfailure ) return
            write(6,'(1x,a14,12x,a6,9x,2(1x,e12.4))') 'H_vap (Direct)'
     &           ,'kJ/mol',value,stdev
         endif
c        --- heat of vaporization (vapor p)
         keycode = AC_H_VAP_VAPOR_P
         ibox = 0
         imolty = 1
c        --- see if this was ever computed
         call twh_averages(lfailure,AVG_EXTRACT_COUNT,keycode,ibox
     &        ,imolty,value,outcount)
         if ( lfailure ) return
         if ( outcount .gt. 0 ) then
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &           ,ibox,imolty,value,nblock)
            if ( lfailure ) return
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &           ,ibox,imolty,stdev,nblock)
            if ( lfailure ) return
            write(6,'(1x,a15,11x,a6,9x,2(1x,e12.4))') 'H_vap (vapor p)'
     &           ,'kJ/mol',value,stdev
         endif
c        --- heat of vaporization (pV/n=RT)
         keycode = AC_H_VAP_PV_NRT
         ibox = 0
         imolty = 1
c        --- see if this was ever computed
         call twh_averages(lfailure,AVG_EXTRACT_COUNT,keycode,ibox
     &        ,imolty,value,outcount)
         if ( lfailure ) return
         if ( outcount .gt. 0 ) then
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &           ,ibox,imolty,value,nblock)
            if ( lfailure ) return
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &           ,ibox,imolty,stdev,nblock)
            if ( lfailure ) return
            write(6,'(1x,a15,11x,a6,9x,2(1x,e12.4))') 'H_vap (pV/n=RT)'
     &           ,'kJ/mol',value,stdev
         endif

c       --- write block averages 
         write(6,*)
         write(6,*) '-----block averages ------'
         do ibox = 1, numboxes
            write(6,1403) ibox,pressurestr
 1403       format('Box: ' ,i4,/,'Block Energy',9x
     &           ,'Density',8x,a6,1x,'Press.',2x,'Mol fracs') 
            do iblock = 1, nblock
               imolty = 0
c              --- first energy term
               keycode = 1
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL,keycode
     &              ,ibox,imolty,baver(1),iblock)
               if ( lfailure ) return
               keycode = AC_SPECIFIC_DENSITY
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL,keycode
     &              ,ibox,imolty,baver(2),iblock)
               if ( lfailure ) return
               if ( pressurestyle .eq. PRES_VIRIAL .or.
     &              pressurestyle .eq. PRES_RADIAL ) then
                  keycode = AC_PRESSURE
               elseif ( pressurestyle .eq. PRES_THERMO ) then
                  keycode = AC_THERMO_PRESSURE
               endif
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL,keycode
     &              ,ibox,imolty,baver(3),iblock)
               if ( lfailure ) return
               keycode = AC_MOL_FRACTION
               do imolty = 1,nmolty
                  call twh_averages(lfailure,AVG_EXTRACT_BLOCK_VAL
     &                 ,keycode,ibox,imolty,baver(3+imolty),iblock)
                  if ( lfailure ) return
               enddo
               write(6,1402) iblock,(baver(zz),zz=1,3+nmolty)
            enddo
         enddo
 1402    format(i5,1x,3(e14.8,1x), NTMAX (f10.8,1x))

c        --- Special output for constructing phase diagrams, etc.
         if ( numboxes .eq. 2 .and. nmolty .eq. 1 ) then
c           --- write out information for constructing the single 
c           --- component vapor-liquid coexistence curve
            if (speden(1) .lt. speden(2)) then
               ivap = 1
               iliq = 2
            else
               ivap = 2
               iliq = 1
            endif
            filename = 'towhee_vlcc'
            call twh_openfile(lfailure,41,filename,'DEFAULT','UNKNOWN'
     &           ,'DEFAULT')
            if ( lfailure ) return
c           --- vapor specific density, vapor std, liquid specific den
c           --- liquid std, vapor pressure, vapor p std
            keycode = AC_SPECIFIC_DENSITY
            imolty = 0
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &           ,ivap,imolty,baver(1),nblock)
            if ( lfailure ) return
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &           ,iliq,imolty,baver(3),nblock)
            if ( lfailure ) return
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &           ,ivap,imolty,baver(2),nblock)
            if ( lfailure ) return
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &           ,iliq,imolty,baver(4),nblock)
            if ( lfailure ) return
            if ( pressurestyle .eq. PRES_VIRIAL .or.
     &           pressurestyle .eq. PRES_RADIAL ) then
c              --- virial or radial pressure
               keycode = AC_PRESSURE
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &              ,ivap,imolty,baver(5),nblock)
               if ( lfailure ) return
               call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &              ,ivap,imolty,baver(6),nblock)
               if ( lfailure ) return
            endif
c           --- thermodynamic pressure
            keycode = AC_THERMO_PRESSURE
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &           ,ivap,imolty,baver(7),nblock)
            if ( lfailure ) return
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &           ,ivap,imolty,baver(8),nblock)
            if ( lfailure ) return

            write(41,1401) temperature,(baver(zz),zz=1,8)
 1401       format(g20.10,8(1x,e12.6))
            close(41)
         endif

c        --- write ostwald values for each moltyp
         filename = 'towhee_ostwald'
         call twh_openfile(lfailure,42,filename,'DEFAULT','UNKNOWN'
     &        ,'DEFAULT')
         if ( lfailure ) return
         gconst = 8.314/(1000*beta)
         do ibox = 1,numboxes-1
            do jbox = ibox+1,numboxes
               write(42,*) 'Boxes ',ibox,' and ',jbox,' Partition K'
     &              ,' std dev   Delta G   std dev'
               do imolty = 1,nmolty
c                 --- baver stores num den(ibox), std(ibox) 
c                 --- num den(jbox) std(jbox)
                  keycode = AC_NUMBER_DENSITY
                  call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG
     &                 ,keycode,ibox,imolty,baver(1),nblock)
                  if ( lfailure ) return
                  call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG
     &                 ,keycode,jbox,imolty,baver(3),nblock)
                  if ( lfailure ) return
                  call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD
     &                 ,keycode,ibox,imolty,baver(2),nblock)
                  if ( lfailure ) return
                  call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD
     &                 ,keycode,jbox,imolty,baver(4),nblock)
                  if ( lfailure ) return
                  if ( baver(1)*baver(3) .ne. 0.0d0 ) then 
                     ostwald = baver(1)/baver(3)
                     stdost  = ostwald*dsqrt( (baver(2)/baver(1))**2
     &                    +  (baver(4)/baver(3))**2 )
                     call twh_nunit(GLB_GET,imolty,nunit)
                     write(42,*) imolty,nunit,ostwald,stdost
     &                    ,-(gconst*log(ostwald)),gconst*stdost/ostwald
                  endif
               enddo
            enddo
         enddo
         close(42)

c        --- information for analyzing pressure vs. density isotherms
         if ( numboxes .eq. 1 .and. ensemble .eq. ENS_NVT ) then
            filename = 'towhee_pvsrho'
            call twh_openfile(lfailure,41,filename,'DEFAULT','UNKNOWN'
     &           ,'DEFAULT')
            if ( lfailure ) return
            imolty = 0
            keycode = AC_SPECIFIC_DENSITY
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &           ,ibox,imolty,baver(1),nblock)
            if ( lfailure ) return
            if ( pressurestyle .eq. PRES_VIRIAL .or.
     &           pressurestyle .eq. PRES_RADIAL ) then
               keycode = AC_PRESSURE
            elseif ( pressurestyle .eq. PRES_THERMO ) then
               keycode = AC_THERMO_PRESSURE
            endif
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_AVG,keycode
     &           ,ibox,imolty,baver(2),nblock)
            if ( lfailure ) return
            call twh_averages(lfailure,AVG_EXTRACT_BLOCK_STD,keycode
     &           ,ibox,imolty,baver(3),nblock)
            if ( lfailure ) return
            write(41,*) baver(1),baver(2),baver(3)
            close(41)
         endif
      endif

      return
      end
