#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2013 Marcus G. Martin                       *
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
      subroutine twh_setclassical(lfailure,numtypes)
c     ******************************************************************
c     * sets up many variables for classical potentials                *
c     * otherwise combinations are only computed if lhere is true      *
c     *                                                                *
c     * split from readtowhee 11-01-2004 by M.G. Martin                *
c     * last modified 07-30-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_CLASSICAL_MIXRULE
#define FUNCTION_GET_LHERE
#define FUNCTION_GET_LOFHERE
#define FUNCTION_GET_NBCOEFF
#define FUNCTION_GET_OFCOEFF
#define FUNCTION_GET_OFSTYLE
#define FUNCTION_GETNBTYPE
#define FUNCTION_LIMITED_DOUBLE
#define FUNCTION_VONEFIVE
#include "functions.h"

c     --- variables passed to/from the subroutine
      logical lfailure
      integer numtypes
c     --- local variables
c     --- logical scalars
      logical lovrlap,lshift,ltailc,lcoulomb
c     --- character string scalars
      character*20 eampairstyle,eamdensstyle,keyname
c     --- integer scalars
      integer imolty,iunit,itype,jtype,ii,ntot,cross,iutemp,ij,cutstyle
      integer cpnumber,nmolty,nunit,nmtype,kmax,nener
c     --- double precision scalars
      double precision totalcharge,qsum,vdum,rijsq,swcut,dist,dtemp
      double precision charge,tpair,eamdens,rcut,rcutsq,rcutin,rcutinsq
      double precision minboxdim
c     --- initialize
      lovrlap = .false.

c     --- retrieve constants
      call twh_cp_number(GLB_GET,cpnumber)
      call twh_ltailc(GLB_GET,ltailc)
      call twh_lcoulomb(GLB_GET,lcoulomb)

c     --- get the initial value of nener
      call twh_nener(GLB_GET,nener)
c     --- set the names of the energy terms
      nener = nener + 1
      keyname = 'Total Classical'
      call twh_vname(GLB_SET,nener,keyname)
      nener = nener + 1
      keyname = 'Inter vdw'
      call twh_vname(GLB_SET,nener,keyname)
      nener = nener + 1
      keyname = 'Angle'
      call twh_vname(GLB_SET,nener,keyname)
      nener = nener + 1
      keyname = 'Torsion'
      call twh_vname(GLB_SET,nener,keyname)
      nener = nener + 1
      keyname = 'Intra vdw'
      call twh_vname(GLB_SET,nener,keyname)
      nener = nener + 1
      keyname = 'External Field'
      call twh_vname(GLB_SET,nener,keyname)
      nener = nener + 1
      keyname = 'Vibration'
      call twh_vname(GLB_SET,nener,keyname)
      nener = nener + 1
      keyname = 'Coulomb'
      call twh_vname(GLB_SET,nener,keyname)
      nener = nener + 1
      keyname = 'Tail vdw'
      call twh_vname(GLB_SET,nener,keyname)
      nener = nener + 1
      keyname = 'Solvation'
      call twh_vname(GLB_SET,nener,keyname)
c     --- set the final value of nener
      call twh_nener(GLB_SET,nener)

c     --- set up the tail correction information (if needed)
      if ( ltailc ) call twh_tail()

c     --- check the charge on molecules and in system
      totalcharge = 0.0d0
      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
c        --- initialize qsum
         qsum = 0.0d0
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            call twh_qqatom(GLB_GET,imolty,iunit,charge)
            if ( .not. lcoulomb .or. charge .eq. 0.0d0 ) then
c               --- no coulmbic interaction for this atom type
               call twh_luse_c(GLB_SET,imolty,iunit,.false.)
            else
c              --- compute the total charge for this molecule
               qsum = qsum + charge
            endif

            call twh_ntype(GLB_GET,imolty,iunit,iutemp)
c           --- set luse logicals for coulomb & VdW
            if ( cpnumber .eq. CP_LENNARD_JONES
     &           .or. cpnumber .eq. CP_SCALED_LENNARD_JONES
     &           ) then
               ii = twh_getnbtype(iutemp,iutemp)
               if ((twh_get_nbcoeff(1,ii) .eq. 0.0d0) .and.
     &              (twh_get_classical_mixrule() .ne. 'Explicit')) then
c                 --- zero energy from the Lennard-Jones term
                  call twh_luse_v(GLB_SET,imolty,iunit,.false.)
               endif
            endif
         enddo

         write(6,1025) imolty,twh_limited_double(qsum,1.0d-10)
 1025    format('Default total charge on molecule ',i3,' is ',f9.5)
c        --- add up total charge in the simulation
         call twh_nmolectyp(GLB_GET,imolty,nmtype)
         totalcharge = totalcharge + qsum*dble(nmtype)
      enddo
      write(6,'(a38,f10.5)') 'Total charge in the simulation system:'
     &     ,twh_limited_double(totalcharge,1.0d-10)
      if ( totalcharge**2 .gt. 1.0d-5 ) then
         write(6,*) 'SETCLASSICAL: total system charge is not zero'
         call twh_ewald_kmax(GLB_GET,1,kmax)
         if ( kmax .ne. 0 ) then
            write(6,*) 'Ewald sum does not work properly for this case'
            lfailure = .true.
            return
         endif
      endif

c     --- output the stretching, bending, and torsional constants
      call twh_writebond(6)
      call twh_writeangle(6)
      call twh_writeangang(6)
      call twh_writetorsion(6)
      call twh_writeimproper(6)
c     --- setup the shift constants for onefive
      cutstyle = -1
      call twh_rcut(GLB_GET,rcut)
      rijsq = rcut*rcut
c     --- retrieve the shift logical
      call twh_lshift(GLB_GET,lshift)
      do itype = 1,TOFMAX
         if ( twh_get_lofhere(itype) ) then
            if ( twh_get_ofstyle(itype) .eq. 1 ) then
c              --- store the modified parameters
               dtemp = twh_get_ofcoeff(itype,1)**6
               call twh_ofcoeff(GLB_SET,itype,3,dtemp)
               dtemp = twh_get_ofcoeff(itype,2)*4.0d0
               call twh_ofcoeff(GLB_SET,itype,4,dtemp)
               if ( lshift ) then
c                 --- setup the shift constant
                  vdum = twh_vonefive(lfailure,rijsq,itype,cutstyle
     &                 ,lovrlap)
                  if ( lfailure ) return
                  if ( lovrlap ) then
                     write(6,*) 
     &                    'SETCLASSICAL: problem setting up one-five'
                     return
                  endif
c                 --- reset the shift
                  call twh_ofcoeff(GLB_SET,itype,0,vdum)
               endif
            elseif ( twh_get_ofstyle(itype) .eq. 2 ) then
c              --- this is not valid with a shifted potential
               if ( lshift ) then
                  write(6,*) 'one-five style of 2 and lshift are not'
     &                 ,' compatable'
                  lfailure = .true.
                  return
               endif
c              --- otherwise nothing to be done here
            else
c              --- unknown one-five style
               write(6,*) 'SETCLASSICAL: unknown one-five style'
               lfailure = .true.
               return
            endif
         endif
      enddo
c     --- output any one-five information
      call twh_writeonefive(6)

c     --- set minboxdim
      if ( cpnumber .eq. CP_HARD_SPHERE .or.
     &     cpnumber .eq. CP_REPULSIVE_SPHERE ) then
c        --- find the largest diameter, or at least 1.0
         minboxdim = 1.0d0
         do itype = 1,numtypes
            do jtype = itype,numtypes
               if ( twh_get_lhere(itype) .and. twh_get_lhere(jtype)
     &              ) then
                  ij = twh_getnbtype(itype,jtype)
                  if ( twh_get_nbcoeff(1,ij) .gt. minboxdim ) then
                     minboxdim = twh_get_nbcoeff(1,ij)
                  endif
               endif
            enddo
         enddo
c        --- double the value to compute minboxdim
         minboxdim = 2.0d0*minboxdim
         call twh_minboxdim(GLB_SET,minboxdim)
      elseif ( cpnumber .eq. CP_STILLINGER_WEBER 
     &        .or. cpnumber .eq. CP_SW_PAIR_ONLY ) then
c        --- stillinger-weber: find the smallest "a" parameter
         rcut = 0.0d0
         do itype = 1,numtypes
            do jtype = itype,numtypes
               if ( twh_get_lhere(itype) .and. twh_get_lhere(jtype)
     &              ) then
                  ij = twh_getnbtype(itype,jtype)
                  swcut = twh_get_nbcoeff(6,ij)*twh_get_nbcoeff(4,ij)
                  if ( swcut .gt. rcut ) then
                     rcut = swcut
                  endif
               endif
            enddo
         enddo
c        --- set the final value of rcut
         call twh_rcut(GLB_SET,rcut)
c        --- set rcutin to rcut
         call twh_rcutin(GLB_SET,rcut)
c        --- use 2*rcut
         minboxdim = 2.0d0*rcut
         call twh_minboxdim(GLB_SET,minboxdim)
      elseif ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD .or.
     &        cpnumber .eq. CP_EAM_PAIR_ONLY ) then
c        --- embedded-atom, find the largest distance in any of
c        --- the splines
         minboxdim = 0.0d0
         do itype = 1,numtypes
            if ( twh_get_lhere(itype) ) then
               if ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD ) then
c                 --- check eam_dens
                  do jtype = 1,numtypes
                     if ( twh_get_lhere(jtype) ) then
                        call twh_eam_dens_style(GLB_GET,jtype,itype
     &                       ,eamdensstyle)
                        if ( eamdensstyle .eq. 'table' ) then
                           call twh_eam_ndens(GLB_GET,jtype,itype,ntot)
                           call twh_eam_dens(GLB_GET,1,jtype,itype,ntot
     &                          ,eamdens)
                           if ( eamdens .gt. minboxdim ) then
                              minboxdim = eamdens
                           endif
                        elseif (eamdensstyle .eq. 'Ackland cubic sum'
     &                          .or. eamdensstyle .eq. 'exponential'
     &                          ) then
                           call twh_eam_ndens(GLB_GET,jtype,itype,ntot)
                           call twh_eam_dens(GLB_GET,1,jtype,itype,ntot
     &                          ,eamdens)
                           if ( eamdens .gt .minboxdim ) then
                              minboxdim = eamdens
                           endif
                        else 
                           write(6,*) 'SETCLASSICAL:unknown'
     &                          ,' eam_dens_style'
                           write(6,*) eamdensstyle
                           lfailure = .true.
                           return
                        endif
                     endif
                  enddo
               endif
c              --- check table_pair
               do jtype = itype,numtypes
                  if ( twh_get_lhere(jtype) ) then
                     cross = twh_getnbtype(itype,jtype)
                     call twh_eam_pair_style(GLB_GET,cross,eampairstyle)
                     if ( eampairstyle .eq. 'table' .or.
     &                    eampairstyle .eq. 'Ackland 3-part' ) then
                        call twh_table_npair(GLB_GET,cross,ntot)
                        call twh_table_pair(GLB_GET,1,cross,ntot,tpair)
                        if ( tpair .gt. minboxdim ) minboxdim = tpair
                     elseif ( eampairstyle .eq. 'Belashchenko 5-8') then
c                       --- the final rho value is the limit (nb 6)
                        call twh_nbcoeff(GLB_GET,6,cross,tpair)
                        if ( tpair .gt. minboxdim ) minboxdim = tpair
                     else
                        write(6,*) 'SETCLASSICAL:unknown eam_pair_style'
                        write(6,*) eampairstyle
                        lfailure = .true.
                        return
                     endif
                  endif
               enddo
            endif
         enddo
c        --- set rcut to the maximum distance if larger than rcut
         call twh_rcut(GLB_GET,rcut)
         if ( minboxdim .gt. rcut ) then
            rcut = minboxdim
            call twh_rcut(GLB_SET,rcut)
            write(6,'(a45,a34)') 
     &           'SETCLASSICAL: adjusted rcut to be the maximum'
     &           ,' distance from tabulated potential'
            write(6,'(a17,1x,f10.5)') 'New rcut value is',rcut
         endif
         if ( .not. lcoulomb ) then
c           --- turn off the dual cutoff
            call twh_ldualcutoff(GLB_SET,.false.)
         endif
c        --- set rcutin to rcut
         call twh_rcut(GLB_GET,rcut)
         call twh_rcutin(GLB_SET,rcut)
c        --- set minboxdim to twice rcut
         minboxdim = 2.0d0*rcut
         call twh_minboxdim(GLB_SET,minboxdim)
      elseif ( cpnumber .eq. CP_SQUARE_WELL .or.
     &        cpnumber .eq. CP_REPULSIVE_WELL ) then
c        --- use twice the largest outer radius, or at least 1.0
c        --- note that nbcoeff contains the square of the diameters 
c        --- at this point as we have already called setmixterms
         minboxdim = 1.0d0
         do itype = 1,numtypes
            do jtype = itype,numtypes
               if ( twh_get_lhere(itype) .and. twh_get_lhere(jtype)
     &              ) then
                  ij = twh_getnbtype(itype,jtype)
                  if ( twh_get_nbcoeff(2,ij) .gt. minboxdim ) then
                     minboxdim = twh_get_nbcoeff(2,ij)
                  endif
               endif
            enddo
         enddo
c        --- minboxdim is twice the square root
         minboxdim = 2.0d0*dsqrt(minboxdim)
         call twh_minboxdim(GLB_SET,minboxdim)
      elseif ( cpnumber .eq. CP_TABULATED_PAIR .or.
     &        cpnumber .eq. CP_MULTIWELL .or.
     &        cpnumber .eq. CP_REPULSIVE_MULTIWELL ) then
c        --- tabulated potential find the largest distance in any of
c        --- the tables
         minboxdim = 0.0d0
         do itype = 1,numtypes
            if ( twh_get_lhere(itype) ) then
c              --- check table_pair
               do jtype = itype,numtypes
                  if ( twh_get_lhere(jtype) ) then
                     cross = twh_getnbtype(itype,jtype)
                     call twh_table_npair(GLB_GET,cross,ntot)
                     call twh_table_pair(GLB_GET,1,cross,ntot,tpair)
                     if ( tpair .gt. minboxdim ) minboxdim = tpair
                  endif
               enddo
            endif
         enddo
c        --- set rcut to the maximum distance
         rcut = minboxdim
         call twh_rcut(GLB_SET,rcut)
         call twh_rcutin(GLB_SET,rcut)
c        --- set minboxdim to twice rcut
         minboxdim = 2.0d0*rcut
         call twh_minboxdim(GLB_SET,minboxdim)
      elseif ( cpnumber .eq. CP_HARD_2580_MULTISTEP .or.
     &        cpnumber .eq. CP_REPULSIVE_2580_MULTISTEP ) then
c        --- find the largest distance of any of the wells
         minboxdim = 0.0d0
         do itype = 1,numtypes
            if ( twh_get_lhere(itype) ) then
               do jtype = itype,numtypes
                  if ( twh_get_lhere(jtype) ) then
                     cross = twh_getnbtype(itype,jtype)
                     dist = dsqrt(twh_get_nbcoeff(11,cross))
                     if ( dist .gt. minboxdim ) minboxdim = dist
                  endif
               enddo
            endif
         enddo
c        --- set rcut to the maximum distance
         rcut = minboxdim
         call twh_rcut(GLB_SET,rcut)
         call twh_rcutin(GLB_SET,rcut)
c        --- set minboxdim to twice rcut
         minboxdim = 2.0d0*rcut
         call twh_minboxdim(GLB_SET,minboxdim)
      else
c        --- use 2*rcut
         call twh_rcut(GLB_GET,rcut)
         minboxdim = 2.0d0*rcut
         call twh_minboxdim(GLB_SET,minboxdim)
      endif

c     --- for those potentials that do not specify an rcut we 
c     --- need to assign a sensible value for use with the additional
c     --- center of mass checking system
      if ( cpnumber .eq. CP_HARD_SPHERE .or.
     &     cpnumber .eq. CP_REPULSIVE_SPHERE .or.
     &     cpnumber .eq. CP_SQUARE_WELL .or.
     &     cpnumber .eq. CP_REPULSIVE_WELL ) then
c        --- use half minboxdim
         call twh_minboxdim(GLB_GET,minboxdim)
         rcut = 0.5d0*minboxdim
         call twh_rcut(GLB_SET,rcut)
c        --- set rcutin to the same value
         call twh_rcutin(GLB_SET,rcut)
      endif

c     --- compute and store the squares of the cutoffs
      call twh_rcut(GLB_GET,rcut)
      rcutsq = rcut*rcut
      call twh_rcutsq(GLB_SET,rcutsq)
      call twh_rcutin(GLB_GET,rcutin)
      rcutinsq = rcutin*rcutin
      call twh_rcutinsq(GLB_SET,rcutinsq)

c     --- set the pressurestyle
      if ( cpnumber .eq. CP_LENNARD_JONES
     &     .or. cpnumber .eq. CP_SCALED_LENNARD_JONES
     &     .or. cpnumber .eq. CP_12_6_PLUS_12_10_H_BOND
     &     .or. cpnumber .eq. CP_12_6_PLUS_SOLVATION
     &     .or. cpnumber .eq. CP_UFF_12_6
     &     .or. cpnumber .eq. CP_9_6
     &     .or. cpnumber .eq. CP_12_9_6
     &     .or. cpnumber .eq. CP_EAM_PAIR_ONLY
     &     .or. cpnumber .eq. CP_EXPONENTIAL_6
     &     .or. cpnumber .eq. CP_EXPONENTIAL_12_6
     &     .or. cpnumber .eq. CP_GORDON_N_6
     &     ) then
         call twh_pressurestyle(GLB_SET,PRES_VIRIAL)
      elseif ( cpnumber .eq. CP_SQUARE_WELL
     &        .or. cpnumber .eq. CP_HARD_SPHERE
     &        .or. cpnumber .eq. CP_HARD_2580_MULTISTEP
     &        ) then
         call twh_pressurestyle(GLB_SET,PRES_RADIAL)
      else
         call twh_pressurestyle(GLB_SET,PRES_THERMO)
      endif

      return
      end

