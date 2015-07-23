#include "preproc.h"
#define SOLVTYPE 17
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
      function twh_veefone(ichain,ibox,cordtype)
c     ******************************************************************
c     * Computes the solvation energy using the effective energy       *
c     * function (EEF1) as descriped in Lazaridis and Karplus,         *
c     * Proteins 35, 133-152 (1999)                                    *
c     *                                                                *
c     * originally written 11-04-2002 by M.G. Martin and C.D. Lorenz   *
c     * last update 08-09-2011 by M.G. Martin                          *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_IJVIB
#define FUNCTION_GET_INVIB
#define FUNCTION_GET_NBCOEFF
#define FUNCTION_VEEFONE
#include "functions.h"
#define DEBUG_VEEFONE 0

c     --- variables passed to/from the function
      integer ichain,cordtype,ibox

c     --- local variables
      integer j,eeftypej,iunit,typei,iexclude,no_exclude,imolty
      integer jexclude,nni_exclude,nnj_exclude,findexclude
      integer iitype,jjtype,nunit

      integer eef_exclude
      dimension eef_exclude(20)

      double precision veef1i,solvent_cutoff,dgisum
      double precision piterm,eefxi,eefyi,eefzi,eefdgi_ref
      double precision eefdgi_free,eefvi,eeflambdai,eefri,eefsum
      double precision eefxj,eefyj,eefzj,eefdx,eefdy,eefdz,eefrij
      double precision eefvj,eeflambdaj,eefrj,eefxij,eefxji,totsum
      double precision eefdgj_ref,eefdgj_free,eefterm1,eefterm2,onepi

c     --- this value from the cutoff was taken from the charmm source
c     --- code 9
      solvent_cutoff = 9.0d0
      totsum = 0.0d0
      twh_veefone = 0.0d0
      dgisum = 0.0d0
#if DEBUG_VEEFONE
      write(6,*) "starting twh_veefone"
#endif

      call twh_constant_pi(GLB_GET,onepi)
      call twh_moltyp(GLB_GET,ichain,imolty)
      call twh_nunit(GLB_GET,imolty,nunit)

      piterm = 4.0d0 * onepi * dsqrt(onepi)
      do iunit = 1,nunit - 1
         call twh_ntype(GLB_GET,imolty,iunit,typei)
c        --- determine the coordinates of atom i
         call twh_coordinates(GLB_GET,cordtype,ichain,iunit
     &        ,eefxi,eefyi,eefzi)
c        --- make exclusion list because nearest and 
c        --- next-nearest neighbor lists are excluded
c        --- from the solvation free energy calcs
         no_exclude = 1
         eef_exclude(no_exclude) = iunit
         do iexclude = 1, twh_get_invib(imolty,iunit)
            no_exclude = no_exclude + 1
            nni_exclude = twh_get_ijvib(imolty,iunit,iexclude)
            eef_exclude(no_exclude) = nni_exclude
            do jexclude = 1,twh_get_invib(imolty,nni_exclude)
               nnj_exclude = twh_get_ijvib(imolty,nni_exclude,jexclude)
               if (nnj_exclude .gt. iunit) then
                  no_exclude = no_exclude + 1
                  eef_exclude(no_exclude) = nnj_exclude
               endif
            enddo
         enddo
c        --- determine the van der Waals radius (defined as 1/2 the 
c        --- distance to the energy minimum in the LJ potential)
c        --- CHECK: nbcoeff(1,type) = sigma
         
         eefsum = 0.0d0
         iitype = twh_getnbtype(typei,typei)
         eefdgi_ref = twh_get_nbcoeff(6,iitype)

         if (typei .le. SOLVTYPE ) then         
            do j = iunit+1, nunit
c           --- determine the coordinates of atom j
               findexclude = 0
               do jexclude = 1,no_exclude
                  if (findexclude .eq. 0) then
                     if (j .eq. eef_exclude(jexclude)) then
                        findexclude = 1
                     endif
                  endif
               enddo
               eefterm1 = 0.0d0
               eefterm2 = 0.0d0
               call twh_ntype(GLB_GET,imolty,j,eeftypej)
               if (eeftypej .le. SOLVTYPE ) then
                  if (findexclude .eq. 0) then
                     call twh_coordinates(GLB_GET,cordtype,ichain,j
     &                    ,eefxj,eefyj,eefzj)
c                    --- determine the distance between atoms i and j
                     call twh_mimage(ibox,eefxi,eefyi,eefzi,eefxj,eefyj
     &                    ,eefzj,eefdx,eefdy,eefdz)
                     eefrij = twh_distance(eefdx,eefdy,eefdz)
                     if (eefrij .lt. solvent_cutoff) then
                        jjtype = twh_getnbtype(eeftypej,eeftypej)
c                       --- determine solvation free energies for atom i
                        eefdgi_free = twh_get_nbcoeff(7,iitype)
c                       --- look up the volume of atom i
                        eefvi = twh_get_nbcoeff(5,iitype)
c                       --- determine the van der Waals radius (defined
c                       --- as 1/2 the distance to the energy minimum
c                       --- in the LJ potential)
c                       --- CHECK: twh_get_nbcoeff(1,type) = sigma
                        eefri = twh_get_nbcoeff(11,iitype)
c                       --- determine the correlation length
                        eeflambdai = twh_get_nbcoeff(10,iitype)
                        if (eeflambdai .eq. 0.0d0) then
                           eeflambdai = 3.5d0
                        endif

c                       --- look up the volume of atom j
                        eefvj = twh_get_nbcoeff(5,jjtype)
c                       --- determine the correlation length
                        eeflambdaj = twh_get_nbcoeff(10,jjtype) 
                        if (eeflambdaj .eq. 0.0d0) then
                           eeflambdaj = 3.5d0
                        endif
c                       --- determine the van der Waals radius (defined
c                       --- as 1/2 the distance to the energy minimum
c                       --- in the LJ potential)
c                       --- CHECK: twh_get_nbcoeff(1,type) = sigma
                        eefrj = twh_get_nbcoeff(11,jjtype)
c                       --- determine the exponential terms: xij and xji
                        eefxij = (eefrij - eefri)/eeflambdai
                        eefxji = (eefrij - eefrj)/eeflambdaj
c                       --- determine solvation free energies for atom j
                        eefdgj_ref = twh_get_nbcoeff(6,jjtype)
                        eefdgj_free = twh_get_nbcoeff(7,jjtype)
                        eefterm1 = ((2.0d0*eefdgi_free)
     &                       /(piterm*eeflambdai*eefrij*eefrij))
     &                       *exp(-1.0d0*eefxij*eefxij)*eefvj
                        eefterm2 = ((2.0d0*eefdgj_free)
     &                       /(piterm*eeflambdaj*eefrij*eefrij))
     &                       *exp(-1.0d0*eefxji*eefxji)*eefvi
                        eefsum = eefsum + (eefterm1 + eefterm2)
                     endif
                  endif
               endif
            enddo
         endif
c        --- calculate delta G_i
         totsum = totsum + eefsum
         dgisum = dgisum + eefdgi_ref
         veef1i = eefdgi_ref - eefsum   
         twh_veefone = twh_veefone + veef1i
#if DEBUG_VEEFONE
         write(6,*) "i,typei,eefvi,eefdgi_ref,eefdgi_free,lambdai",
     &          iunit,typei,eefvi,eefdgi_ref,eefdgi_free,eeflambdai 
         write(6,*) "totsum, dgisum, veef1i, twh_veefone",
     &          totsum, dgisum, veef1i, twh_veefone
#endif
      enddo
      twh_veefone = twh_veefone * CON_KCALTOK

#if DEBUG_VEEFONE
      write(6,*) "twh_veefone returning with value ", twh_veefone
#endif
      return
      end


      subroutine twh_thermocorrect(tempi,tempref)
c     ******************************************************************
c     * Computes temperature corrected values of the solvation         *
c     * free energy for use in EEF1                                    *
c     *                                                                *
c     * originally written 11-04-2002 by M.G. Martin and C.D. Lorenz   *
c     * last update 10-13-2008 by M.G. Martin                          *
c     ******************************************************************
      implicit  none
#define FUNCTION_GET_NBCOEFF
#include "functions.h"
c     --- variables passed to/from the function
      double precision tempi,tempref
c     --- local variables
      integer atype
      double precision ds_ref,dtemp

      do atype = 1, SOLVTYPE
         ds_ref = (twh_get_nbcoeff(8,atype)
     &        - twh_get_nbcoeff(6,atype))/tempref
      
         dtemp = twh_get_nbcoeff(6,atype) - 
     &        (ds_ref * (tempi - tempref)) - 
     &        ((twh_get_nbcoeff(9,atype)/1000.0d0) * tempi * 
     &        log(tempi/tempref)) + 
     &        ((twh_get_nbcoeff(9,atype)/1000.0d0) * 
     &        (tempi - tempref))
         call twh_nbcoeff(GLB_SET,6,atype,dtemp)

         dtemp = twh_get_nbcoeff(7,atype) - 
     &        (ds_ref * (tempi - tempref)) - 
     &        ((twh_get_nbcoeff(9,atype)/1000.0d0) * tempi * 
     &        log(tempi/tempref)) + 
     &        ((twh_get_nbcoeff(9,atype)/1000.0d0) *
     &        (tempi - tempref))
         call twh_nbcoeff(GLB_SET,7,atype,dtemp)
      enddo

      return
      end
