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
      subroutine twh_qqassign(lfailure,imolty,forcefield)
c     ******************************************************************
c     * makes the charge assignment for a molecule we just assembled   *
c     *                                                                *
c     * originally written 12-10-2004 by M.G. Martin                   *
c     * last modified 08-05-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_CBNAME
#define FUNCTION_GET_IJVIB
#define FUNCTION_GET_INVIB
#define FUNCTION_GET_NBNAME
#define FUNCTION_GET_QBASEVALUE
#define FUNCTION_GET_QBIFF
#define FUNCTION_GET_QBINAMES
#define FUNCTION_GET_QBIVALUE
#include "functions.h"
c     --- variables passed to/from the subroutine
      character*(*) forcefield
      logical lfailure
      integer imolty
c     --- local variables
c     --- character string scalars
      character*(FFNAMELEN) viborder,testorder,wildorder
      character*30 charge_assignment
c     --- character string arrays
      character*(FFNAMELEN) names
      dimension names(2)
c     --- integer scalars
      integer param,type,ivib,itype,jtype,iunit,junit,nunit
c     --- double precision scalars
      double precision qadjust,totalcharge,charge

c      write(6,*) 'start QQASSIGN'

      call twh_nunit(GLB_GET,imolty,nunit)

      call twh_charge_assignment(GLB_GET,imolty,charge_assignment)
      if ( charge_assignment .eq. 'bond increment' ) then
c        --- sort through all of the atoms in the molecule and 
c        --- determine the charges using bond increment method
         do iunit = 1,nunit
            call twh_ntype(GLB_GET,imolty,iunit,itype)
c           --- assign the base charge
            charge = twh_get_qbasevalue(itype)
c           --- compute the bond increments
            names(1) = twh_get_cbname(itype)
            do ivib = 1,twh_get_invib(imolty,iunit)
               junit = twh_get_ijvib(imolty,iunit,ivib)
               call twh_ntype(GLB_GET,imolty,junit,jtype)
               names(2) = twh_get_cbname(jtype)
c              --- get the vibration order
               call twh_vibration_order(GLB_GET,iunit,ivib,viborder)
c              --- sort through the list and look for matches
               param = 0
               qadjust = 0.0d0
c              --- on the first pass check for exact bond-increment type
               type = 0
               do while ( param .eq. 0 .and. type .lt. TBIMAX )
                  type = type + 1
                  call twh_biorder(GLB_GET,type,testorder)
                  if ( twh_get_qbiff(type) .eq. forcefield .and.
     &                 viborder .eq. testorder ) then
                     if ( twh_get_qbinames(type,1) .eq. names(1) 
     &                    .and. 
     &                    twh_get_qbinames(type,2) .eq. names(2) ) then
                        param = type
                        qadjust = twh_get_qbivalue(param)
                     endif
                     if ( twh_get_qbinames(type,1) .eq. names(2) 
     &                    .and. 
     &                    twh_get_qbinames(type,2) .eq. names(1) ) then
                        param = type
c                       --- move charge in opposite direction
                        qadjust = -twh_get_qbivalue(param)
                     endif
                  endif
               enddo
c              --- on the second pass check for wild bond-increment type
               wildorder = 'wild'
               type = 0
               do while ( param .eq. 0 .and. type .lt. TBIMAX )
                  type = type + 1
                  call twh_biorder(GLB_GET,type,testorder)
                  if ( twh_get_qbiff(type) .eq. forcefield .and.
     &                 wildorder .eq. testorder ) then
                     if ( twh_get_qbinames(type,1) .eq. names(1) 
     &                    .and. 
     &                    twh_get_qbinames(type,2) .eq. names(2) ) then
                        param = type
                        qadjust = twh_get_qbivalue(param)
                     endif
                     if ( twh_get_qbinames(type,1) .eq. names(2) 
     &                    .and. 
     &                    twh_get_qbinames(type,2) .eq. names(1) ) then
                        param = type
c                       --- move charge in opposite direction
                        qadjust = -twh_get_qbivalue(param)
                     endif
                  endif
               enddo
c              --- search for global wildcard matches
               type = 0
               do while ( param .eq. 0 .and. type .lt. TBIMAX )
                  type = type + 1
                  if ( twh_get_qbiff(type) .eq. forcefield ) then
                     if ( twh_get_qbinames(type,1) .eq. 'wild' 
     &                    .and. 
     &                    twh_get_qbinames(type,2) .eq. 'wild' ) then
                        param = type
                        qadjust = twh_get_qbivalue(param)
                     endif
                  endif
               enddo
               if ( param .eq. 0 ) then
                  write(6,*) 'parameters not found for bond increment'
                  write(6,*) 'atoms ',iunit,' and ',junit,' vib ',ivib
                  write(6,*) 'names [',names(1),'][',names(2),']'
                  lfailure = .true.
                  return
               endif
               charge = charge + qadjust
            enddo
c           --- store the final charge for this atom
            call twh_qqatom(GLB_SET,imolty,iunit,charge)
         enddo
c        --- output the assigned charges
         write(6,'(a35,1x,i4)') 'Charges assigned for Molecule Type:'
     &        ,imolty
         totalcharge = 0.0d0
         do iunit = 1,nunit
            call twh_ntype(GLB_GET,imolty,iunit,itype)
            call twh_qqatom(GLB_GET,imolty,iunit,charge)
            write(6,'(a5,i5,1x,a7,1x,a10,1x,a7,1x,f10.5)') 
     &           'Unit:',iunit,'nbname:',twh_get_nbname(itype)
     &           ,'Charge:',charge
            totalcharge = totalcharge + charge
         enddo
         write(6,'(a,1x,i4,1x,a3,1x,f10.5)') 
     &        'Total charge for Molecule Type:',imolty,'is:'
     &        ,totalcharge
      elseif ( charge_assignment .eq. 'none' ) then
c        --- zero charge on every interaction site
         do iunit = 1,nunit
            charge = 0.0d0
            call twh_qqatom(GLB_SET,imolty,iunit,charge)
         enddo
      else
         write(6,*) 'QQASSIGN: unknown charge_assignment'
         write(6,*) charge_assignment
         lfailure = .true.
         return
      endif
c      write(6,*) 'finished QQASSIGN'

      return
      end
