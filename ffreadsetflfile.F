#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2003-2011 Marcus G. Martin                       *
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
      subroutine twh_ffreadsetflfile(lfailure,nbtype,ffname)
c     ******************************************************************
c     * reads in a setfl force field file and transforms it into       *
c     * corresponding Towhee data structures ao that it can create     *
c     * a Towhee force field file                                      *
c     *                                                                *
c     * originally written 03-27-2003 by M.G. Martin                   *
c     * last modifed 08-05-2011 by M.G. Martin                         *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer nbtype
      character*(FFNAMELEN) ffname
c     --- local variables
      character*2 ename
      character*10 eleflag
      character*30 filename
      integer nrho,nr,itype,ielement,jtype,cross,irho,ir,ipair
c     --- double precision scalars
      double precision drho,dr,rcutall,value,dtemp,dvalue,rcut
c     --- double precision arrays
      double precision tableread
      dimension tableread(MAXTABLE)

      write(6,*) 'Please enter name of the setfl file'
      write(6,*) ' (max. 30 characters) or type exit if you have'
      write(6,*) ' entered this choice by mistake'
      read(5,*) filename

      if ( filename .eq. 'exit' ) then
         lfailure = .true.
         return
      endif


      write(6,*) 'Enter a name for this forcefield (10 characters max)'
      read(5,*) ffname

c     --- this is EAM
      call twh_classical_potential(GLB_SET,'Embedded Atom Method')
      call twh_classical_mixrule(GLB_SET,'Explicit')

      call twh_openfile(lfailure,10,filename,'DEFAULT','OLD','DEFAULT')
      if ( lfailure ) return
c     --- skip the first 3 lines
      read(10,*)
      read(10,*)
      read(10,*)
      read(10,*) nbtype

      if (nbtype .gt. MAXTABTYPE) then
         write(6,*)'error: number of types greater than MAXTABTYPE'
         lfailure = .true.
         return
      endif
c     --- set nbff
      do itype = 1,nbtype
         call twh_nbff(GLB_SET,itype,ffname)
      enddo
c     --- set the styles all to table
      do itype = 1,nbtype
         call twh_eam_dens_style(GLB_SET,itype,itype,'table')
         call twh_eam_embed_style(GLB_SET,itype,'table')
         do jtype = itype,nbtype
            cross = twh_getnbtype(itype,jtype)
            call twh_eam_pair_style(GLB_SET,cross,'table')
         enddo
      enddo

      read(10,9901) nrho, drho, nr, dr, rcutall
c     --- safety check nrho and nr
      if ( nrho .gt. MAXTABLE ) then
         write(6,*) 'nrho exceeds MAXTABLE'
         lfailure = .true.
         return
      endif
      if ( nr .gt. MAXTABLE ) then
         write(6,*) 'nr exceeds MAXTABLE'
         lfailure = .true.
         return
      endif

c     --- translate to Towhee variable
      do itype = 1,nbtype
c        --- set the size of the arrays
         call twh_eam_nembed(GLB_SET,itype,nrho)
         call twh_eam_ndens(GLB_SET,itype,itype,nr)
c        --- set the distances in eam_embed
         do irho = 1,nrho
            value = dble(irho-1)*drho
            call twh_eam_embed(GLB_SET,1,itype,irho,value)
         enddo
c        --- set the distances in eam_dens
         do ir = 1,nr
            value = dble(ir-1)*dr
            call twh_eam_dens(GLB_SET,1,itype,itype,ir,value)
         enddo
         do jtype = 1,nbtype
            cross = twh_getnbtype(itype,jtype)
            call twh_table_npair(GLB_SET,cross,nr)
            do ir = 1,nr
               value = dble(ir-1)*dr
               call twh_table_pair(GLB_SET,1,cross,ir,value)
            enddo
         enddo
      enddo

      rcut = dsqrt(rcutall)
      call twh_rcut(GLB_SET,rcut)
      eleflag = 'getelement'
      do itype = 1,nbtype
         read(10,30) ielement,dtemp
         call twh_mass(GLB_SET,itype,dtemp)
         call twh_getelement(lfailure,eleflag,dtemp,ename,ielement)
         if ( lfailure ) return
         read(10,9902) (tableread(irho),irho=1,nrho)
         do irho = 1,nrho
            call twh_eam_embed(GLB_SET,2,itype,irho,tableread(irho))
         enddo
c        --- convert the units on the embedding function
         do irho = 1,nrho
            call twh_eam_embed(GLB_GET,2,itype,irho,dvalue)
            dvalue = CON_EVTOK*dvalue
            call twh_eam_embed(GLB_SET,2,itype,irho,dvalue)
         enddo
         read(10,9902) (tableread(ir),ir=1,nr)
         do ir = 1,nr
            call twh_eam_dens(GLB_SET,2,itype,itype,ir,tableread(ir))
         enddo
      enddo
      do itype = 1,nbtype
         do jtype= 1,itype
            cross = twh_getnbtype(itype,jtype)
            read(10,9902) (tableread(ir),ir=1,nr)
c           --- need to convert this by dividing out the distance
c           --- this is not valid for the first point, which we are 
c           --- going to remove
            do ir = 2,nr
c              --- set the second index of table_pair
               call twh_table_pair(GLB_GET,1,cross,ir,dtemp)
               dvalue = CON_EVTOK*tableread(ir) / dtemp
               call twh_table_pair(GLB_SET,2,cross,ir-1,dvalue)
c              --- shift the value down one slot for the first index
               call twh_table_pair(GLB_SET,1,cross,ir-1,dtemp)
               if ( itype .eq. 1 .and. jtype .eq. 1 ) then
                  write(12,*) dtemp,dvalue/CON_EVTOK
               elseif ( itype .eq. 2 .and. jtype .eq. 1 ) then
                  write(13,*) dtemp,dvalue/CON_EVTOK
               else
                  write(14,*) dtemp,dvalue/CON_EVTOK
               endif
            enddo
c           --- reduce table_npair by 1
            call twh_table_npair(GLB_GET,cross,ipair)
            ipair = ipair - 1
            call twh_table_npair(GLB_SET,cross,ipair)
         enddo
      enddo

      close(10)

c     --- set all of the cross term densities to the self densities
      do itype = 1,nbtype
         call setcrossdens(itype)
      enddo

 30   format( i5, g15.5 )
 9901 format( i5, e24.16, i5, 2 e24.16 )
 9902 format( 5 e24.16 )

      return
      end

