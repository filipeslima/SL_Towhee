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
      subroutine twh_ffpmf(lfailure,ntypes)
c     ******************************************************************
c     * reads in a series of potentials of mean force and turns them   *
c     * into a forcefield                                              *
c     *                                                                *
c     * originally written 09-23-2004 by M.G. Martin                   *
c     * last modifed 08-04-2011 by M.G. Martin                         *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes
c     --- local variables
      character*2 ename
      character*(FFNAMELEN) ffname,ctemp
      character*15 idname
      character*30 filename
      integer iversion,nversion,ntij,itype,itable,jtype,ipair
      double precision tpairone,tpairtwo
      double precision dtemp

      idname = 'ffpmf'
c     --- currently version 1
      nversion = 1
c     --- tabulated potential
      call twh_classical_potential(GLB_SET,'Tabulated Pair')
c     --- explicit mixrule
      call twh_classical_mixrule(GLB_SET,'Explicit')

      write(6,*) 'Enter a name for this forcefield (10 characters max)'
      read(5,*) ffname

c     --- get ntypes
      write(6,*) 'Enter the number of atoms types in this forcefield'
      read(5,*) ntypes
      if ( ntypes .gt. MAXTABTYPE ) then
         write(6,*)'error: number of types greater than MAXTABTYPE'
         lfailure = .true.
         return
      endif

c     --- set nbff
      do itype = 1,ntypes
         write(6,*) 'Enter a nonbond name for atom type:',itype
         read(5,*) ctemp
         call twh_nbname(GLB_SET,itype,ctemp)
         write(6,*) 'Enter a bonding name for atom type:',itype
         read(5,*) ctemp
         call twh_cbname(GLB_SET,itype,ctemp)
         write(6,*) 'Enter a bend angle name for atom type:',itype
         read(5,*) ctemp
         call twh_caname(GLB_SET,itype,ctemp)
         write(6,*) 'Enter a torsion name for atom type:',itype
         read(5,*) ctemp
         call twh_ctname(GLB_SET,itype,ctemp)
         write(6,*) 'Enter a mass for atom type:',itype
         read(5,*) dtemp
         call twh_mass(GLB_SET,itype,dtemp)
         write(6,*) 'Enter 2 letter element code for atom type:',itype
         read(5,*) ename
         call twh_elementname(GLB_SET,itype,ename)
c        --- run through all cross terms (explicit)
         do jtype = itype,ntypes
            call twh_nbff(GLB_SET,itype,ffname)
            ntij = twh_getnbtype(ntypes,ntypes)
c           --- get the filename
            write(6,*) 'Please enter the PMF file for the potential'
     &           ,' of type:',itype,' with type:',jtype
            read(5,*) filename
            call twh_openfile(lfailure,10,filename,'DEFAULT','OLD'
     &           ,'DEFAULT')
            if ( lfailure ) return

c           --- read in the pmf version
            iversion = twh_read_labeled_integer(lfailure,10
     &           ,'pmf version',.false., idname)
            if ( lfailure ) return
            if ( iversion .ne. nversion ) then
               write(6,*) 'invalid towhee pmf version',iversion
               lfailure = .true.
               return
            endif

c           --- read in the number of datapoints
            ipair = twh_read_labeled_integer(lfailure,10,
     &           'datapoints', .false., idname)
            if ( lfailure ) return
            call twh_table_npair(GLB_SET,ntij,ipair)
            if ( ipair .lt. 1 ) then
               write(6,*) 'datapoints must be positive',ipair
               lfailure = .true.
               return
            endif
            if ( ipair .gt. MAXTABLE ) then
               write(6,*) 'FFPMF: datapoints exceed MAXTABLE'
               write(6,*) 'datapoints',ipair,'MAXTABLE'
     &              ,MAXTABLE
               lfailure = .true.
               return
            endif

c           --- read in the data
            call twh_read_label(lfailure,10, 'data', .false., idname)
            if ( lfailure ) return
            do itable = 1,ipair
               read(10,*) tpairone,tpairtwo
               call twh_table_pair(GLB_SET,1,ntij,itable,tpairone)
               call twh_table_pair(GLB_SET,2,ntij,itable,tpairtwo)
            enddo
            close(10)
            rewind(10)
         enddo
      enddo
      return
      end

