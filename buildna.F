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
      subroutine twh_buildna(lfailure,imolty)
c     ******************************************************************
c     * builds the input file for a nucleic acid sequence starting     *
c     * from the 5'-terminus                                           *
c     * imolty: molecule type number                                   *
c     *                                                                *
c     * originally written 04-10-2003 by C.D. Lorenz                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"
c     ---- variables passed to/from the subroutine
      logical lfailure
      integer imolty
c     --- local variables
      integer hterm,nunit
      logical lautotorsion

      character*2 pepname
      dimension pepname(NUMAX)
      character*3 shortname
      character*10 atomname
      dimension atomname(NUMAX)
      character*10 forcefield
      character*15 idname
      character*20 termflag
      character*40 aminoname

      integer imonomer,iatom,iprev,total,totatoms,itemp,ibondend

c      write(6,*) 'start BUILDNA'
      idname = 'buildna'
      write(6,*) 'Building the nucleic acid input file'
c     --- set lautotorsion
      lautotorsion = .true.
      
      hterm = twh_read_labeled_integer(lfailure,4,'terminus',
     &    .false.,idname)
      if ( lfailure ) return

c     --- error check hterm
      if ( hterm .eq. 3 ) then
         write(6,*) '3 prime end is hydrogen terminated'
         write(6,*) '5 prime end is phosphorated'
      elseif ( hterm .eq. 5 ) then
         write(6,*) '3 prime end is phosphorated'
         write(6,*) '5 prime end is hydrogen terminated'
      else
         write(6,*) 'BUILDNA: hterm must be 3 or 5'
         lfailure = .true.
         return
      endif

      call twh_read_labeled_string(lfailure,4,'forcefield',.false.
     &     ,idname,forcefield)
      if ( lfailure ) return

c     --- check for valid forcefield
      if ( forcefield .eq. 'Charmm27' ) then
         write(6,*) '   using the Charmm27 forcefield'
      else
         write(6,*) 'Invalid forcefield of ',forcefield
         lfailure = .true.
         return
      endif

      call twh_read_label(lfailure,4,'monomername',.false.,idname)
      if ( lfailure ) return
c     --- read in the nucleic acid sequence
      call twh_nunit(GLB_GET,imolty,nunit)
      do imonomer = 1,nunit
         call twh_read_string(lfailure,4, 'monomername', .false.
     &        ,pepname(imonomer))
         if ( lfailure ) return
         write(6,15) 'nucleic acid:',imonomer,'name:',pepname(imonomer)
 15      format(a13,i5,1x,a5,a2)
      enddo

c     --- initialize iprev and total
      iprev = 1
      total = 0
      totatoms = 0
c     --- initialize ibondend
      ibondend = 0

      do imonomer = 1,nunit
c        --- increment the current atom count by the previous atom count
c        --- plus the number of atoms in the previous residue
         iatom = iprev + total
         if ( imonomer .eq. 1 ) then
c           --- build the 5' terminus
            termflag = '5 prime'
         elseif ( imonomer .eq. nunit ) then
c           --- build the 3' terminus
            termflag = '3 prime'
         else
            termflag = 'neither'
         endif
         if ( forcefield .eq. 'Charmm27' ) then
            call twh_polycharmm27na(lfailure,imolty,iatom,termflag
     &           ,total,ibondend,imonomer,pepname,atomname,aminoname
     &           ,shortname,hterm)
            if ( lfailure ) return
         else
            write(6,*) 'BUILDNA: unknown forcefield'
            lfailure = .true.
            return
         endif

         totatoms = totatoms + total

         if ( totatoms .gt. NUMAX ) then
            write(6,*) 'total, totatoms',total,totatoms
            write(6,*) 'totatoms > NUMAX'
            lfailure = .true.
            return
         endif

c        --- assign the shortcut names to the atoms
         do itemp = iatom,iatom+total-1
            call twh_aminoshort(GLB_SET,imolty,itemp,shortname)
            call twh_aminonum(GLB_SET,imolty,itemp,imonomer)
         enddo

c        --- update iprev
         iprev = iatom
         if ( iprev +35 .gt. NUMAX ) then
            write(6,*) 'WARNING: approaching numax in buildna'
            write(6,*) 'total, totatoms',total,totatoms
            write(6,*) 'if the code suddenly crashes then you exceeded'
     &           ,' numax'
         endif
      enddo

      call twh_assemble(lfailure,imolty,totatoms,0,forcefield,atomname
     &     ,'standard','manual',lautotorsion)
      if ( lfailure ) return
c      write(6,*) 'finished BUILDNA'

      return
      end
