#include "preproc.h"
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
      subroutine twh_readlammps(lfailure)
c     ******************************************************************
c     * reads in the lammps input parameters and uses those to set up  *
c     * the towhee_altinp and towhee_coords files.  Obviously, there   *
c     * Lammps and Towhee variables which do not overlap, and for those*
c     * a default value is put in place.                               *
c     * this subroutine is designed as a conversion routine to turn    *
c     * Lammps runs into a form suitable for a Towhee run              *
c     * thanks to Steve Plimpton for helpful discussions about how to  *
c     * set up this file, and for allowing me to take some code from   *
c     * Lammps 2001                                                    *
c     *                                                                *
c     * originally written 02-27-2002 by M.G. Martin                   *
c     * last modified 01-31-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_MASS
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- function declarations local to the readlammps subroutine
      logical match
c     --- local variables
      character*2 ename
      character*10 eleflag
      character*80 str
      logical ldone
      integer m,itype,massflag,natomflag,nbondflag,nangleflag
      integer ndihedflag,nimproflag,noncoeffflag,nbondcoeffflag
      integer nanglecoeffflag,ndihedcoeffflag,nimprocoeffflag
      integer nbondbondflag,nbondangleflag,nmidbondtorsionflag
      integer nendbondtorsionflag,nangletorsionflag
      integer nangleangletorsionflag,nbondbond13flag
      integer nangleangleflag,amberflag
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nimprotypes
     &     ,naatypes,nbitypes,noftypes,nhbondtypes
      integer natoms,nbonds,nangles,ndihedrals,nimpropers,nangleangle
      integer nonstyle,bondstyle,lmpangstyle,dihedstyle,improstyle
      integer ielement,numboxes,total
c     --- double precision scalars
      double precision localmass,dvalue
      double precision box
      dimension box(2,3)
      double precision special
      dimension special(3)
      character*30 filename

c     --- initialize arrays
      call twh_initialize(lfailure)
      if ( lfailure ) return
c     --- set to a single box
      numboxes = 1
c     --- allocate and initialize arrays for numboxes
      call twh_allocate_numboxes(numboxes)
c     --- allocate and initialize arrays for numboxes and MAXMOLTY
      call twh_allocate_maxmolty_numboxes(NTMAX,numboxes)
c     --- allocate for torsions
      call twh_allocate_maxmolty_maxunit_maxnumtor(NTMAX,NUMAX
     &     ,MAXTOR)

c     --- set pmrotate to 1.0
      call twh_pmrotate(GLB_SET,1.0d0)

c     --- initialize special
      special(1) = 0.0d0
      special(2) = 0.0d0
      special(3) = 0.0d0
      call twh_loutlammps(GLB_SET,.true.)

c     --- first scan through the lammps input file to pick out some
c     --- useful information about the force field styles, the cutoffs
c     --- temperature, pressure
      call readlammps_input(lfailure,nonstyle,bondstyle,lmpangstyle
     &     ,dihedstyle,improstyle,amberflag,special)
      if ( lfailure ) return

      write (6,*) 'Reading lammps_data file ...'

      ntypes = 0
      nbondtypes = 0
      nangletypes = 0
      ndihedtypes = 0
      nimprotypes = 0
      naatypes = 0
      nbitypes = 0
      noftypes = 0
      nhbondtypes = 0
      nangleangle = 0
      
c     --- read data file header - fixed format
      filename = 'lammps_data'
      call twh_openfile(lfailure,2,filename,'DEFAULT','OLD','DEFAULT')
      if ( lfailure ) return
      read (2,*)
      read (2,*)
c     --- read the number of atoms
      read (2,*) natoms
      if (natoms .le. 0) stop 'Atoms .le. 0'
c     --- allocate memory for temporary coordinate storage 
      total = 3*natoms
      dvalue = 0.0d0
      call twh_coordtemp(GLB_ALLOC,total,dvalue)
c     --- read the number of bonds
      read (2,*) nbonds
      if (nbonds .lt. 0) stop 'Bonds .lt. 0'
c     --- read the number of angles
      read (2,*) nangles
      if (nangles.lt.0) stop 'Angles .lt. 0'
c     --- read the number of dihedrals
      read (2,*) ndihedrals
      if (ndihedrals.lt.0) stop 'Dihedrals .lt. 0'
c     --- read the number of impropers
      read (2,*) nimpropers
      if (nimpropers.lt.0) stop 'Impropers .lt. 0'
c     --- spacer
      read (2,*)
c     --- read the number of atom types
      read (2,*) ntypes
      if (ntypes.le.0) stop 'Atom types .le. 0'
      if ( ntypes .gt. NNTYPE ) stop 'Atoms types .gt. NNTYPE'
c     --- read the number of bond types
      if (nbonds .gt. 0) then
         read (2,*) nbondtypes
         if ( nbondtypes .le. 0) stop 'Bond types .le. 0'
         if ( nbondtypes .gt. TVIBMAX ) stop 'Bond types .gt. TVIBMAX'
      endif
c     --- read the number of angle types
      if (nangles.gt.0) then
         read (2,*) nangletypes
         if ( nangletypes .le. 0) stop 'Angle types .le. 0'
         if ( nangletypes .gt. TBENMAX ) stop 'Angle types .gt. TBENMAX'
      endif
c     --- read the number of dihedral types
      if (ndihedrals.gt.0) then
         read (2,*) ndihedtypes
         if ( ndihedtypes .le. 0 ) stop 'Dihedral types .le. 0'
         if ( ndihedtypes .gt. TTORMAX )
     &        stop 'Dihedral types .gt. TTORMAX'
      endif
c     --- read in the number of improper types
      if (nimpropers.gt.0) then
         read (2,*) nimprotypes
         if ( nimprotypes .le. 0 ) stop 'Improper types .le. 0'
         if ( nimprotypes .gt. TTORMAX ) 
     &        stop 'Improper types .gt. TTORMAX'
      endif
c     --- spacer
      read (2,*)
c     --- read in the box dimensions
      read (2,*) box(1,1),box(2,1)
      read (2,*) box(1,2),box(2,2)
      read (2,*) box(1,3),box(2,3)
      if (box(1,1) .ge. box(2,1) .or. box(1,2) .ge. box(2,2) .or.
     &     box(1,3) .ge. box(2,3) ) then
         stop 'Illegal box dimensions'
      endif
c     --- initialize hmatrix
      call twh_hmatrix(GLB_INIT,1,3,3,0.0d0)
c     --- set hmatrix
      dvalue = box(2,1)-box(1,1)
      call twh_hmatrix(GLB_SET,1,1,1,dvalue)
      dvalue = box(2,2)-box(1,2)
      call twh_hmatrix(GLB_SET,1,2,2,dvalue)
      dvalue = box(2,3)-box(1,3)
      call twh_hmatrix(GLB_SET,1,3,3,dvalue)

c     --- read remainder of data file - free format
c     --- initialize variables
      massflag = 0
      natomflag = 0
      nbondflag = 0
      nangleflag = 0
      ndihedflag = 0
      nimproflag = 0

      noncoeffflag = 0
      nbondcoeffflag = 0
      nanglecoeffflag = 0
      ndihedcoeffflag = 0
      nimprocoeffflag = 0

      nbondbondflag = 0
      nbondangleflag = 0
      nmidbondtorsionflag = 0
      nendbondtorsionflag = 0
      nangletorsionflag = 0
      nangleangletorsionflag = 0
      nbondbond13flag = 0
      nangleangleflag = 0

      ldone = .false.
      do while (.not. ldone)
c        --- read identifier string
         read (2,*,end=999,err=999)
         read (2,900,end=999,err=999) str
 900     format (a)
         read (2,*,end=999,err=999)
c        --- end of file 
         if (match('All Done',str,m)) then
            ldone = .true.
         elseif (match('Masses',str,m)) then
c           --- read atom masses
            call twh_read_mass(massflag,ntypes)
         elseif (match('Atoms',str,m)) then
c           --- read atoms
            call twh_read_atoms(lfailure,natomflag,natoms,box(1,1)
     &           ,box(1,2),box(1,3))
            if ( lfailure ) return
         elseif (match('Velocities',str,m)) then
c           --- read velocities
            if (natomflag .eq. 0) then
               write(6,*) 'Atoms entry must appear before Velocities'
               stop
            endif
            call twh_read_vels(natoms)
         elseif (match('Bonds',str,m)) then
c           --- read bonds
            if (natomflag .eq. 0) then
               write(6,*) 'Atoms entry must appear before Bonds'
               stop
            endif
            call twh_read_bonds(lfailure,nbondflag,nbonds,nbondtypes)
            if ( lfailure ) return
         elseif (match('Angles',str,m)) then
c           --- read angles
            if (natomflag .eq. 0) then
               write(6,*) 'Atoms entry must appear before Angles'
               stop
            endif
            call twh_read_angles(lfailure,nangleflag,nangles
     &           ,nangletypes)
            if ( lfailure ) return
        elseif (match('Dihedrals',str,m)) then
c         --- read dihedrals
          if (natomflag .eq. 0) then
             write(6,*) 'Atoms entry must appear before Dihedrals'
             stop
          endif
          call twh_read_dihedrals(lfailure,ndihedflag,ndihedrals
     &         ,ndihedtypes)
          if ( lfailure ) return
        elseif (match('Impropers',str,m)) then
c         ---read impropers
           if (natomflag .eq. 0) then
              write(6,*) 'Atoms entry must appear before Impropers'
              stop
           endif
           call twh_read_impropers(lfailure,nimproflag,nimpropers
     &          ,nimprotypes)
           if ( lfailure ) return

        elseif (match('Nonbond Coeffs',str,m)) then
c          --- read nonbond coeffs
           call twh_read_nonbond_coeff(noncoeffflag,nonstyle,ntypes)

        elseif (match('Bond Coeffs',str,m)) then
c          --- read bond coeffs
           call twh_read_bond_coeff(nbondcoeffflag,bondstyle,nbondtypes
     &          ,nbonds)

        elseif (match('Angle Coeffs',str,m)) then
c          --- read angle coeffs
           call twh_read_angle_coeff(nanglecoeffflag,lmpangstyle,nangles
     &          ,nangletypes)

        elseif (match('Dihedral Coeffs',str,m)) then
c          --- read dihedral coeffs
           call twh_read_dihedral_coeff(ndihedcoeffflag,dihedstyle
     &          ,ndihedrals,ndihedtypes)

        elseif (match('Improper Coeffs',str,m)) then
c          --- read improper coeffs
           call twh_read_improper_coeff(nimprocoeffflag,improstyle
     &          ,nimprotypes,nimpropers)

        elseif (match('BondBond Coeffs',str,m)) then
c          --- read bond-bond coeffs for class 2 force field
           call twh_read_bondbond_coeff(nbondbondflag,nangles
     &          ,lmpangstyle,nangletypes)

        elseif (match('BondAngle Coeffs',str,m)) then
c          --- read bond-angle coeffs for class 2 force field
           call twh_read_bondangle_coeff(nbondangleflag,nangles
     &          ,lmpangstyle,nangletypes)

        elseif (match('MiddleBondTorsion Coeffs',str,m)) then
c          --- read class 2 middle-bond-torsion coeffs
           call twh_read_midbondtorsion_coeff(nmidbondtorsionflag
     &          ,ndihedrals,dihedstyle,ndihedtypes)

        elseif (match('EndBondTorsion Coeffs',str,m)) then
c          --- read class 2 end-bond-torsion coeffs
           call twh_read_endbondtorsion_coeff(nendbondtorsionflag
     &          ,ndihedrals,dihedstyle,ndihedtypes)

        elseif (match('AngleTorsion Coeffs',str,m)) then
c          --- read class 2 angle-torsion coeffs
           call twh_read_angletorsion_coeff(nangletorsionflag,ndihedrals
     &          ,dihedstyle,ndihedtypes)

        elseif (match('AngleAngleTorsion Coeffs',str,m)) then
c          --- read class 2 angle-angle-torsion coeffs
           call twh_read_angleangletorsion_coeff(nangleangletorsionflag
     &          ,ndihedrals,dihedstyle,ndihedtypes)

        elseif (match('BondBond13 Coeffs',str,m)) then
c          --- read class 2 bondbond13 coeffs
           call twh_read_bondbond13_coeff(nbondbond13flag,ndihedrals
     &          ,dihedstyle,ndihedtypes)

        elseif (match('AngleAngle Coeffs',str,m)) then
c          --- read class 2 angle-angle coeffs
           call twh_read_angleangle_coeff(nangleangleflag
     &          ,nimpropers,improstyle,nimprotypes,naatypes)

        else
c          --- unknown identifier
           write (6,*) 'UNKNOWN: ',str
           write(6,*) 'Unknown identifier in data file'
           stop
        endif
      enddo
 999  continue

c     --- close lammps_data file
      close (2)

c     --- check that all required options were read in from data file
      if (massflag.eq.0) stop 'No Masses in data file'
c     --- assign the element names based on the masses
      eleflag = 'getelement'
      do itype = 1,ntypes
         localmass = twh_get_mass(itype)
         call twh_getelement(lfailure,eleflag,localmass,ename,ielement)
         if ( lfailure ) return
         call twh_elementname(GLB_SET,itype,ename)
      enddo

      if (natomflag.eq.0) stop 'No Atoms in data file'

      if (nbonds.gt.0.and.nbondflag.eq.0)
     &     stop'No Bonds in data file'
      if (nangles.gt.0.and.nangleflag.eq.0)
     &     stop'No Angles in data file'
      if (ndihedrals.gt.0.and.ndihedflag.eq.0)
     &     stop'No Dihedrals in data file'
      if (nimpropers.gt.0.and.nimproflag.eq.0)
     &     stop'No Impropers in data file'

      if (nbonds.gt.0.and.bondstyle.eq.5.and.nbondcoeffflag.eq.0)
     &     stop'No Bond Coeffs in data file'
      if (nangles.gt.0.and.lmpangstyle.eq.2.and.nanglecoeffflag.eq.0)
     &     stop'No Angle Coeffs in data file'
      if (ndihedrals.gt.0.and.dihedstyle.eq.2.and.ndihedcoeffflag.eq.0)
     &     stop'No Dihedral Coeffs in data file'
      if (nimpropers.gt.0.and.improstyle.eq.3.and.nimprocoeffflag.eq.0)
     &     stop'No Improper Coeffs in data file'

      if (nangles.gt.0.and.lmpangstyle.eq.2.and.nbondbondflag.eq.0)
     &     stop'No BondBond Coeffs in data file'
      if (nangles.gt.0.and.lmpangstyle.eq.2.and.nbondangleflag.eq.0)
     &     stop'No BondAngle Coeffs in data file'

      if (ndihedrals.gt.0.and.dihedstyle.eq.2.and.
     &     nmidbondtorsionflag.eq.0)
     &     stop'No MiddleBondTorsion Coeffs in data file'
      if (ndihedrals.gt.0.and.dihedstyle.eq.2.and.
     &     nendbondtorsionflag.eq.0)
     &     stop'No EndBondTorsion Coeffs in data file'
      if (ndihedrals.gt.0.and.dihedstyle.eq.2.and.
     &     nangletorsionflag.eq.0)
     &     stop'No AngleTorsion Coeffs in data file'
      if (ndihedrals.gt.0.and.dihedstyle.eq.2.and.
     &     nangleangletorsionflag.eq.0)
     &     stop'No AngleAngleTorsion Coeffs in data file'
      if (ndihedrals.gt.0.and.dihedstyle.eq.2.and.
     &     nbondbond13flag.eq.0)
     &     stop'No BondBond13 Coeffs in data file'

      if (nimpropers.gt.0.and.improstyle.eq.3.and.
     &     nangleangleflag.eq.0)
     &     stop'No AngleAngle Coeffs in data file'

c     --- check improstyle, if it is 3 then those impropers are also
c     --- angle angle terms.  LAMMPS defines an angle-angle and a 
c     --- Wilson oop improper for each and just sets a lot of force 
c     --- constants to 0.0 if it just wants one or the other
      if ( improstyle .eq. 3 ) then
         nangleangle = nimpropers
      endif

c     --- create the molecule data structures
      call twh_createmolecule(lfailure,natoms,nbonds,nangles,ndihedrals
     &     ,nimpropers,nangleangle,ndihedtypes)
      if ( lfailure ) return

      call twh_rwforcefield(lfailure,WRITE_FLAG,ntypes,nbondtypes
     &     ,nangletypes,ndihedtypes,nimprotypes,naatypes,nbitypes
     &     ,noftypes,nhbondtypes,1)
      if ( lfailure ) return

      return

      end

c     ******************************************************************
c     * assorted input functions to help this subroutine               *
c     ******************************************************************

      subroutine twh_read_mass(massflag,ntypes)
c     ******************************************************************
c     * read_mass                                                      *
c     * reads in the masses of each atom type in the system            *
c     *                                                                *
c     * rewritten from LAMMPS 2001 on 02-28-2002 by M.G. Martin        *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer massflag,ntypes
c     --- local variables
      integer itype,jtmp
      double precision localmass

      write (6,*) 'Masses ...'
      massflag = 1

      do itype = 1,ntypes
         read (2,*) jtmp,localmass
         call twh_mass(GLB_SET,itype,localmass)
      enddo

      return
      end


      subroutine twh_read_atoms(lfailure,natomflag,natoms
     &     ,xmin,ymin,zmin)
c     ******************************************************************
c     * read_atoms                                                     *
c     * reads in all of the atoms into temporary storage for use later *
c     * when we sort them and put them into the permenant arrays       *
c     *                                                                *
c     * rewritten from LAMMPS 2001 on 02-28-2002 by M.G. Martin        *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer natomflag,natoms
      double precision xmin,ymin,zmin
c     --- local variables
      integer iatom,itmp,imolec,itype
      double precision qtmp,xtmp,ytmp,ztmp
      character*30 filename

      write (6,*) 'Atoms ...'
      natomflag = 1
      filename = 'lammps_atoms'
      call twh_openfile(lfailure,13,filename,'DEFAULT','UNKNOWN'
     &     ,'DEFAULT')
      if ( lfailure ) return

      do iatom = 1,natoms
         read(2,*) itmp,imolec,itype,qtmp,xtmp,ytmp,ztmp
c        --- adjust coordinates to get 0.0 as lower boundary of box
         write(13,*) itmp,itype,qtmp,xtmp-xmin,ytmp-ymin,ztmp-zmin
      enddo
      rewind(13)
      close(13)
      return
      end


      subroutine twh_read_vels(natoms)
c     ******************************************************************
c     * read_vels                                                      *
c     * we do not need the velecities so this just skips past all of   *
c     * that data                                                      *
c     *                                                                *
c     * rewritten from LAMMPS 2001 on 02-28-2002 by M.G. Martin        *
c     * last modified 02-28-2002 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer natoms
c     --- local variables
      integer iatom

      write (6,*) 'Velocities ...'

      do iatom = 1,natoms
         read(2,*)
      enddo
      
      return
      end

      subroutine twh_read_bonds(lfailure,nbondflag,nbonds,nbondtypes)
c     ******************************************************************
c     * read_bonds                                                     *
c     * reads the bond topology and types and stores this to a file    *
c     *                                                                *
c     * rewritten from LAMMPS 2001 on 02-28-2002 by M.G. Martin        *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer nbondflag,nbonds,nbondtypes
c     --- local variables
      integer ibond,itmp,itype,iatom,jatom
      character*30 filename

      write (6,*) 'Bonds ...'
      if (nbonds.eq.0.or.nbondtypes.eq.0) stop'Should be no Bonds entry'
      nbondflag = 1
c     --- check array size
      if ( nbondtypes .gt. TVIBMAX ) then
         write(6,*) 'nbondtypes exceeds TVIBMAX'
         stop
      endif

c     --- open a scratch file for output to avoid the large arrays we
c     --- would need to set up this information
      filename = 'lammps_bonds'
      call twh_openfile(lfailure,13,filename,'DEFAULT','UNKNOWN'
     &     ,'DEFAULT')
      if ( lfailure ) return
      do ibond = 1,nbonds
         read(2,*) itmp,itype,iatom,jatom
         write(13,*) iatom,jatom,itype
      enddo

      rewind(13)
      close(13)

      return
      end

      subroutine twh_read_angles(lfailure,nangleflag,nangles
     &     ,nangletypes)
c     ******************************************************************
c     * read_angles                                                    *
c     * reads in the angle and their types and then stores this to a   *
c     * file                                                           *
c     *                                                                *
c     * rewritten from LAMMPS 2001 on 02-28-2002 by M.G. Martin        *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer nangleflag,nangles,nangletypes
c     --- local variables
      integer iangle,atype,atom1,atom2,atom3,idum
      character*30 filename

      write (6,*) 'Angles ...'
      if (nangles.eq.0.or.nangletypes.eq.0)
     &     stop'Should be no Angles entry'
      nangleflag = 1

      filename = 'lammps_angles'
      call twh_openfile(lfailure,13,filename,'DEFAULT','UNKNOWN'
     &     ,'DEFAULT')
      if ( lfailure ) return

      do iangle = 1,nangles
         read(2,*) idum,atype,atom1,atom2,atom3
         write(13,*) atom1,atom2,atom3,atype
      enddo

      rewind(13)
      close(13)

      return
      end

      subroutine twh_read_dihedrals(lfailure,ndihedflag,ndihedrals
     &     ,ndihedtypes)
c     ******************************************************************
c     * read_dihedral                                                  *
c     * reads in the torsional terms and stores this to a file         *
c     *                                                                *
c     * rewritten from LAMMPS 2001 on 02-28-2002 by M.G. Martin        *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ndihedflag,ndihedrals,ndihedtypes
c     --- local variables
      integer idihedral,ttype,atom1,atom2,atom3,atom4,idum
      character*30 filename

      write(6,*) 'Dihedrals ...'
      if (ndihedrals.eq.0.or.ndihedtypes.eq.0)
     &     stop'Should be no Dihedrals entry'
      ndihedflag = 1

      filename = 'lammps_dihedrals'
      call twh_openfile(lfailure,13,filename,'DEFAULT','UNKNOWN'
     &     ,'DEFAULT')
      if ( lfailure ) return

      do idihedral = 1,ndihedrals
         read(2,*) idum,ttype,atom1,atom2,atom3,atom4
         write(13,*) atom1,atom2,atom3,atom4,ttype
      enddo

      rewind(13)
      close(13)

      return
      end

      subroutine twh_read_impropers(lfailure,nimproflag,nimpropers
     &     ,nimprotypes)
c     ******************************************************************
c     * read_impropers                                                 *
c     * reads in the improper terms and writes them to a file          *
c     *                                                                *
c     * rewritten from LAMMPS 2001 on 02-28-2002 by M.G. Martin        *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer nimproflag,nimpropers,nimprotypes
c     --- local variables
      integer idum,type,atom1,atom2,atom3,atom4,impropers
      character*30 filename

      write (6,*) 'Impropers ...'
      if (nimpropers.eq.0.or.nimprotypes.eq.0)
     &     stop'Should be no Impropers entry'
      nimproflag = 1

      filename = 'lammps_impropers'
      call twh_openfile(lfailure,13,filename,'DEFAULT','UNKNOWN'
     &     ,'DEFAULT')
      if ( lfailure ) return

      do impropers = 1,nimpropers
         read(2,*) idum,type,atom1,atom2,atom3,atom4
         write(13,*) atom1,atom2,atom3,atom4,type
      enddo

      rewind(13)
      close(13)

      return
      end

      subroutine twh_read_nonbond_coeff(noncoeffflag,nonstyle,ntypes)
c     ******************************************************************
c     * read_nonbond_coeff                                             *
c     * reads in the nonbonded force field coefficients and turns them *
c     * into the proper towhee formats                                 *
c     *                                                                *
c     * rewritten from LAMMPS 2001 on 02-28-2002 by M.G. Martin        *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed/from subroutine
      integer noncoeffflag,nonstyle,ntypes
c     --- local variables
      integer i,ii,k,jtmp
      double precision tmp
      dimension tmp(4)
      double precision nbclocal(MINNBCOEFF:MAXNBCOEFF)
      character*(FFNAMELEN) ctemp

      write (6,*) 'Nonbond Coeffs ...'
      if (nonstyle.eq.0)
     &     stop'Nonbond Coeffs entry not allowed for nonstyle 0'
      noncoeffflag = 1

      do i = 1,ntypes
         ii = twh_getnbtype(i,i)
         if (nonstyle.eq.1) then
c           --- lj/cutoff, this is the standard 12-6 LJ in Towhee
            call twh_classical_potential(GLB_SET,'Lennard-Jones')
            read (2,*) jtmp,(tmp(k),k=1,2)
c           --- sigma
            nbclocal(1) = tmp(2)
c           --- epsilon (in kcal)
            nbclocal(2) = tmp(1)*CON_KCALTOK
c           --- set the 1-4 terms to the same values
            nbclocal(3) = nbclocal(1)
            nbclocal(4) = nbclocal(2)
         elseif (nonstyle.eq.2) then
c           --- lj/switch, no complete equivalent in Towhee so just use 
c           --- the standard 12-6 LJ
            call twh_classical_potential(GLB_SET,'Lennard-Jones')
            read (2,*) jtmp,(tmp(k),k=1,2)
c           --- sigma
            nbclocal(1) = tmp(2)
c           --- epsilon (in kcal)
            nbclocal(2) = tmp(1)*CON_KCALTOK
c           --- set the 1-4 terms to the same values
            nbclocal(3) = nbclocal(1)
            nbclocal(4) = nbclocal(2) 
         elseif (nonstyle.eq.3) then
c           --- lj/shift unlike towhee this adjusts the distances instead
c           --- of applying an energy shift.  Just treat as standard 
c           --- 12-6 LJ and ignore the third term
            call twh_classical_potential(GLB_SET,'Lennard-Jones')
            read (2,*) jtmp,(tmp(k),k=1,3)
c           --- sigma
            nbclocal(1) = tmp(2)
c           --- epsilon (in kcal)
            nbclocal(2) = tmp(1)*CON_KCALTOK
c           --- set the 1-4 terms to the same values
            nbclocal(3) = nbclocal(1)
            nbclocal(4) = nbclocal(2)
         elseif (nonstyle.eq.4) then
c           --- soft, not currently supported in Towhee
            stop 'cannot read in Lammps soft potential'
         elseif (nonstyle.eq.5) then
            call twh_classical_potential(GLB_SET,'9-6')
c           --- class1/cutoff, same as LJ 9-6 in Towhee
            read (2,*) jtmp,(tmp(k),k=1,2)
c           --- sigma/Ro
            nbclocal(1) = tmp(2)
c           --- epsilon (kcal)
            nbclocal(2) = tmp(1)*CON_KCALTOK
         elseif (nonstyle.eq.6) then
c           --- lj/charmm, same as Towhee 12-6 LJ where there are
c           --- different 1-4 terms from the regular nonbond
            call twh_classical_potential(GLB_SET,'Lennard-Jones')
            read (2,*) jtmp,(tmp(k),k=1,4)
c           --- sigma
            nbclocal(1) = tmp(2)
c           --- epsilon (in kcal)
            nbclocal(2) = tmp(1)*CON_KCALTOK
c           --- 1-4 sigma
            nbclocal(3) = tmp(4)
c           --- 1-4j epsilon (kcal)
            nbclocal(4) = tmp(3)*CON_KCALTOK
         endif
         call twh_set_nbcoeff(ii,nbclocal)
c        --- there is no support for automatic typing so set all of the
c        --- atom typing variables to "no name" and set the forcefield
c        --- to "none"
         ctemp = 'none'
         call twh_nbff(GLB_SET,i,ctemp)
         ctemp = 'no name'
         call twh_nbname(GLB_SET,i,ctemp)
         call twh_cbname(GLB_SET,i,ctemp)
         call twh_caname(GLB_SET,i,ctemp)
         call twh_ctname(GLB_SET,i,ctemp)
      enddo
      return
      end

      subroutine twh_read_bond_coeff(nbondcoeffflag,bondstyle,nbondtypes
     &     ,nbonds)
c     ******************************************************************
c     * read_bond_coeff                                                *
c     * reads in the bond coefficients                                 *
c     * rewritten from LAMMPS 2001 on 03-01-2002 by M.G. Martin        *
c     * last modified 07-14-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer nbondcoeffflag,bondstyle,nbondtypes,nbonds
c     --- local variables
      integer i,j,jtmp
      double precision tmp
      dimension tmp(5)

      write (6,*) 'Bond Coeffs ...'
      if (nbonds.eq.0.or.bondstyle.eq.0)
     &     stop'Bond Coeffs entry not allowed'
      nbondcoeffflag = 1

      if (bondstyle.eq.1) then
c        --- standard harmonic
         do i = 1,nbondtypes
            read (2,*) jtmp,(tmp(j),j=1,2)
c           --- force constant (kcal)
            call twh_vibcoeff(GLB_SET,i,1,tmp(1)*CON_KCALTOK)
c           --- equilibrium value
            call twh_vibcoeff(GLB_SET,i,0,tmp(2))
c           --- this is towhee style 2
            call twh_vibstyle(GLB_SET,i,BOND_STYLE_HARMONIC)
         enddo
      elseif (bondstyle.eq.2) then
c        --- FENE/standard (not in Towhee)
         stop 'FENE/standard is not in Towhee'
      elseif (bondstyle.eq.3) then
c        --- FENE/shift (not in Towhee)
         stop 'FENE/shift is not in Towhee'
      elseif (bondstyle.eq.4) then
c        --- nonlinear, same as Towhee nonlinear
         do i = 1,nbondtypes
            read (2,*) jtmp,(tmp(j),j=1,3)
c           --- force constant (kcal)
            call twh_vibcoeff(GLB_SET,i,1,tmp(1)*CON_KCALTOK)
c           --- equilibrium value
            call twh_vibcoeff(GLB_SET,i,0,tmp(2))
c           --- lamda squared
            call twh_vibcoeff(GLB_SET,i,2,tmp(3)**2)
c           --- this is towhee style 4
            call twh_vibstyle(GLB_SET,i,BOND_STYLE_NONLINEAR)
         enddo

      elseif (bondstyle.eq.5) then
c        --- Class 2 (Compass quartic in Towhee)
         do i = 1,nbondtypes
            read (2,*) jtmp,(tmp(j),j=1,4)
c           --- equilibrium distance
            call twh_vibcoeff(GLB_SET,i,0,tmp(1))
c           --- square force constant (kcal)
            call twh_vibcoeff(GLB_SET,i,1,tmp(2)*CON_KCALTOK)
c           --- cubic force constant (kcal)
            call twh_vibcoeff(GLB_SET,i,2,tmp(3)*CON_KCALTOK)
c           --- quartic force constant (kcal)
            call twh_vibcoeff(GLB_SET,i,3,tmp(4)*CON_KCALTOK)
c           --- quartic is type 6 in Towhee
            call twh_vibstyle(GLB_SET,i,BOND_STYLE_COMPASS)
         enddo
      else
         write(6,*) 'invalid bondstyle ',bondstyle
         stop
      endif
      
      return
      end

      subroutine twh_read_angle_coeff(nanglecoeffflag,lmpangstyle
     &     ,nangles,nangletypes)
c     ******************************************************************
c     * read_angle_coeff                                               *
c     * reads in the angle coefficients                                *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-01-2002 by M.G. Martin           *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer nanglecoeffflag,lmpangstyle,nangles,nangletypes
c     --- local variables
      integer i,j,jtmp,iset
      double precision tmp
      dimension tmp(4)
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)

      write (6,*) 'Angle Coeffs ...'
      if (nangles.eq.0.or.lmpangstyle.eq.0)
     &     stop'Angle Coeffs entry not allowed'
      nanglecoeffflag = 1

c     --- make sure we have enough room to read in the angles
      if ( nangletypes .gt. TBENMAX ) then
         write(6,*) 'angle types exceeds TBENMAX'
         stop
      endif

      if (lmpangstyle.eq.1) then
c        --- harmonic (also harmonic in Towhee)
         do i = 1,nangletypes
            read (2,*) jtmp,(tmp(j),j=1,2)
c           --- force constant (kcal/radian^2)
            bclocal(1) = tmp(1)*CON_KCALTOK
            iset = 1
            call twh_bencoeff(GLB_SET,i,iset,bclocal(iset))
c           --- equil distance (degrees) - later converted to radians
            bclocal(0) = tmp(2)
            iset = 0
            call twh_bencoeff(GLB_SET,i,iset,bclocal(iset))
c           --- this is anglestyle 1 in Towhee
            call twh_anglestyle(GLB_SET,i,ANG_STYLE_HARMONIC)
         enddo
      elseif (lmpangstyle.eq.2) then
c        --- class 2 This is part of the Compass quartic angle in Towhee
c        --- later will add in the bond-angle terms
c        --- use the style 8 version as the terms are explicitly 
c        --- declared
         do i = 1,nangletypes
            read (2,*) jtmp,(tmp(j),j=1,4)
c           --- equilibrium angle (degrees), later converted to radians
            bclocal(0) = tmp(1)
            iset = 0
            call twh_bencoeff(GLB_SET,i,iset,bclocal(iset))
c           --- square force constant (kcal)
            bclocal(1) = tmp(2)*CON_KCALTOK
            iset = 1
            call twh_bencoeff(GLB_SET,i,iset,bclocal(iset))
c           --- cubic force constant (kcal)
            bclocal(2) = tmp(3)*CON_KCALTOK
            iset = 2
            call twh_bencoeff(GLB_SET,i,iset,bclocal(iset))
c           --- quartic force constant (kcal)
            bclocal(3) = tmp(4)*CON_KCALTOK
            iset = 3
            call twh_bencoeff(GLB_SET,i,iset,bclocal(iset))
c           --- this is anglestyle 8 in Towhee
            call twh_anglestyle(GLB_SET,i,ANG_STYLE_COMPASS_EX)
c           --- set the lbonang term to false until we sort it out
            call twh_lbonang(GLB_SET,i,.false.)
         enddo
      elseif (lmpangstyle.eq.3) then
         do i = 1,nangletypes
c           --- Charmm angles, also charmm with Urey-Bradley in Towhee
            read (2,*) jtmp,(tmp(j),j=1,4)
c           --- equilibrium angle (degrees), later convert to radians
            bclocal(0) = tmp(2)
            iset = 0
            call twh_bencoeff(GLB_SET,i,iset,bclocal(iset))
c           --- angle force constant (kcal)
            bclocal(1) = tmp(1)*CON_KCALTOK
            iset = 1
            call twh_bencoeff(GLB_SET,i,iset,bclocal(iset))
c           --- Urey-Bradley distance
            bclocal(2) = tmp(4)
            iset = 2
            call twh_bencoeff(GLB_SET,i,iset,bclocal(iset))
c           --- UB force constant (kcal)
            bclocal(3) = tmp(3)*CON_KCALTOK
            iset = 3
            call twh_bencoeff(GLB_SET,i,iset,bclocal(iset))
c           --- this is anglestyle 5
            call twh_anglestyle(GLB_SET,i,ANG_STYLE_CHARMM_UB)
         enddo
      endif

      return
      end

      subroutine twh_read_dihedral_coeff(ndihedcoeffflag,dihedstyle
     &     ,ndihedrals,ndihedtypes)
c     ******************************************************************
c     * read_dihedral_coeff                                            *
c     * reads in the dihedral coefficients                             *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-01-2002 by M.G. Martin           *
c     * last modified 04-17-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_ARCCOS
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ndihedcoeffflag,dihedstyle,ndihedrals,ndihedtypes
c     --- local variables
      integer i,j,jtmp
      double precision tmp,onepi
      dimension tmp(6)
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)

      call twh_constant_pi(GLB_GET,onepi)

      write (6,*) 'Dihedral Coeffs ...'
      if (ndihedrals.eq.0.or.dihedstyle.eq.0)
     &     stop'Dihedral Coeffs entry not allowed'
      ndihedcoeffflag = 1
c     --- check for space in the arrays
      if ( ndihedtypes .gt. TTORMAX ) then
         write(6,*) 'torsion types exceeds TTORMAX'
         stop
      endif

      if (dihedstyle.eq.1) then
c        --- LAMMPS harmonic, can be rewritten into the Amber/Charmm/
c        --- Gromos general Cosine series
         do i = 1,ndihedtypes
            read (2,*) jtmp,(tmp(j),j=1,3)
c           --- assign the type for the Amber/Charmm/Gromos Cos series
            call twh_torstyle(GLB_SET,i,3)
c           --- only a single loop
            call twh_ntorloop(GLB_SET,i,1)
c           --- set default 1-4 to false and zero
            call twh_loftor(GLB_SET,i,.false.)
            tclocal(-1) = 0.0d0
c           --- prefactor force constant (kcal)
            tclocal(1) = tmp(1)*CON_KCALTOK
c           --- angle coefficient inside the cosine
            tclocal(2) = tmp(3)
c           --- LAMMPS uses a d constant before the cosine, this is 
c           --- rewritten in Towhee to be an adjustment inside the 
c           --- cosine
            tclocal(3) = -twh_arccos(tmp(2))
            call twh_set_torcoeff(i,tclocal)
         enddo
      elseif (dihedstyle.eq.2) then
c        --- LAMMPS class 2, same as Towhee Class 2/Compass with
c        --- the cross terms and explicit declaration of ordering.
         do i = 1,ndihedtypes
            read (2,*) jtmp,(tmp(j),j=1,6)
c           --- assign the type for explicit class 2
            call twh_torstyle(GLB_SET,i,9)
c           --- by default assume 1-4 interactions with 1.0 scale
            call twh_loftor(GLB_SET,i,.true.)
            tclocal(-1) = 1.0d0
c           --- assign the force constants (kcal)
            tclocal(0) = tmp(1)*CON_KCALTOK
            tclocal(1) = tmp(3)*CON_KCALTOK
            tclocal(2) = tmp(5)*CON_KCALTOK
c           --- assign the shift factors (convert from degrees to rad)
            tclocal(3) = tmp(2)*onepi/180.0d0
            tclocal(4) = tmp(4)*onepi/180.0d0
            tclocal(5) = tmp(6)*onepi/180.0d0
            call twh_set_torcoeff(i,tclocal)
         enddo
      elseif (dihedstyle.eq.3) then
c        --- LAMMPS multiharmonic, same as the Towhee Cosine Power 
c        --- Series
         do i = 1,ndihedtypes
            read (2,*) jtmp,(tmp(j),j=1,5)
            call twh_torstyle(GLB_SET,i,10)
c           --- assume no 1-4 interactions
            call twh_loftor(GLB_SET,i,.false.)
            tclocal(-1) = 0.0d0
c           --- LAMMPS default is a 5 term series, which means a torloop
c           --- of 4 in towhee as towhee starts with 0
            call twh_ntorloop(GLB_SET,i,4)
c           --- assign the 5 force constants (kcal)
            tclocal(0) = tmp(1)*CON_KCALTOK
            tclocal(1) = tmp(2)*CON_KCALTOK
            tclocal(2) = tmp(3)*CON_KCALTOK
            tclocal(3) = tmp(4)*CON_KCALTOK
            tclocal(4) = tmp(5)*CON_KCALTOK
            call twh_set_torcoeff(i,tclocal)
         enddo
      elseif (dihedstyle.eq.4) then
c        --- LAMMPS charmm, same as Amber/Charmm/Gromos in Towhee
         do i = 1,ndihedtypes
            read (2,*) jtmp,(tmp(j),j=1,4)
c           --- Amber/Charmm/Gromos style
            call twh_torstyle(GLB_SET,i,3)
c           --- this potential has 1-4 interactions if tmp(4) is .gt. 0
            if ( tmp(4) .gt. 0.0d0 ) then
               call twh_loftor(GLB_SET,i,.true.)
               tclocal(-1) = tmp(4)
            else
               call twh_loftor(GLB_SET,i,.false.)
            endif
c           --- just one loop
            call twh_ntorloop(GLB_SET,i,1)
c           --- force constant (kcal)
            tclocal(1) = tmp(1)*CON_KCALTOK
c           --- interior multiplier
            tclocal(2) = tmp(2)
c           --- shift (sign different) (convert from degrees to radians)
            tclocal(3) = -tmp(3)*onepi/180.0d0
            call twh_set_torcoeff(i,tclocal)
         enddo
      endif

      return
      end

      subroutine twh_read_improper_coeff(nimprocoeffflag,improstyle
     &     ,nimprotypes,nimpropers)
c     ******************************************************************
c     * read_improper_coeff                                            *
c     * reads the improper torsion coefficients from a file and creates*
c     * a data structure that holds them.                              *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-04-2002 by M.G. Martin           *
c     * last modified 04-17-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_ARCCOS
#include "functions.h" 
c     --- variables passed to/from the subroutine
      integer nimprocoeffflag,improstyle,nimprotypes
      integer nimpropers
c     --- local variables
      integer i,j,jtmp
      double precision tmp,onepi
      dimension tmp(3)
      double precision iclocal
      dimension iclocal(MINIMPCOEFF:MAXIMPCOEFF)

      call twh_constant_pi(GLB_GET,onepi)

      write (6,*) 'Improper Coeffs ...'
      if (nimpropers.eq.0.or.improstyle.eq.0)
     &     stop'Improper Coeffs entry not allowed'
      nimprocoeffflag = 1
c     --- safety check nimprotypes
      if ( nimprotypes .gt. TIMPMAX ) then
         write(6,*) 'improper types exceeds TIMPMAX '
         stop
      endif

      if (improstyle .eq. 1) then
c        --- LAMMPS harmonic, this one is the same as Towhee harmonic
         do i = 1,nimprotypes
            read (2,*) jtmp,(tmp(j),j=1,2)
c           --- this is the normal oop form
            call twh_impform(GLB_SET,i,2)
c           --- harmonic is torsion type 1
            call twh_impstyle(GLB_SET,i,1)
c           --- energy constant (kcal)
            iclocal(0) = tmp(1)*CON_KCALTOK
c           --- angle (degrees, change to radians)
            iclocal(1) = tmp(2)*onepi/180.0d0
c           --- set the improper coefficients
            call twh_set_impcoeff(i,iclocal)
         enddo
      elseif (improstyle.eq.2) then
c        --- LAMMPS cvff, same as Amber/Charmm/Gromos cos series
         do i = 1,nimprotypes
            read (2,*) jtmp,(tmp(j),j=1,3)
c           --- this is the normal oop form
            call twh_impform(GLB_SET,i,2)
c           --- this is torsion type 3
            call twh_impstyle(GLB_SET,i,3)
c           --- prefactor (kcal)
            iclocal(1) = tmp(1)*CON_KCALTOK
c           --- cosine multiplier, move inside the cos
            iclocal(2) = -twh_arccos(tmp(2))
c           --- interior cosine factor
            iclocal(3) = tmp(3)
c           --- set the improper coefficients
            call twh_set_impcoeff(i,iclocal)
         enddo
      elseif (improstyle.eq.3) then
c        --- LAMMPS class 2, same as wilson oop harmonic
         do i = 1,nimprotypes
            read (2,*) jtmp,(tmp(j),j=1,2)
c           --- this is the wilson oop form
            call twh_impform(GLB_SET,i,3)
c           --- wilson harmonic is torsion type 5
            call twh_impstyle(GLB_SET,i,5)
c           --- energy constant (kcal)
            iclocal(0) = tmp(1)*CON_KCALTOK
c           --- angle (degrees, change to radians)
            iclocal(1) = tmp(2)*onepi/180.0d0
c           --- set the improper coefficients
            call twh_set_impcoeff(i,iclocal)
         enddo
      endif
      return
      end

      subroutine twh_read_bondbond_coeff(nbondbondflag,nangles
     &     ,lmpangstyle,nangletypes)
c     ******************************************************************
c     * read_bondbond_coeff                                            *
c     * reads in the bondbond coefficients and sets up the bondbond    *
c     * structures                                                     *
c     * it is assumed that each bondbond type goes with the same angle *
c     * type                                                           *
c     *                                                                *
c     * rewritten from LAMMPS 2001 by M.G. Martin                      *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer nbondbondflag,nangles,lmpangstyle,nangletypes
c     --- local variables
      integer i,j,jtmp,iset
      double precision tmp
      dimension tmp(3)
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)

      write (6,*) 'BondBond Coeffs ...'
      if (nangles.eq.0.or.lmpangstyle.ne.2)
     &     stop'BondBond Coeffs entry not allowed'
      nbondbondflag = 1

      do i = 1,nangletypes
c        --- reset the bondbond logical
         call twh_lbonbon(GLB_SET,i,.true.)
         read (2,*) jtmp,(tmp(j),j=1,3)
c        --- put these into the angles
         call twh_lbonbon(GLB_SET,i,.true.)
c        --- force constant (kcal)
         bclocal(8) = tmp(1)*CON_KCALTOK
         iset = 8
         call twh_bencoeff(GLB_SET,i,iset,bclocal(iset))
c        --- equilibrium distances
         bclocal(9) = tmp(2)
         iset = 9
         call twh_bencoeff(GLB_SET,i,iset,bclocal(iset))
         bclocal(10) = tmp(3)
         iset = 10
         call twh_bencoeff(GLB_SET,i,iset,bclocal(iset))
      enddo

      return
      end

      subroutine twh_read_bondangle_coeff(nbondangleflag,nangles
     &     ,lmpangstyle,nangletypes)
c     ******************************************************************
c     * read_bondangle_coeff                                           *
c     * reads in the bondangle coefficients                            *
c     * assumes the same angle order and same atom order as the angle  *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-04-2002 by M.G. Martin           *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer nbondangleflag,nangles,lmpangstyle,nangletypes
c     --- local variables
      integer i,j,jtmp,iset
      double precision tmp
      dimension tmp(4)
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)

      write (6,*) 'BondAngle Coeffs ...'
      if (nangles.eq.0.or.lmpangstyle.ne.2)
     &     stop'BondAngle Coeffs entry not allowed'
      nbondangleflag = 1

      do i = 1,nangletypes
c        --- reset the bond-angle logical
         call twh_lbonang(GLB_SET,i,.true.)
         read (2,*) jtmp,(tmp(j),j=1,4)
c        --- assign the force constants (kcal)
         bclocal(4) = tmp(1)*CON_KCALTOK
         iset = 4
         call twh_bencoeff(GLB_SET,i,iset,bclocal(iset))
         bclocal(6) = tmp(2)*CON_KCALTOK
         iset = 6
         call twh_bencoeff(GLB_SET,i,iset,bclocal(iset))
c        --- assign the equilibrium distances
         bclocal(5) = tmp(3)
         iset = 5
         call twh_bencoeff(GLB_SET,i,iset,bclocal(iset))
         bclocal(7) = tmp(4)
         iset = 7
         call twh_bencoeff(GLB_SET,i,iset,bclocal(iset))
      enddo

      return
      end

      subroutine twh_read_midbondtorsion_coeff(nmidbondtorsionflag
     &     ,ndihedrals,dihedstyle,ndihedtypes)
c     ******************************************************************
c     * read_midbondtorsion_coeff                                      *
c     * read midbondtorsion coeffs from data file                      *
c     * assumes are in same order as torsions                          *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-04-2002 by M.G. Martin           *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer nmidbondtorsionflag,ndihedrals,dihedstyle,ndihedtypes
c     --- local variables
      integer i,j,jtmp,index
      double precision tclocal,tmp
      dimension tmp(4)

      write (6,*) 'MiddleBondTorsion Coeffs ...'
      if (ndihedrals.eq.0.or.dihedstyle.ne.2) stop
     &     'MiddleBondTorsion Coeffs entry not allowed'
      nmidbondtorsionflag = 1

      do i = 1,ndihedtypes
         read (2,*) jtmp,(tmp(j),j=1,4)
c        --- assign the force constants manually so we do not overwrite
c        --- those assigned for the base part of this torsion
         index = 6
         tclocal = tmp(1)*CON_KCALTOK
         call twh_torcoeff(GLB_SET,i,index,tclocal)
         index = 7
         tclocal = tmp(2)*CON_KCALTOK
         call twh_torcoeff(GLB_SET,i,index,tclocal)
         index = 8
         tclocal = tmp(3)*CON_KCALTOK
         call twh_torcoeff(GLB_SET,i,index,tclocal)
c        --- assign the distance
         index = 9
         tclocal = tmp(4)
         call twh_torcoeff(GLB_SET,i,index,tclocal)
      enddo

      return
      end

      subroutine twh_read_endbondtorsion_coeff(nendbondtorsionflag
     &     ,ndihedrals,dihedstyle,ndihedtypes)
c     ******************************************************************
c     * read_endbondtorsion_coeff                                      *
c     * read endbondtorsion coeffs from data file                      *
c     * assumes exact same atom order as the torsions                  *
c     *                                                                *
c     * rewritten from LAMMPS 2002 03-04-2002 by M.G. Martin           *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer nendbondtorsionflag,ndihedrals,dihedstyle,ndihedtypes
c     --- local variables
      integer i,jtmp,j,index
      double precision tclocal,tmp
      dimension tmp(8)

      write (6,*) 'EndBondTorsion Coeffs ...'
      if (ndihedrals.eq.0.or.dihedstyle.ne.2) stop
     &     'EndBondTorsion Coeffs entry not allowed'
      nendbondtorsionflag = 1

      do i = 1,ndihedtypes
         read (2,*) jtmp,(tmp(j),j=1,8)
c        --- assign the force constants (kcal)
         index = 10
         tclocal = tmp(1)*CON_KCALTOK
         call twh_torcoeff(GLB_SET,i,index,tclocal)
         index = 11
         tclocal = tmp(2)*CON_KCALTOK
         call twh_torcoeff(GLB_SET,i,index,tclocal)
         index = 12
         tclocal = tmp(3)*CON_KCALTOK
         call twh_torcoeff(GLB_SET,i,index,tclocal)
         index = 14
         tclocal = tmp(4)*CON_KCALTOK
         call twh_torcoeff(GLB_SET,i,index,tclocal)
         index = 15
         tclocal = tmp(5)*CON_KCALTOK
         call twh_torcoeff(GLB_SET,i,index,tclocal)
         index = 16
         tclocal = tmp(6)*CON_KCALTOK
         call twh_torcoeff(GLB_SET,i,index,tclocal)
c        --- assign the distances
         index = 13
         tclocal = tmp(7)
         call twh_torcoeff(GLB_SET,i,index,tclocal)
         index = 17
         tclocal = tmp(8)
         call twh_torcoeff(GLB_SET,i,index,tclocal)
      enddo

      return
      end

      subroutine twh_read_angletorsion_coeff(nangletorsionflag
     &     ,ndihedrals,dihedstyle,ndihedtypes)
c     ******************************************************************
c     * read_angletorsion_coeff                                        *
c     * read angletorsion coeffs from data file                        *
c     * assumes same order as in the torsion                           *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-04-2002 by M.G. Martin           *
c     * last modified 06-14-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer nangletorsionflag,ndihedrals,dihedstyle,ndihedtypes
c     --- local variables
      integer i,jtmp,j,index
      double precision tclocal,tmp,onepi
      dimension tmp(8)

      call twh_constant_pi(GLB_GET,onepi)

      write (6,*) 'AngleTorsion Coeffs ...'
      if (ndihedrals.eq.0.or.dihedstyle.ne.2) stop
     &     'AngleTorsion Coeffs entry not allowed'
      nangletorsionflag = 1

      do i = 1,ndihedtypes
         read (2,*) jtmp,(tmp(j),j=1,8)
c        --- assign force constants (kcal)
         index = 18
         tclocal = tmp(1)*CON_KCALTOK
         call twh_torcoeff(GLB_SET,i,index,tclocal)
         index = 19
         tclocal = tmp(2)*CON_KCALTOK
         call twh_torcoeff(GLB_SET,i,index,tclocal)
         index = 20
         tclocal = tmp(3)*CON_KCALTOK
         call twh_torcoeff(GLB_SET,i,index,tclocal)
         index = 21
         tclocal = tmp(4)*CON_KCALTOK
         call twh_torcoeff(GLB_SET,i,index,tclocal)
         index = 22
         tclocal = tmp(5)*CON_KCALTOK
         call twh_torcoeff(GLB_SET,i,index,tclocal)
         index = 23
         tclocal = tmp(6)*CON_KCALTOK
         call twh_torcoeff(GLB_SET,i,index,tclocal)
c        --- assign angles (convert from degrees)
         index = 24
         tclocal = tmp(7)*onepi/180.0d0
         call twh_torcoeff(GLB_SET,i,index,tclocal)
         index = 25
         tclocal = tmp(8)*onepi/180.0d0
         call twh_torcoeff(GLB_SET,i,index,tclocal)
      enddo

      return
      end

      subroutine twh_read_angleangletorsion_coeff(nangleangletorsionflag
     &     ,ndihedrals,dihedstyle,ndihedtypes)
c     ******************************************************************
c     * read_angleangletorsion_coeff                                   *
c     * read angleangletorsion coeffs from data file                   *
c     * assumes same order as the torsions                             *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-04-2002 by M.G. Martin           *
c     * last modified 06-14-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer nangleangletorsionflag,ndihedrals,dihedstyle,ndihedtypes
c     --- local variables
      integer i,jtmp,j,index
      double precision tclocal,tmp,onepi
      dimension tmp(3)

      write (6,*) 'AngleAngleTorsion Coeffs ...'
      if (ndihedrals.eq.0.or.dihedstyle.ne.2) stop
     &     'AngleAngleTorsion Coeffs entry not allowed'
      nangleangletorsionflag = 1

      call twh_constant_pi(GLB_GET,onepi)

      do i = 1,ndihedtypes
         read (2,*) jtmp,(tmp(j),j=1,3)
c        --- assign the force constant (kcal)
         index = 26
         tclocal = tmp(1)*CON_KCALTOK
         call twh_torcoeff(GLB_SET,i,index,tclocal)
c        --- assign the angles (convert from degrees)
         index = 27
         tclocal = tmp(2)*onepi/180.0d0
         call twh_torcoeff(GLB_SET,i,index,tclocal)
         index = 28
         tclocal = tmp(3)*onepi/180.0d0
         call twh_torcoeff(GLB_SET,i,index,tclocal)
      enddo

      return
      end

      subroutine twh_read_bondbond13_coeff(nbondbond13flag,ndihedrals
     &     ,dihedstyle,ndihedtypes)
c     ******************************************************************
c     * read_bondbond13_coeff                                          *
c     * read bondbond13 coeffs from data file                          *
c     * assumes same order as torsions                                 *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-04-2002 by M.G. Martin           *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer nbondbond13flag,ndihedrals,dihedstyle,ndihedtypes
c     --- local variables
      integer i,jtmp,j,index
      double precision tclocal,tmp
      dimension tmp(3)

      write (6,*) 'BondBond13 Coeffs ...'
      if (ndihedrals.eq.0.or.dihedstyle.ne.2) stop
     &     'BondBond13 Coeffs entry not allowed'
      nbondbond13flag = 1

      do i = 1,ndihedtypes
         read (2,*) jtmp,(tmp(j),j=1,3)
c        --- assign the force constant (kcal)
         index = 29
         tclocal = tmp(1)*CON_KCALTOK
         call twh_torcoeff(GLB_SET,i,index,tclocal)
c        --- assign the distances
         index = 30
         tclocal = tmp(2)
         call twh_torcoeff(GLB_SET,i,index,tclocal)
         index = 31
         tclocal = tmp(3)
         call twh_torcoeff(GLB_SET,i,index,tclocal)
      enddo

      return
      end

      subroutine twh_read_angleangle_coeff(nangleangleflag
     &     ,nimpropers,improstyle,nimprotypes,naatypes)
c     ******************************************************************
c     * read_angleangle_coeff                                          *
c     * reads in the angle angle coefficients, unfortunately LAMMPS    *
c     * handles these terms fairly differently than Towhee so it is    *
c     * not so easy to switch between the types                        *
c     * we will stash these variables in the improper torsions until   *
c     * we can figure out what to do with them                         *
c     *                                                                *
c     * rewritten from LAMMPS 2002 03-04-2002 by M.G. Martin           *
c     * last modified 06-14-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer nangleangleflag,nimprotypes,improstyle
      integer nimpropers,naatypes
c     --- local variables
      integer i,jtmp,j
      double precision tmp,onepi
      dimension tmp(6)
      double precision aaclocal
      dimension aaclocal(MINAACOEFF:MAXAACOEFF)

      write (6,*) 'AngleAngle Coeffs ...'
      if (nimpropers.eq.0.or.improstyle.ne.3) stop
     &     'AngleAngle Coeffs entry not allowed'
      nangleangleflag = 1

      call twh_constant_pi(GLB_GET,onepi)

      do i = 1,nimprotypes
         read (2,*) jtmp,(tmp(j),j=1,6)
c        --- make sure we are not going to exceed TAAMAX
         if ( naatypes + 3 .gt. TAAMAX ) then
            write(6,*) 'READLAMMPS: angle-angle types exceeds TAAMAX'
            stop
         endif

c        --- each of these three will be a separate entry in aacoeff
c        --- angles ABC and CBD
         naatypes = naatypes + 1
         call twh_aastyle(GLB_SET,naatypes,2)
         aaclocal(0) = tmp(1)*CON_KCALTOK
         aaclocal(1) = tmp(4)*onepi/180.0d0
         aaclocal(2) = tmp(5)*onepi/180.0d0
         call twh_set_aacoeff(naatypes,aaclocal)
c        ---- angles ABC and ABD
         naatypes = naatypes + 1
         call twh_aastyle(GLB_SET,naatypes,2)
         aaclocal(0) = tmp(2)*CON_KCALTOK
         aaclocal(1) = tmp(4)*onepi/180.0d0
         aaclocal(2) = tmp(6)*onepi/180.0d0
         call twh_set_aacoeff(naatypes,aaclocal)
c        --- angles ABD and CBD
         naatypes = naatypes + 1
         call twh_aastyle(GLB_SET,naatypes,2)
         aaclocal(0) = tmp(3)*CON_KCALTOK
         aaclocal(1) = tmp(6)*onepi/180.0d0
         aaclocal(2) = tmp(5)*onepi/180.0d0
         call twh_set_aacoeff(naatypes,aaclocal)
      enddo

      return
      end

c     ******************************************************************
c     * match                                                          *
c     * this function returns TRUE if str1 matches 1st chars of str2   *
c     *               returns FALSE otherwise                          *
c     * also returns m = loc of next char in str2                      *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-04-2002 by M.G. Martin           *
c     * last modified 05-08-2003 by M.G. Martin                        *
c     ******************************************************************
      function match(str1,str2,m)
      implicit none
c     --- function declarations local to the readlammps subroutine
      logical match
c     --- variables passed to/from this function
      character*(*) str1,str2
      integer m

      match = .FALSE.
      m = len(str1) + 1
      if (len(str1).gt.len(str2)) return
      if (str1.eq.str2(1:len(str1))) match = .TRUE.

      return
      end



c     ******************************************************************
c     * strread                                                        *
c     * return the substring that starts at loc m in str               *
c     * skip initial spaces and tabs                                   *
c     * substring is any chars up to next space or tab                 *
c     * also return m = loc of char after end of substr in str         *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-04-2002 by M.G. Martin           *
c     * last modified 06-17-2005 by M.G. Martin                        *
c     ******************************************************************
      subroutine strread(str,m,substr)
      implicit none
c     --- variables passed to/from the subroutine
      character*(*) str,substr
      integer m
c     --- local variables
      integer n1,n2,nmax
      
      n1 = m
      nmax = len(str)
      do while (n1.le.nmax.and.
     &     (str(n1:n1).eq.' '.or.ichar(str(n1:n1)).eq.9))
         n1 = n1 + 1
      enddo
      n2 = n1
      do while (n2.le.nmax.and.str(n2:n2).ne.' '
     &     .and.ichar(str(n2:n2)).ne.9)
         n2 = n2 + 1
      enddo
      substr = str(n1:n2-1)
      m = n2
      
      return
      end

      subroutine readlammps_input(lfailure,nonstyle,bondstyle
     &     ,lmpangstyle,dihedstyle,improstyle,amberflag,special)
c     ******************************************************************
c     * readlammps_input                                               *
c     * this subroutine reads in some important information from the   *
c     * lammps input file                                              *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-04-2002 by M.G. Martin           *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer nonstyle,bondstyle,lmpangstyle,dihedstyle,improstyle
      integer amberflag
      double precision special
      dimension special(3)
c     --- local variables
      character*200 str
      logical ldone
      integer m,mixflag,mixstyle
      double precision pressure,dielect,temperature
      character*30 filename
c     --- function declarations local to the readlammps subroutine
      logical match

c     --- set the LAMMPS defaults
      nonstyle = 1
      bondstyle = 1
      lmpangstyle = 1
      dihedstyle = 1
      improstyle = 1
      mixflag = 0
c     --- set some towhee defaults
      pressure = 1.0d0
      write(6,*) 'Reading lammps_input'

      filename = 'lammps_input'
      call twh_openfile(lfailure,11,filename,'DEFAULT','OLD','DEFAULT')
      if ( lfailure ) return
      ldone = .false.
      do while (.not. ldone)
        read (11,900,end=999,err=999) str
 900    format (a)
        if (match('mixing style',str,m)) then
           call in_mixing_style(str,m,mixflag,mixstyle)
        elseif (match('nonbond style',str,m)) then
           call in_nonbond_style(str,m,nonstyle,mixflag,mixstyle)
        elseif (match('bond style',str,m)) then
           call in_bond_style(str,m,bondstyle)
        elseif (match('angle style',str,m)) then
           call in_angle_style(str,m,lmpangstyle)
        elseif (match('dihedral style',str,m)) then
           call in_dihedral_style(str,m,dihedstyle)
        elseif (match('improper style',str,m)) then
           call in_improper_style(str,m,improstyle)
        else if (match('special bonds',str,m)) then
          call in_special_bonds(str,m,amberflag,special)
        elseif (match('temp control',str,m)) then
           call in_temp_control(str,m)
        else if (match('create temp',str,m)) then
          call in_create_temp(str,m)
        else if (match('press control',str,m)) then
          call in_press_control(str,m,pressure)
        elseif (match('dielectric',str,m)) then
           call in_dielectric(str,m,dielect)
        endif

      enddo
 999  close(11)
      
c     --- assign the Towhee mixrule
      if (mixstyle .eq. 1) then
         call twh_classical_mixrule(GLB_SET,'Geometric')
      elseif (mixstyle .eq. 2) then
         call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')
      elseif (mixstyle .eq. 3) then
         call twh_classical_mixrule(GLB_SET,'Sixth Power')
      endif
c     --- assign the Towhee temperature
      call twh_temperature(GLB_GET,temperature)
c     --- assign the Towhee express
      pressure = pressure*CON_KPATOSIM
      call twh_express(GLB_SET,pressure)
      return
      end

c     ******************************************************************
c     * in_mixing_style                                                *
c     * reads in the nonbonded mixing style and sets it in Towhee      *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-04-2002 by M.G. Martin           *
c     * last modified 06-17-2005 by M.G. Martin                        *
c     ******************************************************************
      subroutine in_mixing_style(str,m,mixflag,mixstyle)
      implicit none
c     --- function declarations local to the readlammps subroutine
      logical match
c     --- variables passed to/from subroutine
      character*(*) str
      integer m,mixflag,mixstyle
c     --- local variables
      character*16 work
      integer mm

      call strread(str,m,work)
      mixflag = 1
      if (match('geometric',work,mm)) then
        mixstyle = 1
      elseif (match('arithmetic',work,mm)) then
        mixstyle = 2
      elseif (match('sixthpower',work,mm)) then
        mixstyle = 3
      else
        stop 'Bad mixing style parameter'
      endif

      return
      end


c     ******************************************************************
c     * realread                                                       *
c     * return the double precision that starts at loc m in str        *
c     * skip initial spaces and tabs                                   *
c     * scans to end of number marked by space or tab                  *
c     * also return m = loc of char after end of number in str         *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-06-2002 by M.G. Martin           *
c     * last modified 06-17-2005 by M.G. Martin                        *
c     ******************************************************************
      function realread(str,m)
      implicit none
c     --- fucntion declarations local to readlammps
      double precision realread
c     --- variables passed to/from the function
      character*(*) str
      integer m
c     --- local variables
      integer n1,n2,nmax
      
      n1 = m
      nmax = len(str)
      do while (n1.le.nmax.and.
     $     (str(n1:n1).eq.' '.or.ichar(str(n1:n1)).eq.9))
        n1 = n1 + 1
      enddo
      n2 = n1
      do while (n2.le.nmax.and.str(n2:n2).ne.' '
     $     .and.ichar(str(n2:n2)).ne.9)
        n2 = n2 + 1
      enddo
      read (str(n1:n2-1),*) realread
      m = n2

      return
      end


c     ******************************************************************
c     * intread                                                        *
c     * return the integer that starts at loc m in str                 *
c     * skip initial spaces and tabs                                   *
c     * scans to end of number marked by space or tab                  *
c     * also return m = loc of char after end of number in str         *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-06-2002 by M.G. Martin           *
c     * last modified 06-17-2005 by M.G. Martin                        *
c     ******************************************************************
      function intread(str,m)
      implicit none
c     --- function declarations local to readlammps
      integer intread
c     --- variables passed to/from the function
      character*(*) str
      integer m
c     --- local variables
      integer n1,n2,nmax

      n1 = m
      nmax = len(str)
      do while (n1.le.nmax.and.
     $     (str(n1:n1).eq.' '.or.ichar(str(n1:n1)).eq.9))
        n1 = n1 + 1
      enddo
      n2 = n1
      do while (n2.le.nmax.and.str(n2:n2).ne.' '
     $     .and.ichar(str(n2:n2)).ne.9)
        n2 = n2 + 1
      enddo
      read (str(n1:n2-1),*) intread
      m = n2

      return
      end
      
      
c     ******************************************************************
c     * in_nonbond_style                                               *
c     * reads in the LAMMPS nonbond style and uses this to set the     *
c     * Towhee classical_potential also sets LAMMPS mixstyle if not    *
c     * set already                                                    *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-06-2002 by M.G. Martin           *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      subroutine in_nonbond_style(str,m,nonstyle,mixflag,mixstyle)
      implicit none
c     --- function declarations local to readlammps
      logical match
      integer intread
      double precision realread
c     --- variables passed to/from the subroutine
      character*(*) str
      integer m,nonstyle,mixflag,mixstyle
c     --- local variables
      character*16 work
      integer mm,offsetflag
      double precision dummy,rcut,rcutin

      call strread(str,m,work)
      if (match('none',work,mm)) then
         stop 'Towhee cannot use a LAMMPS nonstyle of none'

      elseif (match('lj/cutoff',work,mm)) then
         nonstyle = 1
c        --- set the cutoff
         rcut = realread(str,m)
         call twh_rcut(GLB_SET,rcut)
c        --- set lshift
         offsetflag = intread(str,m)
         if ( offsetflag .eq. 0 ) then
            call twh_lshift(GLB_SET,.false.)
         else
            call twh_lshift(GLB_SET,.true.)
         endif
         if (mixflag .eq. 0) mixstyle = 1
         
      elseif (match('lj/smooth',work,mm)) then
         nonstyle = 2
c        --- set the cutoff
         dummy = realread(str,m)
         rcut = realread(str,m)
         call twh_rcut(GLB_SET,rcut)
         if (mixflag .eq. 0) mixstyle = 1
c        --- set lshift to false
         call twh_lshift(GLB_SET,.false.)
         
      elseif (match('lj/shift',work,mm)) then
         nonstyle = 3
c        --- set the cutoff
         rcut = realread(str,m)
         call twh_rcut(GLB_SET,rcut)
c        --- set lshift
         offsetflag = intread(str,m)
         if ( offsetflag .eq. 0 ) then
            call twh_lshift(GLB_SET,.false.)
         else
            call twh_lshift(GLB_SET,.true.)
         endif
         if (mixflag .eq. 0) mixstyle = 1
         
      elseif (match('soft',work,mm)) then
         nonstyle = 4
         
      elseif (match('class2/cutoff',work,mm)) then
         nonstyle = 5
c        --- set the cutoff
         rcut = realread(str,m)
         call twh_rcut(GLB_SET,rcut)
c        --- set lshift
         offsetflag = intread(str,m)
         if ( offsetflag .eq. 0 ) then
            call twh_lshift(GLB_SET,.false.)
         else
            call twh_lshift(GLB_SET,.true.)
         endif
         if (mixflag .eq. 0) mixstyle = 3
         
      elseif (match('lj/charmm',work,mm)) then
         nonstyle = 6
c        --- set the cutoff
         rcut = realread(str,m)
         call twh_rcut(GLB_SET,rcut)
         dummy = realread(str,m)
         if ( mixflag .eq. 0 ) mixstyle = 2
c        --- set lshift
         call twh_lshift(GLB_SET,.false.)

      else
         stop 'Bad nonbond style parameter'
      endif

c     --- set rcutin
      call twh_rcut(GLB_GET,rcut)
      rcutin = rcut
      call twh_rcutin(GLB_SET,rcutin)
      
      return
      end


c     ******************************************************************
c     * in_bond_style                                                  *
c     * reads in the bond style for all of the bonds.                  *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-06-2002 by M.G. Martin           *
c     * last modified 03-03-2006 by M.G. Martin                        *
c     ******************************************************************
      subroutine in_bond_style(str,m,bondstyle)
      implicit none
c     --- function declarations local to the readlammps subroutine
      logical match
c     --- variables passed to/from the subroutine
      character*(*) str
      integer m,bondstyle
c     --- local variables
      integer mm
      character*16 work

      call strread(str,m,work)
      if (match('none',work,mm)) then
         bondstyle = 0
      elseif (match('harmonic',work,mm)) then
         bondstyle = 1
      elseif (match('fene/standard',work,mm)) then
         bondstyle = 2
      elseif (match('fene/shift',work,mm)) then
         bondstyle = 3
      elseif (match('nonlinear',work,mm)) then
         bondstyle = 4
      elseif (match('class2',work,mm)) then
         bondstyle = 5
      else
         stop 'Bad bond style parameter'
      endif

      return
      end

c     ******************************************************************
c     * in_angle_style                                                 *
c     * reads in the angle style                                       *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-06-2002 by M.G. Martin           *
c     * last modified 03-03-2006 by M.G. Martin                        *
c     ******************************************************************
      subroutine in_angle_style(str,m,lmpangstyle)
      implicit none
c     --- function declarations local to the readlammps subroutine
      logical match
c     --- variables passed to/from the subroutine
      character*(*) str
      integer m,lmpangstyle
c     --- local variables
      integer mm
      character*16 work

      call strread(str,m,work)
      if (match('none',work,mm)) then
        lmpangstyle = 0
      elseif (match('harmonic',work,mm)) then
        lmpangstyle = 1
      elseif (match('class2',work,mm)) then
        lmpangstyle = 2
      elseif (match('charmm',work,mm)) then
        lmpangstyle = 3
      else
        stop 'Bad angle style parameter'
      endif

      return
      end

c     ******************************************************************
c     * in_dihedral_style                                              *
c     * reads in the LAMMPS 2001 dihedral style                        *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-06-2002 by M.G. Martin           *
c     * last modified 03-03-2006 by M.G. Martin                        *
c     ******************************************************************
      subroutine in_dihedral_style(str,m,dihedstyle)
      implicit none
c     --- function declarations local to the readlammps subroutine
      logical match
c     --- variables passed to/from the subroutine
      character*(*) str
      integer m,dihedstyle
c     --- local variables
      integer mm
      character*16 work

      call strread(str,m,work)
      if (match('none',work,mm)) then
        dihedstyle = 0
      elseif (match('harmonic',work,mm)) then
        dihedstyle = 1
      elseif (match('class2',work,mm)) then
        dihedstyle = 2
      elseif (match('multiharmonic',work,mm)) then
        dihedstyle = 3
      elseif (match('charmm',work,mm)) then
        dihedstyle = 4
      else
        stop 'Bad dihedral style parameter'
      endif

      return
      end

c     ******************************************************************
c     * in_improper_style                                              *
c     * reads in the LAMMPS 2001 improper style                        *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-06-2002 by M.G. Martin           *
c     * last modified 03-03-2006 by M.G. Martin                        *
c     ******************************************************************
      subroutine in_improper_style(str,m,improstyle)
      implicit none
c     --- function declarations local to the readlammps subroutine
      logical match
c     --- variables passed to/from the subroutine
      character*(*) str
      integer m,improstyle
c     --- local variables
      integer mm
      character*16 work

      call strread(str,m,work)
      if (match('none',work,mm)) then
        improstyle = 0
      elseif (match('harmonic',work,mm)) then
        improstyle = 1
      elseif (match('cvff',work,mm)) then
        improstyle = 2
      elseif (match('class2',work,mm)) then
        improstyle = 3
      else
        stop 'Bad improper style parameter'
      endif

      return
      end

c     ******************************************************************
c     * in_special_bonds                                               *
c     * reads in the list of special bonds                             *
c     * returns some variables to help determine the 1-4 factors for   *
c     * vdw and coulombic interactions, these are then used later in   *
c     * the torsion parts (read_dihedral_coeff)                        *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-06-2002 by M.G. Martin           *
c     * last modified 06-17-2005 by M.G. Martin                        *
c     ******************************************************************
      subroutine in_special_bonds(str,m,amberflag,special)
      implicit none
c     --- function declarations local to the readlammps subroutine
      logical match
      double precision realread
c     --- variables passed to/from the subroutine
      character*(*) str
      integer m,amberflag
      double precision special
      dimension special(3)
c     --- local variables
      character*16 work
      integer mm

      mm = m
      call strread(str,mm,work)
      if (match('charmm',work,mm)) then
        amberflag = 0
        special(1) = 0.0
        special(2) = 0.0
        special(3) = 0.0
      elseif (match('amber',work,mm)) then
        amberflag = 1
        special(1) = 0.0
        special(2) = 0.0
        special(3) = 1.0d0/1.2d0
      else
        amberflag = 0
        special(1) = realread(str,m)
        special(2) = realread(str,m)
        special(3) = realread(str,m)
        if (special(1) .lt. 0.0.or.special(1) .gt. 1.0.or.
     &       special(2) .lt. 0.0.or.special(2) .gt. 1.0.or.
     &       special(3) .lt. 0.0.or.special(3) .gt. 1.0)
     &       stop 'Bad special bonds parameter'
      endif

      return
      end

c     ******************************************************************
c     * in_temp_control                                                *
c     * extracts the temperature out of the LAMMPS input file          *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-06-2002 by M.G. Martin           *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      subroutine in_temp_control(str,m)
      implicit none
c     --- function declarations local to the readlammps subroutine
      logical match
      double precision realread
c     --- variables passed to/from the subroutine
      character*(*) str
      integer m
c     --- local variables
      integer mm
      double precision dummy,temperature
      character*16 work

      call strread(str,m,work)
      if (match('none',work,mm)) then
c        --- set default temp of 1.0
         temperature = 1.0d0
         call twh_temperature(GLB_SET,temperature)
      elseif (match('rescale',work,mm)) then
         dummy = realread(str,m)
         temperature = realread(str,m)
         call twh_temperature(GLB_SET,temperature)
      elseif (match('replace',work,mm)) then
         dummy = realread(str,m)
         temperature = realread(str,m)
         call twh_temperature(GLB_SET,temperature)
      elseif (match('langevin',work,mm)) then
         dummy = realread(str,m)
         temperature = realread(str,m)
         call twh_temperature(GLB_SET,temperature)
      elseif (match('nose/hoover',work,mm)) then
         dummy = realread(str,m)
         temperature = realread(str,m)
         call twh_temperature(GLB_SET,temperature)
      else
         stop 'Bad temp control parameter'
      endif

      return
      end

      subroutine in_press_control(str,m,pressure)
c     ******************************************************************
c     * in_press_control                                               *
c     * extracts the pressure from lammps                              *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-06-2002 by M.G. Martin           *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- function declarations local to the readlammps subroutine
      logical match
      double precision realread
c     --- variables passed to/from the subroutine
      character*(*) str
      integer m
      double precision pressure,dtemp
c     --- local variables
      integer mm,presscouple
      character*16 work
      double precision dummy

      call strread(str,m,work)

      if (match('none',work,mm)) then
c        --- set default pressure
         pressure = 1.0d0
      elseif (match('nose/hoover',work,mm)) then

        call strread(str,m,work)
        if (match('xyz',work,mm)) then
          presscouple = 0
        elseif (match('xy',work,mm)) then
          presscouple = 1
        elseif (match('yz',work,mm)) then
          presscouple = 2
        elseif (match('xz',work,mm)) then
          presscouple = 3
        elseif (match('aniso',work,mm)) then
          presscouple = 4
        else
          stop 'Bad press control parameter'
        endif

        if (presscouple .eq. 0) then
           dummy = realread(str,m)
           pressure = realread(str,m)
        else
           mm = m
           call strread(str,m,work)
           if (work(1:4) .eq. 'NULL') then
              pressure = 1.0d0
              call strread(str,m,work)
           else
              m = mm
              dummy = realread(str,m)
              pressure = realread(str,m)
           endif
        endif
      else
         stop 'Bad press control parameter'
      endif

c     --- convert from atm to kpa
      pressure = pressure*CON_ATMTOKPA
c     --- convert from kpa to sim
      dtemp = pressure*CON_KPATOSIM
      call twh_express(GLB_SET,dtemp)
      call twh_ensemble(GLB_SET,ENS_NPT)

      return
      end

      subroutine in_dielectric(str,m,dielect)
c     ******************************************************************
c     * in_dielectric                                                  *
c     * reads in the dielectric constant from LAMMPS                   *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-06-2002 by M.G. Martin           *
c     * last modified 06-17-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- fucntion declarations local to readlammps
      double precision realread
c     --- variables passed to/from the subroutine
      character*(*) str
      integer m
      double precision dielect

      dielect = realread(str,m)
      if (dielect .le. 0.0) stop 'Bad dielectric parameter'

      return
      end


      subroutine in_create_temp(str,m)
c     ******************************************************************
c     * in_create_temp                                                 *
c     * extracts the temperature from LAMMPS                           *
c     *                                                                *
c     * rewritten from LAMMPS 2001 03-22-2002 by M.G. Martin           *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- function declarations local to the readlammps subroutine
      logical match
      double precision realread
c     --- variables passed to/from the subroutine
      character*(*) str
      integer m
c     --- local variables
      character*16 work
      integer mm
      double precision temperature

      call strread(str,m,work)
      if (match('uniform',work,mm)) then
         temperature = realread(str,m)
         call twh_temperature(GLB_SET,temperature)
      else if (match('gaussian',work,mm)) then
         temperature = realread(str,m)
         call twh_temperature(GLB_SET,temperature)
      else if (match('velocity',work,mm)) then
c        --- no temperature to extract
      else
         stop 'Bad create temp parameter'
      endif

      return
      end

