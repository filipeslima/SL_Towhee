c     ******************************************************************
c     * SDB: Utilities for handling the sdb database format            *
c     * Copyright (C) 2004-2006 Marcus G. Martin, Peter A. Schultz     *
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
      subroutine SDB_RDSTRUCT( ierr, IDBFILE, IWR,
     $     edata,eatom_eV,ebind_eV,etotal_Ry,
     $     ndim,natm,ntyp, itypa, natmnm,atmnm,typnm, rcell,ratom,
     $     nnote,notes, structnm )
c     ******************************************************************
c     * Purpose: parse training set structure data                     *
c     * Default distance units are Angstrom                            *
c     * Default energy units (Ebinding, Eatomize) are eV, entire cell  *
c     *                                                                *
c     * originally written: Peter A. Schultz, January-2004             *
c     * last modified: 03-03-2006 by M.G. Martin                       *
c     ******************************************************************
      IMPLICIT NONE
c     
c     Input:
c     IWR = user output unit (output supressed if IWR<1)
c     IDBFILE = file unit for data read
      INTEGER    IWR, IDBFILE, IWRX
c     natmnm = length of atom/type names
      INTEGER    natmnm
c
c     Output:
c     ierr = error flag:
c     0 = completed ok
c     <0 = EOF discovered
c     >0 = internal error
      INTEGER    ierr, jerr
c     edata = string recording type of energes obtained
      CHARACTER  edata*(3)
c     eatom_eV  = atomization energy in eV; valid if edata(1:1).ne.blank
c     ebind_eV  = binding energy in eV    ; valid if edata(2:2).ne.blank
c     etotal_Ry = total energy in Ry      ; valid if edata(3:3).ne.blank
      double precision     eatom_eV, ebind_eV, etotal_Ry
c     ndim=dimensionality of system; natm=no. of atoms; ntyp=no. of types
      INTEGER    ndim,natm,ntyp,  idim,iatm
c     rcell(3,3) = unit cell vectors (Cartesian vectors in Angstrom)
c     ratom(3,natm) = atom positions (Cartesian vectors in Angstrom)
      double precision     rcell, ratom
      DIMENSION  rcell(3,3), ratom(3,*)
c     itypa(natm) = the type of atom for each atom
      INTEGER    itypa
      DIMENSION  itypa(*)
c     atmnm(natm)=name for each atom; typnm(natm)=name for each type
      CHARACTER*(*)  atmnm(*),typnm(*)
c     nnote=number of note lines; notes(9) = 80-char note lines
      INTEGER    nnote, inote
      CHARACTER*(80)  notes(*)
      CHARACTER  structnm*(40)
c     
c     Local declarations
c     sectname = local module data name
      CHARACTER*(8)  sectname, endname
      DATA           sectname / 'structur' /
c     Characters denoting tagline, local and external
      CHARACTER  taglocal*1,tagxtrnl*1
      DATA       taglocal,  tagxtrnl   / '$','@' /
c     
      CHARACTER  line*(128), line1*(1)
      EQUIVALENCE  (line,line1)
      CHARACTER  tagname*8,string8*8
c     
      LOGICAL  printon
      INTEGER  ncellrd,natmrd,ntyprd
      INTEGER  i,iline
      INTEGER  lenstr
      INTEGER  ic1,icn,nch, iu1,iun
c     
c     Distance unit conversion:
      double precision  units_r
c     Energy unit conversion:
      double precision  units_e
c     
      double precision  one

c     
c     >>>> EXECUTABLE CODE:
c
      one = 1.0d0
      ierr = 0
c     Set default for print diagnostics:
      printon = .false.
c     And turn off prints if unit is invalid:
      if( IWR.lt.1 )  printon = .false.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     Initialize some defaults
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     Nonsense energy values to trigger problems if inadvertantly used:
      eatom_eV = 9.d99
      ebind_eV = 9.d99
      etotal_Ry = 9.d99
c     Initialize/default defining parameters:
      nnote = 0
      ndim = 3
      ncellrd = 0
      natm = 0
      natmrd = 0
      ntyp = 0
      ntyprd = 0
      edata = '   '
c     Default is Angstrom units:
      units_r = one
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                        Look for next tag
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      iline = 0
 100  continue
c     
      read(IDBFILE,8128,end=913,err=913)  line
      if( printon ) write(IWR,8128)  line
 8128 format(a128)
 8080 format(a80)
      iline = iline + 1
c     
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     On first line, validate type of dataset and get its name
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if( iline.eq.1 )then
c     
c     Check that this is the right kind of data set
         if( line1 .eq. tagxtrnl )then
c     First, check to skip module-specifying tag for this:
c     
            call SDB_STRGET( tagname,8, lenstr, line(2:128) )
            if( tagname .ne. sectname )then
c     The @module designated does not match this module
               ierr = 1
               goto 913
            endif
c     
            iline = 0
            goto 100
         endif
c     
         if( line1 .eq. taglocal )then
c     
            call SDB_STRGET( tagname,8, lenstr, line(2:128) )
            if( tagname .eq. 'start   ' )then
c
               call SDB_STRGET( string8,8, lenstr, line(8:128) )
               if( string8 .eq. sectname )then
c     Correct module for this dataset, get dataset name
                  read(IDBFILE,8128)  line
                  call SDB_STRPARS( line, 1,128, ic1,icn )
                  structnm = ' '
                  if( ic1 .gt. 0 )then
                     nch = icn - ic1 + 1
                     if( nch.gt.40 )then
                        nch = 40
                        icn = ic1 - 1 + 40
                     endif
                     structnm(1:nch) = line(ic1:icn)
                  endif
                  goto 100
c     
               endif
            endif
c     
c     End first line check
         endif
c     
c     We do not have the correct routine for this dataset
         ierr = 2
         goto 913
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     Find and process tags
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if( line1 .eq. tagxtrnl )then
c
c       >>>>> Indicates external module to be processed
c
         call SDB_STRGET( tagname,8, lenstr, line(2:128) )
c     
c     ... and check if we know where to access this module:
         if( tagname .eq. '        ' )then
c     Blank module name not allowed
            ierr = 3
            goto 913
         else
c     Do not recognize module, skip to its end:
            call SDB_RDTOEND( jerr, IDBFILE, IWR, tagname )
            if( jerr.ne.0 )then
               write(6,*) 'SDB_RDSTRUCT: error'
               write(6,*) line1
               ierr = 4
               goto 913
            endif
         endif
         goto 100
c     
c     <<<<< End module tag
c     
      elseif( line1 .eq. taglocal )then
c     
c     >>>>> Indicates local tag to be processed
c     
         call SDB_STRGET( tagname,8, lenstr, line(2:128) )
c     
         if( tagname .eq. 'end     ' )then
c     We have found an "end", verify that it is the correct end
            call SDB_STRGET( endname,8, lenstr, line(5:128) )
            if( endname .eq. sectname )then
c     We are done with this section of data
               goto 900
            else
c     Oops, something is badly wrong, bad section end name
               ierr = 5
               goto 913
            endif
c     
         elseif( tagname(1:5) .eq. 'notes' )then
c     Read the note lines
            nnote = 0
            if( tagname(6:6) .ne. ' ' ) read(tagname(6:6),'(i1)')  nnote
            if( nnote.eq.0 ) nnote = 1
            do  inote=1,nnote
               read(IDBFILE,8080)  notes(inote)
            enddo
c     
         elseif( tagname .eq. 'distance' )then
c     Read the distance unit used in the file
c     
            read(IDBFILE,8080)  line
            call SDB_STRPARS( line,1,128,ic1,icn )
            if( ic1 .lt. 1 )then
c     Cannot have blank string here.
               ierr = 6
               goto 913
            endif
            call SDB_RTOANG( line(ic1:icn), units_r )
c     
         elseif( tagname .eq. 'dimensio' )then
c     Read the system dimensionality
            read(IDBFILE,*)  ndim
            if( ndim.lt.0 .or. ndim.gt.3 )then
c     If given, must be valid entry
               ierr = 9
               goto 913
            endif
c     
         elseif( tagname .eq. 'cell    ' )then
c     Read the primitive cell vectors
            do  idim=1,3
               read(IDBFILE,*)  (rcell(i,idim),i=1,3)
            enddo
            ncellrd = 1
c     
         elseif( tagname .eq. 'number  ' )then
c     Read number of atoms
            read(IDBFILE,*)  natm
            if( natm.lt.0 )then
               ierr = 10
               goto 913
            endif
c     
         elseif( tagname .eq. 'atom    ' )then
c     Read atoms, their types, and position vectors:
            IWRX = -1
            if( printon ) IWRX = IWR
            call SDB_CONFIGRD( IDBFILE,IWRX, natmnm, ntyp,typnm,
     $           natm, atmnm, itypa, ratom )
            natmrd = natm
            ntyprd = ntyp
c     
         elseif( tagname .eq. 'Eatomize' )then
c     Read binding energy (eV) w.r.t. atoms, and E-unit (if given)
            edata(1:1) = 'A'
            read(IDBFILE,8128)  line
            read(line,*)  eatom_eV
c     Fetch energy unit and convert to eV if it is there ...
            call SDB_STRPARS( line,1,128,ic1,icn )
            call SDB_STRPARS( line,icn+1,128, iu1,iun )
            if( iu1.gt.0 )then
               call SDB_ETOEV( line(iu1:iun), units_e )
               eatom_eV = eatom_eV * units_e
            endif
c     
         elseif( tagname .eq. 'Ebinding' )then
c     Read binding energy (eV) w.r.t. standard state, and E-unit
            edata(2:2) = 'B'
            read(IDBFILE,8128)  line
            read(line,*)  ebind_eV
c     Fetch energy unit and convert to eV if it is there ...
            call SDB_STRPARS( line,1,128,ic1,icn )
            call SDB_STRPARS( line,icn+1,128, iu1,iun )
            if( iu1.gt.0 )then
               call SDB_ETOEV( line(iu1:iun), units_e )
               ebind_eV = ebind_eV * units_e
            endif
c     
         elseif( tagname .eq. 'Etotal  ' )then
c     Read total energy (Ry), raw total energy from calculation
            read(IDBFILE,*)  etotal_Ry
            edata(3:3) = 'T'
c     
         else
c     Unknown tag, skip it and look for next tag:
            goto 100
         endif
c     
c     <<<<< End local tags
      else
c     Not a local or xtrnl tag, skip it and look for next tag
         goto 100
      endif
c     
c     Go look for next tag:
      goto 100
c     
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                  Local post-processing of data
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c    Check for data sufficiency, and convert data to standard form
c
 900  continue
      if( ncellrd.eq.0 .or. natmrd.eq.0 .or. ntyprd.eq.0 )then
c     Error: did not pick up a unit cell, atom positions, or types
         ierr = 11
         goto 913
      endif
c     
c  Convert input distances into appropriate units
c
      if( ncellrd .ne. 0 )then
c     Supercell vectors
         do  idim=1,3
            do i=1,3
               rcell(i,idim) = rcell(i,idim) * units_r
            enddo
         enddo
      endif
c     
      if( natmrd .ne. 0 ) then
c     Presuming that atom coordinates are Cartesian
c     No provision (yet) for atoms in lattice/direct units
         do  iatm=1,natm
            do i=1,3
               ratom(i,iatm) = ratom(i,iatm) * units_r
            enddo
         enddo
      endif
c     
c  Go to routine exit
      goto 999
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                       Error conditions
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 913  continue
      if( ierr.eq.0 ) ierr = -1
      goto 999
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                       Exit from routine
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 999  continue
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SDB_ETOEV
c
c
      subroutine SDB_ETOEV( u_name, u_convert )
c---------------------------------------------------------------
c Purpose: energy conversion to eV
c
c Written: Peter A. Schultz, January-2004
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Returns in u_convert the value in eV of the energy unit that is
c    input in u_name.  Use minimum number of characters to
c    identify the unit, in all sensible upper/lower case combos.
c
c Values taken from:
c    CRC Handbook of Chemistry and Physics, 72nd edit. (1991-1992)
c    David R. Lide, Editor-in-Chief
c
      IMPLICIT NONE
c
      CHARACTER  u_name*(*)
      double precision     u_convert
c
      CHARACTER  ustr*128
      double precision  u_target, u_value
c Energy units: e_xx per 1 eV, ei_xx in eV/xx
      double precision  e_eV, e_meV, ei_Ry, ei_H, e_kcalm, e_kJm, e_cm
     &     , ei_kelv
      DATA  e_eV    / 1. d0 /
      DATA  e_meV   / 1. d3 /
      DATA  ei_Ry   / 13.605698 d0 /
      DATA  ei_H    / 27.211396 d0 /
      DATA  e_kcalm / 23.0605 d0 /
      DATA  e_kJm   / 96.4853 d0 /
      DATA  e_cm    / 8065.54 d0 /
      DATA  ei_kelv / 1.16045d-04 /
c
      double precision  one
      DATA  one / 1.d0 /
c
c >>>> EXECUTABLE CODE:
c
c Specify target unit:
      u_target = e_eV
c Transfer into string of guaranteed adequate length
      ustr = ' '
      ustr = u_name
      if(     ustr(1:2) .eq. 'eV'    .or. ustr(1:2) .eq. 'EV'    )then
        u_value = e_eV
      elseif( ustr(1:3) .eq. 'meV'   .or. ustr(1:2) .eq. 'MEV'   )then
        u_value = e_meV
      elseif( ustr(1:2) .eq. 'Ry'    .or. ustr(1:2) .eq. 'RY'    .or.
     $        ustr(1:2) .eq. 'ry'     )then
        u_value = one / ei_Ry
      elseif( ustr(1:1) .eq. 'h'     .or. ustr(1:1) .eq. 'H'     )then
        u_value = one / ei_H
      elseif( ustr(1:4) .eq. 'kcal'  .or. ustr(1:4) .eq. 'KCAL'  )then
        u_value = e_kcalm
      elseif( ustr(1:2) .eq. 'kJ'    .or. ustr(1:2) .eq. 'KJ'    .or.
     $        ustr(1:2) .eq. 'kj'    )then
        u_value = e_kJm
      elseif( ustr(1:2) .eq. 'cm'    .or. ustr(1:2) .eq. 'CM'    )then
        u_value = e_cm
      elseif( ustr(1:2) .eq. 'ke'    .or. ustr(1:2) .eq. 'KE'    .or.
     $        ustr(1:2) .eq. 'Ke'    )then
        u_value = one / ei_kelv
      else
c       Failed to match a unit - leave the data unchanged
        u_target = one
        u_value = one
      endif
c
      u_convert = u_target / u_value
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SDB_RTOANG
c
c
      subroutine SDB_RTOANG( u_name, u_convert )
c---------------------------------------------------------------
c Purpose: distance conversion to Angstrom
c
c Written: Peter A. Schultz, January-2004
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Returns in u_convert the value in A of the distance unit that is
c    input in u_name.  Use minimum number of characters to
c    identify the unit, in all sensible upper/lower case combos.
c
c Values taken from: memory
c
      IMPLICIT NONE
c
      CHARACTER  u_name*(*)
      double precision u_convert
c
      CHARACTER  ustr*128
      double precision  u_target, u_value
c Distance units, r_xx per 1 Angstrom
      double precision  r_A, r_bohr, r_nm
      parameter (r_A=1.0d0,r_bohr=1.8897265d0,r_nm=0.1d0)
c
      double precision  one
      data    one / 1.0d0 /
c
c >>>> EXECUTABLE CODE:
c
c Specify target unit:
      u_target = r_A
c
c Transfer into string of guaranteed adequate length
      ustr = ' '
      ustr = u_name
      if(     ustr(1:1) .eq. 'B'     .or. ustr(1:1) .eq. 'b'     .or.
     $        ustr(1:2) .eq. 'AU'    .or. ustr(1:2) .eq. 'au'    )then
        u_value = r_bohr
      elseif( ustr(1:1) .eq. 'A'     .or. ustr(1:3) .eq. 'ang'   )then
c       Angstrom units (must check for "A" after check "AU" for bohr)
        u_value = r_A
      elseif( ustr(1:2) .eq. 'nm'    .or. ustr(1:2) .eq. 'NM'    )then
        u_value = r_nm
      else
c       Failed to match a unit - leave the data unchanged
c       Other (default) option is Angstrom, this is not
        u_target = one
        u_value = one
      endif
c
      u_convert = u_target / u_value
c
      RETURN
      END
c
c
