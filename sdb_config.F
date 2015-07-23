c     ******************************************************************
c     * SDB: Utilities for handling the sdb database format            *
c     * Copyright (C) 2004-2011 Marcus G. Martin, Peter A. Schultz     *
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
      subroutine SDB_CONFIGRD( IDT,IWR, natmnm, ntyp,typnm,
     $ natm, atmnm, itypa, ratm )
c     ******************************************************************
c     * Purpose: input atomic configuration                            *
c     *                                                                *
c     * originally written: Peter A. Schultz, 26-November-2001         *
c     * last modified 11-04-2004 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer IDT,IWR,natmnm,ntyp,natm,itypa
      DIMENSION  itypa(*)
      double precision ratm
      DIMENSION  ratm(3,natm)
c     input: (type names)
      CHARACTER*(*)  typnm
      DIMENSION  typnm(*)
c     output: atom names, atom types, atom coordinates
      CHARACTER*(*)  atmnm
      DIMENSION  atmnm(*)
c     --- local variables 
      integer iatm

      do iatm = 1,natm
         call SDB_ATOMRD( IDT,IWR, natmnm, ntyp, typnm,
     $        atmnm(iatm), iatm,itypa, ratm(1,iatm) )
      enddo

      return
      end

      subroutine SDB_ATOMRD(IDT, IWR, natmnm, ntyp, typnm, atmnm, iatm
     &     ,itypa, ratm )
c     ******************************************************************
c     * Purpose: input atom coordinates                                *
c     *                                                                *
c     * originally written: Peter A. Schultz, 26-November-2001         *
c     * last modified 11-08-2006 by M.G. Martin                        *
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c     input: (type names)
      CHARACTER*(*)  typnm
      DIMENSION  typnm(*)
c     output: atom names, atom types, atom coordinates
      CHARACTER*(*)  atmnm
      DIMENSION  itypa(*)
      DIMENSION  ratm(3)
c     local declarations:
      PARAMETER  ( nline=80 )
      CHARACTER  line*(nline), strtmp*(nline)
      LOGICAL    printon

      printon = ( IWR .ge. 0 )
c
c     Read the line with the atomic information:
c     
      read(IDT,'(a80)',end=1301,err=1301 )  line
      if( printon ) write(IWR,'(a)' )   line
      nl1 = 1
      nl2 = nline
c     
c     Parse the atom number/name (first string on line):
c     
      call SDB_STRPARS( line, nl1,nl2, n1,n2 )
      if( n1 .eq. 0 ) STOP 'atm-name'
      nmlen = n2 - n1 + 1
      if( nmlen .gt. natmnm ) then
         write(6,*) 'nmlen,natmnm',nmlen,natmnm
         write(6,*) line
         STOP 'longatom'
      endif
      strtmp = ' '
      strtmp(1:nmlen) = line(n1:n2)
      atmnm = strtmp(1:natmnm)
      nl1 = n2 + 1
c     
c     First, check if atoms input as number sequence.
c     If they are, check that atoms are indeed in sequence.
      read( atmnm, * , err=120 ) jatm
      if( jatm .ne. iatm )  STOP 'atom-seq'
c     Shift atom name to right aligned 6-char deep
      atmnm = ' '
      write( atmnm, '(i6)' ) iatm
c     
c     Parse the atom type:
c     
 120  continue
      call SDB_STRPARS( line, nl1,nl2, n1,n2 )
      if( n1 .eq. 0 )  STOP 'no type'
      nmlen = n2 - n1 + 1
      if( nmlen .gt. natmnm ) STOP 'longtype'
      nl1 = n2 + 1
c     
c     And now try and get the atom type:
c     
      itypa(iatm) = 0
      do  ityp=1,ntyp
         strtmp = typnm(ityp)
         call SDB_STRPARS( strtmp, 1,nline, nt1,nt2 )
         ntlen = nt2 - nt1 + 1
         if( ntlen .eq. nmlen .and.
     $        line(n1:n2) .eq. strtmp(nt1:nt2) )then
            itypa(iatm) = ityp
            goto 200
         endif
      enddo

c     Declare a new type, and name it
c     **** WARNING **** ntyp dim not checked to be adequate
      ntyp = ntyp + 1
      if( IWR.gt.0 ) write(IWR,9125) iatm,ntyp,line(n1:n2)
 9125 format(1x,'>>>> Atom#',i4,' - new type#',i2,' named >',a,'<')
      typnm(ntyp) = line(n1:n2)
      itypa(iatm) = ntyp
      
 200  continue
c     
c     Get the atom coordinates:
      read( line(nl1:nl2), * )  ratm
      
      
      RETURN
 1301 STOP 'err-ratm'
      end
