c     $Id: sdb_rdtoend.F,v 1.3 2006/02/09 20:47:38 mattwycz Exp $
c     ******************************************************************
c     * SDB: Utilities for handling the sdb database format            *
c     * Copyright (C) 2004 Marcus G. Martin, Peter A. Schultz          *
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
      subroutine SDB_RDTOEND( ierr, IDBFILE, IWR, sectname )
c     ******************************************************************
c     * Purpose: read database file thru end of designated section     *
c     *                                                                *
c     * originally written 12-05-2003 by Peter A. Schultz              *
c     * last modified 01-27-2004 by M.G. Martin                        *
c     ******************************************************************
      IMPLICIT NONE
c     Input:
c     IWR = user output unit (output supressed if IWR<1)
c     IDBFILE = file unit for data read
      INTEGER    IWR, IDBFILE
c     sectname = name of section to be skipped thru to end
      CHARACTER  sectname*8
c
c     Output:
c     ierr = error flag; if ok returns 0, if error return non-0
      INTEGER    ierr
c     
c     Local declarations
c     character denoting tagline, local and external, and the tag
      CHARACTER  taglocal*1
      DATA       taglocal  / '$' /
c     
      CHARACTER  line*(128), line1*1
      EQUIVALENCE  (line,line1)
      CHARACTER  tagname*8,endname*8
c     
      LOGICAL  printon
      INTEGER  lenstr
c     
c     >>>> EXECUTABLE CODE:
c
      ierr = 0
      printon = .false.
      if( IWR.gt.0 )  printon = .true.
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                        Look for next tag
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c We are looking for the "$end <sectname>" that concludes section
c
 100  continue
c
      read(IDBFILE,9128,end=913,err=913)  line
 9128 format(a128)
      
      if( line1 .eq. taglocal )then
c
c        >>>>> Identify local tag to be parsed
c
         call SDB_STRGET( tagname,8, lenstr, line(2:128) )
         if( tagname .eq. 'end     ' )then
c           We have found an "end", check that the correct end
            call SDB_STRGET( endname,8, lenstr, line(5:128) )
            if( endname .eq. sectname )then
c              We found it!
               goto 999
            endif
         endif
      endif
c     Keep looking
      goto 100
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                       Error conditions
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 913  continue
      if( ierr.eq.0 ) ierr = 1
      goto 999
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                       Exit from routine
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 999  continue
c
      RETURN
      END
