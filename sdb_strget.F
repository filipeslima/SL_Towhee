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
      subroutine SDB_STRGET( strout,nstr, lenstr, inline )
c     ******************************************************************
c     * Purpose: get first string (at most nstr characters) out of a   *
c     * an input line "inline"                                         *
c     *                                                                *
c     * originally written: 12-04-2003 by Peter A. Schultz             *
c     * last modified 03-03-2006 by M.G. Martin                        *
c     ******************************************************************
      IMPLICIT NONE
c     
c     Input declaration:
      INTEGER    nstr
      CHARACTER  strout*(*)
c     Output declarations:
      INTEGER    lenstr
      CHARACTER  inline*(*)
c     
c     Local declarations:
      CHARACTER  line*(128)
c     local typing declarations:
      INTEGER    ic1,icn
c     
c     >>>> EXECUTABLE CODE:
c     
c     Default to blank output string of length 0:
      lenstr = 0
      strout = ' '
c     
c     Set full line to blanks:
      line = ' '
c     Safely transfer inline (unknown length, .le.128)
c     into line (known length=128):
      line = inline
c     
c     Determine bounds of first non-blank contiguous string in line:
      call SDB_STRPARS( line, 1,128, ic1,icn )
c     
      if( ic1 .eq. 0 )then
c     String is blank, return ...
         RETURN
      endif
c     
      lenstr = icn - ic1 + 1
      lenstr = MIN( lenstr, nstr )
      icn = ic1 + lenstr - 1
      strout(1:lenstr) = line(ic1:icn)
c     
      RETURN
      END
