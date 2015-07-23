c     $Id: sdb_strpars.F,v 1.4 2006/02/09 20:47:38 mattwycz Exp $
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
      subroutine SDB_STRPARS( line, nl1,nl2, n1,n2 )
c     ******************************************************************
c     * Purpose: parse bounds of first character string within line    *
c     *                                                                *
c     * originally written 11-26-2001 by Peter A. Schultz              *
c     * last modified 11-04-2004 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      CHARACTER*(*)  line
      integer nl1,nl2,n1,n2
c     --- local variables
      CHARACTER*(1)  blank
      DATA           blank / ' ' /
      integer il,jl
c
c     >>>> EXECUTABLE CODE:
c
      n1 = 0
      n2 = nl2
      do  il=nl1,nl2
         if( line(il:il) .ne. blank )then
c           Found the first non-blank character ...
            n1 = il
            n2 = n1
c           ... and find the last character in string:
            do  jl=n1+1,nl2
               if( line(jl:jl) .eq. blank ) goto 999
               n2 = jl
            enddo
            goto 999
         endif
      enddo
     
 999  continue

      RETURN
      END
