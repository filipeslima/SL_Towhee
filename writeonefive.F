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
      subroutine twh_writeonefive(index)
c     ******************************************************************
c     * writes all of the information about the one-five interactions  *
c     * to the file index                                              *
c     *                                                                *
c     * originally written 12-15-2004 by M.G. Martin                   *
c     * last modified 08-09-2011 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LOFHERE
#define FUNCTION_GET_OFCOEFF
#define FUNCTION_GET_OFNAMES
#define FUNCTION_GET_OFSTYLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer index

c     --- local variables
      logical lfound
      integer type

c     --- see if there are any special one-five types
      lfound = .false.
      do type = 1,TOFMAX
         if ( twh_get_lofhere(type) )  then
            if ( .not. lfound ) then
c              --- this is the first found interaction - write the intro
               write(index,'(a22)') 'Special One-Five Types'
               lfound = .true.
            endif
            if ( twh_get_ofstyle(type) .eq. 1 ) then
               write(index,10) type,twh_get_ofnames(type,1)
     &              ,twh_get_ofnames(type,5)
     &              ,twh_get_ofcoeff(type,1),twh_get_ofcoeff(type,2)
     &              ,twh_get_ofcoeff(type,0)
            elseif ( twh_get_ofstyle(type) .eq. 2 ) then
               write(index,20) type,twh_get_ofnames(type,1)
     &              ,twh_get_ofnames(type,5)
     &              ,twh_get_ofcoeff(type,1)
            else
               write(index,*) 'WRITEONEFIVE: unknown ofstyle '
     &              ,twh_get_ofstyle(type),' for type ',type
               stop
            endif
         endif
      enddo

 10   format('Type:',1x,i3,1x,'Atom(i):',1x,a10,1x,'Atom(j):',1x,a10,1x
     &     ,'Sigma:',1x,f7.3,1x,'Epsilon:',1x,f7.3,1x,'Shift:',1x,f7.4)
 20   format('Type:',1x,i3,1x,'Atom(i):',1x,a10,1x,'Atom(j):',1x,a10,1x
     &     ,'A/r^12:',1x,e10.4)

      return
      end
    
