c     $Id: crossproduct.F,v 1.4 2006/02/09 20:47:32 mattwycz Exp $
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2005 Marcus G. Martin                       *
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
      subroutine twh_crossproduct( xone,yone,zone,xtwo,ytwo,ztwo,xcross
     &     ,ycross,zcross)
c     ******************************************************************
c     * computes the cross product of Vone * Vtwo and returns Vcross   *
c     * xone: x for vector one (passed in)                             *
c     * yone: y for vector one (passed in)                             *
c     * zone: z for vector one (passed in)                             *
c     * xtwo: x for vector two (passed in)                             *
c     * ytwo: y for vector two (passed in)                             *
c     * ztwo: z for vector two (passed in)                             *
c     * xcross: x for the cross product vector (returned)              *
c     * ycross: y for the cross product vector (returned)              *
c     * zcross: z for the cross product vector (returned)              *
c     *                                                                *
c     * originally written 09-21-2001 by M.G. Martin                   *
c     * last modified 06-16-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none

c     --- variables passed to/from the subroutine
      double precision xone,yone,zone,xtwo,ytwo,ztwo,xcross,ycross
     &     ,zcross

      xcross = yone*ztwo - zone*ytwo
      ycross = zone*xtwo - xone*ztwo
      zcross = xone*ytwo - yone*xtwo

      return
      end
