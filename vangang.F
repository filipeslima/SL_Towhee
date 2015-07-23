#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2000-2012 Marcus G. Martin                       *
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
      function twh_vangang(cordstyle,ibox,imolty,iunit,index,angle
     &     ,xcord,ycord,zcord)
c     ******************************************************************
c     * computes the angle-angle cross term for a given angle-angle    *
c     * type and the single angle types plus the two angles            *
c     *                                                                *
c     * originally written 02-24-2000 by M.G. Martin                   *
c     * last modified 08-30-2012 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_ARCCOS
#define FUNCTION_DISTANCE
#define FUNCTION_GET_AACOEFF
#define FUNCTION_GET_AASTYLE
#define FUNCTION_GET_BENCOEFF
#define FUNCTION_GET_ITAA
#define FUNCTION_VANGANG
#include "functions.h"
c     --- variables passed to/from the function
      integer ibox,imolty,iunit,index,cordstyle
      double precision angle,xcord,ycord,zcord
      dimension angle(2)
      dimension xcord(4),ycord(4),zcord(4)
c     --- local variables
      integer type,atypeone,atypetwo,style
      double precision angone,angtwo,xone,yone,zone,xtwo,ytwo,ztwo
      double precision distone,disttwo

c     --- get the angles
      if ( cordstyle .eq. AA_COORD ) then
c        --- compute angles from the cords
c        --- compute the aa1 - aa0 - aa2 angle
         call twh_mimage(ibox,xcord(2),ycord(2),zcord(2)
     &        ,xcord(1),ycord(1),zcord(1),xone,yone,zone )
         distone = twh_distance(xone,yone,zone)
         call twh_mimage(ibox,xcord(3),ycord(3),zcord(3)
     &        ,xcord(1),ycord(1),zcord(1),xtwo,ytwo,ztwo )
         disttwo = twh_distance(xtwo,ytwo,ztwo)
         angone = twh_arccos( (xone*xtwo + yone*ytwo + zone*ztwo)
     &        /(distone*disttwo) )
c        --- compute the aa1 - aa0 - aa3 angle (xone is the same)
         call twh_mimage(ibox,xcord(4),ycord(4),zcord(4)
     &        ,xcord(1),ycord(1),zcord(1),xtwo,ytwo,ztwo )
         disttwo = twh_distance(xtwo,ytwo,ztwo)
         angtwo = twh_arccos( (xone*xtwo + yone*ytwo + zone*ztwo)
     &        /(distone*disttwo) )
      elseif ( cordstyle .eq. AA_ANGLE ) then
c        --- use the angles that were passed to the function
         angone = twh_arccos(angle(1))
         angtwo = twh_arccos(angle(2))
      else
         write(6,*) 'cordstyle ',cordstyle
         stop 'invalid cordstyle in vangang'
      endif
         
c     --- determine angle-angle type
      type = twh_get_itaa(imolty,iunit,index)
      style = twh_get_aastyle(type)

      if ( style .eq. 1 ) then
c        --- Compass angle-angle cross term
c        --- energy = k* (angle(1) - equil(1))*(angle(2) - equil(2))
         call twh_iataa1(GLB_GET,imolty,iunit,index,atypeone)
         call twh_iataa2(GLB_GET,imolty,iunit,index,atypetwo)
         twh_vangang = twh_get_aacoeff(type,0)
     &        *(angone - twh_get_bencoeff(atypeone,0))
     &        *(angtwo - twh_get_bencoeff(atypetwo,0))
      elseif ( style .eq. 2 ) then
c        --- Compass angle-angle term, but without the clever logic
c        --- of storing the angles in the bencoeffs.  
         twh_vangang = twh_get_aacoeff(type,0)
     &        *(angone - twh_get_aacoeff(type,1))
     &        *(angtwo - twh_get_aacoeff(type,2))
      else
         stop 'non defined angle-angle type'
      endif

      return
      end

