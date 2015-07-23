#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2011 Marcus G. Martin                       *
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
      subroutine twh_mimage ( ibox,rxui,ryui,rzui,rxuj,ryuj,rzuj
     &     ,rxuij,ryuij,rzuij )
c     ******************************************************************
c     * this subroutine determines the minimum distance between any two*
c     * vectors of the coordinates through all of the periodic images  *
c     * values of arguments rxuij,ryuij,rzuij are set here             *
c     *                                                                *
c     * completely rewritten 02-14-2001 by M.G. Martin                 *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_GT
#define FUNCTION_CMP_LT
#define FUNCTION_MINBOXLENGTH
#include "functions.h"

c     --- variables passed to/from the subroutine
      integer ibox
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj
      double precision rxuij,ryuij,rzuij
c     --- local variables
      logical lrectbox
      integer idim,iii,jjj,kkk
      double precision arbcord
      dimension arbcord(3)
      double precision hbx,xtemp,ytemp,ztemp
      double precision rijsq,rxuij_temp,ryuij_temp,rzuij_temp
      double precision rijsq_temp,hmvalue,minhbxsq

c      write(6,*) 'start MIMAGE',ibox

c     --- compute the distances
      rxuij = rxui - rxuj
      ryuij = ryui - ryuj
      rzuij = rzui - rzuj

c     --- use periodic boundaries unless this is the ideal gas box
      if ( ibox .ne. 0 ) then
         call twh_lrectbox(GLB_GET,ibox,lrectbox)
         if ( lrectbox ) then
c           --- box is a rectangular parallapiped, no need for fancy
c           --- hmatrix operations
c           --- find the minimum image
            idim = 1
            call twh_hmatrix(GLB_GET,ibox,idim,idim,hmvalue)
            call twh_halfbox(GLB_GET,ibox,idim,hbx)
            do while ( twh_cmp_gt(rxuij, hbx) )
               rxuij = rxuij - hmvalue
            enddo
            do while ( twh_cmp_lt(rxuij, -hbx) )
               rxuij = rxuij + hmvalue
            enddo
            idim = 2
            call twh_hmatrix(GLB_GET,ibox,idim,idim,hmvalue)
            call twh_halfbox(GLB_GET,ibox,idim,hbx)
            do while ( twh_cmp_gt(ryuij, hbx) )
               ryuij = ryuij - hmvalue
            enddo
            do while ( twh_cmp_lt(ryuij, -hbx) )
               ryuij = ryuij + hmvalue
            enddo
            idim = 3
            call twh_hmatrix(GLB_GET,ibox,idim,idim,hmvalue)
            call twh_halfbox(GLB_GET,ibox,idim,hbx)
            do while ( twh_cmp_gt(rzuij, hbx) )
               rzuij = rzuij - hmvalue
            enddo
            do while ( twh_cmp_lt(rzuij, -hbx) )
               rzuij = rzuij + hmvalue
            enddo
         else
c           --- non-rectangular simulation cell
c           --- shift the coordinate system
            call twh_carttoarb(ibox,rxuij,ryuij,rzuij
     &           ,arbcord(1),arbcord(2),arbcord(3))

c     --- this method was used alone prior to version 5.3.0
c     --- it usually finds the minimum image, but it's not guaranteed
c     --- try it first, test, and resort to brute force if needed
c           --- find the 'best guess' minimum image
c           --- put arbitrary coords in [-0.5,0.5]
            do idim = 1,3
               do while ( twh_cmp_gt(arbcord(idim), 0.5d0) )
                  arbcord(idim) = arbcord(idim) - 1.0d0
               enddo
               do while ( twh_cmp_lt(arbcord(idim), -0.5d0) )
                  arbcord(idim) = arbcord(idim) + 1.0d0
               enddo
            enddo
c           --- transform back to cartesian coordinates
            call twh_arbtocart(ibox,arbcord(1),arbcord(2),arbcord(3)
     &           ,rxuij,ryuij,rzuij)

c           --- if minimum image vector length is less than half the
c           --- distance between closest walls of simulation cell
c           --- we must have the correct minimum image
            minhbxsq = twh_minboxlength(ibox) / 2.0d0
            minhbxsq = minhbxsq * minhbxsq
            rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij

            if ( twh_cmp_lt(rijsq,minhbxsq) ) return

c           --- otherwise, try it the hard way to be certain

c           --- put arbitrary coords in [0,1]
            do idim = 1,3
               do while ( twh_cmp_gt(arbcord(idim), 1.0d0) )
                  arbcord(idim) = arbcord(idim) - 1.0d0
               enddo
               do while ( twh_cmp_lt(arbcord(idim), 0.0d0) )
                  arbcord(idim) = arbcord(idim) + 1.0d0
               enddo
            enddo

c           --- check for nearest image in eight cell images
c           --- stopping if image within minhbx is found

            do iii = 0,1
               xtemp = arbcord(1) - dble(iii)
               do jjj = 0,1
                  ytemp = arbcord(2) - dble(jjj)
                  do kkk = 0,1
                     ztemp = arbcord(3) - dble(kkk)
                     call twh_arbtocart(ibox,xtemp,ytemp,ztemp
     &                    ,rxuij_temp,ryuij_temp,rzuij_temp)
                     rijsq_temp = rxuij_temp*rxuij_temp +
     &                    ryuij_temp*ryuij_temp + rzuij_temp*rzuij_temp
                     if ( rijsq_temp .lt. rijsq ) then
                        rijsq = rijsq_temp
                        rxuij = rxuij_temp
                        ryuij = ryuij_temp
                        rzuij = rzuij_temp
                    if ( rijsq .lt. minhbxsq ) return
                     endif
                  enddo
               enddo
            enddo
         endif
      endif

c      write(6,*) 'end MIMAGE'

      return
      end
