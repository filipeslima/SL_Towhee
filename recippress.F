#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2011 Marcus G. Martin                       *
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
      subroutine twh_recippress(ibox,repress)
c     ******************************************************************
c     * calculates the reciprocal space contribution to pressure using *
c     * thermodynamic definition. See J. Chem. Phys. Vol. 109 P2791.   *
c     *                                                                *
c     * originally written in 1998 by Bin Chen                         *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NBOXI
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ibox
      double precision repress
      dimension repress(6)

c     --- local variables
      integer ncount,ichain,ii,imolty,index,nchain,nunit,numvect
      double precision factor,piix,piiy,piiz,xcmi,ycmi,zcmi,arg
      double precision inner,rxui,ryui,rzui,charge
      double precision kx,ky,kz,prefact,calp,ssumr,ssumi

      do index = 1,6
         repress(index) = 0.0d0
      enddo

      call twh_ewald_calp(GLB_GET,ibox,calp)
      call twh_ewald_numvect(GLB_GET,ibox,numvect)
      do ncount = 1,numvect
         call twh_ewald_kvector(GLB_GET,CRD_REAL,ncount,ibox,kx,ky,kz)
         call twh_ewald_prefact(GLB_GET,CRD_REAL,ncount,ibox,prefact)
         call twh_ewald_ssumr(GLB_GET,CRD_REAL,ncount,ibox,ssumr)
         call twh_ewald_ssumi(GLB_GET,CRD_REAL,ncount,ibox,ssumi)
         factor = prefact*(ssumr*ssumr + ssumi*ssumi)
         inner = (1.0d0/(4.0d0*calp*calp)
     &        + 1.0d0/(kx*kx + ky*ky + kz*kz))*2.0d0
         repress(1) = repress(1) + factor*(1.0d0 - inner*kx*kx)
         repress(2) = repress(2) + factor*(1.0d0 - inner*ky*ky)
         repress(3) = repress(3) + factor*(1.0d0 - inner*kz*kz)
         repress(4) = repress(4) + factor*(1.0d0 - inner*kx*ky)
         repress(5) = repress(5) + factor*(1.0d0 - inner*kx*kz)
         repress(6) = repress(6) + factor*(1.0d0 - inner*ky*kz)
      enddo
c     --- the intramolecular part should be substracted
      call twh_nchain(GLB_GET,nchain)
      do ichain = 1, nchain
c        --- check if i is in relevant box 
         if ( twh_get_nboxi(ichain) .eq. ibox ) then
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_com(GLB_GET,CRD_REAL,ichain,xcmi,ycmi,zcmi)
c           --- loop over all beads ii of chain i 
            call twh_nunit(GLB_GET,imolty,nunit)
            do ii = 1, nunit
c              --- compute the vector of the bead to the COM (p)
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,ii
     &              ,rxui,ryui,rzui)
               piix = rxui - xcmi
               piiy = ryui - ycmi
               piiz = rzui - zcmi
             
               call twh_qqatom(GLB_GET,imolty,ii,charge)
               do ncount = 1,numvect
                  call twh_ewald_kvector(GLB_GET,CRD_REAL,ncount,ibox
     &                 ,kx,ky,kz)
                  call twh_ewald_prefact(GLB_GET,CRD_REAL,ncount,ibox
     &                 ,prefact)
                  call twh_ewald_ssumr(GLB_GET,CRD_REAL,ncount,ibox
     &                 ,ssumr)
                  call twh_ewald_ssumi(GLB_GET,CRD_REAL,ncount,ibox
     &                 ,ssumi)
c                 --- compute the dot product of k and r
                  arg = kx*rxui + ky*ryui + kz*rzui
                  factor = prefact*2.0d0
     &                 *(-ssumr*dsin(arg)+ssumi*dcos(arg))*charge
                  repress(1) = repress(1) + factor*kx*piix
                  repress(2) = repress(2) + factor*ky*piiy
                  repress(3) = repress(3) + factor*kz*piiz
c                 --- take the average of the two possible combinations
                  repress(4) = repress(4)
     &                 + factor*0.5d0*(kx*piiy+ky*piix)
                  repress(5) = repress(5)
     &                 + factor*0.5d0*(kx*piiz+kz*piix)
                  repress(6) = repress(6)
     &                 + factor*0.5d0*(ky*piiz+kz*piiy)
               enddo
            enddo
         endif
      enddo

      return
      end
