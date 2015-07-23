#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2000-2014 Marcus G. Martin                       *
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
      function twh_vangle(lfailure,lovrlap,ibox,imolty,iunit,ibend,iflag
     &     ,angcos,distone,disttwo,xcord,ycord,zcord)
c     ******************************************************************
c     * computes the bond bending energy given an angle type and the   *
c     * angle.                                                         *
c     * ibox: the box the molecule is in (used for mimage)             *
c     * imolty: the molecule type                                      *
c     * iunit: the unit of the first atom in the bend                  *
c     * ibend: the index of the bend                                   *
c     * iflag: 0 if we should use the angcos,distone,disttwo info      *
c     *        1 if we should use the xcord,ycord,zcord info           *
c     *       -1 if we want to return the angle, instead of the        *
c     *             energy, from the xcord,ycord,zcord info            *
c     * xcord: xcoordinates of the three atoms                         *
c     * ycord: ycoordinates of the three atoms                         *
c     * zcord: zcoordinates of the three atoms                         *
c     * angcos: the cosine of the bond bending angle                   *
c     * distone: the distance between atoms 1-2                        *
c     * disttwo: the distance between atoms 2-3                        *
c     * lovrlap: returned as true if there is a hard overlap in the    *
c     *          nonbonded part of certain bending potentials          *
c     *          Note: we are a bit tricky as lovrlap is not intialized*
c     *          in this routine, should be intialized by the caller   *
c     *                                                                *
c     * rewritten from old sumup subroutines 01-21-2000 by M.G. Martin *
c     * last modified 07-21-2014 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_ARCCOS
#define FUNCTION_CMP_GT
#define FUNCTION_CMP_LT
#define FUNCTION_DOTPRODUCT
#define FUNCTION_GET_ANGLESTYLE
#define FUNCTION_GET_BENCOEFF
#define FUNCTION_GET_BENPREFACT
#define FUNCTION_GET_IBTBEN1
#define FUNCTION_GET_IBTBEN2
#define FUNCTION_GET_IJBEN3
#define FUNCTION_GET_ITBEN
#define FUNCTION_GET_LBONANG
#define FUNCTION_GET_VIBCOEFF
#define FUNCTION_LAWOFCOSINE
#define FUNCTION_VANGLE
#define FUNCTION_VTWOBODY
#include "functions.h"

#define DEBUG_VANGLE 0
c     --- variables passed to and from the function
      logical lfailure,lovrlap
      integer ibox,imolty,iunit,ibend,iflag
      double precision distone,disttwo,angcos
      double precision xcord,ycord,zcord
      dimension xcord(3),ycord(3),zcord(3)
c     --- local variables
c     --- logical scalars
      logical ldononbond
c     --- integer scalars
      integer style,iutype,jutype,type,junit,iorder,bone,btwo,iallow
c     --- double precision scalars
      double precision dist,distsq,qqa,qqb
      double precision angle,differ,value,prefac,diffsquare,diffcube
      double precision xone,yone,zone,xtwo,ytwo,ztwo,dv,qqfact
#if DEBUG_VANGLE
      double precision radtodeg
#endif
c     --- double precision parameters
      double precision tolerance
      parameter (tolerance = 1d-5)
c     --- double precision arrays
      double precision vf
      dimension vf(MAX_FOREIGN_LAMBDA)

      twh_vangle = 0.0d0
      ldononbond = .false.
      if ( abs(iflag) .eq. 1 ) then
c        --- compute the bond lengths and angles from the coordinates
         call twh_mimage(ibox,xcord(1),ycord(1),zcord(1)
     &        ,xcord(2),ycord(2),zcord(2),xone,yone,zone )
         call twh_mimage(ibox,xcord(3),ycord(3),zcord(3)
     &        ,xcord(2),ycord(2),zcord(2),xtwo,ytwo,ztwo )
         angcos = twh_dotproduct(xone,yone,zone,xtwo,ytwo,ztwo,distone
     &        ,disttwo)
         if ( iflag .eq. -1 ) then
            twh_vangle = twh_arccos(angcos)
            return
         endif
      endif

#if DEBUG_VANGLE
      write(6,*) 'imolty,iunit,ibend:',imolty,iunit,ibend
      angle = twh_arccos(angcos)
      call twh_constant_radtodeg(GLB_GET,radtodeg)
      write(6,*) 'angle:',angle*radtodeg
#endif

c     --- determine type, style and order of the angle
      type = twh_get_itben(imolty,iunit,ibend)
      call twh_orderben(GLB_GET,imolty,iunit,ibend,iorder)

      style = twh_get_anglestyle(type)

      if ( style .eq. ANG_STYLE_RIGID ) then
c        --- rigid angle, make sure it is within tolerance, else overlap
         angle = twh_arccos(angcos)
         value = (angle - twh_get_bencoeff(type,0))**2
         if (twh_cmp_gt(value, twh_get_bencoeff(type,1)) ) then 
c           --- out of tolerance, hard overlap
            lovrlap = .true.
#if DEBUG_VANGLE
            write(6,*) 'vangle angle:',angle*radtodeg
            call twh_constant_radtodeg(GLB_GET,radtodeg)
            write(6,*) 'coeff angle:',twh_get_bencoeff(type,0)*radtodeg
#endif
         else
c           --- in tolerance, no energy
            twh_vangle = 0.0d0
         endif
      elseif ( style .eq. ANG_STYLE_HARMONIC ) then
c        --- standard harmonic bond bending potential
c        --- energy = k* (angle - equil)**2
         angle = twh_arccos(angcos)
         twh_vangle = twh_get_bencoeff(type,1)
     &        *(angle-twh_get_bencoeff(type,0))**2
      elseif ( style .eq. ANG_STYLE_DREIDING ) then
c        --- Dreiding K[ 1 + Cos(angle)]
         twh_vangle = twh_get_bencoeff(type,1)*(1.0d0 + angcos)
      elseif ( style .eq. ANG_STYLE_HARM_COS ) then
c        --- harmonic cosine potential
c        --- energy = k* (cos(angle) - cos(equil))**2
         twh_vangle = twh_get_bencoeff(type,1)
     &        *(angcos - twh_get_bencoeff(type,-1))**2
      elseif ( style .eq. ANG_STYLE_COMPASS_AUTO ) then
c        --- Compass quartic angle with autodetection
         angle = twh_arccos(angcos)
         differ = angle - twh_get_bencoeff(type,0)
         value = differ*differ
         twh_vangle = twh_get_bencoeff(type,1)*(value)
     &        + twh_get_bencoeff(type,2)*(value*differ)
     &        + twh_get_bencoeff(type,3)*(value*value)
         call twh_ntype(GLB_GET,imolty,iunit,iutype)
         if ( twh_get_lbonang(type) ) then
c           --- add in both of the bond-angle cross terms
            bone = twh_get_ibtben1(imolty,iunit,ibend)
            btwo = twh_get_ibtben2(imolty,iunit,ibend)
            if ( iorder .eq. ORD_FORWARD ) then
c              --- force constants match up with bond order
               prefac = twh_get_bencoeff(type,4)*(distone 
     &              - twh_get_vibcoeff(bone,0))
     &              + twh_get_bencoeff(type,5)*(disttwo 
     &              - twh_get_vibcoeff(btwo,0))
            elseif ( iorder .eq. ORD_REVERSE ) then
c              --- force constants match with inverted bond order
               prefac = twh_get_bencoeff(type,5)*(distone 
     &              - twh_get_vibcoeff(bone,0))
     &              + twh_get_bencoeff(type,4)*(disttwo 
     &              - twh_get_vibcoeff(btwo,0))
            else
               write(6,*) 'iutype',iutype,'ang type',type
               write(6,*) 'iorder',iorder
               write(6,*) 'VANGLE: bad type'
               lfailure = .true.
               return
            endif
            twh_vangle = twh_vangle + prefac*differ
         endif
      elseif ( style .eq. ANG_STYLE_CHARMM_UB ) then
c         --- Charmm harmonic with urey-bradley
         angle = twh_arccos(angcos)
         twh_vangle = twh_get_bencoeff(type,1)
     &        *((angle-twh_get_bencoeff(type,0))**2)

c        --- compute the distance between the 1-3 atoms
         dist = twh_lawofcosine(distone,disttwo,angcos)
         twh_vangle = twh_vangle 
     &        + twh_get_bencoeff(type,3)
     &        *((dist-twh_get_bencoeff(type,2))**2)
      elseif ( style .eq. ANG_STYLE_NONBOND ) then
c        --- nonbonded interaction between the 1-3 atoms
         ldononbond = .true.
      elseif ( style .eq. ANG_STYLE_NB_HARM ) then
c        --- harmonic plus nonbond and coulombic between the 1-3 atoms
         angle = twh_arccos(angcos)
         twh_vangle = twh_get_bencoeff(type,1)
     &        *(angle-twh_get_bencoeff(type,0))**2
         ldononbond = .true.
      elseif ( style .eq. ANG_STYLE_COMPASS_EX ) then
c        --- Compass quartic angle with explicit ordering of terms
         angle = twh_arccos(angcos)
         differ = angle - twh_get_bencoeff(type,0)
         value = differ*differ
         twh_vangle = twh_get_bencoeff(type,1)*(value)
     &        + twh_get_bencoeff(type,2)*(value*differ)
     &        + twh_get_bencoeff(type,3)*(value*value)
         call twh_ntype(GLB_GET,imolty,iunit,iutype)
         if ( twh_get_lbonang(type) ) then
c           --- add in both of the bond angle cross terms
            if ( iorder .eq. ORD_FORWARD ) then
c              --- force constants match up with bond order
               prefac = twh_get_bencoeff(type,4)
     &              *(distone - twh_get_bencoeff(type,5))
     &              + twh_get_bencoeff(type,6)
     &              *(disttwo - twh_get_bencoeff(type,7))
            elseif ( iorder .eq. ORD_REVERSE ) then
c              --- force constants match with inverted bond order
               prefac = twh_get_bencoeff(type,6)
     &              *(distone - twh_get_bencoeff(type,7))
     &              + twh_get_bencoeff(type,4)
     &              *(disttwo - twh_get_bencoeff(type,5))
            else
               write(6,*) 'VANGLE: bad iorder value ',iorder
               lfailure = .true.
               return
            endif
            twh_vangle = twh_vangle + prefac*differ
         endif
      elseif ( style .eq. ANG_STYLE_FOURIER ) then
c        --- Fourier expansion with constant minus single term
         angle = twh_arccos(angcos)
         twh_vangle = twh_get_benprefact(imolty,iunit,ibend)
     &        *(1.0d0 + twh_get_bencoeff(type,2)
     &        *dcos( twh_get_bencoeff(type,1)*angle))
      elseif ( style .eq. ANG_STYLE_TTFOURIER ) then
c        --- 3 term Fourier expansion
         angle = twh_arccos(angcos)
         twh_vangle = twh_get_benprefact(imolty,iunit,ibend)
     &        *(twh_get_bencoeff(type,1)
     &        + twh_get_bencoeff(type,2)*angcos
     &        + twh_get_bencoeff(type,3)*dcos(2.0d0*angle) )
      elseif ( style .eq. ANG_STYLE_NONE ) then
c        --- no energy for this term
         twh_vangle = 0.0d0
      elseif ( style .eq. ANG_STYLE_SIXTH_AUTO ) then
c        --- MM2 sixth power angle with autodetection
         angle = twh_arccos(angcos)
         differ = angle - twh_get_bencoeff(type,0)
         value = differ*differ
         twh_vangle = twh_get_bencoeff(type,1)*value*
     &        (1.0d0 + twh_get_bencoeff(type,2)*value*value)
         if ( twh_get_lbonang(type) ) then
c           --- add in both of the bond-angle cross terms
c           --- check based upon only the first letter for this term
            call twh_ntype(GLB_GET,imolty,iunit,iutype)
            if ( iorder .eq. ORD_FORWARD ) then
c              --- force constants match up with bond order
               bone = twh_get_ibtben1(imolty,iunit,ibend)
               btwo = twh_get_ibtben2(imolty,iunit,ibend)
            elseif ( iorder .eq. ORD_REVERSE ) then
c              --- force constants match with inverted bond order
               bone = twh_get_ibtben2(imolty,iunit,ibend)
               btwo = twh_get_ibtben1(imolty,iunit,ibend)
            else
               write(6,*) 'iutype',iutype,'ang type',type
               write(6,*) 'orderben:',iorder
               write(6,*) 'VANGLE: bad type'
               lfailure = .true.
               return
            endif
            twh_vangle = twh_vangle + twh_get_bencoeff(type,4)*differ*(
     &           ( distone - twh_get_vibcoeff(bone,0) ) +
     &           ( disttwo - twh_get_vibcoeff(btwo,0) ) )
         endif
      elseif ( style .eq. ANG_STYLE_SQUARE_WELL ) then
c        --- infinite square well bending angle
c        --- compute the distance between the 1-3 atoms
         dist = twh_lawofcosine(distone,disttwo,angcos)
         if ( twh_cmp_lt(dist, twh_get_bencoeff(type,1)) .or.
     &        twh_cmp_gt(dist, twh_get_bencoeff(type,2)) ) then
            lovrlap = .true.
         else
            twh_vangle = twh_get_bencoeff(type,3)
         endif
      elseif ( style .eq. ANG_STYLE_MULTI_RIGID ) then
c        --- rigid angle with two allowed angles
c        --- make sure it is within tolerance of one of the allowed
c        --- angles, else overlap
         lovrlap = .true.
         angle = twh_arccos(angcos)
         do iallow = 0,1
            value = (angle - twh_get_bencoeff(type,iallow))**2
            if ( .not. twh_cmp_gt(value, twh_get_bencoeff(type,2)) 
     &           ) then 
c              --- in tolerance, no energy
               lovrlap = .false.
               twh_vangle = 0.0d0
            endif
         enddo
      elseif ( style .eq. ANG_STYLE_MMFF ) then
c        --- MMFF cubic plus a bond-angle cross term
         angle = twh_arccos(angcos)
         differ = angle - twh_get_bencoeff(type,0)
         diffsquare = differ*differ
         diffcube = differ*diffsquare
         twh_vangle = twh_get_bencoeff(type,1)*diffsquare
     &        + twh_get_bencoeff(type,2)*diffcube
c        --- cross terms
         call twh_ntype(GLB_GET,imolty,iunit,iutype)
         if ( twh_get_lbonang(type) ) then
c           --- add in both of the bond-angle cross terms
            bone = twh_get_ibtben1(imolty,iunit,ibend)
            btwo = twh_get_ibtben2(imolty,iunit,ibend)
            if ( iorder .eq. ORD_FORWARD ) then
c              --- force constants match up with bond order
               prefac = twh_get_bencoeff(type,3)
     &              *(distone - twh_get_vibcoeff(bone,0))
     &              + twh_get_bencoeff(type,4)
     &              *(disttwo - twh_get_vibcoeff(btwo,0))
            elseif ( iorder .eq. ORD_REVERSE ) then
c              --- force constants match with inverted bond order
               prefac = twh_get_bencoeff(type,4)
     &              *(distone - twh_get_vibcoeff(bone,0))
     &              + twh_get_bencoeff(type,3)
     &              *(disttwo - twh_get_vibcoeff(btwo,0))
            else
               write(6,*) 'iutype',iutype,'ang type',type
               write(6,*) 'benorder:',iorder
               write(6,*) 'VANGLE: bad type'
               lfailure = .true.
               return
            endif
            twh_vangle = twh_vangle + prefac*differ
         endif
      elseif ( style .eq. ANG_STYLE_HARM_COS_OT ) then
c        --- harmonic cosine plus nonbond between the 1-3 atoms
         twh_vangle = twh_get_bencoeff(type,1)
     &        *(angcos - twh_get_bencoeff(type,-1))**2
         ldononbond = .true.
      else
         write(6,*) 'VANGLE: non defined bond angle type'
         lfailure = .true.
         return
      endif

      if ( ldononbond ) then
c        --- compute the distance between the 1-3 atoms
         dist = twh_lawofcosine(distone,disttwo,angcos)
         distsq = dist*dist
c        --- get the nonbonded vdw energy 
         junit = twh_get_ijben3(imolty,iunit,ibend)
         call twh_ntype(GLB_GET,imolty,iunit,iutype)
         call twh_ntype(GLB_GET,imolty,junit,jutype)
c        --- assign the charges
         call twh_qqatom(GLB_GET,imolty,iunit,qqa)
         call twh_qqatom(GLB_GET,imolty,junit,qqb)
         twh_vangle = twh_vangle + twh_vtwobody(lfailure,distsq
     &        ,iutype,jutype,0,lovrlap,vf,dv)
         if ( lfailure ) return
c        --- get the nonbonded coulombic energy (only the real space
c        --- part)
         call twh_constant_qqfact(GLB_GET,qqfact)
         twh_vangle = twh_vangle + qqfact*qqa*qqb/dist
      endif

      return
      end
