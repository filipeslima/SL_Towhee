#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2010 Marcus G. Martin                       *
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
      function twh_vtorsion(lfailure,ibox,cordstyle,imolty,unit_1
     &,torind,xcord,ycord,zcord,phic,loverlap)
c     ******************************************************************
c     * Computes the torsional energy plus cross terms and 1-4 nonbond *
c     * for a particular torsion type.                                 *
c     * ibox: the box the molecule is in - needed for minimum image    *
c     * cordstyle: tells the function what kind of angle we want to    *
c     *            work with and whether it is passed to the function  *
c     *            directly, or should be computed from the coordinates*
c     *            -1: normal dihedral passed as an angle              *
c     *             0: special case, compute the normal dihedral angle *
c     *                from the coordinates and then just return that  *
c     *                angle instead of returning the torsional energy *
c     *                This is used for some of the non-uniform        *
c     *                sampling.                                       *
c     *             1: normal dihedral that must be computed from the  *
c     *                coordinates.  (1-2 x 2-3 * 2-3 x 3-4)           *
c     * imolty: molecule type                                          *
c     * unit_1: unit numbers of the starting atom of the torsion       *
c     * torind: the index of the torsion                               *
c     * x,y,z cord : coordinates of the four atoms in the torsion      *
c     * phic: cosine of the torsional angle                            *
c     * loverlap: true for hard overlap, note it is initialized only in*
c     *  the calling routine.  Do not initialize in any of the vxx     *
c     *  functions                                                     *
c     *                                                                *
c     * originally written prior to 1995                               *
c     * last modified 01-29-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_ARCCOS
#define FUNCTION_CMP_GT
#define FUNCTION_CMP_LT
#define FUNCTION_DISTANCE
#define FUNCTION_EXPON
#define FUNCTION_DOTPRODUCT
#define FUNCTION_GET_BENCOEFF
#define FUNCTION_GET_CTNAME
#define FUNCTION_GET_IJTOR2
#define FUNCTION_GET_IJTOR3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_INVIB
#define FUNCTION_GET_ITTOR
#define FUNCTION_GET_NTORLOOP
#define FUNCTION_GET_NTYPE
#define FUNCTION_GET_TORCOEFF
#define FUNCTION_GET_TORSTRING
#define FUNCTION_GET_TORSTYLE
#define FUNCTION_GET_VIBCOEFF
#define FUNCTION_VCOULOMB
#define FUNCTION_VTORSION
#define FUNCTION_VTWOBODY
#include "functions.h"
#define DEBUG_VTORSION 0
c     --- variables passed to and from function
      logical lfailure,loverlap
      integer imolty,torind,unit_1,cordstyle,ibox
      double precision phic
      double precision xcord,ycord,zcord
      dimension xcord(4),ycord(4),zcord(4)

c     --- local integer scalars
      integer itor,ntii,ntjj,bij,bjk,bkl,aijk,ajkl,i,style
      integer torffnum,itmp
      integer abcord,ivec,unit_2,unit_3,unit_4,index,torofcode
c     --- local integer arrays
      integer atype,btype
      dimension btype(3)
      dimension atype(2)
c     --- local double precision scalars
      double precision phi,rijsq,cosphi,cos2phi,cos3phi,rij
      double precision xaa1,yaa1,zaa1,xaa2,yaa2,zaa2,daa1,daa2
      double precision xcc,ycc,zcc,tcc,qqfact
      double precision torscale,rxij,ryij,rzij
      double precision xone,xtwo,yone,ytwo,zone,ztwo,cosshiftphi
      double precision totbond,dv,onepi,vonefour
c     --- local double precision arrays
      double precision xvec,yvec,zvec
      dimension xvec(4),yvec(4),zvec(4)
      double precision length
      dimension length(4)
      double precision angle
      dimension angle(2)
      double precision vf
      dimension vf(MAX_FOREIGN_LAMBDA)

#if DEBUG_VTORSION
      write(6,*) 'Start VTORSION'
#endif
c     --- initialize vtorsion
      twh_vtorsion = 0.0d0

c     --- determine the torsion forfield number and style
      abcord = iabs(cordstyle)
      if ( abcord .eq. 0 ) then
c        --- just computing an angle and returning, set a dummy ff
         torffnum = 1
      elseif ( abcord .eq. 1 ) then
c        --- regular torsion - use ittor
         torffnum = abs(twh_get_ittor(imolty,unit_1,torind))
      else
         write(6,*) 'invalid cordstyle in VTORSION',cordstyle
         lfailure = .true.
         return
      endif
      
      style = twh_get_torstyle(torffnum)
#if DEBUG_VTORSION
      write(6,*) 'torsion style:',style
      write(6,*) 'cord style:',cordstyle
#endif

      if ( cordstyle .lt. 0 ) then
c        --- use the angle that was passed to the fucntion
         if ( style .eq. 5 ) then
            write(6,*) 'cannot pass angles to VTORSION for Class 2'
            lfailure = .true.
            return
         endif
      elseif ( cordstyle .eq. 0 ) then
c        --- want the standard angle, but need it over the full range
c        --- of -Pi to Pi so we can use this as an adjustment in the
c        --- non-uniform selection of torsional angles.
         do ivec = 1,3
            call twh_mimage(ibox,xcord(ivec+1),ycord(ivec+1)
     &           ,zcord(ivec+1),xcord(ivec),ycord(ivec),zcord(ivec)
     &           ,xvec(ivec),yvec(ivec),zvec(ivec))
         enddo
c        --- we also need a fourth vector for the triple product
         call twh_mimage(ibox,xcord(2),ycord(2),zcord(2)
     &        ,xcord(3),ycord(3),zcord(3),xvec(4),yvec(4),zvec(4))
c        --- calculate cross products 1-2 x 2-3 
         call twh_crossproduct(xvec(1),yvec(1),zvec(1),xvec(2),yvec(2)
     &        ,zvec(2),xaa1,yaa1,zaa1)
c        --- calculate cross products 3-4 x 2-3
         call twh_crossproduct(xvec(3),yvec(3),zvec(3),xvec(2),yvec(2)
     &        ,zvec(2),xaa2,yaa2,zaa2)
c        --- calculate dot product of cross products
         phic = twh_dotproduct(xaa1,yaa1,zaa1,xaa2,yaa2,zaa2,daa1,daa2) 
         phic = twh_arccos(-phic)
c        --- calculate cross product of cross products
         call twh_crossproduct(xaa1,yaa1,zaa1,xaa2,yaa2,zaa2
     &        ,xcc,ycc,zcc)
c        --- calculate scalar triple product
         tcc = xcc*(xvec(4)) + ycc*(yvec(4)) + zcc*(zvec(4))
c        --- compute the angle over a range -pi to +pi
         if ( tcc .gt. 0.0d0 ) phic = -phic
         twh_vtorsion = phic
         return
      elseif ( cordstyle .eq. 1 ) then
c        --- compute the dihedral angle
c        --- defined so that a cis bond is at 0.0 phi
         do ivec = 1,3
            call twh_mimage(ibox,xcord(ivec+1),ycord(ivec+1)
     &           ,zcord(ivec+1),xcord(ivec),ycord(ivec),zcord(ivec)
     &           ,xvec(ivec),yvec(ivec),zvec(ivec))
         enddo
c        --- calculate cross products 1-2 x 2-3 
         call twh_crossproduct(xvec(1),yvec(1),zvec(1),xvec(2),yvec(2)
     &        ,zvec(2),xaa1,yaa1,zaa1)
c        --- calculate cross products 3-4 x 2-3
         call twh_crossproduct(xvec(3),yvec(3),zvec(3),xvec(2),yvec(2)
     &        ,zvec(2),xaa2,yaa2,zaa2)
c        --- calculate dot product of cross products
         phic = -twh_dotproduct(xaa1,yaa1,zaa1,xaa2,yaa2,zaa2,daa1,daa2)
      endif

      if ( style .eq. 1 ) then
c        --- simple harmonic potential 
c        --- used for TraPPE double bonds (cis vs trans) 
         phi = twh_arccos(phic)
         twh_vtorsion = twh_get_torcoeff(torffnum,0)
     &        *(phi-twh_get_torcoeff(torffnum,1))**2.0d0
      elseif ( style .eq. 2  ) then
c        --- OPLS style cosine series 
c        --- used in the OPLS-AA force field
         phi = twh_arccos(phic)
         twh_vtorsion = twh_get_torcoeff(torffnum,1)*(1.d0 + phic) +
     &        twh_get_torcoeff(torffnum,2)*(1.d0 - dcos(2.0d0*phi)) +
     &        twh_get_torcoeff(torffnum,3)*(1.d0 + dcos(3.0d0*phi))
      elseif ( style .eq. 3 ) then
c        --- Gromos/Charmm/Amber cosine potential
         phi = twh_arccos(phic)
         do itor = 1,twh_get_ntorloop(torffnum)
            index = 3*(itor-1)
            twh_vtorsion = twh_vtorsion 
     &           + twh_get_torcoeff(torffnum,index+1)* 
     &           ( 1.0d0 + dcos(twh_get_torcoeff(torffnum,index+2)*phi
     &           - twh_get_torcoeff(torffnum,index+3))
     &           )
         enddo
      elseif ( style .eq. 4 ) then
c        --- Gromos/Charmm/Amber cosine potential
c        --- plus a harmonic term that Charmm would traditionally call
c        --- an improper torsion despite the fact that the bonding 
c        --- pattern is actually that of a regular torsion
         phi = twh_arccos(phic)
         do itor = 1,twh_get_ntorloop(torffnum)
            index = 3*(itor-1)
            twh_vtorsion = twh_vtorsion 
     &           + twh_get_torcoeff(torffnum,index+1)* 
     &           ( 1.0d0 + dcos(twh_get_torcoeff(torffnum,index+2)*phi
     &           - twh_get_torcoeff(torffnum,index+3))
     &           )
         enddo
c        --- reset index
         index = 3*twh_get_ntorloop(torffnum)
c        --- add in the harmonic term
         twh_vtorsion = twh_vtorsion 
     &        + twh_get_torcoeff(torffnum,index+1)
     &        *(phi-twh_get_torcoeff(torffnum,index+2))**2.0d0

      elseif ( style .eq. 5 .or. style .eq. 6) then
c        --- Compass cosine series using the clever atom-based terms 
c        --- to sort out all of the cross terms
c        --- style 5 includes cross terms
c        --- style 6 does not include cross terms
         phi = twh_arccos(phic)
         cosphi = phic
         cos2phi = dcos(2.0d0*phi)
         cos3phi = dcos(3.0d0*phi)

c        --- simple torsion part
         twh_vtorsion = twh_get_torcoeff(torffnum,0)*(1.0d0 - cosphi)
     &        + twh_get_torcoeff(torffnum,1)*(1.0d0 - cos2phi)
     &        + twh_get_torcoeff(torffnum,2)*(1.0d0 - cos3phi)

         if ( style .eq. 5 ) then
c           --- actually need all of the unit numbers to determine the
c           --- order of the torsion
            unit_2 = twh_get_ijtor2(imolty,unit_1,torind)
            unit_3 = twh_get_ijtor3(imolty,unit_1,torind)
            unit_4 = twh_get_ijtor4(imolty,unit_1,torind)

c           --- Compute the cross terms
c           --- Determine whether we are in normal or inverse order
            if ( twh_get_ctname(twh_get_ntype(imolty,unit_1)) 
     &           .eq. twh_get_torstring(torffnum,1) .and.
     &           twh_get_ctname(twh_get_ntype(imolty,unit_2))
     &           .eq. twh_get_torstring(torffnum,2) .and.
     &           twh_get_ctname(twh_get_ntype(imolty,unit_3))
     &           .eq. twh_get_torstring(torffnum,3) .and.
     &           twh_get_ctname(twh_get_ntype(imolty,unit_4))
     &           .eq. twh_get_torstring(torffnum,4) ) then
c              --- normal order
               bij = 1
               bjk = 2
               bkl = 3
               aijk = 1
               ajkl = 2
            elseif ( twh_get_ctname(twh_get_ntype(imolty,unit_4)) 
     &              .eq. twh_get_torstring(torffnum,1) .and.
     &              twh_get_ctname(twh_get_ntype(imolty,unit_3)) 
     &              .eq. twh_get_torstring(torffnum,2) .and.
     &              twh_get_ctname(twh_get_ntype(imolty,unit_2))
     &              .eq. twh_get_torstring(torffnum,3) .and.
     &              twh_get_ctname(twh_get_ntype(imolty,unit_1))
     &              .eq. twh_get_torstring(torffnum,4) ) then
c              --- inverse order
               bij = 3
               bjk = 2
               bkl = 1
               aijk = 2
               ajkl = 1
            else
               write(6,*) 'VTORSION: torffnum',torffnum
               write(6,*) 'ctname'
     &              ,twh_get_ctname(twh_get_ntype(imolty,unit_1))
     &              ,twh_get_ctname(twh_get_ntype(imolty,unit_2))
     &              ,twh_get_ctname(twh_get_ntype(imolty,unit_3))
     &              ,twh_get_ctname(twh_get_ntype(imolty,unit_4))
               write(6,*) 'torstring'
     &              ,(twh_get_torstring(torffnum,i),i=1,4)
               write(6,*) 'torsion does not match up in VTORSION'
               lfailure = .true.
               return
            endif

c           --- compute bond lengths using *vec (already mimaged)
            do itmp = 1,3
               rxij = xvec(itmp)
               ryij = yvec(itmp)
               rzij = zvec(itmp)
               length(itmp) = twh_distance(rxij,ryij,rzij)
            enddo
c           --- compute angles using *vec (already mimaged)
            do itmp = 1,2
               xone = -xvec(itmp)
               xtwo = xvec(itmp+1)
               yone = -yvec(itmp)
               ytwo = yvec(itmp+1)
               zone = -zvec(itmp)
               ztwo = zvec(itmp+1)
               angle(itmp) = twh_dotproduct(xone,yone,zone,xtwo,ytwo
     &              ,ztwo,daa1,daa2)
            enddo
c           --- store the bond and angle types
            call twh_ibttor1(GLB_GET,imolty,unit_1,torind,btype(1))
            call twh_ibttor2(GLB_GET,imolty,unit_1,torind,btype(2))
            call twh_ibttor3(GLB_GET,imolty,unit_1,torind,btype(3))
            call twh_iattor1(GLB_GET,imolty,unit_1,torind,atype(1))
            call twh_iattor2(GLB_GET,imolty,unit_1,torind,atype(2))

c           --- convert angles from cosine to normal
            angle(1) = twh_arccos(angle(1))
            angle(2) = twh_arccos(angle(2))
c           --- bond(ij) - torsion cross
            twh_vtorsion = twh_vtorsion 
     &           + (length(bij) -twh_get_vibcoeff(btype(bij),0))*
     &           (twh_get_torcoeff(torffnum,3)*cosphi 
     &           + twh_get_torcoeff(torffnum,4)*cos2phi
     &           + twh_get_torcoeff(torffnum,5)*cos3phi)
c           --- bond(jk) - torsion cross
     &           + (length(bjk) - twh_get_vibcoeff(btype(bjk),0))*
     &           (twh_get_torcoeff(torffnum,6)*cosphi 
     &           + twh_get_torcoeff(torffnum,7)*cos2phi
     &           + twh_get_torcoeff(torffnum,8)*cos3phi)
c           --- bond(kl) - torsion cross
     &           + (length(bkl) - twh_get_vibcoeff(btype(bkl),0))*
     &           (twh_get_torcoeff(torffnum,9)*cosphi 
     &           + twh_get_torcoeff(torffnum,10)*cos2phi
     &           + twh_get_torcoeff(torffnum,11)*cos3phi)
c           --- angle(ijk) - torsion cross
     &           + (angle(aijk) - twh_get_bencoeff(atype(aijk),0))*
     &           (twh_get_torcoeff(torffnum,12)*cosphi 
     &           + twh_get_torcoeff(torffnum,13)*cos2phi
     &           + twh_get_torcoeff(torffnum,14)*cos3phi)
c           --- angle(jkl) - torsion cross
     &           + (angle(ajkl) - twh_get_bencoeff(atype(ajkl),0))*
     &           (twh_get_torcoeff(torffnum,15)*cosphi 
     &           + twh_get_torcoeff(torffnum,16)*cos2phi
     &           + twh_get_torcoeff(torffnum,17)*cos3phi)
c           --- angle(ijk) - angle(jkl) - torsion cross
     &           + twh_get_torcoeff(torffnum,18)*
     &           (angle(aijk) - twh_get_bencoeff(atype(aijk),0))*
     &           (angle(ajkl) - twh_get_bencoeff(atype(ajkl),0))*cosphi
c           --- bond(ij) - bond(kl) cross
     &           + twh_get_torcoeff(torffnum,19)*
     &           (length(bij) - twh_get_vibcoeff(btype(bij),0))*
     &           (length(bkl) - twh_get_vibcoeff(btype(bkl),0))
         endif

      elseif ( style .eq. 7 ) then
c        --- TraPPE simple cosine function
c        --- trappe is defined so that trans is at 0.0
         call twh_constant_pi(GLB_GET,onepi)
         phi = twh_arccos(phic)-onepi
         twh_vtorsion = twh_get_torcoeff(torffnum,0)
     &        *(1.0d0 - dcos(2.0d0*phi + twh_get_torcoeff(torffnum,1)))

      elseif ( style .eq. 8 ) then
c        --- nonbonded potential only -or absolutely no torsional energy
c        --- any scale factor is accounted for in ffnonbond
c        --- nothing computed here, all down below

      elseif ( style .eq. 9) then
c        --- Compass cosine series plus cross terms using an explicit
c        --- parameter declaration and negative typing to determine 
c        --- the order
         phi = twh_arccos(phic)
         cosshiftphi = dcos(phi - twh_get_torcoeff(torffnum,3))
c        --- simple torsion part
         twh_vtorsion = 
     &        twh_get_torcoeff(torffnum,0)
     &        *(1.0d0 - cosshiftphi )
     &        + twh_get_torcoeff(torffnum,1)
     &        *(1.0d0 - dcos(2.0d0*phi - twh_get_torcoeff(torffnum,4)) )
     &        + twh_get_torcoeff(torffnum,2)
     &        *(1.0d0 - dcos(3.0d0*phi - twh_get_torcoeff(torffnum,5)) )

c        --- the unshifted values of the cosines are used in the cross
c        --- terms
         cosphi = dcos(phi)
         cos2phi = dcos(2.0d0*phi)
         cos3phi = dcos(3.0d0*phi)
c        --- compute bond lengths using *vec (already mimaged)
         do itmp = 1,3
            rxij = xvec(itmp)
            ryij = yvec(itmp)
            rzij = zvec(itmp)
            length(itmp) = twh_distance(rxij,ryij,rzij)
         enddo
c        --- compute angles using *vec (already mimaged)
         do itmp = 1,2
            xone = -xvec(itmp)
            xtwo = xvec(itmp+1)
            yone = -yvec(itmp)
            ytwo = yvec(itmp+1)
            zone = -zvec(itmp)
            ztwo = zvec(itmp+1)
            angle(itmp) = twh_arccos(twh_dotproduct(xone,yone,zone
     &           ,xtwo,ytwo,ztwo,daa1,daa2))
         enddo
c        --- Compute the cross terms
c        --- Determine whether we are in normal or inverse order
         if ( twh_get_ittor(imolty,unit_1,torind) .gt. 0 ) then
c           --- normal order, don't change a thing
            bij = 1
            bjk = 2
            bkl = 3
            aijk = 1
            ajkl = 2
         else
c           --- inverse order, switch all of the bond lengths and the
c           --- angles
            bij = 3
            bjk = 2
            bkl = 1
            aijk = 2
            ajkl = 1
         endif
         twh_vtorsion = twh_vtorsion 
c        --- bond(jk) - torsion cross
     &        + (length(bjk) - twh_get_torcoeff(torffnum,9))*
     &        (twh_get_torcoeff(torffnum,6)*cosphi 
     &        + twh_get_torcoeff(torffnum,7)*cos2phi
     &        + twh_get_torcoeff(torffnum,8)*cos3phi)
c        --- bond(ij) - torsion cross
     &        + (length(bij) - twh_get_torcoeff(torffnum,13))*
     &        (twh_get_torcoeff(torffnum,10)*cosphi 
     &        + twh_get_torcoeff(torffnum,11)*cos2phi
     &        + twh_get_torcoeff(torffnum,12)*cos3phi)
c        --- bond(kl) - torsion cross
     &        + (length(bkl) - twh_get_torcoeff(torffnum,17))*
     &        (twh_get_torcoeff(torffnum,14)*cosphi 
     &        + twh_get_torcoeff(torffnum,15)*cos2phi
     &        + twh_get_torcoeff(torffnum,16)*cos3phi)
c        --- angle(ijk) - torsion cross
     &        + (angle(aijk) - twh_get_torcoeff(torffnum,24))*
     &        (twh_get_torcoeff(torffnum,18)*cosphi 
     &        + twh_get_torcoeff(torffnum,19)*cos2phi
     &        + twh_get_torcoeff(torffnum,20)*cos3phi)
c        --- angle(jkl) - torsion cross
     &        + (angle(ajkl) - twh_get_torcoeff(torffnum,25))*
     &        (twh_get_torcoeff(torffnum,21)*cosphi 
     &        + twh_get_torcoeff(torffnum,22)*cos2phi
     &        + twh_get_torcoeff(torffnum,23)*cos3phi)
c        --- angle(ijk) - angle(jkl) - torsion cross
     &        + twh_get_torcoeff(torffnum,26)*
     &        (angle(aijk) - twh_get_torcoeff(torffnum,27))*
     &        (angle(ajkl) - twh_get_torcoeff(torffnum,28))*cosshiftphi
c        --- bond(ij) - bond(kl) cross
     &        + twh_get_torcoeff(torffnum,29)*
     &        (length(bij) - twh_get_torcoeff(torffnum,30))*
     &        (length(bkl) - twh_get_torcoeff(torffnum,31))

      elseif ( style .eq. 10 ) then
c        --- cosine power series
         do itor = 0,twh_get_ntorloop(torffnum)
            twh_vtorsion = twh_vtorsion 
     &           + twh_get_torcoeff(torffnum,itor)*(phic**itor)
         enddo
      elseif ( style .eq. 11 ) then
c        --- old OPLS style cosine series 
c        --- used in the OPLS-UA force field
         phi = twh_arccos(phic)
         twh_vtorsion = twh_get_torcoeff(torffnum,0) + 
     &        twh_get_torcoeff(torffnum,1)*(1.d0 + phic) +
     &        twh_get_torcoeff(torffnum,2)*(1.d0 - dcos(2.0d0*phi)) +
     &        twh_get_torcoeff(torffnum,3)*(1.d0 + dcos(3.0d0*phi))
      elseif ( style .eq. 12 ) then
c        --- Sum2003 style torsion potential
         phi = twh_arccos(phic)
         do itor = 0,twh_get_ntorloop(torffnum)
            index = itor + 1
            twh_vtorsion = twh_vtorsion 
     &           + twh_get_torcoeff(torffnum,index)* 
     &           ( 1.0d0 - dcos(dble(itor)
     &           *(phi - twh_get_torcoeff(torffnum,0))))
         enddo
      elseif ( style .eq. 13 ) then
c         two term torsional potential used in W.L. Jorgensen;
c         C.J. Swenson; J. Am. Chem. Soc. 107, 569-578 (1985).
         phi = twh_arccos(phic)
         twh_vtorsion = twh_vtorsion 
     &        + twh_get_torcoeff(torffnum,1)*(1.0d0 - phic)
     &        + twh_get_torcoeff(torffnum,2)*(1.0d0 - dcos(2.0d0*phi))
      elseif ( style .eq. 14 ) then
c        --- UFF style 1 - cosine potential
c        --- energy is divided by the total number of torsions that 
c        --- occur across these two atoms
         unit_2 = twh_get_ijtor2(imolty,unit_1,torind)
         unit_3 = twh_get_ijtor3(imolty,unit_1,torind)
         totbond = dble((twh_get_invib(imolty,unit_2)-1)
     &        *(twh_get_invib(imolty,unit_3)-1))
         phi = twh_arccos(phic)
         twh_vtorsion = twh_get_torcoeff(torffnum,1)
     &        *( 1.0d0 - twh_get_torcoeff(torffnum,2)
     &        *dcos( twh_get_torcoeff(torffnum,3)*phi) )/totbond
      elseif ( style .eq. 15 ) then
c         --- DREIDING k[ 1 - Cos ( n (phi - phi_0))]
         phi = twh_arccos(phic)
         twh_vtorsion = twh_get_torcoeff(torffnum,1)*(1.0d0 
     &        - dcos(twh_get_torcoeff(torffnum,2)
     &        *(phi - twh_get_torcoeff(torffnum,3))))
c        --- energy is divided by the total number of torsions that 
c        --- occur across these two atoms
         unit_2 = twh_get_ijtor2(imolty,unit_1,torind)
         unit_3 = twh_get_ijtor3(imolty,unit_1,torind)
         totbond = dble((twh_get_invib(imolty,unit_2)-1)
     &        *(twh_get_invib(imolty,unit_3)-1))
         twh_vtorsion = twh_vtorsion / totbond
      elseif ( style .eq. 16 ) then
c        --- 2-fold cosine
         phi = twh_arccos(phic)
         twh_vtorsion = twh_get_torcoeff(torffnum,1)
     &        *(1.0d0 - dcos(2.0d0*phi))
      elseif ( style .eq. 17 ) then
c        --- Potoff torsions from J. Phys. Chem. B 108 14130-14136 2004
c        --- note that there is a typo in equation 5 of that paper
c        --- the correct functional form is
c        U = c1[1 + cos(phi + f1)] + c2[1 - cos^2(phi)]
         phi = twh_arccos(phic)
         twh_vtorsion = twh_get_torcoeff(torffnum,1)
     &        *(1.0d0 + dcos(phi + twh_get_torcoeff(torffnum,3)) )
     &        + twh_get_torcoeff(torffnum,2)*(1.0d0 - phic*phic)
      elseif ( style .eq. 18 ) then
c        --- square well torsions originally implemented for Elli2002
c        --- can't just use nonbond below as has different parameters
c        --- compute the 1-4 distance squared
         call twh_mimage(ibox,xcord(4),ycord(4),zcord(4)
     &        ,xcord(1),ycord(1),zcord(1),rxij,ryij,rzij)
         rijsq = rxij*rxij + ryij*ryij + rzij*rzij
         rij = dsqrt(rijsq)
         if ( twh_cmp_lt(rij, twh_get_torcoeff(torffnum,1)) .or.
     &        twh_cmp_gt(rij, twh_get_torcoeff(torffnum,2)) ) then
c           --- hard overlap
            loverlap = .true.
         else
c           --- assign default energy
            twh_vtorsion = twh_get_torcoeff(torffnum,3)
         endif
      elseif ( style .eq. 19 ) then
c        --- Gromos/Charmm/Amber cosine potential divided by the 
c        --- total number of torsions sharing these central atoms
         phi = twh_arccos(phic)
         do itor = 1,twh_get_ntorloop(torffnum)
            index = 3*(itor-1)
            twh_vtorsion = twh_vtorsion 
     &           + twh_get_torcoeff(torffnum,index+1)* 
     &           ( 1.0d0 + dcos(twh_get_torcoeff(torffnum,index+2)*phi
     &           - twh_get_torcoeff(torffnum,index+3))
     &           )
         enddo
c        --- divide out by the total number of torsions on these atoms
         unit_2 = twh_get_ijtor2(imolty,unit_1,torind)
         unit_3 = twh_get_ijtor3(imolty,unit_1,torind)
         totbond = dble((twh_get_invib(imolty,unit_2)-1)
     &        *(twh_get_invib(imolty,unit_3)-1))
         twh_vtorsion = twh_vtorsion / totbond
      elseif ( style .eq. 20  ) then
c        --- OPLS-AA Fluorocarbon 4 parameter cosine series 
c        --- used in the OPLS-AA force field.
c        --- J. Phys. Chem. A. v105, 4118-4125 (2001). 
         phi = twh_arccos(phic)
         twh_vtorsion = twh_get_torcoeff(torffnum,1)*(1.d0 + phic) +
     &        twh_get_torcoeff(torffnum,2)*(1.d0 - dcos(2.0d0*phi)) +
     &        twh_get_torcoeff(torffnum,3)*(1.d0 + dcos(3.0d0*phi)) +
     &        twh_get_torcoeff(torffnum,4)*(1.d0 - dcos(4.0d0*phi))
      elseif ( style .eq. 21 ) then
c        --- multiple rigid dihedral
         loverlap = .true.
         phi = twh_arccos(phic)
         do itor = 1,twh_get_ntorloop(torffnum)
            if ( dabs( twh_get_torcoeff(torffnum,itor) - phi) .lt. 
     &           twh_get_torcoeff(torffnum,0) ) then
               loverlap = .false.
               twh_vtorsion = 0.0d0
            endif
         enddo
         call twh_constant_pi(GLB_GET,onepi)
      elseif ( style .eq. 22 ) then
c        --- fluoroalkane from S.T. Cui; J.I. Siepmann; H.D. Cochran;
c        P.T. Cummings; "Intermolecular potentials and vapor-liquid
c        phase equilibria of perfluorinated alkanes",
c        Fluid Phase Equilib. 146 51-61 (1998)
         phi = twh_arccos(phic)
         twh_vtorsion = twh_get_torcoeff(torffnum,0)
     &        + twh_get_torcoeff(torffnum,1)*(1.0d0 - phic)
     &        + twh_get_torcoeff(torffnum,2)*(1.0d0 - dcos(3.0d0*phi))
     &        + twh_get_torcoeff(torffnum,3)*(1.0d0 - phic)**5
     &        + twh_get_torcoeff(torffnum,4)
     &        *twh_expon(-twh_get_torcoeff(torffnum,5)*phi**2)
      elseif ( style .eq. 23 ) then
c        --- simple two term cosine originally from 
c        --- K. Drukker; G. Wu; G. C. Schatz;
c        --- "Model simulations of DNA denaturation dynamics"
c        --- J. Chem. Phys. 114 579-590 (2001)
         phi = twh_arccos(phi)
         twh_vtorsion = twh_get_torcoeff(torffnum,1)
     &        *(1.0d0 - dcos(phi - twh_get_torcoeff(torffnum,0)))
      else
         write(6,*) 'VTORSION: non-defined torsional type'
         lfailure = .true.
         return
      endif

c     --- compute the 1-4 nonbonded terms, if they are used
      call twh_torofcode(GLB_GET,imolty,unit_1,torind,torofcode)
      if ( torofcode .ne. CODE_ZERO ) then
         if ( cordstyle .ne. 1 ) then
            write(6,*) 'cannot pass angle to VTORSION with 1-4'
            lfailure = .true.
            return
         endif
c        --- compute 1-4 vdw interactions
c        --- compute the 1-4 distance squared
         call twh_mimage(ibox,xcord(4),ycord(4),zcord(4)
     &        ,xcord(1),ycord(1),zcord(1),rxij,ryij,rzij)
         rijsq = rxij*rxij + ryij*ryij + rzij*rzij

c        --- need unit number of the opposite atom
         unit_4 = twh_get_ijtor4(imolty,unit_1,torind)
         call twh_ntype(GLB_GET,imolty,unit_1,ntii)
         call twh_ntype(GLB_GET,imolty,unit_4,ntjj)

         vonefour = twh_vtwobody(lfailure,rijsq,ntii,ntjj,3,loverlap
     &        ,vf,dv )
         if ( torofcode .eq. CODE_RESCALE ) then
c           --- get the scaling factor and rescale
            call twh_itscale(GLB_GET,imolty,unit_1,torind,torscale)
            vonefour = torscale*vonefour
         endif
c        --- add into the torsion energy
         twh_vtorsion = twh_vtorsion + vonefour
         if ( lfailure ) return
         if ( loverlap ) return

c        --- compute 1-4 coulombic interactions
         rij = dsqrt(rijsq)

c        --- note that we are only computing the real space portion
c        --- of the coulombic interaction here.  The correction for 
c        --- scaling factors due to the ewald sum is accounted for
c        --- in the eng** routines.  It is done this way for two reasons
c        --- 1) So we don't need any error function info in this step
c        --- 2) So I don't need any special logic to determine 1-4 
c        --- interactions in engatom.

         call twh_constant_qqfact(GLB_GET,qqfact)
         vonefour = twh_get_torcoeff(torffnum,-1)
     &        *qqfact*twh_vcoulomb(ibox,rij,imolty,unit_1,imolty
     &        ,unit_4,.true.,vf,dv,loverlap)
         if ( torofcode .eq. CODE_RESCALE ) then
            vonefour = torscale*vonefour
         endif
         twh_vtorsion = twh_vtorsion + vonefour
         if ( loverlap ) return
      endif

c      write(6,*) 'finished VTORSION'

      return
      end
