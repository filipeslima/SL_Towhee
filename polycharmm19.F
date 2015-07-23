#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2011 Marcus G. Martin                       *
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
      subroutine twh_polycharmm19(lfailure,imolty,iatom,termflag,total
     &,ibondend,ipeptide,iha,icb
     &,pepname,atomname,aminoname,shortname)
c     ******************************************************************
c     * builds the input file for an amino acid sequence starting from *
c     * the N-terminus                                                 *
c     * imolty: molecule type number                                   *
c     *                                                                *
c     * split out of buildprot 12-12-2002 by M.G. Martin               *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
#ifndef DISABLE_POLY
      implicit none
#define FUNCTION_GET_INIMPROP
#include "functions.h"
c     ---- variables passed to/from the subroutine
      logical lfailure
      integer imolty,iatom,termflag,total,ibondend,ipeptide
      integer iha,icb
      character*2 pepname
      dimension pepname(NUMAX)
      character*3 shortname
      character*10 atomname
      dimension atomname(NUMAX)
      character*40 aminoname

c     --- local variables
      logical lsettorsion
      integer shift,heavy,ninsert,hmain,light,itemp,tortype

c      write(6,*) 'start POLYCHARMM19'

c     --- need to explicitly set the torsions for this forcefield
      lsettorsion = .true.
c     --- all torsions we create here will have type 0
      tortype = 0

c     --- determine which amino acid is built this step
c     --- check the special cases first, otherwise this is a generic
c     --- amino acid on a generic backbone.

      call twh_ua_pepname_p0(lfailure,aminoname,shortname,termflag
     &     ,imolty,total,iatom,icb,iha,ibondend,tortype,lsettorsion)
      if ( lfailure ) return

      if ( pepname(ipeptide) .eq. 'p0' .or.
     &     pepname(ipeptide) .eq. 'P0' ) then
c        --- proline requires its own set of special cases
         aminoname = 'proline'
         shortname = 'PRO'
         if ( termflag .eq. -1 ) then
c           --- N-terminal proline
            write(6,*) 'BUILDPROT: Charmm19 does not have a set'
            write(6,*) 'of parameters for N-terminal proline'
            lfailure = .true.
            return
         elseif ( termflag .eq. 0 ) then
c           --- interior proline
c           --- improper torsions
c           --- N - Ca - Cd - prev C
            call twh_inimprop(GLB_INCR,imolty,iatom,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,iatom,itemp,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom,itemp,iatom+6)
            call twh_ijimprop4(GLB_SET,imolty,iatom,itemp,ibondend)
c           --- C - Ca - next N - O
            call twh_inimprop(GLB_INCR,imolty,iatom+2,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,itemp,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,itemp
     &           ,iatom+total)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,itemp,iatom+3)
c           --- Ca - N - C - Cb
            call twh_inimprop(GLB_INCR,imolty,iatom+1,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,iatom+1,itemp,iatom)
            call twh_ijimprop3(GLB_SET,imolty,iatom+1,itemp,iatom+2)
            call twh_ijimprop4(GLB_SET,imolty,iatom+1,itemp,iatom+4)
         elseif ( termflag .eq. 1 ) then
c           --- C-terminal proline
            write(6,*) 'BUILDPROT: Charmm19 does not have a set'
            write(6,*) 'of parameters for C-terminal proline'
            lfailure = .true.
            return
         endif
c        --- update ibondend for proline (same as generic amino acid)
         ibondend = iatom+2
      else
c        --- generic amino acid side chain on normal backbone
c        --- set up the side chain parts that are independent of
c        --- whether this is n-term, c-term, or interior
c        --- also include some non-standard amino acids that still
c        --- have the same backbone

         call twh_ua_pepname_generic(lfailure,termflag,hmain,iatom
     &        ,ninsert)
         if ( lfailure ) return

         if ( pepname(ipeptide) .eq. 'a0' .or.
     &        pepname(ipeptide) .eq. 'A0' ) then
            call twh_ua_pepname_a0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
c           --- assign atomtypes
            atomname(hmain) = 'C(3)'
            call twh_qqatom(GLB_SET,imolty,hmain,0.0d0)
c           --- no side chain torsions
c           --- no extra side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 'r+' .or.
     &           pepname(ipeptide) .eq. 'R+' ) then
            call twh_ua_pepname_r_plus(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.00d0)
            atomname(hmain+2) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.10d0)
            atomname(hmain+3) = 'N1'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.40d0)
            atomname(hmain+4) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.50d0)
            atomname(hmain+5) = 'NC2'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.45d0)
            atomname(hmain+6) = 'NC2'
            call twh_qqatom(GLB_SET,imolty,hmain+6,-0.45d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.30d0)
            atomname(shift+1) = 'HC'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.35d0)
            atomname(shift+2) = 'HC'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.35d0)
            atomname(shift+3) = 'HC'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.35d0)
            atomname(shift+4) = 'HC'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.35d0)
c           --- side chain improper torsions
c           --- Ne - Cd - Cz - He
            call twh_inimprop(GLB_INCR,imolty,hmain+3,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,itemp,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,itemp,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,itemp,shift)
c           --- Cz - Nh1 - Nh2 - Ne
            call twh_inimprop(GLB_INCR,imolty,hmain+4,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,itemp,hmain+5)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,itemp,hmain+6)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,itemp,hmain+3)
         elseif ( pepname(ipeptide) .eq. 'n0' .or.
     &           pepname(ipeptide) .eq. 'N0' ) then
            call twh_ua_pepname_n0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.55d0)
            atomname(hmain+2) = 'O'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.55d0)
            atomname(hmain+3) = 'N2'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.60d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.30d0)
            atomname(shift+1) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.30d0)
c           --- side chain improper torsions
c           --- Cg - Od1 - Nd2 - Cb
            call twh_inimprop(GLB_INCR,imolty,hmain+1,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp,hmain)
c           --- Nd2 - Hd21 - Hd22 - Cg
            call twh_inimprop(GLB_INCR,imolty,hmain+3,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,itemp,shift)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,itemp,shift+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,itemp,hmain+1)
         elseif ( pepname(ipeptide) .eq. 'd-' .or.
     &           pepname(ipeptide) .eq. 'D-' ) then
            call twh_ua_pepname_d_neg(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.16d0)
            atomname(hmain+1) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.36d0)
            atomname(hmain+2) = 'OC'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.60d0)
            atomname(hmain+3) = 'OC'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.60d0)
c           --- side chain improper torsion
c           --- Cg - Od1 - Od2 - Cb
            call twh_inimprop(GLB_INCR,imolty,hmain+1,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp,hmain)
         elseif ( pepname(ipeptide) .eq. 'c0' .or.
     &           pepname(ipeptide) .eq. 'C0' ) then
            call twh_ua_pepname_c0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain,0.19d0)
            atomname(hmain+1) = 'S'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.19d0)
         elseif ( pepname(ipeptide) .eq. 'q0' .or.
     &           pepname(ipeptide) .eq. 'Q0' ) then
            call twh_ua_pepname_q0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.00d0)
            atomname(hmain+2) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.55d0)
            atomname(hmain+3) = 'O'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.55d0)
            atomname(hmain+4) = 'N2'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.60d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.30d0)
            atomname(shift+1) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.30d0)
c           --- side chain improper torsions
c           --- Cd - Oe1 - Ne2 - Cg
            call twh_inimprop(GLB_INCR,imolty,hmain+2,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,itemp,hmain+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,itemp,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,itemp,hmain+1)
c           --- Ne2 - He21 - He22 - Cd
            call twh_inimprop(GLB_INCR,imolty,hmain+4,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,itemp,shift)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,itemp,shift+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,itemp,hmain+2)
         elseif ( pepname(ipeptide) .eq. 'e-' .or.
     &           pepname(ipeptide) .eq. 'E-' ) then
            call twh_ua_pepname_e_neg(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.16d0)
            atomname(hmain+2) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.36d0)
            atomname(hmain+3) = 'OC'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.60d0)
            atomname(hmain+4) = 'OC'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.60d0)
c           --- side chain improper torsion
c           --- Cd - Oe1 - Oe2 - Cg
            call twh_inimprop(GLB_INCR,imolty,hmain+2,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,itemp,hmain+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,itemp,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,itemp,hmain+1)
         elseif ( pepname(ipeptide) .eq. 'g0' .or.
     &           pepname(ipeptide) .eq. 'G0' ) then
            call twh_ua_pepname_g0(lfailure,aminoname,shortname
     &           ,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
         elseif ( pepname(ipeptide) .eq. 'h+' .or.
     &           pepname(ipeptide) .eq. 'H+' ) then
            call twh_ua_pepname_h_plus(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain,0.10d0)
            atomname(hmain+1) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.15d0)
            atomname(hmain+2) = 'N1'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.30d0)
            atomname(hmain+3) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.20d0)
            atomname(hmain+4) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.45d0)
            atomname(hmain+5) = 'N1'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.30d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.35d0)
            atomname(shift+1) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.35d0)
c           --- side chain true improper torsions
c           --- Cg - Nd1 - Cd2 - Cb
            call twh_inimprop(GLB_INCR,imolty,hmain+1,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp,hmain)
c           --- Nd1 - Cg - Ce1 - Hd1
            call twh_inimprop(GLB_INCR,imolty,hmain+2,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,itemp,hmain+1)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,itemp,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,itemp,shift)
c           --- Ne2 - Ce1 - Cd2 - He2
            call twh_inimprop(GLB_INCR,imolty,hmain+5,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,itemp,hmain+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,itemp,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,itemp,shift+1)
         elseif ( pepname(ipeptide) .eq. 'hd' .or.
     &           pepname(ipeptide) .eq. 'HD' .or.
     &           pepname(ipeptide) .eq. 'Hd' ) then
            call twh_ua_pepname_hd(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.10d0)
            atomname(hmain+2) = 'N1'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.40d0)
            atomname(hmain+3) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.10d0)
            atomname(hmain+4) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.30d0)
            atomname(hmain+5) = 'NR'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.40d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.30d0)
c           --- side chain true improper torsions
c           --- Cg - Nd1 - Cd2 - Cb
            call twh_inimprop(GLB_INCR,imolty,hmain+1,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp,hmain)
c           --- Nd1 - Cg - Ce1 - Hd1
            call twh_inimprop(GLB_INCR,imolty,hmain+2,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,itemp,hmain+1)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,itemp,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,itemp,shift)
         elseif ( pepname(ipeptide) .eq. 'he' .or.
     &           pepname(ipeptide) .eq. 'HE' .or.
     &           pepname(ipeptide) .eq. 'He' ) then
            call twh_ua_pepname_he(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.10d0)
            atomname(hmain+2) = 'NR'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.40d0)
            atomname(hmain+3) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.10d0)
            atomname(hmain+4) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.30d0)
            atomname(hmain+5) = 'N1'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.40d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.30d0)
c           --- side chain true improper torsions
c           --- Cg - Nd1 - Cd2 - Cb
            call twh_inimprop(GLB_INCR,imolty,hmain+1,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp,hmain)
c           --- Ne2 - Ce1 - Cd2 - He2
            call twh_inimprop(GLB_INCR,imolty,hmain+5,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,itemp,hmain+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,itemp,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,itemp,shift+1)
         elseif ( pepname(ipeptide) .eq. 'i0' .or.
     &           pepname(ipeptide) .eq. 'I0' ) then
            call twh_ua_pepname_i0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.00d0)
            atomname(hmain+2) = 'C(3)'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.00d0)
            atomname(hmain+3) = 'C(3)'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.00d0)
c           --- side chain improper torsions
c           --- Cb - Cg1 - Cg2 - Ca
            call twh_inimprop(GLB_INCR,imolty,hmain,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain,itemp,hmain+1)
            call twh_ijimprop3(GLB_SET,imolty,hmain,itemp,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain,itemp,iatom+1)
         elseif ( pepname(ipeptide) .eq. 'l0' .or.
     &           pepname(ipeptide) .eq. 'L0' ) then
            call twh_ua_pepname_l0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'CH'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.00d0)
            atomname(hmain+2) = 'C(3)'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.00d0)
            atomname(hmain+3) = 'C(3)'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.00d0)
c           --- additional regular torsions that are not communal
c           --- Ca - Cb - Cg - Cd1
            call twh_create_torsion(imolty,iatom+1,hmain,hmain+1
     &           ,hmain+2,tortype)
c           --- reverse: Cd1 - Cg - Cb - Ca
            call twh_create_torsion(imolty,hmain+2,hmain+1,hmain
     &           ,iatom+1,tortype)
c           --- side chain impropers
c           --- Cg - Cd2 - Cd1 - Cb
            call twh_inimprop(GLB_INCR,imolty,hmain+1,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp,hmain+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp,hmain)
         elseif ( pepname(ipeptide) .eq. 'k+' .or.
     &           pepname(ipeptide) .eq. 'K+' ) then
            call twh_ua_pepname_k_plus(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.00d0)
            atomname(hmain+2) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.00d0)
            atomname(hmain+3) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.25d0)
            atomname(hmain+4) = 'N3'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.30d0)
            atomname(shift) = 'HC'
            call twh_qqatom(GLB_SET,imolty,shift,0.35d0)
            atomname(shift+1) = 'HC'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.35d0)
            atomname(shift+2) = 'HC'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.35d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 'm0' .or.
     &           pepname(ipeptide) .eq. 'M0' ) then
            call twh_ua_pepname_m0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.06d0)
            atomname(hmain+2) = 'S'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.12d0)
            atomname(hmain+3) = 'C(3)'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.06d0)
c           --- no impropers on side chain
         elseif ( pepname(ipeptide) .eq. 'f0' .or.
     &           pepname(ipeptide) .eq. 'F0' ) then
            call twh_ua_pepname_f0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.0d0)
            atomname(hmain+2) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.00d0)
            atomname(hmain+3) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.00d0)
            atomname(hmain+4) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.00d0)
            atomname(hmain+5) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+5,0.00d0)
            atomname(hmain+6) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+6,0.00d0)
c           --- side chain improper torsions
c           --- Cg - Cd1 - Cd2 - Cb
            call twh_inimprop(GLB_INCR,imolty,hmain+1,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp,hmain)
         elseif ( pepname(ipeptide) .eq. 's0' .or.
     &           pepname(ipeptide) .eq. 'S0' ) then
            call twh_ua_pepname_s0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain,0.25d0)
            atomname(hmain+1) = 'OH'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.65d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.40d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 't0' .or.
     &           pepname(ipeptide) .eq. 'T0' ) then
            call twh_ua_pepname_t0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH'
            call twh_qqatom(GLB_SET,imolty,hmain,0.25d0)
            atomname(hmain+1) = 'OH'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.65d0)
            atomname(hmain+2) = 'C(3)'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.00d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.40d0)
c           --- side chain improper torsions
c           --- Cb - Og1 - Cg2 - Ca
            call twh_inimprop(GLB_INCR,imolty,hmain,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain,itemp,hmain+1)
            call twh_ijimprop3(GLB_SET,imolty,hmain,itemp,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain,itemp,iatom+1)
         elseif ( pepname(ipeptide) .eq. 'w0' .or.
     &           pepname(ipeptide) .eq. 'W0' ) then
            call twh_ua_pepname_w0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.03d0)
            atomname(hmain+2) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.06d0)
            atomname(hmain+3) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.10d0)
            atomname(hmain+4) = 'N1'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.36d0)
            atomname(hmain+5) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.04d0)
            atomname(hmain+6) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+6,-0.03d0)
            atomname(hmain+7) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+7,0.00d0)
            atomname(hmain+8) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+8,0.00d0)
            atomname(hmain+9) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+9,0.00d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.30d0)
c           --- side chain improper torsions
c           --- Cg - Cd1 - Cd2 - Cb
            call twh_inimprop(GLB_INCR,imolty,hmain+1,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp,hmain)
c           --- Ne1 - Cd1 - Ce2 - He1
            call twh_inimprop(GLB_INCR,imolty,hmain+4,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,itemp,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,itemp,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,itemp,shift)
         elseif ( pepname(ipeptide) .eq. 'y0' .or.
     &           pepname(ipeptide) .eq. 'Y0' ) then
            call twh_ua_pepname_y0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'C(2)'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.00d0)
            atomname(hmain+2) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.00d0)
            atomname(hmain+3) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.00d0)
            atomname(hmain+4) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.00d0)
            atomname(hmain+5) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+5,0.00d0)
            atomname(hmain+6) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+6,0.25d0)
            atomname(hmain+7) = 'OH'
            call twh_qqatom(GLB_SET,imolty,hmain+7,-0.65d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.40d0)
c           --- side chain improper torsions
c           --- Cg - Cd1 - Cd2 - Cb
            call twh_inimprop(GLB_INCR,imolty,hmain+1,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp,hmain)
c           --- Cz - Ce1 - Ce2 - Oh
            call twh_inimprop(GLB_INCR,imolty,hmain+6,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain+6,itemp,hmain+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+6,itemp,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+6,itemp,hmain+7)
         elseif ( pepname(ipeptide) .eq. 'v0' .or.
     &           pepname(ipeptide) .eq. 'V0' ) then
            call twh_ua_pepname_v0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'C(3)'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.00d0)
            atomname(hmain+2) = 'C(3)'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.00d0)
c           --- side chain improper torsions
c           --- Cb - Cg2 - Cg1 - Ca
            call twh_inimprop(GLB_INCR,imolty,hmain,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,hmain,itemp,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain,itemp,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain,itemp,iatom+1)
         else
            write(6,*) 'invalid peptide name ',ipeptide
     &           ,' ',pepname(ipeptide)
            lfailure = .true.
            return
         endif

         call twh_ua_pepname_main_chain(lfailure,termflag
     &        ,imolty,total,iatom,ibondend,heavy,light,ninsert
     &        ,hmain,lsettorsion,tortype)
         if ( lfailure ) return

         if ( termflag .eq. -1 ) then
c           --- n-terminus main chain
c           --- atom types and charges
            atomname(iatom) =   'N3'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.30d0)
            atomname(iatom+1) = 'CH'
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.25d0)
            atomname(iatom+2) = 'C'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.51d0)
            atomname(iatom+3) = 'O'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.51d0)
c           --- N hydrogen
            atomname(hmain+heavy) = 'HC'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.35d0)
c           --- extra N-terminus hydrogens
            atomname(hmain+heavy+1) = 'HC'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.35d0)
            atomname(hmain+heavy+2) = 'HC'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.35d0)
            if ( heavy .eq. 0 ) then
c              --- rename the main chain
               atomname(iatom+1) = 'C(2)'
            endif
c           --- improper torsions
c           --- backbone C-Ca-N-O
            call twh_inimprop(GLB_INCR,imolty,iatom+2,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,itemp,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,itemp,iatom+total)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,itemp,iatom+3)
            if ( heavy .ne. 0 ) then
c              --- stereochem sidechain Ca - N - C - Cb
               call twh_inimprop(GLB_INCR,imolty,iatom+1,1)
               itemp = twh_get_inimprop(imolty,iatom)
               call twh_ijimprop2(GLB_SET,imolty,iatom+1,itemp,iatom)
               call twh_ijimprop3(GLB_SET,imolty,iatom+1,itemp,iatom+2)
               call twh_ijimprop4(GLB_SET,imolty,iatom+1,itemp,hmain)
            endif
         elseif ( termflag .eq. 0 ) then
c           --- interior segment main chain
c           --- assign atom type and charges
            atomname(iatom) =   'N1'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.35d0)
            atomname(iatom+1) = 'CH'
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.10d0)
            atomname(iatom+2) = 'C'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.55d0)
            atomname(iatom+3) = 'O'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.55d0)
            atomname(hmain+heavy) = 'H'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.25d0)
            if ( heavy .eq. 0 ) then
c              --- rename the main chain
               atomname(iatom+1) = 'C(2)'
            endif
c           --- improper torsions
c           --- backbone N - prev C - Ca - Hn
            call twh_inimprop(GLB_INCR,imolty,iatom,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,iatom,itemp,ibondend)
            call twh_ijimprop3(GLB_SET,imolty,iatom,itemp,iatom+1)
            call twh_ijimprop4(GLB_SET,imolty,iatom,itemp,hmain+heavy)
c           --- backbone C - Ca - next N - O
            call twh_inimprop(GLB_INCR,imolty,iatom+2,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,itemp,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,itemp,iatom+total)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,itemp,iatom+3)
            if ( heavy .ne. 0 ) then
c              --- stereochem side chain Ca - N - C - Cb
               call twh_inimprop(GLB_INCR,imolty,iatom+1,1)
               itemp = twh_get_inimprop(imolty,iatom)
               call twh_ijimprop2(GLB_SET,imolty,iatom+1,itemp,iatom)
               call twh_ijimprop3(GLB_SET,imolty,iatom+1,itemp,iatom+2)
               call twh_ijimprop4(GLB_SET,imolty,iatom+1,itemp,hmain)
            endif
         elseif ( termflag .eq. 1 ) then
c           --- c-terminus main chain
c           --- atom types and charges
            atomname(iatom) =  'N1'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.35d0)
            atomname(iatom+1) = 'CH'
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.10d0)
            atomname(iatom+2) = 'C'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.14d0)
            atomname(iatom+3) = 'OC'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57d0)
c           --- c-terminus oxygen
            atomname(hmain+heavy) = 'OC'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.57d0)
c           --- nitrogen hydrogen
            atomname(hmain+heavy+1) = 'H'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.25d0)
            if ( heavy .eq. 0 ) then
c              --- rename the main chain
               atomname(iatom+1) = 'C(2)'
            endif
c           --- improper torsions
c           --- backbone N - prev C - Ca - Hn
            call twh_inimprop(GLB_INCR,imolty,iatom,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,iatom,itemp,ibondend)
            call twh_ijimprop3(GLB_SET,imolty,iatom,itemp,iatom+1)
            call twh_ijimprop4(GLB_SET,imolty,iatom,itemp,hmain+heavy+1)
c           --- C-terminus C - Ca - O - O
            call twh_inimprop(GLB_INCR,imolty,iatom+2,1)
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,itemp,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,itemp,hmain+heavy)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,itemp,iatom+3)
            if ( heavy .ne. 0 ) then
c              --- stereochem side chain Ca - N - C - Cb
               call twh_inimprop(GLB_INCR,imolty,iatom+1,1)
               itemp = twh_get_inimprop(imolty,iatom)
               call twh_ijimprop2(GLB_SET,imolty,iatom+1,itemp,iatom)
               call twh_ijimprop3(GLB_SET,imolty,iatom+1,itemp,iatom+2)
               call twh_ijimprop4(GLB_SET,imolty,iatom+1,itemp,hmain)
            endif
         else
c           --- invalid amino acid designation
            write(6,*) 'imolty,ipeptide,pepname',imolty,ipeptide
     &           ,pepname(ipeptide),' is not valid'
            lfailure = .true.
            return
         endif
c        --- update ibondend for the typical amino acids
         ibondend = iatom+2
      endif

c      write(6,*) 'finished POLYCHARMM19'

      return
      end
#else
c     Compilation disabled (for speed) with FFLAGS=-DDISABLE_POLY
c     statement in ./configure
      write(6,*) "polycharmm19 disabled (compile flag)"
      lfailure = .true.
      return
      end
#endif
