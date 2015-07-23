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
      subroutine twh_polyc19eef1(lfailure,imolty,iatom,termflag,total
     &,ibondend,ipeptide,iha,icb,pepname,atomname,aminoname,shortname)
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

c      write(6,*) 'start POLYC19EEF1'

      lsettorsion = .true.
      tortype = 0

c     --- determine which amino acid is built this step
c     --- check the special cases first, otherwise this is a generic
c     --- amino acid on a generic backbone.

      if ( pepname(ipeptide) .eq. 'p0' .or.
     &     pepname(ipeptide) .eq. 'P0' ) then
         call twh_ua_pepname_p0(lfailure,aminoname,shortname,termflag
     &        ,imolty,total,iatom,icb,iha,ibondend,tortype,lsettorsion)
         if ( lfailure ) return
         if ( termflag .eq. -1 ) then
c           --- N-terminal proline
            write(6,*) 'BUILDPROT: Charmm19 does not have a set'
            write(6,*) 'of parameters for N-terminal proline'
            lfailure = .true.
            return
         elseif ( termflag .eq. 0 ) then
c           --- interior proline
            atomname(iatom) = 'N'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.20d0)
            atomname(iatom+1) = 'CH1E'
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.10d0)
            atomname(iatom+2) = 'C'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.55d0)
            atomname(iatom+3) = 'O'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.55d0)
            atomname(iatom+4) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,iatom+4,0.00d0)
            atomname(iatom+5) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,iatom+5,0.00d0)
            atomname(iatom+6) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,iatom+6,0.10d0)
c           --- improper torsions
c           --- N - Ca - Cd - prev C
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_inimprop(GLB_SET,imolty,iatom,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,iatom,itemp+1,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom,itemp+1,iatom+6)
            call twh_ijimprop4(GLB_SET,imolty,iatom,itemp+1,ibondend)
c           --- C - Ca - next N - O
            itemp = twh_get_inimprop(imolty,iatom+2)
            call twh_inimprop(GLB_SET,imolty,iatom+2,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,itemp+1,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,itemp+1
     &           ,iatom+total)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,itemp+1,iatom+3)
c           --- Ca - N - C - Cb
            itemp = twh_get_inimprop(imolty,iatom+1)
            call twh_inimprop(GLB_SET,imolty,iatom+1,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+1,1,iatom)
            call twh_ijimprop3(GLB_SET,imolty,iatom+1,1,iatom+2)
            call twh_ijimprop4(GLB_SET,imolty,iatom+1,1,iatom+4)
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
            atomname(hmain) = 'CH3E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.0d0)
c           --- no side chain torsions
c           --- no extra side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 'r+' .or.
     &           pepname(ipeptide) .eq. 'R+' ) then
            call twh_ua_pepname_r_plus(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.00d0)
            atomname(hmain+2) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.10d0)
            atomname(hmain+3) = 'NH1'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.40d0)
            atomname(hmain+4) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.25d0)
            atomname(hmain+5) = 'NC2'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.85d0)
            atomname(hmain+6) = 'NC2'
            call twh_qqatom(GLB_SET,imolty,hmain+6,-0.85d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.15d0)
            atomname(shift+1) = 'HC'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.40d0)
            atomname(shift+2) = 'HC'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.40d0)
            atomname(shift+3) = 'HC'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.40d0)
            atomname(shift+4) = 'HC'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.40d0)
c           --- side chain improper torsions
c           --- Ne - Cd - Cz - He
            itemp = twh_get_inimprop(imolty,hmain+3)
            call twh_inimprop(GLB_SET,imolty,hmain+3,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,itemp+1,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,itemp+1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,itemp+1,shift)
c           --- Cz - Nh1 - Nh2 - Ne
            itemp = twh_get_inimprop(imolty,hmain+4)
            call twh_inimprop(GLB_SET,imolty,hmain+4,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,itemp+1,hmain+5)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,itemp+1,hmain+6)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,itemp+1,hmain+3)
         elseif ( pepname(ipeptide) .eq. 'n0' .or.
     &           pepname(ipeptide) .eq. 'N0' ) then
            call twh_ua_pepname_n0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.55d0)
            atomname(hmain+2) = 'O'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.55d0)
            atomname(hmain+3) = 'NH2'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.60d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.30d0)
            atomname(shift+1) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.30d0)
c           --- side chain improper torsions
c           --- Cg - Od1 - Nd2 - Cb
            itemp = twh_get_inimprop(imolty,hmain+1)
            call twh_inimprop(GLB_SET,imolty,hmain+1,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp+1,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp+1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp+1,hmain)
c           --- Nd2 - Hd21 - Hd22 - Cg
            itemp = twh_get_inimprop(imolty,hmain+3)
            call twh_inimprop(GLB_SET,imolty,hmain+3,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,itemp+1,shift)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,itemp+1,shift+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,itemp+1,hmain+1)
         elseif ( pepname(ipeptide) .eq. 'd-' .or.
     &           pepname(ipeptide) .eq. 'D-' ) then
            call twh_ua_pepname_d_neg(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.15d0)
            atomname(hmain+1) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+1,1.35d0)
            atomname(hmain+2) = 'OC'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.60d0)
            atomname(hmain+3) = 'OC'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.60d0)
c           --- side chain improper torsion
c           --- Cg - Od1 - Od2 - Cb
            itemp = twh_get_inimprop(imolty,hmain+1)
            call twh_inimprop(GLB_SET,imolty,hmain+1,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp+1,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp+1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp+1,hmain)
         elseif ( pepname(ipeptide) .eq. 'c0' .or.
     &           pepname(ipeptide) .eq. 'C0' ) then
            call twh_ua_pepname_c0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.19d0)
            atomname(hmain+1) = 'SH1E'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.19d0)
         elseif ( pepname(ipeptide) .eq. 'q0' .or.
     &           pepname(ipeptide) .eq. 'Q0' ) then
            call twh_ua_pepname_q0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.00d0)
            atomname(hmain+2) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.55d0)
            atomname(hmain+3) = 'O'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.55d0)
            atomname(hmain+4) = 'NH2'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.60d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.30d0)
            atomname(shift+1) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.30d0)
c           --- side chain improper torsions
c           --- Cd - Oe1 - Ne2 - Cg
            itemp = twh_get_inimprop(imolty,hmain+2)
            call twh_inimprop(GLB_SET,imolty,hmain+2,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,itemp+1,hmain+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,itemp+1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,itemp+1,hmain+1)
c           --- Ne2 - He21 - He22 - Cd
            itemp = twh_get_inimprop(imolty,hmain+4)
            call twh_inimprop(GLB_SET,imolty,hmain+4,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,itemp+1,shift)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,itemp+1,shift+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,itemp+1,hmain+2)
         elseif ( pepname(ipeptide) .eq. 'e-' .or.
     &           pepname(ipeptide) .eq. 'E-' ) then
            call twh_ua_pepname_e_neg(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.15d0)
            atomname(hmain+2) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+2,1.35d0)
            atomname(hmain+3) = 'OC'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.60d0)
            atomname(hmain+4) = 'OC'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.60d0)
c           --- side chain improper torsion
c           --- Cd - Oe1 - Oe2 - Cg
            call twh_inimprop(GLB_SET,imolty,hmain+2,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,hmain+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,hmain+1)
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
            atomname(hmain) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.10d0)
            atomname(hmain+1) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.15d0)
            atomname(hmain+2) = 'NH1'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.30d0)
            atomname(hmain+3) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.20d0)
            atomname(hmain+4) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.45d0)
            atomname(hmain+5) = 'NH1'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.30d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.35d0)
            atomname(shift+1) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.35d0)
c           --- side chain true improper torsions
c           --- Cg - Nd1 - Cd2 - Cb
            itemp = twh_get_inimprop(imolty,hmain+1)
            call twh_inimprop(GLB_SET,imolty,hmain+1,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp+1,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp+1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp+1,hmain)
c           --- Nd1 - Cg - Ce1 - Hd1
            itemp = twh_get_inimprop(imolty,hmain+2)
            call twh_inimprop(GLB_SET,imolty,hmain+2,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,itemp+1,hmain+1)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,itemp+1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,itemp+1,shift)
c           --- Ne2 - Ce1 - Cd2 - He2
            itemp = twh_get_inimprop(imolty,hmain+5)
            call twh_inimprop(GLB_SET,imolty,hmain+5,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,itemp+1,hmain+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,itemp+1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,itemp+1,shift+1)
         elseif ( pepname(ipeptide) .eq. 'hd' .or.
     &           pepname(ipeptide) .eq. 'HD' .or.
     &           pepname(ipeptide) .eq. 'Hd' ) then
            call twh_ua_pepname_hd(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.10d0)
            atomname(hmain+2) = 'NH1'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.40d0)
            atomname(hmain+3) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.10d0)
            atomname(hmain+4) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.30d0)
            atomname(hmain+5) = 'NR'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.40d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.30d0)
c           --- side chain true improper torsions
c           --- Cg - Nd1 - Cd2 - Cb
            itemp = twh_get_inimprop(imolty,hmain+1)
            call twh_inimprop(GLB_SET,imolty,hmain+1,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp+1,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp+1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp+1,hmain)
c           --- Nd1 - Cg - Ce1 - Hd1
            itemp = twh_get_inimprop(imolty,hmain+2)
            call twh_inimprop(GLB_SET,imolty,hmain+2,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,itemp+1,hmain+1)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,itemp+1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,itemp+1,shift)
         elseif ( pepname(ipeptide) .eq. 'he' .or.
     &           pepname(ipeptide) .eq. 'HE' .or.
     &           pepname(ipeptide) .eq. 'He' ) then
            call twh_ua_pepname_he(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.10d0)
            atomname(hmain+2) = 'NR'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.40d0)
            atomname(hmain+3) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.10d0)
            atomname(hmain+4) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.30d0)
            atomname(hmain+5) = 'NH1'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.40d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.30d0)
c           --- side chain true improper torsions
c           --- Cg - Nd1 - Cd2 - Cb
            itemp = twh_get_inimprop(imolty,hmain+1)
            call twh_inimprop(GLB_SET,imolty,hmain+1,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp+1,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp+1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp+1,hmain)
c           --- Ne2 - Ce1 - Cd2 - He2
            itemp = twh_get_inimprop(imolty,hmain+5)
            call twh_inimprop(GLB_SET,imolty,hmain+5,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,itemp+1,hmain+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,itemp+1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,itemp+1,shift+1)
         elseif ( pepname(ipeptide) .eq. 'i0' .or.
     &           pepname(ipeptide) .eq. 'I0' ) then
            call twh_ua_pepname_i0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH1E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.00d0)
            atomname(hmain+2) = 'CH3E'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.00d0)
            atomname(hmain+3) = 'CH3E'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.00d0)
c           --- side chain improper torsions
c           --- Cb - Cg1 - Cg2 - Ca
            itemp = twh_get_inimprop(imolty,hmain)
            call twh_inimprop(GLB_SET,imolty,hmain,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain,itemp+1,hmain+1)
            call twh_ijimprop3(GLB_SET,imolty,hmain,itemp+1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain,itemp+1,iatom+1)
         elseif ( pepname(ipeptide) .eq. 'l0' .or.
     &           pepname(ipeptide) .eq. 'L0' ) then
            call twh_ua_pepname_l0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'CH1E'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.00d0)
            atomname(hmain+2) = 'CH3E'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.00d0)
            atomname(hmain+3) = 'CH3E'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.00d0)
c           --- regular torsions that are not communal
c           --- Ca - Cb - Cg - Cd2
            call twh_create_torsion(imolty,iatom+1,hmain,hmain+1
     &           ,hmain+3,tortype)
c           --- reverse: Cd2 - Cg - Cb - Ca
            call twh_create_torsion(imolty,hmain+3,hmain+1,hmain
     &           ,iatom+1,tortype)
c           --- side chain impropers
c           --- Cg - Cd2 - Cd1 - Cb
            itemp = twh_get_inimprop(imolty,hmain+1)
            call twh_inimprop(GLB_SET,imolty,hmain+1,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp+1,hmain+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp+1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp+1,hmain)
         elseif ( pepname(ipeptide) .eq. 'k+' .or.
     &           pepname(ipeptide) .eq. 'K+' ) then
            call twh_ua_pepname_k_plus(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.00d0)
            atomname(hmain+2) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.00d0)
            atomname(hmain+3) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.00d0)
            atomname(hmain+4) = 'NH3'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-1.35d0)
            atomname(shift) = 'HC'
            call twh_qqatom(GLB_SET,imolty,shift,0.45d0)
            atomname(shift+1) = 'HC'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.45d0)
            atomname(shift+2) = 'HC'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.45d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 'm0' .or.
     &           pepname(ipeptide) .eq. 'M0' ) then
            call twh_ua_pepname_m0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.06d0)
            atomname(hmain+2) = 'S'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.12d0)
            atomname(hmain+3) = 'CH3E'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.06d0)
c           --- no impropers on side chain
         elseif ( pepname(ipeptide) .eq. 'f0' .or.
     &           pepname(ipeptide) .eq. 'F0' ) then
            call twh_ua_pepname_f0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.0d0)
            atomname(hmain+2) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.00d0)
            atomname(hmain+3) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.00d0)
            atomname(hmain+4) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.00d0)
            atomname(hmain+5) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+5,0.00d0)
            atomname(hmain+6) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+6,0.00d0)
c           --- side chain improper torsions
c           --- Cg - Cd1 - Cd2 - Cb
            itemp = twh_get_inimprop(imolty,hmain+1)
            call twh_inimprop(GLB_SET,imolty,hmain+1,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp+1,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp+1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp+1,hmain)
         elseif ( pepname(ipeptide) .eq. 's0' .or.
     &           pepname(ipeptide) .eq. 'S0' ) then
            call twh_ua_pepname_s0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.25d0)
            atomname(hmain+1) = 'OH1'
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
            atomname(hmain) = 'CH1E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.25d0)
            atomname(hmain+1) = 'OH1'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.65d0)
            atomname(hmain+2) = 'CH3E'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.00d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.40d0)
c           --- side chain improper torsions
c           --- Cb - Og1 - Cg2 - Ca
            itemp = twh_get_inimprop(imolty,hmain)
            call twh_inimprop(GLB_SET,imolty,hmain,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain,itemp+1,hmain+1)
            call twh_ijimprop3(GLB_SET,imolty,hmain,itemp+1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain,itemp+1,iatom+1)
         elseif ( pepname(ipeptide) .eq. 'w0' .or.
     &           pepname(ipeptide) .eq. 'W0' ) then
            call twh_ua_pepname_w0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.03d0)
            atomname(hmain+2) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.06d0)
            atomname(hmain+3) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.10d0)
            atomname(hmain+4) = 'NH1'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.36d0)
            atomname(hmain+5) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.04d0)
            atomname(hmain+6) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+6,-0.03d0)
            atomname(hmain+7) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+7,0.00d0)
            atomname(hmain+8) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+8,0.00d0)
            atomname(hmain+9) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+9,0.00d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.30d0)
c           --- side chain improper torsions
c           --- Cg - Cd1 - Cd2 - Cb
            itemp = twh_get_inimprop(imolty,hmain+1)
            call twh_inimprop(GLB_SET,imolty,hmain+1,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp+1,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp+1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp+1,hmain)
c           --- Ne1 - Cd1 - Ce2 - He1
            itemp = twh_get_inimprop(imolty,hmain+4)
            call twh_inimprop(GLB_SET,imolty,hmain+4,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,itemp+1,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,itemp+1,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,itemp+1,shift)
         elseif ( pepname(ipeptide) .eq. 'y0' .or.
     &           pepname(ipeptide) .eq. 'Y0' ) then
            call twh_ua_pepname_y0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH2E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.00d0)
            atomname(hmain+2) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.00d0)
            atomname(hmain+3) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.00d0)
            atomname(hmain+4) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.00d0)
            atomname(hmain+5) = 'CR1E'
            call twh_qqatom(GLB_SET,imolty,hmain+5,0.00d0)
            atomname(hmain+6) = 'C'
            call twh_qqatom(GLB_SET,imolty,hmain+6,0.25d0)
            atomname(hmain+7) = 'OH1'
            call twh_qqatom(GLB_SET,imolty,hmain+7,-0.65d0)
            atomname(shift) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift,0.40d0)
c           --- side chain improper torsions
c           --- Cg - Cd1 - Cd2 - Cb
            itemp = twh_get_inimprop(imolty,hmain+1)
            call twh_inimprop(GLB_SET,imolty,hmain+1,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,itemp+1,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,itemp+1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,itemp+1,hmain)
c           --- Cz - Ce1 - Ce2 - Oh
            itemp = twh_get_inimprop(imolty,hmain+6)
            call twh_inimprop(GLB_SET,imolty,hmain+6,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+6,itemp+1,hmain+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+6,itemp+1,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+6,itemp+1,hmain+7)
         elseif ( pepname(ipeptide) .eq. 'v0' .or.
     &           pepname(ipeptide) .eq. 'V0' ) then
            call twh_ua_pepname_v0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &           ,lsettorsion)
            if ( lfailure ) return
            atomname(hmain) = 'CH1E'
            call twh_qqatom(GLB_SET,imolty,hmain,0.00d0)
            atomname(hmain+1) = 'CH3E'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.00d0)
            atomname(hmain+2) = 'CH3E'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.00d0)
c           --- side chain improper torsions
c           --- Cb - Cg2 - Cg1 - Ca
            itemp = twh_get_inimprop(imolty,hmain)
            call twh_inimprop(GLB_SET,imolty,hmain,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,hmain,itemp+1,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain,itemp+1,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain,itemp+1,iatom+1)
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
            atomname(iatom) =   'NH3'
            call twh_qqatom(GLB_SET,imolty,iatom,-1.35d0)
            atomname(iatom+1) = 'CH1E'
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.00d0)
            atomname(iatom+2) = 'C'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.51d0)
            atomname(iatom+3) = 'O'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.51d0)
c           --- N hydrogen
            atomname(hmain+heavy) = 'HC'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.45d0)
c           --- extra N-terminus hydrogens
            atomname(hmain+heavy+1) = 'HC'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.45d0)
            atomname(hmain+heavy+2) = 'HC'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.45d0)
            if ( heavy .eq. 0 ) then
               atomname(iatom+1) = 'CH2E'
            endif
c           --- additional regular torsions that are not communal
c           --- backbone 1H - N - Ca - C
            call twh_create_torsion(imolty,hmain+heavy,iatom,iatom+1
     &           ,iatom+2,tortype)
c           --- reverse backbone 1H - N - Ca - C
            call twh_create_torsion(imolty,iatom+2,iatom+1,iatom
     &           ,hmain+heavy,tortype)
c           --- backbone 2H - N - Ca - C
            call twh_create_torsion(imolty,hmain+heavy+1,iatom,iatom+1
     &           ,iatom+2,tortype)
c           --- reverse backbone 2H - N - Ca - C
            call twh_create_torsion(imolty,iatom+2,iatom+1,iatom
     &           ,hmain+heavy+1,tortype)
c           --- backbone 3H - N - Ca - C
            call twh_create_torsion(imolty,hmain+heavy+2,iatom,iatom+1
     &           ,iatom+2,tortype)
c           --- reverse backbone 3H - N - Ca - C
            call twh_create_torsion(imolty,iatom+2,iatom+1,iatom
     &           ,hmain+heavy+2,tortype)
c           --- improper torsions
c           --- backbone C-Ca-N-O
            itemp = twh_get_inimprop(imolty,iatom+2)
            call twh_inimprop(GLB_SET,imolty,iatom+2,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,itemp+1,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,itemp+1
     &           ,iatom+total)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,itemp+1,iatom+3)
            if ( heavy .ne. 0 ) then
c              --- stereochem sidechain Ca - N - C - Cb
               itemp = twh_get_inimprop(imolty,iatom+1)
               call twh_inimprop(GLB_SET,imolty,iatom+1,itemp+1)
               call twh_ijimprop2(GLB_SET,imolty,iatom+1,itemp+1,iatom)
               call twh_ijimprop3(GLB_SET,imolty,iatom+1,itemp+1
     &              ,iatom+2)
               call twh_ijimprop4(GLB_SET,imolty,iatom+1,itemp+1,hmain)
            endif
         elseif ( termflag .eq. 0 ) then
c           --- interior segment main chain
c           --- assign atom type and charges
            atomname(iatom) =   'NH1'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.35d0)
            atomname(iatom+1) = 'CH1E'
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.10d0)
            atomname(iatom+2) = 'C'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.55d0)
            atomname(iatom+3) = 'O'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.55d0)
            atomname(hmain+heavy) = 'H'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.25d0)
            if ( heavy .eq. 0 ) then
               atomname(iatom+1) = 'CH2E'
            endif
c           --- improper torsions
c           --- backbone N - prev C - Ca - Hn
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_inimprop(GLB_SET,imolty,iatom,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,iatom,itemp+1,ibondend)
            call twh_ijimprop3(GLB_SET,imolty,iatom,itemp+1,iatom+1)
            call twh_ijimprop4(GLB_SET,imolty,iatom,itemp+1,hmain+heavy)
c           --- backbone C - Ca - next N - O
            itemp = twh_get_inimprop(imolty,iatom+2)
            call twh_inimprop(GLB_SET,imolty,iatom+2,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,itemp+1,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,itemp+1
     &           ,iatom+total)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,itemp+1,iatom+3)
            if ( heavy .ne. 0 ) then
c              --- stereochem side chain Ca - N - C - Cb
               itemp = twh_get_inimprop(imolty,iatom+1)
               call twh_inimprop(GLB_SET,imolty,iatom+1,itemp+1)
               call twh_ijimprop2(GLB_SET,imolty,iatom+1,itemp+1,iatom)
               call twh_ijimprop3(GLB_SET,imolty,iatom+1,itemp+1
     &              ,iatom+2)
               call twh_ijimprop4(GLB_SET,imolty,iatom+1,itemp+1,hmain)
            endif
         elseif ( termflag .eq. 1 ) then
c           --- c-terminus main chain
c           --- atom types and charges
            atomname(iatom) =  'NH1'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.35d0)
            atomname(iatom+1) = 'CH1E'
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.10d0)
            atomname(iatom+2) = 'C'
            call twh_qqatom(GLB_SET,imolty,iatom+2,1.20d0)
            atomname(iatom+3) = 'OC'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.60d0)
c           --- c-terminus oxygen
            atomname(hmain+heavy) = 'OC'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.60d0)
c           --- nitrogen hydrogen
            atomname(hmain+heavy+1) = 'H'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.25d0)
            if ( heavy .eq. 0 ) then
               atomname(iatom+1) = 'CH2E'
            endif
c           --- improper torsions
c           --- backbone N - prev C - Ca - Hn
            itemp = twh_get_inimprop(imolty,iatom)
            call twh_inimprop(GLB_SET,imolty,iatom,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,iatom,itemp+1,ibondend)
            call twh_ijimprop3(GLB_SET,imolty,iatom,itemp+1,iatom+1)
            call twh_ijimprop4(GLB_SET,imolty,iatom,itemp+1
     &           ,hmain+heavy+1)
c           --- C-terminus C - Ca - O - O
            itemp = twh_get_inimprop(imolty,iatom+2)
            call twh_inimprop(GLB_SET,imolty,iatom+2,itemp+1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,itemp+1,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,itemp+1
     &           ,hmain+heavy)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,itemp+1,iatom+3)
            if ( heavy .ne. 0 ) then
c              --- stereochem side chain Ca - N - C - Cb
               itemp = twh_get_inimprop(imolty,iatom+1)
               call twh_inimprop(GLB_SET,imolty,iatom+1,itemp+1)
               call twh_ijimprop2(GLB_SET,imolty,iatom+1,itemp+1,iatom)
               call twh_ijimprop3(GLB_SET,imolty,iatom+1,itemp+1
     &              ,iatom+2)
               call twh_ijimprop4(GLB_SET,imolty,iatom+1,itemp+1,hmain)
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

c      write(6,*) 'finished POLYC19eef1'

      return
      end
#else
c     Compilation disabled (for speed) with FFLAGS=-DDISABLE_POLY
c     statement in ./configure
      write(6,*) "polyc19eef1 disabled (compile flag)"
      lfailure = .true.
      return
      end
#endif
