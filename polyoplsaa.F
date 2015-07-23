#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2011 Marcus G. Martin                       *
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
      subroutine twh_polyoplsaa(lfailure,imolty,iatom,termflag,total
     &     ,ibondend,ipeptide,ndisulfide,myaminonum,myatomnum,iha,icb
     &     ,pepname,atomname,aminoname,shortname)
c     ******************************************************************
c     * builds the input file for an amino acid sequence starting from *
c     * the N-terminus for the OPLS-aa force field                     *
c     *                                                                *
c     * Adapted from the Charmm27 builder by A.A. Chen 9-15-2005       *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
#ifndef DISABLE_POLY
      implicit none
c     ---- variables passed to/from the subroutine
      logical lfailure
      integer imolty,iatom,termflag,total,ibondend,ipeptide,ndisulfide
      integer iha,icb
      character*2 pepname
      dimension pepname(NUMAX)
      character*3 shortname
      character*10 atomname
      dimension atomname(NUMAX)
      character*40 aminoname
      integer myaminonum,myatomnum
      dimension myaminonum(NUMAX),myatomnum(NUMAX)

c     --- local variables
      integer shift,heavy,ninsert,hmain,light
      double precision charge

c      write(6,*) 'start POLYOPLSAA'

c     --- determine which amino acid is built this step
c     --- check the special cases first, otherwise this is a generic
c     --- amino acid on a generic backbone.
      if ( pepname(ipeptide) .eq. 'za' ) then
         call twh_aa_pepname_za(lfailure,aminoname,shortname,termflag
     &        ,imolty,total,iatom)
         if ( lfailure ) return
         atomname(iatom) = 'C n'
         call twh_qqatom(GLB_SET,imolty,iatom,0.50d0)
         atomname(iatom+1) = 'CT'
         call twh_qqatom(GLB_SET,imolty,iatom+1,-0.18d0)
         atomname(iatom+2) = 'O'
         call twh_qqatom(GLB_SET,imolty,iatom+2,-0.50d0)
         atomname(iatom+3) = 'HC~1'
         call twh_qqatom(GLB_SET,imolty,iatom+3,0.06d0)
         atomname(iatom+4) = 'HC~1'
         call twh_qqatom(GLB_SET,imolty,iatom+4,0.06d0)
         atomname(iatom+5) = 'HC~1'
         call twh_qqatom(GLB_SET,imolty,iatom+5,0.06d0)
c        --- improper torsion CY-CAY-N(next)-OY
c        --- OPLS type O-Cn-X-Y = type 2.3 for CY
         call twh_inimprop(GLB_SET,imolty,iatom,1)
         call twh_ijimprop2(GLB_SET,imolty,iatom,1,iatom+2)
         call twh_ijimprop3(GLB_SET,imolty,iatom,1,iatom+total)
         call twh_ijimprop4(GLB_SET,imolty,iatom,1,iatom+1)
c        --- update ibondend
         ibondend = iatom
      elseif ( pepname(ipeptide) .eq. 'zm' ) then
         call twh_aa_pepname_zm(lfailure,aminoname,shortname,termflag
     &        ,imolty,total,iatom,ibondend)
         if ( lfailure ) return
         atomname(iatom) = 'N'
         call twh_qqatom(GLB_SET,imolty,iatom,-0.76d0)
         atomname(iatom+1) = 'H'
         call twh_qqatom(GLB_SET,imolty,iatom+1,0.38d0)
         atomname(iatom+2) = 'H'
         call twh_qqatom(GLB_SET,imolty,iatom+2,0.38d0)
c        --- also need to modify the charges and atom types 
c        --- on the previous amino acid backbone
c        --- previous main chain C
         atomname(ibondend) = 'C n'
         call twh_qqatom(GLB_SET,imolty,ibondend,0.50d0)
c        --- previous main chain O
c        --- added atomname to be sure AAC 9/05
         atomname(ibondend+1) = 'O'
         call twh_qqatom(GLB_SET,imolty,ibondend+1,-0.50d0)
c        --- improper torsion on NT-C-HT1-HT2
c        --- OPLS type Z-N-X-Y = type 3.1
         call twh_inimprop(GLB_SET,imolty,iatom,1)
         call twh_ijimprop2(GLB_SET,imolty,iatom,1,iatom+2)
         call twh_ijimprop3(GLB_SET,imolty,iatom,1,ibondend)
         call twh_ijimprop4(GLB_SET,imolty,iatom,1,iatom+1)
c        --- ibondend not relevent for c-terminus
         ibondend = -1
      elseif ( pepname(ipeptide) .eq. 'zn' ) then
         call twh_aa_pepname_zn(lfailure,aminoname,shortname,termflag
     &        ,imolty,total,iatom,ibondend)
         if ( lfailure ) return
         atomname(iatom) = 'N'
         call twh_qqatom(GLB_SET,imolty,iatom,-0.50d0)
         atomname(iatom+1) = 'CT'
         call twh_qqatom(GLB_SET,imolty,iatom+1,0.02d0)
         atomname(iatom+2) = 'H'
         call twh_qqatom(GLB_SET,imolty,iatom+2,0.30d0)
         atomname(iatom+3) = 'HC~1'
         call twh_qqatom(GLB_SET,imolty,iatom+3,0.06d0)
         atomname(iatom+4) = 'HC~1'
         call twh_qqatom(GLB_SET,imolty,iatom+4,0.06d0)
         atomname(iatom+5) = 'HC~1'
         call twh_qqatom(GLB_SET,imolty,iatom+5,0.06d0)
c        --- improper torsion on "backbone" NT-C-CAT-HNT
c        --- OPLS type Z-N-X-Y = type 3.1
         call twh_inimprop(GLB_SET,imolty,iatom,1)
         call twh_ijimprop2(GLB_SET,imolty,iatom,1,iatom+2)
         call twh_ijimprop3(GLB_SET,imolty,iatom,1,iatom+1)
         call twh_ijimprop4(GLB_SET,imolty,iatom,1,ibondend)
c        --- ibondend not relevent for c-terminus
         ibondend = -1
      elseif ( pepname(ipeptide) .eq. 'zf' ) then
         call twh_aa_pepname_zf(lfailure,aminoname,shortname,termflag
     &        ,imolty,total,iatom)
         if ( lfailure ) return
         atomname(iatom) = 'C n'
         call twh_qqatom(GLB_SET,imolty,iatom,0.50d0)
         atomname(iatom+1) = 'O'
         call twh_qqatom(GLB_SET,imolty,iatom+1,-0.50d0)
         atomname(iatom+2) = 'HC~3'
         call twh_qqatom(GLB_SET,imolty,iatom+2,0.00d0)
c        --- improper torsion
c        --- OPLS type O-C-X-Y = type 2.3 for C
         call twh_inimprop(GLB_SET,imolty,iatom,1)
         call twh_ijimprop2(GLB_SET,imolty,iatom,1,iatom+1)
c        --- changed from iatom +3 assumed type AAC
         call twh_ijimprop3(GLB_SET,imolty,iatom,1,iatom+total)
         call twh_ijimprop4(GLB_SET,imolty,iatom,1,iatom+2)
c        --- update ibondend for formaldehyde
         ibondend = iatom
      elseif ( pepname(ipeptide) .eq. 'p0' .or.
     &        pepname(ipeptide) .eq. 'P0' ) then
         call twh_aa_pepname_p0(lfailure,aminoname,shortname,termflag
     &        ,imolty,total,iatom,icb,iha,ibondend)
         if ( lfailure ) return
         if ( termflag .eq. -1 ) then
c           --- N-terminal proline
            atomname(iatom) = 'NT+2'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.24d0)
            atomname(iatom+1) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.25d0)
            atomname(iatom+2) = 'C n'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.50d0)
            atomname(iatom+3) = 'O'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.50d0)
            atomname(iatom+4) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+4,-0.12d0)
            atomname(iatom+5) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+5,-0.12d0)
            atomname(iatom+6) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+6,0.19d0)
            atomname(iatom+7) = 'H'
            call twh_qqatom(GLB_SET,imolty,iatom+7,0.31d0)
            atomname(iatom+8) = 'H'
            call twh_qqatom(GLB_SET,imolty,iatom+8,0.31d0)
            atomname(iatom+9) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+9,0.06d0)
            atomname(iatom+10) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+10,0.06d0)
            atomname(iatom+11) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+11,0.06d0)
            atomname(iatom+12) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+12,0.06d0)
            atomname(iatom+13) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+13,0.06d0)
            atomname(iatom+14) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+14,0.06d0)
            atomname(iatom+15) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+15,0.06d0)
c           --- improper torsions
c           --- OPLS type O-Cn-X-Y = type 2.3 for C
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+3)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,iatom+total)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,iatom+1)
         elseif ( termflag .eq. 0 ) then
c           --- interior proline
            atomname(iatom) = 'N'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.14d0)
            atomname(iatom+1) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.01d0)
            atomname(iatom+2) = 'C n'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.50d0)
            atomname(iatom+3) = 'O'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.50d0)
            atomname(iatom+4) = 'CT '
            call twh_qqatom(GLB_SET,imolty,iatom+4,-0.12d0)
            atomname(iatom+5) = 'CT '
            call twh_qqatom(GLB_SET,imolty,iatom+5,-0.12d0)
            atomname(iatom+6) = 'CT '
            call twh_qqatom(GLB_SET,imolty,iatom+6,-0.05d0)
            atomname(iatom+7) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+7,0.06d0)
            atomname(iatom+8) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+8,0.06d0)
            atomname(iatom+9) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+9,0.06d0)
            atomname(iatom+10) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+10,0.06d0)
            atomname(iatom+11) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+11,0.06d0)
            atomname(iatom+12) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+12,0.06d0)
            atomname(iatom+13) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+13,0.06d0)
c           --- improper torsions
c           --- OPLS type Z-N-X-Y = type 3.1 for N
            call twh_inimprop(GLB_SET,imolty,iatom,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom,1,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom,1,ibondend)
            call twh_ijimprop4(GLB_SET,imolty,iatom,1,iatom+6)
c           --- OPLS type O-Cn-X-Y = type 2.3 for C
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+3)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,iatom+total)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,iatom+1)
         elseif ( termflag .eq. 1 ) then
c           --- C-terminal proline
            atomname(iatom) = 'N'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.14d0)
            atomname(iatom+1) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+1,-0.09d0)
            atomname(iatom+2) = 'C a'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.70d0)
            atomname(iatom+3) = 'O2'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.80d0)
            atomname(iatom+4) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+4,-0.12d0)
            atomname(iatom+5) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+5,-0.12d0)
            atomname(iatom+6) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+6,-0.05d0)
            atomname(iatom+7) = 'O2'
            call twh_qqatom(GLB_SET,imolty,iatom+7,-0.80d0)
            atomname(iatom+8) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+8,0.06d0)
            atomname(iatom+9) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+9,0.06d0)
            atomname(iatom+10) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+10,0.06d0)
            atomname(iatom+11) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+11,0.06d0)
            atomname(iatom+12) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+12,0.06d0)
            atomname(iatom+13) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+13,0.06d0)
            atomname(iatom+14) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,iatom+14,0.06d0)
c           --- improper torsions for proline backbone
c           --- OPLS type Z-N-X-Y = type 3.1 for N
            call twh_inimprop(GLB_SET,imolty,iatom,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom,1,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom,1,ibondend)
            call twh_ijimprop4(GLB_SET,imolty,iatom,1,iatom+6)
c           --- OPLS type O-Ca-X-Y = type 2.1 for C
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+3)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,iatom+7)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,iatom+1)
         endif
c        --- update ibondend for proline (same as generic amino acid)
         ibondend = iatom+2
      else
c        --- generic amino acid side chain on normal backbone
c        --- set up the side chain parts that are independent of
c        --- whether this is n-term, c-term, or interior
c        --- also include some non-standard amino acids that still
c        --- have the same backbone
         
         call twh_aa_pepname_generic(lfailure,termflag,hmain,iatom
     &        ,ninsert)
         if ( lfailure ) return

         if ( pepname(ipeptide) .eq. 'a0' .or.
     &        pepname(ipeptide) .eq. 'A0' ) then
            call twh_aa_pepname_a0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT  '
            call twh_qqatom(GLB_SET,imolty,hmain,-0.18d0)
            atomname(shift) = 'HC~1  '
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1  '
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'HC~1  '
            call twh_qqatom(GLB_SET,imolty,shift+2,0.06d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 'n0' .or.
     &           pepname(ipeptide) .eq. 'N0' ) then
            call twh_aa_pepname_n0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.12d0)
            atomname(hmain+1) = 'C n'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.50d0)
            atomname(hmain+2) = 'O'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.50d0)
            atomname(hmain+3) = 'N'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.76d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.38d0)
            atomname(shift+3) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.38d0)
c           --- side chain improper torsions (if any)
c           --- OPLS type O-Cn-X-Y
            call twh_inimprop(GLB_SET,imolty,hmain+1,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain)
c           --- OPLS type Z-N-X-Y
            call twh_inimprop(GLB_SET,imolty,hmain+3,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,1,hmain+1)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,1,shift+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,1,shift+3)
         elseif ( pepname(ipeptide) .eq. 'c0' .or.
     &           pepname(ipeptide) .eq. 'C0' ) then
            call twh_aa_pepname_c0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,0.06d0)
            atomname(hmain+1) = 'SH'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.435d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'HS'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.255d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 'cs' .or.
     &           pepname(ipeptide) .eq. 'CS' .or.
     &           pepname(ipeptide) .eq. 'cS' .or.
     &           pepname(ipeptide) .eq. 'Cs' ) then
            call twh_aa_pepname_cs(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
c           --- setup information to figure out the disulfide bond
            ndisulfide = ndisulfide + 1
            myaminonum(ndisulfide) = ipeptide
            myatomnum(ndisulfide) = hmain+1
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,0.0975d0)
            atomname(hmain+1) = 'S'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.2175d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
c           --- no side chain improper torsions         
         elseif ( pepname(ipeptide) .eq. 'q0' .or.
     &           pepname(ipeptide) .eq. 'Q0' ) then
            call twh_aa_pepname_q0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.12d0)
            atomname(hmain+1) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.12d0)
            atomname(hmain+2) = 'C n'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.50d0)
            atomname(hmain+3) = 'O'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.50d0)
            atomname(hmain+4) = 'N'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.76d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.06d0)
            atomname(shift+3) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.06d0)
            atomname(shift+4) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.38d0)
            atomname(shift+5) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.38d0)
c           --- side chain improper torsions
c           --- OPLS type O-Cn-X-Y = type 2.3 for CD
            call twh_inimprop(GLB_SET,imolty,hmain+2,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,hmain+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,hmain+4)
c           --- OPLS type Z-N-X-Y = type 3.1 for NE2
            call twh_inimprop(GLB_SET,imolty,hmain+4,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,shift+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,shift+5)
         elseif ( pepname(ipeptide) .eq. 'g0' .or.
     &           pepname(ipeptide) .eq. 'G0' ) then
            call twh_aa_pepname_g0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
c           --- no improper torsions on side chain
         elseif ( pepname(ipeptide) .eq. 'i0' .or.
     &           pepname(ipeptide) .eq. 'I0' ) then
            call twh_aa_pepname_i0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.06d0)
            atomname(hmain+1) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.12d0)
            atomname(hmain+2) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.18d0)
            atomname(hmain+3) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.18d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.06d0)
            atomname(shift+3) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.06d0)
            atomname(shift+4) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.06d0)
            atomname(shift+5) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.06d0)
            atomname(shift+6) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.06d0)
            atomname(shift+7) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+7,0.06d0)
            atomname(shift+8) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+8,0.06d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 'l0' .or.
     &           pepname(ipeptide) .eq. 'L0' ) then
            call twh_aa_pepname_l0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.12d0)
            atomname(hmain+1) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.06d0)
            atomname(hmain+2) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.18d0)
            atomname(hmain+3) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.18d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.06d0)
            atomname(shift+3) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.06d0)
            atomname(shift+4) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.06d0)
            atomname(shift+5) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.06d0)
            atomname(shift+6) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.06d0)
            atomname(shift+7) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+7,0.06d0)
            atomname(shift+8) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+8,0.06d0)
c           --- no side chain impropers
c           --- Note these are the OPLS-aa 1996 parameters, not 2001
c           --- AAC 9/05
         elseif ( pepname(ipeptide) .eq. 'm0' .or.
     &           pepname(ipeptide) .eq. 'M0' ) then
            call twh_aa_pepname_m0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.12d0)
            atomname(hmain+1) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.0975d0)
            atomname(hmain+2) = 'S'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.435d0)
            atomname(hmain+3) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.0375d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.06d0)
            atomname(shift+3) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.06d0)
            atomname(shift+4) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.06d0)
            atomname(shift+5) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.06d0)
            atomname(shift+6) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.06d0)
c           --- no impropers on side chain
         elseif ( pepname(ipeptide) .eq. 'f0' .or.
     &           pepname(ipeptide) .eq. 'F0' ) then
            call twh_aa_pepname_f0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.005d0)
            atomname(hmain+1) = 'CA~'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.115d0)
            atomname(hmain+2) = 'CA~'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.115d0)
            atomname(hmain+3) = 'CA~'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.115d0)
            atomname(hmain+4) = 'CA~'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.115d0)
            atomname(hmain+5) = 'CA~'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.115d0)
            atomname(hmain+6) = 'CA~'
            call twh_qqatom(GLB_SET,imolty,hmain+6,-0.115d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.115d0)
            atomname(shift+3) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.115d0)
            atomname(shift+4) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.115d0)
            atomname(shift+5) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.115d0)
            atomname(shift+6) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.115d0)
c           --- OPLS type Z-Ca~-X-Y = type 1.1 for CD2
            call twh_inimprop(GLB_SET,imolty,hmain+3,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,1,shift+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,1,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,1,hmain+1)
c           --- OPLS type Z-Ca~-X-Y = type 1.1 for CE2
            call twh_inimprop(GLB_SET,imolty,hmain+5,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,1,shift+5)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,1,hmain+6)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,1,hmain+3)
c           --- OPLS type Z-Ca~-X-Y = type 1.1 for CE1
            call twh_inimprop(GLB_SET,imolty,hmain+4,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,shift+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+6)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+2)
c           --- OPLS type Z-Ca~-X-Y = type 1.1 for CD1
            call twh_inimprop(GLB_SET,imolty,hmain+2,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,shift+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,hmain+1)
c           --- OPLS type Z-Ca~-X-Y = type 1.1 for CG
            call twh_inimprop(GLB_SET,imolty,hmain+1,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+2)
c           --- OPLS type Z-Ca~-X-Y = type 1.1 for CZ
            call twh_inimprop(GLB_SET,imolty,hmain+6,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+6,1,shift+6)
            call twh_ijimprop3(GLB_SET,imolty,hmain+6,1,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+6,1,hmain+4)
         elseif ( pepname(ipeptide) .eq. 's0' .or.
     &           pepname(ipeptide) .eq. 'S0' ) then
            call twh_aa_pepname_s0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,0.145d0)
            atomname(hmain+1) = 'OHm'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.683d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'HO'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.418d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 't0' .or.
     &           pepname(ipeptide) .eq. 'T0' ) then
            call twh_aa_pepname_t0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,0.205d0)
            atomname(hmain+1) = 'OHm'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.683d0)
            atomname(hmain+2) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.18d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HO'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.418d0)
            atomname(shift+2) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.06d0)
            atomname(shift+3) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.06d0)
            atomname(shift+4) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.06d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 'w0' .or.
     &           pepname(ipeptide) .eq. 'W0' ) then
            call twh_aa_pepname_w0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.12d0)
            atomname(hmain+1) = 'C*'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.075d0)
            atomname(hmain+2) = 'CW'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.115d0)
            atomname(hmain+3) = 'CB'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.055d0)
            atomname(hmain+4) = 'NA'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.57d0)
            atomname(hmain+5) = 'CN'
            call twh_qqatom(GLB_SET,imolty,hmain+5,0.13d0)
            atomname(hmain+6) = 'CA~'
            call twh_qqatom(GLB_SET,imolty,hmain+6,-0.115d0)
            atomname(hmain+7) = 'CA~'
            call twh_qqatom(GLB_SET,imolty,hmain+7,-0.115d0)
            atomname(hmain+8) = 'CA~'
            call twh_qqatom(GLB_SET,imolty,hmain+8,-0.115d0)
            atomname(hmain+9) = 'CA~'
            call twh_qqatom(GLB_SET,imolty,hmain+9,-0.115d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.115d0)
            atomname(shift+3) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.42d0)
            atomname(shift+4) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.115d0)
            atomname(shift+5) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.115d0)
            atomname(shift+6) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.115d0)
            atomname(shift+7) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+7,0.115d0)
c           --- OPLS type Z-N-X-Y = type 3.1 for NE1
            call twh_inimprop(GLB_SET,imolty,hmain+4,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,shift+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+2)
c           --- OPLS type Z-Ca~-X-Y = type 1.1 for CZ2
            call twh_inimprop(GLB_SET,imolty,hmain+7,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+7,1,shift+5)
            call twh_ijimprop3(GLB_SET,imolty,hmain+7,1,hmain+9)
            call twh_ijimprop4(GLB_SET,imolty,hmain+7,1,hmain+5)
c           --- OPLS type Z-Ca~-X-Y = type 1.1 for CH2
            call twh_inimprop(GLB_SET,imolty,hmain+9,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+9,1,shift+7)
            call twh_ijimprop3(GLB_SET,imolty,hmain+9,1,hmain+8)
            call twh_ijimprop4(GLB_SET,imolty,hmain+9,1,hmain+7)
c           --- OPLS type Z-Ca~-X-Y = type 1.1 for CZ3
            call twh_inimprop(GLB_SET,imolty,hmain+8,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+8,1,shift+6)
            call twh_ijimprop3(GLB_SET,imolty,hmain+8,1,hmain+6)
            call twh_ijimprop4(GLB_SET,imolty,hmain+8,1,hmain+9)
c           --- OPLS type Z-Ca~-X-Y = type 1.1 for CE3
            call twh_inimprop(GLB_SET,imolty,hmain+6,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+6,1,shift+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+6,1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+6,1,hmain+8)
c           --- OPLS type Z-Ca~-X-Y = type 1.1 for CD1
            call twh_inimprop(GLB_SET,imolty,hmain+2,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,shift+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,hmain+1)
c           --- OPLS type Z-C*-X-Y = type 1.5 for CG
            call twh_inimprop(GLB_SET,imolty,hmain+1,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+2)
         elseif ( pepname(ipeptide) .eq. 'y0' .or.
     &           pepname(ipeptide) .eq. 'Y0' ) then
            call twh_aa_pepname_y0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.005d0)
            atomname(hmain+1) = 'CA~'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.115d0)
            atomname(hmain+2) = 'CA~'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.115d0)
            atomname(hmain+3) = 'CA~'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.115d0)
            atomname(hmain+4) = 'CA~'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.115d0)
            atomname(hmain+5) = 'CA~'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.115d0)
            atomname(hmain+6) = 'CA~'
            call twh_qqatom(GLB_SET,imolty,hmain+6,0.150d0)
            atomname(hmain+7) = 'OHp'
            call twh_qqatom(GLB_SET,imolty,hmain+7,-0.585d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.115d0)
            atomname(shift+3) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.115d0)
            atomname(shift+4) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.115d0)
            atomname(shift+5) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.115d0)
            atomname(shift+6) = 'HO'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.435d0)
c           --- OPLS type Z-Ca~-X-Y = type 1.1 for CD2
            call twh_inimprop(GLB_SET,imolty,hmain+3,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,1,shift+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,1,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,1,hmain+1)
c           --- OPLS type Z-Ca~-X-Y = type 1.1 for CE2
            call twh_inimprop(GLB_SET,imolty,hmain+5,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,1,shift+5)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,1,hmain+6)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,1,hmain+3)
c           --- OPLS type Z-Ca~-X-Y = type 1.1 for CE1
            call twh_inimprop(GLB_SET,imolty,hmain+4,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,shift+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+6)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+2)
c           --- OPLS type Z-Ca~-X-Y = type 1.1 for CD1
            call twh_inimprop(GLB_SET,imolty,hmain+2,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,shift+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,hmain+1)
c           --- OPLS type Z-Ca~-X-Y = type 1.1 for CG
            call twh_inimprop(GLB_SET,imolty,hmain+1,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+2)
c           --- OPLS type Z-Ca~-X-Y = type 1.1 for CZ
            call twh_inimprop(GLB_SET,imolty,hmain+6,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+6,1,hmain+7)
            call twh_ijimprop3(GLB_SET,imolty,hmain+6,1,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+6,1,hmain+4)
         elseif ( pepname(ipeptide) .eq. 'v0' .or.
     &           pepname(ipeptide) .eq. 'V0' ) then
            call twh_aa_pepname_v0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.06d0)
            atomname(hmain+1) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.18d0)
            atomname(hmain+2) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.18d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.06d0)
            atomname(shift+3) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.06d0)
            atomname(shift+4) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.06d0)
            atomname(shift+5) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.06d0)
            atomname(shift+6) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.06d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 'd-' .or.
     &           pepname(ipeptide) .eq. 'D-' ) then
            call twh_aa_pepname_d_neg(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.22d0)
            atomname(hmain+1) = 'C a'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.7d0)
            atomname(hmain+2) = 'O2'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.8d0)
            atomname(hmain+3) = 'O2'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.8d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
c           --- side chain improper torsion
c           --- OPLS type O-Ca-X-Y = type 2.1 for CG
            call twh_inimprop(GLB_SET,imolty,hmain+1,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+2)
         elseif ( pepname(ipeptide) .eq. 'e-' .or.
     &           pepname(ipeptide) .eq. 'E-' ) then
            call twh_aa_pepname_e_neg(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.12d0)
            atomname(hmain+1) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.22d0)
            atomname(hmain+2) = 'C a'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.7d0)
            atomname(hmain+3) = 'O2'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.8d0)
            atomname(hmain+4) = 'O2'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.8d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.06d0)
            atomname(shift+3) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.06d0)
c           --- side chain improper torsion
c           --- OPLS type O-Ca-X-Y = type 2.1 for CD
            call twh_inimprop(GLB_SET,imolty,hmain+2,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,hmain+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,hmain+3)
         elseif ( pepname(ipeptide) .eq. 'r+' .or.
     &           pepname(ipeptide) .eq. 'R+' ) then
            call twh_aa_pepname_r_plus(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.12d0)
            atomname(hmain+1) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.05d0)
            atomname(hmain+2) = 'CT '
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.19d0)
            atomname(hmain+3) = 'N2'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.70d0)
            atomname(hmain+4) = 'CA+'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.64d0)
            atomname(hmain+5) = 'N2'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.80d0)
            atomname(hmain+6) = 'N2'
            call twh_qqatom(GLB_SET,imolty,hmain+6,-0.80d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.06d0)
            atomname(shift+3) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.06d0)
            atomname(shift+4) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.06d0)
            atomname(shift+5) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.06d0)
            atomname(shift+6) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.44d0)
            atomname(shift+7) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+7,0.46d0)
            atomname(shift+8) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+8,0.46d0)
            atomname(shift+9) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+9,0.46d0)
            atomname(shift+10) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+10,0.46d0)
c           --- side chain improper torsions
c           --- OPLS type Z-CA+-X-Y = type 1.4 for CZ
            call twh_inimprop(GLB_SET,imolty,hmain+4,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,hmain+6)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+3)
c           --- OPLS type Z-N-X-Y = type 3.1 for NE
            call twh_inimprop(GLB_SET,imolty,hmain+3,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,1,shift+6)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,1,hmain+2)
c           --- OPLS type Z-N-X-Y = type 3.1 for NH1
            call twh_inimprop(GLB_SET,imolty,hmain+5,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,1,shift+8)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,1,shift+7)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,1,hmain+4)
c           --- OPLS type Z-N-X-Y = type 3.1 for NH2
            call twh_inimprop(GLB_SET,imolty,hmain+6,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+6,1,shift+10)
            call twh_ijimprop3(GLB_SET,imolty,hmain+6,1,shift+9)
            call twh_ijimprop4(GLB_SET,imolty,hmain+6,1,hmain+4)
         elseif ( pepname(ipeptide) .eq. 'h+' .or.
     &           pepname(ipeptide) .eq. 'H+' ) then
            call twh_aa_pepname_h_plus(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.05d0)
            atomname(hmain+1) = 'CW'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.215d0)
            atomname(hmain+2) = 'NA'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.54d0)
            atomname(hmain+3) = 'CW'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.215d0)
            atomname(hmain+4) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.385d0)
            atomname(hmain+5) = 'NA'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.54d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.46d0)
            atomname(shift+3) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.115d0)
            atomname(shift+4) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.115d0)
            atomname(shift+5) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.46d0)
c           --- side chain improper torsions
c           --- OPLS type Z-N-X-Y = type 3.1 for ND1
            call twh_inimprop(GLB_SET,imolty,hmain+2,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,shift+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,hmain+1)
c           --- OPLS type Z-N-X-Y = type 3.1 for NE2
            call twh_inimprop(GLB_SET,imolty,hmain+5,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,1,shift+5)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,1,hmain+3)
c           --- OPLS type Z-CW-X-Y = type 1.2 for CD2
            call twh_inimprop(GLB_SET,imolty,hmain+3,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,1,shift+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,1,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,1,hmain+1)
c           --- OPLS type Z-CR-X-Y = type 1.3 for CE1
            call twh_inimprop(GLB_SET,imolty,hmain+4,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,shift+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+2)
c           --- OPLS type Z-CW-X-Y = type 1.2 for CG
            call twh_inimprop(GLB_SET,imolty,hmain+1,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+2)
         elseif ( pepname(ipeptide) .eq. 'hd' .or.
     &           pepname(ipeptide) .eq. 'HD' .or.
     &           pepname(ipeptide) .eq. 'Hd' ) then
            call twh_aa_pepname_hd(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.005d0)
            atomname(hmain+1) = 'CW'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.015d0)
            atomname(hmain+2) = 'NA'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.57d0)
            atomname(hmain+3) = 'CV'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.015d0)
            atomname(hmain+4) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.295d0)
            atomname(hmain+5) = 'NB'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.49d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.42d0)
            atomname(shift+3) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.115d0)
            atomname(shift+4) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.115d0)
c           --- side chain improper torsions
c           --- OPLS type Z-CW-X-Y = type 112 for CG
            call twh_inimprop(GLB_SET,imolty,hmain+1,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+2)
c           --- OPLS type Z-N-X-Y = type 3.1 for ND1
            call twh_inimprop(GLB_SET,imolty,hmain+2,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,shift+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,hmain+1)
c     --- OPLS type Z-CV-X-Y type 1.4 for CD2
            call twh_inimprop(GLB_SET,imolty,hmain+3,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,1,shift+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,1,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,1,hmain+1)
c     --- OPLS type Z-CR-X-Y type 1.3 for CE1
            call twh_inimprop(GLB_SET,imolty,hmain+4,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,shift+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+5)
         elseif ( pepname(ipeptide) .eq. 'he' .or.
     &           pepname(ipeptide) .eq. 'HE' .or.
     &           pepname(ipeptide) .eq. 'He' ) then
            call twh_aa_pepname_he(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.005d0)
            atomname(hmain+1) = 'CV'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.015d0)
            atomname(hmain+2) = 'NB'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.49d0)
            atomname(hmain+3) = 'CW'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.015d0)
            atomname(hmain+4) = 'CR'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.295d0)
            atomname(hmain+5) = 'NA'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.57d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.115d0)
            atomname(shift+3) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.115d0)
            atomname(shift+4) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.42d0)
c           --- side chain improper torsions
c           --- OPLS type Z-CV-X-Y = type 1.4 for CG
            call twh_inimprop(GLB_SET,imolty,hmain+1,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+3)
c           --- OPLS type Z-CW-X-Y = type 1.2 for CD2
            call twh_inimprop(GLB_SET,imolty,hmain+3,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,1,shift+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,1,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,1,hmain+1)
c           --- OPLS type Z-CR-X-Y = type 1.3 for CE1
            call twh_inimprop(GLB_SET,imolty,hmain+4,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,shift+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+5)
c     --- OPLS type Z-N-X-Y = type 3.1 for NE2
            call twh_inimprop(GLB_SET,imolty,hmain+5,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,1,shift+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,1,hmain+3)
         elseif ( pepname(ipeptide) .eq. 'k+' .or.
     &           pepname(ipeptide) .eq. 'K+' ) then
            call twh_aa_pepname_k_plus(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.120d0)
            atomname(hmain+1) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.120d0)
            atomname(hmain+2) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.120d0)
            atomname(hmain+3) = 'CT'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.190d0)
            atomname(hmain+4) = 'NT+2'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.30d0)
            atomname(shift) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift,0.06d0)
            atomname(shift+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.06d0)
            atomname(shift+2) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.06d0)
            atomname(shift+3) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.06d0)
            atomname(shift+4) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.06d0)
            atomname(shift+5) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.06d0)
            atomname(shift+6) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.06d0)
            atomname(shift+7) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,shift+7,0.06d0)
            atomname(shift+8) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+8,0.33d0)
            atomname(shift+9) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+9,0.33d0)
            atomname(shift+10) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+10,0.33d0)
c           --- no side chain improper torsions
         else
            write(6,*) 'invalid peptide name ',ipeptide
     &           ,' ',pepname(ipeptide)
            lfailure = .true.
            return
         endif

         call twh_aa_pepname_main_chain(lfailure,aminoname,termflag
     &        ,imolty,total,iatom,icb,iha,ibondend,heavy,light,ninsert
     &        ,shift,hmain)
         if ( lfailure ) return
         if ( termflag .eq. -1 ) then
c           --- n-terminus main chain
            atomname(iatom) =   'NT+2'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.30d0)
            atomname(iatom+1) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.25d0)
            atomname(iatom+2) = 'C n'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.50d0)
            atomname(iatom+3) = 'O'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.50d0)
c           --- N hydrogen
            atomname(hmain+heavy) = 'H'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.33d0)
c           --- extra N-terminus hydrogens
            atomname(hmain+heavy+1) = 'H'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.33d0)
            atomname(hmain+heavy+2) = 'H'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.33d0)
c           --- C alpha hydrogen
            atomname(hmain+heavy+3) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.06d0)
c           --- improper torsions
c           --- backbone C-Ca-N-O
c           --- OPLS type Z-Cn-X-Y = type 2.3
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+3)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,iatom+total)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,iatom+1)
         elseif ( termflag .eq. 0 ) then
c           --- interior segment main chain
            atomname(iatom) =   'N'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.50d0)
            atomname(iatom+1) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.14d0)
            atomname(iatom+2) = 'C n'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.50d0)
            atomname(iatom+3) = 'O'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.50d0)
            atomname(hmain+heavy) = 'H'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.30d0)
            atomname(hmain+heavy+1) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.06d0)
c           --- improper torsions
c           --- backbone N-C-Ca-Hn 
c           --- OPLS Z-N-X-Y is type 3.1
            call twh_inimprop(GLB_SET,imolty,iatom,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom,1,hmain+heavy)
            call twh_ijimprop3(GLB_SET,imolty,iatom,1,iatom+1)
            call twh_ijimprop4(GLB_SET,imolty,iatom,1,ibondend)
c           --- backbone C-Ca-N-O
c           --- OPLS O-Cn-X-Y is type 2,3
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+3)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,iatom+total)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,iatom+1)
         elseif ( termflag .eq. 1 ) then
c           --- c-terminus main chain
            atomname(iatom) =   'N'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.50d0)
            atomname(iatom+1) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.04d0)
            atomname(iatom+2) = 'C a'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.70d0)
            atomname(iatom+3) = 'O2'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.80d0)
c           --- c-terminus oxygen
            atomname(hmain+heavy) = 'O2'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.80d0)
c           --- main chain hydrogens
            atomname(hmain+heavy+1) = 'H'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.30d0)
            atomname(hmain+heavy+2) = 'HC~1'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.06d0)
c           --- improper torsions
c           --- backbone N-C-Ca-Hn
c           --- OPLS type Z-N-X-Y = type 3.1
            call twh_inimprop(GLB_SET,imolty,iatom,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom,1,hmain+heavy+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom,1,iatom+1)
            call twh_ijimprop4(GLB_SET,imolty,iatom,1,ibondend)
c           --- C-terminus C-Ca-O-O
c           --- OPLS type O-Cn-X-Y = type 2.3
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+3)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,hmain+heavy)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,iatom+1)
         else
c           --- invalid amino acid designation
            write(6,*) 'imolty,ipeptide,pepname',imolty,ipeptide
     &           ,pepname(ipeptide),' is not valid'
            lfailure = .true.
            return
         endif
c        --- update ibondend for the typical amino acids
         ibondend = iatom+2
         if ( aminoname .eq. 'glycine' ) then
c           --- if this was a glycine residue need to reduce
c           --- charge on the main chain C alpha
            atomname(iatom+1) = 'CT'
            call twh_qqatom(GLB_GET,imolty,iatom+1,charge)
            charge = charge - 0.06d0
            call twh_qqatom(GLB_SET,imolty,iatom+1,charge)
         endif
      endif
      
c     write(6,*) 'finished POLYOPLSAA'
      
      return
      end

#else
c     Compilation disabled (for speed) with FFLAGS=-DDISABLE_POLY
c     statement in ./configure
      write(6,*) "polyoplsaa disabled (compile flag)"
      lfailure = .true.
      return
      end
#endif
