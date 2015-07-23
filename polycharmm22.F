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
      subroutine twh_polycharmm22(lfailure,imolty,iatom,termflag,total
     &,ibondend,ipeptide,ndisulfide,myaminonum,myatomnum,iha,icb
     &,pepname,atomname,aminoname,shortname)
c     ******************************************************************
c     * builds the input file for an amino acid sequence starting from *
c     * the N-terminus for the Charmm22 force field                    *
c     *                                                                *
c     * split out of buildprot 12-12-2002 by M.G. Martin               *
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

c      write(6,*) 'start POLYCHARMM22'

c     --- determine which amino acid is built this step
c     --- check the special cases first, otherwise this is a generic
c     --- amino acid on a generic backbone.
      if ( pepname(ipeptide) .eq. 'ze' ) then
c        --- References for this group:
c        --- Ketchem, R. R.; Roux, B.; Cross, T. A.
c        `High-resolution polypeptide structure in a lamerllar
c        phase lipid environment from solid state NMR derived
c        orientational constraints'
c        Structure, 5, 1655-1669 (1997).
c        --- T.B. Woolf personal communication 08-16-2001
c        (M.G. Martin)
         call twh_aa_pepname_ze(lfailure,aminoname,shortname,termflag
     &        ,imolty,total,iatom,ibondend)
         if ( lfailure ) return
         atomname(iatom) = 'NH1'
         call twh_qqatom(GLB_SET,imolty,iatom,-0.47d0)
         atomname(iatom+1) = 'CT2'
         call twh_qqatom(GLB_SET,imolty,iatom+1,-0.02d0)
         atomname(iatom+2) = 'CT2'
         call twh_qqatom(GLB_SET,imolty,iatom+2,0.05d0)
         atomname(iatom+3) = 'OH1'
         call twh_qqatom(GLB_SET,imolty,iatom+3,-0.66d0)
         atomname(iatom+4) = 'H'
         call twh_qqatom(GLB_SET,imolty,iatom+4,0.31d0)
         atomname(iatom+5) = 'HA'
         call twh_qqatom(GLB_SET,imolty,iatom+5,0.09d0)
         atomname(iatom+6) = 'HA'
         call twh_qqatom(GLB_SET,imolty,iatom+6,0.09d0)
         atomname(iatom+7) = 'HA'
         call twh_qqatom(GLB_SET,imolty,iatom+7,0.09d0)
         atomname(iatom+8) = 'HA'
         call twh_qqatom(GLB_SET,imolty,iatom+8,0.09d0)
         atomname(iatom+9) = 'H'
         call twh_qqatom(GLB_SET,imolty,iatom+9,0.43d0)
c        --- improper torsion on "backbone" N-C-Ca-Hn
         call twh_inimprop(GLB_SET,imolty,iatom,1)
         call twh_ijimprop2(GLB_SET,imolty,iatom,1,ibondend)
         call twh_ijimprop3(GLB_SET,imolty,iatom,1,iatom+1)
         call twh_ijimprop4(GLB_SET,imolty,iatom,1,iatom+4)
c        --- ibondend not relevent for c-terminus
         ibondend = -1
      elseif ( pepname(ipeptide) .eq. 'zm' ) then
         call twh_aa_pepname_zm(lfailure,aminoname,shortname,termflag
     &        ,imolty,total,iatom,ibondend)
         if ( lfailure ) return
c        --- References for this group:
c        --- taken from the Charmm22 top_all22_prot_lipid.inp
c        --- as found on the MacKerrel web site
c        --- EXCEPT charges on the H modified to be identical as 
c        --- no method of maintaining the cis and trans is suggested
c        --- by Charmm22 (they use 0.32 and 0.30 instead of 2 0.31)
         atomname(iatom) = 'NH2'
         call twh_qqatom(GLB_SET,imolty,iatom,-0.62d0)
         atomname(iatom+1) = 'H'
         call twh_qqatom(GLB_SET,imolty,iatom+1,0.31d0)
         atomname(iatom+2) = 'H'
         call twh_qqatom(GLB_SET,imolty,iatom+2,0.31d0)
c        --- also need to modify the charges and atom types 
c        --- on the previous amino acid backbone
c        --- previous main chain C
         atomname(ibondend) = 'CC'
         call twh_qqatom(GLB_SET,imolty,ibondend,0.55d0)
c        --- previous main chain O
         call twh_qqatom(GLB_SET,imolty,ibondend+1,-0.55d0)
c        --- improper torsion on NT-C-HT1-HT2
         call twh_inimprop(GLB_SET,imolty,iatom,1)
         call twh_ijimprop2(GLB_SET,imolty,iatom,1,ibondend)
         call twh_ijimprop3(GLB_SET,imolty,iatom,1,iatom+1)
         call twh_ijimprop4(GLB_SET,imolty,iatom,1,iatom+2)
c        --- improper torsion on NT-C-HT2-HT1
         call twh_inimprop(GLB_SET,imolty,iatom,2)
         call twh_ijimprop2(GLB_SET,imolty,iatom,2,ibondend)
         call twh_ijimprop3(GLB_SET,imolty,iatom,2,iatom+2)
         call twh_ijimprop4(GLB_SET,imolty,iatom,2,iatom+1)
c        --- ibondend not relevent for c-terminus
         ibondend = -1
      elseif ( pepname(ipeptide) .eq. 'zn' ) then
         call twh_aa_pepname_zn(lfailure,aminoname,shortname,termflag
     &        ,imolty,total,iatom,ibondend)
         if ( lfailure ) return
c        --- References for this group:
c        --- taken from the Charmm22 top_all22_prot_lipid.inp
c        --- as found on the MacKerrel web site
         atomname(iatom) = 'NH1'
         call twh_qqatom(GLB_SET,imolty,iatom,-0.47d0)
         atomname(iatom+1) = 'CT3'
         call twh_qqatom(GLB_SET,imolty,iatom+1,-0.11d0)
         atomname(iatom+2) = 'H'
         call twh_qqatom(GLB_SET,imolty,iatom+2,0.31d0)
         atomname(iatom+3) = 'HA'
         call twh_qqatom(GLB_SET,imolty,iatom+3,0.09d0)
         atomname(iatom+4) = 'HA'
         call twh_qqatom(GLB_SET,imolty,iatom+4,0.09d0)
         atomname(iatom+5) = 'HA'
         call twh_qqatom(GLB_SET,imolty,iatom+5,0.09d0)
c        --- improper torsion on "backbone" NT-C-CAT-HNT
         call twh_inimprop(GLB_SET,imolty,iatom,1)
         call twh_ijimprop2(GLB_SET,imolty,iatom,1,ibondend)
         call twh_ijimprop3(GLB_SET,imolty,iatom,1,iatom+1)
         call twh_ijimprop4(GLB_SET,imolty,iatom,1,iatom+2)
c        --- ibondend not relevent for c-terminus
         ibondend = -1
      elseif ( pepname(ipeptide) .eq. 'zf' ) then
         call twh_aa_pepname_zf(lfailure,aminoname,shortname,termflag
     &        ,imolty,total,iatom)
         if ( lfailure ) return
c        --- References for this group:
c        --- Ketchem, R. R.; Roux, B.; Cross, T. A.
c        `High-resolution polypeptide structure in a lamerllar
c        phase lipid environment from solid state NMR derived
c        orientational constraints'
c        Structure, 5, 1655-1669 (1997).
c        --- T.B. Woolf personal communication 08-16-2001
c        (M.G. Martin)
         atomname(iatom) = 'C'
         call twh_qqatom(GLB_SET,imolty,iatom,0.42d0)
         atomname(iatom+1) = 'O'
         call twh_qqatom(GLB_SET,imolty,iatom+1,-0.51d0)
         atomname(iatom+2) = 'HA'
         call twh_qqatom(GLB_SET,imolty,iatom+2,0.09d0)
c        --- improper torsion
         call twh_inimprop(GLB_SET,imolty,iatom,1)
         call twh_ijimprop2(GLB_SET,imolty,iatom,1,iatom+2)
         call twh_ijimprop3(GLB_SET,imolty,iatom,1,iatom+3)
         call twh_ijimprop4(GLB_SET,imolty,iatom,1,iatom+1)
c        --- update ibondend
         ibondend = iatom
      elseif ( pepname(ipeptide) .eq. 'za' ) then
c        --- acetylated N-terminus
         call twh_aa_pepname_za(lfailure,aminoname,shortname,termflag
     &        ,imolty,total,iatom)
         if ( lfailure ) return
c        --- References for this group:
c        --- taken from the top_all22_prot_lipid.inp file
c        --- from the MacKerrel web site
         atomname(iatom) = 'C'
         call twh_qqatom(GLB_SET,imolty,iatom,0.51d0)
         atomname(iatom+1) = 'CT3'
         call twh_qqatom(GLB_SET,imolty,iatom+1,-0.27d0)
         atomname(iatom+2) = 'O'
         call twh_qqatom(GLB_SET,imolty,iatom+2,-0.51d0)
         atomname(iatom+3) = 'HA'
         call twh_qqatom(GLB_SET,imolty,iatom+3,0.09d0)
         atomname(iatom+4) = 'HA'
         call twh_qqatom(GLB_SET,imolty,iatom+4,0.09d0)
         atomname(iatom+5) = 'HA'
         call twh_qqatom(GLB_SET,imolty,iatom+5,0.09d0)
c        --- improper torsion CY-CAY-N(next)-OY
         call twh_inimprop(GLB_SET,imolty,iatom,1)
         call twh_ijimprop2(GLB_SET,imolty,iatom,1,iatom+1)
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
            atomname(iatom) = 'NP  '
            call twh_qqatom(GLB_SET,imolty,iatom,-0.07d0)
            atomname(iatom+1) = 'CP1 '
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.16d0)
            atomname(iatom+2) = 'C   '
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.51d0)
            atomname(iatom+3) = 'O   '
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.51d0)
            atomname(iatom+4) = 'CP2 '
            call twh_qqatom(GLB_SET,imolty,iatom+4,-0.18d0)
            atomname(iatom+5) = 'CP2 '
            call twh_qqatom(GLB_SET,imolty,iatom+5,-0.18d0)
            atomname(iatom+6) = 'CP3 '
            call twh_qqatom(GLB_SET,imolty,iatom+6,0.16d0)
            atomname(iatom+7) = 'HC'
            call twh_qqatom(GLB_SET,imolty,iatom+7,0.24d0)
            atomname(iatom+8) = 'HC'
            call twh_qqatom(GLB_SET,imolty,iatom+8,0.24d0)
            atomname(iatom+9) = 'HB'
            call twh_qqatom(GLB_SET,imolty,iatom+9,0.09d0)
            atomname(iatom+10) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+10,0.09d0)
            atomname(iatom+11) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+11,0.09d0)
            atomname(iatom+12) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+12,0.09d0)
            atomname(iatom+13) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+13,0.09d0)
            atomname(iatom+14) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+14,0.09d0)
            atomname(iatom+15) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+15,0.09d0)
c           --- improper torsions
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,iatom+total)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,iatom+3)
         elseif ( termflag .eq. 0 ) then
            atomname(iatom) = 'N   '
            call twh_qqatom(GLB_SET,imolty,iatom,-0.29d0)
            atomname(iatom+1) = 'CP1 '
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.02d0)
            atomname(iatom+2) = 'C   '
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.51d0)
            atomname(iatom+3) = 'O   '
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.51d0)
            atomname(iatom+4) = 'CP2 '
            call twh_qqatom(GLB_SET,imolty,iatom+4,-0.18d0)
            atomname(iatom+5) = 'CP2 '
            call twh_qqatom(GLB_SET,imolty,iatom+5,-0.18d0)
            atomname(iatom+6) = 'CP3 '
            call twh_qqatom(GLB_SET,imolty,iatom+6,0.00d0)
            atomname(iatom+7) = 'HB  '
            call twh_qqatom(GLB_SET,imolty,iatom+7,0.09d0)
            atomname(iatom+8) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+8,0.09d0)
            atomname(iatom+9) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+9,0.09d0)
            atomname(iatom+10) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+10,0.09d0)
            atomname(iatom+11) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+11,0.09d0)
            atomname(iatom+12) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+12,0.09d0)
            atomname(iatom+13) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+13,0.09d0)
c           --- improper torsions
            call twh_inimprop(GLB_SET,imolty,iatom,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom,1,ibondend)
            call twh_ijimprop3(GLB_SET,imolty,iatom,1,iatom+1)
            call twh_ijimprop4(GLB_SET,imolty,iatom,1,iatom+6)
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,iatom+total)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,iatom+3)
         elseif ( termflag .eq. 1 ) then
            atomname(iatom) = 'N   '
            call twh_qqatom(GLB_SET,imolty,iatom,-0.29d0)
            atomname(iatom+1) = 'CP1 '
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.02d0)
            atomname(iatom+2) = 'CC  '
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.34d0)
            atomname(iatom+3) = 'OC  '
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.67d0)
            atomname(iatom+4) = 'CP2 '
            call twh_qqatom(GLB_SET,imolty,iatom+4,-0.18d0)
            atomname(iatom+5) = 'CP2 '
            call twh_qqatom(GLB_SET,imolty,iatom+5,-0.18d0)
            atomname(iatom+6) = 'CP3 '
            call twh_qqatom(GLB_SET,imolty,iatom+6,0.00d0)
            atomname(iatom+7) = 'OC  '
            call twh_qqatom(GLB_SET,imolty,iatom+7,-0.67d0)
            atomname(iatom+8) = 'HB  '
            call twh_qqatom(GLB_SET,imolty,iatom+8,0.09d0)
            atomname(iatom+9) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+9,0.09d0)
            atomname(iatom+10) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+10,0.09d0)
            atomname(iatom+11) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+11,0.09d0)
            atomname(iatom+12) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+12,0.09d0)
            atomname(iatom+13) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+13,0.09d0)
            atomname(iatom+14) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,iatom+14,0.09d0)
c           --- improper torsions
            call twh_inimprop(GLB_SET,imolty,iatom,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom,1,ibondend)
            call twh_ijimprop3(GLB_SET,imolty,iatom,1,iatom+1)
            call twh_ijimprop4(GLB_SET,imolty,iatom,1,iatom+6)
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,iatom+7)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,iatom+3)
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
            aminoname = 'alanine'
            shortname = 'ALA'
            call twh_aa_pepname_a0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT3 '
            call twh_qqatom(GLB_SET,imolty,hmain,-0.27d0)
            atomname(shift) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,shift+2,0.09d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 'n0' .or.
     &           pepname(ipeptide) .eq. 'N0' ) then
            call twh_aa_pepname_n0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.18d0)
            atomname(hmain+1) = 'CC'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.55d0)
            atomname(hmain+2) = 'O'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.55d0)
            atomname(hmain+3) = 'NH2'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.62d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.31d0)
            atomname(shift+3) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.31d0)
c           --- side chain improper torsions
            call twh_inimprop(GLB_SET,imolty,hmain+1,2)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+2)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,2,hmain+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,2,hmain)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,2,hmain+2)
            call twh_inimprop(GLB_SET,imolty,hmain+3,2)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,1,hmain+1)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,1,shift+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,1,shift+3)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,2,hmain+1)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,2,shift+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,2,shift+2)
         elseif ( pepname(ipeptide) .eq. 'c0' .or.
     &           pepname(ipeptide) .eq. 'C0' ) then
            call twh_aa_pepname_c0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.11d0)
            atomname(hmain+1) = 'S'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.23d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'HS'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.16d0)
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
            atomname(hmain) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.10d0)
            atomname(hmain+1) = 'SM'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.08d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 'cp' .or.
     &           pepname(ipeptide) .eq. 'CP' .or.
     &           pepname(ipeptide) .eq. 'cP' .or.
     &           pepname(ipeptide) .eq. 'Cp' ) then
            call twh_aa_pepname_cp(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.11d0)
            atomname(hmain+1) = 'S'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.28d0)
            atomname(hmain+2) = 'CL'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.63d0)
            atomname(hmain+3) = 'CTL2'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.08d0)
            atomname(hmain+4) = 'CTL2'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.18d0)
            atomname(hmain+5) = 'CTL2'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.18d0)
            atomname(hmain+6) = 'CTL2'
            call twh_qqatom(GLB_SET,imolty,hmain+6,-0.18d0)
            atomname(hmain+7) = 'CTL2'
            call twh_qqatom(GLB_SET,imolty,hmain+7,-0.18d0)
            atomname(hmain+8) = 'CTL2'
            call twh_qqatom(GLB_SET,imolty,hmain+8,-0.18d0)
            atomname(hmain+9) = 'CTL2'
            call twh_qqatom(GLB_SET,imolty,hmain+9,-0.18d0)
            atomname(hmain+10) = 'CTL2'
            call twh_qqatom(GLB_SET,imolty,hmain+10,-0.18d0)
            atomname(hmain+11) = 'CTL2'
            call twh_qqatom(GLB_SET,imolty,hmain+11,-0.18d0)
            atomname(hmain+12) = 'CTL2'
            call twh_qqatom(GLB_SET,imolty,hmain+12,-0.18d0)
            atomname(hmain+13) = 'CTL2'
            call twh_qqatom(GLB_SET,imolty,hmain+13,-0.18d0)
            atomname(hmain+14) = 'CTL2'
            call twh_qqatom(GLB_SET,imolty,hmain+14,-0.18d0)
            atomname(hmain+15) = 'CTL2'
            call twh_qqatom(GLB_SET,imolty,hmain+15,-0.18d0)
            atomname(hmain+16) = 'CTL2'
            call twh_qqatom(GLB_SET,imolty,hmain+16,-0.18d0)
            atomname(hmain+17) = 'CTL3'
            call twh_qqatom(GLB_SET,imolty,hmain+17,-0.27d0)
            atomname(hmain+18) = 'OBL'
            call twh_qqatom(GLB_SET,imolty,hmain+18,-0.52d0)
c           --- hydrogens
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.09d0)
            atomname(shift+3) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.09d0)
            atomname(shift+4) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.09d0)
            atomname(shift+5) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.09d0)
            atomname(shift+6) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.09d0)
            atomname(shift+7) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+7,0.09d0)
            atomname(shift+8) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+8,0.09d0)
            atomname(shift+9) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+9,0.09d0)
            atomname(shift+10) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+10,0.09d0)
            atomname(shift+11) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+11,0.09d0)
            atomname(shift+12) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+12,0.09d0)
            atomname(shift+13) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+13,0.09d0)
            atomname(shift+14) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+14,0.09d0)
            atomname(shift+15) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+15,0.09d0)
            atomname(shift+16) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+16,0.09d0)
            atomname(shift+17) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+17,0.09d0)
            atomname(shift+18) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+18,0.09d0)
            atomname(shift+19) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+19,0.09d0)
            atomname(shift+20) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+20,0.09d0)
            atomname(shift+21) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+21,0.09d0)
            atomname(shift+22) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+22,0.09d0)
            atomname(shift+23) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+23,0.09d0)
            atomname(shift+24) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+24,0.09d0)
            atomname(shift+25) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+25,0.09d0)
            atomname(shift+26) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+26,0.09d0)
            atomname(shift+27) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+27,0.09d0)
            atomname(shift+28) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+28,0.09d0)
            atomname(shift+29) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+29,0.09d0)
            atomname(shift+30) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+30,0.09d0)
            atomname(shift+31) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+31,0.09d0)
            atomname(shift+32) = 'HAL'
            call twh_qqatom(GLB_SET,imolty,shift+32,0.09d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 'q0' .or.
     &           pepname(ipeptide) .eq. 'Q0' ) then
            call twh_aa_pepname_q0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.18d0)
            atomname(hmain+1) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.18d0)
            atomname(hmain+2) = 'CC'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.55d0)
            atomname(hmain+3) = 'O'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.55d0)
            atomname(hmain+4) = 'NH2'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.62d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.09d0)
            atomname(shift+3) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.09d0)
            atomname(shift+4) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.31d0)
            atomname(shift+5) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.31d0)
c           --- side chain improper torsions
            call twh_inimprop(GLB_SET,imolty,hmain+4,2)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,shift+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,shift+5)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,2,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,2,shift+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,2,shift+4)
         elseif ( pepname(ipeptide) .eq. 'g0' .or.
     &           pepname(ipeptide) .eq. 'G0' ) then
            call twh_aa_pepname_g0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(shift) = 'HB  '
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
c           --- no improper torsions on side chain
         elseif ( pepname(ipeptide) .eq. 'i0' .or.
     &           pepname(ipeptide) .eq. 'I0' ) then
            call twh_aa_pepname_i0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT1'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.09d0)
            atomname(hmain+1) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.18d0)
            atomname(hmain+2) = 'CT3'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.27d0)
            atomname(hmain+3) = 'CT3'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.27d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.09d0)
            atomname(shift+3) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.09d0)
            atomname(shift+4) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.09d0)
            atomname(shift+5) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.09d0)
            atomname(shift+6) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.09d0)
            atomname(shift+7) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+7,0.09d0)
            atomname(shift+8) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+8,0.09d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 'l0' .or.
     &           pepname(ipeptide) .eq. 'L0' ) then
            call twh_aa_pepname_l0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.18d0)
            atomname(hmain+1) = 'CT1'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.09d0)
            atomname(hmain+2) = 'CT3'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.27d0)
            atomname(hmain+3) = 'CT3'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.27d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.09d0)
            atomname(shift+3) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.09d0)
            atomname(shift+4) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.09d0)
            atomname(shift+5) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.09d0)
            atomname(shift+6) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.09d0)
            atomname(shift+7) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+7,0.09d0)
            atomname(shift+8) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+8,0.09d0)
c           --- no side chain impropers
         elseif ( pepname(ipeptide) .eq. 'm0' .or.
     &           pepname(ipeptide) .eq. 'M0' ) then
            call twh_aa_pepname_m0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.18d0)
            atomname(hmain+1) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.14d0)
            atomname(hmain+2) = 'S'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.09d0)
            atomname(hmain+3) = 'CT3'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.22d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.09d0)
            atomname(shift+3) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.09d0)
            atomname(shift+4) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.09d0)
            atomname(shift+5) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.09d0)
            atomname(shift+6) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.09d0)
c           --- no impropers on side chain
         elseif ( pepname(ipeptide) .eq. 'f0' .or.
     &           pepname(ipeptide) .eq. 'F0' ) then
            call twh_aa_pepname_f0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.18d0)
            atomname(hmain+1) = 'CA'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.0d0)
            atomname(hmain+2) = 'CA'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.115d0)
            atomname(hmain+3) = 'CA'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.115d0)
            atomname(hmain+4) = 'CA'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.115d0)
            atomname(hmain+5) = 'CA'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.115d0)
            atomname(hmain+6) = 'CA'
            call twh_qqatom(GLB_SET,imolty,hmain+6,-0.115d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'HP'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.115d0)
            atomname(shift+3) = 'HP'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.115d0)
            atomname(shift+4) = 'HP'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.115d0)
            atomname(shift+5) = 'HP'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.115d0)
            atomname(shift+6) = 'HP'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.115d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 's0' .or.
     &           pepname(ipeptide) .eq. 'S0' ) then
            call twh_aa_pepname_s0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain,0.05d0)
            atomname(hmain+1) = 'OH1'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.66d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.43d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 't0' .or.
     &           pepname(ipeptide) .eq. 'T0' ) then
            call twh_aa_pepname_t0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT1'
            call twh_qqatom(GLB_SET,imolty,hmain,0.14d0)
            atomname(hmain+1) = 'OH1'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.66d0)
            atomname(hmain+2) = 'CT3'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.27d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.43d0)
            atomname(shift+2) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.09d0)
            atomname(shift+3) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.09d0)
            atomname(shift+4) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.09d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 'w0' .or.
     &           pepname(ipeptide) .eq. 'W0' ) then
            call twh_aa_pepname_w0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.18d0)
            atomname(hmain+1) = 'CY'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.03d0)
            atomname(hmain+2) = 'CA'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.035d0)
            atomname(hmain+3) = 'CPT'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.02d0)
            atomname(hmain+4) = 'NY'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.61d0)
            atomname(hmain+5) = 'CPT'
            call twh_qqatom(GLB_SET,imolty,hmain+5,0.13d0)
            atomname(hmain+6) = 'CA'
            call twh_qqatom(GLB_SET,imolty,hmain+6,-0.115d0)
            atomname(hmain+7) = 'CA'
            call twh_qqatom(GLB_SET,imolty,hmain+7,-0.115d0)
            atomname(hmain+8) = 'CA'
            call twh_qqatom(GLB_SET,imolty,hmain+8,-0.115d0)
            atomname(hmain+9) = 'CA'
            call twh_qqatom(GLB_SET,imolty,hmain+9,-0.115d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'HP'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.115d0)
            atomname(shift+3) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.38d0)
            atomname(shift+4) = 'HP'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.115d0)
            atomname(shift+5) = 'HP'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.115d0)
            atomname(shift+6) = 'HP'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.115d0)
            atomname(shift+7) = 'HP'
            call twh_qqatom(GLB_SET,imolty,shift+7,0.115d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 'y0' .or.
     &           pepname(ipeptide) .eq. 'Y0' ) then
            call twh_aa_pepname_y0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.18d0)
            atomname(hmain+1) = 'CA'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.00d0)
            atomname(hmain+2) = 'CA'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.115d0)
            atomname(hmain+3) = 'CA'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.115d0)
            atomname(hmain+4) = 'CA'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.115d0)
            atomname(hmain+5) = 'CA'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.115d0)
            atomname(hmain+6) = 'CA'
            call twh_qqatom(GLB_SET,imolty,hmain+6,0.11d0)
            atomname(hmain+7) = 'OH1'
            call twh_qqatom(GLB_SET,imolty,hmain+7,-0.54d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'HP'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.115d0)
            atomname(shift+3) = 'HP'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.115d0)
            atomname(shift+4) = 'HP'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.115d0)
            atomname(shift+5) = 'HP'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.115d0)
            atomname(shift+6) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.43d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 'v0' .or.
     &           pepname(ipeptide) .eq. 'V0' ) then
            call twh_aa_pepname_v0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT1'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.09d0)
            atomname(hmain+1) = 'CT3'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.27d0)
            atomname(hmain+2) = 'CT3'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.27d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.09d0)
            atomname(shift+3) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.09d0)
            atomname(shift+4) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.09d0)
            atomname(shift+5) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.09d0)
            atomname(shift+6) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.09d0)
c           --- no side chain improper torsions
         elseif ( pepname(ipeptide) .eq. 'd-' .or.
     &           pepname(ipeptide) .eq. 'D-' ) then
            call twh_aa_pepname_d_neg(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.28d0)
            atomname(hmain+1) = 'CC'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.62d0)
            atomname(hmain+2) = 'OC'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.76d0)
            atomname(hmain+3) = 'OC'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.76d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
c           --- side chain improper torsion
            call twh_inimprop(GLB_SET,imolty,hmain+1,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+3)
         elseif ( pepname(ipeptide) .eq. 'e-' .or.
     &           pepname(ipeptide) .eq. 'E-' ) then
            call twh_aa_pepname_e_neg(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.18d0)
            atomname(hmain+1) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.28d0)
            atomname(hmain+2) = 'CC'
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.62d0)
            atomname(hmain+3) = 'OC'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.76d0)
            atomname(hmain+4) = 'OC'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.76d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.09d0)
            atomname(shift+3) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.09d0)
c           --- side chain improper torsion
            call twh_inimprop(GLB_SET,imolty,hmain+2,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,hmain+1)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,hmain+4)
         elseif ( pepname(ipeptide) .eq. 'r+' .or.
     &           pepname(ipeptide) .eq. 'R+' ) then
            call twh_aa_pepname_r_plus(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT2 '
            call twh_qqatom(GLB_SET,imolty,hmain,-0.18d0)
            atomname(hmain+1) = 'CT2 '
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.18d0)
            atomname(hmain+2) = 'CT2 '
            call twh_qqatom(GLB_SET,imolty,hmain+2,0.20d0)
            atomname(hmain+3) = 'NC2 '
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.70d0)
            atomname(hmain+4) = 'C   '
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.64d0)
            atomname(hmain+5) = 'NC2 '
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.80d0)
            atomname(hmain+6) = 'NC2 '
            call twh_qqatom(GLB_SET,imolty,hmain+6,-0.80d0)
            atomname(shift) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,shift+2,0.09d0)
            atomname(shift+3) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,shift+3,0.09d0)
            atomname(shift+4) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,shift+4,0.09d0)
            atomname(shift+5) = 'HA  '
            call twh_qqatom(GLB_SET,imolty,shift+5,0.09d0)
            atomname(shift+6) = 'HC  '
            call twh_qqatom(GLB_SET,imolty,shift+6,0.44d0)
            atomname(shift+7) = 'HC  '
            call twh_qqatom(GLB_SET,imolty,shift+7,0.46d0)
            atomname(shift+8) = 'HC  '
            call twh_qqatom(GLB_SET,imolty,shift+8,0.46d0)
            atomname(shift+9) = 'HC  '
            call twh_qqatom(GLB_SET,imolty,shift+9,0.46d0)
            atomname(shift+10) = 'HC  '
            call twh_qqatom(GLB_SET,imolty,shift+10,0.46d0)
c           --- side chain improper torsions
            call twh_inimprop(GLB_SET,imolty,hmain+4,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,hmain+5)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+6)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+3)
         elseif ( pepname(ipeptide) .eq. 'h+' .or.
     &           pepname(ipeptide) .eq. 'H+' ) then
            call twh_aa_pepname_h_plus(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.05d0)
            atomname(hmain+1) = 'CPH1'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.19d0)
            atomname(hmain+2) = 'NR3'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.51d0)
            atomname(hmain+3) = 'CPH1'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.19d0)
            atomname(hmain+4) = 'CPH2'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.32d0)
            atomname(hmain+5) = 'NR3'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.51d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.44d0)
            atomname(shift+3) = 'HR1'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.13d0)
            atomname(shift+4) = 'HR2'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.18d0)
            atomname(shift+5) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.44d0)
c           --- side chain improper torsions
            call twh_inimprop(GLB_SET,imolty,hmain+2,2)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,hmain+1)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,shift+2)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,2,hmain+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,2,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,2,shift+2)
            call twh_inimprop(GLB_SET,imolty,hmain+5,2)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,1,hmain+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,1,shift+5)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,2,hmain+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,2,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,2,shift+5)
         elseif ( pepname(ipeptide) .eq. 'hd' .or.
     &           pepname(ipeptide) .eq. 'HD' .or.
     &           pepname(ipeptide) .eq. 'Hd' ) then
            call twh_aa_pepname_hd(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.09d0)
            atomname(hmain+1) = 'CPH1'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.05d0)
            atomname(hmain+2) = 'NR1'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.36d0)
            atomname(hmain+3) = 'CPH1'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.22d0)
            atomname(hmain+4) = 'CPH2'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.25d0)
            atomname(hmain+5) = 'NR2'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.70d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.32d0)
            atomname(shift+3) = 'HR3'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.10d0)
            atomname(shift+4) = 'HR1'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.13d0)
c           --- side chain improper torsions
            call twh_inimprop(GLB_SET,imolty,hmain+2,2)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,hmain+1)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,shift+2)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,2,hmain+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,2,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,2,shift+2)
            call twh_inimprop(GLB_SET,imolty,hmain+3,2)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,1,hmain+1)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,1,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,1,shift+3)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,2,hmain+5)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,2,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,2,shift+3)
            call twh_inimprop(GLB_SET,imolty,hmain+4,2)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,hmain+5)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,shift+4)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,2,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,2,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,2,shift+4)
         elseif ( pepname(ipeptide) .eq. 'he' .or.
     &           pepname(ipeptide) .eq. 'HE' .or.
     &           pepname(ipeptide) .eq. 'He' ) then
            call twh_aa_pepname_he(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.08d0)
            atomname(hmain+1) = 'CPH1'
            call twh_qqatom(GLB_SET,imolty,hmain+1,0.22d0)
            atomname(hmain+2) = 'NR2'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.70d0)
            atomname(hmain+3) = 'CPH1'
            call twh_qqatom(GLB_SET,imolty,hmain+3,-0.05d0)
            atomname(hmain+4) = 'CPH2'
            call twh_qqatom(GLB_SET,imolty,hmain+4,0.25d0)
            atomname(hmain+5) = 'NR1'
            call twh_qqatom(GLB_SET,imolty,hmain+5,-0.36d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'HR3'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.09d0)
            atomname(shift+3) = 'HR1'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.13d0)
            atomname(shift+4) = 'H'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.32d0)
c           --- side chain improper torsions
            call twh_inimprop(GLB_SET,imolty,hmain+3,2)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,1,hmain+1)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,1,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,1,shift+2)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,2,hmain+5)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,2,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,2,shift+2)
            call twh_inimprop(GLB_SET,imolty,hmain+4,2)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,hmain+5)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,shift+3)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,2,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,2,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,2,shift+3)
            call twh_inimprop(GLB_SET,imolty,hmain+5,2)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,1,hmain+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,1,shift+4)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,2,hmain+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,2,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,2,shift+4)
         elseif ( pepname(ipeptide) .eq. 'k+' .or.
     &           pepname(ipeptide) .eq. 'K+' ) then
            call twh_aa_pepname_k_plus(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain,-0.18d0)
            atomname(hmain+1) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain+1,-0.18d0)
            atomname(hmain+2) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain+2,-0.18d0)
            atomname(hmain+3) = 'CT2'
            call twh_qqatom(GLB_SET,imolty,hmain+3,0.21d0)
            atomname(hmain+4) = 'NH3'
            call twh_qqatom(GLB_SET,imolty,hmain+4,-0.30d0)
            atomname(shift) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift,0.09d0)
            atomname(shift+1) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+1,0.09d0)
            atomname(shift+2) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+2,0.09d0)
            atomname(shift+3) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+3,0.09d0)
            atomname(shift+4) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+4,0.09d0)
            atomname(shift+5) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+5,0.09d0)
            atomname(shift+6) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+6,0.05d0)
            atomname(shift+7) = 'HA'
            call twh_qqatom(GLB_SET,imolty,shift+7,0.05d0)
            atomname(shift+8) = 'HC'
            call twh_qqatom(GLB_SET,imolty,shift+8,0.33d0)
            atomname(shift+9) = 'HC'
            call twh_qqatom(GLB_SET,imolty,shift+9,0.33d0)
            atomname(shift+10) = 'HC'
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
            atomname(iatom) =   'NH3'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.30d0)
            atomname(iatom+1) = 'CT1'
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.21d0)
            atomname(iatom+2) = 'C'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.51d0)
            atomname(iatom+3) = 'O'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.51d0)
c           --- N hydrogen
            atomname(hmain+heavy) = 'HC'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.33d0)
c           --- extra N-terminus hydrogens
            atomname(hmain+heavy+1) = 'HC'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.33d0)
            atomname(hmain+heavy+2) = 'HC'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.33d0)
c           --- C alpha hydrogen
            atomname(hmain+heavy+3) = 'HB'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.10d0)
c           --- improper torsions
c           --- backbone C-Ca-N-O
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,iatom+total)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,iatom+3)
         elseif ( termflag .eq. 0 ) then
c           --- interior segment main chain
            atomname(iatom) =   'NH1'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.47d0)
            atomname(iatom+1) = 'CT1'
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.07d0)
            atomname(iatom+2) = 'C'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.51d0)
            atomname(iatom+3) = 'O'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.51d0)
            atomname(hmain+heavy) = 'H'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.31d0)
            atomname(hmain+heavy+1) = 'HB'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.09d0)
c           --- improper torsions
c           --- backbone N-C-Ca-Hn
            call twh_inimprop(GLB_SET,imolty,iatom,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom,1,ibondend)
            call twh_ijimprop3(GLB_SET,imolty,iatom,1,iatom+1)
            call twh_ijimprop4(GLB_SET,imolty,iatom,1,hmain+heavy)
c           --- backbone C-Ca-N-O
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,iatom+total)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,iatom+3)
         elseif ( termflag .eq. 1 ) then
c           --- c-terminus main chain
            atomname(iatom) =   'NH1'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.47d0)
            atomname(iatom+1) = 'CT1'
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.07d0)
            atomname(iatom+2) = 'CC'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.34d0)
            atomname(iatom+3) = 'OC'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.67d0)
c           --- c-terminus oxygen
            atomname(hmain+heavy) = 'OC'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.67d0)
c           --- main chain hydrogens
            atomname(hmain+heavy+1) = 'H'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.31d0)
            atomname(hmain+heavy+2) = 'HB'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.09d0)
c           --- improper torsions
c           --- backbone N-C-Ca-Hn
            call twh_inimprop(GLB_SET,imolty,iatom,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom,1,ibondend)
            call twh_ijimprop3(GLB_SET,imolty,iatom,1,iatom+1)
            call twh_ijimprop4(GLB_SET,imolty,iatom,1,hmain+heavy+1)
c           --- C-terminus C-Ca-O-O
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,hmain+heavy)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,iatom+3)
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
            atomname(iatom+1) = 'CT2'
            call twh_qqatom(GLB_GET,imolty,iatom+1,charge)
            charge = charge - 0.09d0
            call twh_qqatom(GLB_SET,imolty,iatom+1,charge)
         endif
      endif
      
c     write(6,*) 'finished POLYCHARMM22'
      
      return
      end
#else
c     Compilation disabled (for speed) with FFLAGS=-DDISABLE_POLY
c     statement in ./configure
      write(6,*) "polycharmm22 disabled (compile flag)"
      lfailure = .true.
      return
      end
#endif
