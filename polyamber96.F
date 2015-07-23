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
      subroutine twh_polyamber96(lfailure,imolty,iatom,termflag,total
     &,ibondend,ipeptide,ndisulfide,myaminonum,myatomnum,iha,icb
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
c     --- variables passed to/from the subroutine
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

c      write(6,*) 'start POLYAMBER96'
c     --- determine which amino acid is built this step
c     --- check the special cases first, otherwise this is a generic
c     --- amino acid on a generic backbone.
      if ( pepname(ipeptide) .eq. 'za' ) then
c        --- acetylated N-terminus
c        --- set the bond connectivity information
         call twh_aa_pepname_za(lfailure,aminoname,shortname,termflag
     &        ,imolty,total,iatom)
         if ( lfailure ) return
c        --- References for this group:
c        --- taken from amber_all_amino94.in
         atomname(iatom) = 'C'
         call twh_qqatom(GLB_SET,imolty,iatom,0.597200d0)
         atomname(iatom+1) = 'CT'
         call twh_qqatom(GLB_SET,imolty,iatom+1,-0.366200d0)
         atomname(iatom+2) = 'O'
         call twh_qqatom(GLB_SET,imolty,iatom+2,-0.567900d0)
         atomname(iatom+3) = 'HC'
         call twh_qqatom(GLB_SET,imolty,iatom+3,0.112300d0)
         atomname(iatom+4) = 'HC'
         call twh_qqatom(GLB_SET,imolty,iatom+4,0.112300d0)
         atomname(iatom+5) = 'HC'
         call twh_qqatom(GLB_SET,imolty,iatom+5,0.112300d0)
c        --- improper torsion (in amber CH3 +N C O)
         call twh_inimprop(GLB_SET,imolty,iatom,1)
         call twh_ijimprop2(GLB_SET,imolty,iatom,1,iatom+2)
         call twh_ijimprop3(GLB_SET,imolty,iatom,1,iatom+1)
         call twh_ijimprop4(GLB_SET,imolty,iatom,1,iatom+total)
c        --- update ibondend for ACE
         ibondend = iatom
      elseif ( pepname(ipeptide) .eq. 'zm' ) then
c        --- amidated c-terminus
         call twh_aa_pepname_zm(lfailure,aminoname,shortname,termflag
     &        ,imolty,total,iatom,ibondend)
         if ( lfailure ) return
c        --- References for this group:
c        --- taken from all_aminoct94.lib
         atomname(iatom) = 'N'
         call twh_qqatom(GLB_SET,imolty,iatom,-0.4630d0)
         atomname(iatom+1) = 'H'
         call twh_qqatom(GLB_SET,imolty,iatom+1,0.2315d0)
         atomname(iatom+2) = 'H'
         call twh_qqatom(GLB_SET,imolty,iatom+2,0.2315d0)
c        --- improper torsion (in amber -C H1 N H2)
         call twh_inimprop(GLB_SET,imolty,iatom,1)
         call twh_ijimprop2(GLB_SET,imolty,iatom,1,iatom+2)
         call twh_ijimprop3(GLB_SET,imolty,iatom,1,ibondend)
         call twh_ijimprop4(GLB_SET,imolty,iatom,1,iatom+1)
c        --- ibondend not relevent for c-terminus
         ibondend = -1
      elseif ( pepname(ipeptide) .eq. 'p0' .or.
     &        pepname(ipeptide) .eq. 'P0' ) then
c        --- set up connectivity information
         call twh_aa_pepname_p0(lfailure,aminoname,shortname,termflag
     &        ,imolty,total,iatom,icb,iha,ibondend)
         if ( lfailure ) return
c        --- set up force-field specific information
         if ( termflag .eq. -1 ) then
            atomname(iatom) = 'N3'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.20200d0)
            atomname(iatom+1) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.10000d0)
            atomname(iatom+2) = 'C'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.52600d0)
            atomname(iatom+3) = 'O'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.50000d0)
            atomname(iatom+4) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+4,-0.11500d0)
            atomname(iatom+5) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+5,-0.12100d0)
            atomname(iatom+6) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+6,-0.01200d0)
            atomname(iatom+7) = 'H'
            call twh_qqatom(GLB_SET,imolty,iatom+7,0.31200d0)
            atomname(iatom+8) = 'H'
            call twh_qqatom(GLB_SET,imolty,iatom+8,0.31200d0)
            atomname(iatom+9) = 'HP'
            call twh_qqatom(GLB_SET,imolty,iatom+9,0.10000d0)
            atomname(iatom+10) = 'HC'
            call twh_qqatom(GLB_SET,imolty,iatom+10,0.10000d0)
            atomname(iatom+11) = 'HC'
            call twh_qqatom(GLB_SET,imolty,iatom+11,0.10000d0)
            atomname(iatom+12) = 'HC'
            call twh_qqatom(GLB_SET,imolty,iatom+12,0.10000d0)
            atomname(iatom+13) = 'HC'
            call twh_qqatom(GLB_SET,imolty,iatom+13,0.10000d0)
            atomname(iatom+14) = 'HP'
            call twh_qqatom(GLB_SET,imolty,iatom+14,0.10000d0)
            atomname(iatom+15) = 'HP'
            call twh_qqatom(GLB_SET,imolty,iatom+15,0.10000d0)
c           --- improper torsions
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+3)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,iatom+1)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,iatom+total)
         elseif ( termflag .eq. 0 ) then
            atomname(iatom) = 'N'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.25480d0)
            atomname(iatom+1) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+1,-0.02660d0)
            atomname(iatom+2) = 'C'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.58960d0)
            atomname(iatom+3) = 'O'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57480d0)
            atomname(iatom+4) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+4,-0.00700d0)
            atomname(iatom+5) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+5,0.01890d0)
            atomname(iatom+6) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+6,0.01920d0)
            atomname(iatom+7) = 'H1'
            call twh_qqatom(GLB_SET,imolty,iatom+7,0.06410d0)
            atomname(iatom+8) = 'HC'
            call twh_qqatom(GLB_SET,imolty,iatom+8,0.02530d0)
            atomname(iatom+9) = 'HC'
            call twh_qqatom(GLB_SET,imolty,iatom+9,0.02530d0)
            atomname(iatom+10) = 'HC'
            call twh_qqatom(GLB_SET,imolty,iatom+10,0.02130d0)
            atomname(iatom+11) = 'HC'
            call twh_qqatom(GLB_SET,imolty,iatom+11,0.02130d0)
            atomname(iatom+12) = 'H1'
            call twh_qqatom(GLB_SET,imolty,iatom+12,0.03910d0)
            atomname(iatom+13) = 'H1'
            call twh_qqatom(GLB_SET,imolty,iatom+13,0.03910d0)
c           --- improper torsions
            call twh_inimprop(GLB_SET,imolty,iatom,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom,1,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom,1,ibondend)
            call twh_ijimprop4(GLB_SET,imolty,iatom,1,iatom+6)
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+3)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,iatom+1)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,iatom+total)
         elseif ( termflag .eq. 1 ) then
            atomname(iatom) = 'N'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.28020d0)
            atomname(iatom+1) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+1,-0.13360d0)
            atomname(iatom+2) = 'C'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.66310d0)
            atomname(iatom+3) = 'O2'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.76970d0)
            atomname(iatom+4) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+4,-0.05430d0)
            atomname(iatom+5) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+5,0.04660d0)
            atomname(iatom+6) = 'CT'
            call twh_qqatom(GLB_SET,imolty,iatom+6,0.04340d0)
            atomname(iatom+7) = 'O2'
            call twh_qqatom(GLB_SET,imolty,iatom+7,-0.76970d0)
            atomname(iatom+8) = 'H1'
            call twh_qqatom(GLB_SET,imolty,iatom+8,0.07760d0)
            atomname(iatom+9) = 'HC'
            call twh_qqatom(GLB_SET,imolty,iatom+9,0.03810d0)
            atomname(iatom+10) = 'HC'
            call twh_qqatom(GLB_SET,imolty,iatom+10,0.03810d0)
            atomname(iatom+11) = 'HC'
            call twh_qqatom(GLB_SET,imolty,iatom+11,0.01720d0)
            atomname(iatom+12) = 'HC'
            call twh_qqatom(GLB_SET,imolty,iatom+12,0.01720d0)
            atomname(iatom+13) = 'H1'
            call twh_qqatom(GLB_SET,imolty,iatom+13,0.03310d0)
            atomname(iatom+14) = 'H1'
            call twh_qqatom(GLB_SET,imolty,iatom+14,0.03310d0)
c           --- improper torsions
            call twh_inimprop(GLB_SET,imolty,iatom,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom,1,iatom+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom,1,ibondend)
            call twh_ijimprop4(GLB_SET,imolty,iatom,1,iatom+6)
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+3)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,iatom+1)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,iatom+7)
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
c           --- set up atom types and vibrations for side chain
            atomname(hmain) = 'CT'
            atomname(shift) = 'HC'
            atomname(shift+1) = 'HC'
            atomname(shift+2) = 'HC'
c           --- no side chain improper torsions for alanine
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c               --- n-terminus
c               --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.14140d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.09620d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.61630d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57220d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.19970d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.19970d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.19970d0)
c              --- Ca hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.08890d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.05970d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.030d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.030d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.030d0)
            elseif ( termflag .eq. 0 ) then
c               --- main chain
c               --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.41570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.03370d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.59730d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.56790d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27190d0)
c              --- Ca hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.08230d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.18250d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.06030d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.06030d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.06030d0)
            elseif ( termflag .eq. 1 ) then
c               --- c-terminus
c               --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.38210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.17470d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.77310d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.80550d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.80550d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.26810d0)
c              --- Ca hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.10670d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.20930d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.07640d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.07640d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.07640d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'n0' .or.
     &           pepname(ipeptide) .eq. 'N0' ) then
            call twh_aa_pepname_n0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'C'
            atomname(hmain+2) = 'O'
            atomname(hmain+3) = 'N'
            atomname(shift) = 'HC'
            atomname(shift+1) = 'HC'
            atomname(shift+2) = 'H'
            atomname(shift+3) = 'H'
c           --- side chain improper torsions
            call twh_inimprop(GLB_SET,imolty,hmain+1,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+3)
            call twh_inimprop(GLB_SET,imolty,hmain+3,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,1,shift+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,1,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,1,shift+2)
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.18010d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.03680d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.61630d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57220d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.19210d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.19210d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.19210d0)
c               --- Ca hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.12310d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.02830d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.58330d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.57440d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.86340d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.05150d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.05150d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.40970d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.40970d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.41570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.01430d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.59730d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.56790d0)
c               --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27190d0)
c              --- Ca hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.10480d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.20410d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.71300d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.59310d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.91910d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.07970d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.07970d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.41960d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.41960d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.38210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.20800d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.80500d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.81470d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.81470d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.26810d0)
c              --- Ca hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.13580d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.22990d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.71530d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.60100d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.90840d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.10230d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.10230d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.41500d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.41500d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'c0' .or.
     &           pepname(ipeptide) .eq. 'C0' ) then
            call twh_aa_pepname_c0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'SH'
            atomname(shift) = 'H1'
            atomname(shift+1) = 'H1'
            atomname(shift+2) = 'HS'
c           --- no side chain improper torsions
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.13250d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.09270d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.61230d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57130d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.20230d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.20230d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.20230d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.14110d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.11950d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.32980d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.11880d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.11880d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.19750d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.41570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.02130d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.59730d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.56790d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27190d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.11240d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.12310d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.31190d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.11120d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.11120d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.19330d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.38210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.16350d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.74970d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.79810d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.79810d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.26810d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.13960d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.19960d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.31020d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.14370d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.14370d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.20680d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'cs' .or.
     &           pepname(ipeptide) .eq. 'CS' .or.
     &           pepname(ipeptide) .eq. 'cS' .or.
     &           pepname(ipeptide) .eq. 'Cs' ) then
            call twh_aa_pepname_cs(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'S'
            atomname(shift) = 'H1'
            atomname(shift+1) = 'H1'
c           --- setup information to figure out the disulfide bond
            ndisulfide = ndisulfide + 1
            myaminonum(ndisulfide) = ipeptide
            myatomnum(ndisulfide) = hmain+1
c           --- no side chain improper torsions
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.20690d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.10550d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.61230d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57130d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.18150d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.18150d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.18150d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.09220d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.02770d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.09840d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.06800d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.06800d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.41570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.04290d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.59730d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.56790d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27190d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.07660d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.07900d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.10810d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.09100d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.09100d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.38210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.13180d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.76180d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.80410d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.80410d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.26810d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.09380d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.19430d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.05290d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.12280d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.12280d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'q0' .or.
     &           pepname(ipeptide) .eq. 'Q0' ) then
            call twh_aa_pepname_q0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'CT'
            atomname(hmain+2) = 'C'
            atomname(hmain+3) = 'O'
            atomname(hmain+4) = 'N'
            atomname(shift) = 'HC'
            atomname(shift+1) = 'HC'
            atomname(shift+2) = 'HC'
            atomname(shift+3) = 'HC'
            atomname(shift+4) = 'H'
            atomname(shift+5) = 'H'
c           --- side chain improper torsions
            call twh_inimprop(GLB_SET,imolty,hmain+2,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,hmain+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,hmain+4)
            call twh_inimprop(GLB_SET,imolty,hmain+4,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,shift+5)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,shift+4)
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.14930d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.05360d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.61230d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57130d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.19960d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.19960d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.19960d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.10150d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.06510d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.09030d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,0.73540d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.61330d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-1.00310d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.00500d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.00500d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.03310d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.03310d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.44290d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.44290d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.41570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.00310d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.59730d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.56790d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27190d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.08500d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.00360d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.06450d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,0.69510d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.60860d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.94070d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.01710d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.01710d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.03520d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.03520d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.42510d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.42510d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.38210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.22480d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.77750d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.80420d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.80420d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.26810d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.12320d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.06640d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.02100d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,0.70930d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.60980d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.95740d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.04520d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.04520d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.02030d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.02030d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.43040d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.43040d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'g0' .or.
     &           pepname(ipeptide) .eq. 'G0' ) then
            call twh_aa_pepname_g0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
c           --- the atomname for the hydrogen in glycine depends
c           --- on the termflag
c           --- no improper torsions on side chain
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
               atomname(shift) = 'HP'
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.29430d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.01000d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.61630d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57220d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.16420d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.16420d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.16420d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.08950d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,shift,0.08950d0)
            elseif ( termflag .eq. 0 ) then
               atomname(shift) = 'H1'
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.41570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.02520d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.59730d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.56790d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27190d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.06980d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,shift,0.06980d0)
            elseif ( termflag .eq. 1 ) then
               atomname(shift) = 'H1'
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.38210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.24930d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.72310d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.78550d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.78550d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.26810d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.10560d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,shift,0.10560d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'i0' .or.
     &           pepname(ipeptide) .eq. 'I0' ) then
            call twh_aa_pepname_i0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'CT'
            atomname(hmain+2) = 'CT'
            atomname(hmain+3) = 'CT'
            atomname(shift) = 'HC'
            atomname(shift+1) = 'HC'
            atomname(shift+2) = 'HC'
            atomname(shift+3) = 'HC'
            atomname(shift+4) = 'HC'
            atomname(shift+5) = 'HC'
            atomname(shift+6) = 'HC'
            atomname(shift+7) = 'HC'
            atomname(shift+8) = 'HC'
c           --- no side chain improper torsions
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.03110d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.02570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.61230d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57130d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.23290d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.23290d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.23290d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.10310d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.18850d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.03870d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.37200d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.09080d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.02130d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.02010d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.02010d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.09470d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.09470d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.09470d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.02260d0)
               call twh_qqatom(GLB_SET,imolty,shift+7,0.02260d0)
               call twh_qqatom(GLB_SET,imolty,shift+8,0.02260d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.41570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.05970d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.59730d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.56790d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27190d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.08690d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.13030d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.04300d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.32040d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.06600d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.01870d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.02360d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.02360d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.08820d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.08820d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.08820d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.01860d0)
               call twh_qqatom(GLB_SET,imolty,shift+7,0.01860d0)
               call twh_qqatom(GLB_SET,imolty,shift+8,0.01860d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.38210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.31000d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.83430d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.81900d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.81900d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.26810d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.13750d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.03630d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.03230d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.34980d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.06990d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.07660d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.03210d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.03210d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.10210d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.10210d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.10210d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.01960d0)
               call twh_qqatom(GLB_SET,imolty,shift+7,0.01960d0)
               call twh_qqatom(GLB_SET,imolty,shift+8,0.01960d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'l0' .or.
     &           pepname(ipeptide) .eq. 'L0' ) then
            call twh_aa_pepname_l0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'CT'
            atomname(hmain+2) = 'CT'
            atomname(hmain+3) = 'CT'
            atomname(shift) = 'HC'
            atomname(shift+1) = 'HC'
            atomname(shift+2) = 'HC'
            atomname(shift+3) = 'HC'
            atomname(shift+4) = 'HC'
            atomname(shift+5) = 'HC'
            atomname(shift+6) = 'HC'
            atomname(shift+7) = 'HC'
            atomname(shift+8) = 'HC'
c           --- no side chain impropers
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.10100d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.01040d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.61230d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57130d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.21480d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.21480d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.21480d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.10530d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.02440d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.34210d0)
c              --- technically the charges on the two carbons
c              --- below are listed as -0.41060 and -0.41040 but
c              --- it seemed absurd to use different charges on
c              --- symmetric atoms so I took the average
c              --- 03-16-2001 M.G. Martin
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.41050d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.41050d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.02560d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.02560d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,-0.03800d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.09800d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.09800d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.09800d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.09800d0)
               call twh_qqatom(GLB_SET,imolty,shift+7,0.09800d0)
               call twh_qqatom(GLB_SET,imolty,shift+8,0.09800d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.41570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.05180d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.59730d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.56790d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27190d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.09220d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.11020d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.35310d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.41210d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.41210d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.04570d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.04570d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,-0.03610d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.10000d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.10000d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.10000d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.10000d0)
               call twh_qqatom(GLB_SET,imolty,shift+7,0.10000d0)
               call twh_qqatom(GLB_SET,imolty,shift+8,0.10000d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.38210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.28470d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.83260d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.81990d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.81990d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.26810d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.13460d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.24690d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.37060d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.41630d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.41630d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.09740d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.09740d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,-0.03740d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.10380d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.10380d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.10380d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.10380d0)
               call twh_qqatom(GLB_SET,imolty,shift+7,0.10380d0)
               call twh_qqatom(GLB_SET,imolty,shift+8,0.10380d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'm0' .or.
     &           pepname(ipeptide) .eq. 'M0' ) then
            call twh_aa_pepname_m0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'CT'
            atomname(hmain+2) = 'S'
            atomname(hmain+3) = 'CT'
            atomname(shift) = 'HC'
            atomname(shift+1) = 'HC'
            atomname(shift+2) = 'H1'
            atomname(shift+3) = 'H1'
            atomname(shift+4) = 'H1'
            atomname(shift+5) = 'H1'
            atomname(shift+6) = 'H1'
c           --- no impropers on side chain
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.15920d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.02210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.61230d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57130d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.19840d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.19840d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.19840d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.11160d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.08650d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.03340d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.27740d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.03410d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.01250d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.01250d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.02920d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.02920d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.05970d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.05970d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.05970d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.41570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.02370d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.59730d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.56790d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27190d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.08800d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.03420d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.00180d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.27370d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.05360d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.02410d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.02410d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.04400d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.04400d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.06840d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.06840d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.06840d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.38210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.25970d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.80130d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.81050d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.81050d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.26810d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.12770d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.02360d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.04920d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.26920d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.03760d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.04800d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.04800d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.03170d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.03170d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.06250d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.06250d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.06250d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'f0' .or.
     &           pepname(ipeptide) .eq. 'F0' ) then
            call twh_aa_pepname_f0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'CA'
            atomname(hmain+2) = 'CA'
            atomname(hmain+3) = 'CA'
            atomname(hmain+4) = 'CA'
            atomname(hmain+5) = 'CA'
            atomname(hmain+6) = 'CA'
            atomname(shift) = 'HC'
            atomname(shift+1) = 'HC'
            atomname(shift+2) = 'HA'
            atomname(shift+3) = 'HA'
            atomname(shift+4) = 'HA'
            atomname(shift+5) = 'HA'
            atomname(shift+6) = 'HA'
c           --- side chain improper torsions
            call twh_inimprop(GLB_SET,imolty,hmain+1,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+3)
            call twh_inimprop(GLB_SET,imolty,hmain+2,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,shift+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,hmain+4)
            call twh_inimprop(GLB_SET,imolty,hmain+3,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,1,shift+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,1,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,1,hmain+5)
            call twh_inimprop(GLB_SET,imolty,hmain+4,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,shift+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+6)
            call twh_inimprop(GLB_SET,imolty,hmain+5,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,1,shift+5)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,1,hmain+6)
            call twh_inimprop(GLB_SET,imolty,hmain+6,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+6,1,shift+6)
            call twh_ijimprop3(GLB_SET,imolty,hmain+6,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+6,1,hmain+5)
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.17370d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.07330d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.61230d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57130d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.19210d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.19210d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.19210d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.10410d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.03300d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.00310d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.13920d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.13910d0)
c              --- again the amber people distinguish between 
c              --- indistinguishible atoms
c              --- they list slightly different values for ce1 &
c              --- ce2.  I take the average for both as they are
c              --- identical 05-15-2001 M.G. Martin
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.16025d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.16025d0)
               call twh_qqatom(GLB_SET,imolty,hmain+6,-0.12080d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.01040d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.01040d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.13740d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.13740d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.14330d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.14330d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.13290d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.41570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.00240d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.59730d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.56790d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27190d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.09780d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.03430d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.01180d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.12560d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.12560d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.17040d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.17040d0)
               call twh_qqatom(GLB_SET,imolty,hmain+6,-0.10720d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.02950d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.02950d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.13300d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.13300d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.14300d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.14300d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.12970d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.38210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.18250d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.76600d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.80260d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.80260d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.26810d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.10980d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.09590d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.05520d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.13000d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.13000d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.18470d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.18470d0)
               call twh_qqatom(GLB_SET,imolty,hmain+6,-0.09440d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.04430d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.04430d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.14080d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.14080d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.14610d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.14610d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.12800d0)
            endif
         elseif ( pepname(ipeptide) .eq. 's0' .or.
     &           pepname(ipeptide) .eq. 'S0' ) then
            call twh_aa_pepname_s0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'OH'
            atomname(shift) = 'H1'
            atomname(shift+1) = 'H1'
            atomname(shift+2) = 'HO'
c           --- no side chain improper torsions
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.18490d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.05670d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.61630d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57220d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.18980d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.18980d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.18980d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.07820d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.25960d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.67140d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.02730d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.02730d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.42390d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.41570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.02490d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.59730d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.56790d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27190d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.08430d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.21170d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.65460d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.03520d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.03520d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.42750d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.38210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.27220d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.81130d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.81320d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.81320d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.26810d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.13040d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.11230d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.65140d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.08130d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.08130d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.44740d0)
            endif
         elseif ( pepname(ipeptide) .eq. 't0' .or.
     &           pepname(ipeptide) .eq. 'T0' ) then
            call twh_aa_pepname_t0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'OH'
            atomname(hmain+2) = 'CT'
            atomname(shift) = 'H1'
            atomname(shift+1) = 'HO'
            atomname(shift+2) = 'HC'
            atomname(shift+3) = 'HC'
            atomname(shift+4) = 'HC'
c           --- no side chain improper torsions
c            --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.18120d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.00340d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.61630d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57220d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.19340d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.19340d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.19340d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.10870d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.45140d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.67640d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.25540d0)
               call twh_qqatom(GLB_SET,imolty,shift,-0.03230d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.40700d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.06270d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.06270d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.06270d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.41570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.03890d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.59730d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.56790d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27190d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.10070d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.36540d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.67610d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.24380d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.00430d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.41020d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.06420d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.06420d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.06420d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.38210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.24200d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.78100d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.80440d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.80440d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.26810d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.12070d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.30250d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.64960d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.18530d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.00780d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.41190d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.05860d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.05860d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.05860d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'w0' .or.
     &           pepname(ipeptide) .eq. 'W0' ) then
            call twh_aa_pepname_w0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'C*'
            atomname(hmain+2) = 'CW'
            atomname(hmain+3) = 'CB'
            atomname(hmain+4) = 'NA'
            atomname(hmain+5) = 'CN'
            atomname(hmain+6) = 'CA'
            atomname(hmain+7) = 'CA'
            atomname(hmain+8) = 'CA'
            atomname(hmain+9) = 'CA'
            atomname(shift) = 'HC'
            atomname(shift+1) = 'HC'
            atomname(shift+2) = 'H4'
            atomname(shift+3) = 'H'
            atomname(shift+4) = 'HA'
            atomname(shift+5) = 'HA'
            atomname(shift+6) = 'HA'
            atomname(shift+7) = 'HA'
c           --- side chain improper torsions
            call twh_inimprop(GLB_SET,imolty,hmain+1,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+3)
            call twh_inimprop(GLB_SET,imolty,hmain+2,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,shift+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,hmain+4)
            call twh_inimprop(GLB_SET,imolty,hmain+4,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,shift+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+5)
            call twh_inimprop(GLB_SET,imolty,hmain+6,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+6,1,shift+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+6,1,hmain+8)
            call twh_ijimprop4(GLB_SET,imolty,hmain+6,1,hmain+3)
            call twh_inimprop(GLB_SET,imolty,hmain+7,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+7,1,shift+5)
            call twh_ijimprop3(GLB_SET,imolty,hmain+7,1,hmain+5)
            call twh_ijimprop4(GLB_SET,imolty,hmain+7,1,hmain+9)
            call twh_inimprop(GLB_SET,imolty,hmain+8,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+8,1,shift+6)
            call twh_ijimprop3(GLB_SET,imolty,hmain+8,1,hmain+9)
            call twh_ijimprop4(GLB_SET,imolty,hmain+8,1,hmain+6)
            call twh_inimprop(GLB_SET,imolty,hmain+9,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+9,1,shift+7)
            call twh_ijimprop3(GLB_SET,imolty,hmain+9,1,hmain+7)
            call twh_ijimprop4(GLB_SET,imolty,hmain+9,1,hmain+8)
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.19130d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.04210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.61230d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57130d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.18880d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.18880d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.18880d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.11620d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.05430d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.16540d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.17880d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,0.11320d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.34440d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,0.15750d0)
               call twh_qqatom(GLB_SET,imolty,hmain+6,-0.22650d0)
               call twh_qqatom(GLB_SET,imolty,hmain+7,-0.27100d0)
               call twh_qqatom(GLB_SET,imolty,hmain+8,-0.20340d0)
               call twh_qqatom(GLB_SET,imolty,hmain+9,-0.10800d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.02220d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.02220d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.21950d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.34120d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.16460d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.15890d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.14580d0)
               call twh_qqatom(GLB_SET,imolty,shift+7,0.14110d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.41570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.02750d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.59730d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.56790d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27190d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.11230d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.00500d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.14150d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.16380d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,0.12430d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.34180d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,0.13800d0)
               call twh_qqatom(GLB_SET,imolty,hmain+6,-0.23870d0)
               call twh_qqatom(GLB_SET,imolty,hmain+7,-0.26010d0)
               call twh_qqatom(GLB_SET,imolty,hmain+8,-0.19720d0)
               call twh_qqatom(GLB_SET,imolty,hmain+9,-0.11340d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.03390d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.03390d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.20620d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.34120d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.17000d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.15720d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.14470d0)
               call twh_qqatom(GLB_SET,imolty,shift+7,0.14170d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.38210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.20840d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.76580d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.80110d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.80110d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.26810d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.12720d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.07420d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.07960d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.18080d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,0.10780d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.33160d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,0.12220d0)
               call twh_qqatom(GLB_SET,imolty,hmain+6,-0.18370d0)
               call twh_qqatom(GLB_SET,imolty,hmain+7,-0.25940d0)
               call twh_qqatom(GLB_SET,imolty,hmain+8,-0.22870d0)
               call twh_qqatom(GLB_SET,imolty,hmain+9,-0.10200d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.04970d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.04970d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.20430d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.34130d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.14910d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.15670d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.15070d0)
               call twh_qqatom(GLB_SET,imolty,shift+7,0.14010d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'y0' .or.
     &           pepname(ipeptide) .eq. 'Y0' ) then
            call twh_aa_pepname_y0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'CA'
            atomname(hmain+2) = 'CA'
            atomname(hmain+3) = 'CA'
            atomname(hmain+4) = 'CA'
            atomname(hmain+5) = 'CA'
            atomname(hmain+6) = 'C'
            atomname(hmain+7) = 'OH'
            atomname(shift) = 'HC'
            atomname(shift+1) = 'HC'
            atomname(shift+2) = 'HA'
            atomname(shift+3) = 'HA'
            atomname(shift+4) = 'HA'
            atomname(shift+5) = 'HA'
            atomname(shift+6) = 'HO'
c           --- side chain improper torsions
            call twh_inimprop(GLB_SET,imolty,hmain+1,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+3)
            call twh_inimprop(GLB_SET,imolty,hmain+2,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,shift+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,hmain+4)
            call twh_inimprop(GLB_SET,imolty,hmain+3,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,1,shift+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,1,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,1,hmain+5)
            call twh_inimprop(GLB_SET,imolty,hmain+4,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,shift+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+6)
            call twh_inimprop(GLB_SET,imolty,hmain+5,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,1,shift+5)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,1,hmain+6)
            call twh_inimprop(GLB_SET,imolty,hmain+6,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+6,1,hmain+7)
            call twh_ijimprop3(GLB_SET,imolty,hmain+6,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+6,1,hmain+5)
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.19400d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.05700d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.61230d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57130d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.18730d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.18730d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.18730d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.09830d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.06590d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.02050d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.20020d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.20020d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.22390d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.22390d0)
               call twh_qqatom(GLB_SET,imolty,hmain+6,0.31390d0)
               call twh_qqatom(GLB_SET,imolty,hmain+7,-0.55780d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.01020d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.01020d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.17200d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.17200d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.16500d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.16500d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.40010d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.41570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.00140d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.59730d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.56790d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27190d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.08760d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.01520d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.00110d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.19060d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.19060d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.23410d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.23410d0)
               call twh_qqatom(GLB_SET,imolty,hmain+6,0.32260d0)
               call twh_qqatom(GLB_SET,imolty,hmain+7,-0.55790d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.02950d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.02950d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.16990d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.16990d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.16560d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.16560d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.39920d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.38210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.20150d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.78170d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.80700d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.80700d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.26810d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.10920d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.07520d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.02430d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.19220d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.19220d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.24580d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.24580d0)
               call twh_qqatom(GLB_SET,imolty,hmain+6,0.33950d0)
               call twh_qqatom(GLB_SET,imolty,hmain+7,-0.56430d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.04900d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.04900d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.17800d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.17800d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.16730d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.16730d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.40170d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'v0' .or.
     &           pepname(ipeptide) .eq. 'V0' ) then
            call twh_aa_pepname_v0(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'CT'
            atomname(hmain+2) = 'CT'
            atomname(shift) = 'HC'
            atomname(shift+1) = 'HC'
            atomname(shift+2) = 'HC'
            atomname(shift+3) = 'HC'
            atomname(shift+4) = 'HC'
            atomname(shift+5) = 'HC'
            atomname(shift+6) = 'HC'
c           --- no side chain improper torsions
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.05770d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.00540d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.61630d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57220d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.22720d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.22720d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.22720d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.10930d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.31960d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.31290d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.31290d0)
               call twh_qqatom(GLB_SET,imolty,shift,-0.02210d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.07350d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.07350d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.07350d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.07350d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.07350d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.07350d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.41570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.08750d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.59730d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.56790d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27190d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.09690d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.29850d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.31920d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.31920d0)
               call twh_qqatom(GLB_SET,imolty,shift,-0.02970d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.07910d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.07910d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.07910d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.07910d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.07910d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.07910d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.38210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.34380d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.83500d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.81730d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.81730d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.26810d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.14380d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.19400d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.30640d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.30640d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.03080d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.08360d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.08360d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.08360d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.08360d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.08360d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.08360d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'd-' .or.
     &           pepname(ipeptide) .eq. 'D-' ) then
            call twh_aa_pepname_d_neg(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'C'
            atomname(hmain+2) = 'O2'
            atomname(hmain+3) = 'O2'
            atomname(shift) = 'HC'
            atomname(shift+1) = 'HC'
c           --- side chain improper torsion
            call twh_inimprop(GLB_SET,imolty,hmain+1,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+2)
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.07820d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.02920d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.56210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.58890d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.22000d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.22000d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.22000d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.11410d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.02350d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.81940d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.80840d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.80840d0)
               call twh_qqatom(GLB_SET,imolty,shift,-0.01690d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,-0.01690d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.51630d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.03810d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.53660d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.58190d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.29360d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.08800d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.03030d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.79940d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.80140d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.80140d0)
               call twh_qqatom(GLB_SET,imolty,shift,-0.01220d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,-0.01220d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.51920d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.18170d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.72560d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.78870d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.78870d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.30550d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.10460d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.06770d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.88510d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.81620d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.81620d0)
               call twh_qqatom(GLB_SET,imolty,shift,-0.02120d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,-0.02120d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'e-' .or.
     &           pepname(ipeptide) .eq. 'E-' ) then
            call twh_aa_pepname_e_neg(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'CT'
            atomname(hmain+2) = 'C'
            atomname(hmain+3) = 'O2'
            atomname(hmain+4) = 'O2'
            atomname(shift) = 'HC'
            atomname(shift+1) = 'HC'
            atomname(shift+2) = 'HC'
            atomname(shift+3) = 'HC'
c           --- side chain improper torsion
            call twh_inimprop(GLB_SET,imolty,hmain+2,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,hmain+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,hmain+3)
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.00170d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.05880d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.56210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.58890d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.23910d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.23910d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.23910d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.12020d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.09090d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.02360d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,0.80870d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.81890d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.81890d0)
               call twh_qqatom(GLB_SET,imolty,shift,-0.02320d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,-0.02320d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,-0.03150d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,-0.03150d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.51630d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.03970d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.53660d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.58190d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.29360d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.11050d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.05600d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.01360d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,0.80540d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.81880d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.81880d0)
               call twh_qqatom(GLB_SET,imolty,shift,-0.01730d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,-0.01730d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,-0.04250d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,-0.04250d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.51920d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.20590d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.74200d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.79300d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.79300d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.30550d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.13990d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.00710d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.06750d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,0.81830d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.82200d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.82200d0)
               call twh_qqatom(GLB_SET,imolty,shift,-0.00780d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,-0.00780d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,-0.05480d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,-0.05480d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'r+' .or.
     &           pepname(ipeptide) .eq. 'R+' ) then
            call twh_aa_pepname_r_plus(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'CT'
            atomname(hmain+2) = 'CT'
            atomname(hmain+3) = 'N2'
            atomname(hmain+4) = 'CA'
            atomname(hmain+5) = 'N2'
            atomname(hmain+6) = 'N2'
            atomname(shift) = 'HC'
            atomname(shift+1) = 'HC'
            atomname(shift+2) = 'HC'
            atomname(shift+3) = 'HC'
            atomname(shift+4) = 'H1'
            atomname(shift+5) = 'H1'
            atomname(shift+6) = 'H'
            atomname(shift+7) = 'H'
            atomname(shift+8) = 'H'
            atomname(shift+9) = 'H'
            atomname(shift+10) = 'H'
c           --- side chain improper torsions
            call twh_inimprop(GLB_SET,imolty,hmain+3,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,1,shift+6)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,1,hmain+4)
            call twh_inimprop(GLB_SET,imolty,hmain+4,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,hmain+6)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+3)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+5)
            call twh_inimprop(GLB_SET,imolty,hmain+5,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,1,shift+8)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,1,shift+7)
            call twh_inimprop(GLB_SET,imolty,hmain+6,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+6,1,shift+10)
            call twh_ijimprop3(GLB_SET,imolty,hmain+6,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+6,1,shift+9)
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.13050d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.02230d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.72140d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.60130d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.20830d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.20830d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.20830d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.12420d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.01180d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.02360d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,0.09350d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.56500d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,0.82810d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.86930d0)
               call twh_qqatom(GLB_SET,imolty,hmain+6,-0.86930d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.02260d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.02260d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.03090d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.03090d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.05270d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.05270d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.35920d0)
               call twh_qqatom(GLB_SET,imolty,shift+7,0.44940d0)
               call twh_qqatom(GLB_SET,imolty,shift+8,0.44940d0)
               call twh_qqatom(GLB_SET,imolty,shift+9,0.44940d0)
               call twh_qqatom(GLB_SET,imolty,shift+10,0.44940d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.34790d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.26370d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.73410d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.58940d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27470d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.15600d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.00070d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.03900d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,0.04860d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.52950d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,0.80760d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.86270d0)
               call twh_qqatom(GLB_SET,imolty,hmain+6,-0.86270d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.03270d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.03270d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.02850d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.02850d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.06870d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.06870d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.34560d0)
               call twh_qqatom(GLB_SET,imolty,shift+7,0.44780d0)
               call twh_qqatom(GLB_SET,imolty,shift+8,0.44780d0)
               call twh_qqatom(GLB_SET,imolty,shift+9,0.44780d0)
               call twh_qqatom(GLB_SET,imolty,shift+10,0.44780d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.34810d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.30680d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.85570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.82660d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.82660d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.27640d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.14470d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.03740d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.07440d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,0.11140d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.55640d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,0.83680d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.87370d0)
               call twh_qqatom(GLB_SET,imolty,hmain+6,-0.87370d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.03710d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.03710d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.01850d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.01850d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.04680d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.04680d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.34790d0)
               call twh_qqatom(GLB_SET,imolty,shift+7,0.44930d0)
               call twh_qqatom(GLB_SET,imolty,shift+8,0.44930d0)
               call twh_qqatom(GLB_SET,imolty,shift+9,0.44930d0)
               call twh_qqatom(GLB_SET,imolty,shift+10,0.44930d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'h+' .or.
     &           pepname(ipeptide) .eq. 'H+' ) then
            call twh_aa_pepname_h_plus(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'CC'
            atomname(hmain+2) = 'NA'
            atomname(hmain+3) = 'CW'
            atomname(hmain+4) = 'CR'
            atomname(hmain+5) = 'NA'
            atomname(shift) = 'HC'
            atomname(shift+1) = 'HC'
            atomname(shift+2) = 'H'
            atomname(shift+3) = 'H4'
            atomname(shift+4) = 'H5'
            atomname(shift+5) = 'H'
c           --- side chain improper torsions
            call twh_inimprop(GLB_SET,imolty,hmain+1,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+3)
            call twh_inimprop(GLB_SET,imolty,hmain+2,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,shift+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,hmain+4)
            call twh_inimprop(GLB_SET,imolty,hmain+3,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,1,shift+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,1,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,1,hmain+5)
            call twh_inimprop(GLB_SET,imolty,hmain+4,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,shift+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+5)
            call twh_inimprop(GLB_SET,imolty,hmain+5,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,1,shift+5)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,1,hmain+3)
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.25600d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.05810d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.72140d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.60130d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.17040d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.17040d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.17040d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.10470d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.04840d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.02360d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.15100d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.14330d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.00110d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.17390d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.05310d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.05310d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.38210d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.24950d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.26450d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.39210d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.34790d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.13540d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.73410d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.58940d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27470d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.12120d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.04140d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.00120d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.15130d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.11410d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.01700d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.17180d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.08100d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.08100d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.38660d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.23170d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.26810d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.39110d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.34810d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.14450d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.80320d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.81770d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.81770d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.27640d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.11150d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.08000d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.02980d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.15010d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.12560d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.02510d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.16830d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.08680d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.08680d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.38830d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.23360d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.26940d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.39130d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'hd' .or.
     &           pepname(ipeptide) .eq. 'HD' .or.
     &           pepname(ipeptide) .eq. 'Hd' ) then
            call twh_aa_pepname_hd(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'CC'
            atomname(hmain+2) = 'NA'
            atomname(hmain+3) = 'CV'
            atomname(hmain+4) = 'CR'
            atomname(hmain+5) = 'NB'
            atomname(shift) = 'HC'
            atomname(shift+1) = 'HC'
            atomname(shift+2) = 'H'
            atomname(shift+3) = 'H4'
            atomname(shift+4) = 'H5'
c           --- side chain improper torsions
            call twh_inimprop(GLB_SET,imolty,hmain+1,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+3)
            call twh_inimprop(GLB_SET,imolty,hmain+2,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+2,1,shift+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+2,1,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+2,1,hmain+4)
            call twh_inimprop(GLB_SET,imolty,hmain+3,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,1,shift+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,1,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,1,hmain+5)
            call twh_inimprop(GLB_SET,imolty,hmain+4,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,shift+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+5)
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.15420d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.09640d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.61230d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57130d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.19630d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.19630d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.19630d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.09580d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.02590d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.03990d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.38190d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,0.10460d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,0.21270d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.57110d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.02090d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.02090d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.36320d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.12990d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.13850d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.41570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.01880d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.59730d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.56790d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27190d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.08810d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.04620d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.02660d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.38110d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,0.12920d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,0.20570d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.57270d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.04020d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.04020d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.36490d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.11470d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.13920d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.38210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.17390d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.76150d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.80160d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.80160d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.26810d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.11000d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.10460d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.02930d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.38920d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,0.10010d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,0.19250d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.56290d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.05650d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.05650d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.37550d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.12410d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.14180d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'he' .or.
     &           pepname(ipeptide) .eq. 'HE' .or.
     &           pepname(ipeptide) .eq. 'He' ) then
            call twh_aa_pepname_he(lfailure,aminoname,shortname,imolty
     &           ,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'CC'
            atomname(hmain+2) = 'NB'
            atomname(hmain+3) = 'CW'
            atomname(hmain+4) = 'CR'
            atomname(hmain+5) = 'NA'
            atomname(shift) = 'HC'
            atomname(shift+1) = 'HC'
            atomname(shift+2) = 'H4'
            atomname(shift+3) = 'H5'
            atomname(shift+4) = 'H'
c           --- side chain improper torsions
            call twh_inimprop(GLB_SET,imolty,hmain+1,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain)
            call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+3)
            call twh_inimprop(GLB_SET,imolty,hmain+3,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+3,1,shift+2)
            call twh_ijimprop3(GLB_SET,imolty,hmain+3,1,hmain+1)
            call twh_ijimprop4(GLB_SET,imolty,hmain+3,1,hmain+5)
            call twh_inimprop(GLB_SET,imolty,hmain+4,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,shift+3)
            call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+2)
            call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+5)
            call twh_inimprop(GLB_SET,imolty,hmain+5,1)
            call twh_ijimprop2(GLB_SET,imolty,hmain+5,1,shift+4)
            call twh_ijimprop3(GLB_SET,imolty,hmain+5,1,hmain+4)
            call twh_ijimprop4(GLB_SET,imolty,hmain+5,1,hmain+3)
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.14720d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,0.02360d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.61230d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57130d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.20160d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.20160d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.20160d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.13800d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.04890d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.17400d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.55790d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.23490d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,0.18040d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.27810d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.02230d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.02230d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.19630d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.13970d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.33240d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.41570d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.05810d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.59730d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.56790d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27190d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.13600d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.00740d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.18680d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.54320d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.22070d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,0.16350d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.27950d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.03670d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.03670d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.18620d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.14350d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.33390d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.38210d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.26990d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.79160d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.80650d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.80650d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.26810d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.16500d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.10680d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.27240d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.55170d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.25880d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,0.15580d0)
               call twh_qqatom(GLB_SET,imolty,hmain+5,-0.26700d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.06200d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.06200d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.19570d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.14480d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.33190d0)
            endif
         elseif ( pepname(ipeptide) .eq. 'k+' .or.
     &           pepname(ipeptide) .eq. 'K+' ) then
            call twh_aa_pepname_k_plus(lfailure,aminoname,shortname
     &           ,imolty,iatom,hmain,heavy,light,shift,ninsert)
            if ( lfailure ) return
            atomname(hmain) = 'CT'
            atomname(hmain+1) = 'CT'
            atomname(hmain+2) = 'CT'
            atomname(hmain+3) = 'CT'
            atomname(hmain+4) = 'N3'
            atomname(shift) = 'HC'
            atomname(shift+1) = 'HC'
            atomname(shift+2) = 'HC'
            atomname(shift+3) = 'HC'
            atomname(shift+4) = 'HC'
            atomname(shift+5) = 'HC'
            atomname(shift+6) = 'HP'
            atomname(shift+7) = 'HP'
            atomname(shift+8) = 'H'
            atomname(shift+9) = 'H'
            atomname(shift+10) = 'H'
c           --- no side chain improper torsions
c           --- set up the charges (depends on termflag)
            if ( termflag .eq. -1 ) then
c              --- n-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,0.09660d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.00150d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.72140d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.60130d0)
c              --- n-terminus hydrogens
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.21650d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.21650d0)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.21650d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.11800d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,0.02120d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,-0.00480d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.06080d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.01810d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.37640d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.02830d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.02830d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.01210d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.01210d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.06330d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.06330d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.11710d0)
               call twh_qqatom(GLB_SET,imolty,shift+7,0.11710d0)
               call twh_qqatom(GLB_SET,imolty,shift+8,0.33820d0)
               call twh_qqatom(GLB_SET,imolty,shift+9,0.33820d0)
               call twh_qqatom(GLB_SET,imolty,shift+10,0.33820d0)
            elseif ( termflag .eq. 0 ) then
c              --- main chain
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.34790d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.24000d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.73410d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.58940d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.27470d0)
c              --- Ca hydrogen (Ha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.14260d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.00940d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.01870d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.04790d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.01430d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.38540d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.03620d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.03620d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.01030d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.01030d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.06210d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.06210d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.11350d0)
               call twh_qqatom(GLB_SET,imolty,shift+7,0.11350d0)
               call twh_qqatom(GLB_SET,imolty,shift+8,0.34000d0)
               call twh_qqatom(GLB_SET,imolty,shift+9,0.34000d0)
               call twh_qqatom(GLB_SET,imolty,shift+10,0.34000d0)
            elseif ( termflag .eq. 1 ) then
c              --- c-terminus
c              --- N - Ca - C = O
               call twh_qqatom(GLB_SET,imolty,iatom,-0.34810d0)
               call twh_qqatom(GLB_SET,imolty,iatom+1,-0.29030d0)
               call twh_qqatom(GLB_SET,imolty,iatom+2,0.84880d0)
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.82520d0)
c              --- extra oxygen on c-terminus
               call twh_qqatom(GLB_SET,imolty,hmain+heavy,-0.82520d0)
c              --- N hydrogen
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.27640d0)
c              --- Ca hydrogen (Halpha)
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.14380d0)
c              --- side chain
               call twh_qqatom(GLB_SET,imolty,hmain,-0.05380d0)
               call twh_qqatom(GLB_SET,imolty,hmain+1,0.02270d0)
               call twh_qqatom(GLB_SET,imolty,hmain+2,-0.03920d0)
               call twh_qqatom(GLB_SET,imolty,hmain+3,-0.01760d0)
               call twh_qqatom(GLB_SET,imolty,hmain+4,-0.37410d0)
               call twh_qqatom(GLB_SET,imolty,shift,0.04820d0)
               call twh_qqatom(GLB_SET,imolty,shift+1,0.04820d0)
               call twh_qqatom(GLB_SET,imolty,shift+2,0.01340d0)
               call twh_qqatom(GLB_SET,imolty,shift+3,0.01340d0)
               call twh_qqatom(GLB_SET,imolty,shift+4,0.06110d0)
               call twh_qqatom(GLB_SET,imolty,shift+5,0.06110d0)
               call twh_qqatom(GLB_SET,imolty,shift+6,0.11210d0)
               call twh_qqatom(GLB_SET,imolty,shift+7,0.11210d0)
               call twh_qqatom(GLB_SET,imolty,shift+8,0.33740d0)
               call twh_qqatom(GLB_SET,imolty,shift+9,0.33740d0)
               call twh_qqatom(GLB_SET,imolty,shift+10,0.33740d0)
            endif
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
            atomname(iatom) =   'N3'
            atomname(iatom+1) = 'CT'
            atomname(iatom+2) = 'C'
            atomname(iatom+3) = 'O'
c           --- N-terminus hydrogens
            atomname(hmain+heavy) = 'H'
            atomname(hmain+heavy+1) = 'H'
            atomname(hmain+heavy+2) = 'H'
c           --- C alpha hydrogen
            atomname(hmain+heavy+3) = 'HP'
c           --- improper torsions
c           --- backbone C-Ca-N-O (in amber Ca +N C O)
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+3)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,iatom+1)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,iatom+total)
         elseif ( termflag .eq. 0 ) then
            atomname(iatom) =   'N'
            atomname(iatom+1) = 'CT'
            atomname(iatom+2) = 'C'
            atomname(iatom+3) = 'O'
c           --- N hydrogen
            atomname(hmain+heavy) = 'H'
c           --- Calpha hydrogen
            atomname(hmain+heavy+1) = 'H1'
c           --- improper torsions
c           --- backbone N-C-Ca-Hn (in amber -C Ca N Hn)
            call twh_inimprop(GLB_SET,imolty,iatom,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom,1,hmain+heavy)
            call twh_ijimprop3(GLB_SET,imolty,iatom,1,ibondend)
            call twh_ijimprop4(GLB_SET,imolty,iatom,1,iatom+1)
c           --- backbone C-Ca-N-O (in amber Ca +N C O)
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+3)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,iatom+1)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,iatom+total)
         elseif ( termflag .eq. 1 ) then
c           --- c-terminus main chain
c           --- N - Ca - C = O
            atomname(iatom) =   'N'
            atomname(iatom+1) = 'CT'
            atomname(iatom+2) = 'C'
            atomname(iatom+3) = 'O2'
c            --- c-terminus oxygen
            atomname(hmain+heavy) = 'O2'
c           --- N hydrogen
            atomname(hmain+heavy+1) = 'H'
c           --- Calpha hydrogen
            atomname(hmain+heavy+2) = 'H1'
c           --- improper torsions
c           --- backbone N-C-Ca-Hn (in amber -C Ca N Hn)
            call twh_inimprop(GLB_SET,imolty,iatom,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom,1,hmain+heavy+1)
            call twh_ijimprop3(GLB_SET,imolty,iatom,1,ibondend)
            call twh_ijimprop4(GLB_SET,imolty,iatom,1,iatom+1)
c           --- C-terminus C-Ca-O-O (in amber Ca Ox C O)
            call twh_inimprop(GLB_SET,imolty,iatom+2,1)
            call twh_ijimprop2(GLB_SET,imolty,iatom+2,1,iatom+3)
            call twh_ijimprop3(GLB_SET,imolty,iatom+2,1,iatom+1)
            call twh_ijimprop4(GLB_SET,imolty,iatom+2,1,hmain+heavy)
         else
c           --- unknown termflag
            write(6,*) 'POLYAMBER96 unknown termflag:',termflag
            write(6,*) 'imolty,ipeptide,pepname',imolty,ipeptide
     &           ,pepname(ipeptide)
            lfailure = .true.
            return
         endif
c        --- update ibondend for the typical amino acids
         ibondend = iatom+2
      endif
c      write(6,*) 'finished POLYAMBER96'
      return
      end

#else
c     Compilation disabled (for speed) with FFLAGS=-DDISABLE_POLY
c     statement in ./configure
      write(6,*) "polyamber96 disabled (compile flag)"
      lfailure = .true.
      return

      end
#endif
