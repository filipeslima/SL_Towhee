#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2003-2011 Marcus G. Martin                       *
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
      subroutine twh_polycharmm27na(lfailure,imolty,iatom,termflag
     &     ,total,ibondend,inucleic,pepname,atomname,aminoname
     &     ,shortname,hterm)
c     ******************************************************************
c     * builds the input file for a nucleic acid sequence starting from*
c     * the N-terminus for the Charmm27 force field                    *
c     *                                                                *
c     * originally written 04-10-2003 by C.D. Lorenz                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
#ifndef DISABLE_POLY
      implicit none
c     ---- variables passed to/from the subroutine
      logical lfailure
      integer imolty,iatom,total,ibondend,inucleic
      integer hterm
      character*2 pepname
      dimension pepname(NUMAX)
      character*3 shortname
      character*10 atomname
      character*20 termflag
      dimension atomname(NUMAX)
      character*40 aminoname

c     --- local variables
      integer shift,heavy,ninsert,hmain,deoxy,ilinker,light

c     write(6,*) 'start POLYCHARMM27NA'
      if (pepname(inucleic) .eq. 'A' .or.
     &     pepname(inucleic) .eq. 'a' .or.           
     &     pepname(inucleic) .eq. 'G' .or.           
     &     pepname(inucleic) .eq. 'g' .or.           
     &     pepname(inucleic) .eq. 'U' .or.           
     &     pepname(inucleic) .eq. 'u' .or.           
     &     pepname(inucleic) .eq. 'C' .or.           
     &     pepname(inucleic) .eq. 'c') then
         deoxy = 0    
      else
         deoxy = 1
      endif 
      if (termflag .eq. '5 prime' ) then
c        --- 5' terminal end
         if (hterm .eq. 5) then
            if (deoxy .eq. 0) then
               hmain = iatom + 9 
            elseif (deoxy .eq. 1) then
               hmain = iatom + 8
            else
               write(6,*) 'ERROR: unknown type of nucleic acid'
               lfailure = .true.
               return
            endif
            ilinker = iatom+4
         elseif (hterm .eq. 3) then
            if (deoxy .eq. 0) then
               hmain = iatom + 13
            elseif (deoxy .eq. 1) then
               hmain = iatom + 12
            else
               write(6,*) "ERROR: unknown type of nucleic acid"
               lfailure = .true.
               return
            endif
            ilinker = iatom+8
         endif
         ninsert = 8
      elseif (termflag .eq. 'neither' ) then
         if (deoxy .eq. 0) then
            hmain = iatom + 12
         elseif (deoxy .eq. 1) then
            hmain = iatom + 11
         else
            write(6,*) "ERROR:unknown type of nucleic acid"
            lfailure = .true.
            return
         endif
         ninsert = 7
         ilinker = iatom+7
      elseif (termflag .eq. '3 prime' ) then
         if (hterm .eq. 3) then
c           --- 3' terminal end
            if (deoxy .eq. 0) then
               hmain = iatom + 12
            elseif (deoxy .eq. 1) then
               hmain = iatom + 11
            else
               write(6,*) "ERROR:unknown type of nucleic acid"
               lfailure = .true.
               return
            endif  
         elseif (hterm .eq. 5) then
            if (deoxy .eq. 0) then
               hmain = iatom + 16
            elseif (deoxy .eq. 1) then
               hmain = iatom + 15
            else
               write(6,*) "ERROR:unknown type of nucleic acid"
               lfailure = .true.
               return
            endif  
         endif
         ninsert = 8
         ilinker = iatom+7
      else
         write(6,*) 'POLYCHARMM27na: unknown termflag'
         write(6,*) termflag
         lfailure = .true.
         return
      endif   
c     --- determine which nucleic acid is built this step
      if ( pepname(inucleic) .eq. 'G' .or. 
     &     pepname(inucleic) .eq. 'g' .or.
     &     pepname(inucleic) .eq. 'dG' .or.
     &     pepname(inucleic) .eq. 'dg') then
c        --- guanine base
         if (pepname(inucleic) .eq. 'G' .or.
     &        pepname(inucleic) .eq. 'g') then
            aminoname = 'guanosine'
            shortname = 'GUA'
         elseif (pepname(inucleic) .eq. 'dG' .or.
     &           pepname(inucleic) .eq. 'dg') then
            aminoname = 'deoxyguanosine'
            shortname = 'dGU'
         endif
         call twh_aa_nucleic_g(imolty,heavy,light,shift,hmain,ilinker
     &        ,ninsert)
         atomname(hmain) = 'NN2B '
         call twh_qqatom(GLB_SET,imolty,hmain,-0.02d0)
         atomname(hmain+1) = 'CN4 '
         call twh_qqatom(GLB_SET,imolty,hmain+1,0.25d0)
         atomname(hmain+2) = 'NN4 '
         call twh_qqatom(GLB_SET,imolty,hmain+2,-0.60d0)
         atomname(hmain+3) = 'CN5G '
         call twh_qqatom(GLB_SET,imolty,hmain+3,0.00d0)
         atomname(hmain+4) = 'CN1 '
         call twh_qqatom(GLB_SET,imolty,hmain+4,0.54d0)
         atomname(hmain+5) = 'ON1 '
         call twh_qqatom(GLB_SET,imolty,hmain+5,-0.51d0)
         atomname(hmain+6) = 'NN2G '
         call twh_qqatom(GLB_SET,imolty,hmain+6,-0.34d0)
         atomname(hmain+7) = 'CN2 '
         call twh_qqatom(GLB_SET,imolty,hmain+7,0.75d0)
         atomname(hmain+8) = 'NN1 '
         call twh_qqatom(GLB_SET,imolty,hmain+8,-0.68d0)
         atomname(hmain+9) = 'NN3G '
         call twh_qqatom(GLB_SET,imolty,hmain+9,-0.74d0)
         atomname(hmain+10) = 'CN5 '
         call twh_qqatom(GLB_SET,imolty,hmain+10,0.26d0)
         atomname(shift) = 'HN3  '
         call twh_qqatom(GLB_SET,imolty,shift,0.16d0)
         atomname(shift+1) = 'HN2  '
         call twh_qqatom(GLB_SET,imolty,shift+1,0.26d0)
         atomname(shift+2) = 'HN1 '
         call twh_qqatom(GLB_SET,imolty,shift+2,0.32d0)
         atomname(shift+3) = 'HN1 '
         call twh_qqatom(GLB_SET,imolty,shift+3,0.35d0)
c        --- side chain improper torsions
         call twh_inimprop(GLB_SET,imolty,hmain+4,1)
         call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,hmain+6)
         call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+3)
         call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+5)
         call twh_inimprop(GLB_SET,imolty,hmain+7,1)
         call twh_ijimprop2(GLB_SET,imolty,hmain+7,1,hmain+9)
         call twh_ijimprop3(GLB_SET,imolty,hmain+7,1,hmain+6)
         call twh_ijimprop4(GLB_SET,imolty,hmain+7,1,hmain+8)
         call twh_inimprop(GLB_SET,imolty,hmain+8,1)
         call twh_ijimprop2(GLB_SET,imolty,hmain+8,1,shift+2)
         call twh_ijimprop3(GLB_SET,imolty,hmain+8,1,hmain+7)
         call twh_ijimprop4(GLB_SET,imolty,hmain+8,1,shift+3)
      elseif ( pepname(inucleic) .eq. 'A' .or.
     &        pepname(inucleic) .eq. 'a' .or.
     &        pepname(inucleic) .eq. 'dA' .or.
     &        pepname(inucleic) .eq. 'da' ) then
c        --- adenine base
         if (pepname(inucleic) .eq. 'A' .or.
     &         pepname(inucleic) .eq. 'a') then
             aminoname = 'adenosine'
             shortname = 'AMP'
         elseif (pepname(inucleic) .eq. 'dA' .or.
     &             pepname(inucleic) .eq. 'da') then
             aminoname = 'deoxyadenosine'
             shortname = 'dAM'
         endif
         call twh_aa_nucleic_a(imolty,heavy,light,shift,hmain,ilinker
     &        ,ninsert)
         atomname(hmain) = 'NN2 '
         call twh_qqatom(GLB_SET,imolty,hmain,-0.05d0)
         atomname(hmain+1) = 'CN4 '
         call twh_qqatom(GLB_SET,imolty,hmain+1,0.34d0)
         atomname(hmain+2) = 'NN4 '
         call twh_qqatom(GLB_SET,imolty,hmain+2,-0.71d0)
         atomname(hmain+3) = 'CN5 '
         call twh_qqatom(GLB_SET,imolty,hmain+3,0.28d0)
         atomname(hmain+4) = 'CN2 '
         call twh_qqatom(GLB_SET,imolty,hmain+4,0.46d0)
         atomname(hmain+5) = 'NN1 '
         call twh_qqatom(GLB_SET,imolty,hmain+5,-0.77d0)
         atomname(hmain+6) = 'NN3A '
         call twh_qqatom(GLB_SET,imolty,hmain+6,-0.74d0)
         atomname(hmain+7) = 'CN4 '
         call twh_qqatom(GLB_SET,imolty,hmain+7,0.50d0)
         atomname(hmain+8) = 'NN3A '
         call twh_qqatom(GLB_SET,imolty,hmain+8,-0.75d0)
         atomname(hmain+9) = 'CN5 '
         call twh_qqatom(GLB_SET,imolty,hmain+9,0.43d0)
         atomname(shift) = 'HN3  '
         call twh_qqatom(GLB_SET,imolty,shift,0.12d0)
         atomname(shift+1) = 'HN1  '
         call twh_qqatom(GLB_SET,imolty,shift+1,0.38d0)
         atomname(shift+2) = 'HN1 '
         call twh_qqatom(GLB_SET,imolty,shift+2,0.38d0)
         atomname(shift+3) = 'HN3 '
         call twh_qqatom(GLB_SET,imolty,shift+3,0.13d0)
c        --- side chain improper torsions
         call twh_inimprop(GLB_SET,imolty,hmain+4,1)
         call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,hmain+6)
         call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+3)
         call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+5)
         call twh_inimprop(GLB_SET,imolty,hmain+5,1)
         call twh_ijimprop2(GLB_SET,imolty,hmain+5,1,hmain+4)
         call twh_ijimprop3(GLB_SET,imolty,hmain+5,1,shift+1)
         call twh_ijimprop4(GLB_SET,imolty,hmain+5,1,shift+2)
      elseif ( pepname(inucleic) .eq. 'C' .or.
     &        pepname(inucleic) .eq. 'c' .or.
     &        pepname(inucleic) .eq. 'dC' .or.
     &        pepname(inucleic) .eq. 'dc') then
c        --- cytosine base
         if (pepname(inucleic) .eq. 'C' .or.
     &        pepname(inucleic) .eq. 'c') then
            aminoname = 'cytidine'
            shortname = 'CMP'
         elseif (pepname(inucleic) .eq. 'dC' .or.
     &           pepname(inucleic) .eq. 'dc') then
            aminoname = 'deoxycytidine'
            shortname = 'dCM'
         endif    
         call twh_aa_nucleic_c(imolty,heavy,light,shift,hmain,ilinker
     &        ,ninsert)
         atomname(hmain) = 'NN2 '
         call twh_qqatom(GLB_SET,imolty,hmain,-0.13d0)
         atomname(hmain+1) = 'CN1 '
         call twh_qqatom(GLB_SET,imolty,hmain+1,0.52d0)
         atomname(hmain+2) = 'ON1C '
         call twh_qqatom(GLB_SET,imolty,hmain+2,-0.49d0)
         atomname(hmain+3) = 'NN3 '
         call twh_qqatom(GLB_SET,imolty,hmain+3,-0.66d0)
         atomname(hmain+4) = 'CN2 '
         call twh_qqatom(GLB_SET,imolty,hmain+4,0.65d0)
         atomname(hmain+5) = 'NN1 '
         call twh_qqatom(GLB_SET,imolty,hmain+5,-0.75d0)
         atomname(hmain+6) = 'CN3 '
         call twh_qqatom(GLB_SET,imolty,hmain+6,-0.13d0)
         atomname(hmain+7) = 'CN3 '
         call twh_qqatom(GLB_SET,imolty,hmain+7,0.05d0)
         atomname(shift) = 'HN1  '
         call twh_qqatom(GLB_SET,imolty,shift,0.37d0)
         atomname(shift+1) = 'HN1  '
         call twh_qqatom(GLB_SET,imolty,shift+1,0.33d0)
         atomname(shift+2) = 'HN3 '
         call twh_qqatom(GLB_SET,imolty,shift+2,0.07d0)
         atomname(shift+3) = 'HN3 '
         call twh_qqatom(GLB_SET,imolty,shift+3,0.17d0)
c        --- side chain improper torsions
         call twh_inimprop(GLB_SET,imolty,hmain+1,1)
         call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain)
         call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+3)
         call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+2)
         call twh_inimprop(GLB_SET,imolty,hmain+4,1)
         call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,hmain+3)
         call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+6)
         call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+5)
         call twh_inimprop(GLB_SET,imolty,hmain+5,1)
         call twh_ijimprop2(GLB_SET,imolty,hmain+5,1,hmain+4)
         call twh_ijimprop3(GLB_SET,imolty,hmain+5,1,shift)
         call twh_ijimprop4(GLB_SET,imolty,hmain+5,1,shift+1)
      elseif ( pepname(inucleic) .eq. 'dT' .or.
     &        pepname(inucleic) .eq. 'dt' .or. 
     &        pepname(inucleic) .eq. 'T' .or.
     &        pepname(inucleic) .eq. 't') then
c        --- thymine base - only present in DNA
         if (pepname(inucleic) .eq. 'T' .or.
     &        pepname(inucleic) .eq. 't') then
            aminoname = 'thymidine'
            shortname = 'TMP'
         elseif (pepname(inucleic) .eq. 'dT' .or.
     &           pepname(inucleic) .eq. 'dt') then
            aminoname = 'deoxythymidine'
            shortname = 'dTM'
         endif    
         call twh_aa_nucleic_t(imolty,heavy,light,shift,hmain,ilinker
     &        ,ninsert)
         atomname(hmain) = 'NN2B '
         call twh_qqatom(GLB_SET,imolty,hmain,-0.34d0)
         atomname(hmain+1) = 'CN1T '
         call twh_qqatom(GLB_SET,imolty,hmain+1,0.51d0)
         atomname(hmain+2) = 'ON1 '
         call twh_qqatom(GLB_SET,imolty,hmain+2,-0.41d0)
         atomname(hmain+3) = 'NN2U '
         call twh_qqatom(GLB_SET,imolty,hmain+3,-0.46d0)
         atomname(hmain+4) = 'CN1 '
         call twh_qqatom(GLB_SET,imolty,hmain+4,0.50d0)
         atomname(hmain+5) = 'ON1 '
         call twh_qqatom(GLB_SET,imolty,hmain+5,-0.45d0)
         atomname(hmain+6) = 'CN3T '
         call twh_qqatom(GLB_SET,imolty,hmain+6,-0.15d0)
         atomname(hmain+7) = 'CN9 '
         call twh_qqatom(GLB_SET,imolty,hmain+7,-0.11d0)
         atomname(hmain+8) = 'CN3 '
         call twh_qqatom(GLB_SET,imolty,hmain+8,0.17d0)
         atomname(shift) = 'HN2  '
         call twh_qqatom(GLB_SET,imolty,shift,0.36d0)
         atomname(shift+1) = 'HN9  '
         call twh_qqatom(GLB_SET,imolty,shift+1,0.07d0)
         atomname(shift+2) = 'HN9 '
         call twh_qqatom(GLB_SET,imolty,shift+2,0.07d0)
         atomname(shift+3) = 'HN9 '
         call twh_qqatom(GLB_SET,imolty,shift+3,0.07d0)
         atomname(shift+4) = 'HN3 '
         call twh_qqatom(GLB_SET,imolty,shift+4,0.17d0)
c        --- side chain improper torsions
         call twh_inimprop(GLB_SET,imolty,hmain+1,1)
         call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain)
         call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+3)
         call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+2)
         call twh_inimprop(GLB_SET,imolty,hmain+4,1)
         call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,hmain+3)
         call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+6)
         call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+5)
         call twh_inimprop(GLB_SET,imolty,hmain+6,1)
         call twh_ijimprop2(GLB_SET,imolty,hmain+6,1,hmain+4)
         call twh_ijimprop3(GLB_SET,imolty,hmain+6,1,hmain+8)
         call twh_ijimprop4(GLB_SET,imolty,hmain+6,1,hmain+7)
      elseif ( pepname(inucleic) .eq. 'U' .or.
     &        pepname(inucleic) .eq. 'u' .or. 
     &        pepname(inucleic) .eq. 'dU' .or.
     &        pepname(inucleic) .eq. 'du' ) then
c        --- uracil base - only present in RNA
c        --- however, the possibility to use it as a 
c        --- deoxyribonucleic acid has been included
c        --- just in case there is any reason why someone
c        --- might want to use it
         if (pepname(inucleic) .eq. 'U' .or.
     &        pepname(inucleic) .eq. 'u') then
            aminoname = 'uridine'
            shortname = 'UMP'
         elseif (pepname(inucleic) .eq. 'dU' .or.
     &           pepname(inucleic) .eq. 'du') then
            aminoname = 'deoxyuridine'
            shortname = 'dUM'
         endif    
         call twh_aa_nucleic_u(imolty,heavy,light,shift,hmain,ilinker
     &        ,ninsert)
         atomname(hmain) = 'NN2B '
         call twh_qqatom(GLB_SET,imolty,hmain,-0.34d0)
         atomname(hmain+1) = 'CN1T '
         call twh_qqatom(GLB_SET,imolty,hmain+1,0.55d0)
         atomname(hmain+2) = 'ON1 '
         call twh_qqatom(GLB_SET,imolty,hmain+2,-0.45d0)
         atomname(hmain+3) = 'NN2U '
         call twh_qqatom(GLB_SET,imolty,hmain+3,-0.46d0)
         atomname(hmain+4) = 'CN1 '
         call twh_qqatom(GLB_SET,imolty,hmain+4,0.53d0)
         atomname(hmain+5) = 'ON1 '
         call twh_qqatom(GLB_SET,imolty,hmain+5,-0.48d0)
         atomname(hmain+6) = 'CN3 '
         call twh_qqatom(GLB_SET,imolty,hmain+6,-0.15d0)
         atomname(hmain+7) = 'CN3 '
         call twh_qqatom(GLB_SET,imolty,hmain+7,0.20d0)
         atomname(shift) = 'HN2  '
         call twh_qqatom(GLB_SET,imolty,shift,0.36d0)
         atomname(shift+1) = 'HN3  '
         call twh_qqatom(GLB_SET,imolty,shift+1,0.10d0)
         atomname(shift+2) = 'HN3 '
         call twh_qqatom(GLB_SET,imolty,shift+2,0.14d0)
c        --- side chain improper torsions
         call twh_inimprop(GLB_SET,imolty,hmain+1,1)
         call twh_ijimprop2(GLB_SET,imolty,hmain+1,1,hmain)
         call twh_ijimprop3(GLB_SET,imolty,hmain+1,1,hmain+3)
         call twh_ijimprop4(GLB_SET,imolty,hmain+1,1,hmain+2)
         call twh_inimprop(GLB_SET,imolty,hmain+4,1)
         call twh_ijimprop2(GLB_SET,imolty,hmain+4,1,hmain+3)
         call twh_ijimprop3(GLB_SET,imolty,hmain+4,1,hmain+6)
         call twh_ijimprop4(GLB_SET,imolty,hmain+4,1,hmain+5)
c        --- no side chain improper torsions
      elseif (pepname(inucleic) .eq. 'T' .or.
     &        pepname(inucleic) .eq. 't') then
         write(6,*) 'thymine is generally found in the form'
         write(6,*) 'of a deoxyribonucleotide'
         write(6,*) 'Double check that if simulating RNA'
         write(6,*) 'you want to use uracil base (U)' 
         lfailure = .true.
         return
      elseif (pepname(inucleic) .eq. 'dU' .or.
     &        pepname(inucleic) .eq. 'du') then
         write(6,*) 'uracil is generally found in the form'
         write(6,*) 'of a ribonucleotide'
         write(6,*) 'Double check that if simulating DNA'
         write(6,*) 'you want to use thymidine base (dT)'
         lfailure = .true.
         return
      else
         write(6,*) 'invalid nucleotide name ',inucleic
     &        ,' ',pepname(inucleic)
         lfailure = .true.
         return
      endif

      call twh_aa_nucleic_main_chain(lfailure,termflag,hterm,total,heavy
     &     ,ninsert,hmain,light,iatom,shift,deoxy,imolty,ibondend)
      if ( lfailure ) return

      if ( termflag .eq. '5 prime' ) then
         if (hterm .eq. 5) then
c           --- 5'-terminus main chain
            atomname(iatom) = 'ON5'
            call twh_qqatom(GLB_SET,imolty,iatom,-0.66d0)
            atomname(iatom+1) = 'CN8B'
            call twh_qqatom(GLB_SET,imolty,iatom+1,0.05d0)
            atomname(iatom+2) = 'CN7'
            call twh_qqatom(GLB_SET,imolty,iatom+2,0.16d0)
            if (deoxy .eq. 0) then
               atomname(iatom+3) = 'ON6B'
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.50d0)
            else if (deoxy .eq. 1) then
               atomname(iatom+3) = 'ON6'
               call twh_qqatom(GLB_SET,imolty,iatom+3,-0.50d0)
            end if
            atomname(iatom+4) = 'CN7B'
            call twh_qqatom(GLB_SET,imolty,iatom+4,0.16d0)
            if (deoxy .eq. 0) then
               atomname(iatom+5) = 'CN7B'
               call twh_qqatom(GLB_SET,imolty,iatom+5,0.14d0)
            elseif (deoxy .eq. 1) then
               atomname(iatom+5) = 'CN8'
               call twh_qqatom(GLB_SET,imolty,iatom+5,-0.18d0)
            else
               write(6,*) "Illegal value of deoxy",deoxy
               lfailure = .true.
               return
            endif
            atomname(iatom+6) = 'CN7'
            call twh_qqatom(GLB_SET,imolty,iatom+6,0.01d0)
            atomname(iatom+7) =   'ON2'
            call twh_qqatom(GLB_SET,imolty,iatom+7,-0.57d0)
            if (deoxy .eq. 0) then
               atomname(iatom+8) = 'ON5'
               call twh_qqatom(GLB_SET,imolty,iatom+8,-0.66d0)
            endif
            atomname(hmain+heavy) = 'HN5'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.43d0)
            atomname(hmain+heavy+1) = 'HN8'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.09d0)
            atomname(hmain+heavy+2) = 'HN8'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.09d0)
            atomname(hmain+heavy+3) = 'HN7'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.09d0)
            atomname(hmain+heavy+4) = 'HN7'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+4,0.09d0)
            if (deoxy .eq. 0) then
               atomname(hmain+heavy+6) = 'HN5'
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+6,0.43d0)
               atomname(hmain+heavy+5) = 'HN7'
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+5,0.09d0)
            elseif (deoxy .eq. 1) then
               atomname(hmain+heavy+6) = 'HN8'
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+6,0.09d0)
               atomname(hmain+heavy+5) = 'HN8'
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+5,0.09d0)
            endif
            atomname(hmain+heavy+7) = 'HN7'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+7,0.09d0)
         elseif (hterm .eq. 3) then 
c           --- 3'-terminus main chain
            atomname(iatom) = 'P'
            call twh_qqatom(GLB_SET,imolty,iatom,1.50d0)
            atomname(iatom+1) = 'ON4'
            call twh_qqatom(GLB_SET,imolty,iatom+1,-0.68d0)
            atomname(iatom+2) = 'ON3'
            call twh_qqatom(GLB_SET,imolty,iatom+2,-0.82d0)
            atomname(iatom+3) = 'ON3'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.82d0)
            atomname(iatom+4) = 'ON2'
            call twh_qqatom(GLB_SET,imolty,iatom+4,-0.62d0)
            atomname(iatom+5) = 'CN8B'
            call twh_qqatom(GLB_SET,imolty,iatom+5,-0.08d0)
            atomname(iatom+6) = 'CN7'
            call twh_qqatom(GLB_SET,imolty,iatom+6,0.16d0)
            if (deoxy .eq. 0) then
               atomname(iatom+7) = 'ON6B'
               call twh_qqatom(GLB_SET,imolty,iatom+7,-0.50d0)
            else if (deoxy .eq. 1) then
               atomname(iatom+7) = 'ON6'
               call twh_qqatom(GLB_SET,imolty,iatom+7,-0.50d0)
            end if
            atomname(iatom+8) = 'CN7B'
            call twh_qqatom(GLB_SET,imolty,iatom+8,0.16d0)
            if (deoxy .eq. 0) then
               atomname(iatom+9) = 'CN7B'
               call twh_qqatom(GLB_SET,imolty,iatom+9,0.14d0)
            elseif (deoxy .eq. 1) then
               atomname(iatom+9) = 'CN8'
               call twh_qqatom(GLB_SET,imolty,iatom+9,-0.18d0)
            else
               write(6,*) "Illegal value of deoxy",deoxy
               lfailure = .true.
               return
            endif
            atomname(iatom+10) = 'CN7'
            call twh_qqatom(GLB_SET,imolty,iatom+10,0.01d0)
            atomname(iatom+11) =   'ON2'
            call twh_qqatom(GLB_SET,imolty,iatom+11,-0.57d0)
            if (deoxy .eq. 0) then
               atomname(iatom+12) = 'ON5'
               call twh_qqatom(GLB_SET,imolty,iatom+12,-0.66d0)
            endif
            atomname(hmain+heavy) = 'HN4'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.34d0)
            atomname(hmain+heavy+1) = 'HN8'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.09d0)
            atomname(hmain+heavy+2) = 'HN8'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.09d0)
            atomname(hmain+heavy+3) = 'HN7'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.09d0)
            atomname(hmain+heavy+4) = 'HN7'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+4,0.09d0)
            if (deoxy .eq. 0) then
               atomname(hmain+heavy+6) = 'HN5'
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+6,0.43d0)
               atomname(hmain+heavy+5) = 'HN7'
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+5,0.09d0)
            elseif (deoxy .eq. 1) then
               atomname(hmain+heavy+6) = 'HN8'
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+6,0.09d0)
               atomname(hmain+heavy+5) = 'HN8'
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+5,0.09d0)
            endif
            atomname(hmain+heavy+7) = 'HN7'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+7,0.09d0)
         endif
      elseif ( termflag .eq. 'neither' ) then
c        --- interior segment main chain
         atomname(iatom) =   'P'
         call twh_qqatom(GLB_SET,imolty,iatom,1.50d0)
         atomname(iatom+1) = 'ON3'
         call twh_qqatom(GLB_SET,imolty,iatom+1,-0.78d0)
         atomname(iatom+2) = 'ON3'
         call twh_qqatom(GLB_SET,imolty,iatom+2,-0.78d0)
         atomname(iatom+3) = 'ON2'
         call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57d0)
         atomname(iatom+4) = 'CN8B'
         call twh_qqatom(GLB_SET,imolty,iatom+4,-0.08d0)
         atomname(iatom+5) = 'CN7'
         call twh_qqatom(GLB_SET,imolty,iatom+5,0.16d0)
         if (deoxy .eq. 0) then
            atomname(iatom+6) = 'ON6B'
            call twh_qqatom(GLB_SET,imolty,iatom+6,-0.50d0)
         else if (deoxy .eq. 1) then
            atomname(iatom+6) = 'ON6'
            call twh_qqatom(GLB_SET,imolty,iatom+6,-0.50d0)
         end if     
         atomname(iatom+7) = 'CN7B'
         call twh_qqatom(GLB_SET,imolty,iatom+7,0.16d0)
         if (deoxy .eq. 0) then
             atomname(iatom+8) = 'CN7B'
             call twh_qqatom(GLB_SET,imolty,iatom+8,0.14d0)
         elseif (deoxy .eq. 1) then
             atomname(iatom+8) = 'CN8'
             call twh_qqatom(GLB_SET,imolty,iatom+8,-0.18d0)
         endif
         atomname(iatom+9) = 'CN7'
         call twh_qqatom(GLB_SET,imolty,iatom+9,0.01d0)
         atomname(iatom+10) = 'ON2'
         call twh_qqatom(GLB_SET,imolty,iatom+10,-0.57d0)
         if (deoxy .eq. 0) then
             atomname(iatom+11) = 'ON5'
             call twh_qqatom(GLB_SET,imolty,iatom+11,-0.66d0)
         endif
         atomname(hmain+heavy) = 'HN8'
         call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.09d0)
         atomname(hmain+heavy+1) = 'HN8'
         call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.09d0)
         atomname(hmain+heavy+2) = 'HN7'
         call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.09d0)
         atomname(hmain+heavy+3) = 'HN7'
         call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.09d0)
         if (deoxy .eq. 0) then         
            atomname(hmain+heavy+4) = 'HN7'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+4,0.09d0)
         elseif (deoxy .eq. 1) then
            atomname(hmain+heavy+4) = 'HN8'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+4,0.09d0)
         endif
         if (deoxy .eq. 0) then
            atomname(hmain+heavy+5) = 'HN5'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+5,0.43d0)
         elseif (deoxy .eq. 1) then
            atomname(hmain+heavy+5) = 'HN8'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+5,0.09d0)
         endif
         atomname(hmain+heavy+6) = 'HN7'
         call twh_qqatom(GLB_SET,imolty,hmain+heavy+6,0.09d0)
      elseif ( termflag .eq. '3 prime' ) then
         if (hterm .eq. 3) then
c           --- 3'-terminus main chain
            atomname(iatom) = 'P'
            call twh_qqatom(GLB_SET,imolty,iatom,1.50d0)
            atomname(iatom+1) = 'ON3'
            call twh_qqatom(GLB_SET,imolty,iatom+1,-0.78d0)
            atomname(iatom+2) = 'ON3'
            call twh_qqatom(GLB_SET,imolty,iatom+2,-0.78d0)
            atomname(iatom+3) = 'ON2'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57d0)
            atomname(iatom+4) = 'CN8B'
            call twh_qqatom(GLB_SET,imolty,iatom+4,-0.08d0)
            atomname(iatom+5) = 'CN7'
            call twh_qqatom(GLB_SET,imolty,iatom+5,0.16d0)
            if (deoxy .eq. 0) then
               atomname(iatom+6) = 'ON6B'
               call twh_qqatom(GLB_SET,imolty,iatom+6,-0.50d0)
            else if (deoxy .eq. 1) then
               atomname(iatom+6) = 'ON6'
               call twh_qqatom(GLB_SET,imolty,iatom+6,-0.50d0)
            endif            
            atomname(iatom+7) = 'CN7B'
            call twh_qqatom(GLB_SET,imolty,iatom+7,0.16d0)
            if (deoxy .eq. 0) then
               atomname(iatom+8) = 'CN7B'
               call twh_qqatom(GLB_SET,imolty,iatom+8,0.14d0)
            elseif (deoxy .eq. 1) then
               atomname(iatom+8) = 'CN8'
               call twh_qqatom(GLB_SET,imolty,iatom+8,-0.18d0)
            endif
            atomname(iatom+9) = 'CN7'
            call twh_qqatom(GLB_SET,imolty,iatom+9,0.14d0)
            atomname(iatom+10) =   'ON5'
            call twh_qqatom(GLB_SET,imolty,iatom+10,-0.66d0)
            if (deoxy .eq. 0) then
               atomname(iatom+11) = 'ON5'
               call twh_qqatom(GLB_SET,imolty,iatom+11,-0.66d0)
            endif
            atomname(hmain+heavy) = 'HN8'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.09d0)
            atomname(hmain+heavy+1) = 'HN8'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.09d0)
            atomname(hmain+heavy+2) = 'HN7'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.09d0)
            atomname(hmain+heavy+3) = 'HN7'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.09d0)
            if (deoxy .eq. 0) then
               atomname(hmain+heavy+5) = 'HN5'
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+5,0.43d0)
               atomname(hmain+heavy+4) = 'HN7'
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+4,0.09d0)
            elseif (deoxy .eq. 1) then
               atomname(hmain+heavy+5) = 'HN8'
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+5,0.09d0)
               atomname(hmain+heavy+4) = 'HN8'
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+4,0.09d0)
            endif
            atomname(hmain+heavy+6) = 'HN7'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+6,0.09d0)
            atomname(hmain+heavy+7) = 'HN5'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+7,0.43d0)
         elseif (hterm .eq. 5) then
c           --- 5'-terminus main chain
            atomname(iatom) = 'P'
            call twh_qqatom(GLB_SET,imolty,iatom,1.50d0)
            atomname(iatom+1) = 'ON3'
            call twh_qqatom(GLB_SET,imolty,iatom+1,-0.78d0)
            atomname(iatom+2) = 'ON3'
            call twh_qqatom(GLB_SET,imolty,iatom+2,-0.78d0)
            atomname(iatom+3) = 'ON2'
            call twh_qqatom(GLB_SET,imolty,iatom+3,-0.57d0)
            atomname(iatom+4) = 'CN8B'
            call twh_qqatom(GLB_SET,imolty,iatom+4,-0.08d0)
            atomname(iatom+5) = 'CN7'
            call twh_qqatom(GLB_SET,imolty,iatom+5,0.16d0)
            if (deoxy .eq. 0) then
               atomname(iatom+6) = 'ON6B'
               call twh_qqatom(GLB_SET,imolty,iatom+6,-0.50d0)
            else if (deoxy .eq. 1) then
               atomname(iatom+6) = 'ON6'
               call twh_qqatom(GLB_SET,imolty,iatom+6,-0.50d0)
            endif            
            atomname(iatom+7) = 'CN7B'
            call twh_qqatom(GLB_SET,imolty,iatom+7,0.16d0)
            if (deoxy .eq. 0) then
               atomname(iatom+8) = 'CN7B'
               call twh_qqatom(GLB_SET,imolty,iatom+8,0.14d0)
            elseif (deoxy .eq. 1) then
               atomname(iatom+8) = 'CN8'
               call twh_qqatom(GLB_SET,imolty,iatom+8,-0.18d0)
            endif
            atomname(iatom+9) = 'CN7'
            call twh_qqatom(GLB_SET,imolty,iatom+9,0.01d0)
            atomname(iatom+10) =   'ON2'
            call twh_qqatom(GLB_SET,imolty,iatom+10,-0.62d0)
            atomname(iatom+11) = 'P'
            call twh_qqatom(GLB_SET,imolty,iatom+11,1.50d0)
            atomname(iatom+12) = 'ON3'
            call twh_qqatom(GLB_SET,imolty,iatom+12,-0.82d0)
            atomname(iatom+13) = 'ON3'
            call twh_qqatom(GLB_SET,imolty,iatom+13,-0.82d0)
            atomname(iatom+14) = 'ON4'
            call twh_qqatom(GLB_SET,imolty,iatom+14,-0.68d0)
            if (deoxy .eq. 0) then
               atomname(iatom+15) = 'ON5'
               call twh_qqatom(GLB_SET,imolty,iatom+15,-0.66d0)
            endif
            atomname(hmain+heavy) = 'HN8'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy,0.09d0)
            atomname(hmain+heavy+1) = 'HN8'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+1,0.09d0)
            atomname(hmain+heavy+2) = 'HN7'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+2,0.09d0)
            atomname(hmain+heavy+3) = 'HN7'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+3,0.09d0)
            if (deoxy .eq. 0) then
               atomname(hmain+heavy+5) = 'HN5'
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+5,0.43d0)
               atomname(hmain+heavy+4) = 'HN7'
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+4,0.09d0)
            elseif (deoxy .eq. 1) then
               atomname(hmain+heavy+5) = 'HN8'
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+5,0.09d0)
               atomname(hmain+heavy+4) = 'HN8'
               call twh_qqatom(GLB_SET,imolty,hmain+heavy+4,0.09d0)
            endif
            atomname(hmain+heavy+6) = 'HN7'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+6,0.09d0)
            atomname(hmain+heavy+7) = 'HN4'
            call twh_qqatom(GLB_SET,imolty,hmain+heavy+7,0.34d0)
         endif
      else

c     --- invalid nucleic acid designation
         write(6,*) 'imolty,inucleic,pepname',imolty,inucleic
     &        ,pepname(inucleic),' is not valid'
         lfailure = .true.
         return
      endif
c     --- update ibondend for the typical nucleic acids
      if (termflag .eq. '5 prime' ) then
         if (hterm .eq. 5) then
            ibondend = iatom+7
         elseif (hterm .eq. 3) then
            ibondend = iatom+11
         endif
      else 
         ibondend = iatom+10
      endif
      
c     write(6,*) 'finished POLYCHARMM27NA'
      
      return
      end
#else
c     Compilation disabled (for speed) with FFLAGS=-DDISABLE_POLY
c     statement in ./configure
      write(6,*) "polycharmm27na disabled (compile flag)"
      lfailure = .true.
      return
      end
#endif
