#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2006-2011 Marcus G. Martin                       *
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

c     ******************************************************************
c     * contains the shared data for monomers used in the polypeptide  *
c     * builder and others                                             *
c     *                                                                *
c     * originally written 03-25-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      subroutine twh_aa_nucleic_main_chain(lfailure,termflag,hterm
     &     ,total,heavy,ninsert,hmain,light,iatom,shift,deoxy,imolty
     &     ,ibondend)
c     ******************************************************************
c     * sets up the communal bond information for the nucleic acid main*
c     * chain portion of the monomers for the nucleic acid builder     *
c     *                                                                *
c     * split from polycharmm27na 03-27-2006 by M.G. Martin            *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*20 termflag
      logical lfailure
      integer hterm,total,heavy,ninsert,hmain,light,iatom,shift,deoxy
      integer imolty,ibondend

c     --- generic nucleic acid
      total = heavy+ninsert+hmain+light-iatom

      if ( termflag .eq. '5 prime' ) then
         if (hterm .eq. 5) then
c           --- 5'-terminus main chain
c           --- vibrations
            call twh_invib(GLB_SET,imolty,iatom,2)
            call twh_ijvib(GLB_SET,imolty,iatom,1,hmain+heavy)
            call twh_ijvib(GLB_SET,imolty,iatom,2,iatom+1)
            call twh_invib(GLB_SET,imolty,iatom+1,4)
            call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
            call twh_ijvib(GLB_SET,imolty,iatom+1,2,iatom+2)
            call twh_ijvib(GLB_SET,imolty,iatom+1,3,hmain+heavy+1)
            call twh_ijvib(GLB_SET,imolty,iatom+1,4,hmain+heavy+2)
            call twh_invib(GLB_SET,imolty,iatom+2,4)
            call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom+1)
            call twh_ijvib(GLB_SET,imolty,iatom+2,2,iatom+3)
            call twh_ijvib(GLB_SET,imolty,iatom+2,3,iatom+6)
            call twh_ijvib(GLB_SET,imolty,iatom+2,4,hmain+heavy+3)
            call twh_invib(GLB_SET,imolty,iatom+3,2)
            call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom+2)
            call twh_ijvib(GLB_SET,imolty,iatom+3,2,iatom+4)
            call twh_invib(GLB_SET,imolty,iatom+4,4)
            call twh_ijvib(GLB_SET,imolty,iatom+4,1,iatom+3)
            call twh_ijvib(GLB_SET,imolty,iatom+4,2,iatom+5)
            call twh_ijvib(GLB_SET,imolty,iatom+4,3,hmain+heavy+4)
            if ( heavy .ne. 0 ) then
               call twh_ijvib(GLB_SET,imolty,iatom+4,4,hmain)
            else
               call twh_ijvib(GLB_SET,imolty,iatom+4,4,shift)
            endif
            call twh_invib(GLB_SET,imolty,iatom+5,4)
            call twh_ijvib(GLB_SET,imolty,iatom+5,1,iatom+4)
            call twh_ijvib(GLB_SET,imolty,iatom+5,2,iatom+6)
            call twh_ijvib(GLB_SET,imolty,iatom+5,3,hmain+heavy+5)
c           --- if doing deoxidized nucleic acids then don't have O2' 
            if (deoxy .eq. 0) then          
               call twh_ijvib(GLB_SET,imolty,iatom+5,4,iatom+8)
            elseif (deoxy .eq. 1) then
               call twh_ijvib(GLB_SET,imolty,iatom+5,4,hmain+heavy+6)
            else   
               write(6,*) "Illegal value of deoxy",deoxy
               lfailure = .true.
               return
            endif
            call twh_invib(GLB_SET,imolty,iatom+6,4)
            call twh_ijvib(GLB_SET,imolty,iatom+6,1,iatom+5)
            call twh_ijvib(GLB_SET,imolty,iatom+6,2,iatom+2)
            call twh_ijvib(GLB_SET,imolty,iatom+6,3,iatom+7)
            call twh_ijvib(GLB_SET,imolty,iatom+6,4,hmain+heavy+7)
            call twh_invib(GLB_SET,imolty,iatom+7,2)
            call twh_ijvib(GLB_SET,imolty,iatom+7,1,iatom+6)
            call twh_ijvib(GLB_SET,imolty,iatom+7,2,iatom+total)
c           --- if doing deoxidized nucleic acids then don't have O2' 
            if (deoxy .eq. 0) then
               call twh_invib(GLB_SET,imolty,iatom+8,2)
               call twh_ijvib(GLB_SET,imolty,iatom+8,1,iatom+5)
               call twh_ijvib(GLB_SET,imolty,iatom+8,2,hmain+heavy+6)
            endif                          
            call twh_invib(GLB_SET,imolty,hmain+heavy,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy,1,iatom)
            call twh_invib(GLB_SET,imolty,hmain+heavy+1,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+1,1,iatom+1)
            call twh_invib(GLB_SET,imolty,hmain+heavy+2,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+2,1,iatom+1)
            call twh_invib(GLB_SET,imolty,hmain+heavy+3,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+3,1,iatom+2)
            call twh_invib(GLB_SET,imolty,hmain+heavy+4,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+4,1,iatom+4)
            call twh_invib(GLB_SET,imolty,hmain+heavy+5,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+5,1,iatom+5)
            call twh_invib(GLB_SET,imolty,hmain+heavy+6,1)
            if (deoxy .eq. 0) then
               call twh_ijvib(GLB_SET,imolty,hmain+heavy+6,1,iatom+8)
            elseif (deoxy .eq. 1) then
               call twh_ijvib(GLB_SET,imolty,hmain+heavy+6,1,iatom+5)
            else
               write(6,*) "Illegal value of deoxy",deoxy
               lfailure = .true.
               return
            endif
            call twh_invib(GLB_SET,imolty,hmain+heavy+7,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+7,1,iatom+6)
c           --- pdb names
            call twh_pdbname(GLB_SET,imolty,iatom," O5'")
            call twh_pdbname(GLB_SET,imolty,iatom+1," C5'")
            call twh_pdbname(GLB_SET,imolty,iatom+2," C4'")
            call twh_pdbname(GLB_SET,imolty,iatom+3," O4'")
            call twh_pdbname(GLB_SET,imolty,iatom+4," C1'")
            call twh_pdbname(GLB_SET,imolty,iatom+5," C2'")
            call twh_pdbname(GLB_SET,imolty,iatom+6," C3'")
            call twh_pdbname(GLB_SET,imolty,iatom+7," O3'")
            if (deoxy .eq. 0) then
               call twh_pdbname(GLB_SET,imolty,iatom+8," O2'")
            endif
            call twh_pdbname(GLB_SET,imolty,hmain+heavy," H5T")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+1," H5'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+2,"H5''")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+3," H4'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+4," H1'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+5,"H2''")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+6," H2'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+7," H3'")
         elseif (hterm .eq. 3) then 
c           --- 3'-terminus main chain
c           --- vibrations
            call twh_invib(GLB_SET,imolty,iatom,4)
            call twh_ijvib(GLB_SET,imolty,iatom,1,iatom+1)
            call twh_ijvib(GLB_SET,imolty,iatom,2,iatom+2)
            call twh_ijvib(GLB_SET,imolty,iatom,3,iatom+3)
            call twh_ijvib(GLB_SET,imolty,iatom,4,iatom+4)
            call twh_invib(GLB_SET,imolty,iatom+1,2)
            call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
            call twh_ijvib(GLB_SET,imolty,iatom+1,2,hmain+heavy)
            call twh_invib(GLB_SET,imolty,iatom+2,1)
            call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom)
            call twh_invib(GLB_SET,imolty,iatom+3,1)
            call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom)
            call twh_invib(GLB_SET,imolty,iatom+4,2)
            call twh_ijvib(GLB_SET,imolty,iatom+4,1,iatom)
            call twh_ijvib(GLB_SET,imolty,iatom+4,2,iatom+5)
            call twh_invib(GLB_SET,imolty,iatom+5,4)
            call twh_ijvib(GLB_SET,imolty,iatom+5,1,iatom+4)
            call twh_ijvib(GLB_SET,imolty,iatom+5,2,iatom+6)
            call twh_ijvib(GLB_SET,imolty,iatom+5,3,hmain+heavy+1)
            call twh_ijvib(GLB_SET,imolty,iatom+5,4,hmain+heavy+2)
            call twh_invib(GLB_SET,imolty,iatom+6,4)
            call twh_ijvib(GLB_SET,imolty,iatom+6,1,iatom+5)
            call twh_ijvib(GLB_SET,imolty,iatom+6,2,iatom+7)
            call twh_ijvib(GLB_SET,imolty,iatom+6,3,iatom+10)
            call twh_ijvib(GLB_SET,imolty,iatom+6,4,hmain+heavy+3)
            call twh_invib(GLB_SET,imolty,iatom+7,2)
            call twh_ijvib(GLB_SET,imolty,iatom+7,1,iatom+6)
            call twh_ijvib(GLB_SET,imolty,iatom+7,2,iatom+8)
            call twh_invib(GLB_SET,imolty,iatom+8,4)
            call twh_ijvib(GLB_SET,imolty,iatom+8,1,iatom+7)
            call twh_ijvib(GLB_SET,imolty,iatom+8,2,iatom+9)
            call twh_ijvib(GLB_SET,imolty,iatom+8,3,hmain+heavy+4)
            if ( heavy .ne. 0 ) then
               call twh_ijvib(GLB_SET,imolty,iatom+8,4,hmain)
            else
               call twh_ijvib(GLB_SET,imolty,iatom+8,4,shift)
            endif
            call twh_invib(GLB_SET,imolty,iatom+9,4)
            call twh_ijvib(GLB_SET,imolty,iatom+9,1,iatom+8)
            call twh_ijvib(GLB_SET,imolty,iatom+9,2,iatom+10)
            call twh_ijvib(GLB_SET,imolty,iatom+9,3,hmain+heavy+5)
c           --- if doing deoxidized nucleic acids then don't have O2' 
            if (deoxy .eq. 0) then          
               call twh_ijvib(GLB_SET,imolty,iatom+9,4,iatom+12)
            elseif (deoxy .eq. 1) then
               call twh_ijvib(GLB_SET,imolty,iatom+9,4,hmain+heavy+6)
            else   
               write(6,*) "Illegal value of deoxy",deoxy
               lfailure = .true.
               return
            endif
            call twh_invib(GLB_SET,imolty,iatom+10,4)
            call twh_ijvib(GLB_SET,imolty,iatom+10,1,iatom+9)
            call twh_ijvib(GLB_SET,imolty,iatom+10,2,iatom+6)
            call twh_ijvib(GLB_SET,imolty,iatom+10,3,iatom+11)
            call twh_ijvib(GLB_SET,imolty,iatom+10,4,hmain+heavy+7)
            call twh_invib(GLB_SET,imolty,iatom+11,2)
            call twh_ijvib(GLB_SET,imolty,iatom+11,1,iatom+10)
            call twh_ijvib(GLB_SET,imolty,iatom+11,2,iatom+total)
c           --- if doing deoxidized nucleic acids then don't have O2' 
            if (deoxy .eq. 0) then
               call twh_invib(GLB_SET,imolty,iatom+12,2)
               call twh_ijvib(GLB_SET,imolty,iatom+12,1,iatom+9)
               call twh_ijvib(GLB_SET,imolty,iatom+12,2,hmain+heavy+6)
            endif                          
            call twh_invib(GLB_SET,imolty,hmain+heavy,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy,1,iatom+1)
            call twh_invib(GLB_SET,imolty,hmain+heavy+1,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+1,1,iatom+5)
            call twh_invib(GLB_SET,imolty,hmain+heavy+2,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+2,1,iatom+5)
            call twh_invib(GLB_SET,imolty,hmain+heavy+3,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+3,1,iatom+6)
            call twh_invib(GLB_SET,imolty,hmain+heavy+4,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+4,1,iatom+8)
            call twh_invib(GLB_SET,imolty,hmain+heavy+5,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+5,1,iatom+9)
            call twh_invib(GLB_SET,imolty,hmain+heavy+6,1)
            if (deoxy .eq. 0) then
               call twh_ijvib(GLB_SET,imolty,hmain+heavy+6,1,iatom+12)
            elseif (deoxy .eq. 1) then
               call twh_ijvib(GLB_SET,imolty,hmain+heavy+6,1,iatom+9)
            else
               write(6,*) "Illegal value of deoxy",deoxy
               lfailure = .true.
               return
            endif
            call twh_invib(GLB_SET,imolty,hmain+heavy+7,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+7,1,iatom+10)
c           --- pdb names
            call twh_pdbname(GLB_SET,imolty,iatom," P")
            call twh_pdbname(GLB_SET,imolty,iatom+1," O5T")
            call twh_pdbname(GLB_SET,imolty,iatom+2," O1P")
            call twh_pdbname(GLB_SET,imolty,iatom+3," O2P")
            call twh_pdbname(GLB_SET,imolty,iatom+4," O5'")
            call twh_pdbname(GLB_SET,imolty,iatom+5," C5'")
            call twh_pdbname(GLB_SET,imolty,iatom+6," C4'")
            call twh_pdbname(GLB_SET,imolty,iatom+7," O4'")
            call twh_pdbname(GLB_SET,imolty,iatom+8," C1'")
            call twh_pdbname(GLB_SET,imolty,iatom+9," C2'")
            call twh_pdbname(GLB_SET,imolty,iatom+10," C3'")
            call twh_pdbname(GLB_SET,imolty,iatom+11," O3'")
            if (deoxy .eq. 0) then
               call twh_pdbname(GLB_SET,imolty,iatom+12," O2'")
            endif
            call twh_pdbname(GLB_SET,imolty,hmain+heavy," H5T")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+1," H5'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+2,"H5''")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+3," H4'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+4," H1'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+5,"H2''")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+6," H2'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+7," H3'")
         endif
      elseif ( termflag .eq. 'neither' ) then
c        --- interior segment main chain
c        --- vibrations
         call twh_invib(GLB_SET,imolty,iatom,4)
         call twh_ijvib(GLB_SET,imolty,iatom,1,iatom+3)
         call twh_ijvib(GLB_SET,imolty,iatom,2,ibondend)
         call twh_ijvib(GLB_SET,imolty,iatom,3,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom,4,iatom+2)
         call twh_invib(GLB_SET,imolty,iatom+1,1)
         call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
         call twh_invib(GLB_SET,imolty,iatom+2,1)
         call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom)
         call twh_invib(GLB_SET,imolty,iatom+3,2)
         call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom)
         call twh_ijvib(GLB_SET,imolty,iatom+3,2,iatom+4)
         call twh_invib(GLB_SET,imolty,iatom+4,4)
         call twh_ijvib(GLB_SET,imolty,iatom+4,1,iatom+3)
         call twh_ijvib(GLB_SET,imolty,iatom+4,2,iatom+5)
         call twh_ijvib(GLB_SET,imolty,iatom+4,3,hmain+heavy)
         call twh_ijvib(GLB_SET,imolty,iatom+4,4,hmain+heavy+1)
         call twh_invib(GLB_SET,imolty,iatom+5,4)
         call twh_ijvib(GLB_SET,imolty,iatom+5,1,iatom+4)
         call twh_ijvib(GLB_SET,imolty,iatom+5,2,iatom+6)
         call twh_ijvib(GLB_SET,imolty,iatom+5,3,iatom+9)
         call twh_ijvib(GLB_SET,imolty,iatom+5,4,hmain+heavy+2)
         call twh_invib(GLB_SET,imolty,iatom+6,2)
         call twh_ijvib(GLB_SET,imolty,iatom+6,1,iatom+7)
         call twh_ijvib(GLB_SET,imolty,iatom+6,2,iatom+5)
         call twh_invib(GLB_SET,imolty,iatom+7,4)
         call twh_ijvib(GLB_SET,imolty,iatom+7,1,iatom+6)
         call twh_ijvib(GLB_SET,imolty,iatom+7,2,iatom+8)
         call twh_ijvib(GLB_SET,imolty,iatom+7,3,hmain+heavy+3)
         if ( heavy .ne. 0 ) then
             call twh_ijvib(GLB_SET,imolty,iatom+7,4,hmain)
         else
             call twh_ijvib(GLB_SET,imolty,iatom+7,4,shift)
         endif
         call twh_invib(GLB_SET,imolty,iatom+8,4)
         call twh_ijvib(GLB_SET,imolty,iatom+8,1,iatom+7)
         call twh_ijvib(GLB_SET,imolty,iatom+8,2,iatom+9)
         call twh_ijvib(GLB_SET,imolty,iatom+8,3,hmain+heavy+4)
c        --- if doing deoxidized nucleic acids then don't have O2' 
         if (deoxy .eq. 0) then
            call twh_ijvib(GLB_SET,imolty,iatom+8,4,iatom+11)
         elseif (deoxy .eq. 1) then
            call twh_ijvib(GLB_SET,imolty,iatom+8,4,hmain+heavy+5)
         else
            write(6,*) "Illegal value of deoxy",deoxy
            lfailure = .true.
            return
         endif
         call twh_invib(GLB_SET,imolty,iatom+9,4)
         call twh_ijvib(GLB_SET,imolty,iatom+9,1,iatom+10)
         call twh_ijvib(GLB_SET,imolty,iatom+9,2,iatom+8)
         call twh_ijvib(GLB_SET,imolty,iatom+9,3,iatom+5)
         call twh_ijvib(GLB_SET,imolty,iatom+9,4,hmain+heavy+6)
         call twh_invib(GLB_SET,imolty,iatom+10,2)
         call twh_ijvib(GLB_SET,imolty,iatom+10,1,iatom+9)
         call twh_ijvib(GLB_SET,imolty,iatom+10,2,iatom + total)
c        --- if doing deoxidized nucleic acids then don't have O2' 
         if (deoxy .eq. 0) then
             call twh_invib(GLB_SET,imolty,iatom+11,2)
             call twh_ijvib(GLB_SET,imolty,iatom+11,1,iatom+8)
             call twh_ijvib(GLB_SET,imolty,iatom+11,2,hmain+heavy+5)
         endif                          
         call twh_invib(GLB_SET,imolty,hmain+heavy,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy,1,iatom+4)
         call twh_invib(GLB_SET,imolty,hmain+heavy+1,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy+1,1,iatom+4)
         call twh_invib(GLB_SET,imolty,hmain+heavy+2,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy+2,1,iatom+5)
         call twh_invib(GLB_SET,imolty,hmain+heavy+3,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy+3,1,iatom+7)
         call twh_invib(GLB_SET,imolty,hmain+heavy+4,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy+4,1,iatom+8)
         call twh_invib(GLB_SET,imolty,hmain+heavy+5,1)
         if (deoxy .eq. 0) then
             call twh_ijvib(GLB_SET,imolty,hmain+heavy+5,1,iatom+11)
         elseif (deoxy .eq. 1) then 
             call twh_ijvib(GLB_SET,imolty,hmain+heavy+5,1,iatom+8)
         endif
         call twh_invib(GLB_SET,imolty,hmain+heavy+6,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy+6,1,iatom+9)
c        --- pdb names
         call twh_pdbname(GLB_SET,imolty,iatom," P")
         call twh_pdbname(GLB_SET,imolty,iatom+1," O1P")
         call twh_pdbname(GLB_SET,imolty,iatom+2," O2P")
         call twh_pdbname(GLB_SET,imolty,iatom+3," O5'")
         call twh_pdbname(GLB_SET,imolty,iatom+4," C5'")
         call twh_pdbname(GLB_SET,imolty,iatom+5," C4'")
         call twh_pdbname(GLB_SET,imolty,iatom+6," O4'")
         call twh_pdbname(GLB_SET,imolty,iatom+7," C1'")
         call twh_pdbname(GLB_SET,imolty,iatom+8," C2'")
         call twh_pdbname(GLB_SET,imolty,iatom+9," C3'")
         call twh_pdbname(GLB_SET,imolty,iatom+10," O3'")
         if (deoxy .eq. 0) then
            call twh_pdbname(GLB_SET,imolty,iatom+11," O2'")
         endif
         call twh_pdbname(GLB_SET,imolty,hmain+heavy," H5'")
         call twh_pdbname(GLB_SET,imolty,hmain+heavy+1,"H5''")
         call twh_pdbname(GLB_SET,imolty,hmain+heavy+2," H4'")
         call twh_pdbname(GLB_SET,imolty,hmain+heavy+3," H1'")
         call twh_pdbname(GLB_SET,imolty,hmain+heavy+4,"H2''")
         call twh_pdbname(GLB_SET,imolty,hmain+heavy+5," H2'")
         call twh_pdbname(GLB_SET,imolty,hmain+heavy+6," H3'")
      elseif ( termflag .eq. '3 prime' ) then
         if (hterm .eq. 3) then
c           --- 3'-terminus main chain
c           --- vibrations
            call twh_invib(GLB_SET,imolty,iatom,4)
            call twh_ijvib(GLB_SET,imolty,iatom,1,iatom+1)
            call twh_ijvib(GLB_SET,imolty,iatom,2,iatom+2)
            call twh_ijvib(GLB_SET,imolty,iatom,3,iatom+3)
            call twh_ijvib(GLB_SET,imolty,iatom,4,ibondend)
            call twh_invib(GLB_SET,imolty,iatom+1,1)
            call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
            call twh_invib(GLB_SET,imolty,iatom+2,1)
            call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom)
            call twh_invib(GLB_SET,imolty,iatom+3,2)
            call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom)
            call twh_ijvib(GLB_SET,imolty,iatom+3,2,iatom+4)
            call twh_invib(GLB_SET,imolty,iatom+4,4)
            call twh_ijvib(GLB_SET,imolty,iatom+4,1,iatom+3)
            call twh_ijvib(GLB_SET,imolty,iatom+4,2,iatom+5)
            call twh_ijvib(GLB_SET,imolty,iatom+4,3,hmain+heavy)
            call twh_ijvib(GLB_SET,imolty,iatom+4,4,hmain+heavy+1)
            call twh_invib(GLB_SET,imolty,iatom+5,4)
            call twh_ijvib(GLB_SET,imolty,iatom+5,1,iatom+4)
            call twh_ijvib(GLB_SET,imolty,iatom+5,2,iatom+6)
            call twh_ijvib(GLB_SET,imolty,iatom+5,3,iatom+9)
            call twh_ijvib(GLB_SET,imolty,iatom+5,4,hmain+heavy+2)
            call twh_invib(GLB_SET,imolty,iatom+6,2)
            call twh_ijvib(GLB_SET,imolty,iatom+6,1,iatom+5)
            call twh_ijvib(GLB_SET,imolty,iatom+6,2,iatom+7)
            call twh_invib(GLB_SET,imolty,iatom+7,4)
            call twh_ijvib(GLB_SET,imolty,iatom+7,1,iatom+6)
            call twh_ijvib(GLB_SET,imolty,iatom+7,2,iatom+8)
            call twh_ijvib(GLB_SET,imolty,iatom+7,3,hmain+heavy+3)
            if ( heavy .ne. 0 ) then
               call twh_ijvib(GLB_SET,imolty,iatom+7,4,hmain)
            else
               call twh_ijvib(GLB_SET,imolty,iatom+7,4,shift)
            endif
            call twh_invib(GLB_SET,imolty,iatom+8,4)
            call twh_ijvib(GLB_SET,imolty,iatom+8,1,iatom+7)
            call twh_ijvib(GLB_SET,imolty,iatom+8,2,iatom+9)
            call twh_ijvib(GLB_SET,imolty,iatom+8,3,hmain+heavy+4)
c           --- if doing deoxidized nucleic acids then don't have O2' 
            if (deoxy .eq. 0) then
               call twh_ijvib(GLB_SET,imolty,iatom+8,4,iatom+11)
            elseif (deoxy .eq. 1) then
               call twh_ijvib(GLB_SET,imolty,iatom+8,4,hmain+heavy+5)
            else 
               write(6,*) "Illegal value of deoxy",deoxy
               lfailure = .true.
               return
            endif
            call twh_invib(GLB_SET,imolty,iatom+9,4)
            call twh_ijvib(GLB_SET,imolty,iatom+9,1,iatom+8)
            call twh_ijvib(GLB_SET,imolty,iatom+9,2,iatom+5)
            call twh_ijvib(GLB_SET,imolty,iatom+9,3,iatom+10)
            call twh_ijvib(GLB_SET,imolty,iatom+9,4,hmain+heavy+6)
            call twh_invib(GLB_SET,imolty,iatom+10,2)
            call twh_ijvib(GLB_SET,imolty,iatom+10,1,iatom+9)
            call twh_ijvib(GLB_SET,imolty,iatom+10,2,hmain+heavy+7)
c           --- if doing deoxidized nucleic acids then don't have O2' 
            if (deoxy .eq. 0) then
               call twh_invib(GLB_SET,imolty,iatom+11,2)
               call twh_ijvib(GLB_SET,imolty,iatom+11,1,iatom+8)
               call twh_ijvib(GLB_SET,imolty,iatom+11,2,hmain+heavy+5)
            endif                          
            call twh_invib(GLB_SET,imolty,hmain+heavy,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy,1,iatom+4)
            call twh_invib(GLB_SET,imolty,hmain+heavy+1,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+1,1,iatom+4)
            call twh_invib(GLB_SET,imolty,hmain+heavy+2,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+2,1,iatom+5)
            call twh_invib(GLB_SET,imolty,hmain+heavy+3,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+3,1,iatom+7)
            call twh_invib(GLB_SET,imolty,hmain+heavy+4,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+4,1,iatom+8)
            call twh_invib(GLB_SET,imolty,hmain+heavy+5,1)
            if (deoxy .eq. 0) then
               call twh_ijvib(GLB_SET,imolty,hmain+heavy+5,1,iatom+11)
            elseif (deoxy .eq. 1) then 
               call twh_ijvib(GLB_SET,imolty,hmain+heavy+5,1,iatom+8)
            endif
            call twh_invib(GLB_SET,imolty,hmain+heavy+6,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+6,1,iatom+9)
            call twh_invib(GLB_SET,imolty,hmain+heavy+7,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+7,1,iatom+10)
c           --- pdb names
            call twh_pdbname(GLB_SET,imolty,iatom," P")
            call twh_pdbname(GLB_SET,imolty,iatom+1," O1P")
            call twh_pdbname(GLB_SET,imolty,iatom+2," O2P")
            call twh_pdbname(GLB_SET,imolty,iatom+3," O5'")
            call twh_pdbname(GLB_SET,imolty,iatom+4," C5'")
            call twh_pdbname(GLB_SET,imolty,iatom+5," C4'")
            call twh_pdbname(GLB_SET,imolty,iatom+6," O4'")
            call twh_pdbname(GLB_SET,imolty,iatom+7," C1'")
            call twh_pdbname(GLB_SET,imolty,iatom+8," C2'")
            call twh_pdbname(GLB_SET,imolty,iatom+9," C3'")
            call twh_pdbname(GLB_SET,imolty,iatom+10," O3'")
            if (deoxy .eq. 0) then
               call twh_pdbname(GLB_SET,imolty,iatom+11," O2'")
            endif
            call twh_pdbname(GLB_SET,imolty,hmain+heavy," H5'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+1,"H5''")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+2," H4'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+3," H1'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+4,"H2''")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+5," H2'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+6," H3'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+7," H3T")
         elseif (hterm .eq. 5) then
c           --- 5'-terminus main chain
c           --- vibrations
            call twh_invib(GLB_SET,imolty,iatom,4)
            call twh_ijvib(GLB_SET,imolty,iatom,1,iatom+1)
            call twh_ijvib(GLB_SET,imolty,iatom,2,iatom+2)
            call twh_ijvib(GLB_SET,imolty,iatom,3,iatom+3)
            call twh_ijvib(GLB_SET,imolty,iatom,4,ibondend)
            call twh_invib(GLB_SET,imolty,iatom+1,1)
            call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
            call twh_invib(GLB_SET,imolty,iatom+2,1)
            call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom)
            call twh_invib(GLB_SET,imolty,iatom+3,2)
            call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom)
            call twh_ijvib(GLB_SET,imolty,iatom+3,2,iatom+4)
            call twh_invib(GLB_SET,imolty,iatom+4,4)
            call twh_ijvib(GLB_SET,imolty,iatom+4,1,iatom+3)
            call twh_ijvib(GLB_SET,imolty,iatom+4,2,iatom+5)
            call twh_ijvib(GLB_SET,imolty,iatom+4,3,hmain+heavy)
            call twh_ijvib(GLB_SET,imolty,iatom+4,4,hmain+heavy+1)
            call twh_invib(GLB_SET,imolty,iatom+5,4)
            call twh_ijvib(GLB_SET,imolty,iatom+5,1,iatom+4)
            call twh_ijvib(GLB_SET,imolty,iatom+5,2,iatom+6)
            call twh_ijvib(GLB_SET,imolty,iatom+5,3,iatom+9)
            call twh_ijvib(GLB_SET,imolty,iatom+5,4,hmain+heavy+2)
            call twh_invib(GLB_SET,imolty,iatom+6,2)
            call twh_ijvib(GLB_SET,imolty,iatom+6,1,iatom+5)
            call twh_ijvib(GLB_SET,imolty,iatom+6,2,iatom+7)
            call twh_invib(GLB_SET,imolty,iatom+7,4)
            call twh_ijvib(GLB_SET,imolty,iatom+7,1,iatom+6)
            call twh_ijvib(GLB_SET,imolty,iatom+7,2,iatom+8)
            call twh_ijvib(GLB_SET,imolty,iatom+7,3,hmain+heavy+3)
            if ( heavy .ne. 0 ) then
               call twh_ijvib(GLB_SET,imolty,iatom+7,4,hmain)
            else
               call twh_ijvib(GLB_SET,imolty,iatom+7,4,shift)
            endif
            call twh_invib(GLB_SET,imolty,iatom+8,4)
            call twh_ijvib(GLB_SET,imolty,iatom+8,1,iatom+7)
            call twh_ijvib(GLB_SET,imolty,iatom+8,2,iatom+9)
            call twh_ijvib(GLB_SET,imolty,iatom+8,3,hmain+heavy+4)
c           --- if doing deoxidized nucleic acids then don't have O2' 
            if (deoxy .eq. 0) then
               call twh_ijvib(GLB_SET,imolty,iatom+8,4,iatom+15)
            elseif (deoxy .eq. 1) then
               call twh_ijvib(GLB_SET,imolty,iatom+8,4,hmain+heavy+5)
            else 
               write(6,*) "Illegal value of deoxy",deoxy
               lfailure = .true.
               return
            endif
            call twh_invib(GLB_SET,imolty,iatom+9,4)
            call twh_ijvib(GLB_SET,imolty,iatom+9,1,iatom+8)
            call twh_ijvib(GLB_SET,imolty,iatom+9,2,iatom+10)
            call twh_ijvib(GLB_SET,imolty,iatom+9,3,iatom+5)
            call twh_ijvib(GLB_SET,imolty,iatom+9,4,hmain+heavy+6)
            call twh_invib(GLB_SET,imolty,iatom+10,2)
            call twh_ijvib(GLB_SET,imolty,iatom+10,1,iatom+9)
            call twh_ijvib(GLB_SET,imolty,iatom+10,2,iatom+11)
            call twh_invib(GLB_SET,imolty,iatom+11,4)
            call twh_ijvib(GLB_SET,imolty,iatom+11,1,iatom+10)
            call twh_ijvib(GLB_SET,imolty,iatom+11,2,iatom+12)
            call twh_ijvib(GLB_SET,imolty,iatom+11,3,iatom+13)
            call twh_ijvib(GLB_SET,imolty,iatom+11,4,iatom+14)
            call twh_invib(GLB_SET,imolty,iatom+12,1)
            call twh_ijvib(GLB_SET,imolty,iatom+12,1,iatom+11)
            call twh_invib(GLB_SET,imolty,iatom+13,1)
            call twh_ijvib(GLB_SET,imolty,iatom+13,1,iatom+11)
            call twh_invib(GLB_SET,imolty,iatom+14,2)
            call twh_ijvib(GLB_SET,imolty,iatom+14,1,iatom+11)
            call twh_ijvib(GLB_SET,imolty,iatom+14,2,hmain+heavy+7)
c           --- if doing deoxidized nucleic acids then don't have O2' 
            if (deoxy .eq. 0) then
               call twh_invib(GLB_SET,imolty,iatom+15,2)
               call twh_ijvib(GLB_SET,imolty,iatom+15,1,iatom+8)
               call twh_ijvib(GLB_SET,imolty,iatom+15,2,hmain+heavy+5)
            endif                          
            call twh_invib(GLB_SET,imolty,hmain+heavy,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy,1,iatom+4)
            call twh_invib(GLB_SET,imolty,hmain+heavy+1,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+1,1,iatom+4)
            call twh_invib(GLB_SET,imolty,hmain+heavy+2,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+2,1,iatom+5)
            call twh_invib(GLB_SET,imolty,hmain+heavy+3,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+3,1,iatom+7)
            call twh_invib(GLB_SET,imolty,hmain+heavy+4,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+4,1,iatom+8)
            call twh_invib(GLB_SET,imolty,hmain+heavy+5,1)
            if (deoxy .eq. 0) then
               call twh_ijvib(GLB_SET,imolty,hmain+heavy+5,1,iatom+15)
            elseif (deoxy .eq. 1) then 
               call twh_ijvib(GLB_SET,imolty,hmain+heavy+5,1,iatom+8)
            endif
            call twh_invib(GLB_SET,imolty,hmain+heavy+6,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+6,1,iatom+9)
            call twh_invib(GLB_SET,imolty,hmain+heavy+7,1)
            call twh_ijvib(GLB_SET,imolty,hmain+heavy+7,1,iatom+14)
c           --- pdb names
            call twh_pdbname(GLB_SET,imolty,iatom," P")
            call twh_pdbname(GLB_SET,imolty,iatom+1," O1P")
            call twh_pdbname(GLB_SET,imolty,iatom+2," O2P")
            call twh_pdbname(GLB_SET,imolty,iatom+3," O5'")
            call twh_pdbname(GLB_SET,imolty,iatom+4," C5'")
            call twh_pdbname(GLB_SET,imolty,iatom+5," C4'")
            call twh_pdbname(GLB_SET,imolty,iatom+6," O4'")
            call twh_pdbname(GLB_SET,imolty,iatom+7," C1'")
            call twh_pdbname(GLB_SET,imolty,iatom+8," C2'")
            call twh_pdbname(GLB_SET,imolty,iatom+9," C3'")
            call twh_pdbname(GLB_SET,imolty,iatom+10," O3'")
            call twh_pdbname(GLB_SET,imolty,iatom+11," P3")
            call twh_pdbname(GLB_SET,imolty,iatom+12,"O1P3")
            call twh_pdbname(GLB_SET,imolty,iatom+13,"O2P3")
            call twh_pdbname(GLB_SET,imolty,iatom+14," O3T")
            if (deoxy .eq. 0) then
               call twh_pdbname(GLB_SET,imolty,iatom+15," O2'")
            endif
            call twh_pdbname(GLB_SET,imolty,hmain+heavy," H5'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+1,"H5''")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+2," H4'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+3," H1'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+4,"H2''")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+5," H2'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+6," H3'")
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+7," H3T")
         endif
      else
c        --- invalid termflag
         write(6,*) 'AA_NUCLEIC_MAIN_CHAIN: invalid termflag'
         write(6,*) termflag
         lfailure = .true.
         return
      endif

      return
      end


      subroutine twh_aa_nucleic_a(imolty,heavy,light,shift,hmain,ilinker
     &     ,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the nucleic acid 'a' *
c     * monomer in the nucleic acid builder                            *
c     * a: oxy or deoxy adenosine                                      *
c     * this monomers makes use of the nucleic acid main chain builder *
c     *                                                                *
c     * split from polycharmm27na 03-27-2006 by M.G. Martin            *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer imolty,heavy,light,shift,hmain,ilinker,ninsert

c     --- all-atom geometry
      heavy = 10
      shift = hmain+heavy+ninsert
      light = 4
c     --- vibrations
c     --- heavy atoms
      call twh_invib(GLB_SET,imolty,hmain,3)
      call twh_ijvib(GLB_SET,imolty,hmain,1,ilinker)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,hmain+9)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,shift)
      call twh_invib(GLB_SET,imolty,hmain+2,2)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+3,3)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,hmain+9)
      call twh_invib(GLB_SET,imolty,hmain+4,3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,hmain+6)
      call twh_invib(GLB_SET,imolty,hmain+5,3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+5,2,shift+1)
      call twh_ijvib(GLB_SET,imolty,hmain+5,3,shift+2)
      call twh_invib(GLB_SET,imolty,hmain+6,2)
      call twh_ijvib(GLB_SET,imolty,hmain+6,1,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+6,2,hmain+7)
      call twh_invib(GLB_SET,imolty,hmain+7,3)
      call twh_ijvib(GLB_SET,imolty,hmain+7,1,hmain+6)
      call twh_ijvib(GLB_SET,imolty,hmain+7,2,hmain+8)
      call twh_ijvib(GLB_SET,imolty,hmain+7,3,shift+3)
      call twh_invib(GLB_SET,imolty,hmain+8,2)
      call twh_ijvib(GLB_SET,imolty,hmain+8,1,hmain+7)
      call twh_ijvib(GLB_SET,imolty,hmain+8,2,hmain+9)
      call twh_invib(GLB_SET,imolty,hmain+9,3)
      call twh_ijvib(GLB_SET,imolty,hmain+9,1,hmain+8)
      call twh_ijvib(GLB_SET,imolty,hmain+9,2,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+9,3,hmain)
c     --- hydrogens
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain+5)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+5)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+7)
c     --- assign pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' N9')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' C8')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' N7')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' C5')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' C6')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' N6')
      call twh_pdbname(GLB_SET,imolty,hmain+6,' N1')
      call twh_pdbname(GLB_SET,imolty,hmain+7,' C2')
      call twh_pdbname(GLB_SET,imolty,hmain+8,' N3')
      call twh_pdbname(GLB_SET,imolty,hmain+9,' C4')
      call twh_pdbname(GLB_SET,imolty,shift,' H8')
      call twh_pdbname(GLB_SET,imolty,shift+1,'1H6')
      call twh_pdbname(GLB_SET,imolty,shift+2,'2H6')
      call twh_pdbname(GLB_SET,imolty,shift+3,' H2')

      return
      end

      subroutine twh_aa_nucleic_c(imolty,heavy,light,shift,hmain,ilinker
     &     ,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the nucleic acid 'c' *
c     * monomer in the nucleic acid builder                            *
c     * c: oxy or deoxy cytidine                                      *
c     * this monomers makes use of the nucleic acid main chain builder *
c     *                                                                *
c     * split from polycharmm27na 03-27-2006 by M.G. Martin            *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer imolty,heavy,light,shift,hmain,ilinker,ninsert

c     --- all-atom geometry
      heavy = 8
      shift = hmain+heavy+ninsert
      light = 4
c     --- vibrations
c     --- heavy atoms
      call twh_invib(GLB_SET,imolty,hmain,3)
      call twh_ijvib(GLB_SET,imolty,hmain,1,ilinker)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,hmain+7)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+3,2)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+4)
      call twh_invib(GLB_SET,imolty,hmain+4,3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,hmain+6)
      call twh_invib(GLB_SET,imolty,hmain+5,3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+5,2,shift)
      call twh_ijvib(GLB_SET,imolty,hmain+5,3,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+6,3)
      call twh_ijvib(GLB_SET,imolty,hmain+6,1,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+6,2,hmain+7)
      call twh_ijvib(GLB_SET,imolty,hmain+6,3,shift+2)
      call twh_invib(GLB_SET,imolty,hmain+7,3)
      call twh_ijvib(GLB_SET,imolty,hmain+7,1,hmain+6)
      call twh_ijvib(GLB_SET,imolty,hmain+7,2,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+7,3,shift+3)
c     --- hydrogens
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain+5)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain+5)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+6)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+7)
c     --- assign pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' N1')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' C2')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' O2')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' N3')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' C4')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' N4')
      call twh_pdbname(GLB_SET,imolty,hmain+6,' C5')
      call twh_pdbname(GLB_SET,imolty,hmain+7,' C6')
      call twh_pdbname(GLB_SET,imolty,shift,'1H4')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2H4')
      call twh_pdbname(GLB_SET,imolty,shift+2,' H5')
      call twh_pdbname(GLB_SET,imolty,shift+3,' H6')

      return
      end

      subroutine twh_aa_nucleic_g(imolty,heavy,light,shift,hmain,ilinker
     &     ,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the nucleic acid 'g' *
c     * monomer in the nucleic acid builder                            *
c     * g: oxy or deoxy guanosine                                      *
c     * this monomers makes use of the nucleic acid main chain builder *
c     *                                                                *
c     * split from polycharmm27na 03-27-2006 by M.G. Martin            *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer imolty,heavy,light,shift,hmain,ilinker,ninsert

c     --- all-atom geometry
      heavy = 11
      shift = hmain+heavy+ninsert
      light = 4
c     --- vibrations
c     --- heavy atoms
      call twh_invib(GLB_SET,imolty,hmain,3)
      call twh_ijvib(GLB_SET,imolty,hmain,1,ilinker)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,hmain+10)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,shift)
      call twh_invib(GLB_SET,imolty,hmain+2,2)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+3,3)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,hmain+10)
      call twh_invib(GLB_SET,imolty,hmain+4,3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,hmain+6)
      call twh_invib(GLB_SET,imolty,hmain+5,1)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+4)
      call twh_invib(GLB_SET,imolty,hmain+6,3)
      call twh_ijvib(GLB_SET,imolty,hmain+6,1,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+6,2,hmain+7)
      call twh_ijvib(GLB_SET,imolty,hmain+6,3,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+7,3)
      call twh_ijvib(GLB_SET,imolty,hmain+7,1,hmain+6)
      call twh_ijvib(GLB_SET,imolty,hmain+7,2,hmain+8)
      call twh_ijvib(GLB_SET,imolty,hmain+7,3,hmain+9)
      call twh_invib(GLB_SET,imolty,hmain+8,3)
      call twh_ijvib(GLB_SET,imolty,hmain+8,1,hmain+7)
      call twh_ijvib(GLB_SET,imolty,hmain+8,2,shift+2)
      call twh_ijvib(GLB_SET,imolty,hmain+8,3,shift+3)
      call twh_invib(GLB_SET,imolty,hmain+9,2)
      call twh_ijvib(GLB_SET,imolty,hmain+9,1,hmain+7)
      call twh_ijvib(GLB_SET,imolty,hmain+9,2,hmain+10)
      call twh_invib(GLB_SET,imolty,hmain+10,3)
      call twh_ijvib(GLB_SET,imolty,hmain+10,1,hmain+9)
      call twh_ijvib(GLB_SET,imolty,hmain+10,2,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+10,3,hmain)
c     --- hydrogens
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain+6)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+8)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+8)
c     --- assign pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' N9')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' C8')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' N7')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' C5')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' C6')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' O6')
      call twh_pdbname(GLB_SET,imolty,hmain+6,' N1')
      call twh_pdbname(GLB_SET,imolty,hmain+7,' C2')
      call twh_pdbname(GLB_SET,imolty,hmain+8,' N2')
      call twh_pdbname(GLB_SET,imolty,hmain+9,' N3')
      call twh_pdbname(GLB_SET,imolty,hmain+10,' C4')
      call twh_pdbname(GLB_SET,imolty,shift,' H8')
      call twh_pdbname(GLB_SET,imolty,shift+1,' H1')
      call twh_pdbname(GLB_SET,imolty,shift+2,'1H2')
      call twh_pdbname(GLB_SET,imolty,shift+3,'2H2')

      return
      end

      subroutine twh_aa_nucleic_t(imolty,heavy,light,shift,hmain,ilinker
     &     ,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the nucleic acid 't' *
c     * monomer in the nucleic acid builder                            *
c     * t: oxy or deoxy thymidine                                      *
c     * this monomers makes use of the nucleic acid main chain builder *
c     *                                                                *
c     * split from polycharmm27na 03-27-2006 by M.G. Martin            *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer imolty,heavy,light,shift,hmain,ilinker,ninsert

      heavy = 9
      shift = hmain+heavy+ninsert
      light = 5
c     --- vibrations
c     --- heavy atoms
      call twh_invib(GLB_SET,imolty,hmain,3)
      call twh_ijvib(GLB_SET,imolty,hmain,1,ilinker)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,hmain+8)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+3,3)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,shift)
      call twh_invib(GLB_SET,imolty,hmain+4,3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,hmain+6)
      call twh_invib(GLB_SET,imolty,hmain+5,1)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+4)
      call twh_invib(GLB_SET,imolty,hmain+6,3)
      call twh_ijvib(GLB_SET,imolty,hmain+6,1,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+6,2,hmain+7)
      call twh_ijvib(GLB_SET,imolty,hmain+6,3,hmain+8)
      call twh_invib(GLB_SET,imolty,hmain+7,4)
      call twh_ijvib(GLB_SET,imolty,hmain+7,1,hmain+6)
      call twh_ijvib(GLB_SET,imolty,hmain+7,2,shift+1)
      call twh_ijvib(GLB_SET,imolty,hmain+7,3,shift+2)
      call twh_ijvib(GLB_SET,imolty,hmain+7,4,shift+3)
      call twh_invib(GLB_SET,imolty,hmain+8,3)
      call twh_ijvib(GLB_SET,imolty,hmain+8,1,hmain+6)
      call twh_ijvib(GLB_SET,imolty,hmain+8,2,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+8,3,shift+4)
c     --- hydrogens
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain+7)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+7)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+7)
      call twh_invib(GLB_SET,imolty,shift+4,1)
      call twh_ijvib(GLB_SET,imolty,shift+4,1,hmain+8)
c     --- assign pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' N1')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' C2')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' O2')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' N3')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' C4')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' O4')
      call twh_pdbname(GLB_SET,imolty,hmain+6,' C5')
      call twh_pdbname(GLB_SET,imolty,hmain+7,'C5M')
      call twh_pdbname(GLB_SET,imolty,hmain+8,' C6')
      call twh_pdbname(GLB_SET,imolty,shift,' H3')
      call twh_pdbname(GLB_SET,imolty,shift+1,'1H5')
      call twh_pdbname(GLB_SET,imolty,shift+2,'2H5')
      call twh_pdbname(GLB_SET,imolty,shift+3,'3H5')
      call twh_pdbname(GLB_SET,imolty,shift+4,' H6')

      return
      end

      subroutine twh_aa_nucleic_u(imolty,heavy,light,shift,hmain,ilinker
     &     ,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the nucleic acid 'u' *
c     * monomer in the nucleic acid builder                            *
c     * u: oxy or deoxy uridine                                        *
c     * this monomers makes use of the nucleic acid main chain builder *
c     *                                                                *
c     * split from polycharmm27na 03-27-2006 by M.G. Martin            *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer imolty,heavy,light,shift,hmain,ilinker,ninsert

      heavy = 8
      shift = hmain+heavy+ninsert
      light = 3
c     --- vibrations
c     --- heavy atoms
      call twh_invib(GLB_SET,imolty,hmain,3)
      call twh_ijvib(GLB_SET,imolty,hmain,1,ilinker)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,hmain+7)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+3,3)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,shift)
      call twh_invib(GLB_SET,imolty,hmain+4,3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,hmain+6)
      call twh_invib(GLB_SET,imolty,hmain+5,1)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+4)
      call twh_invib(GLB_SET,imolty,hmain+6,3)
      call twh_ijvib(GLB_SET,imolty,hmain+6,1,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+6,2,hmain+7)
      call twh_ijvib(GLB_SET,imolty,hmain+6,3,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+7,3)
      call twh_ijvib(GLB_SET,imolty,hmain+7,1,hmain+6)
      call twh_ijvib(GLB_SET,imolty,hmain+7,2,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+7,3,shift+2)
c     --- hydrogens
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain+6)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+7)
c     --- assign pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' N1')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' C2')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' O2')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' N3')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' C4')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' O4')
      call twh_pdbname(GLB_SET,imolty,hmain+6,' C5')
      call twh_pdbname(GLB_SET,imolty,hmain+7,' C6')
      call twh_pdbname(GLB_SET,imolty,shift,' H3')
      call twh_pdbname(GLB_SET,imolty,shift+1,' H5')
      call twh_pdbname(GLB_SET,imolty,shift+2,' H6')

      return
      end

      subroutine twh_aa_pepname_generic(lfailure,termflag,hmain,iatom
     &     ,ninsert)
c     ******************************************************************
c     * sets up the generic information for most all-atom peptides     *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer termflag,hmain,iatom,ninsert

      if ( termflag .eq. -1 ) then
c        --- n-terminus
         hmain = iatom+4
c        --- all atom force field
         ninsert = 4
      elseif ( termflag .eq. 0 ) then
c        --- main chain
         hmain = iatom+4
c        --- all atom force field
         ninsert = 2
      elseif ( termflag .eq. 1 ) then
c        --- c-terminus
         hmain = iatom+4
c        --- all atom force field
         ninsert = 3
      else
         write(6,*) 'AA_PEPNAME_GENERIC: unknown termflag:',termflag
         lfailure = .true.
         return
      endif

      return
      end

      subroutine twh_aa_pepname_main_chain(lfailure,aminoname,termflag
     &     ,imolty,total,iatom,icb,iha,ibondend,heavy,light,ninsert
     &     ,shift,hmain)
c     ******************************************************************
c     * sets up the main chain generic information for most all-atom   *
c     * peptides                                                       *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*40 aminoname
      integer termflag,imolty,total,iatom,icb,iha,ibondend,heavy,shift
      integer hmain,ninsert,light

c     --- generic amino acid setup
      total = heavy+ninsert+hmain+light
      call twh_monomer_check(lfailure,total)
      if ( lfailure ) return
c     --- subtract iatom to get the number of new atoms for this monomer
      total = total - iatom

      if ( termflag .eq. -1 ) then
c        --- vibrations
c        --- set iha
         iha = hmain+heavy+3
c        --- n-terminus main chain
         call twh_invib(GLB_SET,imolty,iatom,4)
         call twh_ijvib(GLB_SET,imolty,iatom,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom,2,hmain+heavy)
         call twh_ijvib(GLB_SET,imolty,iatom,3,hmain+heavy+1)
         call twh_ijvib(GLB_SET,imolty,iatom,4,hmain+heavy+2)
         call twh_invib(GLB_SET,imolty,iatom+1,4)
         call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
         call twh_ijvib(GLB_SET,imolty,iatom+1,2,iatom+2)
         if ( heavy .ne. 0 ) then
            call twh_ijvib(GLB_SET,imolty,iatom+1,3,hmain)
c           --- set icb to the first atom of the chain
            icb = hmain
         else
            call twh_ijvib(GLB_SET,imolty,iatom+1,3,shift)
c           --- set icb to the hydrogen
            icb = shift
         endif
         call twh_ijvib(GLB_SET,imolty,iatom+1,4,hmain+heavy+3)
         call twh_invib(GLB_SET,imolty,iatom+2,3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom+2,2,iatom+3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,3,iatom+total)
         call twh_invib(GLB_SET,imolty,iatom+3,1)
         call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom+2)
         call twh_invib(GLB_SET,imolty,hmain+heavy,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy,1,iatom)
         call twh_invib(GLB_SET,imolty,hmain+heavy+1,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy+1,1,iatom)
         call twh_invib(GLB_SET,imolty,hmain+heavy+2,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy+2,1,iatom)
         call twh_invib(GLB_SET,imolty,hmain+heavy+3,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy+3,1,iatom+1)
c        --- pdb names
         call twh_pdbname(GLB_SET,imolty,iatom,' N')
         call twh_pdbname(GLB_SET,imolty,iatom+1,' CA')
         call twh_pdbname(GLB_SET,imolty,iatom+2,' C')
         call twh_pdbname(GLB_SET,imolty,iatom+3,' O')
         call twh_pdbname(GLB_SET,imolty,hmain+heavy,'1H')
         call twh_pdbname(GLB_SET,imolty,hmain+heavy+1,'2H')
         call twh_pdbname(GLB_SET,imolty,hmain+heavy+2,'3H')
         if ( aminoname .eq. 'glycine' ) then
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+3,'1HA')
         else
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+3,' HA')
         endif
      elseif ( termflag .eq. 0 ) then
c        --- interior segment main chain
c        --- set iha
         iha = hmain+heavy+1
c        --- vibrations
         call twh_invib(GLB_SET,imolty,iatom,3)
         call twh_ijvib(GLB_SET,imolty,iatom,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom,2,hmain+heavy)
         call twh_ijvib(GLB_SET,imolty,iatom,3,ibondend)
         call twh_invib(GLB_SET,imolty,iatom+1,4)
         call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
         call twh_ijvib(GLB_SET,imolty,iatom+1,2,iatom+2)
         if ( heavy .ne. 0 ) then
            call twh_ijvib(GLB_SET,imolty,iatom+1,3,hmain)
c           --- set icb
            icb = hmain
         else
            call twh_ijvib(GLB_SET,imolty,iatom+1,3,shift)
c           --- set icb
            icb = shift
         endif
         call twh_ijvib(GLB_SET,imolty,iatom+1,4,hmain+heavy+1)
         call twh_invib(GLB_SET,imolty,iatom+2,3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom+2,2,iatom+3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,3,iatom+total)
         call twh_invib(GLB_SET,imolty,iatom+3,1)
         call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom+2)
         call twh_invib(GLB_SET,imolty,hmain+heavy,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy,1,iatom)
         call twh_invib(GLB_SET,imolty,hmain+heavy+1,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy+1,1,iatom+1)
c        --- pdb names
         call twh_pdbname(GLB_SET,imolty,iatom,' N')
         call twh_pdbname(GLB_SET,imolty,iatom+1,' CA')
         call twh_pdbname(GLB_SET,imolty,iatom+2,' C')
         call twh_pdbname(GLB_SET,imolty,iatom+3,' O')
         call twh_pdbname(GLB_SET,imolty,hmain+heavy,' H')
         if ( aminoname .eq. 'glycine' ) then
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+1,'1HA')
         else
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+1,' HA')
         endif
      elseif ( termflag .eq. 1 ) then
c        --- c-terminus main chain
c        --- set iha
         iha = hmain+heavy+2
c        --- vibrations
         call twh_invib(GLB_SET,imolty,iatom,3)
         call twh_ijvib(GLB_SET,imolty,iatom,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom,2,hmain+heavy+1)
         call twh_ijvib(GLB_SET,imolty,iatom,3,ibondend)
         call twh_invib(GLB_SET,imolty,iatom+1,4)
         call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
         call twh_ijvib(GLB_SET,imolty,iatom+1,2,iatom+2)
         if ( heavy .ne. 0 ) then
            call twh_ijvib(GLB_SET,imolty,iatom+1,3,hmain)
c           --- set icb
            icb = hmain
         else
            call twh_ijvib(GLB_SET,imolty,iatom+1,3,shift)
c           --- set icb
            icb = shift
         endif
         call twh_ijvib(GLB_SET,imolty,iatom+1,4,hmain+heavy+2)
         call twh_invib(GLB_SET,imolty,iatom+2,3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom+2,2,iatom+3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,3,hmain+heavy)
         call twh_invib(GLB_SET,imolty,iatom+3,1)
         call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom+2)
         call twh_invib(GLB_SET,imolty,hmain+heavy,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy,1,iatom+2)
         call twh_invib(GLB_SET,imolty,hmain+heavy+1,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy+1,1,iatom)
         call twh_invib(GLB_SET,imolty,hmain+heavy+2,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy+2,1,iatom+1)
c        --- pdb names
         call twh_pdbname(GLB_SET,imolty,iatom,' N')
         call twh_pdbname(GLB_SET,imolty,iatom+1,' CA')
         call twh_pdbname(GLB_SET,imolty,iatom+2,' C')
         call twh_pdbname(GLB_SET,imolty,iatom+3,' O')
         call twh_pdbname(GLB_SET,imolty,hmain+heavy,' OXT')
         call twh_pdbname(GLB_SET,imolty,hmain+heavy+1,' H')
         if ( aminoname .eq. 'glycine' ) then
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+2,'1HA')
         else
            call twh_pdbname(GLB_SET,imolty,hmain+heavy+2,' HA')
         endif
      else
         write(6,*) 'AA_PEPNAME_MAIN_CHAIN unknown termflag:',termflag
         lfailure = .true.
         return
      endif

      return
      end

      subroutine twh_aa_pepname_a0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'a0' monomer in  *
c     * the polypeptide builder                                        *
c     * a0: alanine                                                    *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'alanine'
      shortname = 'ALA'
c     --- alanine side chain
c     --- all-atom geometry
      heavy = 1
      shift = hmain+heavy+ninsert
      light = 3
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
c     --- heavy atoms
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift+1)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+2)
c     --- hydrogens
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain)
c     --- assign pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')
      call twh_pdbname(GLB_SET,imolty,shift+2,'3HB')

      return
      end

      subroutine twh_aa_pepname_c0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'c0' monomer in  *
c     * the polypeptide builder                                        *
c     * c0: cysteine with H bonded to the S                            *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'cysteine with hydrogen'
      shortname = 'CYS'
c     --- cysteine side chain
      heavy = 2
      shift = hmain+heavy+ninsert
      light = 3
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,shift+2)
c     --- hydrogens
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+1)
c     --- pdb atom names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
c      call twh_pdbname(GLB_SET,imolty,hmain+1,' SG')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')
      call twh_pdbname(GLB_SET,imolty,shift+2,' HG')

      return
      end

      subroutine twh_aa_pepname_cp(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'cp' monomer in  *
c     * the polypeptide builder                                        *
c     * cp: cysteine bonded to palmitate                               *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polycharmm22 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'cysteine bonded to palmitate'
      shortname = 'CYP'
c     --- cysteine side chain
      heavy = 19
      shift = hmain+heavy+ninsert
      light = 33
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_invib(GLB_SET,imolty,hmain+2,3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,hmain+18)
      call twh_invib(GLB_SET,imolty,hmain+3,4)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,shift+2)
      call twh_ijvib(GLB_SET,imolty,hmain+3,4,shift+3)
      call twh_invib(GLB_SET,imolty,hmain+4,4)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,shift+4)
      call twh_ijvib(GLB_SET,imolty,hmain+4,4,shift+5)
      call twh_invib(GLB_SET,imolty,hmain+5,4)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+5,2,hmain+6)
      call twh_ijvib(GLB_SET,imolty,hmain+5,3,shift+6)
      call twh_ijvib(GLB_SET,imolty,hmain+5,4,shift+7)
      call twh_invib(GLB_SET,imolty,hmain+6,4)
      call twh_ijvib(GLB_SET,imolty,hmain+6,1,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+6,2,hmain+7)
      call twh_ijvib(GLB_SET,imolty,hmain+6,3,shift+8)
      call twh_ijvib(GLB_SET,imolty,hmain+6,4,shift+9)
      call twh_invib(GLB_SET,imolty,hmain+7,4)
      call twh_ijvib(GLB_SET,imolty,hmain+7,1,hmain+6)
      call twh_ijvib(GLB_SET,imolty,hmain+7,2,hmain+8)
      call twh_ijvib(GLB_SET,imolty,hmain+7,3,shift+10)
      call twh_ijvib(GLB_SET,imolty,hmain+7,4,shift+11)
      call twh_invib(GLB_SET,imolty,hmain+8,4)
      call twh_ijvib(GLB_SET,imolty,hmain+8,1,hmain+7)
      call twh_ijvib(GLB_SET,imolty,hmain+8,2,hmain+9)
      call twh_ijvib(GLB_SET,imolty,hmain+8,3,shift+12)
      call twh_ijvib(GLB_SET,imolty,hmain+8,4,shift+13)
      call twh_invib(GLB_SET,imolty,hmain+9,4)
      call twh_ijvib(GLB_SET,imolty,hmain+9,1,hmain+8)
      call twh_ijvib(GLB_SET,imolty,hmain+9,2,hmain+10)
      call twh_ijvib(GLB_SET,imolty,hmain+9,3,shift+14)
      call twh_ijvib(GLB_SET,imolty,hmain+9,4,shift+15)
      call twh_invib(GLB_SET,imolty,hmain+10,4)
      call twh_ijvib(GLB_SET,imolty,hmain+10,1,hmain+9)
      call twh_ijvib(GLB_SET,imolty,hmain+10,2,hmain+11)
      call twh_ijvib(GLB_SET,imolty,hmain+10,3,shift+16)
      call twh_ijvib(GLB_SET,imolty,hmain+10,4,shift+17)
      call twh_invib(GLB_SET,imolty,hmain+11,4)
      call twh_ijvib(GLB_SET,imolty,hmain+11,1,hmain+10)
      call twh_ijvib(GLB_SET,imolty,hmain+11,2,hmain+12)
      call twh_ijvib(GLB_SET,imolty,hmain+11,3,shift+18)
      call twh_ijvib(GLB_SET,imolty,hmain+11,4,shift+19)
      call twh_invib(GLB_SET,imolty,hmain+12,4)
      call twh_ijvib(GLB_SET,imolty,hmain+12,1,hmain+11)
      call twh_ijvib(GLB_SET,imolty,hmain+12,2,hmain+13)
      call twh_ijvib(GLB_SET,imolty,hmain+12,3,shift+20)
      call twh_ijvib(GLB_SET,imolty,hmain+12,4,shift+21)
      call twh_invib(GLB_SET,imolty,hmain+13,4)
      call twh_ijvib(GLB_SET,imolty,hmain+13,1,hmain+12)
      call twh_ijvib(GLB_SET,imolty,hmain+13,2,hmain+14)
      call twh_ijvib(GLB_SET,imolty,hmain+13,3,shift+22)
      call twh_ijvib(GLB_SET,imolty,hmain+13,4,shift+23)
      call twh_invib(GLB_SET,imolty,hmain+14,4)
      call twh_ijvib(GLB_SET,imolty,hmain+14,1,hmain+13)
      call twh_ijvib(GLB_SET,imolty,hmain+14,2,hmain+15)
      call twh_ijvib(GLB_SET,imolty,hmain+14,3,shift+24)
      call twh_ijvib(GLB_SET,imolty,hmain+14,4,shift+25)
      call twh_invib(GLB_SET,imolty,hmain+15,4)
      call twh_ijvib(GLB_SET,imolty,hmain+15,1,hmain+14)
      call twh_ijvib(GLB_SET,imolty,hmain+15,2,hmain+16)
      call twh_ijvib(GLB_SET,imolty,hmain+15,3,shift+26)
      call twh_ijvib(GLB_SET,imolty,hmain+15,4,shift+27)
      call twh_invib(GLB_SET,imolty,hmain+16,4)
      call twh_ijvib(GLB_SET,imolty,hmain+16,1,hmain+15)
      call twh_ijvib(GLB_SET,imolty,hmain+16,2,hmain+17)
      call twh_ijvib(GLB_SET,imolty,hmain+16,3,shift+28)
      call twh_ijvib(GLB_SET,imolty,hmain+16,4,shift+29)
      call twh_invib(GLB_SET,imolty,hmain+17,4)
      call twh_ijvib(GLB_SET,imolty,hmain+17,1,hmain+16)
      call twh_ijvib(GLB_SET,imolty,hmain+17,2,shift+30)
      call twh_ijvib(GLB_SET,imolty,hmain+17,3,shift+31)
      call twh_ijvib(GLB_SET,imolty,hmain+17,4,shift+32)
      call twh_invib(GLB_SET,imolty,hmain+18,1)
      call twh_ijvib(GLB_SET,imolty,hmain+18,1,hmain+2)
c     --- hydrogens
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+4,1)
      call twh_ijvib(GLB_SET,imolty,shift+4,1,hmain+4)
      call twh_invib(GLB_SET,imolty,shift+5,1)
      call twh_ijvib(GLB_SET,imolty,shift+5,1,hmain+4)
      call twh_invib(GLB_SET,imolty,shift+6,1)
      call twh_ijvib(GLB_SET,imolty,shift+6,1,hmain+5)
      call twh_invib(GLB_SET,imolty,shift+7,1)
      call twh_ijvib(GLB_SET,imolty,shift+7,1,hmain+5)
      call twh_invib(GLB_SET,imolty,shift+8,1)
      call twh_ijvib(GLB_SET,imolty,shift+8,1,hmain+6)
      call twh_invib(GLB_SET,imolty,shift+9,1)
      call twh_ijvib(GLB_SET,imolty,shift+9,1,hmain+6)
      call twh_invib(GLB_SET,imolty,shift+10,1)
      call twh_ijvib(GLB_SET,imolty,shift+10,1,hmain+7)
      call twh_invib(GLB_SET,imolty,shift+11,1)
      call twh_ijvib(GLB_SET,imolty,shift+11,1,hmain+7)
      call twh_invib(GLB_SET,imolty,shift+12,1)
      call twh_ijvib(GLB_SET,imolty,shift+12,1,hmain+8)
      call twh_invib(GLB_SET,imolty,shift+13,1)
      call twh_ijvib(GLB_SET,imolty,shift+13,1,hmain+8)
      call twh_invib(GLB_SET,imolty,shift+14,1)
      call twh_ijvib(GLB_SET,imolty,shift+14,1,hmain+9)
      call twh_invib(GLB_SET,imolty,shift+15,1)
      call twh_ijvib(GLB_SET,imolty,shift+15,1,hmain+9)
      call twh_invib(GLB_SET,imolty,shift+16,1)
      call twh_ijvib(GLB_SET,imolty,shift+16,1,hmain+10)
      call twh_invib(GLB_SET,imolty,shift+17,1)
      call twh_ijvib(GLB_SET,imolty,shift+17,1,hmain+10)
      call twh_invib(GLB_SET,imolty,shift+18,1)
      call twh_ijvib(GLB_SET,imolty,shift+18,1,hmain+11)
      call twh_invib(GLB_SET,imolty,shift+19,1)
      call twh_ijvib(GLB_SET,imolty,shift+19,1,hmain+11)
      call twh_invib(GLB_SET,imolty,shift+20,1)
      call twh_ijvib(GLB_SET,imolty,shift+20,1,hmain+12)
      call twh_invib(GLB_SET,imolty,shift+21,1)
      call twh_ijvib(GLB_SET,imolty,shift+21,1,hmain+12)
      call twh_invib(GLB_SET,imolty,shift+22,1)
      call twh_ijvib(GLB_SET,imolty,shift+22,1,hmain+13)
      call twh_invib(GLB_SET,imolty,shift+23,1)
      call twh_ijvib(GLB_SET,imolty,shift+23,1,hmain+13)
      call twh_invib(GLB_SET,imolty,shift+24,1)
      call twh_ijvib(GLB_SET,imolty,shift+24,1,hmain+14)
      call twh_invib(GLB_SET,imolty,shift+25,1)
      call twh_ijvib(GLB_SET,imolty,shift+25,1,hmain+14)
      call twh_invib(GLB_SET,imolty,shift+26,1)
      call twh_ijvib(GLB_SET,imolty,shift+26,1,hmain+15)
      call twh_invib(GLB_SET,imolty,shift+27,1)
      call twh_ijvib(GLB_SET,imolty,shift+27,1,hmain+15)
      call twh_invib(GLB_SET,imolty,shift+28,1)
      call twh_ijvib(GLB_SET,imolty,shift+28,1,hmain+16)
      call twh_invib(GLB_SET,imolty,shift+29,1)
      call twh_ijvib(GLB_SET,imolty,shift+29,1,hmain+16)
      call twh_invib(GLB_SET,imolty,shift+30,1)
      call twh_ijvib(GLB_SET,imolty,shift+30,1,hmain+17)
      call twh_invib(GLB_SET,imolty,shift+31,1)
      call twh_ijvib(GLB_SET,imolty,shift+31,1,hmain+17)
      call twh_invib(GLB_SET,imolty,shift+32,1)
      call twh_ijvib(GLB_SET,imolty,shift+32,1,hmain+17)
c     --- pdb atom names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' SG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' C1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' C2')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' C3')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' C4')
      call twh_pdbname(GLB_SET,imolty,hmain+6,' C5')
      call twh_pdbname(GLB_SET,imolty,hmain+7,' C6')
      call twh_pdbname(GLB_SET,imolty,hmain+8,' C7')
      call twh_pdbname(GLB_SET,imolty,hmain+9,' C8')
      call twh_pdbname(GLB_SET,imolty,hmain+10,' C9')
      call twh_pdbname(GLB_SET,imolty,hmain+11,' C10')
      call twh_pdbname(GLB_SET,imolty,hmain+12,' C11')
      call twh_pdbname(GLB_SET,imolty,hmain+13,' C12')
      call twh_pdbname(GLB_SET,imolty,hmain+14,' C13')
      call twh_pdbname(GLB_SET,imolty,hmain+15,' C14')
      call twh_pdbname(GLB_SET,imolty,hmain+16,' C15')
      call twh_pdbname(GLB_SET,imolty,hmain+17,' C16')
      call twh_pdbname(GLB_SET,imolty,hmain+18,' O1')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')
      call twh_pdbname(GLB_SET,imolty,shift+2,'1H2')
      call twh_pdbname(GLB_SET,imolty,shift+3,'2H2')
      call twh_pdbname(GLB_SET,imolty,shift+4,'1H3')
      call twh_pdbname(GLB_SET,imolty,shift+5,'2H3')
      call twh_pdbname(GLB_SET,imolty,shift+6,'1H4')
      call twh_pdbname(GLB_SET,imolty,shift+7,'2H4')
      call twh_pdbname(GLB_SET,imolty,shift+8,'1H5')
      call twh_pdbname(GLB_SET,imolty,shift+9,'2H5')
      call twh_pdbname(GLB_SET,imolty,shift+10,'1H6')
      call twh_pdbname(GLB_SET,imolty,shift+11,'2H6')
      call twh_pdbname(GLB_SET,imolty,shift+12,'1H7')
      call twh_pdbname(GLB_SET,imolty,shift+13,'2H7')
      call twh_pdbname(GLB_SET,imolty,shift+14,'1H8')
      call twh_pdbname(GLB_SET,imolty,shift+15,'2H8')
      call twh_pdbname(GLB_SET,imolty,shift+16,'1H9')
      call twh_pdbname(GLB_SET,imolty,shift+17,'2H9')
      call twh_pdbname(GLB_SET,imolty,shift+18,'1H10')
      call twh_pdbname(GLB_SET,imolty,shift+19,'2H10')
      call twh_pdbname(GLB_SET,imolty,shift+20,'1H11')
      call twh_pdbname(GLB_SET,imolty,shift+21,'2H11')
      call twh_pdbname(GLB_SET,imolty,shift+22,'1H12')
      call twh_pdbname(GLB_SET,imolty,shift+23,'2H12')
      call twh_pdbname(GLB_SET,imolty,shift+24,'1H13')
      call twh_pdbname(GLB_SET,imolty,shift+25,'2H13')
      call twh_pdbname(GLB_SET,imolty,shift+26,'1H14')
      call twh_pdbname(GLB_SET,imolty,shift+27,'2H14')
      call twh_pdbname(GLB_SET,imolty,shift+28,'1H15')
      call twh_pdbname(GLB_SET,imolty,shift+29,'2H15')
      call twh_pdbname(GLB_SET,imolty,shift+30,'1H16')
      call twh_pdbname(GLB_SET,imolty,shift+31,'2H16')
      call twh_pdbname(GLB_SET,imolty,shift+32,'3H16')

      return
      end


      subroutine twh_aa_pepname_cs(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'cs' monomer in  *
c     * the polypeptide builder                                        *
c     * cs: cysteine bonded to another cysteine                        *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'cysteine in disulfide bond'
      shortname = 'CYD'
c     --- cysteine side chain
      heavy = 2
      shift = hmain+heavy+ninsert
      light = 2
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,-1)
c     --- hydrogens
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
c     --- pdb atom names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' SG')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')

      return
      end

      subroutine twh_aa_pepname_d_neg(lfailure,aminoname,shortname
     &     ,imolty,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'd-' monomer in  *
c     * the polypeptide builder                                        *
c     * d-: aspartic acid deprotonated                                 *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'aspartic acid deprotonated'
      shortname = 'ASP'
c     --- asparagine side chain
      heavy = 4
      shift = hmain+heavy+ninsert
      light = 2
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+3,1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' OD1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' OD2')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')

      return
      end

      subroutine twh_aa_pepname_e_neg(lfailure,aminoname,shortname
     &     ,imolty,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'e-' monomer in  *
c     * the polypeptide builder                                        *
c     * e-: glutamic acid deprotonated                                 *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'glutamic acid deprotonated'
      shortname = 'GLU'
c     --- glutamic acid side chain
      heavy = 5
      shift = hmain+heavy+ninsert
      light = 4
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,4)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,shift+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,4,shift+3)
      call twh_invib(GLB_SET,imolty,hmain+2,3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,hmain+4)
      call twh_invib(GLB_SET,imolty,hmain+3,1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+2)
      call twh_invib(GLB_SET,imolty,hmain+4,1)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+1)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CD')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' OE1')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' OE2')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')
      call twh_pdbname(GLB_SET,imolty,shift+2,'1HG')
      call twh_pdbname(GLB_SET,imolty,shift+3,'2HG')

      return
      end

      subroutine twh_aa_pepname_f0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'f0' monomer in  *
c     * the polypeptide builder                                        *
c     * f0: phenylalanine                                              *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'phenylalanine'
      shortname = 'PHE'
c     --- phenylalanine side chain
      heavy = 7
      shift = hmain+heavy+ninsert
      light = 7
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,shift+2)
      call twh_invib(GLB_SET,imolty,hmain+3,3)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,shift+3)
      call twh_invib(GLB_SET,imolty,hmain+4,3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+6)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,shift+4)
      call twh_invib(GLB_SET,imolty,hmain+5,3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,2,hmain+6)
      call twh_ijvib(GLB_SET,imolty,hmain+5,3,shift+5)
      call twh_invib(GLB_SET,imolty,hmain+6,3)
      call twh_ijvib(GLB_SET,imolty,hmain+6,1,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+6,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+6,3,shift+6)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+4,1)
      call twh_ijvib(GLB_SET,imolty,shift+4,1,hmain+4)
      call twh_invib(GLB_SET,imolty,shift+5,1)
      call twh_ijvib(GLB_SET,imolty,shift+5,1,hmain+5)
      call twh_invib(GLB_SET,imolty,shift+6,1)
      call twh_ijvib(GLB_SET,imolty,shift+6,1,hmain+6)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CD1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CD2')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' CE1')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' CE2')
      call twh_pdbname(GLB_SET,imolty,hmain+6,' CZ')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')
      call twh_pdbname(GLB_SET,imolty,shift+2,' HD1')
      call twh_pdbname(GLB_SET,imolty,shift+3,' HD2')
      call twh_pdbname(GLB_SET,imolty,shift+4,' HE1')
      call twh_pdbname(GLB_SET,imolty,shift+5,' HE2')
      call twh_pdbname(GLB_SET,imolty,shift+6,' HZ')

      return
      end

      subroutine twh_aa_pepname_g0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'g0' monomer in  *
c     * the polypeptide builder                                        *
c     * g0: glycine monomer                                            *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'glycine'
      shortname = 'GLY'
c     --- side chain glycine
c     --- vibrations
      heavy = 0
      shift = hmain+heavy+ninsert
      light = 1
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,iatom+1)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,shift,'2HA')

      return
      end

      subroutine twh_aa_pepname_hd(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'hd' monomer in  *
c     * the polypeptide builder                                        *
c     * hd: proton on delta histidine                                  *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert
      
      aminoname = 'histidine neutral only Nd protonated'
      shortname = 'HIS'
c     --- histidine side chain
      heavy = 6
      shift = hmain+heavy+ninsert
      light = 5
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,shift+2)
      call twh_invib(GLB_SET,imolty,hmain+3,3)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,shift+3)
      call twh_invib(GLB_SET,imolty,hmain+4,3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,shift+4)
      call twh_invib(GLB_SET,imolty,hmain+5,2)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,2,hmain+4)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+4,1)
      call twh_ijvib(GLB_SET,imolty,shift+4,1,hmain+4)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' ND1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CD2')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' CE1')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' NE2')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')
      call twh_pdbname(GLB_SET,imolty,shift+2,' HD1')
      call twh_pdbname(GLB_SET,imolty,shift+3,' HD2')
      call twh_pdbname(GLB_SET,imolty,shift+4,' HE1')

      return
      end

      subroutine twh_aa_pepname_he(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'he' monomer in  *
c     * the polypeptide builder                                        *
c     * he: proton on epsilon histidine                                *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert
      
      aminoname = 'histidine neutral only Ne protonated'
      shortname = 'HIS'
c     --- histidine side chain
      heavy = 6
      shift = hmain+heavy+ninsert
      light = 5
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,2)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+4)
      call twh_invib(GLB_SET,imolty,hmain+3,3)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,shift+2)
      call twh_invib(GLB_SET,imolty,hmain+4,3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,shift+3)
      call twh_invib(GLB_SET,imolty,hmain+5,3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+5,3,shift+4)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+4)
      call twh_invib(GLB_SET,imolty,shift+4,1)
      call twh_ijvib(GLB_SET,imolty,shift+4,1,hmain+5)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' ND1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CD2')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' CE1')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' NE2')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')
      call twh_pdbname(GLB_SET,imolty,shift+2,' HD2')
      call twh_pdbname(GLB_SET,imolty,shift+3,' HE1')
      call twh_pdbname(GLB_SET,imolty,shift+4,' HE2')

      return
      end

      subroutine twh_aa_pepname_h_plus(lfailure,aminoname,shortname
     &     ,imolty,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'h+' monomer in  *
c     * the polypeptide builder                                        *
c     * h+: doubly protonated histidine                                *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'histidine both N protonated'
      shortname = 'HIS'
c     --- histidine side chain
      heavy = 6
      shift = hmain+heavy+ninsert
      light = 6
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,shift+2)
      call twh_invib(GLB_SET,imolty,hmain+3,3)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,shift+3)
      call twh_invib(GLB_SET,imolty,hmain+4,3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,shift+4)
      call twh_invib(GLB_SET,imolty,hmain+5,3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+5,3,shift+5)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+4,1)
      call twh_ijvib(GLB_SET,imolty,shift+4,1,hmain+4)
      call twh_invib(GLB_SET,imolty,shift+5,1)
      call twh_ijvib(GLB_SET,imolty,shift+5,1,hmain+5)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' ND1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CD2')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' CE1')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' NE2')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')
      call twh_pdbname(GLB_SET,imolty,shift+2,' HD1')
      call twh_pdbname(GLB_SET,imolty,shift+3,' HD2')
      call twh_pdbname(GLB_SET,imolty,shift+4,' HE1')
      call twh_pdbname(GLB_SET,imolty,shift+5,' HE2')

      return
      end

      subroutine twh_aa_pepname_i0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'i0' monomer in  *
c     * the polypeptide builder                                        *
c     * i0: isoleucine monomer                                         *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'isoleucine'
      shortname = 'ILE'
c     --- side chain isoleucine
      heavy = 4
      shift = hmain+heavy+ninsert
      light = 9
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift)
      call twh_invib(GLB_SET,imolty,hmain+1,4)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,shift+1)
      call twh_ijvib(GLB_SET,imolty,hmain+1,4,shift+2)
      call twh_invib(GLB_SET,imolty,hmain+2,4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,shift+3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,shift+4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,4,shift+5)
      call twh_invib(GLB_SET,imolty,hmain+3,4)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,shift+6)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,shift+7)
      call twh_ijvib(GLB_SET,imolty,hmain+3,4,shift+8)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+4,1)
      call twh_ijvib(GLB_SET,imolty,shift+4,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+5,1)
      call twh_ijvib(GLB_SET,imolty,shift+5,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+6,1)
      call twh_ijvib(GLB_SET,imolty,shift+6,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+7,1)
      call twh_ijvib(GLB_SET,imolty,shift+7,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+8,1)
      call twh_ijvib(GLB_SET,imolty,shift+8,1,hmain+3)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG1')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CG2')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CD1')
      call twh_pdbname(GLB_SET,imolty,shift,' HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'1HG1')
      call twh_pdbname(GLB_SET,imolty,shift+2,'2HG1')
      call twh_pdbname(GLB_SET,imolty,shift+3,'1HG2')
      call twh_pdbname(GLB_SET,imolty,shift+4,'2HG2')
      call twh_pdbname(GLB_SET,imolty,shift+5,'3HG2')
      call twh_pdbname(GLB_SET,imolty,shift+6,'1HD')
      call twh_pdbname(GLB_SET,imolty,shift+7,'2HD')
      call twh_pdbname(GLB_SET,imolty,shift+8,'3HD')

      return
      end

      subroutine twh_aa_pepname_k_plus(lfailure,aminoname,shortname
     &     ,imolty,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'k+' monomer in  *
c     * the polypeptide builder                                        *
c     * k+: lysine protonated                                          *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'lysine protonated'
      shortname = 'LYS'
c     --- lysine side chain
      heavy = 5
      shift = hmain+heavy+ninsert
      light = 11
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,4)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,shift+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,4,shift+3)
      call twh_invib(GLB_SET,imolty,hmain+2,4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,shift+4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,4,shift+5)
      call twh_invib(GLB_SET,imolty,hmain+3,4)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,shift+6)
      call twh_ijvib(GLB_SET,imolty,hmain+3,4,shift+7)
      call twh_invib(GLB_SET,imolty,hmain+4,4)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,shift+8)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,shift+9)
      call twh_ijvib(GLB_SET,imolty,hmain+4,4,shift+10)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+4,1)
      call twh_ijvib(GLB_SET,imolty,shift+4,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+5,1)
      call twh_ijvib(GLB_SET,imolty,shift+5,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+6,1)
      call twh_ijvib(GLB_SET,imolty,shift+6,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+7,1)
      call twh_ijvib(GLB_SET,imolty,shift+7,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+8,1)
      call twh_ijvib(GLB_SET,imolty,shift+8,1,hmain+4)
      call twh_invib(GLB_SET,imolty,shift+9,1)
      call twh_ijvib(GLB_SET,imolty,shift+9,1,hmain+4)
      call twh_invib(GLB_SET,imolty,shift+10,1)
      call twh_ijvib(GLB_SET,imolty,shift+10,1,hmain+4)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CD')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CE')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' NZ')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')
      call twh_pdbname(GLB_SET,imolty,shift+2,'1HG')
      call twh_pdbname(GLB_SET,imolty,shift+3,'2HG')
      call twh_pdbname(GLB_SET,imolty,shift+4,'1HD')
      call twh_pdbname(GLB_SET,imolty,shift+5,'2HD')
      call twh_pdbname(GLB_SET,imolty,shift+6,'1HE')
      call twh_pdbname(GLB_SET,imolty,shift+7,'2HE')
      call twh_pdbname(GLB_SET,imolty,shift+8,'1HZ')
      call twh_pdbname(GLB_SET,imolty,shift+9,'2HZ')
      call twh_pdbname(GLB_SET,imolty,shift+10,'3HZ')

      return
      end

      subroutine twh_aa_pepname_l0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'l0' monomer in  *
c     * the polypeptide builder                                        *
c     * l0: leucine monomer                                            *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'leucine'
      shortname = 'LEU'
c     --- side chain leucine
      heavy = 4
      shift = hmain+heavy+ninsert
      light = 9
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,4)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,4,shift+2)
      call twh_invib(GLB_SET,imolty,hmain+2,4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,shift+3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,shift+4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,4,shift+5)
      call twh_invib(GLB_SET,imolty,hmain+3,4)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,shift+6)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,shift+7)
      call twh_ijvib(GLB_SET,imolty,hmain+3,4,shift+8)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+4,1)
      call twh_ijvib(GLB_SET,imolty,shift+4,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+5,1)
      call twh_ijvib(GLB_SET,imolty,shift+5,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+6,1)
      call twh_ijvib(GLB_SET,imolty,shift+6,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+7,1)
      call twh_ijvib(GLB_SET,imolty,shift+7,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+8,1)
      call twh_ijvib(GLB_SET,imolty,shift+8,1,hmain+3)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CD1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CD2')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')
      call twh_pdbname(GLB_SET,imolty,shift+2,' HG')
      call twh_pdbname(GLB_SET,imolty,shift+3,'1HD1')
      call twh_pdbname(GLB_SET,imolty,shift+4,'2HD1')
      call twh_pdbname(GLB_SET,imolty,shift+5,'3HD1')
      call twh_pdbname(GLB_SET,imolty,shift+6,'1HD2')
      call twh_pdbname(GLB_SET,imolty,shift+7,'2HD2')
      call twh_pdbname(GLB_SET,imolty,shift+8,'3HD2')

      return
      end

      subroutine twh_aa_pepname_m0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'm0' monomer in  *
c     * the polypeptide builder                                        *
c     * m0: methionine monomer                                         *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'methionine'
      shortname = 'MET'
c     --- methionine side chain
      heavy = 4
      shift = hmain+heavy+ninsert
      light = 7
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,4)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,shift+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,4,shift+3)
      call twh_invib(GLB_SET,imolty,hmain+2,2)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+3,4)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,shift+4)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,shift+5)
      call twh_ijvib(GLB_SET,imolty,hmain+3,4,shift+6)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+4,1)
      call twh_ijvib(GLB_SET,imolty,shift+4,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+5,1)
      call twh_ijvib(GLB_SET,imolty,shift+5,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+6,1)
      call twh_ijvib(GLB_SET,imolty,shift+6,1,hmain+3)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' SD')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CE')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')
      call twh_pdbname(GLB_SET,imolty,shift+2,'1HG')
      call twh_pdbname(GLB_SET,imolty,shift+3,'2HG')
      call twh_pdbname(GLB_SET,imolty,shift+4,'1HE')
      call twh_pdbname(GLB_SET,imolty,shift+5,'2HE')
      call twh_pdbname(GLB_SET,imolty,shift+6,'3HE')

      return
      end

      subroutine twh_aa_pepname_n0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'n0' monomer in  *
c     * the polypeptide builder                                        *
c     * n0: asparagine                                                 *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'asparagine'
      shortname = 'ASN'
c     --- asparagine side chain
c     --- set up atom types and vibrations for side chain
      heavy = 4
      shift = hmain+heavy+ninsert
      light = 4
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
c     --- heavy atoms
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+3,3)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,shift+2)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,shift+3)
c     --- hydrogens
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+3)
c     --- set pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' OD1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' ND2')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')
      call twh_pdbname(GLB_SET,imolty,shift+2,'1HD2')
      call twh_pdbname(GLB_SET,imolty,shift+3,'2HD2')

      return
      end

      subroutine twh_aa_pepname_p0(lfailure,aminoname,shortname,termflag
     &     ,imolty,total,iatom,icb,iha,ibondend)
c     ******************************************************************
c     * sets up the communal bond information for the 'p0' monomer in  *
c     * the polypeptide builder                                        *
c     * p0: proline                                                    *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer termflag,imolty,total,iatom,icb,iha,ibondend
c     --- local variables

c     --- proline requires its own set of special cases
      aminoname = 'proline'
      shortname = 'PRO'
      if ( termflag .eq. -1 ) then
c        --- N-terminal proline
c        --- set icb
         icb = iatom+4
c        --- set iha
         iha = iatom+9
         total = 16
         call twh_polysafetycheck(lfailure,iatom,total)
         if ( lfailure ) return
c        --- vibrations
         call twh_invib(GLB_SET,imolty,iatom,4)
         call twh_ijvib(GLB_SET,imolty,iatom,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom,2,iatom+6)
         call twh_ijvib(GLB_SET,imolty,iatom,3,iatom+7)
         call twh_ijvib(GLB_SET,imolty,iatom,4,iatom+8)
         call twh_invib(GLB_SET,imolty,iatom+1,4)
         call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
         call twh_ijvib(GLB_SET,imolty,iatom+1,2,iatom+2)
         call twh_ijvib(GLB_SET,imolty,iatom+1,3,iatom+4)
         call twh_ijvib(GLB_SET,imolty,iatom+1,4,iatom+9)
         call twh_invib(GLB_SET,imolty,iatom+2,3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom+2,2,iatom+3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,3,iatom+total)
         call twh_invib(GLB_SET,imolty,iatom+3,1)
         call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom+2)
         call twh_invib(GLB_SET,imolty,iatom+4,4)
         call twh_ijvib(GLB_SET,imolty,iatom+4,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom+4,2,iatom+5)
         call twh_ijvib(GLB_SET,imolty,iatom+4,3,iatom+10)
         call twh_ijvib(GLB_SET,imolty,iatom+4,4,iatom+11)
         call twh_invib(GLB_SET,imolty,iatom+5,4)
         call twh_ijvib(GLB_SET,imolty,iatom+5,1,iatom+4)
         call twh_ijvib(GLB_SET,imolty,iatom+5,2,iatom+6)
         call twh_ijvib(GLB_SET,imolty,iatom+5,3,iatom+12)
         call twh_ijvib(GLB_SET,imolty,iatom+5,4,iatom+13)
         call twh_invib(GLB_SET,imolty,iatom+6,4)
         call twh_ijvib(GLB_SET,imolty,iatom+6,1,iatom)
         call twh_ijvib(GLB_SET,imolty,iatom+6,2,iatom+5)
         call twh_ijvib(GLB_SET,imolty,iatom+6,3,iatom+14)
         call twh_ijvib(GLB_SET,imolty,iatom+6,4,iatom+15)
         call twh_invib(GLB_SET,imolty,iatom+7,1)
         call twh_ijvib(GLB_SET,imolty,iatom+7,1,iatom)
         call twh_invib(GLB_SET,imolty,iatom+8,1)
         call twh_ijvib(GLB_SET,imolty,iatom+8,1,iatom)
         call twh_invib(GLB_SET,imolty,iatom+9,1)
         call twh_ijvib(GLB_SET,imolty,iatom+9,1,iatom+1)
         call twh_invib(GLB_SET,imolty,iatom+10,1)
         call twh_ijvib(GLB_SET,imolty,iatom+10,1,iatom+4)
         call twh_invib(GLB_SET,imolty,iatom+11,1)
         call twh_ijvib(GLB_SET,imolty,iatom+11,1,iatom+4)
         call twh_invib(GLB_SET,imolty,iatom+12,1)
         call twh_ijvib(GLB_SET,imolty,iatom+12,1,iatom+5)
         call twh_invib(GLB_SET,imolty,iatom+13,1)
         call twh_ijvib(GLB_SET,imolty,iatom+13,1,iatom+5)
         call twh_invib(GLB_SET,imolty,iatom+14,1)
         call twh_ijvib(GLB_SET,imolty,iatom+14,1,iatom+6)
         call twh_invib(GLB_SET,imolty,iatom+15,1)
         call twh_ijvib(GLB_SET,imolty,iatom+15,1,iatom+6)
c        --- pdb names
         call twh_pdbname(GLB_SET,imolty,iatom,' N')
         call twh_pdbname(GLB_SET,imolty,iatom+1,' CA')
         call twh_pdbname(GLB_SET,imolty,iatom+2,' C')
         call twh_pdbname(GLB_SET,imolty,iatom+3,' O')
         call twh_pdbname(GLB_SET,imolty,iatom+4,' CB')
         call twh_pdbname(GLB_SET,imolty,iatom+5,' CG')
         call twh_pdbname(GLB_SET,imolty,iatom+6,' CD')
         call twh_pdbname(GLB_SET,imolty,iatom+7,'1H')
         call twh_pdbname(GLB_SET,imolty,iatom+8,'2H')
         call twh_pdbname(GLB_SET,imolty,iatom+9,' HA')
         call twh_pdbname(GLB_SET,imolty,iatom+10,'1HB')
         call twh_pdbname(GLB_SET,imolty,iatom+11,'2HB')
         call twh_pdbname(GLB_SET,imolty,iatom+12,'1HG')
         call twh_pdbname(GLB_SET,imolty,iatom+13,'2HG')
         call twh_pdbname(GLB_SET,imolty,iatom+14,'1HD')
         call twh_pdbname(GLB_SET,imolty,iatom+15,'2HD')
      elseif ( termflag .eq. 0 ) then
c        --- interior proline
c        --- set icb
         icb = iatom+4
c        --- set iha
         iha = iatom+7
         total = 14
         call twh_polysafetycheck(lfailure,iatom,total)
         if ( lfailure ) return
c        --- vibrations
         call twh_invib(GLB_SET,imolty,iatom,3)
         call twh_ijvib(GLB_SET,imolty,iatom,1,ibondend)
         call twh_ijvib(GLB_SET,imolty,iatom,2,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom,3,iatom+6)
         call twh_invib(GLB_SET,imolty,iatom+1,4)
         call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
         call twh_ijvib(GLB_SET,imolty,iatom+1,2,iatom+2)
         call twh_ijvib(GLB_SET,imolty,iatom+1,3,iatom+4)
         call twh_ijvib(GLB_SET,imolty,iatom+1,4,iatom+7)
         call twh_invib(GLB_SET,imolty,iatom+2,3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom+2,2,iatom+3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,3,iatom+total)
         call twh_invib(GLB_SET,imolty,iatom+3,1)
         call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom+2)
         call twh_invib(GLB_SET,imolty,iatom+4,4)
         call twh_ijvib(GLB_SET,imolty,iatom+4,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom+4,2,iatom+5)
         call twh_ijvib(GLB_SET,imolty,iatom+4,3,iatom+8)
         call twh_ijvib(GLB_SET,imolty,iatom+4,4,iatom+9)
         call twh_invib(GLB_SET,imolty,iatom+5,4)
         call twh_ijvib(GLB_SET,imolty,iatom+5,1,iatom+4)
         call twh_ijvib(GLB_SET,imolty,iatom+5,2,iatom+6)
         call twh_ijvib(GLB_SET,imolty,iatom+5,3,iatom+10)
         call twh_ijvib(GLB_SET,imolty,iatom+5,4,iatom+11)
         call twh_invib(GLB_SET,imolty,iatom+6,4)
         call twh_ijvib(GLB_SET,imolty,iatom+6,1,iatom)
         call twh_ijvib(GLB_SET,imolty,iatom+6,2,iatom+5)
         call twh_ijvib(GLB_SET,imolty,iatom+6,3,iatom+12)
         call twh_ijvib(GLB_SET,imolty,iatom+6,4,iatom+13)
         call twh_invib(GLB_SET,imolty,iatom+7,1)
         call twh_ijvib(GLB_SET,imolty,iatom+7,1,iatom+1)
         call twh_invib(GLB_SET,imolty,iatom+8,1)
         call twh_ijvib(GLB_SET,imolty,iatom+8,1,iatom+4)
         call twh_invib(GLB_SET,imolty,iatom+9,1)
         call twh_ijvib(GLB_SET,imolty,iatom+9,1,iatom+4)
         call twh_invib(GLB_SET,imolty,iatom+10,1)
         call twh_ijvib(GLB_SET,imolty,iatom+10,1,iatom+5)
         call twh_invib(GLB_SET,imolty,iatom+11,1)
         call twh_ijvib(GLB_SET,imolty,iatom+11,1,iatom+5)
         call twh_invib(GLB_SET,imolty,iatom+12,1)
         call twh_ijvib(GLB_SET,imolty,iatom+12,1,iatom+6)
         call twh_invib(GLB_SET,imolty,iatom+13,1)
         call twh_ijvib(GLB_SET,imolty,iatom+13,1,iatom+6)
c        --- pdb names
         call twh_pdbname(GLB_SET,imolty,iatom,' N')
         call twh_pdbname(GLB_SET,imolty,iatom+1,' CA')
         call twh_pdbname(GLB_SET,imolty,iatom+2,' C')
         call twh_pdbname(GLB_SET,imolty,iatom+3,' O')
         call twh_pdbname(GLB_SET,imolty,iatom+4,' CB')
         call twh_pdbname(GLB_SET,imolty,iatom+5,' CG')
         call twh_pdbname(GLB_SET,imolty,iatom+6,' CD')
         call twh_pdbname(GLB_SET,imolty,iatom+7,' HA')
         call twh_pdbname(GLB_SET,imolty,iatom+8,'1HB')
         call twh_pdbname(GLB_SET,imolty,iatom+9,'2HB')
         call twh_pdbname(GLB_SET,imolty,iatom+10,'1HG')
         call twh_pdbname(GLB_SET,imolty,iatom+11,'2HG')
         call twh_pdbname(GLB_SET,imolty,iatom+12,'1HD')
         call twh_pdbname(GLB_SET,imolty,iatom+13,'2HD')
      elseif ( termflag .eq. 1 ) then
c        --- C-terminal proline
c        --- set icb
         icb = iatom+4
c        --- set iha
         iha = iatom+8
         total = 15
         call twh_polysafetycheck(lfailure,iatom,total)
         if ( lfailure ) return
c        --- vibrations
         call twh_invib(GLB_SET,imolty,iatom,3)
         call twh_ijvib(GLB_SET,imolty,iatom,1,ibondend)
         call twh_ijvib(GLB_SET,imolty,iatom,2,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom,3,iatom+6)
         call twh_invib(GLB_SET,imolty,iatom+1,4)
         call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
         call twh_ijvib(GLB_SET,imolty,iatom+1,2,iatom+2)
         call twh_ijvib(GLB_SET,imolty,iatom+1,3,iatom+4)
         call twh_ijvib(GLB_SET,imolty,iatom+1,4,iatom+8)
         call twh_invib(GLB_SET,imolty,iatom+2,3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom+2,2,iatom+3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,3,iatom+7)
         call twh_invib(GLB_SET,imolty,iatom+3,1)
         call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom+2)
         call twh_invib(GLB_SET,imolty,iatom+4,4)
         call twh_ijvib(GLB_SET,imolty,iatom+4,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom+4,2,iatom+5)
         call twh_ijvib(GLB_SET,imolty,iatom+4,3,iatom+9)
         call twh_ijvib(GLB_SET,imolty,iatom+4,4,iatom+10)
         call twh_invib(GLB_SET,imolty,iatom+5,4)
         call twh_ijvib(GLB_SET,imolty,iatom+5,1,iatom+4)
         call twh_ijvib(GLB_SET,imolty,iatom+5,2,iatom+6)
         call twh_ijvib(GLB_SET,imolty,iatom+5,3,iatom+11)
         call twh_ijvib(GLB_SET,imolty,iatom+5,4,iatom+12)
         call twh_invib(GLB_SET,imolty,iatom+6,4)
         call twh_ijvib(GLB_SET,imolty,iatom+6,1,iatom)
         call twh_ijvib(GLB_SET,imolty,iatom+6,2,iatom+5)
         call twh_ijvib(GLB_SET,imolty,iatom+6,3,iatom+13)
         call twh_ijvib(GLB_SET,imolty,iatom+6,4,iatom+14)
         call twh_invib(GLB_SET,imolty,iatom+7,1)
         call twh_ijvib(GLB_SET,imolty,iatom+7,1,iatom+2)
         call twh_invib(GLB_SET,imolty,iatom+8,1)
         call twh_ijvib(GLB_SET,imolty,iatom+8,1,iatom+1)
         call twh_invib(GLB_SET,imolty,iatom+9,1)
         call twh_ijvib(GLB_SET,imolty,iatom+9,1,iatom+4)
         call twh_invib(GLB_SET,imolty,iatom+10,1)
         call twh_ijvib(GLB_SET,imolty,iatom+10,1,iatom+4)
         call twh_invib(GLB_SET,imolty,iatom+11,1)
         call twh_ijvib(GLB_SET,imolty,iatom+11,1,iatom+5)
         call twh_invib(GLB_SET,imolty,iatom+12,1)
         call twh_ijvib(GLB_SET,imolty,iatom+12,1,iatom+5)
         call twh_invib(GLB_SET,imolty,iatom+13,1)
         call twh_ijvib(GLB_SET,imolty,iatom+13,1,iatom+6)
         call twh_invib(GLB_SET,imolty,iatom+14,1)
         call twh_ijvib(GLB_SET,imolty,iatom+14,1,iatom+6)
c        --- pdb names
         call twh_pdbname(GLB_SET,imolty,iatom,' N')
         call twh_pdbname(GLB_SET,imolty,iatom+1,' CA')
         call twh_pdbname(GLB_SET,imolty,iatom+2,' C')
         call twh_pdbname(GLB_SET,imolty,iatom+3,' O')
         call twh_pdbname(GLB_SET,imolty,iatom+4,' CB')
         call twh_pdbname(GLB_SET,imolty,iatom+5,' CG')
         call twh_pdbname(GLB_SET,imolty,iatom+6,' CD')
         call twh_pdbname(GLB_SET,imolty,iatom+7,' OXT')
         call twh_pdbname(GLB_SET,imolty,iatom+8,' HA')
         call twh_pdbname(GLB_SET,imolty,iatom+9,'1HB')
         call twh_pdbname(GLB_SET,imolty,iatom+10,'2HB')
         call twh_pdbname(GLB_SET,imolty,iatom+11,'1HG')
         call twh_pdbname(GLB_SET,imolty,iatom+12,'2HG')
         call twh_pdbname(GLB_SET,imolty,iatom+13,'1HD')
         call twh_pdbname(GLB_SET,imolty,iatom+14,'2HD')
      else
         write(6,*) 'AA_PEPNAME_P0: invalid termflag:',termflag
         lfailure = .true.
      endif
  
      return
      end

      subroutine twh_aa_pepname_q0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'q0' monomer in  *
c     * the polypeptide builder                                        *
c     * q0: glutamine monomer                                          *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'glutamine'
      shortname = 'GLN'
c     --- glutamine side chain
      heavy = 5
      shift = hmain+heavy+ninsert
      light = 6
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,4)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,shift+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,4,shift+3)
      call twh_invib(GLB_SET,imolty,hmain+2,3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,hmain+4)
      call twh_invib(GLB_SET,imolty,hmain+3,1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+2)
      call twh_invib(GLB_SET,imolty,hmain+4,3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,shift+4)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,shift+5)
c     --- hydrogens
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+4,1)
      call twh_ijvib(GLB_SET,imolty,shift+4,1,hmain+4)
      call twh_invib(GLB_SET,imolty,shift+5,1)
      call twh_ijvib(GLB_SET,imolty,shift+5,1,hmain+4)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CD')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' OE1')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' NE2')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')
      call twh_pdbname(GLB_SET,imolty,shift+2,'1HG')
      call twh_pdbname(GLB_SET,imolty,shift+3,'2HG')
      call twh_pdbname(GLB_SET,imolty,shift+4,'1HE2')
      call twh_pdbname(GLB_SET,imolty,shift+5,'2HE2')

      return
      end

      subroutine twh_aa_pepname_r_plus(lfailure,aminoname,shortname
     &     ,imolty,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'r+' monomer in  *
c     * the polypeptide builder                                        *
c     * r+: protonated arginine                                        *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'arginine protonated'
      shortname = 'ARG'
c     --- arginine side chain
      heavy = 7
      shift = hmain+heavy+ninsert
      light = 11
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,4)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,shift+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,4,shift+3)
      call twh_invib(GLB_SET,imolty,hmain+2,4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,shift+4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,4,shift+5)
      call twh_invib(GLB_SET,imolty,hmain+3,3)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,shift+6)
      call twh_invib(GLB_SET,imolty,hmain+4,3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,hmain+6)
      call twh_invib(GLB_SET,imolty,hmain+5,3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+5,2,shift+7)
      call twh_ijvib(GLB_SET,imolty,hmain+5,3,shift+8)
      call twh_invib(GLB_SET,imolty,hmain+6,3)
      call twh_ijvib(GLB_SET,imolty,hmain+6,1,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+6,2,shift+9)
      call twh_ijvib(GLB_SET,imolty,hmain+6,3,shift+10)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+4,1)
      call twh_ijvib(GLB_SET,imolty,shift+4,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+5,1)
      call twh_ijvib(GLB_SET,imolty,shift+5,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+6,1)
      call twh_ijvib(GLB_SET,imolty,shift+6,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+7,1)
      call twh_ijvib(GLB_SET,imolty,shift+7,1,hmain+5)
      call twh_invib(GLB_SET,imolty,shift+8,1)
      call twh_ijvib(GLB_SET,imolty,shift+8,1,hmain+5)
      call twh_invib(GLB_SET,imolty,shift+9,1)
      call twh_ijvib(GLB_SET,imolty,shift+9,1,hmain+6)
      call twh_invib(GLB_SET,imolty,shift+10,1)
      call twh_ijvib(GLB_SET,imolty,shift+10,1,hmain+6)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CD')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' NE')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' CZ')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' NH1')
      call twh_pdbname(GLB_SET,imolty,hmain+6,' NH2')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')
      call twh_pdbname(GLB_SET,imolty,shift+2,'1HG')
      call twh_pdbname(GLB_SET,imolty,shift+3,'2HG')
      call twh_pdbname(GLB_SET,imolty,shift+4,'1HD')
      call twh_pdbname(GLB_SET,imolty,shift+5,'2HD')
      call twh_pdbname(GLB_SET,imolty,shift+6,' HE')
      call twh_pdbname(GLB_SET,imolty,shift+7,'1HH1')
      call twh_pdbname(GLB_SET,imolty,shift+8,'2HH1')
      call twh_pdbname(GLB_SET,imolty,shift+9,'1HH2')
      call twh_pdbname(GLB_SET,imolty,shift+10,'2HH2')

      return
      end

      subroutine twh_aa_pepname_s0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 's0' monomer in  *
c     * the polypeptide builder                                        *
c     * s0: serine                                                     *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'serine'
      shortname = 'SER'
c     --- serine side chain
      heavy = 2
      shift = hmain+heavy+ninsert
      light = 3
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,shift+2)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+1)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' OG')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')
      call twh_pdbname(GLB_SET,imolty,shift+2,' HG')

      return
      end

      subroutine twh_aa_pepname_t0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 't0' monomer in  *
c     * the polypeptide builder                                        *
c     * t0: threonine                                                  *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'threonine'
      shortname = 'THR'
c     --- threonine side chain
      heavy = 3
      shift = hmain+heavy+ninsert
      light = 5
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift)
      call twh_invib(GLB_SET,imolty,hmain+1,2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+2,4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,shift+2)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,shift+3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,4,shift+4)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+4,1)
      call twh_ijvib(GLB_SET,imolty,shift+4,1,hmain+2)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' OG1')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CG2')
      call twh_pdbname(GLB_SET,imolty,shift,' HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,' HG1')
      call twh_pdbname(GLB_SET,imolty,shift+2,'1HG2')
      call twh_pdbname(GLB_SET,imolty,shift+3,'2HG2')
      call twh_pdbname(GLB_SET,imolty,shift+4,'3HG2')

      return
      end

      subroutine twh_aa_pepname_v0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'v0' monomer in  *
c     * the polypeptide builder                                        *
c     * v0: valine                                                     *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'valine'
      shortname = 'VAL'
c     --- valine side chain
      heavy = 3
      shift = hmain+heavy+ninsert
      light = 7
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift)
      call twh_invib(GLB_SET,imolty,hmain+1,4)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,shift+1)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,shift+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,4,shift+3)
      call twh_invib(GLB_SET,imolty,hmain+2,4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,shift+4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,shift+5)
      call twh_ijvib(GLB_SET,imolty,hmain+2,4,shift+6)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+1)
      call twh_invib(GLB_SET,imolty,shift+4,1)
      call twh_ijvib(GLB_SET,imolty,shift+4,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+5,1)
      call twh_ijvib(GLB_SET,imolty,shift+5,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+6,1)
      call twh_ijvib(GLB_SET,imolty,shift+6,1,hmain+2)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG1')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CG2')
      call twh_pdbname(GLB_SET,imolty,shift,' HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'1HG1')
      call twh_pdbname(GLB_SET,imolty,shift+2,'2HG1')
      call twh_pdbname(GLB_SET,imolty,shift+3,'3HG1')
      call twh_pdbname(GLB_SET,imolty,shift+4,'1HG2')
      call twh_pdbname(GLB_SET,imolty,shift+5,'2HG2')
      call twh_pdbname(GLB_SET,imolty,shift+6,'3HG2')

      return
      end

      subroutine twh_aa_pepname_w0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'w0' monomer in  *
c     * the polypeptide builder                                        *
c     * w0: tryptophan monomer                                         *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'tryptophan'
      shortname = 'TRP'
c     --- tryptophan side chain
      heavy = 10
      shift = hmain+heavy+ninsert
      light = 8
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,shift+2)
      call twh_invib(GLB_SET,imolty,hmain+3,3)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,hmain+6)
      call twh_invib(GLB_SET,imolty,hmain+4,3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,shift+3)
      call twh_invib(GLB_SET,imolty,hmain+5,3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+5,3,hmain+7)
      call twh_invib(GLB_SET,imolty,hmain+6,3)
      call twh_ijvib(GLB_SET,imolty,hmain+6,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+6,2,hmain+8)
      call twh_ijvib(GLB_SET,imolty,hmain+6,3,shift+4)
      call twh_invib(GLB_SET,imolty,hmain+7,3)
      call twh_ijvib(GLB_SET,imolty,hmain+7,1,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+7,2,hmain+9)
      call twh_ijvib(GLB_SET,imolty,hmain+7,3,shift+5)
      call twh_invib(GLB_SET,imolty,hmain+8,3)
      call twh_ijvib(GLB_SET,imolty,hmain+8,1,hmain+6)
      call twh_ijvib(GLB_SET,imolty,hmain+8,2,hmain+9)
      call twh_ijvib(GLB_SET,imolty,hmain+8,3,shift+6)
      call twh_invib(GLB_SET,imolty,hmain+9,3)
      call twh_ijvib(GLB_SET,imolty,hmain+9,1,hmain+7)
      call twh_ijvib(GLB_SET,imolty,hmain+9,2,hmain+8)
      call twh_ijvib(GLB_SET,imolty,hmain+9,3,shift+7)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+4)
      call twh_invib(GLB_SET,imolty,shift+4,1)
      call twh_ijvib(GLB_SET,imolty,shift+4,1,hmain+6)
      call twh_invib(GLB_SET,imolty,shift+5,1)
      call twh_ijvib(GLB_SET,imolty,shift+5,1,hmain+7)
      call twh_invib(GLB_SET,imolty,shift+6,1)
      call twh_ijvib(GLB_SET,imolty,shift+6,1,hmain+8)
      call twh_invib(GLB_SET,imolty,shift+7,1)
      call twh_ijvib(GLB_SET,imolty,shift+7,1,hmain+9)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CD1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CD2')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' NE1')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' CE2')
      call twh_pdbname(GLB_SET,imolty,hmain+6,' CE3')
      call twh_pdbname(GLB_SET,imolty,hmain+7,' CZ2')
      call twh_pdbname(GLB_SET,imolty,hmain+8,' CZ3')
      call twh_pdbname(GLB_SET,imolty,hmain+9,' CH2')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')
      call twh_pdbname(GLB_SET,imolty,shift+2,' HD1')
      call twh_pdbname(GLB_SET,imolty,shift+3,' HE1')
      call twh_pdbname(GLB_SET,imolty,shift+4,' HE3')
      call twh_pdbname(GLB_SET,imolty,shift+5,' HZ2')
      call twh_pdbname(GLB_SET,imolty,shift+6,' HZ3')
      call twh_pdbname(GLB_SET,imolty,shift+7,' HH2')

      return
      end

      subroutine twh_aa_pepname_y0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the 'y0' monomer in  *
c     * the polypeptide builder                                        *
c     * y0: tyrosine                                                   *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'tyrosine'
      shortname = 'TYR'
c     --- tyrosine side chain
      heavy = 8
      shift = hmain+heavy+ninsert
      light = 7
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,4)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,shift)
      call twh_ijvib(GLB_SET,imolty,hmain,4,shift+1)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,shift+2)
      call twh_invib(GLB_SET,imolty,hmain+3,3)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,shift+3)
      call twh_invib(GLB_SET,imolty,hmain+4,3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+6)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,shift+4)
      call twh_invib(GLB_SET,imolty,hmain+5,3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,2,hmain+6)
      call twh_ijvib(GLB_SET,imolty,hmain+5,3,shift+5)
      call twh_invib(GLB_SET,imolty,hmain+6,3)
      call twh_ijvib(GLB_SET,imolty,hmain+6,1,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+6,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+6,3,hmain+7)
      call twh_invib(GLB_SET,imolty,hmain+7,2)
      call twh_ijvib(GLB_SET,imolty,hmain+7,1,hmain+6)
      call twh_ijvib(GLB_SET,imolty,hmain+7,2,shift+6)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+4,1)
      call twh_ijvib(GLB_SET,imolty,shift+4,1,hmain+4)
      call twh_invib(GLB_SET,imolty,shift+5,1)
      call twh_ijvib(GLB_SET,imolty,shift+5,1,hmain+5)
      call twh_invib(GLB_SET,imolty,shift+6,1)
      call twh_ijvib(GLB_SET,imolty,shift+6,1,hmain+7)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CD1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CD2')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' CE1')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' CE2')
      call twh_pdbname(GLB_SET,imolty,hmain+6,' CZ')
      call twh_pdbname(GLB_SET,imolty,hmain+7,' OH')
      call twh_pdbname(GLB_SET,imolty,shift,'1HB')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HB')
      call twh_pdbname(GLB_SET,imolty,shift+2,' HD1')
      call twh_pdbname(GLB_SET,imolty,shift+3,' HD2')
      call twh_pdbname(GLB_SET,imolty,shift+4,' HE1')
      call twh_pdbname(GLB_SET,imolty,shift+5,' HE2')
      call twh_pdbname(GLB_SET,imolty,shift+6,' HH')
      
      return
      end

      subroutine twh_aa_pepname_za(lfailure,aminoname,shortname,termflag
     &     ,imolty,total,iatom)
c     ******************************************************************
c     * sets up the communal bond information for the 'za' monomer in  *
c     * the polypeptide builder                                        *
c     * za: acetylated N-terminus                                      *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer termflag,imolty,total,iatom

      if ( termflag .ne. -1 ) then
c        --- problem, can only have ace on N-term end
         write(6,*) 'AA_PEPNAME_ZA problem, termflag:',termflag
         write(6,*) 'imolty',imolty
         write(6,*) 'can only have acetylated cap on the N-term'
     &        ,' end of the polypeptide'
         lfailure = .true.
         return
      endif

      aminoname = 'acetylated cap on N-terminus'
      shortname = 'ACE'
      total = 6
      call twh_polysafetycheck(lfailure,iatom,total)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,iatom,3)
      call twh_ijvib(GLB_SET,imolty,iatom,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,iatom,2,iatom+2)
c     --- this bond is with the next group
      call twh_ijvib(GLB_SET,imolty,iatom,3,iatom+total)
      call twh_invib(GLB_SET,imolty,iatom+1,4)
      call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
      call twh_ijvib(GLB_SET,imolty,iatom+1,2,iatom+3)
      call twh_ijvib(GLB_SET,imolty,iatom+1,3,iatom+4)
      call twh_ijvib(GLB_SET,imolty,iatom+1,4,iatom+5)
      call twh_invib(GLB_SET,imolty,iatom+2,1)
      call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom)
      call twh_invib(GLB_SET,imolty,iatom+3,1)
      call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom+1)
      call twh_invib(GLB_SET,imolty,iatom+4,1)
      call twh_ijvib(GLB_SET,imolty,iatom+4,1,iatom+1)
      call twh_invib(GLB_SET,imolty,iatom+5,1)
      call twh_ijvib(GLB_SET,imolty,iatom+5,1,iatom+1)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,iatom,' CY')
      call twh_pdbname(GLB_SET,imolty,iatom+1,' CAY')
      call twh_pdbname(GLB_SET,imolty,iatom+2,' OY')
      call twh_pdbname(GLB_SET,imolty,iatom+3,'1HY')
      call twh_pdbname(GLB_SET,imolty,iatom+4,'2HY')
      call twh_pdbname(GLB_SET,imolty,iatom+5,'3HY')
  
      return
      end

      subroutine twh_aa_pepname_ze(lfailure,aminoname,shortname,termflag
     &     ,imolty,total,iatom,ibondend)
c     ******************************************************************
c     * sets up the communal bond information for the 'ze' monomer in  *
c     * the polypeptide builder                                        *
c     * ze: ethanolamine cap on C-terminus                             *
c     *                                                                *
c     * split from polycharmm22 03-27-2006 by M.G. Martin              *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer termflag,imolty,total,iatom,ibondend

c     --- ethanolamine cap on the c-terminus
      if ( termflag .ne. 1 ) then
c        --- problem, can only have ethanolamine on C-term end
         write(6,*) 'AA_PEPNAME_ZE problem: imolty:',imolty
         write(6,*) 'can only have ethanolamine cap on the C-term'
     &        ,' end of the polypeptide'
         lfailure = .true.
         return
      endif
      aminoname = 'ethanolamine cap on C-terminus'
      shortname = 'ETA'
      total = 10
      call twh_polysafetycheck(lfailure,iatom,total)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,iatom,3)
      call twh_ijvib(GLB_SET,imolty,iatom,1,ibondend)
      call twh_ijvib(GLB_SET,imolty,iatom,2,iatom+1)
      call twh_ijvib(GLB_SET,imolty,iatom,3,iatom+4)
      call twh_invib(GLB_SET,imolty,iatom+1,4)
      call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
      call twh_ijvib(GLB_SET,imolty,iatom+1,2,iatom+2)
      call twh_ijvib(GLB_SET,imolty,iatom+1,3,iatom+5)
      call twh_ijvib(GLB_SET,imolty,iatom+1,4,iatom+6)
      call twh_invib(GLB_SET,imolty,iatom+2,4)
      call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,iatom+2,2,iatom+3)
      call twh_ijvib(GLB_SET,imolty,iatom+2,3,iatom+7)
      call twh_ijvib(GLB_SET,imolty,iatom+2,4,iatom+8)
      call twh_invib(GLB_SET,imolty,iatom+3,2)
      call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom+2)
      call twh_ijvib(GLB_SET,imolty,iatom+3,2,iatom+9)
      call twh_invib(GLB_SET,imolty,iatom+4,1)
      call twh_ijvib(GLB_SET,imolty,iatom+4,1,iatom)
      call twh_invib(GLB_SET,imolty,iatom+5,1)
      call twh_ijvib(GLB_SET,imolty,iatom+5,1,iatom+1)
      call twh_invib(GLB_SET,imolty,iatom+6,1)
      call twh_ijvib(GLB_SET,imolty,iatom+6,1,iatom+1)
      call twh_invib(GLB_SET,imolty,iatom+7,1)
      call twh_ijvib(GLB_SET,imolty,iatom+7,1,iatom+2)
      call twh_invib(GLB_SET,imolty,iatom+8,1)
      call twh_ijvib(GLB_SET,imolty,iatom+8,1,iatom+2)
      call twh_invib(GLB_SET,imolty,iatom+9,1)
      call twh_ijvib(GLB_SET,imolty,iatom+9,1,iatom+3)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,iatom,' N')
      call twh_pdbname(GLB_SET,imolty,iatom+1,' C2')
      call twh_pdbname(GLB_SET,imolty,iatom+2,' C1')
      call twh_pdbname(GLB_SET,imolty,iatom+3,' O')
      call twh_pdbname(GLB_SET,imolty,iatom+4,'1HN')
      call twh_pdbname(GLB_SET,imolty,iatom+5,'1H2')
      call twh_pdbname(GLB_SET,imolty,iatom+6,'2H2')
      call twh_pdbname(GLB_SET,imolty,iatom+7,'1H1')
      call twh_pdbname(GLB_SET,imolty,iatom+8,'2H2')
      call twh_pdbname(GLB_SET,imolty,iatom+9,' HO')
  
      return
      end

      subroutine twh_aa_pepname_zf(lfailure,aminoname,shortname,termflag
     &     ,imolty,total,iatom)
c     ******************************************************************
c     * sets up the communal bond information for the 'zf' monomer in  *
c     * the polypeptide builder                                        *
c     * zf: formaldahyde cap on N-terminus                             *
c     *                                                                *
c     * split from polycharmm22 03-27-2006 by M.G. Martin              *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer termflag,imolty,total,iatom

c     --- formaldahyde cap on N-terminal end of peptide
      if ( termflag .ne. -1 ) then
c        --- problem, can only have formaldehyde on N-term end
         write(6,*) 'AA_PEPNAME_ZF problem, imolty:',imolty
         write(6,*) 'can only have formaldehyde cap on the N-term'
     &        ,' end of the polypeptide'
         lfailure = .true.
         return
      endif
      aminoname = 'formaldehyde cap on N-terminus'
      shortname = 'FOR'
      total = 3
      call twh_polysafetycheck(lfailure,iatom,total)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,iatom,3)
      call twh_ijvib(GLB_SET,imolty,iatom,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,iatom,2,iatom+2)
c     --- this bond is with the next group
      call twh_ijvib(GLB_SET,imolty,iatom,3,iatom+total)
      call twh_invib(GLB_SET,imolty,iatom+1,1)
      call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
      call twh_invib(GLB_SET,imolty,iatom+2,1)
      call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,iatom,' C')
      call twh_pdbname(GLB_SET,imolty,iatom+1,' O')
      call twh_pdbname(GLB_SET,imolty,iatom+2,' H')
  
      return
      end

      subroutine twh_aa_pepname_zm(lfailure,aminoname,shortname,termflag
     &     ,imolty,total,iatom,ibondend)
c     ******************************************************************
c     * sets up the communal bond information for the 'zm' monomer in  *
c     * the polypeptide builder                                        *
c     * zm: amidated C-terminus                                        *
c     *                                                                *
c     * split from polyamber96 03-25-2006 by M.G. Martin               *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer termflag,imolty,total,iatom,ibondend
c     --- local variables
      if ( termflag .ne. 1 ) then
c        --- problem, can only be on C-term end
         write(6,*) 'AA_PEPNAME_ZM problem, termflag:',termflag
         write(6,*) 'imolty',imolty
         write(6,*) 'can only have amidated cap on the C-term'
     &        ,' end of the polypeptide'
         lfailure = .true.
         return
      endif
      aminoname = 'Amidated term'
      shortname = 'NHE'
      total = 3
      call twh_polysafetycheck(lfailure,iatom,total)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,iatom,3)
      call twh_ijvib(GLB_SET,imolty,iatom,1,ibondend)
      call twh_ijvib(GLB_SET,imolty,iatom,2,iatom+1)
      call twh_ijvib(GLB_SET,imolty,iatom,3,iatom+2)
      call twh_invib(GLB_SET,imolty,iatom+1,1)
      call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
      call twh_invib(GLB_SET,imolty,iatom+2,1)
      call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,iatom,' NT')
      call twh_pdbname(GLB_SET,imolty,iatom+1,'1HT')
      call twh_pdbname(GLB_SET,imolty,iatom+2,'2HT')
      
      return
      end

      subroutine twh_aa_pepname_zn(lfailure,aminoname,shortname,termflag
     &     ,imolty,total,iatom,ibondend)
c     ******************************************************************
c     * sets up the communal bond information for the 'zn' monomer in  *
c     * the polypeptide builder                                        *
c     * zn: N-methylamide cap on C-terminus                            *
c     *                                                                *
c     * split from polycharmm22 03-27-2006 by M.G. Martin              *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer termflag,imolty,total,iatom,ibondend

c     --- N-methylamide cap on the c-terminus
      if ( termflag .ne. 1 ) then
c        --- problem, can only have N-methylamide on C-term end
         write(6,*) 'AA_PEPNAME_ZN problem imolty:',imolty
         write(6,*)'can only have N-methylamide cap on the C-term'
     &        ,' end of the polypeptide'
         lfailure = .true.
         return
      endif
      aminoname = 'N-methylamide'
      shortname = 'NMA'
      total = 6
      call twh_polysafetycheck(lfailure,iatom,total)
      if ( lfailure ) return
c     --- vibrations
      call twh_invib(GLB_SET,imolty,iatom,3)
      call twh_ijvib(GLB_SET,imolty,iatom,1,ibondend)
      call twh_ijvib(GLB_SET,imolty,iatom,2,iatom+1)
      call twh_ijvib(GLB_SET,imolty,iatom,3,iatom+2)
      call twh_invib(GLB_SET,imolty,iatom+1,4)
      call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
      call twh_ijvib(GLB_SET,imolty,iatom+1,2,iatom+3)
      call twh_ijvib(GLB_SET,imolty,iatom+1,3,iatom+4)
      call twh_ijvib(GLB_SET,imolty,iatom+1,4,iatom+5)
      call twh_invib(GLB_SET,imolty,iatom+2,1)
      call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom)
      call twh_invib(GLB_SET,imolty,iatom+3,1)
      call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom+1)
      call twh_invib(GLB_SET,imolty,iatom+4,1)
      call twh_ijvib(GLB_SET,imolty,iatom+4,1,iatom+1)
      call twh_invib(GLB_SET,imolty,iatom+5,1)
      call twh_ijvib(GLB_SET,imolty,iatom+5,1,iatom+1)
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,iatom,' NT')
      call twh_pdbname(GLB_SET,imolty,iatom+1,' CAT')
      call twh_pdbname(GLB_SET,imolty,iatom+2,' HNT')
      call twh_pdbname(GLB_SET,imolty,iatom+3,'1HT')
      call twh_pdbname(GLB_SET,imolty,iatom+4,'2HT')
      call twh_pdbname(GLB_SET,imolty,iatom+5,'3HT')

      return
      end

      subroutine twh_ua_pepname_generic(lfailure,termflag,hmain,iatom
     &     ,ninsert)
c     ******************************************************************
c     * sets up the generic information for most united-atom peptides  *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer termflag,hmain,iatom,ninsert

      if ( termflag .eq. -1 ) then
c        --- n-terminus
         hmain = iatom+4
c        --- lump H on to C
         ninsert = 3
      elseif ( termflag .eq. 0 ) then
c        --- main chain
         hmain = iatom+4
c        --- lump H on to C
         ninsert = 1
      elseif ( termflag .eq. 1 ) then
c        --- c-terminus
         hmain = iatom+4
c        --- lump H on to C
         ninsert = 2
      else
         write(6,*) 'invalid termflag'
         lfailure = .true.
         return
      endif

      return
      end

      subroutine twh_ua_pepname_main_chain(lfailure,termflag
     &     ,imolty,total,iatom,ibondend,heavy,light,ninsert
     &     ,hmain,lsettorsion,tortype)
c     ******************************************************************
c     * sets up the main chain generic information for most all-atom   *
c     * peptides                                                       *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      integer termflag,imolty,total,iatom,ibondend,heavy
      integer hmain,ninsert,light,tortype

c     --- generic amino acid setup
      total = heavy+ninsert+hmain+light
      call twh_monomer_check(lfailure,total)
      if ( lfailure ) return
c     --- subtract iatom to get the number of new atoms for this monomer
      total = total - iatom

      if ( termflag .eq. -1 ) then
c        --- n-terminus main chain
c        --- pdb names
         call twh_pdbname(GLB_SET,imolty,iatom,' N')
         call twh_pdbname(GLB_SET,imolty,iatom+1,' CA')
         call twh_pdbname(GLB_SET,imolty,iatom+2,' C')
         call twh_pdbname(GLB_SET,imolty,iatom+3,' O')
         call twh_pdbname(GLB_SET,imolty,hmain+heavy,'1H')
         call twh_pdbname(GLB_SET,imolty,hmain+heavy+1,'2H')
         call twh_pdbname(GLB_SET,imolty,hmain+heavy+2,'3H')
c        --- vibrations
         call twh_invib(GLB_SET,imolty,iatom,4)
         call twh_ijvib(GLB_SET,imolty,iatom,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom,2,hmain+heavy)
         call twh_ijvib(GLB_SET,imolty,iatom,3,hmain+heavy+1)
         call twh_ijvib(GLB_SET,imolty,iatom,4,hmain+heavy+2)
         if ( heavy .ne. 0 ) then
c           --- there is a side chain
            call twh_invib(GLB_SET,imolty,iatom+1,3)
            call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
            call twh_ijvib(GLB_SET,imolty,iatom+1,2,iatom+2)
            call twh_ijvib(GLB_SET,imolty,iatom+1,3,hmain)
         else
c           --- no side chain (glycine)
            call twh_invib(GLB_SET,imolty,iatom+1,2)
            call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
            call twh_ijvib(GLB_SET,imolty,iatom+1,2,iatom+2)
         endif
         call twh_invib(GLB_SET,imolty,iatom+2,3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom+2,2,iatom+3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,3,iatom+total)
         call twh_invib(GLB_SET,imolty,iatom+3,1)
         call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom+2)
         call twh_invib(GLB_SET,imolty,hmain+heavy,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy,1,iatom)
         call twh_invib(GLB_SET,imolty,hmain+heavy+1,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy+1,1,iatom)
         call twh_invib(GLB_SET,imolty,hmain+heavy+2,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy+2,1,iatom)
         if ( lsettorsion ) then
c           --- torsions
c           --- backbone N - Ca - C - next N
            call twh_create_torsion(imolty,iatom,iatom+1,iatom+2
     &           ,iatom+total,tortype)
c           --- backbone Ca - C - next N - next Ca
            call twh_create_torsion(imolty,iatom+1,iatom+2,iatom+total
     &           ,iatom+total+1,tortype)
c           --- reverse backbone C - next N - next Ca - next C
            call twh_create_torsion(imolty,iatom+2,iatom+total
     &           ,iatom+total+1,iatom+total+2,tortype)
         endif
      elseif ( termflag .eq. 0 ) then
c        --- interior segment main chain
c        --- pdb names
         call twh_pdbname(GLB_SET,imolty,iatom,' N')
         call twh_pdbname(GLB_SET,imolty,iatom+1,' CA')
         call twh_pdbname(GLB_SET,imolty,iatom+2,' C')
         call twh_pdbname(GLB_SET,imolty,iatom+3,' O')
         call twh_pdbname(GLB_SET,imolty,hmain+heavy,' H')
c        --- vibrations
         call twh_invib(GLB_SET,imolty,iatom,3)
         call twh_ijvib(GLB_SET,imolty,iatom,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom,2,hmain+heavy)
         call twh_ijvib(GLB_SET,imolty,iatom,3,ibondend)
         if ( heavy .ne. 0 ) then
c           --- there is a side chain
            call twh_invib(GLB_SET,imolty,iatom+1,3)
            call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
            call twh_ijvib(GLB_SET,imolty,iatom+1,2,iatom+2)
            call twh_ijvib(GLB_SET,imolty,iatom+1,3,hmain)
         else
c           --- glycine only has H which is lumped in this ff
            call twh_invib(GLB_SET,imolty,iatom+1,2)
            call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
            call twh_ijvib(GLB_SET,imolty,iatom+1,2,iatom+2)
         endif
         call twh_invib(GLB_SET,imolty,iatom+2,3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom+2,2,iatom+3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,3,iatom+total)
         call twh_invib(GLB_SET,imolty,iatom+3,1)
         call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom+2)
         call twh_invib(GLB_SET,imolty,hmain+heavy,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy,1,iatom)
         if ( lsettorsion ) then
c           --- torsions
c           --- backbone N - Ca - C - next N
            call twh_create_torsion(imolty,iatom,iatom+1,iatom+2
     &           ,iatom+total,tortype)
c           --- reverse backbone N - prev C - prev Ca - prev N
            call twh_create_torsion(imolty,iatom,ibondend,ibondend-1
     &           ,ibondend-2,tortype)
c           --- backbone Ca - C - next N - next Ca
            call twh_create_torsion(imolty,iatom+1,iatom+2,iatom+total
     &           ,iatom+total+1,tortype)
c           --- reverse backbone Ca - N - prev C - prev Ca
            call twh_create_torsion(imolty,iatom+1,iatom,ibondend
     &           ,ibondend-1,tortype)
c           ---- backbone C - Ca - N - prev C
            call twh_create_torsion(imolty,iatom+2,iatom+1,iatom
     &           ,ibondend,tortype)
c           --- reverse backbone C - next N - next Ca - next C
            call twh_create_torsion(imolty,iatom+2,iatom+total
     &           ,iatom+total+1,iatom+total+2,tortype)
         endif
      elseif ( termflag .eq. 1 ) then
c        --- c-terminus main chain
c        --- pdb names
         call twh_pdbname(GLB_SET,imolty,iatom,' N')
         call twh_pdbname(GLB_SET,imolty,iatom+1,' CA')
         call twh_pdbname(GLB_SET,imolty,iatom+2,' C')
         call twh_pdbname(GLB_SET,imolty,iatom+3,' O')
         call twh_pdbname(GLB_SET,imolty,hmain+heavy,' OXT')
         call twh_pdbname(GLB_SET,imolty,hmain+heavy+1,' H')
c        --- vibrations
         call twh_invib(GLB_SET,imolty,iatom,3)
         call twh_ijvib(GLB_SET,imolty,iatom,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom,2,hmain+heavy+1)
         call twh_ijvib(GLB_SET,imolty,iatom,3,ibondend)
         if ( heavy .ne. 0 ) then
            call twh_invib(GLB_SET,imolty,iatom+1,3)
            call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
            call twh_ijvib(GLB_SET,imolty,iatom+1,2,iatom+2)
            call twh_ijvib(GLB_SET,imolty,iatom+1,3,hmain)
         else
            call twh_invib(GLB_SET,imolty,iatom+1,2)
            call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
            call twh_ijvib(GLB_SET,imolty,iatom+1,2,iatom+2)
         endif
         call twh_invib(GLB_SET,imolty,iatom+2,3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom+2,2,iatom+3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,3,hmain+heavy)
         call twh_invib(GLB_SET,imolty,iatom+3,1)
         call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom+2)
         call twh_invib(GLB_SET,imolty,hmain+heavy,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy,1,iatom+2)
         call twh_invib(GLB_SET,imolty,hmain+heavy+1,1)
         call twh_ijvib(GLB_SET,imolty,hmain+heavy+1,1,iatom)
         if ( lsettorsion ) then
c           --- torsions
c           --- reverse backbone N - prev C - prev Ca - prev N
            call twh_create_torsion(imolty,iatom,ibondend,ibondend-1
     &           ,ibondend-2,tortype)
c           --- reverse backbone Ca - N - prev C - prev Ca
            call twh_create_torsion(imolty,iatom+1,iatom,ibondend
     &           ,ibondend-1,tortype)
c           ---- backbone C - Ca - N - prev C
            call twh_create_torsion(imolty,iatom+2,iatom+1,iatom
     &           ,ibondend,tortype)
         endif
      else
c        --- invalid termflag
         write(6,*) 'UA_MAIN_CHAIN unknown termflag'
         write(6,*) termflag
         lfailure = .true.
         return
      endif

      return
      end

      subroutine twh_ua_pepname_a0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'a0' *
c     * monomer in the polypeptide builder                             *
c     * a0: alanine                                                    *
c     * this monomer makes use of the generic pepname features         *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert

      aminoname = 'alanine'
      shortname = 'ALA'
c     --- alanine side chain
c     --- united-atom geometry
      heavy = 1
      shift = hmain+heavy+ninsert
      light = 0
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- vibrations
c     --- heavy atoms
      call twh_invib(GLB_SET,imolty,hmain,1)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
c     --- assign pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')

      return
      end

      subroutine twh_ua_pepname_c0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'c0' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype

      aminoname = 'cysteine with hydrogen'
      shortname = 'CYS'
c     --- cysteine side chain
      heavy = 2
      shift = hmain+heavy+ninsert
      light = 0
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb atom names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' SG')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,2)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+1,1)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      if ( lsettorsion ) then
c        --- side chain torsions
c        --- N - Ca - Cb - Sg
         call twh_create_torsion(imolty,iatom,iatom+1,hmain
     &        ,hmain+1,tortype)
c        --- reverse: Sg - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
      endif

      return
      end

      subroutine twh_ua_pepname_d_neg(lfailure,aminoname,shortname
     &     ,imolty,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'd-' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype

      aminoname = 'aspartic acid deprotonated'
      shortname = 'ASP'
c     --- asparagine side chain
      heavy = 4
      shift = hmain+heavy+ninsert
      light = 0
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' OD1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' OD2')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,2)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+3,1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      if ( lsettorsion ) then
c        --- side chain torsions
c        --- N - Ca - Cb - Cg
         call twh_create_torsion(imolty,iatom,iatom+1,hmain,hmain+1
     &        ,tortype)
c        --- reverse: Cg - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
c        --- Ca - Cb - Cg - Od1
         call twh_create_torsion(imolty,iatom+1,hmain,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Od1 - Cg - Cb - Ca
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain
     &        ,iatom+1,tortype)
      endif

      return
      end

      subroutine twh_ua_pepname_e_neg(lfailure,aminoname,shortname
     &     ,imolty,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'e-' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype

      aminoname = 'glutamic acid deprotonated'
      shortname = 'GLU'
c     --- glutamic acid side chain
      heavy = 5
      shift = hmain+heavy+ninsert
      light = 0
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CD')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' OE1')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' OE2')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,2)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+1,2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_invib(GLB_SET,imolty,hmain+2,3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,hmain+4)
      call twh_invib(GLB_SET,imolty,hmain+3,1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+2)
      call twh_invib(GLB_SET,imolty,hmain+4,1)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+2)
      if ( lsettorsion ) then
c        --- side chain torsions
c        --- N - Ca - Cb - Cg
         call twh_create_torsion(imolty,iatom,iatom+1,hmain
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
c        --- Ca - Cb - Cg - Cd
         call twh_create_torsion(imolty,iatom+1,hmain,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Cd - Cg - Cb - Ca
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain
     &        ,iatom+1,tortype)
c        --- Cb - Cg - Cd - Oe1
         call twh_create_torsion(imolty,hmain,hmain+1,hmain+2
     &        ,hmain+3,tortype)
c        --- reverse: Oe1 - Cd - Cg - Cb
         call twh_create_torsion(imolty,hmain+3,hmain+2,hmain+1
     &        ,hmain,tortype)
      endif

      return
      end

      subroutine twh_ua_pepname_f0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'f0' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype

      aminoname = 'phenylalanine'
      shortname = 'PHE'
c     --- phenylalanine side chain
      heavy = 7
      shift = hmain+heavy+ninsert
      light = 0
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CD1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CD2')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' CE1')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' CE2')
      call twh_pdbname(GLB_SET,imolty,hmain+6,' CZ')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,2)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,2)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+4)
      call twh_invib(GLB_SET,imolty,hmain+3,2)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+5)
      call twh_invib(GLB_SET,imolty,hmain+4,2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+6)
      call twh_invib(GLB_SET,imolty,hmain+5,2)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,2,hmain+6)
      call twh_invib(GLB_SET,imolty,hmain+6,2)
      call twh_ijvib(GLB_SET,imolty,hmain+6,1,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+6,2,hmain+5)
      if ( lsettorsion ) then
c        --- side chain torsions
c        --- N - Ca - Cb - Cg
         call twh_create_torsion(imolty,iatom,iatom+1,hmain
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
c        --- Ca - Cb - Cg - Cd1
         call twh_create_torsion(imolty,iatom+1,hmain,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Cd1 - Cg - Cb - Ca
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain
     &        ,iatom+1,tortype)
c        --- regular torsions listed as improper torsions in charmm
c        --- Cg - Cd1 - Ce1 - Cz
         call twh_create_torsion(imolty,hmain+1,hmain+2,hmain+4
     &        ,hmain+6,tortype)
c        --- reverse: Cz - Ce1 - Cd1 - Cg
         call twh_create_torsion(imolty,hmain+6,hmain+4,hmain+2
     &        ,hmain+1,tortype)
c        --- Cd1 - Ce1 - Cz - Ce2
         call twh_create_torsion(imolty,hmain+2,hmain+4,hmain+6
     &        ,hmain+5,tortype)
c        --- reverse: Ce2 - Cz - Ce1 - Cd1
         call twh_create_torsion(imolty,hmain+5,hmain+6,hmain+4
     &        ,hmain+2,tortype)
c        --- Ce1 - Cz - Ce2 - Cd2
         call twh_create_torsion(imolty,hmain+4,hmain+6,hmain+5
     &        ,hmain+3,tortype)
c        --- reverse: Cd2 - Ce2 - Cz - Ce1
         call twh_create_torsion(imolty,hmain+3,hmain+5,hmain+6
     &        ,hmain+4,tortype)
c        --- Cz - Ce2 - Cd2 - Cg
         call twh_create_torsion(imolty,hmain+6,hmain+5,hmain+3
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cd2 - Ce2 - Cz
         call twh_create_torsion(imolty,hmain+1,hmain+3,hmain+5
     &        ,hmain+6,tortype)
c        --- Ce2 - Cd2 - Cg - Cd1
         call twh_create_torsion(imolty,hmain+5,hmain+3,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Cd1 - Cg - Cd2 - Ce2
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain+3
     &        ,hmain+5,tortype)
c        --- Cd2 - Cg - Cd1 - Ce1
         call twh_create_torsion(imolty,hmain+3,hmain+1,hmain+2
     &        ,hmain+4,tortype)
c        --- reverse: Ce1 - Cd1 - Cg - Cd2
         call twh_create_torsion(imolty,hmain+4,hmain+2,hmain+1
     &        ,hmain+3,tortype)
      endif

      return
      end

      subroutine twh_ua_pepname_g0(lfailure,aminoname,shortname
     &     ,hmain,heavy,light,shift,ninsert)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'g0' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*3 shortname
      character*40 aminoname
      integer hmain,heavy,light,shift,ninsert

      aminoname = 'glycine'
      shortname = 'GLY'
c     --- side chain glycine has no atoms in united-atom approach
      heavy = 0
      shift = hmain+heavy+ninsert
      light = 0
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return

      return
      end

      subroutine twh_ua_pepname_hd(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'hd' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype
      
      aminoname = 'histidine neutral only Nd protonated'
      shortname = 'HIS'
c     --- histidine side chain
      heavy = 6
      shift = hmain+heavy+ninsert
      light = 1
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' ND1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CD2')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' CE1')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' NE2')
      call twh_pdbname(GLB_SET,imolty,shift,' HD1')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,2)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,shift)
      call twh_invib(GLB_SET,imolty,hmain+3,2)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+5)
      call twh_invib(GLB_SET,imolty,hmain+4,2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+5)
      call twh_invib(GLB_SET,imolty,hmain+5,2)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,2,hmain+4)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain+2)
      if ( lsettorsion ) then
c        --- torsions that are listed in charmm19 as torsions
c        --- N - Ca - Cb - Cg
         call twh_create_torsion(imolty,iatom,iatom+1,hmain
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
c        --- Ca - Cb - Cg - Nd1
         call twh_create_torsion(imolty,iatom+1,hmain,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Nd1 - Cg - Cb - Ca
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain
     &        ,iatom+1,tortype)
c        --- torsions that are listed in charmm19 as impropers
c        --- Cg - Nd1 - Ce1 - Ne2
         call twh_create_torsion(imolty,hmain+1,hmain+2,hmain+4
     &        ,hmain+5,tortype)
c        --- reverse: Ne2 - Ce1 - Nd1 - Cg
         call twh_create_torsion(imolty,hmain+5,hmain+4,hmain+2
     &        ,hmain+1,tortype)
c        --- Nd1 - Ce1 - Ne2 - Cd2
         call twh_create_torsion(imolty,hmain+2,hmain+4,hmain+5
     &        ,hmain+3,tortype)
c        --- reverse: Cd2 - Ne2 - Ce1 - Nd1
         call twh_create_torsion(imolty,hmain+3,hmain+5,hmain+4
     &        ,hmain+2,tortype)
c        --- Ce1 - Ne2 - Cd2 - Cg
         call twh_create_torsion(imolty,hmain+4,hmain+5,hmain+3
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cd2 - Ne2 - Ce1
         call twh_create_torsion(imolty,hmain+1,hmain+3,hmain+5
     &        ,hmain+4,tortype)
c        --- Ne2 - Cd2 - Cg - Nd1
         call twh_create_torsion(imolty,hmain+5,hmain+3,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Nd1 - Cg - Cd2 - Ne2
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain+3
     &        ,hmain+5,tortype)
c        --- Cd2 - Cg - Nd1 - Ce1
         call twh_create_torsion(imolty,hmain+3,hmain+1,hmain+2
     &        ,hmain+4,tortype)
c        --- reverse: Ce1 - Nd1 - Cg - Cd2
         call twh_create_torsion(imolty,hmain+4,hmain+2,hmain+1
     &        ,hmain+3,tortype)
      endif

      return
      end

      subroutine twh_ua_pepname_he(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'he' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype
      
      aminoname = 'histidine neutral only Ne protonated'
      shortname = 'HIS'
c     --- histidine side chain
      heavy = 6
      shift = hmain+heavy+ninsert
      light = 1
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' ND1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CD2')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' CE1')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' NE2')
      call twh_pdbname(GLB_SET,imolty,shift,' HE2')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,2)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,2)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+4)
      call twh_invib(GLB_SET,imolty,hmain+3,2)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+5)
      call twh_invib(GLB_SET,imolty,hmain+4,2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+5)
      call twh_invib(GLB_SET,imolty,hmain+5,3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+5,3,shift)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain+5)
      if ( lsettorsion ) then
c        --- torsions that are listed in charmm19 as torsions
c        --- N - Ca - Cb - Cg
         call twh_create_torsion(imolty,iatom,iatom+1,hmain
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
c        --- Ca - Cb - Cg - Nd1
         call twh_create_torsion(imolty,iatom+1,hmain,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Nd1 - Cg - Cb - Ca
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain
     &        ,iatom+1,tortype)
c        --- torsions that are listed in charmm19 as impropers
c        --- Cg - Nd1 - Ce1 - Ne2
         call twh_create_torsion(imolty,hmain+1,hmain+2,hmain+4
     &        ,hmain+5,tortype)
c        --- reverse: Ne2 - Ce1 - Nd1 - Cg
         call twh_create_torsion(imolty,hmain+5,hmain+4,hmain+2
     &        ,hmain+1,tortype)
c        --- Nd1 - Ce1 - Ne2 - Cd2
         call twh_create_torsion(imolty,hmain+2,hmain+4,hmain+5
     &        ,hmain+3,tortype)
c        --- reverse: Cd2 - Ne2 - Ce1 - Nd1
         call twh_create_torsion(imolty,hmain+3,hmain+5,hmain+4
     &        ,hmain+2,tortype)
c        --- Ce1 - Ne2 - Cd2 - Cg
         call twh_create_torsion(imolty,hmain+4,hmain+5,hmain+3
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cd2 - Ne2 - Ce1
         call twh_create_torsion(imolty,hmain+1,hmain+3,hmain+5
     &        ,hmain+4,tortype)
c        --- Ne2 - Cd2 - Cg - Nd1
         call twh_create_torsion(imolty,hmain+5,hmain+3,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Nd1 - Cg - Cd2 - Ne2
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain+3
     &        ,hmain+5,tortype)
c        --- Cd2 - Cg - Nd1 - Ce1
         call twh_create_torsion(imolty,hmain+3,hmain+1,hmain+2
     &        ,hmain+4,tortype)
c        --- reverse: Ce1 - Nd1 - Cg - Cd2
         call twh_create_torsion(imolty,hmain+4,hmain+2,hmain+1
     &        ,hmain+3,tortype)
      endif

      return
      end

      subroutine twh_ua_pepname_h_plus(lfailure,aminoname,shortname
     &     ,imolty,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'h+' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype
      
      aminoname = 'histidine both N protonated'
      shortname = 'HIS'
c     --- histidine side chain with both N protonated
      heavy = 6
      shift = hmain+heavy+ninsert
      light = 2
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' ND1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CD2')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' CE1')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' NE2')
      call twh_pdbname(GLB_SET,imolty,shift,' HD1')
      call twh_pdbname(GLB_SET,imolty,shift+1,' HE2')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,2)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,shift)
      call twh_invib(GLB_SET,imolty,hmain+3,2)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+5)
      call twh_invib(GLB_SET,imolty,hmain+4,2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+5)
      call twh_invib(GLB_SET,imolty,hmain+5,3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+5,3,shift+1)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain+2)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain+5)
      if ( lsettorsion ) then
c        --- torsions that are listed in charmm19 as torsions
c        --- N - Ca - Cb - Cg
         call twh_create_torsion(imolty,iatom,iatom+1,hmain
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
c        --- Ca - Cb - Cg - Nd1
         call twh_create_torsion(imolty,iatom+1,hmain,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Nd1 - Cg - Cb - Ca
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain
     &        ,iatom+1,tortype)
c        --- torsions that are listed in charmm19 as impropers
c        --- Cg - Nd1 - Ce1 - Ne2
         call twh_create_torsion(imolty,hmain+1,hmain+2,hmain+4
     &        ,hmain+5,tortype)
c        --- reverse: Ne2 - Ce1 - Nd1 - Cg
         call twh_create_torsion(imolty,hmain+5,hmain+4,hmain+2
     &        ,hmain+1,tortype)
c        --- Nd1 - Ce1 - Ne2 - Cd2
         call twh_create_torsion(imolty,hmain+2,hmain+4,hmain+5
     &        ,hmain+3,tortype)
c        --- reverse: Cd2 - Ne2 - Ce1 - Nd1
         call twh_create_torsion(imolty,hmain+3,hmain+5,hmain+4
     &        ,hmain+2,tortype)
c        --- Ce1 - Ne2 - Cd2 - Cg
         call twh_create_torsion(imolty,hmain+4,hmain+5,hmain+3
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cd2 - Ne2 - Ce1
         call twh_create_torsion(imolty,hmain+1,hmain+3,hmain+5
     &        ,hmain+4,tortype)
c        --- Ne2 - Cd2 - Cg - Nd1
         call twh_create_torsion(imolty,hmain+5,hmain+3,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Nd1 - Cg - Cd2 - Ne2
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain+3
     &        ,hmain+5,tortype)
c        --- Cd2 - Cg - Nd1 - Ce1
         call twh_create_torsion(imolty,hmain+3,hmain+1,hmain+2
     &        ,hmain+4,tortype)
c        --- reverse: Ce1 - Nd1 - Cg - Cd2
         call twh_create_torsion(imolty,hmain+4,hmain+2,hmain+1
     &        ,hmain+3,tortype)
      endif

      return
      end

      subroutine twh_ua_pepname_i0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'i0' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype

      aminoname = 'isoleucine'
      shortname = 'ILE'
c     --- side chain isoleucine
      heavy = 4
      shift = hmain+heavy+ninsert
      light = 0
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG1')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CG2')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CD1')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,3)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,hmain+2)
      call twh_invib(GLB_SET,imolty,hmain+1,2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain)
      call twh_invib(GLB_SET,imolty,hmain+3,1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      if ( lsettorsion ) then
c        --- torsions
c        --- N - Ca - Cb - Cg1
         call twh_create_torsion(imolty,iatom,iatom+1,hmain
     &        ,hmain+1,tortype)
c        --- reverse: Cg1 - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
c        --- Ca - Cb - Cg1 - Cd
         call twh_create_torsion(imolty,iatom+1,hmain,hmain+1
     &        ,hmain+3,tortype)
c        --- reverse: Cd - Cg1 - Cb - Ca
         call twh_create_torsion(imolty,hmain+3,hmain+1,hmain
     &        ,iatom+1,tortype)
      endif

      return
      end

      subroutine twh_ua_pepname_k_plus(lfailure,aminoname,shortname
     &     ,imolty,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'k+' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype

      aminoname = 'lysine protonated'
      shortname = 'LYS'
c     --- lysine side chain
      heavy = 5
      shift = hmain+heavy+ninsert
      light = 3
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CD')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CE')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' NZ')
      call twh_pdbname(GLB_SET,imolty,shift,'1HZ')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HZ')
      call twh_pdbname(GLB_SET,imolty,shift+2,'3HZ')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,2)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+1,2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_invib(GLB_SET,imolty,hmain+2,2)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+3,2)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+4)
      call twh_invib(GLB_SET,imolty,hmain+4,4)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,shift)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,shift+1)
      call twh_ijvib(GLB_SET,imolty,hmain+4,4,shift+2)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain+4)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain+4)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+4)
      if ( lsettorsion ) then
c        --- torsions
c        --- N - Ca - Cb - Cg
         call twh_create_torsion(imolty,iatom,iatom+1,hmain
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
c        --- Ca - Cb - Cg - Cd
         call twh_create_torsion(imolty,iatom+1,hmain,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Cd - Cg - Cb - Ca
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain
     &        ,iatom+1,tortype)
c        --- Cb - Cg - Cd - Ce
         call twh_create_torsion(imolty,hmain,hmain+1,hmain+2
     &        ,hmain+3,tortype)
c        --- reverse: Ce - Cd - Cg - Cb
         call twh_create_torsion(imolty,hmain+3,hmain+2,hmain+1
     &        ,hmain,tortype)
c        --- Cg - Cd - Ce - Nz
         call twh_create_torsion(imolty,hmain+1,hmain+2,hmain+3
     &        ,hmain+4,tortype)
c        --- reverse: Nz - Ce - Cd - Cg
         call twh_create_torsion(imolty,hmain+4,hmain+3,hmain+2
     &        ,hmain+1,tortype)
c        --- Cd - Ce - Nz - Hz1
         call twh_create_torsion(imolty,hmain+2,hmain+3,hmain+4
     &        ,shift,tortype)
c        --- reverse: Hz1 - Nz - Ce - Cd
         call twh_create_torsion(imolty,shift,hmain+4,hmain+3
     &        ,hmain+2,tortype)
c        --- Cd - Ce - Nz - Hz2
         call twh_create_torsion(imolty,hmain+2,hmain+3,hmain+4
     &        ,shift+1,tortype)
c        --- reverse: Hz2 - Nz - Ce - Cd
         call twh_create_torsion(imolty,shift+1,hmain+4,hmain+3
     &        ,hmain+2,tortype)
c        --- Cd - Ce - Nz - Hz3
         call twh_create_torsion(imolty,hmain+2,hmain+3,hmain+4
     &        ,shift+2,tortype)
c        --- reverse: Hz3 - Nz - Ce - Cd
         call twh_create_torsion(imolty,shift+2,hmain+4,hmain+3
     &        ,hmain+2,tortype)
      endif

      return
      end

      subroutine twh_ua_pepname_l0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'l0' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype

      aminoname = 'leucine'
      shortname = 'LEU'
c     --- side chain leucine
      heavy = 4
      shift = hmain+heavy+ninsert
      light = 0
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CD1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CD2')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,2)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+3,1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      if ( lsettorsion ) then
c        --- torsions
c        --- N - Ca - Cb - Cg
         call twh_create_torsion(imolty,iatom,iatom+1,hmain
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
      endif

      return
      end

      subroutine twh_ua_pepname_m0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'm0' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype

      aminoname = 'methionine'
      shortname = 'MET'
c     --- methionine side chain
      heavy = 4
      shift = hmain+heavy+ninsert
      light = 0
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' SD')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CE')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,2)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+1,2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_invib(GLB_SET,imolty,hmain+2,2)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+3,1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+2)
      if ( lsettorsion ) then
c        --- side chain torsions
c        --- N - Ca - Cb - Cg
         call twh_create_torsion(imolty,iatom,iatom+1,hmain
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
c        --- Ca - Cb - Cg - Sd
         call twh_create_torsion(imolty,iatom+1,hmain,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Sd - Cg - Cb - Ca
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain
     &        ,iatom+1,tortype)
c        --- Cb - Cg - Sd - Ce
         call twh_create_torsion(imolty,hmain,hmain+1,hmain+2
     &        ,hmain+3,tortype)
c        --- reverse: Ce - Sd - Cg - Cb
         call twh_create_torsion(imolty,hmain+3,hmain+2,hmain+1
     &        ,hmain,tortype)
      endif

      return
      end

      subroutine twh_ua_pepname_n0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'n0' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype

      aminoname = 'asparagine'
      shortname = 'ASN'
c     --- asparagine side chain
      heavy = 4
      shift = hmain+heavy+ninsert
      light = 2
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- set pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' OD1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' ND2')
      call twh_pdbname(GLB_SET,imolty,shift,'1HD2')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HD2')
c     --- vibrations
c     --- heavy atoms
      call twh_invib(GLB_SET,imolty,hmain,2)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+3,3)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,shift)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,shift+1)
c     --- hydrogens
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain+3)
      if ( lsettorsion ) then
c        --- side chain torsions
c        --- N - Ca - Cb - Cg
         call twh_create_torsion(imolty,iatom,iatom+1,hmain
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
c        --- Ca - Cb - Cg - Od1
         call twh_create_torsion(imolty,iatom+1,hmain,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Od1 - Cg - Cb - Ca
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain
     &        ,iatom+1,tortype)
c        --- Cb - Cg - Nd2 - Hd21
         call twh_create_torsion(imolty,hmain,hmain+1,hmain+3
     &        ,shift,tortype)
c        --- reverse: Hd21 - Nd2 - Cg - Cb
         call twh_create_torsion(imolty,shift,hmain+3,hmain+1
     &        ,hmain,tortype)
      endif

      return
      end


      subroutine twh_ua_pepname_p0(lfailure,aminoname,shortname,termflag
     &     ,imolty,total,iatom,icb,iha,ibondend,tortype,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'p0' *
c     * monomer in the polypeptide builder                             *
c     * p0: proline                                                    *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer termflag,imolty,total,iatom,icb,iha,ibondend,tortype

c     --- proline requires its own set of special cases
      aminoname = 'proline'
      shortname = 'PRO'
      if ( termflag .eq. -1 ) then
c        --- N-terminal proline
         write(6,*) 'UA_PEPNAME_P0: lacking a set '
         write(6,*) 'of parameters for N-terminal proline'
         lfailure = .true.
         return
      elseif ( termflag .eq. 0 ) then
c        --- interior proline
c        --- set icb
         icb = iatom+4
c        --- set iha
         iha = iatom
         total = 7
c        --- pdb names
         call twh_pdbname(GLB_SET,imolty,iatom,' N')
         call twh_pdbname(GLB_SET,imolty,iatom+1,' CA')
         call twh_pdbname(GLB_SET,imolty,iatom+2,' C')
         call twh_pdbname(GLB_SET,imolty,iatom+3,' O')
         call twh_pdbname(GLB_SET,imolty,iatom+4,' CB')
         call twh_pdbname(GLB_SET,imolty,iatom+5,' CG')
         call twh_pdbname(GLB_SET,imolty,iatom+6,' CD')
c        --- vibrations
         call twh_invib(GLB_SET,imolty,iatom,3)
         call twh_ijvib(GLB_SET,imolty,iatom,1,ibondend)
         call twh_ijvib(GLB_SET,imolty,iatom,2,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom,3,iatom+6)
         call twh_invib(GLB_SET,imolty,iatom+1,3)
         call twh_ijvib(GLB_SET,imolty,iatom+1,1,iatom)
         call twh_ijvib(GLB_SET,imolty,iatom+1,2,iatom+2)
         call twh_ijvib(GLB_SET,imolty,iatom+1,3,iatom+4)
         call twh_invib(GLB_SET,imolty,iatom+2,3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom+2,2,iatom+3)
         call twh_ijvib(GLB_SET,imolty,iatom+2,3,iatom+total)
         call twh_invib(GLB_SET,imolty,iatom+3,1)
         call twh_ijvib(GLB_SET,imolty,iatom+3,1,iatom+2)
         call twh_invib(GLB_SET,imolty,iatom+4,2)
         call twh_ijvib(GLB_SET,imolty,iatom+4,1,iatom+1)
         call twh_ijvib(GLB_SET,imolty,iatom+4,2,iatom+5)
         call twh_invib(GLB_SET,imolty,iatom+5,2)
         call twh_ijvib(GLB_SET,imolty,iatom+5,1,iatom+4)
         call twh_ijvib(GLB_SET,imolty,iatom+5,2,iatom+6)
         call twh_invib(GLB_SET,imolty,iatom+6,2)
         call twh_ijvib(GLB_SET,imolty,iatom+6,1,iatom)
         call twh_ijvib(GLB_SET,imolty,iatom+6,2,iatom+5)
         if ( lsettorsion ) then
c           --- all of the torsions for proline
c           --- N - Ca - Cb - Cg
            call twh_create_torsion(imolty,iatom,iatom+1,iatom+4,iatom+5
     &        ,tortype)
c           --- reverse: Cg - Cb - Ca - N
            call twh_create_torsion(imolty,iatom+5,iatom+4,iatom+1,iatom
     &           ,tortype)
c           --- Ca - Cb - Cg - Cd
            call twh_create_torsion(imolty,iatom+1,iatom+4,iatom+5
     &           ,iatom+6,tortype)
c           --- reverse: Cd - Cg - Cb - Ca
            call twh_create_torsion(imolty,iatom+6,iatom+5,iatom+4
     &           ,iatom+1,tortype)
c           --- Cb - Cg - Cd - N
            call twh_create_torsion(imolty,iatom+4,iatom+5,iatom+6
     &           ,iatom,tortype)
c           --- reverse: N - Cd - Cg - Cb
            call twh_create_torsion(imolty,iatom,iatom+6,iatom+5,iatom+4
     &           ,tortype)
c           --- Cg - Cd - N - Ca
            call twh_create_torsion(imolty,iatom+5,iatom+6,iatom,iatom+1
     &           ,tortype)
c           --- reverse: Ca - N - Cd - Cg
            call twh_create_torsion(imolty,iatom+1,iatom,iatom+6,iatom+5
     &           ,tortype)
c           --- main chain torsions
c           --- backbone N - Ca - C - next N
            call twh_create_torsion(imolty,iatom,iatom+1,iatom+2
     &           ,iatom+total,tortype)
c           --- reverse backbone N - prev C - prev Ca - prev N
            call twh_create_torsion(imolty,iatom,ibondend,ibondend-1
     &           ,ibondend-2,tortype)
c           --- backbone Ca - C - next N - next Ca
            call twh_create_torsion(imolty,iatom+1,iatom+2,iatom+total
     &           ,iatom+total+1,tortype)
c           --- reverse backbone Ca - N - prev C - prev Ca
            call twh_create_torsion(imolty,iatom+1,iatom,ibondend
     &           ,ibondend-1,tortype)
c           ---- backbone C - Ca - N - prev C
            call twh_create_torsion(imolty,iatom+2,iatom+1,iatom
     &           ,ibondend,tortype)
c           --- reverse backbone C - next N - next Ca - next C
            call twh_create_torsion(imolty,iatom+2,iatom+total
     &           ,iatom+total+1,iatom+total+2,tortype)
         endif
      elseif ( termflag .eq. 1 ) then
c        --- C-terminal proline
         write(6,*) 'UA_PEPNAME_P0: lacking a set '
         write(6,*) 'of parameters for N-terminal proline'
         lfailure = .true.
         return
      endif

      return
      end

      subroutine twh_ua_pepname_q0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'q0' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype

      aminoname = 'glutamine'
      shortname = 'GLN'
c     --- glutamine side chain
      heavy = 5
      shift = hmain+heavy+ninsert
      light = 2
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CD')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' OE1')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' NE2')
      call twh_pdbname(GLB_SET,imolty,shift,'1HE2')
      call twh_pdbname(GLB_SET,imolty,shift+1,'2HE2')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,2)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+1,2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_invib(GLB_SET,imolty,hmain+2,3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+2,3,hmain+4)
      call twh_invib(GLB_SET,imolty,hmain+3,1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+2)
      call twh_invib(GLB_SET,imolty,hmain+4,3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,shift)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,shift+1)
c     --- hydrogens
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain+4)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain+4)
      if ( lsettorsion ) then
c        --- side chain torsions
c        --- N - Ca - Cb - Cg
         call twh_create_torsion(imolty,iatom,iatom+1,hmain
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
c        --- Ca - Cb - Cg - Cd
         call twh_create_torsion(imolty,iatom+1,hmain,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Cd - Cg - Cb - Ca
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain
     &        ,iatom+1,tortype)
c        --- Cb - Cg - Cd - Oe1
         call twh_create_torsion(imolty,hmain,hmain+1,hmain+2
     &        ,hmain+3,tortype)
c        --- reverse: Oe1 - Cd - Cg - Cb
         call twh_create_torsion(imolty,hmain+3,hmain+2,hmain+1
     &        ,hmain,tortype)
c        --- Cg - Cd - Ne2 - He21
         call twh_create_torsion(imolty,hmain+1,hmain+2,hmain+4
     &        ,shift,tortype)
c        --- reverse: He21 - Ne2 - Cd - Cg
         call twh_create_torsion(imolty,shift,hmain+4,hmain+2
     &        ,hmain+1,tortype)
      endif

      return
      end

      subroutine twh_ua_pepname_r_plus(lfailure,aminoname,shortname
     &     ,imolty,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'r+' *
c     * monomer in the polypeptide builder                             *
c     * r+: arginine protonated                                        *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype

      aminoname = 'arginine protonated'
      shortname = 'ARG'
c     --- arginine side chain
      heavy = 7
      shift = hmain+heavy+ninsert
      light = 5
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CD')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' NE')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' CZ')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' NH1')
      call twh_pdbname(GLB_SET,imolty,hmain+6,' NH2')
      call twh_pdbname(GLB_SET,imolty,shift,' HE')
      call twh_pdbname(GLB_SET,imolty,shift+1,'1HH1')
      call twh_pdbname(GLB_SET,imolty,shift+2,'2HH1')
      call twh_pdbname(GLB_SET,imolty,shift+3,'1HH2')
      call twh_pdbname(GLB_SET,imolty,shift+4,'2HH2')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,2)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+1,2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_invib(GLB_SET,imolty,hmain+2,2)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+3,3)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,shift)
      call twh_invib(GLB_SET,imolty,hmain+4,3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,hmain+6)
      call twh_invib(GLB_SET,imolty,hmain+5,3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+5,2,shift+1)
      call twh_ijvib(GLB_SET,imolty,hmain+5,3,shift+2)
      call twh_invib(GLB_SET,imolty,hmain+6,3)
      call twh_ijvib(GLB_SET,imolty,hmain+6,1,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+6,2,shift+3)
      call twh_ijvib(GLB_SET,imolty,hmain+6,3,shift+4)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain+3)
      call twh_invib(GLB_SET,imolty,shift+1,1)
      call twh_ijvib(GLB_SET,imolty,shift+1,1,hmain+5)
      call twh_invib(GLB_SET,imolty,shift+2,1)
      call twh_ijvib(GLB_SET,imolty,shift+2,1,hmain+5)
      call twh_invib(GLB_SET,imolty,shift+3,1)
      call twh_ijvib(GLB_SET,imolty,shift+3,1,hmain+6)
      call twh_invib(GLB_SET,imolty,shift+4,1)
      call twh_ijvib(GLB_SET,imolty,shift+4,1,hmain+6)
      if ( lsettorsion ) then
c        --- side chain torsions
c        --- N - Ca - Cb - Cg
         call twh_create_torsion(imolty,iatom,iatom+1,hmain,hmain+1
     &        ,tortype)
c        --- reverse: Cg - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1,iatom
     &        ,tortype)
c        --- Ca - Cb - Cg - Cd
         call twh_create_torsion(imolty,iatom+1,hmain,hmain+1,hmain+2
     &        ,tortype)
c        --- reverse: Cd - Cg - Cb - Ca
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain,iatom+1
     &        ,tortype)
c        --- Cb - Cg - Cd - Ne
         call twh_create_torsion(imolty,hmain,hmain+1,hmain+2,hmain+3
     &        ,tortype)
c        --- reverse: Ne - Cd - Cg - Cb
         call twh_create_torsion(imolty,hmain+3,hmain+2,hmain+1
     &        ,hmain,tortype)
c        --- Cg - Cd - Ne - Cz
         call twh_create_torsion(imolty,hmain+1,hmain+2,hmain+3
     &        ,hmain+4,tortype)
c        --- reverse: Cz - Ne - Cd - Cg
         call twh_create_torsion(imolty,hmain+4,hmain+3,hmain+2
     &        ,hmain+1,tortype)
c        --- Cd - Ne - Cz - Nh1
         call twh_create_torsion(imolty,hmain+2,hmain+3,hmain+4
     &        ,hmain+5,tortype)
c        --- reverse: Nh1 - Cz - Ne - Cd
         call twh_create_torsion(imolty,hmain+5,hmain+4,hmain+3
     &        ,hmain+2,tortype)
c        --- Ne - Cz - Nh1 - Hh11
         call twh_create_torsion(imolty,hmain+3,hmain+4,hmain+5
     &        ,shift+1,tortype)
c        --- reverse: Hh11 - Nh1 - Cz - Ne
         call twh_create_torsion(imolty,shift+1,hmain+5,hmain+4
     &        ,hmain+3,tortype)
c        --- Ne - Cz - Nh1 - Hh12
         call twh_create_torsion(imolty,hmain+3,hmain+4,hmain+5
     &        ,shift+2,tortype)
c        --- reverse: Hh12 - Nh1 - Cz - Ne
         call twh_create_torsion(imolty,shift+2,hmain+5,hmain+4
     &        ,hmain+3,tortype)
c        --- Ne - Cz - Nh2 - Hh21
         call twh_create_torsion(imolty,hmain+3,hmain+4,hmain+6
     &        ,shift+3,tortype)
c        --- reverse: Hh21 - Nh2 - Cz - Ne
         call twh_create_torsion(imolty,shift+3,hmain+6,hmain+4
     &        ,hmain+3,tortype)
c        --- Ne - Cz - Nh2 - Hh22
         call twh_create_torsion(imolty,hmain+3,hmain+4,hmain+6
     &        ,shift+4,tortype)
c        --- reverse: Hh22 - Nh2 - Cz - Ne
         call twh_create_torsion(imolty,shift+4,hmain+6,hmain+4
     &        ,hmain+3,tortype)
      endif

      return
      end

      subroutine twh_ua_pepname_s0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 's0' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype

      aminoname = 'serine'
      shortname = 'SER'
c     --- serine side chain
      heavy = 2
      shift = hmain+heavy+ninsert
      light = 1
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' OG')
      call twh_pdbname(GLB_SET,imolty,shift,' HG')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,2)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+1,2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,shift)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain+1)
      if ( lsettorsion ) then
c        --- side chain torsions
c        --- N - Ca - Cb - Og
         call twh_create_torsion(imolty,iatom,iatom+1,hmain
     &        ,hmain+1,tortype)
c        --- reverse: Og - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
c        --- Ca - Cb - Og - Hg
         call twh_create_torsion(imolty,iatom+1,hmain,hmain+1
     &        ,shift,tortype)
c        --- reverse: Hg - Og - Cb - Ca
         call twh_create_torsion(imolty,shift,hmain+1,hmain
     &        ,iatom+1,tortype)
      endif

      return
      end

      subroutine twh_ua_pepname_t0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 't0' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype

      aminoname = 'threonine'
      shortname = 'THR'
c     --- threonine side chain
      heavy = 3
      shift = hmain+heavy+ninsert
      light = 1
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' OG1')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CG2')
      call twh_pdbname(GLB_SET,imolty,shift,' HG1')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,3)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,hmain+2)
      call twh_invib(GLB_SET,imolty,hmain+1,2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,shift)
      call twh_invib(GLB_SET,imolty,hmain+2,1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain+1)
      if ( lsettorsion ) then
c        --- side chain torsions
c        --- N - Ca - Cb - Og1
         call twh_create_torsion(imolty,iatom,iatom+1,hmain
     &        ,hmain+1,tortype)
c        --- reverse: Og1 - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
c        --- Ca - Cb - Og1 - Hg1
         call twh_create_torsion(imolty,iatom+1,hmain,hmain+1
     &        ,shift,tortype)
c        --- reverse: Hg1 - Og1 - Cb - Ca
         call twh_create_torsion(imolty,shift,hmain+1,hmain
     &        ,iatom+1,tortype)
      endif

      return
      end

      subroutine twh_ua_pepname_v0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'v0' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype

      aminoname = 'valine'
      shortname = 'VAL'
c     --- valine side chain
      heavy = 3
      shift = hmain+heavy+ninsert
      light = 0
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG1')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CG2')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,3)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain,3,hmain+2)
      call twh_invib(GLB_SET,imolty,hmain+1,1)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_invib(GLB_SET,imolty,hmain+2,1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain)
      if ( lsettorsion ) then
c        --- side chain torsions
c        --- N - Ca - Cb - Cg1
         call twh_create_torsion(imolty,iatom,iatom+1,hmain
     &        ,hmain+1,tortype)
c        --- reverse: Cg1 - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
      endif
      
      return
      end

      subroutine twh_ua_pepname_w0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'w0' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype

      aminoname = 'tryptophan'
      shortname = 'TRP'
c     --- tryptophan side chain
      heavy = 10
      shift = hmain+heavy+ninsert
      light = 1
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CD1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CD2')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' NE1')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' CE2')
      call twh_pdbname(GLB_SET,imolty,hmain+6,' CE3')
      call twh_pdbname(GLB_SET,imolty,hmain+7,' CZ2')
      call twh_pdbname(GLB_SET,imolty,hmain+8,' CZ3')
      call twh_pdbname(GLB_SET,imolty,hmain+9,' CH2')
      call twh_pdbname(GLB_SET,imolty,shift,' HE1')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,2)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,2)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+4)
      call twh_invib(GLB_SET,imolty,hmain+3,3)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+3,3,hmain+6)
      call twh_invib(GLB_SET,imolty,hmain+4,3)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+4,3,shift)
      call twh_invib(GLB_SET,imolty,hmain+5,3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,2,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+5,3,hmain+7)
      call twh_invib(GLB_SET,imolty,hmain+6,2)
      call twh_ijvib(GLB_SET,imolty,hmain+6,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+6,2,hmain+8)
      call twh_invib(GLB_SET,imolty,hmain+7,2)
      call twh_ijvib(GLB_SET,imolty,hmain+7,1,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+7,2,hmain+9)
      call twh_invib(GLB_SET,imolty,hmain+8,2)
      call twh_ijvib(GLB_SET,imolty,hmain+8,1,hmain+6)
      call twh_ijvib(GLB_SET,imolty,hmain+8,2,hmain+9)
      call twh_invib(GLB_SET,imolty,hmain+9,2)
      call twh_ijvib(GLB_SET,imolty,hmain+9,1,hmain+7)
      call twh_ijvib(GLB_SET,imolty,hmain+9,2,hmain+8)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain+4)
      if ( lsettorsion ) then
c        --- could not implement all of the torsions as two of them
c        --- are this wild ring spanning business between atoms that
c        --- are not bonded together.  The following Charmm19 torsion
c        --- were not implemented
c        --- Cd1 - Ce2 - Cd2 - Cz3
c        --- Cd1 - Cd2 - Ce2 - Ch2
c        --- side chain torsions
c        --- N - Ca - Cb - Cg
         call twh_create_torsion(imolty,iatom,iatom+1,hmain
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
c        --- Ca - Cb - Cg - Cd1
         call twh_create_torsion(imolty,iatom+1,hmain,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Cd1 - Cg - Cb - Ca
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain
     &        ,iatom+1,tortype)
c        --- Cz2 - Ce2 - Cd2 - Cg
         call twh_create_torsion(imolty,hmain+7,hmain+5,hmain+3
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cd2 - Ce2 - Cz2
         call twh_create_torsion(imolty,hmain+1,hmain+3,hmain+5
     &        ,hmain+7,tortype)
c        --- Ce3 - Cd2 - Ce2 - Ne1
         call twh_create_torsion(imolty,hmain+6,hmain+3,hmain+5
     &        ,hmain+4,tortype)
c        --- reverse: Ne1 - Ce2 - Cd2 - Ce3
         call twh_create_torsion(imolty,hmain+4,hmain+5,hmain+3
     &        ,hmain+6,tortype)
c        --- side chain torsions listed as impropers in charmm
c        --- Cd2 - Ce2 - Cz2 - Ch2
         call twh_create_torsion(imolty,hmain+3,hmain+5,hmain+7
     &        ,hmain+9,tortype)
c        --- reverse: Ch2 - Cz2 - Ce2 - Cd2
         call twh_create_torsion(imolty,hmain+9,hmain+7,hmain+5
     &        ,hmain+3,tortype)
c        --- Ce2 - Cz2 - Ch2 - Cz3
         call twh_create_torsion(imolty,hmain+5,hmain+7,hmain+9
     &        ,hmain+8,tortype)
c        --- reverse: Cz3 - Ch2 - Cz2 - Ce2
         call twh_create_torsion(imolty,hmain+8,hmain+9,hmain+7
     &        ,hmain+5,tortype)
c        --- Cz2 - Ch2 - Cz3 - Ce3
         call twh_create_torsion(imolty,hmain+7,hmain+9,hmain+8
     &        ,hmain+6,tortype)
c        --- reverse: Ce3 - Cz3 - Ch2 - Cz2
         call twh_create_torsion(imolty,hmain+6,hmain+8,hmain+9
     &        ,hmain+7,tortype)
c        --- Ch2 - Cz3 - Ce3 - Cd2
         call twh_create_torsion(imolty,hmain+9,hmain+8,hmain+6
     &        ,hmain+3,tortype)
c        --- reverse: Cd2 - Ce3 - Cz3 - Ch2
         call twh_create_torsion(imolty,hmain+3,hmain+6,hmain+8
     &        ,hmain+9,tortype)
c        --- Cz3 - Ce3 - Cd2 - Ce2
         call twh_create_torsion(imolty,hmain+8,hmain+6,hmain+3
     &        ,hmain+5,tortype)
c        --- reverse: Ce2 - Cd2 - Ce3 - Cz3
         call twh_create_torsion(imolty,hmain+5,hmain+3,hmain+6
     &        ,hmain+8,tortype)
c        --- Ce3 - Cd2 - Ce2 - Cz2
         call twh_create_torsion(imolty,hmain+6,hmain+3,hmain+5
     &        ,hmain+7,tortype)
c        --- reverse: Cz2 - Ce2 - Cd2 - Ce3
         call twh_create_torsion(imolty,hmain+7,hmain+5,hmain+3
     &        ,hmain+6,tortype)
c        --- Cg - Cd1 - Ne1 - Ce2
         call twh_create_torsion(imolty,hmain+1,hmain+2,hmain+4
     &        ,hmain+5,tortype)
c        --- reverse: Ce2 - Ne1 - Cd1 - Cg
         call twh_create_torsion(imolty,hmain+5,hmain+4,hmain+2
     &        ,hmain+1,tortype)
c        --- Cd1 - Ne1 - Ce2 - Cd2
         call twh_create_torsion(imolty,hmain+2,hmain+4,hmain+5
     &        ,hmain+3,tortype)
c        --- reverse: Cd2 - Ce2 - Ne1 - Cd1
         call twh_create_torsion(imolty,hmain+3,hmain+5,hmain+4
     &        ,hmain+2,tortype)
c        --- Ne1 - Ce2 - Cd2 - Cg
         call twh_create_torsion(imolty,hmain+4,hmain+5,hmain+3
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cd2 - Ce2 - Ne1
         call twh_create_torsion(imolty,hmain+1,hmain+3,hmain+5
     &        ,hmain+4,tortype)
c        --- Ce2 - Cd2 - Cg - Cd1
         call twh_create_torsion(imolty,hmain+5,hmain+3,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Cd1 - Cg - Cd2 - Ce2
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain+3
     &        ,hmain+5,tortype)
c        --- Cd2 - Cg - Cd1 - Ne1
         call twh_create_torsion(imolty,hmain+3,hmain+1,hmain+2
     &        ,hmain+4,tortype)
c        --- reverse: Ne1 - Cd1 - Cg - Cd2
         call twh_create_torsion(imolty,hmain+4,hmain+2,hmain+1
     &        ,hmain+3,tortype)
      endif

      return
      end

      subroutine twh_ua_pepname_y0(lfailure,aminoname,shortname,imolty
     &     ,iatom,hmain,heavy,light,shift,ninsert,tortype
     &     ,lsettorsion)
c     ******************************************************************
c     * sets up the communal bond information for the united-atom 'y0' *
c     * monomer in the polypeptide builder                             *
c     *                                                                *
c     * split from polycharmm19 03-27-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lsettorsion
      character*3 shortname
      character*40 aminoname
      integer imolty,iatom,hmain,heavy,light,shift,ninsert,tortype

      aminoname = 'tyrosine'
      shortname = 'TYR'
c     --- tyrosine side chain
      heavy = 8
      shift = hmain+heavy+ninsert
      light = 1
      call twh_monomer_check(lfailure,shift+light)
      if ( lfailure ) return
c     --- pdb names
      call twh_pdbname(GLB_SET,imolty,hmain,' CB')
      call twh_pdbname(GLB_SET,imolty,hmain+1,' CG')
      call twh_pdbname(GLB_SET,imolty,hmain+2,' CD1')
      call twh_pdbname(GLB_SET,imolty,hmain+3,' CD2')
      call twh_pdbname(GLB_SET,imolty,hmain+4,' CE1')
      call twh_pdbname(GLB_SET,imolty,hmain+5,' CE2')
      call twh_pdbname(GLB_SET,imolty,hmain+6,' CZ')
      call twh_pdbname(GLB_SET,imolty,hmain+7,' OH')
      call twh_pdbname(GLB_SET,imolty,shift,' HH')
c     --- vibrations
      call twh_invib(GLB_SET,imolty,hmain,2)
      call twh_ijvib(GLB_SET,imolty,hmain,1,iatom+1)
      call twh_ijvib(GLB_SET,imolty,hmain,2,hmain+1)
      call twh_invib(GLB_SET,imolty,hmain+1,3)
      call twh_ijvib(GLB_SET,imolty,hmain+1,1,hmain)
      call twh_ijvib(GLB_SET,imolty,hmain+1,2,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+1,3,hmain+3)
      call twh_invib(GLB_SET,imolty,hmain+2,2)
      call twh_ijvib(GLB_SET,imolty,hmain+2,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+2,2,hmain+4)
      call twh_invib(GLB_SET,imolty,hmain+3,2)
      call twh_ijvib(GLB_SET,imolty,hmain+3,1,hmain+1)
      call twh_ijvib(GLB_SET,imolty,hmain+3,2,hmain+5)
      call twh_invib(GLB_SET,imolty,hmain+4,2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,1,hmain+2)
      call twh_ijvib(GLB_SET,imolty,hmain+4,2,hmain+6)
      call twh_invib(GLB_SET,imolty,hmain+5,2)
      call twh_ijvib(GLB_SET,imolty,hmain+5,1,hmain+3)
      call twh_ijvib(GLB_SET,imolty,hmain+5,2,hmain+6)
      call twh_invib(GLB_SET,imolty,hmain+6,3)
      call twh_ijvib(GLB_SET,imolty,hmain+6,1,hmain+4)
      call twh_ijvib(GLB_SET,imolty,hmain+6,2,hmain+5)
      call twh_ijvib(GLB_SET,imolty,hmain+6,3,hmain+7)
      call twh_invib(GLB_SET,imolty,hmain+7,2)
      call twh_ijvib(GLB_SET,imolty,hmain+7,1,hmain+6)
      call twh_ijvib(GLB_SET,imolty,hmain+7,2,shift)
      call twh_invib(GLB_SET,imolty,shift,1)
      call twh_ijvib(GLB_SET,imolty,shift,1,hmain+7)
      if ( lsettorsion ) then
c        --- side chain torsions
c        --- N - Ca - Cb - Cg
         call twh_create_torsion(imolty,iatom,iatom+1,hmain
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cb - Ca - N
         call twh_create_torsion(imolty,hmain+1,hmain,iatom+1
     &        ,iatom,tortype)
c        --- Ca - Cb - Cg - Cd1
         call twh_create_torsion(imolty,iatom+1,hmain,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Cd1 - Cg - Cb - Ca
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain
     &        ,iatom+1,tortype)
c        --- Ce2 - Cz - Oh - Hh
         call twh_create_torsion(imolty,hmain+5,hmain+6,hmain+7
     &     ,shift,tortype)
c        --- reverse: Hh - Oh - Cz - Ce2
         call twh_create_torsion(imolty,shift,hmain+7,hmain+6
     &        ,hmain+5,tortype)
c        --- side chain torsions listed in charmm as impropers
c        --- Cg - Cd1 - Ce1 - Cz
         call twh_create_torsion(imolty,hmain+1,hmain+2,hmain+4
     &        ,hmain+6,tortype)
c        --- reverse: Cz - Ce1 - Cd1 - Cg
         call twh_create_torsion(imolty,hmain+6,hmain+4,hmain+2
     &        ,hmain+1,tortype)
c        --- Cd1 - Ce1 - Cz - Ce2
         call twh_create_torsion(imolty,hmain+2,hmain+4,hmain+6
     &        ,hmain+5,tortype)
c        --- reverse: Ce2 - Cz - Ce1 - Cd1
         call twh_create_torsion(imolty,hmain+5,hmain+6,hmain+4
     &        ,hmain+2,tortype)
c        --- Ce1 - Cz - Ce2 - Cd2
         call twh_create_torsion(imolty,hmain+4,hmain+6,hmain+5
     &        ,hmain+3,tortype)
c        --- reverse: Cd2 - Ce2 - Cz - Ce1
         call twh_create_torsion(imolty,hmain+3,hmain+5,hmain+6
     &        ,hmain+4,tortype)
c        --- Cz - Ce2 - Cd2 - Cg
         call twh_create_torsion(imolty,hmain+6,hmain+5,hmain+3
     &        ,hmain+1,tortype)
c        --- reverse: Cg - Cd2 - Ce2 - Cz
         call twh_create_torsion(imolty,hmain+1,hmain+3,hmain+5
     &        ,hmain+6,tortype)
c        --- Ce2 - Cd2 - Cg - Cd1
         call twh_create_torsion(imolty,hmain+5,hmain+3,hmain+1
     &        ,hmain+2,tortype)
c        --- reverse: Cd1 - Cg - Cd2 - Ce2
         call twh_create_torsion(imolty,hmain+2,hmain+1,hmain+3
     &        ,hmain+5,tortype)
c        --- Cd2 - Cg - Cd1 - Ce1
         call twh_create_torsion(imolty,hmain+3,hmain+1,hmain+2
     &        ,hmain+4,tortype)
c        --- reverse: Ce1 - Cd1 - Cg - Cd2
         call twh_create_torsion(imolty,hmain+4,hmain+2,hmain+1
     &        ,hmain+3,tortype)
      endif

      return
      end

      subroutine twh_monomer_check(lfailure,value)
c     ******************************************************************
c     * determines whether there is enough room left in the arrays for *
c     * the atoms required by this monomer                             *
c     *                                                                *
c     * originally written 06-28-2006 by M.G. Martin                   *
c     * last modified 06-28-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer value

      if ( value .gt. NUMAX ) then
         write(6,*) 'MONOMERS: error - request for atoms exceeds NUMAX'
         write(6,*) 'requested atom total:',value
         write(6,*) 'NUMAX:',NUMAX
         lfailure = .true.
      endif
      return
      end
