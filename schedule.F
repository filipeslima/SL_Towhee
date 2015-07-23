#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2009 Marcus G. Martin                       *
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
      subroutine twh_schedule(lfailure,numunits,imolty,index,iutry
     &     ,moveinfo,iprev,movetype)
c     ******************************************************************
c     * computes the growth shedule for CBMC type moves                *
c     * this is the heart of the logic that performs the CB moves      *
c     * You will never find a more wretched hive of scum and villany.  *
c     * we must be cautious ...                                        *
c     *                                                                *
c     * numunits: the number of units in the molecule (passed in)      *
c     * imolty: the molecule type (passed in)                          *
c     * index: the number of cb growth steps (computed and passed back)*
c     * iutry: the bead from which the growth begins IF this is a      *
c     *        cbswap move, in which case it is (passed in) and also   *
c     *        used again later in those routines so don't change it.  *
c     *        if this is a backbone move then it is not passed in, but*
c     *        is passed back out with the value of growfrom(1) for    *
c     *        some bookkeeping in cbregrow                            *
c     * moveinfo: movetype == 0 this is the total number of atoms that *
c     *          would be grown using this growth path (passed out)    *
c     *          movetype == 4 this is the torsion number used for the *
c     *          backbone regrowth move (passed out)                   *
c     *          movetype ==6 this is the aminonum of the key CA       *
c     *          (passed out)                                          *
c     * iprev:   if movetype= 0 then use this as the vibration number  *
c     *          of the existing neighbor and compute the growth path  *
c     *          (passed in)                                           *
c     *          if movetype= 3 (passed in) to tell schedule which atom*
c     *          is the previous atom in the growth.  Only valid for   *
c     *          the (currently unsupported) cbswitch move             *
c     * movetype: tells schedule what kind of cb move we are doing     *
c     *        0: determine the number of atoms that would be grown    *
c     *           for a given iutry using iprev as the previous bead   *
c     *        1: standard cbregrow which selects a bead and then      *
c     *           either regrows the entire molecule, or regrows in all*
c     *           directions except for one from the chosen bead       *
c     *        2: cbswap move, starts from bead iutry and grows in all *
c     *           directions                                           *
c     *        3: cbswitch move, starts from bead iutry and grows in   *
c     *           all directions, except for those atoms in direction  *
c     *           iprev                                                *
c     *        4: backbone cbregrowth move.  currently only works for  *
c     *           polypeptides.  Selects an atom that is on the protein*
c     *           backbone, selects another backbone atom that is      *
c     *           connected to the first by a torsion, and then regrows*
c     *           the atoms inbetween                                  *
c     *        5: regrows all atoms except for a list that remain in   *
c     *           the same positions. Originally designed to construct *
c     *           DNA helices from the P positions                     *
c     *        6: regrows a particular sidechain, or occasionally      *
c     *           multiple bonded sidechains, on a polypeptide         *
c     *                                                                *
c     * DECODER for the growth logic                                   *
c     * lexsched is true if the bead exists at that time of the growth *
c     *     it is false when a bead has not yet been grown this time   *
c     * growfrom(index) gives the unit number that is to be grown from *
c     *     for the index step                                         *
c     * growprev(index) gives the bead that exists connected to        *
c     *     growfrom(index)                                            *
c     * grownum(index) gives the number of beads to be grown from      *
c     *     growfrom(index)                                            *
c     * growlist(index,count) gives the unit numbers for the beads     *
c     *     grown from growfrom(index), get count from grownum(index)  *
c     * grownbon(index,count) is the number of atoms that are being    *
c     *     grown this step, or which already exist, which are bonded  *
c     *     to the atom growlist(index,count).  This is used for fecmbc*
c     * growbond(index,count,num) is the list of atoms which are bonded*
c     *     to growlist(index,count) and already exist, or are being   *
c     *     grown during this step.  This is used for fecmbc           *
c     * growxtra(index,count+1) is the list of atoms which are bonded  *
c     *     to growfrom(index) and are not being grown this step, and  *
c     *     are not the growprev(index).  This is used for fecbmc      *
c     *     growxtra(index,1) is the number of xtra atoms              *
c     *                                                                *
c     * originally taken from the old config and swap subroutines and  *
c     * modified by  M.G. Martin 1-29-1998                             *
c     * last modified 06-24-2009 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_ELEMENTNAME
#define FUNCTION_GET_GROWCOUNT
#define FUNCTION_GET_GROWVALIDCOUNT
#define FUNCTION_GET_GROWVALIDLIST
#define FUNCTION_GET_IJTOR2
#define FUNCTION_GET_IJTOR3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_IJVIB
#define FUNCTION_GET_INTOR
#define FUNCTION_GET_INVIB
#define FUNCTION_GET_NTYPE
#define FUNCTION_RANDOM
#include "functions.h"
#define DEBUG_SCHEDULE 0
c     --- variables passed to/from the subroutine
      logical lfailure
      integer numunits,imolty,index,iutry,iprev,movetype,moveinfo
c     --- local variables
c     --- character string scalars
      character*3 aminoshort
      character*4 pdbname
c     --- logical scalars
      logical lfound,lhook,ldone,lbacka,lmida,lmidb,lbackb,lvalid
      logical lbackbone,lbacktwo,lbackthr,lexschedule
c     --- logical arrays
      logical lstart
      dimension lstart(NUMAX)
c     --- integer scalars
      integer icbu,iut,invtry,iu,ibead,count,ivib,itry,atomptr,ixtra
      integer temp_count,zz,outer_num,iufrom,outer_try,itemp
      integer ind_start,ibackbone,iubacka,iubackb,iumida,iumidb
      integer itor,ind_other,iother,vibcount,isearch,ifound,icurrent
      integer nca,iunit,target,test_index,junit,nvibs,black_index
      integer nunit,testring,growfrom,growprev,growtest
      integer grownuma,grownumb,grownbon,xtracount,nmaxcbmc
c     --- integer arrays
      integer viblist
      dimension viblist(NNBOND)
      integer tempback
      dimension tempback(MAXTOR)
      integer temp_store
      dimension temp_store(NNBOND)
      integer outer_sites,outer_prev
      dimension outer_sites(NUMAX),outer_prev(NUMAX)
      integer blacklist
      dimension blacklist(NUMAX)
c     --- double precision scalars
      double precision rvalue,pmall

c      write(6,*) 'start SCHEDULE'

c     --- initialize lexsched, lstart
c     --- initialize grownbon 
      call twh_grownbon(GLB_INIT,numunits,NNBOND,0)
      do iu = 1,numunits
         call twh_lexsched(GLB_SET,iu,.true.)
         lstart(iu) = .false.
         call twh_growxtra(GLB_SET,iu,1,0)
      enddo

c     --- retrieve some constants
      call twh_nunit(GLB_GET,imolty,nunit)

c     --- initialize temp_count
      temp_count = 0
      if ( movetype .eq. 0 ) then
c        --- compute the growth pattern and return the number of atoms
c        --- that would be grown
         call twh_inring(GLB_GET,imolty,iutry,testring)
         if ( testring .ne. 0 ) then
c           --- using current logic, entire molecule regrown
            moveinfo = numunits - 1
            return
         endif
         call twh_growfrom(GLB_SET,1,iutry)
         lstart(iutry) = .true.
         index = 1
         invtry = twh_get_invib(imolty,iutry)
c        --- safety check
         if ( invtry .lt. iprev ) then
            write(6,*) 'SCHEDULE: movetype 0 passed in an iprev'
            write(6,*) 'that is larger than invib',iprev,invtry
            lfailure = .true.
            return
         endif
         if ( invtry .eq. 1 ) then
c           --- not enough neighbors so entire molecule will be regrown 
            moveinfo = numunits - 1
            return
         endif
         iut = twh_get_ijvib(imolty,iutry,iprev)
         count = 0
         do ivib = 1,invtry
            if ( ivib .ne. iprev ) then 
               iut = twh_get_ijvib(imolty,iutry,ivib)
c              --- add this branch to the list of those to grow 
               temp_count = temp_count + 1
               temp_store(temp_count) = iut
            endif
         enddo

c        --- put these atoms into the data structure
         do zz = temp_count,1,-1
c           --- choose grow bead randomly from temp_store
            rvalue = dble(zz)* twh_random()
            itry = idint( rvalue ) + 1
            iut = temp_store(itry)
            count = count + 1
            call twh_growlist(GLB_SET,1,count,iut)
c           --- update temp_store for next iteration
            temp_store(itry) = temp_store(zz)
         enddo
         call twh_grownum(GLB_SET,1,count)

c        ---end the part that is specific for movetype 0

      elseif ( movetype .eq. 1 ) then 
c        -- select the first bead to grow from
c        --- choose an atom at random from the growvalidlist
         rvalue = twh_random() *dble(twh_get_growvalidcount(imolty))
         atomptr = idint(rvalue) + 1
         iutry = twh_get_growvalidlist(imolty,atomptr)
         call twh_growfrom(GLB_SET,1,iutry)
         lstart(iutry) = .true.
         invtry = twh_get_invib(imolty,iutry)
         index = 1
         if ( invtry .eq. 0 ) then
c           --- problem, cannot do cbregrow move on a 1 bead molecule
            write(6,*) 'SCHEDULE'
            write(6,*) 'molecule type ',imolty,' only has 1 atom.'
            write(6,*) 'cannot do CBMC regrow for one unit molecule'
            lfailure = .true.
            return
         elseif ( invtry .eq. 1 ) then
c           --- regrow entire molecule
c           --- set up the grow data structures
            call twh_growprev(GLB_SET,1,0)
            call twh_grownum(GLB_SET,1,invtry)
            ivib = invtry
            iut = twh_get_ijvib(imolty,iutry,ivib)
            call twh_growlist(GLB_SET,1,ivib,iut)
         else
c           --- at a branch point, decide how many branches to regrow
c           --- regrow all of legal branches 
            count = 0
            do ivib = 1,invtry
               iut = twh_get_ijvib(imolty,iutry,ivib)
c              --- add this branch to the list of those to grow 
               temp_count = temp_count + 1
               temp_store(temp_count) = iut
            enddo

c           --- put these atoms into the data structure
            do zz = temp_count,1,-1
c              --- choose grow bead randomly from temp_store
               rvalue = dble(zz)* twh_random()
               itry = idint( rvalue ) + 1
               iut = temp_store(itry)
               count = count + 1
               call twh_growlist(GLB_SET,1,count,iut)
c              --- update temp_store for next iteration
               temp_store(itry) = temp_store(zz)
            enddo

            call twh_grownum(GLB_SET,1,count)
            if ( count .ne. invtry ) then
               write(6,*) 'SCHEDULE: count .ne. invtry'
               lfailure = .true.
               return
            endif
c           --- determine whether we really wanted to regrow in all 
c           --- directions or in all directions except for one
            call twh_inring(GLB_GET,imolty,iutry,testring)
            call twh_pmall(GLB_GET,imolty,pmall)
            if ( testring .eq. 0 .and. twh_random() .gt. pmall ) then
c              --- not pmall so select a previous bead
c              --- initialize viblist
               do ivib = 1,invtry
                  viblist(ivib) = ivib
               enddo
               vibcount = invtry
c              --- select a valid atom from the list (check nmaxcbmc)
               lfound = .false.
               do while ( .not. lfound .and. vibcount .ne. 0 )
                  rvalue = twh_random()*dble(vibcount)
                  atomptr = idint( rvalue ) + 1
                  ivib = viblist(atomptr)
                  call twh_nmaxcbmc(GLB_GET,imolty,nmaxcbmc)
                  if ( twh_get_growcount(imolty,iutry,ivib) 
     &                 .le. nmaxcbmc ) then
c                    --- found a valid target
                     lfound = .true.
                     iut = twh_get_ijvib(imolty,iutry,ivib)
                  else
c                    --- try again
                     viblist(atomptr) = viblist(vibcount)
                     vibcount = vibcount - 1
                  endif
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'SCHEDULE: neighbor atom not found'
                  lfailure = .true.
                  return
               endif
c              --- find the atom we want
               lfound = .false.
               ivib = 0
               do while ( .not. lfound .and. ivib .lt. invtry )
                  ivib = ivib + 1
                  call twh_growlist(GLB_GET,1,ivib,iu)
                  if ( iu .eq. iut ) lfound = .true.
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'SCHEDULE: vibration not found'
                  lfailure = .true.
                  return
               endif
               call twh_growprev(GLB_SET,1,iu)
c              --- replace this unit with the last unit in growlist
               call twh_growlist(GLB_GET,1,invtry,growtest)
               call twh_growlist(GLB_SET,1,ivib,growtest)
c              --- reduce numgrow by 1
               call twh_grownum(GLB_DECR,1,1)
            else
c              --- we regrew all branches, no previous bead
               call twh_growprev(GLB_SET,1,0)
            endif
         endif   
c        ---end the part that is specific for cbregrow

      elseif ( movetype .eq. 2 ) then
c        --- begin the part that is specific for cbswap
c        --- iutry is the first bead inserted - need to grow its 
c        ---- neighbors
         call twh_growfrom(GLB_SET,1,iutry)
         lstart(iutry) = .true.

         invtry = twh_get_invib(imolty,iutry)
         call twh_growprev(GLB_SET,1,0)
         if ( invtry .eq. 0 ) then
c           --- Bead iutry is the only bead to be grown
            index = 0
         else
c           --- grow all of the beads connected to bead iutry
            index = 1
            call twh_grownum(GLB_SET,1,invtry)
            do ivib = 1,invtry
               iut = twh_get_ijvib(imolty,iutry,ivib)
c              --- grow these branches
               temp_count = temp_count + 1
               temp_store(temp_count) = iut
            enddo

            count = 0
            do zz = temp_count,1,-1
c              --- choose grow bead randomly from temp_store
               rvalue = dble(zz)* twh_random()
               itry = idint( rvalue ) + 1
               iut = temp_store(itry)
               count = count + 1
               call twh_growlist(GLB_SET,1,count,iut)
c              --- update temp_store for next iteration
               temp_store(itry) = temp_store(zz)
            enddo
         endif
c        ---end the part that is specific for cbswap

      elseif ( movetype .eq. 3 ) then
c        --- begin part that is specific for cbswitch 
         if ( iutry .eq. 0 ) then
c           --- no beads to be regrown via cbmc
            index = 0
            return
         endif
         call twh_growfrom(GLB_SET,1,iutry)
         lstart(iutry) = .true.
         invtry = twh_get_invib(imolty,iutry)
         call twh_growprev(GLB_SET,1,iprev)

         if ( invtry .eq. 0 ) then
c           --- Bead 1 is the only bead to be grown
            index = 0
         else

c           --- grow all (except iprev) of the beads connected to bead 1
            index = 1
            do ivib = 1,invtry
               iut = twh_get_ijvib(imolty,iutry,ivib)
               if ( iut .ne. iprev ) then
                  temp_count = temp_count + 1
                  temp_store(temp_count) = iut
               endif
            enddo

            count = 0
            do zz = temp_count,1,-1
c              --- choose grow bead randomly from temp_store
               rvalue = dble(zz)* twh_random()
               itry = idint( rvalue ) + 1
               iut = temp_store(itry)
               count = count + 1
               call twh_growlist(GLB_SET,1,count,iut)
c              --- update temp_store for next iteration
               temp_store(itry) = temp_store(zz)
            enddo
            call twh_grownum(GLB_SET,1,count)
         endif
c        --- end part that is specific for cbswitch
      elseif ( movetype .eq. 4 ) then
c        --- begin part that is specific for backbone cbregrow
c        --- select a backbone atom at random
         call twh_backlist(GLB_GET,imolty,0,ibackbone)
c        --- safety check
         if ( ibackbone .eq. 0 ) then
            write(6,*) 'cannot perform a backbone move on a molecule '
            write(6,*) 'that does not have a protein backbone '
            write(6,*) 'molecule type ',imolty
            lfailure = .true.
            return
         endif
         rvalue = twh_random()*dble(ibackbone)
         ibackbone = idint(rvalue) + 1
         call twh_backlist(GLB_GET,imolty,ibackbone,iubacka)

c        --- count the number of torsion partners who are also backbone
c        --- atoms
         ibackbone = 0
         do itor = 1,twh_get_intor(imolty,iubacka)
            iumida = twh_get_ijtor2(imolty,iubacka,itor)
            iumidb = twh_get_ijtor3(imolty,iubacka,itor)
            iubackb = twh_get_ijtor4(imolty,iubacka,itor)
            call twh_lbackbone(GLB_GET,imolty,iumida,lbackbone)
            call twh_lbackbone(GLB_GET,imolty,iumidb,lbacktwo)
            call twh_lbackbone(GLB_GET,imolty,iubackb,lbackthr)
            if ( lbackbone .and. lbacktwo .and. lbackthr ) then
c              --- found another backbone torsion
               ibackbone = ibackbone + 1
               tempback(ibackbone) = itor
            endif
         enddo
c        --- better have at least one other backbone atom to choose from
         if ( ibackbone .eq. 0 ) then
            write(6,*) 'no other backbone atoms in schedule '
            write(6,*) 'imolty iubacka ',imolty,iubacka
            lfailure = .true.
            return
         endif

c        --- choose one of the backbone torsions from the list
         rvalue = twh_random()*dble(ibackbone)
         itor = tempback( idint(rvalue) + 1 )
c        --- assign the atoms
         iumida = twh_get_ijtor2(imolty,iubacka,itor)
         iumidb = twh_get_ijtor3(imolty,iubacka,itor)
         iubackb = twh_get_ijtor4(imolty,iubacka,itor)

c        --- make sure we have not just selected a backbone that will
c        --- regrow in such a way that a starting grow bead and the 
c        --- previous bead fill the N and CA slots of a proline
c        --- that case causes the algorithm to break down due to the
c        --- side chain cyling back on the backbone
         lbacka = .false.
         lmida = .false.
         lmidb = .false.
         lbackb = .false.
         call twh_aminoshort(GLB_GET,imolty,iubacka,aminoshort)
         if ( aminoshort .eq. 'PRO' ) then
            call twh_pdbname(GLB_GET,imolty,iubacka,pdbname)
            if ( pdbname .eq. ' N' .or. 
     &           pdbname .eq. ' CA' ) lbacka = .true.
         endif
         call twh_aminoshort(GLB_GET,imolty,iumida,aminoshort)
         if ( aminoshort .eq. 'PRO' ) then
            call twh_pdbname(GLB_GET,imolty,iumida,pdbname)
            if ( pdbname .eq. ' N' .or. 
     &           pdbname .eq. ' CA' ) lmida = .true.
         endif
         call twh_aminoshort(GLB_GET,imolty,iumidb,aminoshort)
         if ( aminoshort .eq. 'PRO' ) then
            call twh_pdbname(GLB_GET,imolty,iumidb,pdbname)
            if ( pdbname .eq. ' N' .or. 
     &           pdbname .eq. ' CA' ) lmidb = .true.
         endif
         call twh_aminoshort(GLB_GET,imolty,iubackb,aminoshort)
         if ( aminoshort .eq. 'PRO' ) then
            call twh_pdbname(GLB_GET,imolty,iubackb,pdbname)
            if ( pdbname .eq. ' N' .or. 
     &           pdbname .eq. ' CA' ) lbackb = .true.
         endif

         ldone = .false.
         if ( lbackb .and. .not. lmidb ) then
c           --- we will end up with a grow and prev bead that are 
c           --- both on the ring of proline
c           --- shift the atoms one space away from the b end
            ivib = 0
            do while ( ivib .lt. twh_get_invib(imolty,iubacka) 
     &           .and. .not. ldone)
               ivib = ivib + 1
               iu = twh_get_ijvib(imolty,iubacka,ivib)
               call twh_lbackbone(GLB_GET,imolty,iu,lbackbone)
               if ( iu .ne. iumida .and. lbackbone ) then
c                 --- this is the new bead a
                  iubackb = iumidb
                  iumidb = iumida
                  iumida = iubacka
                  iubacka = iu
                  ldone = .true.
               endif
            enddo
            if ( .not. ldone ) then 
c              --- could be a terminal proline, go back through and
c              --- take any atom that is not a carbon
               ivib = 0
               do while ( ivib .lt. twh_get_invib(imolty,iubacka) 
     &              .and. .not. ldone)
                  ivib = ivib + 1
                  iu = twh_get_ijvib(imolty,iubacka,ivib)
                  if ( iu .ne. iumida 
     &                 .and. 
     &                 twh_get_elementname(twh_get_ntype(imolty,iu))
     &                 .ne. ' C' ) then
c                    --- this is the new bead a
                     iubackb = iumidb
                     iumidb = iumida
                     iumida = iubacka
                     iubacka = iu
                     ldone = .true.
                  endif
               enddo
            endif
c           --- safety check
            if ( .not. ldone ) then
               write(6,*) 'SCHEDULE: problem with proline'
               lfailure = .true.
               return
            endif         
         elseif ( lbacka .and. .not. lmida ) then
c           --- we will end up with a grow and prev bead that are 
c           --- both on the ring of proline
c           --- shift the atoms one space away from the a end
            ldone = .false.
            ivib = 0
            do while ( ivib .lt. twh_get_invib(imolty,iubackb) 
     &           .and. .not. ldone)
               ivib = ivib + 1
               iu = twh_get_ijvib(imolty,iubackb,ivib)
               call twh_lbackbone(GLB_GET,imolty,iu,lbackbone)
               if ( iu .ne. iumidb .and. lbackbone ) then
c                 --- this is the new bead b
                  iubacka = iumida
                  iumida = iumidb
                  iumidb = iubackb
                  iubackb = iu
                  ldone = .true.
               endif
            enddo
            if ( .not. ldone ) then
c              --- could be a terminal proline, go back through and
c              --- take any atom that is not a carbon
               ivib = 0
               do while ( ivib .lt. twh_get_invib(imolty,iubackb) 
     &              .and. .not. ldone)
                  ivib = ivib + 1
                  iu = twh_get_ijvib(imolty,iubackb,ivib)
                  if ( iu .ne. iumidb 
     &                 .and. 
     &                 twh_get_elementname(twh_get_ntype(imolty,iu))
     &                 .ne. ' C' ) then
c                    --- this is the new bead b
                     iubacka = iumida
                     iumida = iumidb
                     iumidb = iubackb
                     iubackb = iu
                     ldone = .true.
                  endif
               enddo
            endif
c           --- safety check
            if ( .not. ldone ) then 
               write(6,*) 'SCHEDULE: problem with proline'
               lfailure = .true.
               return
            endif
         endif

         if ( ldone ) then
c            --- find the new torsion number for this
            itor = 0
            do while ( itor .lt. twh_get_intor(imolty,iubacka) 
     &           .and. ldone )
               itor = itor + 1
               if ( twh_get_ijtor2(imolty,iubacka,itor) .eq. iumida 
     &              .and. twh_get_ijtor3(imolty,iubacka,itor) 
     &              .eq. iumidb
     &              .and. twh_get_ijtor4(imolty,iubacka,itor) 
     &              .eq. iubackb
     &              ) ldone = .false.
            enddo
         endif

c        --- set up info for iubacka
         index = 1
         call twh_growfrom(GLB_SET,index,iubacka)
         lstart(iubacka) = .true.
c        --- see how many neighboring atoms are both on the backbone
c        --- and not already iumida
         ibackbone = 0
         do ivib = 1,twh_get_invib(imolty,iubacka)
            iut = twh_get_ijvib(imolty,iubacka,ivib)
            call twh_lbackbone(GLB_GET,imolty,iut,lbackbone)
            if ( lbackbone .and. iut .ne. iumida ) then
               ibackbone = ibackbone + 1
               tempback(ibackbone) = ivib
            endif
         enddo
         if ( ibackbone .eq. 0 ) then
c           --- must be one of the termini, choose one of the atoms 
c           --- (other than iumida) to remain in place
            temp_count = 0
            do ivib = 1,twh_get_invib(imolty,iubacka)
               iut = twh_get_ijvib(imolty,iubacka,ivib)
               if ( iut .ne. iumida ) then
                  temp_count = temp_count + 1
                  temp_store(temp_count) = iut
               endif
            enddo
            if ( temp_count .eq. 0 ) then
c              --- there are no other atoms bonded to this atom
               call twh_growprev(GLB_SET,index,0)
            else
c              --- select one of the neighboring atoms at random
               rvalue = twh_random()*dble(temp_count)
               growprev = temp_store(idint(rvalue)+1)
               call twh_growprev(GLB_SET,index,growprev)
            endif
         else
c           --- choose one of the backbone atoms as prev
            rvalue = twh_random()*dble(ibackbone)
            growprev = twh_get_ijvib(imolty,iubacka
     &           ,tempback( idint(rvalue) + 1))
            call twh_growprev(GLB_SET,index,growprev)
         endif
c        --- set up info for iubackb
         index = 2
         call twh_growfrom(GLB_SET,index,iubackb)
         lstart(iubackb) = .true.
c        --- see how many neighboring atoms are both on the backbone
c        --- and not already iumida
         ibackbone = 0
         do ivib = 1,twh_get_invib(imolty,iubackb)
            iut = twh_get_ijvib(imolty,iubackb,ivib)
            call twh_lbackbone(GLB_GET,imolty,iut,lbackbone)
            if ( lbackbone .and. iut .ne. iumidb ) then
               ibackbone = ibackbone + 1
               tempback(ibackbone) = ivib
            endif
         enddo
         if ( ibackbone .eq. 0 ) then
c           --- must be one of the termini, choose one of the atoms 
c           --- (other than iumidb) to remain in place
            temp_count = 0
            do ivib = 1,twh_get_invib(imolty,iubackb)
               iut = twh_get_ijvib(imolty,iubackb,ivib)
               if ( iut .ne. iumidb ) then
                  temp_count = temp_count + 1
                  temp_store(temp_count) = iut
               endif
            enddo
            if ( temp_count .eq. 0 ) then
c              --- there are no other atoms bonded to this atom
               call twh_growprev(GLB_SET,index,0)
            else
c              --- select one of the neighboring atoms at random
               rvalue = twh_random()*dble(temp_count)
               growprev = temp_store(idint(rvalue)+1)
               call twh_growprev(GLB_SET,index,growprev)
            endif
         else
c           --- choose one of the atoms as prev
            rvalue = twh_random()*dble(ibackbone)
            growprev = twh_get_ijvib(imolty,iubackb
     &           ,tempback( idint(rvalue) + 1))
            call twh_growprev(GLB_SET,index,growprev)
         endif

c        --- grow all (except iprev) of the beads bonded to start beads
         do ind_start = 1,2
            if ( ind_start .eq. 1 ) then
               iutry = iubacka
            else
               iutry = iubackb
            endif
            temp_count = 0
            do ivib = 1,twh_get_invib(imolty,iutry)
               iut = twh_get_ijvib(imolty,iutry,ivib)
               call twh_growprev(GLB_GET,ind_start,growprev)
               if ( iut .ne. growprev ) then
                  temp_count = temp_count + 1
                  temp_store(temp_count) = iut
               endif
            enddo

            count = 0
            do zz = temp_count,1,-1
c              --- choose grow bead randomly from temp_store
               rvalue = dble(zz)* twh_random()
               itry = idint(rvalue) + 1
               iut = temp_store(itry)
               count = count + 1
               call twh_growlist(GLB_SET,ind_start,count,iut)
c              --- update temp_store for next iteration
               temp_store(itry) = temp_store(zz)
            enddo
            call twh_grownum(GLB_SET,ind_start,count)
         enddo
c        --- bit of trickery here to store the torsion number for 
c        --- use in cbregrow
         iutry = iubacka
         moveinfo = itor
c        --- end part that is specific for backbone cbregrow
      elseif ( movetype .eq. 5 ) then
c        --- regrowth of all atoms, except for certain specified ones
c        --- that remain in the same position.  Originally created to 
c        --- construct hexlical DNA from just the P positions
         call twh_cbmc_atomname_tot(GLB_GET,imolty,index)
         do ind_start = 1,index
            call twh_cbmc_atomname_list(GLB_GET,imolty,ind_start,iu)
            call twh_growfrom(GLB_SET,ind_start,iu)
c           --- set the grownum to the number of vibrations
            grownuma = twh_get_invib(imolty,iu)
            call twh_grownum(GLB_SET,ind_start,grownuma)
c           --- set growlist (don't bother with random for now)
            do ivib = 1,twh_get_invib(imolty,iu)
               iut = twh_get_ijvib(imolty,iu,ivib)
               call twh_growlist(GLB_SET,ind_start,ivib,iut)
            enddo
c           --- set lstart
            lstart(iu) = .true.
c           --- set growprev to zero
            call twh_growprev(GLB_SET,ind_start,0)
         enddo
      elseif ( movetype .eq. 6 ) then
c        --- regrowth of a single sidechain on a polypeptide.  Can 
c        --- regrow more than one sidechain if they are bonded together
c        --- the most common instance of this is disulfide bonds

c        --- reinitialize lstart to be true for all atoms
         do iunit = 1,nunit
            lstart(iunit) = .true.
         enddo

c        --- compile a list of 'CA' main chain atoms
         nca = 0
         do iunit = 1,nunit
            call twh_pdbname(GLB_GET,imolty,iunit,pdbname)
            if ( pdbname .eq. ' CA ' ) then
               nca = nca + 1
               blacklist(nca) = iunit
            endif
         enddo

         if ( nca .eq. 0 ) then
            write(6,*) 'SCHEDULE: cannot perform a cbmc sidechain move'
            write(6,*) 'on a molecule that does not have any CA atoms'
            write(6,*) 'molecule type:',imolty
            lfailure = .true.
            return
         endif

c        --- select a CA atom at random from the list
         rvalue = twh_random()*dble(nca)
         nca = idint(rvalue) + 1
         target = blacklist(nca)
         moveinfo = target
c        --- compile a list of all atoms that start from the selected
c        --- CA and can be grown unit hitting the main chain again.
c        --- we reuse the blacklist variable, but it now has a different
c        --- meaning from above
         black_index = 0
         test_index = 0
         ldone = .false.
         do while ( .not. ldone )
            do ivib = 1,twh_get_invib(imolty,target)
               junit = twh_get_ijvib(imolty,target,ivib)
               lvalid = .true.
               call twh_pdbname(GLB_GET,imolty,junit,pdbname)
               if ( pdbname .eq. ' N  ' .or.
     &              pdbname .eq. ' CA ' .or.
     &              pdbname .eq. ' C  ' .or.
     &              pdbname .eq. ' O  ' ) then
c                 --- main chain atom
                  lvalid = .false.
               else
c                 --- make sure it is not already included in the list
                  if ( .not. lstart(junit) ) lvalid = .false.
               endif
               if ( lvalid ) then
c                 --- valid atom found for further growth
                  black_index = black_index + 1
                  blacklist(black_index) = junit
c                 --- set lstart
                  lstart(junit) = .false.
               endif
            enddo
            test_index = test_index + 1
            if ( test_index .gt. black_index ) then
               ldone = .true.
            else
               target = blacklist(test_index)
            endif
         enddo

c        --- now we have the list of all the atoms on the blacklist so
c        --- go through the entire molecule and select all the atoms
c        --- that are not on the blacklist as growth initiation sites
         index = nunit - black_index
         iu = 0
         do ind_start = 1,index
c           --- skip the units that are on the blacklist
            ldone = .false.
            do while ( .not. ldone )
               iu = iu + 1
               lvalid = .true.
               do test_index = 1,black_index
                  if ( iu .eq. blacklist(test_index) ) then
                     lvalid = .false.
                  endif
               enddo
               if ( lvalid ) then
c                 --- we are done
                  ldone = .true.
               endif
               if ( iu .gt. nunit ) then
                  write(6,*) 'SCHEDULE: logic failure in sidechain'
                  lfailure = .true.
                  return
               endif
            enddo
            call twh_growfrom(GLB_SET,ind_start,iu)
c           --- set the grownum to the number of vibrations
            grownuma = twh_get_invib(imolty,iu)
            call twh_grownum(GLB_SET,ind_start,grownuma)
c           --- set growlist in random order
            nvibs = twh_get_invib(imolty,iu)
            do ivib = 1,nvibs
               temp_store(ivib) = twh_get_ijvib(imolty,iu,ivib)
            enddo
c           --- randomize the order in the growlist
            count = 0
            do while ( nvibs .ne. 0 )
               count = count + 1
               rvalue = twh_random()*dble(nvibs)
               ivib = idint(rvalue)+1
               growtest = temp_store(ivib)
               call twh_growlist(GLB_SET,ind_start,count,growtest)
               temp_store(ivib) = temp_store(nvibs)
               nvibs = nvibs - 1
            enddo
c           --- set lstart
            lstart(iu) = .true.
c           --- set growprev to zero
            call twh_growprev(GLB_SET,ind_start,0)
         enddo
      else
c        --- non-existant move type
         write(6,*) 'schedule movetype ',movetype
         write(6,*) 'non-valid move type'
         lfailure = .true.
         return
      endif

c     --- from here on down all movetypes are the same
c     --- set up list of "outer" beads that have further growth sites

c     --- set up the inital logic and then feed this into the 
c     --- recursive loop until all of the logic is complete
c     --- note that some of the logic for the initial setup is 
c     --- different from the logic that occurs in the main loop 
c     --- this is not a mistake, but is required to handle cyclics
      ibead = 0
      outer_num = 0

      do ind_start = 1,index
         call twh_growfrom(GLB_GET,ind_start,iufrom)
c        --- cycle through the atoms that we could potentially grow
c        --- from this atom and make sure that they are not other
c        --- starter atoms
         icbu = 1
         call twh_grownum(GLB_GET,ind_start,grownuma)
         do while ( icbu .le. grownuma )
            call twh_growlist(GLB_GET,ind_start,icbu,iu)
            lfound = .false.
            do isearch = 1,index
               if ( isearch .ne. ind_start ) then
                  call twh_growfrom(GLB_GET,isearch,growfrom)
                  if ( iu .eq. growfrom ) then
c                    --- this is actually another starter atom
                     lfound = .true.
                     ifound = isearch
                  endif
               endif
            enddo
            if ( lfound ) then
               call twh_growprev(GLB_GET,ind_start,growprev)
               if ( growprev .eq. 0 ) then
c                 --- set growprev to this value.
                  call twh_growlist(GLB_GET,ind_start,icbu,growprev)
                  call twh_growprev(GLB_SET,ind_start,growprev)
               else
c                 --- we already have a growprev so this goes in the
c                 --- growextra array instead
c                 --- increment the growxtra array
                  call twh_growxtra(GLB_INCR,ind_start,1,1)
c                 --- add the unit number into the growxtra array
                  call twh_growxtra(GLB_GET,ind_start,1,itemp)
                  itemp = itemp + 1
                  call twh_growlist(GLB_GET,ind_start,icbu,growtest)
                  call twh_growxtra(GLB_SET,ind_start,itemp,growtest)
               endif
c              --- store the current total number of neighbors to be 
c              --- grown
               call twh_grownum(GLB_GET,ind_start,icurrent)
c              --- remove the atom found to already exist from the 
c              --- grownum and replace with the final one
               call twh_growlist(GLB_GET,ind_start,icurrent,growtest)
               call twh_growlist(GLB_SET,ind_start,icbu,growtest)
               call twh_grownum(GLB_DECR,ind_start,1)
c              --- fetch the new value of grownum(ind_start)
               call twh_grownum(GLB_GET,ind_start,grownuma)
            else
c              --- valid growth site, increment icbu
               icbu = icbu + 1
            endif
         enddo
      enddo

c     --- loop through again and remove any of the original start atoms
c     --- that now has 0 grownum
      ind_start = 1
      do while ( ind_start .le. index )
         call twh_grownum(GLB_GET,ind_start,grownuma)
         if ( grownuma .eq. 0 ) then
c           --- remove and replace with index
            call twh_growfrom(GLB_GET,index,growfrom)
            call twh_growfrom(GLB_SET,ind_start,growfrom)
            call twh_growprev(GLB_GET,index,growprev)
            call twh_growprev(GLB_SET,ind_start,growprev)
            call twh_grownum(GLB_GET,index,grownumb)
            call twh_grownum(GLB_SET,ind_start,grownumb)
            do icbu = 1,grownumb
               call twh_growlist(GLB_GET,index,icbu,growtest)
               call twh_growlist(GLB_SET,ind_start,icbu,growtest)
            enddo
            call twh_growxtra(GLB_GET,index,1,xtracount)
            call twh_growxtra(GLB_SET,ind_start,1,xtracount)
            do ixtra = 2,1+xtracount
               call twh_growxtra(GLB_GET,index,ixtra,itemp)
               call twh_growxtra(GLB_SET,ind_start,ixtra,itemp)
            enddo
c           --- decrement index
            index = index - 1
         else
c           --- nothing to do here, move on to the next atom in the list
            ind_start = ind_start + 1
         endif
      enddo

c     --- set lexsched to false for all of the atoms that are still
c     --- listed in the grownum arrays for the current index
      do ind_start = 1,index
         call twh_grownum(GLB_GET,ind_start,grownuma)
         do icbu = 1,grownuma
            call twh_growlist(GLB_GET,ind_start,icbu,iu)
            call twh_lexsched(GLB_SET,iu,.false.)
         enddo
      enddo

      do ind_start = 1,index
         call twh_growfrom(GLB_GET,ind_start,iufrom)
c        --- cycle through the atoms that we are growing from iufrom 
c        --- and look for those which will require a further growth step
         call twh_grownum(GLB_GET,ind_start,grownuma)
         do icbu = 1,grownuma
c           --- increment counter ibead for total number of beads grown
            ibead = ibead + 1
c           --- determine whether this bead has any non-grown neighbors
            call twh_growlist(GLB_GET,ind_start,icbu,iu)
            lfound = .false.
            do ivib = 1,twh_get_invib(imolty,iu)
               iut = twh_get_ijvib(imolty,iu,ivib)
               if ( iut .ne. iufrom ) then
                  lhook = .false.
                  if ( lstart(iut) ) then
                     lhook = .true.
                  else
c                    --- make sure we didn't just grow these atoms in
c                    --- a previous index value at the start.  Only 
c                    --- happens when we are doing a fixed-endpoint
                     do ind_other = 1,ind_start-1
                        call twh_grownum(GLB_GET,ind_other,grownumb)
                        do iother = 1,grownumb
                           call twh_growlist(GLB_GET,ind_other,iother
     &                          ,growtest)
                           if ( growtest .eq. iut ) lhook = .true.
                        enddo
                     enddo
                  endif
                  if ( lhook ) then
c                    --- must be a ring, or we are hooking up to one of
c                    --- the starting sites
                     call twh_grownbon(GLB_INCR,ind_start,icbu,1)
                     call twh_grownbon(GLB_GET,ind_start,icbu,grownbon)
                     call twh_growbond(GLB_SET,ind_start,icbu,grownbon
     &                    ,iut)
                  else
c                    --- this atom has more atoms attached to it that 
c                    --- need to be grown
                     lfound = .true.
                  endif
               endif
            enddo
         
            if ( lfound ) then
c              --- add one to the stack of outer_sites
c              --- note that this just makes this a potential outer 
c              --- growth site, which is definately a site if this is a 
c              --- non-cyclic molecule with a regular regrowth, but may
c              --- no longer be a relevent site when it is time to 
c              --- assign it below in case of a cyclic molecule, or for
c              --- a fixed-endpoint regrowth
               outer_num = outer_num + 1
               outer_sites(outer_num) = iu
               outer_prev(outer_num) = iufrom
            endif
         enddo
      enddo

c     --- begin while loop to grow all outer beads until done
      do while ( outer_num .gt. 0 ) 
c        --- choose one site randomly from the stack
         rvalue =  dble(outer_num)* twh_random()
         outer_try = int( rvalue ) + 1
c        --- increment index (note will decrement index later if this
c        --- is not a valid growth site)
         index = index + 1
c        --- initialize growxtra
         call twh_growxtra(GLB_SET,index,1,0)
c        --- assign growfrom and growprev for this index
         call twh_growfrom(GLB_SET,index,outer_sites(outer_try))
         growprev = outer_prev(outer_try)
         call twh_growprev(GLB_SET,index,growprev)
c        --- temporary assignments of from 
         call twh_growfrom(GLB_GET,index,iufrom)

c        --- determine whether this is still a valid growth site
c        --- a site may become invalid if it is a cyclic molecule
c        --- if it is valid then start filling up the data structures
         temp_count = 0
         do ivib = 1, twh_get_invib(imolty,iufrom)
            iut = twh_get_ijvib(imolty,iufrom,ivib)
            call twh_growprev(GLB_GET,index,growprev)
            if ( iut .ne. growprev ) then
               call twh_lexsched(GLB_GET,iut,lexschedule)
               if ( lexschedule .and. .not. lstart(iut ) ) then
c                 --- add to the list of beads to be grown from iu
                  temp_count = temp_count + 1
                  temp_store(temp_count) = iut
               else
c                 --- we are growing from an atom which is already 
c                 --- bonded to other atoms that we are not regrowing
c                 --- other than the normal previous atom
c                 --- this happens when using fecbmc
                  call twh_growxtra(GLB_INCR,index,1,1)
                  call twh_growxtra(GLB_GET,index,1,itemp)
                  itemp = itemp+1
                  call twh_growxtra(GLB_SET,index,itemp,iut)
               endif
            endif
         enddo

c        --- assign the grow beads in random order
         count = 0
         do zz = temp_count, 1, -1
            rvalue = dble(zz)* twh_random()
            itry = int( rvalue ) + 1
            iut = temp_store(itry)
            count = count + 1
c           --- assign growlist for current index and count
            call twh_growlist(GLB_SET,index,count,iut)
c           --- update logic
            call twh_lexsched(GLB_SET,iut,.false.)
c           --- update temp_store for next iteration
            temp_store(itry) = temp_store(zz)
         enddo
c        --- assign grownum for this index
         call twh_grownum(GLB_SET,index,count)

c        --- determine whether the atoms grown in this step are 
c        --- bonded to any other atoms.  If so then they either need
c        --- to be added to the list of target sites, or there is 
c        --- some special bonded interaction we need to account for 
c        --- in this growth step
         call twh_grownum(GLB_GET,index,grownuma)
         do icbu = 1,grownuma
c           --- increment counter ibead for total number of beads grown
            ibead = ibead + 1
c           --- determine whether this bead has any non-grown neighbors
            call twh_growlist(GLB_GET,index,icbu,iu)

            lfound = .false.
            do ivib = 1, twh_get_invib(imolty,iu)
               iut = twh_get_ijvib(imolty,iu,ivib)
c              --- check against the current from bead (instead of prev)
               if ( iut .ne. iufrom ) then
                  call twh_lexsched(GLB_GET,iut,lexschedule)
                  if ( .not. lexschedule .or. lstart(iut) ) then
c                    --- must be a ring, or we are hooking up to 
c                    --- one of the starting end sites
                     call twh_grownbon(GLB_INCR,index,icbu,1)
                     call twh_grownbon(GLB_GET,index,icbu,grownbon)
                     call twh_growbond(GLB_SET,index,icbu,grownbon,iut)
                  elseif ( .not. lfound ) then
c                    --- this atom has more atoms attached to it that
c                    --- need to be grown, and it has not yet been 
c                    --- accounted for
c                    --- add one to the stack of outer_sites
                     outer_num = outer_num + 1
                     outer_sites(outer_num) = iu
c                    --- the current from is the next prev
                     call twh_growfrom(GLB_GET,index
     &                    ,outer_prev(outer_num))
c                    --- make sure we don't add this atom to the list
c                    --- again
                     lfound = .true.
                  endif
               endif
            enddo
         enddo

c        --- update list of "outer" beads
c        --- remove bead that was just a growth target from outer list
c        --- note that this bead is removed whether it actually was a
c        --- growth target or not
         outer_sites(outer_try) = outer_sites(outer_num)
         outer_prev(outer_try) = outer_prev(outer_num)
         outer_num = outer_num - 1
c        --- decrement index if this was not actually a growth target
         if ( count .eq. 0 ) index = index - 1
      enddo

c     --- set moveinfo to ibead if this is movetype 0
      if ( movetype .eq. 0 ) moveinfo = ibead

#if DEBUG_SCHEDULE
      write(6,*) 'movetype',movetype
      write(6,*) 'index',index
      write(6,*) 'numunits',numunits
      do ibead = 1,index
         write(6,*) 'ibead',ibead
         call twh_growfrom(GLB_GET,ibead,growfrom)
         write(6,*) 'growfrom(ibead)',growfrom
         write(6,*) 'lstart(growfrom(ibead))',lstart(growfrom)
         call twh_growprev(GLB_GET,ibead,growprev)
         write(6,*) 'growprev(ibead)',growprev
         call twh_grownum(GLB_GET,ibead,grownuma)
         write(6,*) 'grownum(ibead)',grownuma
         call twh_growxtra(GLB_GET,ibead,1,xtracount)
         do itemp = 1,xtracount
            call twh_growxtra(GLB_GET,ibead,itemp+1,temp_store(itemp))
         enddo
         write(6,*) 'growxtra num:',xtracount,' vals:'
     &        ,(temp_store(itemp),itemp=1,xtracount)
         do count = 1,grownuma
            call twh_growlist(GLB_GET,ibead,count,growtest)
            write(6,*) 'count,growlist(ibead,count)',count,growtest
            call twh_grownbon(GLB_GET,ibead,count,grownbon)
            do itemp = 1,grownbon
               call twh_growbond(GLB_GET,ibead,count,itemp
     &              ,temp_store(itemp))
            enddo
            write(6,*) 'grownbon,growbond',grownbon
     &           ,(temp_store(itemp),itemp=1,grownbon)
         enddo
      enddo
#endif
      
c      write(6,*) 'end SCHEDULE'
      return
      end
