#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2007 Marcus G. Martin                       *
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
      subroutine twh_assemble(lfailure,imolty,totatoms,spectorcount
     &     ,forcefield,atomname,match_style,improper_assignment
     &     ,lautotorsion)
c     ******************************************************************
c     * assembles the input information suitable for the code          *
c     * imolty: the molecule type                                      *
c     * totatoms: the total number of atoms in the molecule            *
c     * forcefield: text names of the forcefield we are using          *
c     * atomname: text names of each atom in the molecule using the    *
c     *           appropriate names for the forcefield                 *
c     * match_style: way of finding matches                            *
c     *              'standard' exact matches followed by wildcards    *
c     *              'MMFF checkdown' exact matches for nonbond and    *
c     *                 bond and a two character checkdown for angle   *
c     *                 and torsions                                   *
c     *              'MMFF bond checkdown' same as MMFF checkdown      *
c     *                 except the nonbond is matched to the bond      *
c     *                 values instead of the normal nonbond values    *
c     *              'stereoselective' similar to standard, but with   *
c     *                 some additional logic to distinguish between   *
c     *                 torsions made up of the same atoms, but with   *
c     *                 different stereochemistries (ie cis vs. trans) *
c     * lautotorsion: true if we are autodetecting all torsions, false *
c     *               if some torsions have already been determined    *
c     *                                                                *
c     * originally written 02-26-2001 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_AAFF
#define FUNCTION_GET_AANAMES
#define FUNCTION_GET_CANAME
#define FUNCTION_GET_CBNAME
#define FUNCTION_GET_CTNAME
#define FUNCTION_GET_IJBEN2
#define FUNCTION_GET_IJBEN3
#define FUNCTION_GET_IJTOR2
#define FUNCTION_GET_IJTOR3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_IMPFF
#define FUNCTION_GET_INBEN
#define FUNCTION_GET_INIMPROP
#define FUNCTION_GET_INOF
#define FUNCTION_GET_INTOR
#define FUNCTION_GET_ITIMPROP
#define FUNCTION_GET_ITTOR
#define FUNCTION_GET_NAASAME
#define FUNCTION_GET_NBFF
#define FUNCTION_GET_NBNAME
#define FUNCTION_GET_NTYPE
#define FUNCTION_GET_OFFF
#define FUNCTION_GET_OFNAMES
#include "functions.h"
c     --- debug variable
#define DEBUG_ASSEMBLE 0
c     --- variables passed to/from the subroutine
      logical lfailure,lautotorsion
      character*(FFNAMELEN) forcefield
      character*(FFNAMELEN) atomname
      dimension atomname(NUMAX)
      character*(*) match_style,improper_assignment
      integer totatoms,imolty,spectorcount
c     --- local variables
c     --- logical scalars
      logical lfound,langang
c     --- character string scalars
      character*(FFNAMELEN) iname,jname,kname,lname
      character*15 specialtor
c     --- character string arrays
      character*(FFNAMELEN) names
      dimension names(5)
      character*15 spectornames
      dimension spectornames(spectorcount)
c     --- integer scalars
      integer param,type,isame,iof,iorder
      integer iatom,jatom,katom,latom,matom,ivib,ibend,itor,jvib
      integer iaa,iindex,kindex,lindex
      integer iatom_nvib,jatom_nvib,katom_nvib,latom_nvib,itest,jtest
c     --- integer arrays
      integer partners
      dimension partners(3)
      integer specialtorsion
      dimension specialtorsion(spectorcount,4)

c      write(6,*) 'start ASSEMBLE'

c     --- determine atom types
      do iatom = 1,totatoms
         names(1) = atomname(iatom)
#if DEBUG_ASSEMBLE
      write(6,*) 'iatom,atomname:',iatom,names(1)
#endif
         param = 0
         type = 0
         do while ( param .eq. 0 .and. type .lt. NNTYPE)
            type = type + 1
            if ( forcefield .eq. twh_get_nbff(type) ) then
               if ( match_style .eq. 'standard'
     &              .or. match_style .eq. 'MMFF checkdown' 
     &              .or. match_style .eq. 'stereoselective' ) then
c                 --- match to the nonbonded name
                  if ( twh_get_nbname(type) .eq. names(1) ) then
                     param = type
                  endif
               elseif ( match_style .eq. 'MMFF bond checkdown' ) then
c                 --- match to the bond name
                  if ( twh_get_cbname(type) .eq. names(1) ) then
                     param = type
                  endif
               else
                  write(6,*) 'ASSEMBLE: unknown match_style for nonbond'
                  write(6,*) match_style
                  lfailure = .true.
                  return
               endif
            endif
         enddo
         if ( param .eq. 0 ) then
            write(6,*) 'parameters not found for atom'
            write(6,*) 'iatom',iatom,' name'
            write(6,*) names(1)
            write(6,*) 'forcefield'
            write(6,*) forcefield
            lfailure = .true.
            return
         endif
         call twh_ntype(GLB_SET,imolty,iatom,param)
      enddo

c     --- determine vibration types
c     --- currently independent of match_style
      do iatom = 1,totatoms
         call twh_invib(GLB_GET,imolty,iatom,iatom_nvib)
         do ivib = 1,iatom_nvib
            call twh_ijvib(GLB_GET,imolty,iatom,ivib,jatom)
            call twh_find_bond_type(lfailure,totatoms,imolty,iatom
     &           ,jatom,ivib,param,forcefield)
            if ( lfailure ) return
            call twh_itvib(GLB_SET,imolty,iatom,ivib,param)
         enddo
      enddo

c     --- determine bending angles and types
      do iatom = 1,totatoms
         call twh_invib(GLB_GET,imolty,iatom,iatom_nvib)
         do ivib = 1,iatom_nvib
            call twh_ijvib(GLB_GET,imolty,iatom,ivib,jatom)
            call twh_invib(GLB_GET,imolty,jatom,jatom_nvib)
            do jvib = 1,jatom_nvib
               call twh_ijvib(GLB_GET,imolty,jatom,jvib,katom)
               if ( iatom .ne. katom ) then
                  iorder = ORD_NULL
c                 --- initialize param (only done once)
                  param = 0
c                 --- we found a bending angle
                  call twh_inben(GLB_INCR,imolty,iatom,1)
                  ibend = twh_get_inben(imolty,iatom)
                  call twh_ijben2(GLB_SET,imolty,iatom,ibend,jatom)
                  call twh_ijben3(GLB_SET,imolty,iatom,ibend,katom)
                  if ( match_style .eq. 'standard' 
     &                 .or. match_style .eq. 'stereoselective' ) then
                     call twh_angle_match_standard(lfailure,imolty,iatom
     &                    ,jatom,katom,iorder,param,forcefield)
                     if ( lfailure ) return
                  elseif ( match_style .eq. 'MMFF checkdown' ) then
                     call twh_angle_match_mmff(lfailure,imolty,iatom
     &                    ,jatom,katom,iorder,param,forcefield)
                     if ( lfailure ) return
                  else
                     write(6,*) 'ASSEMBLE: unknown match_style'
                     write(6,*) match_style
                     lfailure = .true.
                     return
                  endif

                  call twh_orderben(GLB_SET,imolty,iatom,ibend,iorder)
                  call twh_itben(GLB_SET,imolty,iatom,ibend,param)
               endif
            enddo
         enddo
      enddo

c     --- read in additional information about the special torsions
c     --- if needed
      if ( spectorcount .gt. 0 ) then
         do itest = 1,spectorcount
            read(4,*) (specialtorsion(itest,itor),itor=1,4)
     &           ,spectornames(itest)
            write(6,*) (specialtorsion(itest,itor),itor=1,4)
     &           ,spectornames(itest)
         enddo
      endif

c     --- determine regular torsion angles
      do iatom = 1,totatoms
         do ibend = 1,twh_get_inben(imolty,iatom)
            jatom = twh_get_ijben2(imolty,iatom,ibend)
            katom = twh_get_ijben3(imolty,iatom,ibend)
            call twh_invib(GLB_GET,imolty,katom,katom_nvib)
            do ivib = 1,katom_nvib
               call twh_ijvib(GLB_GET,imolty,katom,ivib,latom)
               if ( latom .ne. jatom .and. iatom .ne. latom ) then
c                 --- we found a valid regular torsion
                  if ( lautotorsion ) then
c                    --- add this to the list of automatically 
c                    --- detected torsions
                     call twh_intor(GLB_INCR,imolty,iatom,1)
                     itor = twh_get_intor(imolty,iatom)
                     call twh_ijtor2(GLB_SET,imolty,iatom,itor,jatom)
                     call twh_ijtor3(GLB_SET,imolty,iatom,itor,katom)
                     call twh_ijtor4(GLB_SET,imolty,iatom,itor,latom)
                  else
c                    --- see if this torsion already is specified
                     lfound = .false.
                     do itor = 1,twh_get_intor(imolty,iatom)
                        if ( twh_get_ijtor2(imolty,iatom,itor) 
     &                       .eq. jatom
     &                       .and. twh_get_ijtor3(imolty,iatom,itor) 
     &                       .eq. katom
     &                       .and. twh_get_ijtor4(imolty,iatom,itor) 
     &                       .eq. latom ) lfound = .true.
                     enddo
                     if ( .not. lfound ) then
c                       --- add to the end of the list
                        call twh_intor(GLB_INCR,imolty,iatom,1)
                        itor = twh_get_intor(imolty,iatom)
                        call twh_ijtor2(GLB_SET,imolty,iatom,itor,jatom)
                        call twh_ijtor3(GLB_SET,imolty,iatom,itor,katom)
                        call twh_ijtor4(GLB_SET,imolty,iatom,itor,latom)
c                       --- set special type flag
                        call twh_ittor(GLB_SET,imolty,iatom,itor,-20)
                     endif
                  endif
               endif
            enddo
         enddo
      enddo

      do iatom = 1,totatoms
         do itor = 1,twh_get_intor(imolty,iatom)
            jatom = twh_get_ijtor2(imolty,iatom,itor)
            katom = twh_get_ijtor3(imolty,iatom,itor)
            latom = twh_get_ijtor4(imolty,iatom,itor)
c           --- set param to the torsion type, generally this is 0
c           --- by initialization, but allows for clever trickery
            param = twh_get_ittor(imolty,iatom,itor)

            if ( match_style .eq. 'standard' ) then
               specialtor = 'wild'
               call twh_torsion_match_standard(lfailure,imolty,iatom
     &              ,jatom,katom,latom,param,forcefield,specialtor)
               if ( lfailure ) return
            elseif ( match_style .eq. 'MMFF checkdown' ) then
               call twh_torsion_match_mmff(lfailure,imolty,iatom
     &              ,jatom,katom,latom,param,forcefield)
               if ( lfailure ) return
            elseif ( match_style .eq. 'stereoselective' ) then
c              --- check for a match in the special torsion list
               specialtor = 'wild'
               itest = 0
               do while ( itest .lt. spectorcount .and.
     &              specialtor .eq. 'wild' )
                  itest = itest + 1
                  if (
     &                 ( iatom .eq. specialtorsion(itest,1) .and.
     &                 jatom .eq. specialtorsion(itest,2) .and.
     &                 katom .eq. specialtorsion(itest,3) .and.
     &                 latom .eq. specialtorsion(itest,4) )
     &                 .or.
     &                 ( iatom .eq. specialtorsion(itest,4) .and.
     &                 jatom .eq. specialtorsion(itest,3) .and.
     &                 katom .eq. specialtorsion(itest,2) .and.
     &                 latom .eq. specialtorsion(itest,1) )
     &                 ) then
                     specialtor = spectornames(itest)
                  endif
               enddo
               call twh_torsion_match_standard(lfailure,imolty,iatom
     &              ,jatom,katom,latom,param,forcefield,specialtor)
               if ( lfailure ) return
            else
               write(6,*) 'ASSEMBLE: invalid match_style for torsions'
               write(6,*) match_style
               lfailure = .true.
               return
            endif
            call twh_ittor(GLB_SET,imolty,iatom,itor,param)
         enddo
      enddo
c     --- find out if this forcefield contains angle-angle terms
      langang = .false.
      do iaa = 1,TAAMAX
         if ( forcefield .eq. twh_get_aaff(iaa) ) langang = .true.
      enddo
      
c     --- determine the angle-angle terms and their types
      do jatom = 1,totatoms
         call twh_invib(GLB_GET,imolty,jatom,jatom_nvib)
         if ( jatom_nvib .ge. 4 .and. langang ) then
            iaa = 0
            do kindex = 1,jatom_nvib
               call twh_ijvib(GLB_GET,imolty,jatom,kindex,katom)
               do iindex = 1,jatom_nvib
                  call twh_ijvib(GLB_GET,imolty,jatom,iindex,iatom)
                  do lindex = iindex+1,jatom_nvib
                     call twh_ijvib(GLB_GET,imolty,jatom,lindex,latom)
                     if ( iindex .ne. kindex .and. lindex .ne. kindex)
     &                    then
c                       --- found a valid angle-angle term
                        iaa = iaa + 1
                        call twh_ijaa0(GLB_SET,imolty,jatom,iaa,katom)
                        call twh_ijaa1(GLB_SET,imolty,jatom,iaa,iatom)
                        call twh_ijaa2(GLB_SET,imolty,jatom,iaa,latom)
c                       --- send the compass angle names for this type
c                       --- use the i,j,k,l order of compass
                        names(1) = twh_get_caname(
     &                       twh_get_ntype(imolty,iatom))
                        names(2) = twh_get_caname(
     &                       twh_get_ntype(imolty,jatom))
                        names(3) = twh_get_caname(
     &                       twh_get_ntype(imolty,katom))
                        names(4) = twh_get_caname(
     &                       twh_get_ntype(imolty,latom))
c                       --- match the parameters
c                       --- positions 1 and 4 are interchangable, 
c                       --- but 2 and 3 must be in this exact order
                        param = 0
                        type = 0
                        do while ( param .eq. 0 .and. type .lt. TAAMAX )
                           type = type + 1
                           if ( forcefield .eq. twh_get_aaff(type)
     &                          ) then
                              do isame = 1,twh_get_naasame(type)
                                 if ( ( names(1)
     &                                .eq. twh_get_aanames(type,isame,1)
     &                                .and. names(2) 
     &                                .eq. twh_get_aanames(type,isame,2)
     &                                .and. names(3) 
     &                                .eq. twh_get_aanames(type,isame,3)
     &                                .and. names(4) 
     &                                .eq. twh_get_aanames(type,isame,4)
     &                                ) .or. (
     &                                names(4) 
     &                                .eq. twh_get_aanames(type,isame,1)
     &                                .and. names(2) 
     &                                .eq. twh_get_aanames(type,isame,2)
     &                                .and. names(3) 
     &                                .eq. twh_get_aanames(type,isame,3)
     &                                .and. names(1) 
     &                                .eq. twh_get_aanames(type,isame,4)
     &                                ) ) then
                                    param = type
                                 endif
                              enddo
                           endif
                        enddo

                        if ( param .eq. 0  ) then
                           write(6,*) 'parameters not found for ang-ang'
                           write(6,*) 'atoms ',iatom,jatom,katom,latom
                           write(6,*) 'names:'
                           write(6,'(4(a10,/))') names(1),names(2)
     &                          ,names(3),names(4)
                           lfailure = .true.
                           return
                        endif
                        call twh_itaa(GLB_SET,imolty,jatom,iaa,param)
                     endif
                  enddo
               enddo
            enddo
c           --- safety check the number of angle-angle terms
            if ( iaa .gt. MAXAA ) then
               write(6,*) 'too many angle-angle terms in assemble'
               write(6,*) 'jatom,iaa',jatom,iaa
               lfailure = .true.
               return
            else
               call twh_inaa(GLB_SET,imolty,jatom,iaa)
            endif
         else 
c           --- no angle-angle terms
            call twh_inaa(GLB_SET,imolty,jatom,0)
         endif
      enddo

      if ( improper_assignment .eq. 'none' .or.
     &     improper_assignment .eq. 'manual' ) then
c        --- nothing to do here, but not an error
      elseif ( improper_assignment .eq. '3-bond: single improper' ) then
         do iatom = 1,totatoms
            call twh_invib(GLB_GET,imolty,iatom,iatom_nvib)
            if ( iatom_nvib .eq. 3 ) then
c              --- improper torsion found
c              --- increment improper torsion count
               call twh_inimprop(GLB_GET,imolty,iatom,itor)
               itor = itor + 1
               call twh_inimprop(GLB_SET,imolty,iatom,itor)
               do ivib = 1,iatom_nvib
                  call twh_ijvib(GLB_GET,imolty,iatom,ivib,jatom)
                  partners(ivib) = jatom
               enddo
c              --- put the types in order
               do ivib = 2,iatom_nvib
                  do jvib = 1,ivib-1
                     itest = partners(ivib)
                     jtest = partners(jvib)
                     if ( twh_get_ntype(imolty,itest) 
     &                    .lt. twh_get_ntype(imolty,jtest) ) then
c                       --- switch the atoms
                        partners(ivib) = jtest
                        partners(jvib) = itest
                     endif
                  enddo
               enddo
c              --- assign the improper partners
               jatom = partners(1)
               call twh_ijimprop2(GLB_SET,imolty,iatom,itor,jatom)
               katom = partners(2)
               call twh_ijimprop3(GLB_SET,imolty,iatom,itor,katom)
               latom = partners(3)
               call twh_ijimprop4(GLB_SET,imolty,iatom,itor,latom)
            endif
         enddo
      else
         write(6,*) 'ASSEMBLE: unknown improper_assignment'
         write(6,*) improper_assignment
         lfailure = .true.
         return
      endif

c     --- determine improper torsion types
      do iatom = 1,totatoms
         do itor = 1,twh_get_inimprop(imolty,iatom)
            call twh_ijimprop2(GLB_GET,imolty,iatom,itor,jatom)
            call twh_ijimprop3(GLB_GET,imolty,iatom,itor,katom)
            call twh_ijimprop4(GLB_GET,imolty,iatom,itor,latom)
c           --- use the atomnames
            iname = twh_get_ctname(twh_get_ntype(imolty,iatom))
            jname = twh_get_ctname(twh_get_ntype(imolty,jatom))
            kname = twh_get_ctname(twh_get_ntype(imolty,katom))
            lname = twh_get_ctname(twh_get_ntype(imolty,latom))
c           --- order of atoms has some exceptions
            if ( forcefield .eq. 'Amber96' ) then
c              --- need to find the correct type for this interaction
c              --- put the atoms back into the amber ordering
               names(1) = kname
               names(2) = lname
               names(3) = iname
               names(4) = jname
            elseif ( forcefield .eq. 'OPLS-aa' .or.
     &              forcefield .eq. 'OPLS-2001' .or.
     &              forcefield .eq. 'MMFF94' ) then
c              --- put the names into the opls and MMFF ordering
               names(1) = jname
               names(2) = iname
               names(3) = kname
               names(4) = lname
            else
c              --- send the atomnames in the standard oop order
               names(1) = iname
               names(2) = jname
               names(3) = kname
               names(4) = lname
            endif
c           --- set the param to the itimprop, this will be zero in 
c           --- all cases except those where a default value has already
c           --- been specified (notably the stereochemistry enforcer 
c           --- for proteins)
            param = twh_get_itimprop(imolty,iatom,itor)
c           --- special hack for the mgmstereochem enforcer if param -1
            if ( param .eq. -1 ) then
               type = 0
               do while ( param .eq. -1 .and. type .lt. TIMPMAX )
                  type = type + 1
                  if ( twh_get_impff(type) .eq. 'mgmstereo' ) then
                     param = type
                  endif
               enddo
c              --- safety check to make sure we have included the 
c              --- required forcefield file
               if ( param .eq. -1 ) then
                  write(6,*) 'ASSEMBLE: you need to include the '
     &                 ,'mgmstereo forcefield'
                  write(6,*) 'to enforce sterochemistry'
                  lfailure = .true.
                  return
               endif
            endif

            if ( match_style .eq. 'standard' ) then
               call twh_improper_match_standard(names,param,forcefield)
            elseif ( match_style .eq. 'MMFF checkdown' ) then
               call twh_improper_match_mmff(names,param,forcefield)
            else
               write(6,*) 'ASSEMBLE: unknown improper match_style'
               write(6,*) match_style
               lfailure = .true.
               return
            endif

            if ( param .le. 0 ) then
               write(6,*) 'parameters not found for improper'
               write(6,*) 'atoms',iatom,jatom,katom,latom
               write(6,*) 'names:'
               write(6,'(4(a10,/))') names(1),names(2),names(3),names(4)
               lfailure = .true.
               return
            endif
            call twh_itimprop(GLB_SET,imolty,iatom,itor,param)
         enddo
      enddo

c     --- determine whether we have any special one-five interactions
      do iatom = 1,totatoms
         do itor = 1,twh_get_intor(imolty,iatom)
            jatom = twh_get_ijtor2(imolty,iatom,itor)
            katom = twh_get_ijtor3(imolty,iatom,itor)
            latom = twh_get_ijtor4(imolty,iatom,itor)
            call twh_invib(GLB_GET,imolty,latom,latom_nvib)
            do ivib = 1,latom_nvib
               call twh_ijvib(GLB_GET,imolty,latom,ivib,matom)
               if ( matom .ne. katom ) then
c                 --- potential one-five interaction
                  names(1) = twh_get_nbname(twh_get_ntype(imolty,iatom))
                  names(2) = twh_get_nbname(twh_get_ntype(imolty,jatom))
                  names(3) = twh_get_nbname(twh_get_ntype(imolty,katom))
                  names(4) = twh_get_nbname(twh_get_ntype(imolty,latom))
                  names(5) = twh_get_nbname(twh_get_ntype(imolty,matom))
c                 --- initialize param
                  param = 0
c                 --- check for exact matches
                  type = 0
                  do while ( param .eq. 0 .and. type .lt. TOFMAX)
                     type = type + 1
                     if ( twh_get_offf(type) .eq. forcefield ) then
                        if ( (
     &                       twh_get_ofnames(type,1) .eq. names(1) .and.
     &                       twh_get_ofnames(type,2) .eq. names(2) .and.
     &                       twh_get_ofnames(type,3) .eq. names(3) .and.
     &                       twh_get_ofnames(type,4) .eq. names(4) .and.
     &                       twh_get_ofnames(type,5) .eq. names(5)
     &                       ) .or. (
     &                       twh_get_ofnames(type,1) .eq. names(5) .and.
     &                       twh_get_ofnames(type,2) .eq. names(4) .and.
     &                       twh_get_ofnames(type,3) .eq. names(3) .and.
     &                       twh_get_ofnames(type,4) .eq. names(2) .and.
     &                       twh_get_ofnames(type,5) .eq. names(1)
     &                       ) ) then
                           param = type
                        endif
                     endif
                  enddo
c                 --- check for exact end matches and triple wild center
                  type = 0
                  do while ( param .eq. 0 .and. type .lt. TOFMAX)
                     type = type + 1
                     if ( twh_get_offf(type) .eq. forcefield ) then
                        if ( (
     &                       twh_get_ofnames(type,1) .eq. names(1) .and.
     &                       twh_get_ofnames(type,2) .eq. 'wild' .and.
     &                       twh_get_ofnames(type,3) .eq. 'wild' .and.
     &                       twh_get_ofnames(type,4) .eq. 'wild' .and.
     &                       twh_get_ofnames(type,5) .eq. names(5)
     &                       ) .or. (
     &                       twh_get_ofnames(type,1) .eq. names(5) .and.
     &                       twh_get_ofnames(type,2) .eq. 'wild' .and.
     &                       twh_get_ofnames(type,3) .eq. 'wild' .and.
     &                       twh_get_ofnames(type,4) .eq. 'wild' .and.
     &                       twh_get_ofnames(type,5) .eq. names(1)
     &                       ) ) then
                           param = type
                        endif
                     endif
                  enddo
c                 --- we only count this if we found a non-zero param
                  if ( param .ne. 0 ) then
c                    --- increment inof
                     call twh_inof(GLB_INCR,imolty,iatom,1)
                     iof = twh_get_inof(imolty,iatom)
                     call twh_ijof5(GLB_SET,imolty,iatom,iof,matom)
                     call twh_itof(GLB_SET,imolty,iatom,iof,param)
                  endif
               endif
            enddo
         enddo
      enddo

c     --- reset nunit
      call twh_nunit(GLB_SET,imolty,totatoms)

c      write(6,*) 'finished ASSEMBLE'

      return
      end


      subroutine twh_find_bond_type(lfailure,totatoms,imolty,iatom
     &     ,jatom,ivib,param,forcefield)
c     ******************************************************************
c     * finds the bond type by searching for a name match              *
c     *                                                                *
c     * split from assemble 05-02-2006 by M.G. Martin                  *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_CBNAME
#define FUNCTION_GET_NTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer totatoms,imolty,iatom,jatom,ivib,param
      character*(FFNAMELEN) forcefield
c     --- local variables
c     --- character string scalars
      character*(FFNAMELEN) vorder,testorder
c     --- character string arrays
      character*(FFNAMELEN) names
      dimension names(2)
c     --- logical scalars
      logical lfound
c     --- integer scalars
      integer jvib,nvib,itest
      
c     --- safety check
      if ( jatom .gt. totatoms ) then
         write(6,*) 'ASSEMBLE: bond partner invalid'
         write(6,*) 'iatom,ivib,jatom',iatom,ivib,jatom
         write(6,*) 'totatoms',totatoms
         write(6,*) 'It is likely that you have either set nunit'
     &        ,' too small'
         write(6,*) 'or entered an invalid bond partner'
         lfailure = .true.
         return
      endif
c     --- use the bonding names
      names(1) = twh_get_cbname(twh_get_ntype(imolty,iatom))
      names(2) = twh_get_cbname(twh_get_ntype(imolty,jatom))
c     --- check the vibration orders
      call twh_vibration_order(GLB_GET,iatom,ivib,vorder)
c     --- find the reverse match
      jvib = 0
      call twh_invib(GLB_GET,imolty,jatom,nvib)
      lfound = .false.
      do while ( jvib .lt. nvib .and. (.not. lfound) )
         jvib = jvib + 1
         call twh_ijvib(GLB_GET,imolty,jatom,jvib,itest)
         if ( itest .eq. iatom ) lfound = .true.
      enddo
      if ( .not. lfound ) then
         write(6,*) 'FIND_BOND_TYPE: reverse bond not found'
         write(6,*) 'iatom,jatom,ivib:',iatom,jatom,ivib
         lfailure = .true.
         return
      endif
      call twh_vibration_order(GLB_GET,jatom,jvib,testorder)
      if ( vorder .ne. testorder ) then
         write(6,*) 'FIND_BOND_TYPE: vibration orders not symmetric'
         write(6,*) 'iatom,ivib,vibration order:',iatom,ivib,vorder
         write(6,*) 'jatom,jvib,vibration order:',jatom,jvib,testorder
         lfailure = .true.
         return
      endif

      param = 0

c     --- first pass through using actual vorder
      call twh_bond_matching(param,vorder,forcefield,names)
c     --- check for success
      if ( param .ne. 0 ) return

c     --- do a second pass with a general vibration order if we haven't
c     --- already
      if ( vorder .ne. 'any' ) then
         vorder = 'any'
         call twh_bond_matching(param,vorder,forcefield,names)
         if ( param .ne. 0 ) return
      endif

c     --- only get here if a match was not found
      write(6,*) 'parameters not found for bond'
      write(6,*) 'atoms:',iatom,' and:',jatom,' vib:',ivib
      write(6,*) 'names:'
      write(6,'(2(a10,/))') names(1),names(2)
      lfailure = .true.

      return
      end

      subroutine twh_bond_matching(param,vorder,forcefield,names)
c     ******************************************************************
c     * looks for a vibname pair that matches the forcefield, names,   *
c     * and vibration order                                            *
c     *                                                                *
c     * split from find_bond_type 12-10-2006 by M.G. Martin            *
c     * last modified 08-11-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NVSAME
#define FUNCTION_GET_VIBFF
#define FUNCTION_GET_VIBNAMES
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer param
      character*(*) vorder,forcefield
      character*(*) names
      dimension names(2)
c     --- local variables
c     --- character scalars
      character*(FFNAMELEN) awild,bwild,ctest
c     --- integer scalars
      integer itype,inum

c     --- sort through the list and look for complete matches
      do itype = 1,TVIBMAX
         if ( twh_get_vibff(itype) .eq. forcefield ) then
            do inum = 1,twh_get_nvsame(itype)
               if ( (twh_get_vibnames(itype,inum,1).eq. names(1) 
     &              .and. 
     &              twh_get_vibnames(itype,inum,2).eq.names(2))
     &              .or. 
     &              (twh_get_vibnames(itype,inum,1).eq.names(2) 
     &              .and. 
     &              twh_get_vibnames(itype,inum,2).eq.names(1)))
     &              then
c                 --- check the vibration order
                  if ( vorder .eq. 'any' ) then
                     param = itype
                     return
                  else
c                    --- get the vibration order from the forcefield
                     call twh_viborder(GLB_GET,itype,ctest)
                     if ( ctest .eq. vorder ) then
                        param = itype
                        return
                     endif
                  endif
               endif
            enddo
         endif
      enddo
c     --- no complete matches found - look for letter wild matches
      awild = ' wild'
      bwild = ' wild'
      do itype = 1,TVIBMAX
         if ( twh_get_vibff(itype) .eq. forcefield ) then
            do inum = 1,twh_get_nvsame(itype)
               awild(1:1) = names(1)(1:1)
               bwild(1:1) = names(2)(1:1)
               if ( (
     &              twh_get_vibnames(itype,inum,1) .eq. awild
     &              .and. 
     &              twh_get_vibnames(itype,inum,2).eq. bwild 
     &              )
     &              .or. 
     &              ( 
     &              twh_get_vibnames(itype,inum,1) .eq. bwild
     &              .and. 
     &              twh_get_vibnames(itype,inum,2) .eq. awild 
     &              )
     &              ) then
                  param = itype
                  return
               endif
            enddo
         endif
      enddo

      return
      end

      subroutine twh_angle_match_standard(lfailure,imolty,iatom,jatom
     &     ,katom,iorder,param,forcefield)
c     *****************************************************************
c     * performs a search of the angle parameters and attempts to find*
c     * the best match according to the standard procedure            *
c     * this starts by looking for exact matches, and then runs       *
c     * through a series of letter-wild cards and the full wildcards  *
c     *                                                               *
c     * split from assemble 12-15-2006 by M.G. Martin                 *
c     * last modified 10-13-2008 by M.G. Martin                       *
c     *****************************************************************
      implicit none
#define FUNCTION_GET_BENDFF
#define FUNCTION_GET_BENDNAMES
#define FUNCTION_GET_CANAME
#define FUNCTION_GET_NBSAME
#define FUNCTION_GET_NTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty,iatom,jatom,katom,iorder,param
      character*(*) forcefield
c     --- local variables
c     --- logical scalars
      logical lfailure
c     --- character string scalars
      character*(FFNAMELEN) awild,zwild,cwild
c     --- character string arrays
      character*(FFNAMELEN) names
      dimension names(3)
c     --- integer scalars
      integer type,inum

c     --- use the angle names for this type
      names(1) = twh_get_caname(twh_get_ntype(imolty,iatom))
      names(2) = twh_get_caname(twh_get_ntype(imolty,jatom))
      names(3) = twh_get_caname(twh_get_ntype(imolty,katom))
c     --- initialize type (done for each kind of search)
      type = 0
c     --- find explicit matches
      do while( param .eq. 0 .and. type .lt. TBENMAX )
         type = type + 1
         if ( forcefield .eq. twh_get_bendff(type) ) then
            do inum = 1,twh_get_nbsame(type)
c              --- need center match
               if ( twh_get_bendnames(type,inum,2) .eq.names(2) )then
c                 --- find exact end matches
                  if ( (twh_get_bendnames(type,inum,1) .eq. names(1)
     &                 .and.
     &                 twh_get_bendnames(type,inum,3).eq.names(3))
     &                 ) then
                     iorder = ORD_FORWARD
                  elseif ( (twh_get_bendnames(type,inum,1) .eq. names(3)
     &                    .and.
     &                    twh_get_bendnames(type,inum,3) .eq. names(1) )
     &                    ) then
                     iorder = ORD_REVERSE
                  endif
                  if ( iorder .ne. ORD_NULL ) then
                     param = type
                  endif
               endif
            enddo
         endif
      enddo
      
c     --- find single letter-wild card matches
c     --- initialize type (done for each kind of search)
      type = 0
c     --- initialize the awild and zwild variables so
c     --- that they end in wild
      awild = ' wild'
      zwild = ' wild'
      do while( param .eq. 0 .and. type .lt. TBENMAX )
         type = type + 1
         if ( forcefield .eq. twh_get_bendff(type) ) then
            do inum = 1,twh_get_nbsame(type)
c              --- assign letter wilds
               awild(1:1) = names(1)(1:1)
               zwild(1:1) = names(3)(1:1)
c              --- need center match
               if ( twh_get_bendnames(type,inum,2)
     &              .eq.names(2))then
c                 --- find exact end matches
                  if ( 
c                      --- awild - match - match
     &                 ( 
     &                 twh_get_bendnames(type,inum,1) .eq. awild
     &                 .and. 
     &                 twh_get_bendnames(type,inum,3) .eq. names(3) 
     &                 )
     &                 .or.
c                      --- match - match - zwild
     &                 ( 
     &                 twh_get_bendnames(type,inum,1) .eq. names(1)
     &                 .and. 
     &                 twh_get_bendnames(type,inum,3) .eq. zwild
     &                 ) ) then
                  iorder = ORD_FORWARD
               elseif ( 
c                      --- zwild - match - match
     &                 (
     &                 twh_get_bendnames(type,inum,1) .eq. zwild
     &                 .and. 
     &                 twh_get_bendnames(type,inum,3) .eq. names(1)
     &                 )
     &                 .or.
c                      --- match - match - awild
     &                 (
     &                 twh_get_bendnames(type,inum,1) .eq. names(3)
     &                 .and. 
     &                 twh_get_bendnames(type,inum,3) .eq. awild
     &                 ) ) then
                  iorder = ORD_REVERSE
               endif
               if ( iorder .ne. ORD_NULL ) then
                  param = type
               endif
            endif
         enddo
      endif
      enddo
      
c     --- find single global-wild card matches
c     --- initialize type (done for each kind of search)
      type = 0
      do while( param .eq. 0 .and. type .lt. TBENMAX )
         type = type + 1
         if ( forcefield .eq. twh_get_bendff(type) ) then
            do inum = 1,twh_get_nbsame(type)
c              --- need center match
               if ( twh_get_bendnames(type,inum,2) .eq. names(2) )then
c                 --- find exact end matches
                  if ( 
c                      --- wild - match - match
     &                 ( 
     &                 twh_get_bendnames(type,inum,1) .eq. 'wild'
     &                 .and. 
     &                 twh_get_bendnames(type,inum,3) .eq. names(3) 
     &                 )
     &                 .or.
c                      --- match - match - wild
     &                 ( 
     &                 twh_get_bendnames(type,inum,1) .eq. names(1)
     &                 .and. 
     &                 twh_get_bendnames(type,inum,3) .eq. 'wild'
     &                 ) ) then
                  iorder = ORD_FORWARD
                  param = type
               endif
            endif
         enddo
      endif
      enddo

c     --- find double end letter-wild card matches
c     --- initialize type (done for each kind of search)
      type = 0
c     --- initialize the awild and zwild variables so
c     --- that they end in wild
      awild = ' wild'
      zwild = ' wild'
      do while( param .eq. 0 .and. type .lt. TBENMAX )
         type = type + 1
         if ( forcefield .eq. twh_get_bendff(type) ) then
            do inum = 1,twh_get_nbsame(type)
c              --- assign letter wilds
               awild(1:1) = names(1)(1:1)
               zwild(1:1) = names(3)(1:1)
c              --- need center match
               if ( twh_get_bendnames(type,inum,2) .eq.names(2) )then
c                 --- find exact end matches
                  if ( 
c                      --- awild - match - zwild
     &                 ( 
     &                 twh_get_bendnames(type,inum,1) .eq. awild
     &                 .and. 
     &                 twh_get_bendnames(type,inum,3) .eq. zwild
     &                 ) ) then
                  iorder = ORD_FORWARD
               elseif (
c                      --- zwild - match - awild
     &                 (
     &                 twh_get_bendnames(type,inum,1) .eq. zwild
     &                 .and. 
     &                 twh_get_bendnames(type,inum,3) .eq. awild
     &                 ) ) then
                  iorder = ORD_REVERSE
               endif
               if ( iorder .ne. ORD_NULL ) then
                  param = type
               endif
            endif
         enddo
      endif
      enddo
      
c     --- find d end-center letter-wild card matches
c     --- initialize type (done for each kind of search)
      type = 0
c     --- initialize the awild and zwild variables so
c     --- that they end in wild
      awild = ' wild'
      cwild = ' wild'
      zwild = ' wild'
      do while( param .eq. 0 .and. type .lt. TBENMAX )
         type = type + 1
         if ( forcefield .eq. twh_get_bendff(type) ) then
            do inum = 1,twh_get_nbsame(type)
c              --- assign letter wilds
               awild(1:1) = names(1)(1:1)
               cwild(1:1) = names(2)(1:1)
               zwild(1:1) = names(3)(1:1)
c              --- need center match
               if ( twh_get_bendnames(type,inum,2) .eq. cwild ) then
c                 --- find end matches
                  if ( 
c                      --- awild - cwild - match
     &                 ( 
     &                 twh_get_bendnames(type,inum,1) .eq. awild
     &                 .and. 
     &                 twh_get_bendnames(type,inum,3) .eq. names(3)
     &                 )
     &                 .or.
c                      --- match - cwild - zwild
     &                 ( 
     &                 twh_get_bendnames(type,inum,1) .eq. names(1)
     &                 .and. 
     &                 twh_get_bendnames(type,inum,3) .eq. zwild
     &                 ) ) then
                  iorder = ORD_FORWARD
               elseif (
c                      --- zwild - cwild - match
     &                 (
     &                 twh_get_bendnames(type,inum,1) .eq. zwild
     &                 .and. 
     &                 twh_get_bendnames(type,inum,3) .eq. names(1)
     &                 ) 
     &                 .or.
c                      --- match - cwild - awild
     &                 (
     &                 twh_get_bendnames(type,inum,1) .eq. names(3)
     &                 .and. 
     &                 twh_get_bendnames(type,inum,3) .eq. awild
     &                 ) 
     &                 ) then
                  iorder = ORD_REVERSE
               endif 
               if ( iorder .ne. ORD_NULL ) then
                  param = type
               endif
            endif
         enddo
      endif
      enddo
      
c     --- find double global-wild card matches
c     --- initialize type (done for each kind of search)
      type = 0
      do while( param .eq. 0 .and. type .lt. TBENMAX )
         type = type + 1
         if ( forcefield .eq. twh_get_bendff(type) ) then
            do inum = 1,twh_get_nbsame(type)
c              --- need center match
               if ( twh_get_bendnames(type,inum,2) .eq. names(2))then
c                 --- wild - match - wild
                  if ( 
     &                 twh_get_bendnames(type,inum,1) .eq. 'wild'
     &                 .and. 
     &                 twh_get_bendnames(type,inum,3) .eq. 'wild'
     &                 ) then
                     iorder = ORD_FORWARD
                     param = type
                  endif
               endif
            enddo
         endif
      enddo
      
c     --- find triple global-wild card matches
c     --- initialize type (done for each kind of search)
      type = 0
      do while( param .eq. 0 .and. type .lt. TBENMAX )
         type = type + 1
         if ( forcefield .eq. twh_get_bendff(type) ) then
            do inum = 1,twh_get_nbsame(type)
c              --- need all to be wild
               if ( 
     &              twh_get_bendnames(type,inum,1) .eq. 'wild'
     &              .and. 
     &              twh_get_bendnames(type,inum,2) .eq. 'wild'
     &              .and. 
     &              twh_get_bendnames(type,inum,3) .eq. 'wild'
     &              ) then
                  iorder = ORD_FORWARD
                  param = type
               endif
            enddo
         endif
      enddo
      
      if ( param .eq. 0 ) then
         write(6,*) 'parameters not found for angle'
         write(6,*) 'molecule type:',imolty
         write(6,*) 'atoms ',iatom,jatom,katom
         write(6,*) 'names:'
         write(6,'(3(a10,/))') names(1),names(2),names(3)
         lfailure = .true.
         return
      endif

      return
      end

      subroutine twh_angle_match_mmff(lfailure,imolty,iatom,jatom,katom
     &     ,iorder,param,forcefield)
c     ******************************************************************
c     * searches for angle matches using the MMFF checkdown procedure  *
c     * that runs through a series of two-letter character codes and   *
c     * takes the first iteration that provides a complete match       *
c     *                                                                *
c     * originally written 12-15-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_BENDFF
#define FUNCTION_GET_BENDNAMES
#define FUNCTION_GET_CANAME
#define FUNCTION_GET_NBSAME
#define FUNCTION_GET_NTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer imolty,iatom,jatom,katom,iorder,param
      character*(FFNAMELEN) forcefield
c     --- local variables
c     --- integer scalars
      integer type,index,inum,count,istart,iend,angletype,batype
c     --- character string scalars
      character*(FFNAMELEN) testname
      character*15 forwardorder,reverseorder,testorder
      character*15 altforward,altreverse
c     --- character string arrays
      character*(FFNAMELEN) names
      dimension names(3)
      character*2 key,testkey
      dimension key(3),testkey(3)

c     --- get the forward mmff angle type
      call twh_mmff_angle_type(lfailure,imolty,iatom,jatom,katom
     &     ,angletype,batype)
      if ( lfailure ) return
      write(forwardorder,99) 'MMFF sb:',batype,'a:',angletype
c     --- get the reverse mmff angle type
      call twh_mmff_angle_type(lfailure,imolty,katom,jatom,iatom
     &     ,angletype,batype)
      if ( lfailure ) return
      write(reverseorder,99) 'MMFF sb:',batype,'a:',angletype
 99   format(a8,i2,1x,a2,i1)

c     --- create the alternate forward and reverse with wilcards for
c     --- the stretch-bend
      altforward = forwardorder
      altforward(9:10) = '**'
      altreverse = reverseorder
      altreverse(9:10) = '**'

c     --- use the angle names for this type
      names(1) = twh_get_caname(twh_get_ntype(imolty,iatom))
      names(2) = twh_get_caname(twh_get_ntype(imolty,jatom))
      names(3) = twh_get_caname(twh_get_ntype(imolty,katom))

      do index = 1,5
         istart = 2*index-1
         iend = 2*index
         do count = 1,3
            key(count) = names(count)(istart:iend)
         enddo

         do type = 1,TBENMAX
c           --- see if the forcefield name matches
            if ( twh_get_bendff(type) .eq. forcefield ) then
c              --- get the bendnames
               do inum = 1,twh_get_nbsame(type)
                  do count = 1,3
                     testname = twh_get_bendnames(type,inum,count)
                     testkey(count) = testname(1:2)
                  enddo
c                 --- check for a forward match
                  if ( key(1) .eq. testkey(1) .and.
     &                 key(2) .eq. testkey(2) .and.
     &                 key(3) .eq. testkey(3) ) then
c                    --- check the angleorder
                     call twh_angleorder(GLB_GET,type,testorder)
                     if ( testorder .eq. forwardorder
     &                    .or. testorder .eq. altforward ) then
c                       --- forward match
                        iorder = ORD_FORWARD
                        param = type
                        return
                     endif
                  endif
c                 --- check for a reverse match
                  if ( key(1) .eq. testkey(3) .and.
     &                 key(2) .eq. testkey(2) .and.
     &                 key(3) .eq. testkey(1) ) then
c                    --- check the angleorder
                     call twh_angleorder(GLB_GET,type,testorder)
                     if ( testorder .eq. reverseorder
     &                    .or. testorder .eq. altreverse ) then
c                       --- reverse match
                        iorder = ORD_REVERSE
                        param = type
                        return
                     endif
                  endif
               enddo
            endif
         enddo
      enddo

c     --- if we got this far then a match was not found
      write(6,*) 'parameters not found for angle'
      write(6,*) 'molecule type:',imolty
      write(6,*) 'atoms ',iatom,jatom,katom
      write(6,*) 'names:'
      write(6,'(3(a10,/))') names(1),names(2),names(3)
      write(6,*) 'forward and reverse orders'
      write(6,*) forwardorder
      write(6,*) reverseorder
      write(6,*) 'alternate forward and reverse orders'
      write(6,*) altforward
      write(6,*) altreverse
      lfailure = .true.

      return
      end

      subroutine twh_torsion_match_standard(lfailure,imolty,iatom
     &     ,jatom,katom,latom,param,forcefield,specialtor)
c     ******************************************************************
c     * finds torsion matches using the standard method of staring     *
c     * with exact matches and then looking for various wildcards      *
c     *                                                                *
c     * split from assemble 12-15-2006 by M.G. Martin                  *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_CTNAME
#define FUNCTION_GET_NTSAME
#define FUNCTION_GET_NTYPE
#define FUNCTION_GET_TORFF
#define FUNCTION_GET_TORNAMES
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer imolty,iatom,jatom,katom,latom,param
      character*(FFNAMELEN) forcefield
      character*(*) specialtor
c     --- local variables
c     --- character string scalars
      character*(FFNAMELEN) awild,zwild,cwild,dwild
      character*15 tororder
c     --- character string arrays
      character*(FFNAMELEN) names
      dimension names(4)
c     --- integer scalars
      integer type,inum

      if ( param .eq. -20 ) then
c        --- we want the all null torsion
         names(1) = 'null'
         names(2) = 'null'
         names(3) = 'null'
         names(4) = 'null'
c        --- reset param
         param = 0
      else
c        --- use the torsion names
         names(1) = twh_get_ctname(twh_get_ntype(imolty,iatom))
         names(2) = twh_get_ctname(twh_get_ntype(imolty,jatom))
         names(3) = twh_get_ctname(twh_get_ntype(imolty,katom))
         names(4) = twh_get_ctname(twh_get_ntype(imolty,latom))
      endif
      
c     --- run through the list of torsions and check for 
c     --- exact matches without wildcards
      type = 0
      do while ( param .eq. 0 .and. type .lt. TTORMAX )
         type = type + 1
         if ( twh_get_torff(type) .eq. forcefield ) then
            do inum = 1,twh_get_ntsame(type)
               if ( (twh_get_tornames(type,inum,1) .eq. names(1)
     &              .and. 
     &              twh_get_tornames(type,inum,2) .eq. names(2)
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(3)
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. names(4))
     &              .or.
     &              (twh_get_tornames(type,inum,1) .eq. names(4)
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(3)
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(2)
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. names(1)) 
     &              )then
c                 --- check torsion order
                  call twh_torsionorder(GLB_GET,type,tororder)
                  if ( tororder .eq. specialtor ) then
                     param = type
                  endif
               endif
            enddo
         endif
      enddo
c     --- look for single letter-wild matches
c     --- these have a match of the first letter of names
c     --- and tornames, while letters 2:5 of tornames are 
c     --- the characters 'wild'
      type = 0
c     --- initialize the awild and zwild variables so
c     --- that they end in wild
      awild = ' wild'
      zwild = ' wild'
      do while (param .eq. 0 .and. type .lt. TTORMAX)
         type = type + 1
         if ( twh_get_torff(type) .eq. forcefield ) then
            do inum = 1,twh_get_ntsame(type)
               awild(1:1) = names(1)(1:1)
               zwild(1:1) = names(4)(1:1)
               if ( ( twh_get_tornames(type,inum,1) .eq. awild 
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(2)
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(3)
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. names(4))
     &              .or.
     &              ( twh_get_tornames(type,inum,1) .eq. zwild
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(3)
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(2)
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. names(1))
     &              .or.
     &              (twh_get_tornames(type,inum,4) .eq. zwild
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(3)
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(2)
     &              .and.
     &              twh_get_tornames(type,inum,1) .eq. names(1))
     &              .or.
     &              (twh_get_tornames(type,inum,4) .eq. awild 
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(2)
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(3)
     &              .and.
     &              twh_get_tornames(type,inum,1) .eq. names(4))
     &              ) then
c                 --- check torsion order
                  call twh_torsionorder(GLB_GET,type,tororder)
                  if ( tororder .eq. specialtor ) then
                     param = type
                  endif
               endif
            enddo
         endif
      enddo

c     --- if nothing yet found then run through again looking for
c     --- single global-wildcard matches
      type = 0
      do while ( param .eq. 0 .and. type .lt. TTORMAX )
         type = type + 1
         if ( twh_get_torff(type) .eq. forcefield ) then
            do inum = 1,twh_get_ntsame(type)
               if ( (twh_get_tornames(type,inum,1) .eq. 'wild'
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(2)
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(3)
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. names(4))
     &              .or.
     &              (twh_get_tornames(type,inum,1) .eq. 'wild'
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(3)
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(2)
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. names(1))
     &              .or.
     &              (twh_get_tornames(type,inum,1) .eq. names(1)
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(2)
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(3)
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. 'wild')
     &              .or.
     &              (twh_get_tornames(type,inum,1) .eq. names(4)
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(3)
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(2)
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. 'wild')
     &              )then
c                 --- check torsion order
                  call twh_torsionorder(GLB_GET,type,tororder)
                  if ( tororder .eq. specialtor ) then
                     param = type
                  endif
               endif
            enddo
         endif
      enddo

c     --- look for double letter-wild matches on either end
      type = 0
      do while (param .eq. 0 .and. type .lt. TTORMAX)
         type = type + 1
         if ( twh_get_torff(type) .eq. forcefield ) then
            do inum = 1,twh_get_ntsame(type)
c              --- forward direction
               awild(1:1) = names(1)(1:1)
               zwild(1:1) = names(4)(1:1)
               if ( (twh_get_tornames(type,inum,1) .eq. awild
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(2)
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(3)
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. zwild 
     &              ) .or. ( 
     &              twh_get_tornames(type,inum,1) .eq. zwild
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(3)
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(2)
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. awild )
     &              ) then
c                 --- check torsion order
                  call twh_torsionorder(GLB_GET,type,tororder)
                  if ( tororder .eq. specialtor ) then
                     param = type
                  endif
               endif
            enddo
         endif
      enddo

c     --- look for double letter-wild matches bond on front end
      type = 0
      do while (param .eq. 0 .and. type .lt. TTORMAX)
         type = type + 1
         if ( twh_get_torff(type) .eq. forcefield ) then
            do inum = 1,twh_get_ntsame(type)
               awild(1:1) = names(1)(1:1)
               zwild(1:1) = names(2)(1:1)
               if ( ( twh_get_tornames(type,inum,1) .eq. awild
     &              .and. twh_get_tornames(type,inum,2) .eq. zwild
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(3)
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. names(4) )
     &              .or.
     &              ( twh_get_tornames(type,inum,1) .eq. names(4)
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(3)
     &              .and. twh_get_tornames(type,inum,3) .eq. zwild
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. awild ) 
     &              ) then
c                 --- check torsion order
                  call twh_torsionorder(GLB_GET,type,tororder)
                  if ( tororder .eq. specialtor ) then
                     param = type
                  endif
               endif
            enddo
         endif
      enddo

c     --- look for double letter-wild matches bond on back end
      type = 0
      do while (param .eq. 0 .and. type .lt. TTORMAX)
         type = type + 1
         if ( twh_get_torff(type) .eq. forcefield ) then
            do inum = 1,twh_get_ntsame(type)
               awild(1:1) = names(3)(1:1)
               zwild(1:1) = names(4)(1:1)
               if ( ( twh_get_tornames(type,inum,1) .eq. names(1)
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(2)
     &              .and. twh_get_tornames(type,inum,3) .eq. awild
     &              .and. twh_get_tornames(type,inum,4) .eq. zwild
     &              ) .or. (
     &              twh_get_tornames(type,inum,1) .eq. zwild
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. awild
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(2)
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. names(1) ) 
     &              ) then
c                 --- check torsion order
                  call twh_torsionorder(GLB_GET,type,tororder)
                  if ( tororder .eq. specialtor ) then
                     param = type
                  endif
               endif
            enddo
         endif
      enddo
      
c     --- if no matches were found then run 
c     --- through again looking for double global-wildcard matches
c     --- with the double wilds at the two ends
      type = 0
      do while ( param .eq. 0 .and. type .lt. TTORMAX )
         type = type + 1
         if ( twh_get_torff(type) .eq. forcefield ) then
            do inum = 1,twh_get_ntsame(type)
               if ( (twh_get_tornames(type,inum,1) .eq. 'wild'
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(2)
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(3)
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. 'wild'
     &              ) .or. (
     &              twh_get_tornames(type,inum,1) .eq. 'wild'
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(3)
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(2)
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. 'wild')
     &              ) then
c                 --- check torsion order
                  call twh_torsionorder(GLB_GET,type,tororder)
                  if ( tororder .eq. specialtor ) then
                     param = type
                  endif
               endif
            enddo
         endif
      enddo

c     --- if no matches were found then run 
c     --- through again looking for double global-wildcard matches
c     --- with the double wilds in the middle
      type = 0
      do while ( param .eq. 0 .and. type .lt. TTORMAX )
         type = type + 1
         if ( twh_get_torff(type) .eq. forcefield ) then
            do inum = 1,twh_get_ntsame(type)
               if ( (twh_get_tornames(type,inum,1) .eq. names(1)
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. 'wild'
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. 'wild'
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. names(4))
     &              .or.
     &              (twh_get_tornames(type,inum,1) .eq. names(4)
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. 'wild'
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. 'wild'
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. names(1)))
     &              then
c                 --- check torsion order
                  call twh_torsionorder(GLB_GET,type,tororder)
                  if ( tororder .eq. specialtor ) then
                     param = type
                  endif
               endif
            enddo
         endif
      enddo

c     --- if no matches were found then run 
c     --- through again looking for double global-wildcard matches
c     --- with the double wilds at one end
      type = 0
      do while ( param .eq. 0 .and. type .lt. TTORMAX )
         type = type + 1
         if ( twh_get_torff(type) .eq. forcefield ) then
            do inum = 1,twh_get_ntsame(type)
               if ( (twh_get_tornames(type,inum,1) .eq. 'wild'
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. 'wild'
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(3)
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. names(4)
     &              ) .or. (
     &              twh_get_tornames(type,inum,1) .eq. names(4)
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(3)
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. 'wild'
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. 'wild')
     &              ) then
c                 --- check torsion order
                  call twh_torsionorder(GLB_GET,type,tororder)
                  if ( tororder .eq. specialtor ) then
                     param = type
                  endif
               endif
            enddo
         endif
      enddo
      
c     --- look for triple letter-wild matches where the complete
c     --- match is one of the central atoms
      type = 0
c     --- setup the letter wilds
      awild = ' wild'
      cwild = ' wild'
      dwild = ' wild'
      zwild = ' wild'
      do while ( param .eq. 0 .and. type .lt. TTORMAX )
         type = type + 1
         if ( twh_get_torff(type) .eq. forcefield ) then
            do inum = 1,twh_get_ntsame(type)
               awild(1:1) = names(1)(1:1)
               cwild(1:1) = names(2)(1:1)
               dwild(1:1) = names(3)(1:1)
               zwild(1:1) = names(4)(1:1)
               if (
     &              ( twh_get_tornames(type,inum,1) .eq. awild
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(2)
     &              .and. twh_get_tornames(type,inum,3) .eq. dwild
     &              .and. twh_get_tornames(type,inum,4) .eq. zwild
     &              ) .or. (
     &              twh_get_tornames(type,inum,1) .eq. awild
     &              .and. twh_get_tornames(type,inum,2) .eq. cwild
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(3)
     &              .and. twh_get_tornames(type,inum,4) .eq. zwild
     &              ) .or. (
     &              twh_get_tornames(type,inum,4) .eq. awild
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. names(2)
     &              .and. twh_get_tornames(type,inum,2) .eq. dwild
     &              .and. twh_get_tornames(type,inum,1) .eq. zwild
     &              ) .or. (
     &              twh_get_tornames(type,inum,4) .eq. awild
     &              .and. twh_get_tornames(type,inum,3) .eq. cwild
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. names(3)
     &              .and. twh_get_tornames(type,inum,1) .eq. zwild
     &              )
     &              ) then
c                 --- check torsion order
                  call twh_torsionorder(GLB_GET,type,tororder)
                  if ( tororder .eq. specialtor ) then
                     param = type
                  endif
               endif
            enddo
         endif
      enddo

c     --- if no matches were found then run through again and 
c     --- see if this forcefield has a global-wildcard
      type = 0
      do while ( param .eq. 0 .and. type .lt. TTORMAX )
         type = type + 1
         if ( twh_get_torff(type) .eq. forcefield ) then
            do inum = 1,twh_get_ntsame(type)
               if ( twh_get_tornames(type,inum,1) .eq. 'wild'
     &              .and.
     &              twh_get_tornames(type,inum,2) .eq. 'wild'
     &              .and.
     &              twh_get_tornames(type,inum,3) .eq. 'wild'
     &              .and.
     &              twh_get_tornames(type,inum,4) .eq. 'wild'
     &              ) then
c                 --- check torsion order
                  call twh_torsionorder(GLB_GET,type,tororder)
                  if ( tororder .eq. specialtor ) then
                     param = type
                  endif
               endif
            enddo
         endif
      enddo
      
      if ( param .eq. 0  ) then
         write(6,*) 'parameters not found for torsion'
         write(6,*) 'atoms ',iatom,jatom,katom,latom
         write(6,*) 'names:'
         write(6,'(4(a10,/))') names(1),names(2),names(3),names(4)
         lfailure = .true.
         return
      endif

      return
      end

      subroutine twh_torsion_match_mmff(lfailure,imolty,iatom
     &     ,jatom,katom,latom,param,forcefield)
c     ******************************************************************
c     * uses the MMFF checkdown procedure to find a torsion match      *
c     *                                                                *
c     * originally written 12-15-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_CTNAME
#define FUNCTION_GET_NTSAME
#define FUNCTION_GET_NTYPE
#define FUNCTION_GET_TORFF
#define FUNCTION_GET_TORNAMES
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer imolty,iatom,jatom,katom,latom,param
      character*(FFNAMELEN) forcefield

c     --- local variables
c     --- integer scalars
      integer type,index,inum,count,istart,iend
c     --- character string scalars
      character*2 wildcard
      character*(FFNAMELEN) testname
      character*15 forwardorder,reverseorder,testorder
c     --- character string arrays
      character*(FFNAMELEN) names
      dimension names(4)
      character*2 key,testkey
      dimension key(4),testkey(4)

      call twh_mmff_torsion_type(lfailure,imolty,iatom,jatom,katom,latom
     &     ,type)
      if ( lfailure ) return
      write(forwardorder,55) type
      call twh_mmff_torsion_type(lfailure,imolty,latom,katom,jatom,iatom
     &     ,type)
      if ( lfailure ) return
      write(reverseorder,55) type
 55   format('MMFF t:',i1)

c     --- use the torsion names for this type
      names(1) = twh_get_ctname(twh_get_ntype(imolty,iatom))
      names(2) = twh_get_ctname(twh_get_ntype(imolty,jatom))
      names(3) = twh_get_ctname(twh_get_ntype(imolty,katom))
      names(4) = twh_get_ctname(twh_get_ntype(imolty,latom))

      do index = 1,5
         istart = 2*index-1
         iend = 2*index
         do count = 1,4
            key(count) = names(count)(istart:iend)
         enddo

         do type = 1,TTORMAX
c           --- see if the forcefield name matches
            if ( twh_get_torff(type) .eq. forcefield ) then
c              --- get the tornames
               do inum = 1,twh_get_ntsame(type)
                  do count = 1,4
                     testname = twh_get_tornames(type,inum,count)
                     testkey(count) = testname(1:2)
                  enddo
c                 --- check for a forward match
                  if ( key(1) .eq. testkey(1) .and.
     &                 key(2) .eq. testkey(2) .and.
     &                 key(3) .eq. testkey(3) .and.
     &                 key(4) .eq. testkey(4) ) then
c                    --- check torsion order
                     call twh_torsionorder(GLB_GET,type,testorder)
                     if ( testorder .eq. forwardorder ) then
c                       --- forward match
                        param = type
                        return
                     endif
                  endif
c                 --- check for a reverse match
                  if ( key(1) .eq. testkey(4) .and.
     &                 key(2) .eq. testkey(3) .and.
     &                 key(3) .eq. testkey(2) .and.
     &                 key(4) .eq. testkey(1) ) then
c                    --- check torsion order
                     call twh_torsionorder(GLB_GET,type,testorder)
                     if ( testorder .eq. reverseorder ) then
c                       --- reverse match
                        param = type
                        return
                     endif
                  endif
               enddo
            endif
         enddo
      enddo

      index = 1
      istart = 2*index-1
      iend = 2*index
      do count = 1,4
         key(count) = names(count)(istart:iend)
      enddo

      wildcard = ' 0'
      do type = 1,TTORMAX
c        --- see if the forcefield name matches
         if ( twh_get_torff(type) .eq. forcefield ) then
c            --- get the tornames
            do inum = 1,twh_get_ntsame(type)
               do count = 1,4
                  testname = twh_get_tornames(type,inum,count)
                  testkey(count) = testname(1:2)
               enddo
c              --- check for a first position single-wild forward match
               if ( wildcard .eq. testkey(1) .and.
     &              key(2) .eq. testkey(2) .and.
     &              key(3) .eq. testkey(3) .and.
     &              key(4) .eq. testkey(4) ) then
c                 --- forward match
                  param = type
                  return
               endif
c              --- check for a final position single-wild forward match
               if ( key(1) .eq. testkey(1) .and.
     &              key(2) .eq. testkey(2) .and.
     &              key(3) .eq. testkey(3) .and.
     &              wildcard .eq. testkey(4) ) then
c                 --- forward match
                  param = type
                  return
               endif
c              --- check for a first position single-wild reverse match
               if ( wildcard .eq. testkey(4) .and.
     &              key(2) .eq. testkey(3) .and.
     &              key(3) .eq. testkey(2) .and.
     &              key(4) .eq. testkey(1) ) then
c                 --- reverse match
                  param = type
                  return
               endif
c              --- check for a final position single-wild reverse match
               if ( key(1) .eq. testkey(4) .and.
     &              key(2) .eq. testkey(3) .and.
     &              key(3) .eq. testkey(2) .and.
     &              wildcard .eq. testkey(1) ) then
c                 --- reverse match
                  param = type
                  return
               endif
            enddo
         endif
      enddo

c     --- check for double-wildcard matches
      do type = 1,TTORMAX
c        --- see if the forcefield name matches
         if ( twh_get_torff(type) .eq. forcefield ) then
c            --- get the tornames
            do inum = 1,twh_get_ntsame(type)
               do count = 1,4
                  testname = twh_get_tornames(type,inum,count)
                  testkey(count) = testname(1:2)
               enddo
c              --- check for a double-wild forward match
               if ( wildcard .eq. testkey(1) .and.
     &              key(2) .eq. testkey(2) .and.
     &              key(3) .eq. testkey(3) .and.
     &              wildcard .eq. testkey(4) ) then
c                 --- forward match
                  param = type
                  return
               endif
c              --- check for a double-wild reverse match
               if ( wildcard .eq. testkey(4) .and.
     &              key(2) .eq. testkey(3) .and.
     &              key(3) .eq. testkey(2) .and.
     &              wildcard .eq. testkey(1) ) then
c                 --- reverse match
                  param = type
                  return
               endif
            enddo
         endif
      enddo

c     --- if we got this far then a match was not found
      write(6,*) 'parameters not found for torsion'
      write(6,*) 'molecule type:',imolty
      write(6,*) 'atoms ',iatom,jatom,katom,latom
      write(6,*) 'names:'
      write(6,'(4(a10,/))') names(1),names(2),names(3),names(4)
      write(6,*) 'forward and reverse order'
      write(6,*) forwardorder
      write(6,*) reverseorder

      lfailure = .true.

      return
      end


      subroutine twh_mmff_bond_type(lfailure,imolty,iatom,jatom,bij)
c     ******************************************************************
c     * determines the integer code for the MMFF bond order            *
c     *                                                                *
c     * originally written 12-15-2006 by M.G. Martin                   *
c     * last modified 12-15-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer imolty,iatom,jatom,bij
c     --- local variables
c     --- integer scalars
      integer ivib,invib,itest,type
c     --- character string scalars
      character*(FFNAMELEN) order

c     --- search from the second atom
      call twh_invib(GLB_GET,imolty,jatom,invib)
      do ivib = 1,invib
         call twh_ijvib(GLB_GET,imolty,jatom,ivib,itest)
         if ( itest .eq. iatom ) then
c           --- match found: get the vibtype
            call twh_itvib(GLB_GET,imolty,jatom,ivib,type)
c           --- get the viborder for this type
            call twh_viborder(GLB_GET,type,order)
            if ( order .eq. 'single' ) then
               bij = 1
               return
            elseif ( order .eq. 'wild' ) then
               bij = 0
               return
            else
               write(6,*) 'MMFF_BOND_TYPE: unknown viborder'
               write(6,*) order
               lfailure = .true.
               return
            endif
         endif
      enddo

c     --- if we made it this far there is a problem
      write(6,*) 'MMFF_BOND_TYPE: bond partner not found'
      write(6,*) 'imolty,iatom,jatom',imolty,iatom,jatom
      lfailure = .true.
      return
      end

      subroutine twh_mmff_angle_type(lfailure,imolty,iatom,jatom
     &     ,katom,angletype,batype)
c     ******************************************************************
c     * determine the mmff angle, and stretch-bend types               *
c     *                                                                *
c     * angle type rules                                               *
c     * 0: The angle i-j-k is a "normal" bond angle                    *
c     * 1: Either bond i-j or bond j-k has a bond type of 1            *
c     * 2: Both i-j and j-k have bond types of 1; the sum is 2.        *
c     * 3: The angle occurs in a three-membered ring                   *
c     * 4: The angle occurs in a four-membered ring                    *
c     * 5: Is in a three-membered ring and the sum of the bond types   *
c     *    is 1                                                        *
c     * 6: Is in a three-membered ring and the sum of the bond types   *
c     *    is 2                                                        *
c     * 7: Is in a four-membered ring and the sum of the bond types    *
c     *    is 1                                                        *
c     * 8: Is in a four-membered ring and the sum of the bond types    *
c     *    is 2                                                        *
c     *                                                                *
c     * bond-angle type rules                                          *
c     * 0: angle:0 bij: 0 bjk: 0                                       *
c     * 1: angle:1 bij: 1 bjk: 0                                       *
c     * 2: angle:1 bij: 0 bjk: 1                                       *
c     * 3: angle:2 bij: 1 bjk: 1                                       *
c     * 4: angle:4 bij: 0 bjk: 0                                       *
c     * 5: angle:3 bij: 0 bjk: 0                                       *
c     * 6: angle:5 bij:	1 bjk: 0                                       *
c     * 7: angle:5 bij: 0 bjk: 1                                       *
c     * 8: angle:6 bij: 1 bjk: 1                                       *
c     * 9: angle:7 bij:	1 bjk: 0                                       *
c     * 10: angle:7 bij: 0 bjk: 1                                      *
c     * 11: angle:8 bij: 1 bkj: 1                                      *
c     *                                                                *
c     * originally written 12-15-2006 by M.G. Martin                   *
c     * last modified 12-15-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer imolty,iatom,jatom,katom,angletype,batype
c     --- local variables
c     --- integer scalars
      integer invib,ivib,bondsum,itest,ibend,inbend,jtest,bij,bjk

c     --- compute the sum of the bond types
      bondsum = 0
      call twh_mmff_bond_type(lfailure,imolty,iatom,jatom,bij)
      if ( lfailure ) return
      call twh_mmff_bond_type(lfailure,imolty,jatom,katom,bjk)
      if ( lfailure ) return
      bondsum = bij + bjk

c     --- determine if this is a four-membered ring
      call twh_inben(GLB_GET,imolty,katom,inbend)
      do ibend = 1,inbend
         call twh_ijben3(GLB_GET,imolty,katom,ibend,itest)
         if ( itest .eq. iatom ) then
c           --- make sure this is a different angle
            call twh_ijben2(GLB_GET,imolty,katom,ibend,jtest)
            if ( jtest .ne. jatom ) then
c              --- four-membered ring
               if ( bondsum .eq. 2 ) then
c                 --- four-membered ring with a bondsum of 2
                  angletype = 8
                  batype = 11
               elseif ( bondsum .eq. 1 ) then
c                 --- four-membered ring with a bondsum of 1
                  angletype = 7
                  if ( bij .eq. 1 ) then
                     batype = 9
                  else
                     batype = 10
                  endif
               else
c                 --- four-membered ring with a bondsum of 0
                  angletype = 4
                  batype = 4
               endif
               return
            endif
         endif
      enddo

c     --- determine if this is a three-membered ring
      call twh_invib(GLB_GET,imolty,katom,invib)
      do ivib = 1,invib
         call twh_ijvib(GLB_GET,imolty,katom,ivib,itest)
         if ( itest .eq. iatom ) then
c           --- three-membered ring
            if ( bondsum .eq. 2 ) then
c              --- three-membered ring with a bondsum of 2
               angletype = 6
               batype = 8
            elseif ( bondsum .eq. 1 ) then
c              --- three-membered ring with a bondsum of 1
               angletype = 5
               if ( bij .eq. 1 ) then
                  batype = 6
               else
                  batype = 7
               endif
            else
c              --- three-membered ring with a bondsum of 0
               angletype = 3
               batype = 5
            endif
            return
         endif
      enddo

c     --- not a three or four membered ring
      angletype = bondsum
      if ( angletype .eq. 2 ) then
         batype = 3
      elseif ( angletype .eq. 1 ) then
         if ( bij .eq. 1 ) then
            batype = 1
         else
            batype = 2
         endif
      else
         batype = 0
      endif

      return
      end

      subroutine twh_mmff_torsion_type(lfailure,imolty,iatom,jatom,katom
     &     ,latom,type)
c     ******************************************************************
c     * determines the MMFF torsion type                               *
c     *                                                                *
c     * 0: if not one of the other values                              *
c     * 1: if the j-k bond has a bond type index of 1                  *
c     * 2: if j-k is 0, but i-j and/or k-l is 1                        *
c     * 4: when all four atoms make up a the same 4-membered ring      *
c     * 5: when all four atoms make up a 5-membered ring that is not   *
c     *    aromatic and contains no unsaturation                       *
c     * last modified 10-14-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer imolty,iatom,jatom,katom,latom,type
c     --- local variables
c     --- logical scalars
      logical lvalid
c     --- character string scalars
      character*(FFNAMELEN) order
c     --- integer scalars
      integer invib,ivib,itest,inbend,ibend,matom,bij,bjk,bkl

c     --- check for four-membered ring
      call twh_invib(GLB_GET,imolty,latom,invib)
      do ivib = 1,invib
         call twh_ijvib(GLB_GET,imolty,latom,ivib,itest)
         if ( itest .eq. iatom ) then
c           --- four-membered ring
            type = 4
            return
         endif
      enddo

c     --- check for five-membered ring
      call twh_inben(GLB_GET,imolty,latom,inbend)
      do ibend = 1,inbend
         call twh_ijben3(GLB_GET,imolty,latom,ibend,itest)
         if ( itest .eq. iatom ) then
c           --- five-membered ring
            call twh_ijben2(GLB_GET,imolty,latom,ibend,matom)
c           --- make sure it is all single-bonds
            lvalid = .true.
            call twh_bond_order(lfailure,imolty,iatom,jatom,order)
            if ( lfailure ) return
            if ( order .ne. 'single' ) lvalid = .false.
            call twh_bond_order(lfailure,imolty,jatom,katom,order)
            if ( lfailure ) return
            if ( order .ne. 'single' ) lvalid = .false.
            call twh_bond_order(lfailure,imolty,katom,latom,order)
            if ( lfailure ) return
            if ( order .ne. 'single' ) lvalid = .false.
            call twh_bond_order(lfailure,imolty,latom,matom,order)
            if ( lfailure ) return
            if ( order .ne. 'single' ) lvalid = .false.
            call twh_bond_order(lfailure,imolty,matom,iatom,order)
            if ( lfailure ) return
            if ( order .ne. 'single' ) lvalid = .false.
            if ( lvalid ) then
               type = 5
               return
            endif
         endif
      enddo

c     --- check the j-k bond index from the forcefield
      call twh_mmff_bond_type(lfailure,imolty,jatom,katom,bjk)
      if ( lfailure ) return
      if ( bjk .eq. 1 ) then
         type = 1
         return
      endif

      if ( bjk .eq. 0 ) then
c        --- check bij and bkl to see if either of them are 1
         call twh_mmff_bond_type(lfailure,imolty,iatom,jatom,bij)
         if ( lfailure ) return
         call twh_mmff_bond_type(lfailure,imolty,katom,latom,bkl)
         if ( lfailure ) return
         if ( bij .eq. 1 .or. bkl .eq. 1 ) then
            type = 2
            return
         endif
      endif

c     --- if we got this far then it is a default type of 0
      type = 0

      return
      end

      subroutine twh_bond_order(lfailure,imolty,testa,testb,order)
c     ******************************************************************
c     * returns the vibration order between two atoms in a molecule as *
c     * specified in the towhee_input setup                            *
c     *                                                                *
c     * originally written 12-16-2006 by M.G. Martin                   *
c     * last modified 12-16-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer imolty,testa,testb
      character*(*) order
c     --- local variables
      integer invib,ivib,itest

      call twh_invib(GLB_GET,imolty,testa,invib)
      do ivib = 1,invib
         call twh_ijvib(GLB_GET,imolty,testa,ivib,itest)
         if ( itest .eq. testb ) then
c           --- get vibration order
            call twh_vibration_order(GLB_GET,testa,ivib,order)
            return
         endif
      enddo

c     --- if we got this far a match was not found
      write(6,*) 'BOND_ORDER: match not found'
      write(6,*) 'imolty,testa,testb',imolty,testa,testb
      lfailure = .true.
      return
      end

      subroutine twh_improper_match_standard(names,param,forcefield)
c     ******************************************************************
c     * matches up the improper parameters to determine the improper   *
c     * type                                                           *
c     *                                                                *
c     * split from assemble 12-20-2006 by M.G. Martin                  *
c     * last modified 12-20-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IMPFF
#define FUNCTION_GET_IMPNAMES
#define FUNCTION_GET_NIMPSAME
#include "functions.h"
c     --- variables passed to/from the subroutine
      character*(FFNAMELEN) names
      dimension names(4)
      integer param
      character*(FFNAMELEN) forcefield
c     --- local variables
      integer type,inum

c     --- check for exact matches
      type = 0
      do while ( param .eq. 0 .and. type .lt. TIMPMAX)
         type = type + 1
         if ( twh_get_impff(type) .eq. forcefield ) then
            do inum = 1,twh_get_nimpsame(type)
               if ( (twh_get_impnames(type,inum,1) .eq. names(1)
     &              .and. 
     &              twh_get_impnames(type,inum,2).eq.names(2)
     &              .and.
     &              twh_get_impnames(type,inum,3).eq.names(3)
     &              .and.
     &              twh_get_impnames(type,inum,4).eq.names(4)
     &              ) .or. (
     &              twh_get_impnames(type,inum,1).eq.names(4)
     &              .and.
     &              twh_get_impnames(type,inum,2).eq.names(3)
     &              .and.
     &              twh_get_impnames(type,inum,3).eq.names(2)
     &              .and.
     &              twh_get_impnames(type,inum,4).eq.names(1)
     &              ) ) then
                  param = type
               endif
            enddo
         endif
      enddo
c     --- check for single wildcard matches in the first spot
      type = 0
      do while ( param .eq. 0 .and. type .lt. TIMPMAX)
         type = type + 1
         if ( twh_get_impff(type) .eq. forcefield ) then
            do inum = 1,twh_get_nimpsame(type)
               if ( twh_get_impnames(type,inum,1) .eq. 'wild' 
     &              .and.
     &              twh_get_impnames(type,inum,2) .eq. names(2)
     &              .and.
     &              twh_get_impnames(type,inum,3) .eq. names(3)
     &              .and.
     &              twh_get_impnames(type,inum,4) .eq. names(4))
     &              then
                  param = type
               endif
            enddo
         endif
      enddo
c     --- check for double wildcard matches in the first 2 spots
      type = 0
      do while ( param .eq. 0 .and. type .lt. TIMPMAX )
         type = type + 1
         if ( twh_get_impff(type) .eq. forcefield ) then
            do inum = 1,twh_get_nimpsame(type)
               if ( twh_get_impnames(type,inum,1) .eq. 'wild' 
     &              .and.
     &              twh_get_impnames(type,inum,2) .eq. 'wild'
     &              .and.
     &              twh_get_impnames(type,inum,3) .eq. names(3)
     &              .and.
     &              twh_get_impnames(type,inum,4) .eq. names(4)
     &              ) then
                  param = type
               endif
            enddo
         endif
      enddo
c     --- check for double wildcard matches in places 2 and 3
      type = 0
      do while ( param .eq. 0 .and. type .lt. TIMPMAX)
         type = type + 1
         if ( twh_get_impff(type) .eq. forcefield ) then
            do inum = 1,twh_get_nimpsame(type)
               if ( twh_get_impnames(type,inum,1) .eq. names(1) 
     &              .and.
     &              twh_get_impnames(type,inum,2) .eq. 'wild'
     &              .and.
     &              twh_get_impnames(type,inum,3) .eq. 'wild'
     &              .and.
     &              twh_get_impnames(type,inum,4).eq.names(4)
     &              )then
                  param = type
               endif
            enddo
         endif
      enddo
c     --- check for double wildcard matches in places 3 and 4
      type = 0
      do while ( param .eq. 0 .and. type .lt. TIMPMAX)
         type = type + 1
         if ( twh_get_impff(type) .eq. forcefield ) then
            do inum = 1,twh_get_nimpsame(type)
               if ( twh_get_impnames(type,inum,1) .eq. names(1) 
     &              .and.
     &              twh_get_impnames(type,inum,2) .eq. names(2)
     &              .and.
     &              twh_get_impnames(type,inum,3) .eq. 'wild'
     &              .and.
     &              twh_get_impnames(type,inum,4).eq. 'wild' 
     &              )then
                  param = type
               endif
            enddo
         endif
      enddo
c     --- check for double wildcard matches in places 2 and 4
      type = 0
      do while ( param .eq. 0 .and. type .lt. TIMPMAX)
         type = type + 1
         if ( twh_get_impff(type) .eq. forcefield ) then
            do inum = 1,twh_get_nimpsame(type)
               if ( twh_get_impnames(type,inum,1) .eq. names(1) 
     &              .and.
     &              twh_get_impnames(type,inum,2) .eq. 'wild'
     &              .and.
     &              twh_get_impnames(type,inum,3) .eq. names(3)
     &              .and.
     &              twh_get_impnames(type,inum,4).eq. 'wild' 
     &              )then
                  param = type
               endif
            enddo
         endif
      enddo
c     --- check for double wildcard matches in places 2 and 3
c     --- combined with a backwards search
      type = 0
      do while ( param .eq. 0 .and. type .lt. TIMPMAX)
         type = type + 1
         if ( twh_get_impff(type) .eq. forcefield ) then
            do inum = 1,twh_get_nimpsame(type)
               if ( twh_get_impnames(type,inum,4) .eq. names(1) 
     &              .and.
     &              twh_get_impnames(type,inum,3) .eq. 'wild'
     &              .and.
     &              twh_get_impnames(type,inum,2) .eq. 'wild'
     &              .and.
     &              twh_get_impnames(type,inum,1).eq.names(4)
     &              )then
                  param = type
               endif
            enddo
         endif
      enddo
c     --- check for triple wildcard matches in all but spot 2
      type = 0
      do while ( param .eq. 0 .and. type .lt. TIMPMAX)
         type = type + 1
         if ( twh_get_impff(type) .eq. forcefield ) then
            do inum = 1,twh_get_nimpsame(type)
               if ( twh_get_impnames(type,inum,1) .eq. 'wild'
     &              .and.
     &              twh_get_impnames(type,inum,2).eq.names(2)
     &              .and.
     &              twh_get_impnames(type,inum,3) .eq. 'wild'
     &              .and.
     &              twh_get_impnames(type,inum,4).eq. 'wild'
     &              )then
                  param = type
               endif
            enddo
         endif
      enddo
c     --- check for triple wildcard matches in all but spot 1
      type = 0
      do while ( param .eq. 0 .and. type .lt. TIMPMAX)
         type = type + 1
         if ( twh_get_impff(type) .eq. forcefield ) then
            do inum = 1,twh_get_nimpsame(type)
               if ( twh_get_impnames(type,inum,1) .eq. names(1)
     &              .and.
     &              twh_get_impnames(type,inum,2).eq. 'wild'
     &              .and.
     &              twh_get_impnames(type,inum,3) .eq. 'wild'
     &              .and.
     &              twh_get_impnames(type,inum,4).eq. 'wild'
     &              )then
                  param = type
               endif
            enddo
         endif
      enddo

      return
      end

      subroutine twh_improper_match_mmff(names,param,forcefield)
c     ******************************************************************
c     * computes the improper types using the MMFF checkdown procedure *
c     *                                                                *
c     * originally written 12-20-2006 by M.G. Martin                   *
c     * last modified 12-20-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IMPFF
#define FUNCTION_GET_IMPNAMES
#define FUNCTION_GET_NIMPSAME
#include "functions.h"
c     --- variables passed to/from the subroutine
      character*(FFNAMELEN) names
      dimension names(4)
      integer param
      character*(FFNAMELEN) forcefield
c     --- local variables
      integer index,istart,iend,count,type,inum
c     --- character string scalars
      character*2 wildcard
      character*(FFNAMELEN) testname
c     --- character string arrays
      character*2 key,testkey
      dimension key(4),testkey(4)

c     --- names have already been put in the proper order for MMFF 
c     --- so we just need to check for matches, but using the checkdown
c     --- procedure

      do index = 1,5
         istart = 2*index-1
         iend = 2*index
         do count = 1,4
            key(count) = names(count)(istart:iend)
         enddo

         do type = 1,TIMPMAX
c           --- see if the forcefield name matches
            if ( twh_get_impff(type) .eq. forcefield ) then
c              --- get the tornames
               do inum = 1,twh_get_nimpsame(type)
                  do count = 1,4
                     testname = twh_get_impnames(type,inum,count)
                     testkey(count) = testname(1:2)
                  enddo
c                 --- check for a forward match
                  if ( key(1) .eq. testkey(1) .and.
     &                 key(2) .eq. testkey(2) .and.
     &                 key(3) .eq. testkey(3) .and.
     &                 key(4) .eq. testkey(4) ) then
c                    --- forward match
                     param = type
                     return
                  endif
               enddo
            endif
         enddo
      enddo

c     --- if that did not work then check for triple wildcard matches
c     --- on all except for jatom, just using the first key
      index = 1
      istart = 2*index-1
      iend = 2*index
      do count = 1,4
         key(count) = names(count)(istart:iend)
      enddo

      wildcard = ' 0'
      do type = 1,TIMPMAX
c        --- see if the forcefield name matches
         if ( twh_get_impff(type) .eq. forcefield ) then
c           --- get the tornames
            do inum = 1,twh_get_nimpsame(type)
               do count = 1,4
                  testname = twh_get_impnames(type,inum,count)
                  testkey(count) = testname(1:2)
               enddo
c              --- check for a triple-wild match
               if ( wildcard .eq. testkey(1) .and.
     &              key(2) .eq. testkey(2) .and.
     &              wildcard .eq. testkey(3) .and.
     &              wildcard .eq. testkey(4) ) then
c                 --- match
                  param = type
                  return
               endif
            enddo
         endif
      enddo

      return
      end
