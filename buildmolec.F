#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2011 Marcus G. Martin                       *
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
      subroutine twh_buildmolec(lfailure,imolty,lpdbnames,input_style
     &     ,forcefield,imp_assign,match_style)
c     ******************************************************************
c     * builds the towhee_input file for a molecule given only the     *
c     * forcefield the atom names, and the set of vibrations           *
c     * imolty: the molecule type number                               *
c     *                                                                *
c     * originally written 03-07-2001 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_READ_INTEGER
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"
c     --- debug toggle
#define DEBUG_BUILDMOLEC 0
c     ---- variables passed to/from the subroutine
      logical lfailure,lpdbnames
      integer imolty
      character*(*) input_style
      character*(*) forcefield
      character*(*) imp_assign
      character*(*) match_style
c     --- local variables
c     --- character string scalars
      character*(FFNAMELEN) order
      character*3 aminoshort
      character*4 pdbname
      character*15 idname
      character*30 charge_assignment
c     --- character string arrays
      character*10 atomname
      dimension atomname(NUMAX)
c     --- logical scalars
      logical lautotorsion
c     --- integer scalars
      integer iatom,totatoms,ivib,itemp,imp,numvib,imptype,aminonum
      integer implist2,implist3,implist4,inimp,spectorcount,nunit
c     --- integer arrays
      integer viblist
      dimension viblist(NNBOND)
c     --- double precision scalars
      double precision qtemp
      
#if DEBUG_BUILDMOLED
      write(6,*) 'start BUILDMOLEC'
#endif

      idname = 'buildmolec'
      write(6,'(a45,1x,i5)') 
     &     'Building the input file for molecule type:',imolty

      call twh_charge_assignment(GLB_GET,imolty,charge_assignment)
c     --- read in the atom sequence
      call twh_nunit(GLB_GET,imolty,nunit)
      do iatom = 1,nunit
         call twh_read_label_formatted(lfailure,4,'unit',.false.
     &        ,idname, '(a4)')
         if ( lfailure ) return
         if ( charge_assignment .eq. 'manual' ) then
            read(4,*) itemp,atomname(iatom),qtemp
            call twh_qqatom(GLB_SET,imolty,iatom,qtemp)
            write(6,'(a5,i5,1x,a5,a10,1x,a7,f10.5)')
     &           'unit:',itemp,'name:',atomname(iatom),'charge:'
     &           ,qtemp
         else
            read(4,*) itemp,atomname(iatom)
            write(6,'(a5,i5,1x,a5,a10)')
     &           'unit:',itemp,'name:',atomname(iatom)
         endif
         if ( itemp .ne. iatom ) then
            write(6,*) 'iunit and itemp do not match in towhee_input'
            write(6,*) 'imolty,iunit,itemp',imolty,iatom,itemp
            lfailure = .true.
            return
         endif

         if ( lpdbnames ) then
c                 --- also read in the pdb name
            call twh_read_pdbname(lfailure,4
     &           ,'pdbname aminonum aminoshort',
     &         pdbname, aminonum,aminoshort)
            if ( lfailure ) return
            call twh_pdbname(GLB_SET,imolty,iatom,pdbname)
            call twh_aminonum(GLB_SET,imolty,iatom,aminonum)
            call twh_aminoshort(GLB_SET,imolty,iatom,aminoshort)
            write(6,37) 'pdbname:',pdbname
     &           ,'aminonum:',aminonum,'aminoshort:',aminoshort
37               format(a8,1x,a4,1x,a9,1x,i10,1x,a11,1x,a3)
         endif

         numvib = twh_read_labeled_integer(lfailure,4,
     &      'vibration', .false., idname)
         if ( lfailure ) return
         call twh_invib(GLB_SET,imolty,iatom,numvib)
c        --- safety check
         if ( numvib .gt. NNBOND ) then
            write(6,*) 'invib exceeds NNBOND'
            write(6,*) 'imolty,iatom,invib',imolty,iatom,numvib
            lfailure = .true.
            return
         endif
         if ( numvib .ne. 0 ) then
            if ( input_style .eq. 'basic connectivity map' ) then
c              --- all bond partners on one line without vibration order
               read(4,*) (viblist(ivib),ivib=1,numvib)
               do ivib = 1,numvib
                  call twh_ijvib(GLB_SET,imolty,iatom,ivib
     &                 ,viblist(ivib))
               enddo
            elseif ( input_style .eq. 'advanced connectivity map' ) then
c              --- bond partners and vibration orders on each line
               do ivib = 1,numvib
                  read(4,*) viblist(ivib),order
                  call twh_ijvib(GLB_SET,imolty,iatom,ivib
     &                 ,viblist(ivib))
                  call twh_vibration_order(GLB_SET,iatom,ivib,order)
               enddo
            endif
               
         endif

         if ( imp_assign .eq. 'manual' ) then
c           --- read in improper information
            call twh_read_label_formatted(lfailure,4,'improper',.false.
     &           ,idname, '(a8)')
            if ( lfailure ) return
            inimp = twh_read_integer(lfailure,4, 'improper',.false.)
            if ( lfailure ) return
            call twh_inimprop(GLB_SET,imolty,iatom,inimp)
c           --- safety check
            if ( inimp .gt. MAXIMPROP ) then
               write(6,*) 'imolty,iatom,inimp',imolty,iatom,inimp
               write(6,*) 'too many improper torsions'
               lfailure = .true.
               return
            endif
            do imp = 1,inimp
               read(4,*) implist2,implist3,implist4,imptype
               call twh_itimprop(GLB_SET,imolty,iatom,imp,imptype)
               call twh_ijimprop2(GLB_SET,imolty,iatom,imp,implist2)
               call twh_ijimprop3(GLB_SET,imolty,iatom,imp,implist3)
               call twh_ijimprop4(GLB_SET,imolty,iatom,imp,implist4)
            enddo
         endif
      enddo

      spectorcount = 0
      if ( match_style .eq. 'stereoselective' ) then
c        --- get the special torsion count
         spectorcount = twh_read_labeled_integer(lfailure,4,
     &        'special_torsion_count', .true., idname)
         if ( lfailure ) return
      endif

      call twh_nunit(GLB_GET,imolty,totatoms)
      lautotorsion = .true.
      call twh_assemble(lfailure,imolty,totatoms,spectorcount,forcefield
     &     ,atomname,match_style,imp_assign,lautotorsion)
      if ( lfailure ) return
c     --- check to see if we need to assign the charges
      if ( charge_assignment .ne. 'manual' ) then
         call twh_qqassign(lfailure,imolty,forcefield)
         if ( lfailure ) return
      endif

#if DEBUG_BUILDMOLEC
      write(6,*) 'finished BUILDMOLEC'
#endif

      return
      end
