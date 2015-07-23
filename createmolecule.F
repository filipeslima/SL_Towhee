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
      subroutine twh_createmolecule(lfailure,natoms,nbonds,nangles
     &     ,ndihedrals,nimpropers,nangleangle,ndihedtypes)
c     ******************************************************************
c     * reads in the lammps connectivity information that has already  *
c     * been parsed and output to several files.  Uses this information*
c     * to construct all of the molecule connectivity templates needed *
c     * for a Towhee run                                               *
c     *                                                                *
c     * originally written 03-18-2002 by M.G. Martin                   *
c     * last modified 11-21-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJIMPROP2
#define FUNCTION_GET_IJIMPROP3
#define FUNCTION_GET_IJIMPROP4
#define FUNCTION_GET_IJTOR2
#define FUNCTION_GET_IJTOR3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_INTOR
#define FUNCTION_GET_INIMPROP
#define FUNCTION_GET_ITIMPROP
#define FUNCTION_GET_ITTOR
#include "functions.h"

c     --- variables passed to/from the subroutine
      logical lfailure
      integer natoms,nbonds,nangles,ndihedrals,nimpropers,nangleangle
      integer ndihedtypes
c     --- local variables
      logical lopen,lmatch,lfound,lsame
      dimension lopen(NTMAX)
      character*2 numbers
      character*30 filename
      integer iatom,atomtype,mybonds,bondtype,bondpart,ibond,base
      integer inix,iniy,iniz
      integer imp2,imp3,imp4
      dimension bondtype(NNBOND),bondpart(NNBOND)
      double precision charge,xtmp,ytmp,ztmp,checkx,checky,checkz,s
      integer iangle,angtype,angpart,myangles,iunit,ibox,nmtype
      dimension angtype(MAXBEND)
      dimension angpart(MAXBEND,2)
      integer mydiheds,dihtype,dihpart,itor,iimp,ilast
      dimension dihtype(MAXTOR)
      dimension dihpart(MAXTOR,3)
      integer myimprops,impropertype,imppart
      dimension impropertype(MAXIMPROP)
      dimension imppart(MAXIMPROP,3)
      integer ninvert,listinvert,current,lmpnum,itmp,imolty,isearch
      dimension listinvert(NUMAX)
      integer listatom,type,index,zero,one,two,newmolty
      integer myangang,myaatype,aapart
      dimension myaatype(MAXAA)
      dimension aapart(MAXAA,3)
      integer isplit,tortype,imptype,imatch,inew,iother,totmolecs,itemp
      integer iorder,nmolty

c     Floating point comparisons vary between architectures.  
c     See http://www.lahey.com/float.htm for discussion
      double precision epsilon
      parameter (epsilon = 1d-8)

      write(6,*) 'Creating molecule data structures'

c     --- set xcordtmp to -1.0 for all atoms to indicate that none of 
c     --- them have been analyzed for data structures
      do iatom = 1,natoms
         index = 3*(iatom-1) + 1
         call twh_coordinates(GLB_FIX,CRD_TEMP,0,index
     &        ,-1.0d0,0.0d0,0.0d0)
      enddo
c     --- set lopen to false to indicate we have not opened any of 
c     --- the lammps coordinate files
      filename = 'lammps_coords_xx'
      base = 80
      do imolty = 1,NTMAX
         lopen(imolty) = .false.
      enddo
c     --- initialize nmolectyp
      call twh_nmolectyp(GLB_INIT,NTMAX,0)
      call twh_initmol(GLB_INIT,MAXBOX,NTMAX,0)
c     --- only have one box from lammps, set ibox to 1
      call twh_numboxes(GLB_SET,1)
      ibox = 1
c     --- initialize listatom and nmolty
      listatom = 1
      nmolty = 0
      do while ( listatom .le. natoms )
         write(6,*) 'listatom',listatom
         index = 3*(listatom-1) + 1
         call twh_coordinates(GLB_CHK,CRD_TEMP,0,index
     &        ,checkx,checky,checkz)
         if ( checkx .lt. 0.0d0 ) then
c           --- start with this atom to build a molecule data structure
            nmolty = nmolty + 1
            if ( nmolty .gt. NTMAX ) then
c              --- problem
               write(6,*) 'maximum molecule types (NTMAX) exceeded '
               write(6,*) 'increase NTMAX'
               lfailure = .true.
               return
            endif
c           --- initialize the current working atom
            current = 1
c           --- initialize the inverse list
            ninvert = 1
            listinvert(1) = listatom
            do while ( current .le. ninvert )
c              --- get the lammps atom number for the current atom
               lmpnum = listinvert(current)
c              --- set the atomtype and charge for this atom
               call twh_getatominfo(lfailure,lmpnum,atomtype,charge
     &              ,xtmp,ytmp,ztmp)
               if ( lfailure ) return
               call twh_qqatom(GLB_SET,nmolty,current,charge)
               if ( charge .ne. 0.0d0) then
                  call twh_lcoulomb(GLB_SET,.true.)
                  call twh_lewald(GLB_SET,.true.)
               endif
               call twh_ntype(GLB_SET,nmolty,current,atomtype)

c              --- find all of the atoms bonded to this atom
               call twh_getbondinfo(lfailure,nbonds,lmpnum,mybonds
     &              ,bondtype,bondpart)
               if ( lfailure ) return
c              --- assign the number of bonds
               call twh_invib(GLB_SET,nmolty,current,mybonds)
c              --- sort through the bonds and see if they already have a
c              --- number, if not then assign them one and add to the 
c              --- list
               do ibond = 1,mybonds
c                 --- assign type
                  call twh_itvib(GLB_SET,nmolty,current,ibond
     &                 ,bondtype(ibond))
c                 --- assign bond partner
                  isearch = bondpart(ibond)
                  call twh_findatom(lfailure,isearch,itmp,ninvert
     &                 ,listinvert)
                  if ( lfailure ) return
                  call twh_ijvib(GLB_SET,nmolty,current,ibond,itmp)
               enddo

c              --- find all of the angles that have this atom at one end
               call twh_getangleinfo(lfailure,nangles,lmpnum,myangles
     &              ,angtype,angpart)
               if ( lfailure ) return
c              --- assign the number of angles
               call twh_inben(GLB_SET,nmolty,current,myangles)
c              --- sort through the angles, assign types, add to lists
               do iangle = 1,myangles
c                 --- set the type
                  if ( angtype(iangle) .lt. 0 ) then
                     angtype(iangle) = -angtype(iangle)
                     iorder = ORD_REVERSE
                  else
                     iorder = ORD_FORWARD
                  endif
                  write(6,*) 'iangle,angtype:',iangle,angtype(iangle)
                  call twh_itben(GLB_SET,nmolty,current,iangle
     &                 ,angtype(iangle))
                  call twh_orderben(GLB_SET,nmolty,current,iangle
     &                 ,iorder)
c                 --- assign the middle atom
                  isearch = angpart(iangle,1)
                  call twh_findatom(lfailure,isearch,itmp,ninvert
     &                 ,listinvert)
                  if ( lfailure ) return
                  call twh_ijben2(GLB_SET,nmolty,current,iangle,itmp)
c                 --- assign other end atom
                  isearch = angpart(iangle,2)
                  call twh_findatom(lfailure,isearch,itmp,ninvert
     &                 ,listinvert)
                  if ( lfailure ) return
                  call twh_ijben3(GLB_SET,nmolty,current,iangle,itmp)
               enddo

c              --- find all of the torsions that have this at one end
               call twh_getdihedralinfo(lfailure,ndihedrals,lmpnum
     &              ,mydiheds,dihtype,dihpart)
               if ( lfailure ) return
c              --- assign the number of dihedrals
               call twh_intor(GLB_SET,nmolty,current,mydiheds)
c              --- sort through dihedrals, assign types, add to lists
               do itor = 1,mydiheds
c                 --- set the type
                  call twh_ittor(GLB_SET,nmolty,current,itor
     &                 ,dihtype(itor))
c                 --- assign torsion atom 2
                  isearch = dihpart(itor,1)
                  call twh_findatom(lfailure,isearch,itmp,ninvert
     &                 ,listinvert)
                  if ( lfailure ) return
                  call twh_ijtor2(GLB_SET,nmolty,current,itor,itmp)
c                 --- assign torsion atom 3
                  isearch = dihpart(itor,2)
                  call twh_findatom(lfailure,isearch,itmp,ninvert
     &                 ,listinvert)
                  if ( lfailure ) return
                  call twh_ijtor3(GLB_SET,nmolty,current,itor,itmp)
c                 --- assign torsion atom 4
                  isearch = dihpart(itor,3)
                  call twh_findatom(lfailure,isearch,itmp,ninvert
     &                 ,listinvert)
                  if ( lfailure ) return
                  call twh_ijtor4(GLB_SET,nmolty,current,itor,itmp)
               enddo

c              --- find all of the impropers that have this atom at 
c              --- the first position
               call twh_getimproperinfo(lfailure,nimpropers,lmpnum
     &              ,myimprops,impropertype,imppart)
               if ( lfailure ) return
c              --- assign the number of impropers
               call twh_inimprop(GLB_SET,nmolty,current,myimprops)
c              --- sort through impropers, assign types, add to list
               do iimp = 1,myimprops
c                 --- set the type
                  call twh_itimprop(GLB_SET,nmolty,current,iimp
     &                 ,impropertype(iimp))
c                 --- assign improper atom 2
                  isearch = imppart(iimp,1)
                  call twh_findatom(lfailure,isearch,itmp,ninvert
     &                 ,listinvert)
                  if ( lfailure ) return
                  call twh_ijimprop2(GLB_SET,nmolty,current,iimp,itmp)
c                 --- assign improper atom 3
                  isearch = imppart(iimp,2)
                  call twh_findatom(lfailure,isearch,itmp,ninvert
     &                 ,listinvert)
                  if ( lfailure ) return
                  call twh_ijimprop3(GLB_SET,nmolty,current,iimp,itmp)
c                 --- assign improper atom 4
                  isearch = imppart(iimp,3)
                  call twh_findatom(lfailure,isearch,itmp,ninvert
     &                 ,listinvert)
                  if ( lfailure ) return
                  call twh_ijimprop4(GLB_SET,nmolty,current,iimp,itmp)
               enddo

c              --- find all of the angle-angle terms that have this 
c              --- atom at the second position
               call twh_getangleangleinfo(lfailure,nangleangle,lmpnum
     &              ,myangang,myaatype,aapart)
               if ( lfailure ) return
c              --- assign the number of angle-angle terms
               call twh_inaa(GLB_SET,nmolty,current,3*myangang)
c              --- sort through impropers, assign types, add to list
               do itor = 1,myangang
c                 --- need to divide the LAMMPS angle-angle terms into
c                 --- three different angle-angle terms for use in 
c                 --- towhee
                  do isplit = 1,3
                     if ( isplit .eq. 1 ) then
                        zero = 2
                        one = 1
                        two = 3
                     elseif ( isplit .eq. 2 ) then
                        zero = 1
                        one = 2
                        two = 3
                     elseif ( isplit .eq. 3 ) then
                        zero = 3
                        one = 1
                        two = 2
                     else
                        write(6,*) 'CREATEMOLECULE: unknown isplit'
     &                       ,isplit
                        lfailure = .true.
                        return
                     endif
                     index = 3*(itor-1) + isplit
                     type = 3*(myaatype(itor)-1) + isplit
c                    --- set the type
                     call twh_itaa(GLB_SET,nmolty,current,index,type)
c                    --- assign angle-angle atom 0
                     isearch = aapart(itor,zero)
                     call twh_findatom(lfailure,isearch,itmp,ninvert
     &                    ,listinvert)
                     if ( lfailure ) return
                     call twh_ijaa0(GLB_SET,nmolty,current,index,itmp)
c                    --- assign angle-angle atom 1
                     isearch = aapart(itor,one)
                     call twh_findatom(lfailure,isearch,itmp,ninvert
     &                    ,listinvert)
                     if ( lfailure ) return
                     call twh_ijaa1(GLB_SET,nmolty,current,index,itmp)
c                    --- assign angle-angle atom 2
                     isearch = aapart(itor,two)
                     call twh_findatom(lfailure,isearch,itmp,ninvert
     &                    ,listinvert)
                     if ( lfailure ) return
                     call twh_ijaa2(GLB_SET,nmolty,current,index,itmp)
                  enddo
               enddo
c              --- set current to the next atom
               current = current + 1
            enddo
c           --- set the number of units
            call twh_nunit(GLB_SET,nmolty,ninvert)
            call twh_nmaxcbmc(GLB_SET,nmolty,ninvert)

c           --- check to see if any of those improper torsions are 
c           --- actually across regular torsion atoms.  This occurs
c           --- with the LAMMPS implementation of the Charmm force 
c           --- field
            do current = 1,ninvert
               iimp = 1
               imatch = 0
               do while ( iimp .le. twh_get_inimprop(nmolty,current) )
                  lmatch = .false.
                  do itor = 1,twh_get_intor(nmolty,current)
                     if ( twh_get_ijimprop2(nmolty,current,iimp) .eq. 
     &                    twh_get_ijtor2(nmolty,current,itor)
     &                    .and.  
     &                    twh_get_ijimprop3(nmolty,current,iimp) .eq. 
     &                    twh_get_ijtor3(nmolty,current,itor)
     &                    .and.  
     &                    twh_get_ijimprop4(nmolty,current,iimp) .eq. 
     &                    twh_get_ijtor4(nmolty,current,itor) ) then
c                       --- this improper is actually an existing 
c                       --- regular torsion
                        if ( lmatch ) then
c                          --- strangely matches two torsions
                           write(6,*) 'CREATEMOLECUE: 2 torsion matches'
                           lfailure = .true.
                           return
                        endif
c                       --- flag this for adjustment
                        lmatch = .true.
                        imatch = itor
                     endif
                  enddo

                  if ( lmatch ) then
c                    --- adjust the types and numbers
                     tortype = twh_get_ittor(nmolty,current,imatch)
                     imptype = twh_get_itimprop(nmolty,current,iimp)
c                    --- combine to make a new torsion type, see if 
c                    --- this particular combination has been 
c                    --- added, if not they create it
                     
                     call twh_combineregandimp(lfailure,inew,tortype
     &                    ,imptype,ndihedtypes)
                     if ( lfailure ) return
c                    --- reassign torsion type
                     call twh_ittor(GLB_SET,nmolty,current,imatch,inew)
c                    --- change the torsion type of the symmetric 
c                    --- partner
                     iother = twh_get_ijtor4(nmolty,current,imatch)
                     lfound = .false.
                     do itor = 1,twh_get_intor(nmolty,iother)
c                       --- see if this is symmetric with current
                        if ( twh_get_ijtor2(nmolty,iother,itor) .eq. 
     &                       twh_get_ijtor3(nmolty,current,imatch) .and.
     &                       twh_get_ijtor3(nmolty,iother,itor) .eq. 
     &                       twh_get_ijtor2(nmolty,current,imatch) .and.
     &                       twh_get_ijtor4(nmolty,iother,itor) .eq. 
     &                       current ) then
c                          --- found the symmetric torsion
                           if ( lfound ) then
c                             --- should only find one sym torsion
                              write(6,*) 'CREATEMOLECULE: two symmtric'
     &                             ,' torsion found'
                              lfailure = .true.
                              return
                           endif
                           lfound = .true.
c                          --- assign the new torsion type
c                          --- make sure we maintain the sign of the
c                          --- other torsion as that matters for class2
                           if ( twh_get_ittor(nmolty,iother,itor) 
     &                          .gt. 0 ) then
                              call twh_ittor(GLB_SET,nmolty,iother,itor
     &                             ,inew)
                           else
                              call twh_ittor(GLB_SET,nmolty,iother,itor
     &                             ,-inew)
                           endif
                        endif
                     enddo
                     if ( .not. lfound ) then
                        write(6,*) 'CREATEMOLECULE: symmetric torsion '
     &                       ,'not found'
                        lfailure = .true.
                        return
                     endif

c                    --- switch improper information with the last
c                    --- improper
                     ilast = twh_get_inimprop(nmolty,current)
                     itemp = twh_get_itimprop(nmolty,current,ilast) 
                     call twh_itimprop(GLB_SET,nmolty,current,iimp
     &                    ,itemp)
                     imp2 = twh_get_ijimprop2(nmolty,current,ilast)
                     call twh_ijimprop2(GLB_SET,nmolty,current,iimp
     &                    ,imp2)
                     imp3 = twh_get_ijimprop3(nmolty,current,ilast)
                     call twh_ijimprop3(GLB_SET,nmolty,current,iimp
     &                    ,imp3)
                     imp4 = twh_get_ijimprop4(nmolty,current,ilast)
                     call twh_ijimprop4(GLB_SET,nmolty,current,iimp
     &                    ,imp4)
c                    --- decrease inimp, but don't change iimp
                     call twh_inimprop(GLB_DECR,nmolty,current,1)
                  else
c                    --- increment iimp
                     iimp = iimp + 1
                  endif
               enddo
            enddo

c           --- reset xcordtmp for all units in this molecule
            do current = 1,ninvert
               isearch = listinvert(current)
               index = 3*(isearch-1) + 1
             call twh_coordinates(GLB_FIX,CRD_TEMP,0,index
     &              ,1.0d0,0.0d0,0.0d0)
            enddo
c           --- see if this "new" molecule type is actually the same
c           --- as one we have already generated
            newmolty = nmolty
            imolty = 1
            lsame = .false.
            do while ( imolty .lt. nmolty .and. .not. lsame )
               call twh_comparestruc(imolty,nmolty,lsame)
               if ( lsame ) then
                  newmolty = imolty
                  nmolty = nmolty - 1
               endif
               imolty = imolty + 1
            enddo

c           --- increment initmol and nmolectyp
            call twh_initmol(GLB_INCR,ibox,newmolty,1)
            call twh_nmolectyp(GLB_INCR,newmolty,1)
c           --- see if we have opened a file to output the coords for
c           --- this type yet
            if ( .not. lopen(newmolty) ) then
c              --- open a file to place these coordinates
               write(numbers,'(i2.2)') newmolty
               filename(15:16) = numbers
               call twh_openfile(lfailure,base+newmolty,filename
     &              ,'DEFAULT','UNKNOWN','DEFAULT')
               if ( lfailure ) return
c              --- reset lopen
               lopen(newmolty) = .true.
            endif
c           --- output the coordinates
            do iunit = 1,ninvert
c              --- get the lammps unit number for this atom
               isearch = listinvert(iunit)
c              --- get the coordinates
               call twh_getatominfo(lfailure,isearch,atomtype,charge
     &              ,xtmp,ytmp,ztmp)
               if ( lfailure ) return
c              --- output the coordinates
               write(base+newmolty,*) xtmp,ytmp,ztmp
            enddo
         endif
c        --- increment listatom
         listatom = listatom + 1
      enddo
c     --- store the final value of nmolty
      call twh_nmolty(GLB_SET,nmolty)
c     --- close the lammps_coord files
      do imolty = 1,nmolty
         close(base+imolty)
      enddo
c     --- set inix, iniy and iniz based on the total number of molecules
      totmolecs = 0
      do imolty = 1,nmolty
         call twh_nmolectyp(GLB_GET,imolty,nmtype)
         totmolecs = totmolecs + nmtype
      enddo

c     if s = totmolecs**(1/3) is an integer, then box size = s
c     Otherwise, box size = next integer larger than s 
      s = totmolecs**(1.0d0/3.0d0)
      if (abs(s - int(s)) .lt. epsilon) then
          inix = int(s)
          iniy = int(s)
          iniz = int(s)
      else
          inix = int(s) + 1
          iniy = int(s) + 1
          iniz = int(s) + 1
      endif 
      call twh_inixyz(GLB_SET,1,inix,iniy,iniz)

c     --- output the towhee_input header
      call twh_writetowhee(lfailure,.false.,0)
      if ( lfailure ) return

      do imolty = 1,nmolty
c        --- output the data structure
         call twh_writetowhee(lfailure,.false.,imolty)
         if ( lfailure ) return
      enddo

      return
      end

      subroutine twh_getatominfo(lfailure,lam_id,type,charge
     &     ,xtmp,ytmp,ztmp)
c     ******************************************************************
c     * searches through the lammps_atoms file to find the atom with   *
c     * LAMMPS identifier of lam_id and then returns the LAMMPS atom   *
c     * type, charge, and positions                                    *
c     *                                                                *
c     * originally written 03-18-2002 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer lam_id,type
      double precision charge,xtmp,ytmp,ztmp
c     --- local variables
      integer iatom,itmp
      character*30 filename

      filename = 'lammps_atoms'
      call twh_openfile(lfailure,13,filename,'DEFAULT','OLD','DEFAULT')
      if ( lfailure ) return
      do iatom = 1,lam_id
         read(13,*) itmp,type,charge,xtmp,ytmp,ztmp
      enddo
c     --- the data in the variables is now that for lam_id so close file
c     --- and return
      rewind(13)
      close(13)
      return
      end

      subroutine twh_getbondinfo(lfailure,nbonds,lam_id,nmybonds,types
     &     ,partners)
c     ******************************************************************
c     * searches through the lammps_bonds file to find the bond which  *
c     * contain the atom lam_id.  It then makes a list of those bonds  *
c     * and returns.  There is a safety check to make sure we don't    *
c     * have more than NNBOND bonds from any atom                      *
c     *                                                                *
c     * originally written 03-18-2002 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer nbonds,lam_id,nmybonds,types,partners
      dimension types(NNBOND),partners(NNBOND)
c     --- local variables
      logical lnew
      integer ibond,atoma,atomb,btype,other
      character*30 filename

      nmybonds = 0

      filename = 'lammps_bonds'
      call twh_openfile(lfailure,13,filename,'DEFAULT','OLD','DEFAULT')
      if ( lfailure ) return
      do ibond = 1,nbonds
         lnew = .false.
         read(13,*) atoma,atomb,btype
         if ( atoma .eq. lam_id ) then
            lnew = .true.
            other = atomb
         elseif ( atomb .eq. lam_id ) then
            lnew = .true.
            other = atoma
         endif

         if ( lnew ) then
c           --- found a new bond, add it to the data structure
            nmybonds = nmybonds + 1
            if ( nmybonds .gt. NNBOND ) then
               write(6,*) 'too many bonds to atom ',lam_id
               write(6,*) 'you need to increase NNBOND'
               lfailure = .true.
               return
            endif
            types(nmybonds) = btype
            partners(nmybonds) = other
         endif
      enddo
c     --- the data in the variables is now that for lam_id so close file
c     --- and return
      rewind(13)
      close(13)
      return
      end

      subroutine twh_getangleinfo(lfailure,nangles,lam_id,myangles,types
     &     ,partners)
c     ******************************************************************
c     * searches through the lammps_angles file to find the angles that*
c     * contain the atom lam_id at one of the ends of the angle.  Then *
c     * makes a list of those angles and returns.                      *
c     * There is a safety check to make sure we don't have more than   *
c     * MAXBEND angles from any atom                                   *
c     *                                                                *
c     * originally written 03-18-2002 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_ANGLESTYLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer nangles,lam_id,myangles,types,partners
      dimension types(MAXBEND)
      dimension partners(MAXBEND,2)
c     --- local variables
      logical lnew
      integer iangle,atoma,atomb,atomc,atype,othera,otherb,iorder
      character*30 filename

      myangles = 0

      filename = 'lammps_angles'
      call twh_openfile(lfailure,13,filename,'DEFAULT','OLD','DEFAULT')
      if ( lfailure ) return
      do iangle = 1,nangles
         lnew = .false.
         read(13,*) atoma,atomb,atomc,atype
         if ( atoma .eq. lam_id ) then
            lnew = .true.
            iorder = 1
            othera = atomb
            otherb = atomc
         elseif ( atomc .eq. lam_id ) then
            iorder = -1
            lnew = .true.
            othera = atomb
            otherb = atoma
         endif

         if ( lnew ) then
c           --- found a new angle, add it to the data structure
            myangles = myangles + 1
            if ( myangles .gt. MAXBEND ) then
               write(6,*) 'too many angles to atom ',lam_id
               write(6,*) 'you need to increase NNBOND'
               lfailure = .true.
               return
            endif
            if ( twh_get_anglestyle(atype) .eq. ANG_STYLE_COMPASS_EX 
     &           .and. iorder .eq. -1 ) then
c              --- this is an explicitly declared class 2 angle that 
c              --- originated from the opposite atom.  This gets a 
c              --- negative type number
               types(myangles) = -atype
            else
c              --- normal positive type numbering system
               types(myangles) = atype
            endif
            partners(myangles,1) = othera
            partners(myangles,2) = otherb
         endif
      enddo
c     --- the data in the variables is now that for lam_id so close file
c     --- and return
      rewind(13)
      close(13)
      return
      end

      subroutine twh_getdihedralinfo(lfailure,ndihedrals,lam_id,mydiheds
     &     ,types,partners)
c     ******************************************************************
c     * searches through the lammps_dihedrals file to find the         *
c     * dihedrals that contain the atom lam_id at one of the ends      *
c     * Then makes a list of those dihedrals and returns.              *
c     * There is a safety check to make sure we don't have more than   *
c     * MAXTOR dihedrals from any atom                                 *
c     *                                                                *
c     * originally written 03-18-2002 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_TORSTYLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ndihedrals,lam_id,mydiheds,types,partners
      dimension types(MAXTOR)
      dimension partners(MAXTOR,3)
c     --- local variables
      logical lnew
      integer itor,atoma,atomb,atomc,atomd,ttype,othera,otherb,otherc
      integer iorder
      character*30 filename

      mydiheds = 0

      filename = 'lammps_dihedrals'
      call twh_openfile(lfailure,13,filename,'DEFAULT','OLD','DEFAULT')
      if ( lfailure ) return
      do itor = 1,ndihedrals
         lnew = .false.
         read(13,*) atoma,atomb,atomc,atomd,ttype
         if ( atoma .eq. lam_id ) then
            lnew = .true.
            othera = atomb
            otherb = atomc
            otherc = atomd
            iorder = 1
         elseif ( atomd .eq. lam_id ) then
            lnew = .true.
            othera = atomc
            otherb = atomb
            otherc = atoma
            iorder = -1
         endif

         if ( lnew ) then
c           --- found a new dihedral, add it to the data structure
            mydiheds = mydiheds + 1
            if ( mydiheds .gt. MAXTOR ) then
               write(6,*) 'too many dihedrals to atom ',lam_id
               write(6,*) 'you need to increase NNBOND'
               lfailure = .true.
               return
            endif
            if ( twh_get_torstyle(ttype) .eq. 9 
     &           .and. iorder .eq. -1 ) then
c              --- this is the explicitly declared class 2 torsion 
c              --- which is not symmetric so we declare a negative 
c              --- type when going from the other direction
               types(mydiheds) = -ttype
            else
c              --- either regular torsion or right order
               types(mydiheds) = ttype
            endif
            partners(mydiheds,1) = othera
            partners(mydiheds,2) = otherb
            partners(mydiheds,3) = otherc
         endif
      enddo
c     --- the data in the variables is now that for lam_id so close file
c     --- and return
      rewind(13)
      close(13)
      return
      end

      subroutine twh_getimproperinfo(lfailure,nimpropers,lam_id
     &     ,myimprops,types,partners)
c     ******************************************************************
c     * searches through the lammps_impropers file to find the         *
c     * impropers that contain the atom lam_id at the first end        *
c     * Then makes a list of those impropers and returns.              *
c     * There is a safety check to make sure we don't have more than   *
c     * MAXIMPROP impropers from any atom                              *
c     *                                                                *
c     * originally written 03-18-2002 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IMPCOEFF
#define FUNCTION_GET_IMPFORM
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer nimpropers,lam_id,myimprops,types,partners
      dimension types(MAXIMPROP)
      dimension partners(MAXIMPROP,3)
c     --- local variables
      logical lnew
      integer itor,atoma,atomb,atomc,atomd,ttype,othera,otherb,otherc
      character*30 filename
      
      myimprops = 0
c     --- there is no file to open if nimpropers is 0
      if ( nimpropers .eq. 0 ) return

c     --- otherwise open the file and read in the improper information
      filename = 'lammps_impropers'
      call twh_openfile(lfailure,13,filename,'DEFAULT','OLD','DEFAULT')
      if ( lfailure ) return
      do itor = 1,nimpropers
         lnew = .false.
         read(13,*) atoma,atomb,atomc,atomd,ttype
c        --- for maximum difficulty, it appears that the harmonic and 
c        --- class 2 force fields use different defaults for which atom
c        --- is central.  Check the style of this type to figure out 
c        --- which it is.  If the style is zero then this is a null 
c        --- improper which is actually an angle-angle term
         if ( twh_get_impform(ttype) .eq. 2 ) then
c           --- defined from the first atom
            if ( atoma .eq. lam_id ) then
               lnew = .true.
               othera = atomb
               otherb = atomc
               otherc = atomd
            endif
         elseif ( twh_get_impform(ttype) .eq. 3 ) then
c           --- class 2, defined from the second atom
            if ( atomb .eq. lam_id ) then
c              --- don't bother with the "null" improper torsions
               if ( twh_get_impcoeff(ttype,0) .gt. 0.0d0 ) then
                  lnew = .true.
                  othera = atoma
                  otherb = atomc
                  otherc = atomd
               endif
            endif
         endif

         if ( lnew ) then
c           --- found a new improper, add it to the data structure
            myimprops = myimprops + 1
            if ( myimprops .gt. MAXIMPROP ) then
               write(6,*) 'myimprops, maximprop '
     &              ,myimprops,MAXIMPROP
               write(6,*) 'too many impropers to atom ',lam_id
               write(6,*) 'you need to increase NNBOND'
               lfailure = .true.
               return
            endif
            types(myimprops) = ttype
            partners(myimprops,1) = othera
            partners(myimprops,2) = otherb
            partners(myimprops,3) = otherc
         endif
      enddo
c     --- the data in the variables is now that for lam_id so close file
c     --- and return
      rewind(13)
      close(13)
      return
      end

      subroutine twh_getangleangleinfo(lfailure,nangleangle,lam_id
     &     ,myangang,types,partners)
c     ******************************************************************
c     * searches through the lammps_impropers file to find the         *
c     * angle-angle terms (LAMMPS treats them as impropers) that       *
c     * contain the atom lam_id at the second position                 *
c     * Then makes a list of those angle-angles and returns.           *
c     * There is a safety check to make sure we don't have more than   *
c     * MAXAA terms from any atom                                      *
c     *                                                                *
c     * originally written 03-22-2002 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_AACOEFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer nangleangle,lam_id,myangang,types,partners
      dimension types(MAXAA)
      dimension partners(MAXAA,3)
c     --- local variables
      logical lnew
      integer itor,atoma,atomb,atomc,atomd,ttype,othera,otherb,otherc
      integer index
      character*30 filename

      myangang = 0
c     --- there is nothing to do if nangleangle is zero
      if ( nangleangle .eq. 0 ) return

c     --- otherwise open the file and search the impropers
      filename = 'lammps_impropers'
      call twh_openfile(lfailure,13,filename,'DEFAULT','OLD','DEFAULT')
      if ( lfailure ) return
      do itor = 1,nangleangle
         lnew = .false.
         read(13,*) atoma,atomb,atomc,atomd,ttype
         if ( atomb .eq. lam_id ) then
c           --- make sure we have some nonzero prefactors
            index = 3*(ttype-1) 
            if ( (twh_get_aacoeff(index+1,0) .ne. 0.0d0) 
     &           .or. (twh_get_aacoeff(index+2,0) .ne. 0.0d0 ) 
     &           .or. ( twh_get_aacoeff(index+3,0) .ne. 0.0d0) ) then
               lnew = .true.
               othera = atoma
               otherb = atomc
               otherc = atomd
            endif
         endif

         if ( lnew ) then
c           --- found a new angle-angle, add it to the data structure
            myangang = myangang + 1
            if (myangang.gt. MAXAA )then
               write(6,*) 'too many angle-angles to atom ',lam_id
               write(6,*) 'you need to increase NNBOND'
               lfailure = .true.
               return
            endif
            types(myangang) = ttype
            partners(myangang,1) = othera
            partners(myangang,2) = otherb
            partners(myangang,3) = otherc
         endif
      enddo
c     --- the data in the variables is now that for lam_id so close file
c     --- and return
      rewind(13)
      close(13)
      return
      end

      subroutine twh_findatom(lfailure,isearch,itmp,ninvert,listinvert)
c     ******************************************************************
c     * this subroutine searches through the invertlist to see if      *
c     * isearch is on that list.  If so, then it returns the position  *
c     * if not, then it adds it to the list and returns the position   *
c     *                                                                *
c     * originally written 03-19-2002 by M.G. Martin                   *
c     * last modified 05-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer isearch,itmp,ninvert,listinvert
      dimension listinvert(NUMAX)
c     --- local variables
      logical lfound

      lfound = .false.
      itmp = 0
      do while ( itmp .lt. ninvert .and. (.not. lfound) )
         itmp = itmp + 1
         if ( isearch .eq. listinvert(itmp)) lfound = .true.
      enddo
c     --- see if we found the atom
      if ( .not. lfound ) then
c        --- we did not find the atom, add it to the end of
c        --- the list
         ninvert = ninvert + 1
         if ( ninvert .gt. NUMAX ) then
            write(6,*) 'CREATEMOLECULE: ninvert > NUMAX'
            lfailure = .true.
            return
         endif
c        --- store the lammps atom number
         listinvert(ninvert) = isearch
c        --- set itmp to ninvert
         itmp = ninvert
      endif
      return
      end

      subroutine twh_combineregandimp(lfailure,inew,tortype,imptype
     &     ,ndihedtypes)
c     ******************************************************************
c     * this subroutine combines a regular torsion and an improper     *
c     * torsion to make a new regular torsion.  First it checks to see *
c     * if this new combined torsion already exists, otherwise it adds *
c     * a new torsion                                                  *
c     *                                                                *
c     * originally written 03-22-2002 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IMPCOEFF
#define FUNCTION_GET_IMPSTYLE
#define FUNCTION_GET_NTORLOOP
#define FUNCTION_GET_TORCOEFF
#define FUNCTION_GET_TORSTYLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer inew,tortype,imptype,ndihedtypes
c     --- local variables
      logical lmatch
      integer ittype,newstyle,index,itor,tstyle,istyle,itemp
      double precision toler,dtemp
      parameter (toler = 0.001d0)

c     --- determine how to combine this regular and improper
      tstyle = twh_get_torstyle(tortype)
      istyle = twh_get_impstyle(imptype)

      if ( tstyle .eq. 3 .and. istyle .eq. 1 ) then
c        --- we can handle this
         newstyle = 4
      else
         write(6,*) 'CREATEMOLECULE: cannot combine regular torsion '
     &        ,'style ',tstyle,' with improper style ',istyle
         lfailure = .true.
         return
      endif
      
      do ittype = 1,ndihedtypes
         lmatch = .true.
c        --- see if this dihedral has the new style
         if ( twh_get_torstyle(ittype) .eq. newstyle ) then
c           --- assume a match, check coefficients to see if any don't
c           --- match
            if ( newstyle .eq. 4 ) then
c              --- check the regular torsion part
               if ( twh_get_ntorloop(ittype) 
     &              .ne. twh_get_ntorloop(tortype) ) then
                  lmatch = .false.
               else
c                 --- check all the regular torsion coefficients
                  do itor = 1,3*twh_get_ntorloop(tortype)
                     if ( (twh_get_torcoeff(tortype,itor)
     &                    -twh_get_torcoeff(ittype,itor))**2 
     &                    .gt. toler ) lmatch = .false.
                  enddo
c                 --- check the harmonic improper part
                  index = 3*twh_get_ntorloop(tortype)+1
                  do itor = 0,1
                     if ( (twh_get_impcoeff(imptype,itor)
     &                    -twh_get_torcoeff(ittype,index+itor))**2
     &                    .gt. toler ) lmatch = .false.
                  enddo
               endif
            endif
         else
c           --- did not find a match
            lmatch = .false.
         endif
         if ( lmatch ) then
c           --- found a match, return the torsion type
            inew = ittype
            goto 200
         endif
      enddo

c     --- if we got here then we did not find a match, make a new
c     --- torsion type
      ndihedtypes = ndihedtypes + 1
      call twh_torstyle(GLB_SET,ndihedtypes,newstyle)
      if ( newstyle .eq. 4 ) then
         itemp = twh_get_ntorloop(tortype)
         call twh_ntorloop(GLB_SET,ndihedtypes,itemp)
c        --- regular torsion part
         do itor = 1,3*twh_get_ntorloop(tortype)
            dtemp = twh_get_torcoeff(tortype,itor)
            call twh_torcoeff(GLB_SET,ndihedtypes,itor,dtemp)
         enddo
c        --- harmonic part
         index = 3*twh_get_ntorloop(tortype)+1
         do itor = 0,1
            dtemp = twh_get_impcoeff(imptype,itor)
            call twh_torcoeff(GLB_SET,ndihedtypes,index+itor,dtemp)
         enddo
      endif
c     --- set inew 
      inew = ndihedtypes

 200  return
      end
