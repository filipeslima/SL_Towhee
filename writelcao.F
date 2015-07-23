#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2008 Marcus G. Martin                       *
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
      subroutine twh_writelcao(lfailure,ibox,iflag,nchainmove,chaininfo
     &     ,natoms,gsenergy)
c     ******************************************************************
c     * creates the lcao.in file and resets the lcao.stat file in      *
c     * preparation for computing the total system energy using the    *
c     * quantum code lcao (quest)                                      *
c     *                                                                *
c     * itrial: if negative then we are computing energy for full      *
c     *         molecules (engtotal, or typical engmolec)              *
c     *         if 0 then we are computing partial atoms from engmolec *
c     *         if positive then we have partial atoms from engatom    *
c     *                                                                *
c     * originally written 11-05-2004 by M.G. Martin                   *
c     * last modified 10-13-2008 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_ELEMENTNAME
#include "functions.h"
#define DEBUG_WRITELCAO 0
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ibox,iflag,natoms,nchainmove
      integer chaininfo
      dimension chaininfo(nchainmove,3)
      double precision gsenergy
c     --- local variables
      logical lvalid
      character*20 twentychar
      character*30 filename
      integer imolty,ichain,iunit,idim,jdim,itype,iatom,ilist,testchain
      integer ichainmove,testcrd,nmolty,nunit,ncmt,lcao_atomtypes
c     --- double precision scalars
      double precision xtemp,ytemp,ztemp,dvalue,temperature
      double precision lcao_gridmultiplier,lcao_kgridproduct
c     --- double precision arrays
      double precision vector
      dimension vector(3)
      double precision hmtemp
      dimension hmtemp(3)

#if DEBUG_WRITELCAO
      write(6,*) 'starting WRITELCAO'
#endif

c     --- initialize the ground state energy
      gsenergy = 0.0d0

c     --- clean out the lcao.stat file
      filename = 'lcao.stat'
      call twh_openfile(lfailure,71,filename,'DEFAULT','UNKNOWN'
     &     ,'DEFAULT')
      if ( lfailure ) return
      write(71,'(a10)') 'END STATUS'
      rewind(71)
      close(71)

c     --- output lcao.in
      filename = 'lcao.in'
      call twh_openfile(lfailure,72,filename,'DEFAULT','UNKNOWN'
     &     ,'DEFAULT')
      if ( lfailure ) return

c     --- output the header information
      write(72,'(a8)') 'do setup'
      write(72,'(a8)') 'do iters'
      write(72,'(a8)') 'no force'
      write(72,'(a8)') 'no relax'
      write(72,'(a11)') 'setup data:'
      write(72,'(a10)') 'functional'
      call twh_lcao_functional(GLB_GET,twentychar)
      write(72,'(a20)') twentychar
      write(72,'(a9)') 'dimension'
      write(72,'(a1)') '3'
      write(72,'(a5)') 'scale'
      write(72,'(a9)') '1.8897265'
      write(72,'(a26)') 'primitive lattice vectors:'
      do idim = 1,3
         do jdim = 1,3
            call twh_hmatrix(GLB_GET,ibox,idim,jdim,hmtemp(jdim))
         enddo
         write(72,*) (hmtemp(jdim),jdim=1,3)
         vector(idim) = 0.0d0
         do jdim = 1,3
            call twh_hmatrix(GLB_GET,ibox,idim,jdim,dvalue)
            vector(idim) = vector(idim) + dvalue**2
         enddo
         vector(idim) = dsqrt(vector(idim))
      enddo
      write(72,'(a15)') 'grid dimensions'
      call twh_lcao_gridmultiplier(GLB_GET,lcao_gridmultiplier)
      write(72,'(3(i10,1x))') 
     &     (idint(vector(idim)*lcao_gridmultiplier)+1,idim=1,3)

      write(72,'(a10)') 'atom types'
      call twh_lcao_atomtypes(GLB_GET,lcao_atomtypes)
      write(72,*) lcao_atomtypes
      write(72,'(a9)') 'atom file'
      do itype = 1,lcao_atomtypes
         call twh_lcao_filename(GLB_GET,itype,twentychar)
         write(72,'(a20)') twentychar
      enddo

c     --- figure out the number of atoms
      natoms = 0
      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
c        --- count all of the molecules that are not ichain
         call twh_ncmt(GLB_GET,ibox,imolty,ncmt)
         do ilist = 1,ncmt
            call twh_chainlist(GLB_GET,ilist,ibox,imolty,testchain)
            lvalid = .true.
            ichainmove = 0
            do while ( ichainmove .lt. nchainmove .and. lvalid )
               ichainmove = ichainmove + 1
               if ( testchain .eq. chaininfo(ichainmove,3) ) then
                  lvalid = .false.
               endif
            enddo
            if ( lvalid ) then
               call twh_nunit(GLB_GET,imolty,nunit)
               natoms = natoms + nunit
            endif
         enddo
      enddo
#if DEBUG_WRITELCAO
      write(6,*) 'natoms counting those that are not special',natoms
#endif
c     --- account for the special molecules (if appropriate)
      do ichainmove = 1,nchainmove
         testcrd = chaininfo(ichainmove,iflag)
         testchain = chaininfo(ichainmove,3)
         if ( testcrd .ne. CRD_NONE ) then
            call twh_moltyp(GLB_GET,testchain,imolty)
c           --- count these atoms
            call twh_nunit(GLB_GET,imolty,nunit)
            natoms = natoms + nunit
         endif
      enddo
      write(72,'(a)') 'number of atoms'
      write(72,*) natoms
#if DEBUG_WRITELCAO
      write(6,*) 'total natoms:',natoms
#endif
c     --- output the coordinates
      write(72,'(a27)') 'atom, type, position vector'
      iatom = 0
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         call twh_ncmt(GLB_GET,ibox,imolty,ncmt)
         do ilist = 1,ncmt
c           --- assign ichain from the chainlist
            call twh_chainlist(GLB_GET,ilist,ibox,imolty,ichain)
            lvalid = .true.
            ichainmove = 0
            do while ( ichainmove .lt. nchainmove .and. lvalid )
               ichainmove = ichainmove + 1
               if ( ichain .eq. chaininfo(ichainmove,3) ) then
                  lvalid = .false.
               endif
            enddo
            
            if ( lvalid ) then
               do iunit = 1,nunit
                  call twh_ntype(GLB_GET,imolty,iunit,itype)
c                 --- retrive from coordinate system
                  call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &                 ,xtemp,ytemp,ztemp)
                  iatom = iatom + 1
                  write(72,111) iatom,twh_get_elementname(itype)
     &                 ,xtemp,ytemp,ztemp
c                 --- increment the ground state energy
                  call twh_quantum_gsenergy(GLB_GET,imolty,dvalue)
                  gsenergy = gsenergy + dvalue
               enddo
            endif
         enddo
      enddo

c     --- now write out the atom positions for the special chains
      do ichainmove = 1,nchainmove
         ichain = chaininfo(ichainmove,3)
         testcrd = chaininfo(ichainmove,iflag)
         if ( testcrd .ne. CRD_NONE ) then
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
c              --- retrive from coordinate system
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
               iatom = iatom + 1
               write(72,111) iatom,twh_get_elementname(itype)
     &              ,xtemp,ytemp,ztemp
c              --- increment the ground state energy
               call twh_quantum_gsenergy(GLB_GET,imolty,dvalue)
               gsenergy = gsenergy + dvalue
            enddo
         endif
      enddo
 111  format(i10,1x,a2,1x,3(f15.8,1x))
c     --- output the kgrid
      write(72,'(a5)') 'kgrid'
      call twh_lcao_kgridproduct(GLB_GET,lcao_kgridproduct)
      write(72,'(3(i10,1x))') (idint(lcao_kgridproduct/vector(idim))+1
     &     ,idim=1,3)

      write(72,'(a20)') 'end setup phase data'
      write(72,'(a20)') 'run phase input data'
      write(72,'(a11)') 'temperature'
c     --- convert temperature to rydberg
      call twh_temperature(GLB_GET,temperature)
      write(72,*) temperature/(CON_RYDBERGTOEV*CON_EVTOK)
      
      rewind(72)
      close(72)

c     --- safety check
      if ( iatom .ne. natoms ) then
         write(6,*) 'WRITELCAO: mismatch of estimated atoms and '
     &        ,'observed atoms'
         write(6,*) 'iatom:',iatom
         write(6,*) 'natoms:',natoms
         stop
      endif

#if DEBUG_WRITELCAO
      write(6,*) 'done with writelcao'
#endif
      return
      end

