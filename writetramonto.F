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
      subroutine twh_writetramonto(lfailure,iflag,activebox,nchainmove
     &     ,chaininfo)
c     ******************************************************************
c     * creates a dft_surfaces file and a dft_decoder file for use in  *
c     * the Tramonto classical density functional theory program.      *
c     * this is done only for box one in the simulation                *
c     *                                                                *
c     * originally written 07-11-2001 by M.G. Martin                   *
c     * last modified 08-12-2011 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define DEBUG_WRITETRAMONTO 0
c     --- variables passed to/from the subroutine
      logical lfailure
      integer iflag,activebox,nchainmove
      integer chaininfo
      dimension chaininfo(nchainmove,3)
c     --- local variables
      logical lvalid
      integer imolty,ichain,iunit,count,type,ilist,ichainmove,testcrd
      integer isolvtype,nmolty,nunit,ncmt
c     --- integer array
      integer atomcode
      dimension atomcode(NNTYPE)
      character*30 filename

#if DEBUG_WRITETRAMONTO
      write(6,*) 'starting writetramonto'
#endif

c     --- initialize atomcode
      do type = 1, NNTYPE
         atomcode(type) = -1
      enddo
c     --- initialize count
      count = -1

c     --- open dft_surfaces.dat and dft_decode.dat
      call twh_isolvtype(GLB_GET,isolvtype)
      if ( isolvtype .ne. SOLV_TRAMONTO ) then
         filename = 'dft_decode.dat'
         call twh_openfile(lfailure,71,filename,'DEFAULT','UNKNOWN'
     &        ,'DEFAULT')
         if ( lfailure ) then
            write(6,*) 'WRITETRAMONTO: failure to open file'
            return
         endif
      endif
      filename = 'dft_surfaces.dat'
      call twh_openfile(lfailure,72,filename,'DEFAULT','UNKNOWN'
     &     ,'DEFAULT')
      if ( lfailure ) then
         write(6,*) 'WRITETRAMONTO: failure to open file'
         return
      endif
c     --- output the information for the non-special chains
      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         call twh_ncmt(GLB_GET,activebox,imolty,ncmt)
         do ilist = 1,ncmt
            call twh_chainlist(GLB_GET,ilist,activebox,imolty,ichain)
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
                  call twh_writetramonto_atom(lfailure,ichain,imolty
     &                 ,iunit,count,CRD_REAL,atomcode)
                  if ( lfailure ) then
                     write(6,*) 'WRITETRAMONTO: failure in tram_atom'
                     return
                  endif
               enddo
            endif
         enddo
      enddo
c     --- output the information for the special chains
      do ichainmove = 1,nchainmove
         testcrd = chaininfo(ichainmove,3)
         if ( testcrd .ne. CRD_NONE ) then
            ichain = chaininfo(ichainmove,iflag)
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
               call twh_writetramonto_atom(lfailure,ichain,imolty,iunit
     &              ,count,testcrd,atomcode)
               if ( lfailure ) then
                  write(6,*) 'WRITETRAMNTO: failure in tram atom spec'
                  return
               endif
            enddo
         endif
      enddo

      close(71)
      rewind(72)
      close(72)

#if DEBUG_WRITETRAMONTO
      write(6,*) 'done with writetramonto'
#endif
      return
      end


      subroutine twh_writetramonto_atom(lfailure,ichain,imolty,iunit
     &     ,count,testcrd,atomcode)
c     ******************************************************************
c     * writes out the information for a single atom in the tramonto   *
c     * format                                                         *
c     *                                                                *
c     * originally written 08-04-2006 by M.G. Martin                   *
c     * last modified 08-12-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_CLASSICAL_POTENTIAL
#define FUNCTION_GET_NBCOEFF
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_QQATOM
#include "functions.h"
c     ---- variables passed to/from this subroutine
      logical lfailure
      integer ichain,imolty,iunit,count,testcrd
      integer atomcode
      dimension atomcode(NNTYPE)
c     --- local variables
      logical lcoulomb
      integer type,ij,isolvtype
      double precision xtemp,ytemp,ztemp

c     --- we are creating the initial file for tramonto
      call twh_ntype(GLB_GET,imolty,iunit,type)
      if ( atomcode(type) .eq. -1 ) then
c        --- we do not yet have an atom of this type, need 
c        --- output to dft_decode.dat
         count = count + 1
         atomcode(type) = count
         call twh_isolvtype(GLB_GET,isolvtype)
         if ( isolvtype .ne. SOLV_TRAMONTO ) then
            if ( twh_get_classical_potential() .eq. 'Lennard-Jones' .or.
     &           twh_get_classical_potential() .eq. 'UFF 12-6' ) then
c              --- LJ 12-6 output
               ij = twh_getnbtype(type,type)
               write(71,*) atomcode(type),twh_get_nbcoeff(1,ij)
     &              ,twh_get_nbcoeff(2,ij)
            else
c              --- not coded up for tramonto yet
               write(6,*) 'WRITETRAMONTO_ATOM: classical_potential '
     &              ,twh_get_classical_potential()
     &              ,'not coded up for tramonto output.'
               lfailure = .true.
               return
            endif
         endif
      endif
      call twh_coordinates(GLB_GET,testcrd,ichain,iunit,xtemp,ytemp
     &     ,ztemp)
c     --- write to dft_surf.dat
      call twh_lcoulomb(GLB_GET,lcoulomb)
      if ( lcoulomb ) then
         write(72,22) atomcode(type),ichain-1,xtemp,ytemp,ztemp
     &        ,twh_get_qqatom(imolty,iunit)
      else
         write(72,22) atomcode(type),ichain-1,xtemp,ytemp,ztemp,0.0d0
      endif
 22   format(i6,1x,i6,1x,4(f10.5,1x))

      return
      end
