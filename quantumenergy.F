#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2011 Marcus G. Martin                       *
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
      subroutine twh_quantumenergy(lfailure,ibox,iflag,nchainmove
     &     ,chaininfo,vquantum)
c     ******************************************************************
c     * computed the energy using quantum mechanics                    *
c     *                                                                *
c     * originally written 11-04-2004 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- debug toggle
#define DEBUG_QUANTUMENERGY 0
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ibox,iflag,nchainmove,chaininfo
      dimension chaininfo(nchainmove,3)
      double precision vquantum
c     --- local variables
      integer natoms
      character*20 externalcode
      double precision gsenergy

      gsenergy = 0.0d0
      natoms = 0
      call twh_external_code(GLB_GET,externalcode)
      if ( externalcode .eq. 'lcao' ) then
c        --- output the lcao.in file with the coordinates
         call twh_writelcao(lfailure,ibox,iflag,nchainmove,chaininfo
     &        ,natoms,gsenergy)
         if ( lfailure ) return
      endif

      if ( natoms .eq. 0 ) then
         vquantum = 0.0d0
      else
         call twh_quantumwrapper(lfailure,vquantum)
         if ( lfailure ) return
c        --- subtract the base atom energy
         vquantum = vquantum - gsenergy
#if DEBUG_QUANTUMENERGY
         write(6,*) 'vquantum:',vquantum
#endif
      endif

      return
      end

      subroutine twh_quantumwrapper(lfailure,vquantum)
c     ******************************************************************
c     * handles the call to the quantum mechanics codes that are       *
c     * compiled as libraries                                          *
c     *                                                                *
c     * originally written 02-13-2005 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#ifdef USEMPI
#include "mpif.h"
#endif
c     --- variables passed to/from the subroutine
      logical lfailure
      double precision vquantum
c     --- local variables
#ifdef USEMPI
      integer imessage,ierr
#endif
#ifdef USELCAO
      logical lsubflag
#endif
      character*(20) externalcode

#ifdef USEMPI
      imessage = 1
      call MPI_BCAST(imessage,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c     --- check error condition
      if ( ierr .eq. 1 ) then
         write(6,*) 'QUANTUMWRAPPER: error with broadcast'
         lfailure = .true.
         return
      endif
#endif

      vquantum = 0.0d0
      if ( externalcode .eq. 'lcao' ) then
#ifdef USELCAO
         lsubflag = .true.
         call subquest(vquantum,lsubflag)
c        --- convert into Kelvins
         vquantum = vquantum*CON_RYDBERGTOEV*CON_EVTOK
#else
         write(6,*) 'QUANTUMENERGY: LCAO not enabled'
         lfailure = .true.
         return
#endif
      else
         write(6,*) 'QUANTUMENERGY: unknown externalcode'
         write(6,*) externalcode
         lfailure = .true.
         return
      endif

      return
      end

      subroutine twhquantumloop()
c     ******************************************************************
c     * place where all of the processers that are not in charge wait  *
c     * for the signal to go into the quantum code                     *
c     *                                                                *
c     * originally written 02-13-2005 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#ifdef USEMPI
#include "mpif.h"
#endif
c     --- local variables
      logical lsubflag
      integer imessage
#ifdef USEMPI
      integer ierr,myrank
#endif
#ifdef USELCAO
      double precision vdummy
#endif
#ifdef USEMPI
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
#endif
      imessage = -1
      lsubflag = .true.
      do while ( imessage .ne. 0 ) 
c        --- get instructions from the master
#ifdef USEMPI
         call MPI_BCAST(imessage,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
c        --- execute instructions
         if ( imessage .eq. 1 ) then
#ifdef USELCAO
            vdummy = 0.0d0
            call subquest(vdummy,lsubflag)
#endif
         endif
      enddo
      return
      end

      subroutine twhendquantum()
c     ******************************************************************
c     * the master node comes here to send the terminate message to the*
c     * workers                                                        *
c     *                                                                *
c     * originally written 02-15-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#ifdef USEMPI
#include "mpif.h"
#endif
c     --- local variables
      integer imessage
#ifdef USEMPI
      integer ierr
#endif

      imessage = 0
#ifdef USEMPI
      call  MPI_BCAST(imessage,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      return
      end
