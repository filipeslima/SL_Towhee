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
      subroutine twh_readquantum(lfailure,rwflag,ifile)
c     ******************************************************************
c     * reads or writes all input parameters from towhee_input that are*
c     * related to the lcao external_code                              *
c     * rwflag is either READ_FLAG or WRITE_FLAG                       *
c     *                                                                *
c     * originally written 11-01-2004 by M.G. Martin                   *
c     * last modifed 08-04-2011 by M.G. Martin                         *
c     ******************************************************************
      implicit none
#define FUNCTION_READ_LABELED_FLOAT
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"

c     --- variables passed to/from the subroutine
      logical lfailure
      integer ifile,rwflag
c     --- local variables
      integer itype,lcao_atomtypes
      character*15 idname,rwformat
      character*20 twentychar
      character*50 goodstring
      double precision dvalue

c     --- check rwflag
      if ( rwflag .ne. READ_FLAG .and. rwflag .ne. WRITE_FLAG ) then
         write(6,*) 'READQUANTUM: unknown rwflag'
         write(6,*) rwflag
         lfailure = .true.
         return
      endif
c     --- set subroutine idname
      idname = 'readquantum'

c     --- make sure we have compiled the code properly for LCAO
#ifndef USELCAO
      write(6,*) 'You must build with the LCAO configure option'
      write(6,*) 'in order to use this functionallity'
      write(6,*) './configure --enable-lcao'
      lfailure = .true.
      return
#endif
         
      goodstring = 'lcao_functional'
      if ( rwflag .eq. READ_FLAG ) then
         call twh_read_labeled_string(lfailure,ifile,goodstring,.true.
     &        , idname,twentychar)
         if ( lfailure ) return
         call twh_lcao_functional(GLB_SET,twentychar)
      else
         rwformat = '(a15)'
         write(ifile,rwformat) goodstring
         call twh_lcao_functional(GLB_GET,twentychar)
         write(ifile,'(a1,a20,a1)') char(39),twentychar,char(39)
      endif
         
      goodstring = 'lcao_atomtypes'
      if ( rwflag .eq. READ_FLAG ) then
         lcao_atomtypes = twh_read_labeled_integer(lfailure,ifile
     &        ,goodstring,.true., idname)
         if ( lfailure ) return
         call twh_lcao_atomtypes(GLB_SET,lcao_atomtypes)
c        --- safty check
         if ( lcao_atomtypes .gt. NNTYPE ) then
            write(6,*) 'READQUANTUM: lcao_atomtypes > NNTYPE'
            lfailure = .true.
            return
         endif
      else
         rwformat = '(a14)'
         write(ifile,rwformat) goodstring
         call twh_lcao_atomtypes(GLB_GET,lcao_atomtypes)
         write(ifile,*) lcao_atomtypes
      endif
      
      goodstring = 'lcao_filename'
      if ( rwflag .eq. READ_FLAG ) then
         call twh_read_label(lfailure,ifile,goodstring, .false., idname)
         if ( lfailure ) return
         do itype = 1,lcao_atomtypes
            call twh_read_string(lfailure,ifile,goodstring,.false.
     &           ,twentychar)
            if ( lfailure ) return
            call twh_lcao_filename(GLB_SET,itype,twentychar)
            write(6,'(a14,i3,1x,a1,1x,a50)') 'lcao_filename:'
     &           ,itype,'=',twentychar
         enddo
      else
         rwformat = '(a13)'
         write(ifile,rwformat) goodstring
         do itype = 1,lcao_atomtypes
            call twh_lcao_filename(GLB_GET,itype,twentychar)
            write(ifile,'(a1,a50,a1)') char(39),twentychar,char(39)
         enddo
      endif

      goodstring = 'lcao_gridmultiplier'
      if ( rwflag .eq. READ_FLAG ) then
         dvalue = twh_read_labeled_float(lfailure,ifile,goodstring
     &        ,.true.  , idname)
         if ( lfailure ) return
         call twh_lcao_gridmultiplier(GLB_SET,dvalue)
      else
         rwformat = '(a19)'
         write(ifile,rwformat) goodstring
         call twh_lcao_gridmultiplier(GLB_GET,dvalue)
         write(ifile,*) dvalue
      endif

      goodstring = 'lcao_kgridproduct'
      if ( rwflag .eq. READ_FLAG ) then
         dvalue = twh_read_labeled_float(lfailure,ifile,goodstring
     &        ,.true.  , idname)
         if ( lfailure ) return
         call twh_lcao_kgridproduct(GLB_SET,dvalue)
      else
         rwformat = '(a17)'
         write(ifile,rwformat) goodstring
         call twh_lcao_kgridproduct(GLB_GET,dvalue)
         write(ifile,*) dvalue
      endif

      return
      end
