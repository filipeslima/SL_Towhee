#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2007-2012 Marcus G. Martin                       *
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
      subroutine twh_rwexternal(lfailure,rwflag,ifile)
c     ******************************************************************
c     * reads or writes all input parameters from towhee_input that are*
c     * related to the external potential style                        *
c     * rwflag is either READ_FLAG or WRITE_FLAG                       *
c     *                                                                *
c     * rewritten from readquantum 01-12-2007 by M.G. Martin           *
c     * last modifed 01-10-2012 by M.G. Martin                         *
c     ******************************************************************
#ifdef USEOPENKIM
#include "KIMstatus.h"
      use KIMservice
#endif
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ifile,rwflag
c     --- local variables
      character*15 idname,rwformat
      character*20 external_code
      character*50 goodstring,kim_model
c     --- cray pointers
#ifdef USEOPENKIM
      integer ierr
      integer(kind=kim_intptr) pkim
#endif

c     --- check rwflag
      if ( rwflag .ne. READ_FLAG .and. rwflag .ne. WRITE_FLAG ) then
         write(6,*) 'RWEXTERNAL: unknown rwflag'
         write(6,*) rwflag
         lfailure = .true.
         return
      endif
c     --- set subroutine idname
      idname = 'rwexternal'

      goodstring = 'external_code'
      if ( rwflag .eq. READ_FLAG ) then
         call twh_read_labeled_string(lfailure,ifile,goodstring,.true.
     &        ,idname,external_code)
         if ( lfailure ) return
         call twh_external_code(GLB_SET,external_code)
      else
         rwformat = '(a11)'
         write(ifile,rwformat) goodstring
         call twh_external_code(GLB_GET,external_code)
         write(ifile,'(a1,a20,a1)') char(39),external_code,char(39)
      endif

      if ( external_code .eq. 'lcao' ) then
         call twh_readquantum(lfailure,rwflag,ifile)
         if ( lfailure ) return
      elseif ( external_code .eq. 'KIM' ) then
         goodstring = 'kim_model'
         if ( rwflag .eq. READ_FLAG ) then
            call twh_read_labeled_string(lfailure,ifile,goodstring
     &           ,.true.,idname,kim_model)
            if ( lfailure ) return
            call twh_kim_model(GLB_SET,kim_model)
         else
            call twh_kim_model(GLB_GET,kim_model)
            rwformat = '(a9)'
            write(ifile,rwformat) goodstring
            write(ifile,'(a1,a,a1)') kim_model
         endif
         
#ifdef USEOPENKIM
         ierr = kim_api_init_f(pkim,'test_towhee',kim_model)
         call twh_pointerkim(GLB_SET,pkim)
         if ( ierr .lt. KIM_STATUS_OK ) then
            call kim_api_report_error_f(__LINE__,__FILE__
     &           ,"kim_api_init_f",ierr)
            lfailure = .true.
            return
         endif
#endif
      else
         write(6,*) 'RWEXTERNAL: invalid external_code'
         lfailure = .true.
         return
      endif

      return
      end
