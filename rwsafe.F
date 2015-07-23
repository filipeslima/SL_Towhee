#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C)2005-2011 Marcus G. Martin                        *
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
      subroutine twh_rwsafe(lfailure,rwflag)
c     ******************************************************************
c     * reads or writes all of the self-adapting fixed endpoint bias   *
c     * information                                                    *
c     * rwflag is either READ_FLAG or WRITE_FLAG                       *
c     *                                                                *
c     * originally written 10-13-2005 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_INBEN
#define FUNCTION_GET_INTOR
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer rwflag
c     --- local variables
c     --- character string scalars
      character*10 bondstyle
      character*15 idname
      character*30 filename,rwformat
      character*50 goodstring,localstyle,twobbias_style,twobbias_init
      character*50 threebbias_init,threebbias_style
c     --- logical scalars
      logical lreturn
c     --- integer scalars
      integer rwfile,currentversion,iversion,itest,iunit,imolty
      integer numbend,ibend,pstart,pend,iparam,numdihed,idihed
      integer nmolty,nunit
c     --- double precision arrays
      double precision dpararray
      dimension dpararray(MAXTWOBOND)
      double precision dthreearray
      dimension dthreearray(MAXTHREEBOND)

c     --- see if we need to do anything here
      lreturn = .true.
      if ( rwflag .eq. READ_FLAG ) then
         call twh_two_bond_bias_initial_value(GLB_GET,twobbias_init)
         if ( twobbias_init .eq. 'file' ) then
            lreturn = .false.
         endif
         call twh_three_bond_bias_initial_value(GLB_GET,threebbias_init)
         if ( threebbias_init .eq. 'file' ) then
            lreturn = .false.
         endif
      elseif (rwflag .eq. WRITE_FLAG) then
c        --- alway write
         lreturn = .false.
      else
         write(6,*) "Unknown rwflag in RWSAFE: ", rwflag
         stop
      endif
      if ( lreturn ) return

c     --- set subroutine idname
      idname = 'rwsafe'
c     --- set the file number
      rwfile = 91
c     --- set the current version of towhee_ff
      currentversion = 1
c     --- open the forcefield file
      if ( rwflag .eq. WRITE_FLAG ) then
c        --- writing a safe file
         filename = 'towhee_safe_final'
      else
c        --- reading a forcefield file
         filename = 'towhee_safe_initial'
      endif
      call twh_openfile(lfailure,rwfile,filename,'DEFAULT','UNKNOWN'
     &     ,'DEFAULT')
      if ( lfailure ) return
c     --- set bondstyle to both for the initial reads
      bondstyle = 'both'
c     --- read the towhee_safe version
      goodstring = 'towhee SAFE Version'
      rwformat = '(a19)'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(rwfile,rwformat) goodstring
         write(rwfile,*) currentversion
c        --- set iversion to the currentversion
         iversion = currentversion
      else
         iversion = twh_read_labeled_integer(lfailure,rwfile,
     &         goodstring, .true.,idname)
         if ( lfailure ) return
         if ( iversion .gt. currentversion ) then
            write(6,*) 'invalid towhee_safe version ',iversion
            stop
         endif
      endif

c     --- number of molecule types
      goodstring = 'Number of molecule types'
      rwformat = '(a24)'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(rwfile,rwformat) goodstring
         call twh_nmolty(GLB_GET,nmolty)
         write(rwfile,*) nmolty
      else
         itest =  twh_read_labeled_integer(lfailure,rwfile,
     &         goodstring, .true.,idname)
         if ( lfailure ) return
         call twh_nmolty(GLB_GET,nmolty)
         if ( itest .ne. nmolty ) then
            call twh_rwsafeproblem(bondstyle)
            close(rwfile)
            return
         endif
      endif

      do imolty = 1,nmolty
         goodstring = 'Number of units'
         rwformat = '(a15)'
         call twh_nunit(GLB_GET,imolty,nunit)
         if ( rwflag .eq. WRITE_FLAG ) then
            write(rwfile,rwformat) goodstring
            write(rwfile,*) nunit
         else
            itest =  twh_read_labeled_integer(lfailure,rwfile,
     &            goodstring, .true.,idname)
            if ( lfailure ) return
            if ( itest .ne. nunit ) then
               call twh_rwsafeproblem(bondstyle)
               close(rwfile)
               return
            endif
         endif
      enddo

c     --- two_bond biasing style
      call twh_two_bond_fe_bias_style(GLB_GET,twobbias_style)
      goodstring = 'Two bond biasing style'
      rwformat = '(a22)'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(rwfile,rwformat) goodstring
         write(rwfile,'(a50)') twobbias_style
         localstyle = twobbias_style
      else
         call twh_read_labeled_string(lfailure,rwfile,
     &         goodstring, .true.,idname,localstyle)
         if ( lfailure ) return
      endif

c     --- set some values based upon the biasing
      if ( localstyle .eq. 
     &     'self adapting gaussian using 1-3 distance' ) then
         pstart = 1
         pend = 2
      elseif ( localstyle .eq. 'none' 
     &        .or. localstyle .eq. 'autofit gaussian'
     &        .or. localstyle .eq. 'analytic Boltzmann using angles'
     &        .or. 
     &        localstyle .eq. 'analytic Boltzmann dihedral energy sum'
     &        ) then
c        --- nothing to input or output
         pstart = 0
         pend = 0
      else
         write(6,*) 'RWSAFE: unknown two bond biasing style'
         bondstyle = 'both'
         call twh_rwsafeproblem(bondstyle)
         close(rwfile)
         return
      endif

      if ( pstart .ne. 0 ) then
         do imolty = 1,nmolty
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
               goodstring = 'Unit number and bend count'
               rwformat = '(a26)'
               if ( rwflag .eq. WRITE_FLAG ) then
                  write(rwfile,rwformat) goodstring
                  write(rwfile,*) iunit,twh_get_inben(imolty,iunit)
c                 --- set numbend
                  numbend = twh_get_inben(imolty,iunit)
               else
                  call twh_read_label(lfailure,rwfile,goodstring
     &                 ,.false.,idname)
                  if ( lfailure ) return
                  read(rwfile,*) itest,numbend
               endif
               do ibend = 1,numbend
                  if ( rwflag .eq. WRITE_FLAG ) then
                     do iparam = pstart,pend
                        call twh_two_bond_bias_info(GLB_GET,imolty
     &                       ,iunit,ibend,iparam,dpararray(iparam))
                     enddo
                     write(rwfile,*) (dpararray(iparam)
     &                    ,iparam=pstart,pend)
                  else
                     read(rwfile,*) (dpararray(iparam)
     &                    ,iparam=pstart,pend)
                     do iparam = pstart,pend
                        call twh_two_bond_bias_info(GLB_SET,imolty
     &                       ,iunit,ibend,iparam,dpararray(iparam))
                     enddo
                  endif
               enddo
            enddo
         enddo
      endif

      if ( localstyle .ne. twobbias_style ) then
c        --- the information we just read is useless, reset but keep
c        --- reading from the file for the three bond information
         bondstyle = 'two'
         call twh_rwsafeproblem(bondstyle)
      endif

c     --- three bond biasing
      call twh_three_bond_fe_bias_style(GLB_GET,threebbias_style)
      goodstring = 'Three bond biasing style'
      rwformat = '(a24)'
      if ( rwflag .eq. WRITE_FLAG ) then
         write(rwfile,rwformat) goodstring
         write(rwfile,'(a50)') threebbias_style
         localstyle = threebbias_style
      else
         call twh_read_labeled_string(lfailure,rwfile,
     &         goodstring, .true.,idname,localstyle)
         if ( lfailure ) return
      endif
c     --- set some values based on the localstyle
      if ( localstyle .eq. 'self adapting gaussian using 1-4 distance'
     &     ) then
         pstart = 1
         pend = 2
      elseif ( localstyle .eq. 'none' .or.
     &        localstyle .eq. 'analytic using max and min 2-4 distance'
     &        .or. localstyle 
     &        .eq. 'autofit gaussian using max and min 2-4 distance'
     &        ) then
         pstart = 0
         pend = 0
      else
         write(6,*) 'RWSAFE: unknown three bond biasing style'
         bondstyle = 'three'
         call twh_rwsafeproblem(bondstyle)
         close(rwfile)
         return
      endif

      if ( pstart .ne. 0 ) then
         do imolty = 1,nmolty
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
               goodstring = 'Unit number and dihedral count'
               rwformat = '(a30)'
               if ( rwflag .eq. WRITE_FLAG ) then
                  write(rwfile,rwformat) goodstring
                  write(rwfile,*) iunit,twh_get_intor(imolty,iunit)
                  numdihed = twh_get_intor(imolty,iunit)
               else
                  call twh_read_label(lfailure,rwfile,goodstring
     &                 ,.false.,idname)
                  if ( lfailure ) return
                  read(rwfile,*) itest,numdihed
               endif
               do idihed = 1,numdihed
                  if ( rwflag .eq. WRITE_FLAG ) then
                     do iparam = pstart,pend
                        call twh_three_bond_bias_info(GLB_GET,imolty
     &                       ,iunit,idihed,iparam,dthreearray(iparam))
                     enddo
                     write(rwfile,*) (dthreearray(iparam)
     &                    ,iparam=pstart,pend)
                  else
                     read(rwfile,*) (dthreearray(iparam)
     &                    ,iparam=pstart,pend)
                     do iparam = pstart,pend
                        call twh_three_bond_bias_info(GLB_SET,imolty
     &                       ,iunit,idihed,iparam,dthreearray(iparam))
                     enddo
                  endif
               enddo
            enddo
         enddo
      endif

      if ( localstyle .ne. threebbias_style ) then
c        --- the stuff we read is not useful.  Reset to autofit gaussian
         bondstyle = 'three'
         call twh_rwsafeproblem(bondstyle)
      endif

      close(rwfile)

      return
      end

      subroutine twh_rwsafeproblem(bondstyle)
c     ******************************************************************
c     * outputs a brief error message and resets the initial value of  *
c     * the two bond bias to autofit gaussian so we can use that for   *
c     * the initial safe parameters                                    *
c     *                                                                *
c     * originally written 10-13-2005 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*10 bondstyle

      write(6,*) 'RWSAFE: problem with towhee_safe_initial'
      write(6,*) 'resetting to use autofit gaussian instead'
      if ( bondstyle .eq. 'two' .or. bondstyle .eq. 'both' ) then
         call twh_two_bond_bias_initial_value(GLB_SET
     &        ,'autofit gaussian')
      endif
      if ( bondstyle .eq. 'three' .or. bondstyle .eq. 'both' ) then
         call twh_three_bond_bias_initial_value(GLB_SET
     &        ,'autofit gaussian')
      endif

      return
      end
