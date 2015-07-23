#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2008 Marcus G. Martin                       *
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
      subroutine twh_writeruntime(press,nnn)
c     ******************************************************************
c     * outputs some runtime information to the standard output        *
c     *                                                                *
c     * split from old pupdate 10-13-2005 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_DERIVATIVE
#define FUNCTION_GET_NCMT
#include "functions.h"
c     --- variables passed to and from the subroutine
      integer nnn,l
      double precision press
      dimension press(MAXBOX)
      double precision v(MAX_FOREIGN_LAMBDA)
c     --- local variables
      logical scp_eval_derivative
      integer imolty,ibox,nmolty,scp_num_foreign_lambda,potentialstyle
      integer numboxes
      double precision energy,boxvolume

c     --- retrieve constants
      call twh_potentialstyle(GLB_GET,potentialstyle)
      call twh_nmolty(GLB_GET,nmolty)
      call twh_numboxes(GLB_GET,numboxes)

c     --- write out runtime information
      if ( potentialstyle .eq. POT_INTERNAL ) then
         call twh_boxvclassic(GLB_GET,1,energy)
      elseif ( potentialstyle .eq. POT_EXTERNAL ) then
         call twh_boxvquantum(GLB_GET,1,energy)
      endif
      call twh_scp_eval_derivative(GLB_GET,scp_eval_derivative)
      if ( scp_eval_derivative ) then
c        --- for du/d[lambda], do only first box.  
c        --- This can be changed if need be.
         call twh_boxvolume(GLB_GET,1,boxvolume)
         write(6,110) nnn, 1,energy,boxvolume,press(1),
     &         scp_get_derivative(1,LAMBDA_LJ),
     &         scp_get_derivative(1,LAMBDA_C),
     &         (twh_get_ncmt(1,imolty),imolty=1,nmolty)
 110     format(i10,1x,'B:',i2,1x,e15.8,1x,e12.5,1x,f12.1,1x,
     &        e12.5,1x,e12.5,NTMAX(i5))
      else
         call twh_boxvolume(GLB_GET,1,boxvolume)
         write(6,111) nnn, 1,energy,boxvolume,press(1)
     &        ,(twh_get_ncmt(1,imolty),imolty=1,nmolty)
 111     format(i10,1x,'B:',i2,1x,e11.4,1x,e11.4,1x,f12.1,1x,
     &        NTMAX(i5))
      endif

c     write the foreign energies on a separate line, box 1 only.  Ugly...      
      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
      if (scp_num_foreign_lambda .ne. 0) then
         call scp_get_foreign_energies(1,v)
         write(6,115) (v(l),l=1,scp_num_foreign_lambda)
115      format(' fu:', MAX_FOREIGN_LAMBDA(e14.6))
      endif


      do ibox = 2,numboxes
         if ( potentialstyle .eq. POT_INTERNAL ) then
            call twh_boxvclassic(GLB_GET,ibox,energy)
         elseif ( potentialstyle .eq. POT_EXTERNAL ) then
            call twh_boxvquantum(GLB_GET,ibox,energy)
         endif
         call twh_boxvolume(GLB_GET,ibox,boxvolume)
         write(6,112) ibox,energy,boxvolume
     &        ,press(ibox),(twh_get_ncmt(ibox,imolty),imolty=1,nmolty)
      enddo
 112  format(11x,'B:',i2,1x,e11.4,1x,e11.4,1x,f12.1,1x,NTMAX(i5))
      
      return
      end


      subroutine twh_writeruntime_header()
c     ******************************************************************
c     * Writes header information for the twh_writeruntime call        *
c     *                                                                *
c     * split from old pupdate 10-13-2005 by M.G. Martin               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- local variables
      logical scp_eval_derivative
      character*6 ctemp
      character*8 stepstyle
      integer scp_num_foreign_lambda

      call twh_stepstyle(GLB_GET,stepstyle)
      if ( stepstyle .eq. 'cycles' ) then
         ctemp = 'Cycle'
      elseif ( stepstyle .eq. 'moves' ) then
         ctemp = 'Move '
      else
c        --- unknown ctemp
         ctemp = '?????'
      endif

      call twh_scp_eval_derivative(GLB_GET,scp_eval_derivative)
      if (scp_eval_derivative) then 
         write(6,'(a5,6x,a3,3x,a10,2x,a12,1x,a12,1x,a9,1x,a8,1x,a7)') 
     &        ctemp,'Box','Energy [K]','Volume [A^3]'
     &        ,'Press. [kPa]','dU/dl_lj','dU/dl_c','Molecules'
      else
         write(6,'(a5,6x,a3,3x,a10,2x,a12,1x,a12,1x,a9)') 
     &        ctemp,'Box','Energy [K]','Volume [A^3]'
     &        ,'Press. [kPa]','Molecules'
      endif

c  very ugly.
      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
      if (scp_num_foreign_lambda .ne. 0) then
         write(6,*) '   Foreign energies, by lambda index'
      endif

      return
      end

