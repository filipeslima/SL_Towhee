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
      subroutine twh_comparestruc(imolty,jmolty,lsame)
c     ******************************************************************
c     * this subroutine checks to see if two molecules have identical  *
c     * molecule templates                                             *
c     *                                                                *
c     * originally written 07-03-2002 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJAA0
#define FUNCTION_GET_IJAA1
#define FUNCTION_GET_IJAA2
#define FUNCTION_GET_IJBEN2
#define FUNCTION_GET_IJBEN3
#define FUNCTION_GET_IJIMPROP2
#define FUNCTION_GET_IJIMPROP3
#define FUNCTION_GET_IJIMPROP4
#define FUNCTION_GET_IJTOR2
#define FUNCTION_GET_IJTOR3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_IJVIB
#define FUNCTION_GET_INAA
#define FUNCTION_GET_INBEN
#define FUNCTION_GET_INIMPROP
#define FUNCTION_GET_INTOR
#define FUNCTION_GET_INVIB
#define FUNCTION_GET_ITAA
#define FUNCTION_GET_ITBEN
#define FUNCTION_GET_ITIMPROP
#define FUNCTION_GET_ITTOR
#define FUNCTION_GET_ITVIB
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lsame
      integer imolty,jmolty,nmaxi,nmaxj

c     --- local variables
      integer unit,ivib,iben,itor,iaa,iunit,junit

c     --- initialize lsame
      lsame = .true.

      call twh_nunit(GLB_GET,imolty,iunit)
      call twh_nunit(GLB_GET,jmolty,junit)
      if ( iunit .ne. junit ) then
         lsame = .false.
         return
      endif

      call twh_nmaxcbmc(GLB_GET,imolty,nmaxi)
      call twh_nmaxcbmc(GLB_GET,jmolty,nmaxj)
      if ( nmaxi .ne. nmaxj ) then
         lsame = .false.
         return
      endif

      do unit = 1,iunit
c        --- vibrations
         if ( twh_get_invib(imolty,unit) 
     &        .ne. twh_get_invib(jmolty,unit) ) then
            lsame = .false.
            return
         endif
         do ivib = 1,twh_get_invib(imolty,unit)
            if ( twh_get_ijvib(imolty,unit,ivib) 
     &           .ne. twh_get_ijvib(jmolty,unit,ivib) 
     &           .or. twh_get_itvib(imolty,unit,ivib)
     &           .ne.twh_get_itvib(imolty,unit,ivib)
     &           ) then
               lsame = .false.
               return
            endif
         enddo
c        --- angles
         if ( twh_get_inben(imolty,unit) 
     &        .ne. twh_get_inben(jmolty,unit) ) then
            lsame = .false.
            return
         endif
         do iben = 1,twh_get_inben(imolty,unit)
            if ( twh_get_ijben2(imolty,unit,iben) 
     &           .ne. twh_get_ijben2(jmolty,unit,iben)
     &           .or. ( twh_get_ijben3(imolty,unit,iben)
     &           .ne. twh_get_ijben3(jmolty,unit,iben))
     &           .or. ( twh_get_itben(imolty,unit,iben)
     &           .ne. twh_get_itben(jmolty,unit,iben)) ) then
               lsame = .false.
               return
            endif
         enddo
c        --- torsions
         if ( twh_get_intor(imolty,unit) 
     &        .ne. twh_get_intor(jmolty,unit) ) then
            lsame = .false.
            return
         endif
         do itor = 1,twh_get_intor(imolty,unit)
            if ( twh_get_ijtor2(imolty,unit,itor) 
     &           .ne. twh_get_ijtor2(jmolty,unit,itor)
     &           .or. ( twh_get_ijtor3(imolty,unit,itor) 
     &           .ne. twh_get_ijtor3(jmolty,unit,itor) )
     &           .or. ( twh_get_ijtor4(imolty,unit,itor) 
     &           .ne. twh_get_ijtor4(jmolty,unit,itor) )
     &           .or. ( twh_get_ittor(imolty,unit,itor) 
     &           .ne. twh_get_ittor(jmolty,unit,itor) ) ) then
               lsame = .false.
               return
            endif
         enddo
c        --- improper torsions
         if ( twh_get_inimprop(imolty,unit) 
     &        .ne. twh_get_inimprop(jmolty,unit) ) then
            lsame = .false.
            return
         endif
         do itor = 1,twh_get_inimprop(imolty,unit)
            if ( twh_get_ijimprop2(imolty,unit,itor) 
     &           .ne. twh_get_ijimprop2(jmolty,unit,itor)
     &           .or. ( twh_get_ijimprop3(imolty,unit,itor) 
     &           .ne. twh_get_ijimprop3(jmolty,unit,itor) )
     &           .or. ( twh_get_ijimprop4(imolty,unit,itor) 
     &           .ne. twh_get_ijimprop4(jmolty,unit,itor) )
     &           .or. ( twh_get_itimprop(imolty,unit,itor) 
     &           .ne. twh_get_itimprop(jmolty,unit,itor) ) ) then
               lsame = .false.
               return
            endif
         enddo
c        --- angle-angle
         if ( twh_get_inaa(imolty,unit) 
     &        .ne. twh_get_inaa(jmolty,unit) ) then
            lsame = .false.
            return
         endif
         do iaa = 1,twh_get_inaa(imolty,unit)
            if ( twh_get_ijaa0(imolty,unit,iaa) 
     &           .ne. twh_get_ijaa0(jmolty,unit,iaa)
     &           .or. 
     &           twh_get_ijaa1(imolty,unit,iaa) 
     &           .ne. twh_get_ijaa1(jmolty,unit,iaa)
     &           .or. 
     &           twh_get_ijaa2(imolty,unit,iaa) 
     &           .ne. twh_get_ijaa2(jmolty,unit,iaa)
     &           .or. twh_get_itaa(imolty,unit,iaa) 
     &           .ne. twh_get_itaa(jmolty,unit,iaa)
     &           ) then
               lsame = .false.
               return
            endif
         enddo
      enddo

      return
      end

