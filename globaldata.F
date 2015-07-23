#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2006-2014 Marcus G. Martin                       *
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

c     ******************************************************************
c     * this file contains all of the globally accessed data structures*
c     *                                                                *
c     * originally written 03-20-2006 by M.G. Martin                   *
c     * last modified 07-19-2014 by M.G. Martin                        *
c     ******************************************************************

      subroutine twh_aabead(flag,iaa,iatom,ivalue)
c     ******************************************************************
c     * atoms involved in angle-angle terms during a CBMC growth       *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modiied 10-05-2008 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,iaa,iatom,ivalue
c     --- storage variable
      integer aabead
      dimension aabead(NNBOND*MAXAA,4)
      save aabead

      if ( flag .eq. GLB_GET ) then
         ivalue = aabead(iaa,iatom)
      elseif ( flag .eq. GLB_SET ) then
         aabead(iaa,iatom) = ivalue
      else
         write(6,*) 'AABEAD: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_aacoeff(flag,itype,icoeff,dvalue)
c     ******************************************************************
c     * the coefficients for each angle-angle type in the simulation   *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modiied 04-13-2006 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,icoeff
      double precision dvalue
c     --- local variables
      integer ittemp,ictemp
c     --- storage variable
      double precision aacoeff
      dimension aacoeff(TAAMAX,MINAACOEFF:MAXAACOEFF)
      save aacoeff

      if ( flag .eq. GLB_SET ) then
         aacoeff(itype,icoeff) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = aacoeff(itype,icoeff)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            do ictemp = MINAACOEFF,icoeff
               aacoeff(ittemp,ictemp) = dvalue
            enddo
         enddo
      else
         write(6,*) 'AACOEFF: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_aaff(flag,itype,cvalue)
c     ******************************************************************
c     * the force field names for each angle-angle type in the         *
c     * simulation                                                     *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modiied 04-13-2006 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype
      character*(*) cvalue
c     --- local variables
      integer ittemp
c     --- storage variable
      character*(FFNAMELEN) aaff
      dimension aaff(TAAMAX)
      save aaff

      if ( flag .eq. GLB_SET ) then
         aaff(itype) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = aaff(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            aaff(itype) = cvalue
         enddo
      else
         write(6,*) 'AAFF: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_aaind(flag,iaa,ivalue)
c     ******************************************************************
c     * index of the angle-angle terms used during CBMC moves          *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modiied 10-05-2008 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,iaa,ivalue
c     --- storage variable
      integer aaind
      dimension aaind(NNBOND*MAXAA)
      save aaind

      if ( flag .eq. GLB_GET ) then
         ivalue = aaind(iaa)
      elseif ( flag .eq. GLB_SET ) then
         aaind(iaa) = ivalue
      else
         write(6,*) 'AAIND: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_aanames(flag,itype,isame,index,cvalue)
c     ******************************************************************
c     * the atom names for each angle-angle type in the simulation     *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modiied 04-13-2006 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,isame,index
      character*(*) cvalue
c     --- local variables
      integer ittemp,istemp,iitemp
c     --- storage variable
      character*(FFNAMELEN) aanames
      dimension aanames(TAAMAX,MAXSAMEAA,4)
      save aanames

      if ( flag .eq. GLB_SET ) then
         aanames(itype,isame,index) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = aanames(itype,isame,index)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            do istemp = 1,isame
               do iitemp = 1,index
                  aanames(ittemp,istemp,iitemp) = cvalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'AANAMES: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_aastyle(flag,itype,ivalue)
c     ******************************************************************
c     * the style of each angle-angle type in the simulation           *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modiied 04-13-2006 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer ittemp
c     --- storage variable
      integer aastyle
      dimension aastyle(TAAMAX)
      save aastyle

      if ( flag .eq. GLB_SET ) then
         aastyle(itype) = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = aastyle(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            aastyle(ittemp) = ivalue
         enddo
      else
         write(6,*) 'AASTYLE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_aminonum(flag,imolty,iunit,ivalue)
c     ******************************************************************
c     * integer count of the amino acid each atom is in, mainly used   *
c     * for output                                                     *
c     *                                                                *
c     * originally written 10-03-2008 by M.G. Martin                   *
c     * last modified 10-03-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ivalue
c     --- local variables
      integer imtemp,iutemp
c     --- information variables
      integer aminonum
      dimension aminonum(NTMAX,NUMAX)
      save aminonum

      if ( flag .eq. GLB_SET ) then
         aminonum(imolty,iunit) = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = aminonum(imolty,iunit)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               aminonum(imtemp,iutemp) = ivalue
            enddo
         enddo
      else
         write(6,*) 'AMINONUM: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_aminoshort(flag,imolty,iunit,cvalue)
c     ******************************************************************
c     * 3 character code used to generate pdb files, especially the    *
c     * amino acid portions                                            *
c     *                                                                *
c     * originally written 10-02-2008 by M.G. Martin                   *
c     * last modified 10-02-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit
      character*(*) cvalue
c     --- information variables
      character*3 aminoshort
      dimension aminoshort(NTMAX,NUMAX)
      save aminoshort

      if ( flag .eq. GLB_SET ) then
         aminoshort(imolty,iunit) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = aminoshort(imolty,iunit)
      else
         write(6,*) 'AMINOSHORT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_angleorder(flag,itype,cvalue)
c     ******************************************************************
c     * the angle order for each angle type for matching in assemble   *
c     *                                                                *
c     * originally written 12-12-2006 by M.G. Martin                   *
c     * last modified 12-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype
      character*(*) cvalue
c     --- local variables
      integer ittemp
c     --- information variables
      character*(15) angleorder
      dimension angleorder(TBENMAX)
      save angleorder

      if ( flag .eq. GLB_SET ) then
         angleorder(itype) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = angleorder(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            angleorder(ittemp) = cvalue
         enddo
      else
         write(6,*) 'ANGLEORDER: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_anglestyle(flag,itype,ivalue)
c     ******************************************************************
c     * the integer angle style for each type of angle in the          *
c     * simulation                                                     *
c     *                                                                *
c     * originally written 04-04-2006 by M.G. Martin                   *
c     * last modiied 04-04-2006 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer itemp
c     --- storage variable
      integer anglestyle
      dimension anglestyle(TBENMAX)
      save anglestyle

      if ( flag .eq. GLB_GET ) then
         ivalue = anglestyle(itype)
      elseif ( flag .eq. GLB_SET ) then
         anglestyle(itype) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do itemp = 1,itype
            anglestyle(itemp) = ivalue
         enddo
      else
         write(6,*) 'ANGLESTYLE: invalid global flag',flag
         stop
      endif

      return
      end


      subroutine twh_arbcmo(flag,ichain,avalue,bvalue,cvalue)
c     ******************************************************************
c     * the center of mass coordinates stored in arbitrary units for   *
c     * use during moves that scale all molecules by their center of   *
c     * mass (e.g. volume moves)                                       *
c     *                                                                *
c     * originally written 09-24-2008 by M.G. Martin                   *
c     * last modiied 06-23-2009 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ichain
      double precision avalue,bvalue,cvalue
c     --- local variables
      integer newindex
      double precision vector
      dimension vector(3)

      newindex = 3*(ichain-1)
      if ( flag .eq. GLB_GET ) then
         call twh_arbcmofield(GLB_GET_TRIPLE,newindex,vector)
         avalue = vector(1)
         bvalue = vector(2)
         cvalue = vector(3)
      elseif ( flag .eq. GLB_SET ) then
         vector(1) = avalue
         vector(2) = bvalue
         vector(3) = cvalue
         call twh_arbcmofield(GLB_SET_TRIPLE,newindex,vector)
      else
         write(6,*) 'ARBCMO: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_avb1rad(flag,dvalue)
c     ******************************************************************
c     * the radius for the aggregation-volume-bias move type 1         *
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modiied 12-24-2007 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      double precision dvalue
c     --- storage variable
      double precision avb1rad
      save avb1rad

      if ( flag .eq. GLB_GET ) then
         dvalue = avb1rad
      elseif ( flag .eq. GLB_SET ) then
         avb1rad = dvalue
      else
         write(6,*) 'AVB1RAD: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_avb2rad(flag,dvalue)
c     ******************************************************************
c     * the radius for the aggregation-volume-bias move type 2         *
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modiied 12-24-2007 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      double precision dvalue
c     --- storage variable
      double precision avb2rad
      save avb2rad

      if ( flag .eq. GLB_GET ) then
         dvalue = avb2rad
      elseif ( flag .eq. GLB_SET ) then
         avb2rad = dvalue
      else
         write(6,*) 'AVB2RAD: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_avb3rad(flag,dvalue)
c     ******************************************************************
c     * the radius for the aggregation-volume-bias move type 3         *
c     *                                                                *
c     * originally written 02-07-2008 by M.G. Martin                   *
c     * last modiied 02-07-2008 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      double precision dvalue
c     --- storage variable
      double precision avb3rad
      save avb3rad

      if ( flag .eq. GLB_GET ) then
         dvalue = avb3rad
      elseif ( flag .eq. GLB_SET ) then
         avb3rad = dvalue
      else
         write(6,*) 'AVB3RAD: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_baavb1(flag,imolty,jmolty,ibox,imove,icount)
c     ******************************************************************
c     * running count of accepted aggregation-volume-bias type 1 moves *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,jmolty,ibox,imove,icount
c     --- local variables
      integer imtemp,intemp,ibtemp,ivtemp
c     --- storage variable
      integer baavb1
      dimension baavb1(NTMAX,NTMAX,MAXBOX,0:3)
      save baavb1

      if ( flag .eq. GLB_INCR ) then
         baavb1(imolty,jmolty,ibox,imove) 
     &        = baavb1(imolty,jmolty,ibox,imove) + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = baavb1(imolty,jmolty,ibox,imove)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do intemp = 1,jmolty
               do ibtemp = 1,ibox
                  do ivtemp = 0,imove
                     baavb1(imtemp,intemp,ibtemp,ivtemp) = icount
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'BAAVB1: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_baavb2(flag,imolty,jmolty,ibox,imove,icount)
c     ******************************************************************
c     * running count of accepted aggregation-volume-bias type 2 moves *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,jmolty,ibox,imove,icount
c     --- local variables
      integer imtemp,intemp,ibtemp,ivtemp
c     --- storage variable
      integer baavb2
      dimension baavb2(NTMAX,NTMAX,MAXBOX,0:3)
      save baavb2

      if ( flag .eq. GLB_INCR ) then
         baavb2(imolty,jmolty,ibox,imove) 
     &        = baavb2(imolty,jmolty,ibox,imove) + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = baavb2(imolty,jmolty,ibox,imove)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do intemp = 1,jmolty
               do ibtemp = 1,ibox
                  do ivtemp = 0,imove
                     baavb2(imtemp,intemp,ibtemp,ivtemp) = icount
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'BAAVB2: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_baavb3(flag,imolty,jmolty,ibox,icount)
c     ******************************************************************
c     * running count of accepted aggregation-volume-bias type 3 moves *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,jmolty,ibox,icount
c     --- local variables
      integer imtemp,intemp,ibtemp
c     --- storage variable
      integer baavb3
      dimension baavb3(NTMAX,NTMAX,MAXBOX)
      save baavb3

      if ( flag .eq. GLB_INCR ) then
         baavb3(imolty,jmolty,ibox) = baavb3(imolty,jmolty,ibox) 
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = baavb3(imolty,jmolty,ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do intemp = 1,jmolty
               do ibtemp = 1,ibox
                  baavb3(imtemp,intemp,ibtemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BAAVB3: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_baback(flag,imolty,iunit,itor,ibox,icount)
c     ******************************************************************
c     * running count of accepted backbone configurational-bias        *
c     * regrowth moves since the last update                           *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibox,itor,icount
c     --- local variables
      integer imtemp,iutemp,ittemp,ibtemp
c     --- storage variable
      integer baback
      dimension baback(NTMAX,NUMAX,MAXTOR,MAXBOX)
      save baback

      if ( flag .eq. GLB_INCR ) then
         baback(imolty,iunit,itor,ibox) = baback(imolty,iunit,itor,ibox)
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = baback(imolty,iunit,itor,ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ittemp = 1,itor
                  do ibtemp = 1,ibox
                     baback(imtemp,iutemp,ittemp,ibtemp) = 0
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'BABACK: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bacb(flag,imolty,iunit,ibox,icount)
c     ******************************************************************
c     * running count of accepted configurational-bias regrowth moves  *
c     * since the last update                                          *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibox,icount
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer bacb
      dimension bacb(NTMAX,NUMAX,MAXBOX)
      save bacb

      if ( flag .eq. GLB_INCR ) then
         bacb(imolty,iunit,ibox) = bacb(imolty,iunit,ibox) + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bacb(imolty,iunit,ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibox
                  bacb(imtemp,iutemp,ibtemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BACB: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bacbside(flag,imolty,iamino,ibox,icount)
c     ******************************************************************
c     * running count of accepted configurational-bias side-chain      *
c     * regrowth moves                                                 *
c     *                                                                *
c     * originally written 06-14-2006 by M.G. Martin                   *
c     * last modified 06-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iamino,ibox,icount
c     --- local variables
      integer imtemp,iatemp,ibtemp
c     --- storage variable
      integer bacbside
      dimension bacbside(NTMAX,NUMAX,MAXBOX)
      save bacbside

      if ( flag .eq. GLB_INCR ) then
         bacbside(imolty,iamino,ibox) = bacbside(imolty,iamino,ibox)
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bacbside(imolty,iamino,ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iatemp = 1,iamino
               do ibtemp = 1,ibox
                  bacbside(imtemp,iatemp,ibtemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BACBSIDE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bacbswap(flag,imolty,inbox,outbox,icount)
c     ******************************************************************
c     * running count of accepted configurational-bias swap moves      *
c     * and also configurational-bias grand-canonical moves (involving *
c     * box 0) since the last update                                   *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,inbox,outbox,icount
c     --- local variables
      integer imtemp,iouttemp,iintemp
c     --- storage variable
      integer bacbswap
      dimension bacbswap(NTMAX,0:MAXBOX,0:MAXBOX)
      save bacbswap

      if ( flag .eq. GLB_INCR ) then
         bacbswap(imolty,inbox,outbox) = bacbswap(imolty,inbox,outbox)
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bacbswap(imolty,inbox,outbox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iintemp = 0,inbox
               do iouttemp = 0,outbox
                  bacbswap(imtemp,iintemp,iouttemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BACBSWAP: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_backlist(flag,imolty,iunit,ivalue)
c     ******************************************************************
c     * integer code to identify the backbone atoms                    *
c     *                                                                *
c     * originally written 10-03-2008 by M.G. Martin                   *
c     * last modified 10-03-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ivalue
c     --- storage variable
      integer backlist
      dimension backlist(NTMAX,0:NUMAX)
      save backlist

      if ( flag .eq. GLB_GET ) then
         ivalue = backlist(imolty,iunit)
      elseif ( flag .eq. GLB_SET ) then
         backlist(imolty,iunit) = ivalue
      else
         write(6,*) 'BACKLIST: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_backupfreq(flag,ivalue)
c     ******************************************************************
c     * frequency for outputting a backup restart file                 *
c     *                                                                *
c     * originally written 03-23-2006 by MAW                           *
c     * last modified 03-23-2006 by MAW                                *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer backupfreq
      save backupfreq

      if ( flag .eq. GLB_GET ) then
         ivalue = backupfreq
      elseif ( flag .eq. GLB_SET ) then
         backupfreq = ivalue
      else
         write(6,*) 'BACKUPFREQ: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_baconrot(flag,imolty,iunit,ibox,icount)
c     ******************************************************************
c     * running count of accepted concerted rotation moves since the   *
c     * last update                                                    *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibox,icount
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer baconrot
      dimension baconrot(NTMAX,NUMAX,MAXBOX)
      save baconrot

      if ( flag .eq. GLB_INCR ) then
         baconrot(imolty,iunit,ibox) = baconrot(imolty,iunit,ibox)
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = baconrot(imolty,iunit,ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibox
                  baconrot(imtemp,iutemp,ibtemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BACONROT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bacrback(flag,imolty,iunit,ibox,icount)
c     ******************************************************************
c     * running count of accepted backbone concerted rotation moves    *
c     * since the last update                                          *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibox,icount
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer bacrback
      dimension bacrback(NTMAX,NUMAX,MAXBOX)
      save bacrback

      if ( flag .eq. GLB_INCR ) then
         bacrback(imolty,iunit,ibox) = bacrback(imolty,iunit,ibox)
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bacrback(imolty,iunit,ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibox
                  bacrback(imtemp,iutemp,ibtemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BACRBACK: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bapivbin(flag,ibox,ibin,icount)
c     ******************************************************************
c     * running count of accepted pivot moves in each angle region bin *
c     *                                                                *
c     * originally written 03-21-2006 by M.G. Martin                   *
c     * last modified 03-21-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox,ibin,icount
c     --- local variables
      integer iboxtemp,ibintemp
c     --- storage variable
      integer bapivbin
      dimension bapivbin(MAXBOX,MAXPIVOTBIN)
      save bapivbin

      if ( flag .eq. GLB_INCR ) then
         bapivbin(ibox,ibin) = bapivbin(ibox,ibin) + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bapivbin(ibox,ibin)
      elseif ( flag .eq. GLB_INIT ) then
         do iboxtemp = 1,ibox
            do ibintemp = 1,ibin
               bapivbin(iboxtemp,ibintemp) = icount
            enddo
         enddo
      else
         write(6,*) 'BAPIVBIN: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bapivot(flag,imolty,iunit,ivib,ibox,icount)
c     ******************************************************************
c     * running count of accepted pivot moves                          *
c     *                                                                *
c     * originally written 03-21-2006 by M.G. Martin                   *
c     * last modified 03-21-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ivib,ibox,icount
c     --- local variables
      integer itemp,iutemp,ivtemp,ibtemp
c     --- storage variable
      integer bapivot
      dimension bapivot(NTMAX,NUMAX,NNBOND,MAXBOX)
      save bapivot

      if ( flag .eq. GLB_INCR ) then
         bapivot(imolty,iunit,ivib,ibox) 
     &        = bapivot(imolty,iunit,ivib,ibox) + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bapivot(imolty,iunit,ivib,ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do itemp = 1,imolty
            do iutemp = 1,iunit
               do ivtemp = 1,ivib
                  do ibtemp = 1,ibox
                     bapivot(itemp,iutemp,ivtemp,ibtemp) = icount
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'BAPIVOT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_barbswap(flag,imolty,inbox,outbox,icount)
c     ******************************************************************
c     * running count of accepted rotational-bias 2-box swap moves     *
c     * since the last update                                          *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,inbox,outbox,icount
c     --- local variables
      integer imtemp,iintemp,iouttemp
c     --- storage variable
      integer barbswap
      dimension barbswap(NTMAX,MAXBOX,MAXBOX)
      save barbswap

      if ( flag .eq. GLB_INCR ) then
         barbswap(imolty,inbox,outbox) = barbswap(imolty,inbox,outbox)
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = barbswap(imolty,inbox,outbox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iintemp = 1,inbox
               do iouttemp = 1,outbox
                  barbswap(imtemp,iintemp,iouttemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BARBSWAP: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bashift(flag,imove,ibox,arbindex,icount)
c     ******************************************************************
c     * running count of accepted shift moves since the last update    *
c     *                                                                *
c     * originally written 03-21-2006 by M.G. Martin                   *
c     * last modified 03-21-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imove,ibox,arbindex,icount
c     --- local variables
      integer itemp,ibtemp,iatemp
c     --- storage variable
      integer bashift
      dimension bashift(2,MAXBOX,3)
      save bashift

      if ( flag .eq. GLB_INCR ) then
         bashift(imove,ibox,arbindex) = bashift(imove,ibox,arbindex)
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bashift(imove,ibox,arbindex)
      elseif ( flag .eq. GLB_INIT ) then
         do itemp = 1,imove
            do ibtemp = 1,ibox
               do iatemp = 1,arbindex
                  bashift(itemp,ibtemp,iatemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BASHIFT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bavol(flag,imove,icount)
c     ******************************************************************
c     * running count of accepted volume moves since the last update   *
c     *                                                                *
c     * originally written 03-21-2006 by M.G. Martin                   *
c     * last modified 03-21-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imove,icount
c     --- local variables
      integer itemp
c     --- storage variable
      integer bavol
      dimension bavol(MAXPBOX)
      save bavol

      if ( flag .eq. GLB_INCR ) then
         bavol(imove) = bavol(imove) + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bavol(imove)
      elseif ( flag .eq. GLB_INIT ) then
         do itemp = 1,imove
            bavol(itemp) = icount
         enddo
      else
         write(6,*) 'BAVOL: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bencoeff(flag,itype,icoeff,dvalue)
c     ******************************************************************
c     * the coefficients that determine the equilibrium angle and the  *
c     * bending energy                                                 *
c     *                                                                *
c     * originally written 04-04-2006 by M.G. Martin                   *
c     * last modified 12-11-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,icoeff
      double precision dvalue
c     --- local variables
      integer ittemp,ictemp,postype
c     --- storage variable
      double precision bencoeff
      dimension bencoeff(TBENMAX,MINBENDCOEFF:MAXBENDCOEFF)
      save bencoeff

c     --- the type can be negative to signal the order, make sure it
c     --- is positive here
      postype = iabs(itype)

      if ( flag .eq. GLB_SET ) then
         bencoeff(postype,icoeff) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = bencoeff(postype,icoeff)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,postype
            do ictemp = MINBENDCOEFF,icoeff
               bencoeff(ittemp,ictemp) = dvalue
            enddo
         enddo
      else
         write(6,*) 'BENCOEFF: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bend_a_ideal_fraction(flag,dvalue)
c     ******************************************************************
c     * ideal fraction for the bending A distributions using certain   *
c     * CBMC growth options                                            *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      double precision dvalue
c     --- storage variable
      double precision bend_a_ideal_fraction
      save bend_a_ideal_fraction

      if ( flag .eq. GLB_GET ) then
         dvalue = bend_a_ideal_fraction
      elseif ( flag .eq. GLB_SET ) then
         bend_a_ideal_fraction = dvalue
      else
         write(6,*) 'BEND_A_IDEAL_FRACTION: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bend_a_sdev_multiplier(flag,dvalue)
c     ******************************************************************
c     * multiplier for the bending A standard deviation using certain  *
c     * CBMC growth options                                            *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      double precision dvalue
c     --- storage variable
      double precision bend_a_sdev_multiplier
      save bend_a_sdev_multiplier

      if ( flag .eq. GLB_GET ) then
         dvalue = bend_a_sdev_multiplier
      elseif ( flag .eq. GLB_SET ) then
         bend_a_sdev_multiplier = dvalue
      else
         write(6,*) 'BEND_A_SDEV_MULTIPLIER: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bend_b_ideal_fraction(flag,dvalue)
c     ******************************************************************
c     * fraction of the bending B distribution that is ideal for use   *
c     * with certain CBMC growth options                               *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      double precision dvalue
c     --- storage variable
      double precision bend_b_ideal_fraction
      save bend_b_ideal_fraction

      if ( flag .eq. GLB_GET ) then
         dvalue = bend_b_ideal_fraction
      elseif ( flag .eq. GLB_SET ) then
         bend_b_ideal_fraction = dvalue
      else
         write(6,*) 'BEND_B_IDEAL_FRACTION: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bend_b_sdev_multiplier(flag,dvalue)
c     ******************************************************************
c     * multiplier for the bending B standard deviation using certain  *
c     * CBMC growth options                                            *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      double precision dvalue
c     --- storage variable
      double precision bend_b_sdev_multiplier
      save bend_b_sdev_multiplier

      if ( flag .eq. GLB_GET ) then
         dvalue = bend_b_sdev_multiplier
      elseif ( flag .eq. GLB_SET ) then
         bend_b_sdev_multiplier = dvalue
      else
         write(6,*) 'BEND_B_SDEV_MULTIPLIER: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bendang(flag,ibond,dvalue)
c     ******************************************************************
c     * bending angles for use during a CBMC regrowth move             *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modified 10-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibond
      double precision dvalue
c     --- storage variable
      double precision bendang
      dimension bendang(NNBOND)
      save bendang

      if ( flag .eq. GLB_GET ) then
         dvalue = bendang(ibond)
      elseif ( flag .eq. GLB_SET ) then
         bendang(ibond) = dvalue
      else
         write(6,*) 'BENDANG: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bendff(flag,itype,cvalue)
c     ******************************************************************
c     * the name of the forcefield for each bending angle type for use *
c     * with the automatical parameter assignment features             *
c     *                                                                *
c     * originally written 04-05-2006 by M.G. Martin                   *
c     * last modified 04-05-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype
      character*(FFNAMELEN) cvalue
c     --- local variables
      integer ittemp
c     --- storage variable
      character*(FFNAMELEN) bendff
      dimension bendff(TBENMAX)
      save bendff

      if ( flag .eq. GLB_SET ) then
         bendff(itype) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = bendff(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            bendff(ittemp) = cvalue
         enddo
      else
         write(6,*) 'BENDFF: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bendnames(flag,itype,isame,index,cvalue)
c     ******************************************************************
c     * the atom names for the triplets that assign the bending terms  *
c     *                                                                *
c     * originally written 04-05-2006 by M.G. Martin                   *
c     * last modified 04-05-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,isame,index
      character*(FFNAMELEN) cvalue
c     --- local variables
      integer ittemp,istemp,iitemp
c     --- storage variable
      character*(FFNAMELEN) bendnames
      dimension bendnames(TBENMAX,MAXSAMEB,MAXBNAME)
      save bendnames

      if ( flag .eq. GLB_SET ) then
         bendnames(itype,isame,index) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = bendnames(itype,isame,index)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            do istemp = 1,isame
               do iitemp = 1,index
                  bendnames(ittemp,istemp,iitemp) = cvalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'BENDNAMES: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_benprefact(flag,imolty,iunit,ibend,dvalue)
c     ******************************************************************
c     * the angle prefactor for each angle in the molecule             *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibend
      double precision dvalue
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      double precision benprefact
      dimension benprefact(NTMAX,NUMAX,MAXBEND)
      save benprefact

      if ( flag .eq. GLB_SET ) then
         benprefact(imolty,iunit,ibend) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = benprefact(imolty,iunit,ibend)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibend
                  benprefact(imtemp,iutemp,ibtemp) = dvalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'BENPREFACT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_beta(flag,dvalue)
c     ******************************************************************
c     * the inverse temperature                                        *
c     *                                                                *
c     * originally written 11-02-2007 by M.G. Martin                   *
c     * last modified 11-02-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      double precision dvalue
c     --- storage variable
      double precision beta
      save beta

      if ( flag .eq. GLB_GET ) then
         dvalue = beta
      elseif ( flag .eq. GLB_SET ) then
         beta = dvalue
      else
         write(6,*) 'BETA: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_biorder(flag,itype,cvalue)
c     ******************************************************************
c     * the bond order for each bond increment type for matching in    *
c     * assemble                                                       *
c     *                                                                *
c     * originally written 12-20-2006 by M.G. Martin                   *
c     * last modified 12-20-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype
      character*(*) cvalue
c     --- local variables
      integer ittemp
c     --- information variables
      character*(FFNAMELEN) biorder
      dimension biorder(TVIBMAX)
      save biorder

      if ( flag .eq. GLB_SET ) then
         biorder(itype) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = biorder(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            biorder(ittemp) = cvalue
         enddo
      else
         write(6,*) 'BIORDER: invalid global flag',flag
         stop
      endif

      return
      end


      subroutine twh_blocksize(flag,ivalue)
c     ******************************************************************
c     * size of the blocks for the block averaging                     *
c     *                                                                *
c     * originally written 12-09-2007 by M.G. Martin                   *
c     * last modified 12-09-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- information variables
      integer blocksize
      save blocksize

      if ( flag .eq. GLB_GET ) then
         ivalue = blocksize
      elseif ( flag .eq. GLB_SET ) then
         blocksize = ivalue
      else
         write(6,*) 'BLOCKSIZE: invalid global flag',flag
         stop
      endif

      return
      end


      subroutine twh_bnavb1(flag,imolty,jmolty,ibox,imove,icount)
c     ******************************************************************
c     * running count of attempted aggregation-volume-bias type 1 moves*
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,jmolty,ibox,imove,icount
c     --- local variables
      integer imtemp,intemp,ibtemp,ivtemp
c     --- storage variable
      integer bnavb1
      dimension bnavb1(NTMAX,NTMAX,MAXBOX,-2:3)
      save bnavb1

      if ( flag .eq. GLB_INCR ) then
         bnavb1(imolty,jmolty,ibox,imove) 
     &        = bnavb1(imolty,jmolty,ibox,imove) + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bnavb1(imolty,jmolty,ibox,imove)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do intemp = 1,jmolty
               do ibtemp = 1,ibox
                  do ivtemp = -2,imove
                     bnavb1(imtemp,intemp,ibtemp,ivtemp) = icount
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'BNAVB1: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bnavb2(flag,imolty,jmolty,ibox,imove,icount)
c     ******************************************************************
c     * running count of attempted aggregation-volume-bias type 2 moves*
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,jmolty,ibox,imove,icount
c     --- local variables
      integer imtemp,intemp,ibtemp,ivtemp
c     --- storage variable
      integer bnavb2
      dimension bnavb2(NTMAX,NTMAX,MAXBOX,-2:3)
      save bnavb2

      if ( flag .eq. GLB_INCR ) then
         bnavb2(imolty,jmolty,ibox,imove) 
     &        = bnavb2(imolty,jmolty,ibox,imove) + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bnavb2(imolty,jmolty,ibox,imove)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do intemp = 1,jmolty
               do ibtemp = 1,ibox
                  do ivtemp = -2,imove
                     bnavb2(imtemp,intemp,ibtemp,ivtemp) = icount
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'BNAVB2: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bnavb3(flag,imolty,jmolty,ibox,icount)
c     ******************************************************************
c     * running count of attempted aggregation-volume-bias type 3 moves*
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,jmolty,ibox,icount
c     --- local variables
      integer imtemp,intemp,ibtemp
c     --- storage variable
      integer bnavb3
      dimension bnavb3(NTMAX,NTMAX,MAXBOX)
      save bnavb3

      if ( flag .eq. GLB_INCR ) then
         bnavb3(imolty,jmolty,ibox) = bnavb3(imolty,jmolty,ibox) 
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bnavb3(imolty,jmolty,ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do intemp = 1,jmolty
               do ibtemp = 1,ibox
                  bnavb3(imtemp,intemp,ibtemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BNAVB3: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bnback(flag,imolty,iunit,itor,ibox,icount)
c     ******************************************************************
c     * running count of attempted backbone configurational-bias       *
c     * regrowth moves since the last update                           *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibox,itor,icount
c     --- local variables
      integer imtemp,iutemp,ittemp,ibtemp
c     --- storage variable
      integer bnback
      dimension bnback(NTMAX,NUMAX,MAXTOR,MAXBOX)
      save bnback

      if ( flag .eq. GLB_INCR ) then
         bnback(imolty,iunit,itor,ibox) = bnback(imolty,iunit,itor,ibox)
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bnback(imolty,iunit,itor,ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ittemp = 1,itor
                  do ibtemp = 1,ibox
                     bnback(imtemp,iutemp,ittemp,ibtemp) = 0
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'BNBACK: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bncb(flag,imolty,iunit,ibox,icount)
c     ******************************************************************
c     * running count of attempted configurational-bias regrowth moves *
c     * since the last update                                          *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibox,icount
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer bncb
      dimension bncb(NTMAX,NUMAX,MAXBOX)
      save bncb

      if ( flag .eq. GLB_INCR ) then
         bncb(imolty,iunit,ibox) = bncb(imolty,iunit,ibox) + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bncb(imolty,iunit,ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibox
                  bncb(imtemp,iutemp,ibtemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BNCB: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bncbside(flag,imolty,iamino,ibox,icount)
c     ******************************************************************
c     * running count of attempted configurational-bias side chain     *
c     * regrowth moves                                                 *
c     *                                                                *
c     * originally written 06-13-2006 by M.G. Martin                   *
c     * last modified 06-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iamino,ibox,icount
c     --- local variables
      integer imtemp,iatemp,ibtemp
c     --- storage variable
      integer bncbside
      dimension bncbside(NTMAX,NUMAX,MAXBOX)
      save bncbside

      if ( flag .eq. GLB_INCR ) then
         bncbside(imolty,iamino,ibox) = bncbside(imolty,iamino,ibox)
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bncbside(imolty,iamino,ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iatemp = 1,iamino
               do ibtemp = 1,ibox
                  bncbside(imtemp,iatemp,ibtemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BNCBSIDE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bncbswap(flag,imolty,inbox,outbox,icount)
c     ******************************************************************
c     * running count of attempted configurational-bias swap moves     *
c     * and configurational-bias grand canonical insert/delete moves   *
c     * since the last update                                          *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,inbox,outbox,icount
c     --- local variables
      integer imtemp,iouttemp,iintemp
c     --- storage variable
      integer bncbswap
      dimension bncbswap(NTMAX,0:MAXBOX,0:MAXBOX)
      save bncbswap

      if ( flag .eq. GLB_INCR ) then
         bncbswap(imolty,inbox,outbox) = bncbswap(imolty,inbox,outbox)
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bncbswap(imolty,inbox,outbox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iintemp = 0,inbox
               do iouttemp = 0,outbox
                  bncbswap(imtemp,iintemp,iouttemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BNCBSWAP: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bnconrot(flag,imolty,iunit,ibox,icount)
c     ******************************************************************
c     * running count of attempted concerted rotation moves since the  *
c     * last update                                                    *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibox,icount
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer bnconrot
      dimension bnconrot(NTMAX,NUMAX,MAXBOX)
      save bnconrot

      if ( flag .eq. GLB_INCR ) then
         bnconrot(imolty,iunit,ibox) = bnconrot(imolty,iunit,ibox)
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bnconrot(imolty,iunit,ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibox
                  bnconrot(imtemp,iutemp,ibtemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BNCONROT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bncrback(flag,imolty,iunit,ibox,icount)
c     ******************************************************************
c     * running count of attempted backbone concerted rotation moves   *
c     * since the last update                                          *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibox,icount
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer bncrback
      dimension bncrback(NTMAX,NUMAX,MAXBOX)
      save bncrback

      if ( flag .eq. GLB_INCR ) then
         bncrback(imolty,iunit,ibox) = bncrback(imolty,iunit,ibox)
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bncrback(imolty,iunit,ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibox
                  bncrback(imtemp,iutemp,ibtemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BNCRBACK: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bnpivbin(flag,ibox,ibin,icount)
c     ******************************************************************
c     * running count of attempted pivot moves in each angle region bin*
c     *                                                                *
c     * originally written 03-21-2006 by M.G. Martin                   *
c     * last modified 03-21-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox,ibin,icount
c     --- local variables
      integer iboxtemp,ibintemp
c     --- storage variable
      integer bnpivbin
      dimension bnpivbin(MAXBOX,MAXPIVOTBIN)
      save bnpivbin

      if ( flag .eq. GLB_INCR ) then
         bnpivbin(ibox,ibin) = bnpivbin(ibox,ibin) + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bnpivbin(ibox,ibin)
      elseif ( flag .eq. GLB_INIT ) then
         do iboxtemp = 1,ibox
            do ibintemp = 1,ibin
               bnpivbin(iboxtemp,ibintemp) = icount
            enddo
         enddo
      else
         write(6,*) 'BNPIVBIN: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bnpivot(flag,imolty,iunit,ivib,ibox,icount)
c     ******************************************************************
c     * running count of attempted pivot moves                         *
c     *                                                                *
c     * originally written 03-21-2006 by M.G. Martin                   *
c     * last modified 03-21-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ivib,ibox,icount
c     --- local variables
      integer itemp,iutemp,ivtemp,ibtemp
c     --- storage variable
      integer bnpivot
      dimension bnpivot(NTMAX,NUMAX,NNBOND,MAXBOX)
      save bnpivot

      if ( flag .eq. GLB_INCR ) then
         bnpivot(imolty,iunit,ivib,ibox) 
     &        = bnpivot(imolty,iunit,ivib,ibox) + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bnpivot(imolty,iunit,ivib,ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do itemp = 1,imolty
            do iutemp = 1,iunit
               do ivtemp = 1,ivib
                  do ibtemp = 1,ibox
                     bnpivot(itemp,iutemp,ivtemp,ibtemp) = icount
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'BNPIVOT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bnrbswap(flag,imolty,inbox,outbox,icount)
c     ******************************************************************
c     * running count of attempted rotational-bias 2-box swap moves    *
c     * since the last update                                          *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,inbox,outbox,icount
c     --- local variables
      integer imtemp,iintemp,iouttemp
c     --- storage variable
      integer bnrbswap
      dimension bnrbswap(NTMAX,MAXBOX,MAXBOX)
      save bnrbswap

      if ( flag .eq. GLB_INCR ) then
         bnrbswap(imolty,inbox,outbox) = bnrbswap(imolty,inbox,outbox)
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bnrbswap(imolty,inbox,outbox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iintemp = 1,inbox
               do iouttemp = 1,outbox
                  bnrbswap(imtemp,iintemp,iouttemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BNRBSWAP: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bnshift(flag,imove,ibox,arbindex,icount)
c     ******************************************************************
c     * running count of attempted shift moves since the last update   *
c     *                                                                *
c     * originally written 03-21-2006 by M.G. Martin                   *
c     * last modified 03-21-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imove,ibox,arbindex,icount
c     --- local variables
      integer itemp,ibtemp,iatemp
c     --- storage variable
      integer bnshift
      dimension bnshift(2,MAXBOX,3)
      save bnshift

      if ( flag .eq. GLB_INCR ) then
         bnshift(imove,ibox,arbindex) = bnshift(imove,ibox,arbindex)
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bnshift(imove,ibox,arbindex)
      elseif ( flag .eq. GLB_INIT ) then
         do itemp = 1,imove
            do ibtemp = 1,ibox
               do iatemp = 1,arbindex
                  bnshift(itemp,ibtemp,iatemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BNSHIFT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bnvol(flag,imove,icount)
c     ******************************************************************
c     * running count of attempted volume moves since the last update  *
c     *                                                                *
c     * originally written 03-21-2006 by M.G. Martin                   *
c     * last modified 03-21-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imove,icount
c     --- local variables
      integer itemp
c     --- storage variable
      integer bnvol
      dimension bnvol(MAXPBOX)
      save bnvol

      if ( flag .eq. GLB_INCR ) then
         bnvol(imove) = bnvol(imove) + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bnvol(imove)
      elseif ( flag .eq. GLB_INIT ) then
         do itemp = 1,imove
            bnvol(itemp) = icount
         enddo
      else
         write(6,*) 'BNVOL: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bond_sdev_multiplier(flag,dvalue)
c     ******************************************************************
c     * multiplier to the computed bond standard deviation for use     *
c     * with certain CBMC options                                      *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      double precision dvalue
c     --- storage variable
      double precision bond_sdev_multiplier
      save bond_sdev_multiplier

      if ( flag .eq. GLB_GET ) then
         dvalue = bond_sdev_multiplier
      elseif ( flag .eq. GLB_SET ) then
         bond_sdev_multiplier = dvalue
      else
         write(6,*) 'BOND_SDEV_MULTIPLIER: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bondlen(flag,ibond,dvalue)
c     ******************************************************************
c     * bondlengths of atoms currently being grown in a CBMC move      *
c     *                                                                *
c     * originally written 04-12-2006 by M.G. Martin                   *
c     * last modified 04-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibond
      double precision dvalue
c     --- storage variable
      double precision bondlen
      dimension bondlen(NNBOND)
      save bondlen

      if ( flag .eq. GLB_GET ) then
         dvalue = bondlen(ibond)
      elseif ( flag .eq. GLB_SET ) then
         bondlen(ibond) = dvalue
      else
         write(6,*) 'BONDLEN: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bondpatt(flag,itype,cvalue)
c     ******************************************************************
c     * the bond pattern of each atom type in the simulation for use   *
c     * with some of the CBMC growth styles                            *
c     *                                                                *
c     * originally written 04-12-2006 by M.G. Martin                   *
c     * last modified 04-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype
      character*(*) cvalue
c     --- local variables
      integer ittemp
c     --- storage variable
      character*5 bondpatt
      dimension bondpatt(NNTYPE)
      save bondpatt

      if ( flag .eq. GLB_SET ) then
         bondpatt(itype) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = bondpatt(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            bondpatt(ittemp) = cvalue
         enddo
      else
         write(6,*) 'BONDPATT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_box_number_density(flag,ibox,dvalue)
c     ******************************************************************
c     * the initial number density in each simulation box for the      *
c     * purpose of setting up the appropriate initial volumes          *
c     *                                                                *
c     * originally written 12-09-2007 by M.G. Martin                   *
c     * last modified 12-09-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox
      double precision dvalue
c     --- storage variable
      double precision box_number_density
      dimension box_number_density(MAXBOX)
      save box_number_density

      if ( flag .eq. GLB_GET ) then
         dvalue = box_number_density(ibox)
      elseif ( flag .eq. GLB_SET ) then
         box_number_density(ibox) = dvalue
      else
         write(6,*) 'BOX_NUMBER_DENSITY: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_boxvangle(flag,ibox,dvalue)
c     ******************************************************************
c     * the classical angle and class 2 angle-angle energies that are  *
c     * not expressed as part of the dihedral in each simulation box   *
c     *                                                                *
c     * originally written 11-07-2007 by M.G. Martin                   *
c     * last modified 11-08-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox
      double precision dvalue
c     --- storage variable
      double precision boxvangle
      dimension boxvangle(MAXBOX)
      save boxvangle

      if ( flag .eq. GLB_GET ) then
         dvalue = boxvangle(ibox)
      elseif ( flag .eq. GLB_SET ) then
         boxvangle(ibox) = dvalue
      elseif ( flag .eq. GLB_INCR ) then
         boxvangle(ibox) = boxvangle(ibox) + dvalue
      else
         write(6,*) 'BOXVANGLE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_boxvolume(flag,ibox,dvalue)
c     ******************************************************************
c     * the current volume in each simulation box                      *
c     *                                                                *
c     * originally written 10-06-2008 by M.G. Martin                   *
c     * last modified 10-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox
      double precision dvalue
c     --- storage variable
      double precision boxvolume
      dimension boxvolume(MAXBOX)
      save boxvolume

      if ( flag .eq. GLB_GET ) then
         dvalue = boxvolume(ibox)
      elseif ( flag .eq. GLB_SET ) then
         boxvolume(ibox) = dvalue
      else
         write(6,*) 'BOXVOLUME: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_boxvbond(flag,ibox,dvalue)
c     ******************************************************************
c     * the classical bond energies in each simulation box             *
c     *                                                                *
c     * originally written 11-07-2007 by M.G. Martin                   *
c     * last modified 11-08-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox
      double precision dvalue
c     --- storage variable
      double precision boxvbond
      dimension boxvbond(MAXBOX)
      save boxvbond

      if ( flag .eq. GLB_GET ) then
         dvalue = boxvbond(ibox)
      elseif ( flag .eq. GLB_SET ) then
         boxvbond(ibox) = dvalue
      elseif ( flag .eq. GLB_INCR ) then
         boxvbond(ibox) = boxvbond(ibox) + dvalue
      else
         write(6,*) 'BOXVBOND: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_boxvclassic(flag,ibox,dvalue)
c     ******************************************************************
c     * the total classical energy in each simulation box              *
c     *                                                                *
c     * originally written 11-05-2007 by M.G. Martin                   *
c     * last modified 11-08-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox
      double precision dvalue
c     --- storage variable
      double precision boxvclassic
      dimension boxvclassic(MAXBOX)
      save boxvclassic

      if ( flag .eq. GLB_GET ) then
         dvalue = boxvclassic(ibox)
      elseif ( flag .eq. GLB_SET ) then
         boxvclassic(ibox) = dvalue
      elseif ( flag .eq. GLB_INCR ) then
         boxvclassic(ibox) = boxvclassic(ibox) + dvalue
      else
         write(6,*) 'BOXVCLASSIC: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_boxvelect(flag,ibox,dvalue)
c     ******************************************************************
c     * the total electrostatic energy in each simulation box          *
c     *                                                                *
c     * originally written 11-07-2007 by M.G. Martin                   *
c     * last modified 11-08-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox
      double precision dvalue
c     --- storage variable
      double precision boxvelect
      dimension boxvelect(MAXBOX)
      save boxvelect

      if ( flag .eq. GLB_GET ) then
         dvalue = boxvelect(ibox)
      elseif ( flag .eq. GLB_SET ) then
         boxvelect(ibox) = dvalue
      elseif ( flag .eq. GLB_INCR ) then
         boxvelect(ibox) = boxvelect(ibox) + dvalue
      else
         write(6,*) 'BOXVELECT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_boxvfield(flag,ibox,dvalue)
c     ******************************************************************
c     * the total external field energy in each simulation box         *
c     *                                                                *
c     * originally written 11-07-2007 by M.G. Martin                   *
c     * last modified 11-08-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox
      double precision dvalue
c     --- storage variable
      double precision boxvfield
      dimension boxvfield(MAXBOX)
      save boxvfield

      if ( flag .eq. GLB_GET ) then
         dvalue = boxvfield(ibox)
      elseif ( flag .eq. GLB_SET ) then
         boxvfield(ibox) = dvalue
      elseif ( flag .eq. GLB_INCR ) then
         boxvfield(ibox) = boxvfield(ibox) + dvalue
      else
         write(6,*) 'BOXVFIELD: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_boxvinter(flag,ibox,dvalue)
c     ******************************************************************
c     * the classical intermolecular energy (non-coulombic) in each    *
c     * simulation box                                                 *
c     *                                                                *
c     * originally written 11-07-2007 by M.G. Martin                   *
c     * last modified 11-08-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox
      double precision dvalue
c     --- storage variable
      double precision boxvinter
      dimension boxvinter(MAXBOX)
      save boxvinter

      if ( flag .eq. GLB_GET ) then
         dvalue = boxvinter(ibox)
      elseif ( flag .eq. GLB_SET ) then
         boxvinter(ibox) = dvalue
      elseif ( flag .eq. GLB_INCR ) then
         boxvinter(ibox) = boxvinter(ibox) + dvalue
      else
         write(6,*) 'BOXVINTER: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_boxvintra(flag,ibox,dvalue)
c     ******************************************************************
c     * the classical intramolecular energy (non-coulombic) in each    *
c     * simulation box                                                 *
c     *                                                                *
c     * originally written 11-07-2007 by M.G. Martin                   *
c     * last modified 11-08-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox
      double precision dvalue
c     --- storage variable
      double precision boxvintra
      dimension boxvintra(MAXBOX)
      save boxvintra

      if ( flag .eq. GLB_GET ) then
         dvalue = boxvintra(ibox)
      elseif ( flag .eq. GLB_SET ) then
         boxvintra(ibox) = dvalue
      elseif ( flag .eq. GLB_INCR ) then
         boxvintra(ibox) = boxvintra(ibox) + dvalue
      else
         write(6,*) 'BOXVINTRA: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_boxvquantum(flag,ibox,dvalue)
c     ******************************************************************
c     * the total quantum energy in each simulation box                *
c     *                                                                *
c     * originally written 11-08-2007 by M.G. Martin                   *
c     * last modified 11-08-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox
      double precision dvalue
c     --- storage variable
      double precision boxvquantum
      dimension boxvquantum(MAXBOX)
      save boxvquantum

      if ( flag .eq. GLB_GET ) then
         dvalue = boxvquantum(ibox)
      elseif ( flag .eq. GLB_SET ) then
         boxvquantum(ibox) = dvalue
      elseif ( flag .eq. GLB_INCR ) then
         boxvquantum(ibox) = boxvquantum(ibox) + dvalue
      else
         write(6,*) 'BOXVQUANTUM: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_boxvsolvation(flag,ibox,dvalue)
c     ******************************************************************
c     * the total solvation energy in each simulation box              *
c     *                                                                *
c     * originally written 11-08-2007 by M.G. Martin                   *
c     * last modified 11-08-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox
      double precision dvalue
c     --- storage variable
      double precision boxvsolvation
      dimension boxvsolvation(MAXBOX)
      save boxvsolvation

      if ( flag .eq. GLB_GET ) then
         dvalue = boxvsolvation(ibox)
      elseif ( flag .eq. GLB_SET ) then
         boxvsolvation(ibox) = dvalue
      elseif ( flag .eq. GLB_INCR ) then
         boxvsolvation(ibox) = boxvsolvation(ibox) + dvalue
      else
         write(6,*) 'BOXVSOLVATION: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_boxvtail(flag,ibox,dvalue)
c     ******************************************************************
c     * the classical intermolecular tail correction energy in each    *
c     * simulation box                                                 *
c     *                                                                *
c     * originally written 11-07-2007 by M.G. Martin                   *
c     * last modified 11-08-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox
      double precision dvalue
c     --- storage variable
      double precision boxvtail
      dimension boxvtail(MAXBOX)
      save boxvtail

      if ( flag .eq. GLB_GET ) then
         dvalue = boxvtail(ibox)
      elseif ( flag .eq. GLB_SET ) then
         boxvtail(ibox) = dvalue
      elseif ( flag .eq. GLB_INCR ) then
         boxvtail(ibox) = boxvtail(ibox) + dvalue
      else
         write(6,*) 'BOXVTAIL: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_boxvtorsion(flag,ibox,dvalue)
c     ******************************************************************
c     * the classical dihedral/torsion energies in each box            *
c     *                                                                *
c     * originally written 11-07-2007 by M.G. Martin                   *
c     * last modified 11-08-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox
      double precision dvalue
c     --- storage variable
      double precision boxvtorsion
      dimension boxvtorsion(MAXBOX)
      save boxvtorsion

      if ( flag .eq. GLB_GET ) then
         dvalue = boxvtorsion(ibox)
      elseif ( flag .eq. GLB_SET ) then
         boxvtorsion(ibox) = dvalue
      elseif ( flag .eq. GLB_INCR ) then
         boxvtorsion(ibox) = boxvtorsion(ibox) + dvalue
      else
         write(6,*) 'BOXVTORSION: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bsback(flag,imolty,iunit,itor,ibox,icount)
c     ******************************************************************
c     * running count of successfully grown backbone configurational   *
c     * bias regrowth moves since the last update                      *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibox,itor,icount
c     --- local variables
      integer imtemp,iutemp,ittemp,ibtemp
c     --- storage variable
      integer bsback
      dimension bsback(NTMAX,NUMAX,MAXTOR,MAXBOX)
      save bsback

      if ( flag .eq. GLB_INCR ) then
         bsback(imolty,iunit,itor,ibox) = bsback(imolty,iunit,itor,ibox)
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bsback(imolty,iunit,itor,ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ittemp = 1,itor
                  do ibtemp = 1,ibox
                     bsback(imtemp,iutemp,ittemp,ibtemp) = 0
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'BSBACK: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bscb(flag,imolty,iunit,ibox,icount)
c     ******************************************************************
c     * running count of successfully grown configurational-bias       *
c     * regrowth moves since the last update                           *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibox,icount
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer bscb
      dimension bscb(NTMAX,NUMAX,MAXBOX)
      save bscb

      if ( flag .eq. GLB_INCR ) then
         bscb(imolty,iunit,ibox) = bscb(imolty,iunit,ibox) + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bscb(imolty,iunit,ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibox
                  bscb(imtemp,iutemp,ibtemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BSCB: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bscbside(flag,imolty,iamino,ibox,icount)
c     ******************************************************************
c     * running count of successfully grown configurational-bias       *
c     * side-chain regrowth moves                                      *
c     *                                                                *
c     * originally written 06-14-2006 by M.G. Martin                   *
c     * last modified 06-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iamino,ibox,icount
c     --- local variables
      integer imtemp,iatemp,ibtemp
c     --- storage variable
      integer bscbside
      dimension bscbside(NTMAX,NUMAX,MAXBOX)
      save bscbside

      if ( flag .eq. GLB_INCR ) then
         bscbside(imolty,iamino,ibox) = bscbside(imolty,iamino,ibox)
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bscbside(imolty,iamino,ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iatemp = 1,iamino
               do ibtemp = 1,ibox
                  bscbside(imtemp,iatemp,ibtemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BSCBSIDE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_bscbswap(flag,imolty,inbox,outbox,icount)
c     ******************************************************************
c     * running count of successfully grown configurational-bias       *
c     * 2-box swap or grand canoncial insertion moves                  *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-23-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,inbox,outbox,icount
c     --- local variables
      integer imtemp,iintemp,iouttemp
c     --- storage variable
      integer bscbswap
      dimension bscbswap(NTMAX,0:MAXBOX,0:MAXBOX)
      save bscbswap

      if ( flag .eq. GLB_INCR ) then
         bscbswap(imolty,inbox,outbox) = bscbswap(imolty,inbox,outbox) 
     &        + icount
      elseif ( flag .eq. GLB_GET ) then
         icount = bscbswap(imolty,inbox,outbox)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iintemp = 0,inbox
               do iouttemp = 0,outbox
                  bscbswap(imtemp,iintemp,iouttemp) = icount
               enddo
            enddo
         enddo
      else
         write(6,*) 'BSCBSWAP: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_c_element(flag,ivalue)
c     ******************************************************************
c     * TMMC collection matrix element                                 *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer c_element
      save c_element

      if ( flag .eq. GLB_GET ) then
         ivalue = c_element
      elseif ( flag .eq. GLB_SET ) then
         c_element = ivalue
      else
         write(6,*) 'C_ELEMENT unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_c_matrix_freq(flag,ivalue)
c     ******************************************************************
c     * TMMC frequency to write c to file                              *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer c_matrix_freq
      save c_matrix_freq

      if ( flag .eq. GLB_GET ) then
         ivalue = c_matrix_freq
      elseif ( flag .eq. GLB_SET ) then
         c_matrix_freq = ivalue
      else
         write(6,*) 'C_MATRIX_FREQ unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_caname(flag,itype,cvalue)
c     ******************************************************************
c     * the name to use when matching angle types for each type of     *
c     * nonbonded interaction in the system                            *
c     *                                                                *
c     * originally written 04-12-2006 by M.G. Martin                   *
c     * last modified 04-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype
      character*(*) cvalue
c     --- local variables
      integer ittemp
c     --- storage variable
      character*(FFNAMELEN) caname
      dimension caname(NNTYPE)
      save caname

      if ( flag .eq. GLB_GET ) then
         cvalue = caname(itype)
      elseif ( flag .eq. GLB_SET ) then
         caname(itype) = cvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            caname(ittemp) = cvalue
         enddo
      else
         write(6,*) 'CANAME unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbbonddist(flag,imolty,iunit,ivib,ivalue)
c     ******************************************************************
c     * distribution to use when generating each bond in a cbmc move   *
c     *                                                                *
c     * originally written 07-25-2010 by M.G. Martin                   *
c     * last modified 07-16-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ivib,ivalue
c     --- local variables
      integer imtemp,iutemp,ivtemp
c     --- information variables
      integer cbbonddist
      dimension cbbonddist(NTMAX,NUMAX,NNBOND)
      save cbbonddist

      if ( flag .eq. GLB_SET ) then
         cbbonddist(imolty,iunit,ivib) = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = cbbonddist(imolty,iunit,ivib)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ivtemp = 1,ivib
                  cbbonddist(imtemp,iutemp,ivtemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'CBBONDDIST: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_cbmc_analysis(flag,cvalue)
c     ******************************************************************
c     * character string describing the CBMC analysis level            *
c     *                                                                *
c     * originally written 07-16-2014 by M.G. Martin                   *
c     * last modified 07-16-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*30 cbmc_analysis
      save cbmc_analysis

      if ( flag .eq. GLB_GET ) then
         cvalue = cbmc_analysis
      elseif ( flag .eq. GLB_SET ) then
         cbmc_analysis = cvalue
      else
         write(6,*) 'CBMC_ANALYSIS unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmc_atomname_list(flag,imolty,iunit,ivalue)
c     ******************************************************************
c     * atomlist used for special CBMC moves                           *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ivalue
c     --- storage variable
      integer cbmc_atomname_list
      dimension cbmc_atomname_list(NTMAX,NUMAX)
      save cbmc_atomname_list

      if ( flag .eq. GLB_GET ) then
         ivalue = cbmc_atomname_list(imolty,iunit)
      elseif ( flag .eq. GLB_SET ) then
         cbmc_atomname_list(imolty,iunit) = ivalue
      else
         write(6,*) 'CBMC_ATOMNAME_LIST unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmc_atomname_tot(flag,imolty,ivalue)
c     ******************************************************************
c     * atomname used for special CBMC moves                           *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,ivalue
c     --- storage variable
      integer cbmc_atomname_tot
      dimension cbmc_atomname_tot(NTMAX)
      save cbmc_atomname_tot

      if ( flag .eq. GLB_GET ) then
         ivalue = cbmc_atomname_tot(imolty)
      elseif ( flag .eq. GLB_SET ) then
         cbmc_atomname_tot(imolty) = ivalue
      else
         write(6,*) 'CBMC_ATOMNAME_TOT unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmc_bend_generation(flag,cvalue)
c     ******************************************************************
c     * the trial distribution to use when generating trials for the   *
c     * angle terms in a CBMC growth                                   *
c     *                                                                *
c     * originally written 10-11-2008 by M.G. Martin                   *
c     * last modified 10-11-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*30 cbmc_bend_generation
      save cbmc_bend_generation

      if ( flag .eq. GLB_GET ) then
         cvalue = cbmc_bend_generation
      elseif ( flag .eq. GLB_SET ) then
         cbmc_bend_generation = cvalue
      else
         write(6,*) 'CBMC_BEND_GENERATION unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmc_bendgen_code(flag,ivalue)
c     ******************************************************************
c     * integer code for the character string cbmc_bend_generation in  *
c     * order to speed up matching algorithms                          *
c     *                                                                *
c     * originally written 10-14-2008 by M.G. Martin                   *
c     * last modified 10-14-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer cbmc_bendgen_code
      save cbmc_bendgen_code

      if ( flag .eq. GLB_GET ) then
         ivalue = cbmc_bendgen_code
      elseif ( flag .eq. GLB_SET ) then
         cbmc_bendgen_code = ivalue
      else
         write(6,*) 'CBMC_BENDGEN_CODE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmc_bond_generation(flag,cvalue)
c     ******************************************************************
c     * the trial distribution to use when generating trials for the   *
c     * bond length terms in a CBMC growth                             *
c     *                                                                *
c     * originally written 10-11-2008 by M.G. Martin                   *
c     * last modified 10-11-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*30 cbmc_bond_generation
      save cbmc_bond_generation

      if ( flag .eq. GLB_GET ) then
         cvalue = cbmc_bond_generation
      elseif ( flag .eq. GLB_SET ) then
         cbmc_bond_generation = cvalue
      else
         write(6,*) 'CBMC_BOND_GENERATION unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmc_bondgen_code(flag,ivalue)
c     ******************************************************************
c     * integer code for the cbmc_bond_generation character string used*
c     * to improve matching speeds                                     *
c     *                                                                *
c     * originally written 10-14-2008 by M.G. Martin                   *
c     * last modified 10-14-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer cbmc_bondgen_code
      save cbmc_bondgen_code

      if ( flag .eq. GLB_GET ) then
         ivalue = cbmc_bondgen_code
      elseif ( flag .eq. GLB_SET ) then
         cbmc_bondgen_code = ivalue
      else
         write(6,*) 'CBMC_BONDGEN_CODE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmc_dihedral_generation(flag,cvalue)
c     ******************************************************************
c     * the trial distribution to use when generating trials for the   *
c     * dihedral terms in a CBMC growth                                *
c     *                                                                *
c     * originally written 10-11-2008 by M.G. Martin                   *
c     * last modified 10-11-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*30 cbmc_dihedral_generation
      save cbmc_dihedral_generation

      if ( flag .eq. GLB_GET ) then
         cvalue = cbmc_dihedral_generation
      elseif ( flag .eq. GLB_SET ) then
         cbmc_dihedral_generation = cvalue
      else
         write(6,*) 'CBMC_DIHEDRAL_GENERATION unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmc_dihedgen_code(flag,ivalue)
c     ******************************************************************
c     * integer code to improve match speed for the                    *
c     * cbmc_dihedral_generation variable                              *
c     *                                                                *
c     * originally written 10-14-2008 by M.G. Martin                   *
c     * last modified 10-14-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer cbmc_dihedgen_code
      save cbmc_dihedgen_code

      if ( flag .eq. GLB_GET ) then
         ivalue = cbmc_dihedgen_code
      elseif ( flag .eq. GLB_SET ) then
         cbmc_dihedgen_code = ivalue
      else
         write(6,*) 'CBMC_DIHEDGEN_CODE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmc_fit_strategy(flag,cvalue)
c     ******************************************************************
c     * character string describing the CBMC strategy for the lower    *
c     * level terms to use when fitting higher level terms.            *
c     * e.g. what bond length to use when fitting angle distributions  *
c     *                                                                *
c     * originally written 07-19-2014 by M.G. Martin                   *
c     * last modified 07-19-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*30 cbmc_fit_strategy
      save cbmc_fit_strategy

      if ( flag .eq. GLB_GET ) then
         cvalue = cbmc_fit_strategy
      elseif ( flag .eq. GLB_SET ) then
         cbmc_fit_strategy = cvalue
      else
         write(6,*) 'CBMC_FIT_STRATEGY unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmc_nb_one_gencode(flag,ibox,itype,ivalue)
c     ******************************************************************
c     * the integer code for the cbmc_nb_one_generation variable       *
c     *                                                                *
c     * originally written 10-30-2008 by M.G. Martin                   *
c     * last modified 10-30-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox,itype,ivalue
c     --- local variables
      integer ibtemp,ittemp
c     --- storage variable
      integer cbmc_nb_one_gencode
      dimension cbmc_nb_one_gencode(MAXBOX,NNTYPE)
      save cbmc_nb_one_gencode

      if ( flag .eq. GLB_GET ) then
         ivalue = cbmc_nb_one_gencode(ibox,itype)
      elseif ( flag .eq. GLB_SET ) then
         cbmc_nb_one_gencode(ibox,itype) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 1,ibox
            do ittemp = 1,itype
               cbmc_nb_one_gencode(ibtemp,ittemp) = ivalue
            enddo
         enddo
      else
         write(6,*) 'CBMC_NB_ONE_GENCODE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmc_nb_one_generation(flag,ibox,itype,cvalue)
c     ******************************************************************
c     * the trial distribution to use when generation trials for the   *
c     * first atom inserted into a simulation box for each molecule    *
c     * type and simulation box in the simulation                      *
c     *                                                                *
c     * originally written 06-07-2006 by M.G. Martin                   *
c     * last modified 06-07-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox,itype
      character*(*) cvalue
c     --- local variables
      integer ibtemp,ittemp
c     --- storage variable
      character*30 cbmc_nb_one_generation
      dimension cbmc_nb_one_generation(MAXBOX,NNTYPE)
      save cbmc_nb_one_generation

      if ( flag .eq. GLB_GET ) then
         cvalue = cbmc_nb_one_generation(ibox,itype)
      elseif ( flag .eq. GLB_SET ) then
         cbmc_nb_one_generation(ibox,itype) = cvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 1,ibox
            do ittemp = 1,itype
               cbmc_nb_one_generation(ibtemp,ittemp) = cvalue
            enddo
         enddo
      else
         write(6,*) 'CBMC_NB_ONE_GENERATION unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbname(flag,itype,cvalue)
c     ******************************************************************
c     * the name to use when matching bond types for each type of      *
c     * nonbonded interaction in the system                            *
c     *                                                                *
c     * originally written 04-12-2006 by M.G. Martin                   *
c     * last modified 04-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype
      character*(*) cvalue
c     --- local variables
      integer ittemp
c     --- storage variable
      character*(FFNAMELEN) cbname
      dimension cbname(NNTYPE)
      save cbname

      if ( flag .eq. GLB_GET ) then
         cvalue = cbname(itype)
      elseif ( flag .eq. GLB_SET ) then
         cbname(itype) = cvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            cbname(ittemp) = cvalue
         enddo
      else
         write(6,*) 'CBNAME unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmc_setting_style(flag,cvalue)
c     ******************************************************************
c     * character string describing the CBMC setting style             *
c     *                                                                *
c     * originally written 10-11-2008 by M.G. Martin                   *
c     * last modified 10-11-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*30 cbmc_setting_style
      save cbmc_setting_style

      if ( flag .eq. GLB_GET ) then
         cvalue = cbmc_setting_style
      elseif ( flag .eq. GLB_SET ) then
         cbmc_setting_style = cvalue
      else
         write(6,*) 'CBMC_SETTING_STYLE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmcaalocl(flag,iatom,icount,ivalue)
c     ******************************************************************
c     * the extra angle-angle terms required from the CBMC angle part  *
c     * but not fully determined until the CBMC dihedral part          *
c     *                                                                *
c     * originally written 12-03-2007 by M.G. Martin                   *
c     * last modified 12-03-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,iatom,icount,ivalue
c     --- storage variable
      integer cbmcaalocl
      dimension cbmcaalocl(0:2,MAXAA)
      save cbmcaalocl

      if ( flag .eq. GLB_GET ) then
         ivalue = cbmcaalocl(iatom,icount)
      elseif ( flag .eq. GLB_SET ) then
         cbmcaalocl(iatom,icount) = ivalue
      else
         write(6,*) 'CBMCAALOCL unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmcaatot(flag,ibond,icount,ivalue)
c     ******************************************************************
c     * the extra angle-angle terms required from the CBMC angle part  *
c     * but not fully determined until the CBMC dihedral part          *
c     *                                                                *
c     * originally written 12-03-2007 by M.G. Martin                   *
c     * last modified 12-03-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibond,icount,ivalue
c     --- storage variable
      integer cbmcaatot
      dimension cbmcaatot(NNBOND,3)
      save cbmcaatot

      if ( flag .eq. GLB_GET ) then
         ivalue = cbmcaatot(ibond,icount)
      elseif ( flag .eq. GLB_SET ) then
         cbmcaatot(ibond,icount) = ivalue
      elseif ( flag .eq. GLB_INCR ) then
         cbmcaatot(ibond,icount) = cbmcaatot(ibond,icount) + ivalue
      else
         write(6,*) 'CBMCAATOT unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmcaatype(flag,ibond,icount,iaa,ivalue)
c     ******************************************************************
c     * the extra angle-angle terms required from the CBMC angle part  *
c     * but not fully determined until the CBMC dihedral part          *
c     *                                                                *
c     * originally written 12-03-2007 by M.G. Martin                   *
c     * last modified 12-03-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibond,icount,iaa,ivalue
c     --- storage variable
      integer cbmcaatype
      dimension cbmcaatype(NNBOND,3,MAXAA)
      save cbmcaatype

      if ( flag .eq. GLB_GET ) then
         ivalue = cbmcaatype(ibond,icount,iaa)
      elseif ( flag .eq. GLB_SET ) then
         cbmcaatype(ibond,icount,iaa) = ivalue
      else
         write(6,*) 'CBMCAATYPE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmcbendcount(flag,imolty,iunit,ibend,jbend,ivalue)
c     ******************************************************************
c     * count of the number of peaks needed to describe a bending      *
c     * angle distribution using certain CBMC settings                 *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibend,jbend,ivalue
c     --- local variables
      integer imtemp,iutemp,ibtemp,ijtemp
c     --- storage variable
      integer cbmcbendcount
      dimension cbmcbendcount(NTMAX,NUMAX,MAXBEND,MAXBEND)
      save cbmcbendcount

      if ( flag .eq. GLB_GET ) then
         ivalue = cbmcbendcount(imolty,iunit,ibend,jbend)
      elseif ( flag .eq. GLB_SET ) then
         cbmcbendcount(imolty,iunit,ibend,jbend) = ivalue
      elseif ( flag .eq. GLB_INCR ) then
         cbmcbendcount(imolty,iunit,ibend,jbend) =
     &        cbmcbendcount(imolty,iunit,ibend,jbend) + ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibend
                  do ijtemp = 1,jbend
                     cbmcbendcount(imtemp,iutemp,ibtemp,ijtemp) = ivalue
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'CBMCBENDCOUNT unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmcbenddist(flag,imolty,iunit,ibend,jbend,ivalue)
c     ******************************************************************
c     * integer code for the arbitrary trial distribution form to use  *
c     * when generating this angle either for step A or step B         *
c     *                                                                *
c     * originally written 07-14-2014 by M.G. Martin                   *
c     * last modified 07-14-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibend,jbend,ivalue
c     --- local variables
      integer imtemp,iutemp,ibtemp,ijtemp
c     --- storage variable
      integer cbmcbenddist
      dimension cbmcbenddist(NTMAX,NUMAX,MAXBEND,MAXBEND)
      save cbmcbenddist

      if ( flag .eq. GLB_GET ) then
         ivalue = cbmcbenddist(imolty,iunit,ibend,jbend)
      elseif ( flag .eq. GLB_SET ) then
         cbmcbenddist(imolty,iunit,ibend,jbend) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibend
                  do ijtemp = 1,jbend
                     cbmcbenddist(imtemp,iutemp,ibtemp,ijtemp) 
     &                    = DIST_NULL
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'CBMCBENDDIST unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmcbendinfo(flag,imolty,iunit,ibend,jbend,ipeak
     &     ,iparam,dvalue)
c     ******************************************************************
c     * information required to create the bending angle distributions *
c     * for certain CBMC options                                       *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibend,jbend,ipeak,iparam
      double precision dvalue
c     --- local variables
      integer imtemp,iutemp,ibtemp,ijtemp,iptemp,ipartemp
c     --- storage variable
      double precision cbmcbendinfo
      dimension cbmcbendinfo(NTMAX,NUMAX,MAXBEND,MAXBEND
     &,CBMC_MAXBENDPEAKS,CBMCMAXBEND)
      save cbmcbendinfo

      if ( flag .eq. GLB_GET ) then
         dvalue = cbmcbendinfo(imolty,iunit,ibend,jbend,ipeak,iparam)
      elseif ( flag .eq. GLB_SET ) then
         cbmcbendinfo(imolty,iunit,ibend,jbend,ipeak,iparam) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibend
                  do ijtemp = 1,jbend
                     do iptemp = 1,ipeak
                        do ipartemp = 1,iparam
                           cbmcbendinfo(imtemp,iutemp,ibtemp,ijtemp
     &                          ,iptemp,ipartemp) = dvalue
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'CBMCBENDINFO unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmcbondinfo(flag,imolty,iunit,ibond,iparam,dvalue)
c     ******************************************************************
c     * information required to describe the functions used to         *
c     * generation trial bond distributions                            *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 07-18-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibond,iparam
      double precision dvalue
c     --- local variables
      integer imtemp,iutemp,ibtemp,iptemp
c     --- storage variable
      double precision cbmcbondinfo
      dimension cbmcbondinfo(NTMAX,NUMAX,NNBOND,0:CBMCMAXBOND)
      save cbmcbondinfo

      if ( flag .eq. GLB_GET ) then
         dvalue = cbmcbondinfo(imolty,iunit,ibond,iparam)
      elseif ( flag .eq. GLB_SET ) then
         cbmcbondinfo(imolty,iunit,ibond,iparam) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibond
                  do iptemp = 0,iparam
                     cbmcbondinfo(imtemp,iutemp,ibtemp,iptemp) = dvalue
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'CBMCBONDINFO unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmcimproper(flag,iatom,icount,ivalue)
c     ******************************************************************
c     * the extra improper torsions required from the CBMC angle part  *
c     * but not fully determined until the CBMC dihedral part          *
c     *                                                                *
c     * originally written 12-02-2007 by M.G. Martin                   *
c     * last modified 12-02-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,iatom,icount,ivalue
c     --- storage variable
      integer cbmcimproper
      dimension cbmcimproper(2:4,MAXIMPROP)
      save cbmcimproper

      if ( flag .eq. GLB_GET ) then
         ivalue = cbmcimproper(iatom,icount)
      elseif ( flag .eq. GLB_SET ) then
         cbmcimproper(iatom,icount) = ivalue
      else
         write(6,*) 'CBMCIMPROPER unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmcrimptype(flag,icount,ivalue)
c     ******************************************************************
c     * information about improper torsions that need to be generated  *
c     * in a given step of the CBMC move                               *
c     * the 0 index contains the total number of imptypes              *
c     *                                                                *
c     * originally written 12-03-2007 by M.G. Martin                   *
c     * last modified 12-03-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,icount,ivalue
c     --- storage variable
      integer cbmcrimptype
      dimension cbmcrimptype(0:MAXIMPROP)
      save cbmcrimptype

      if ( flag .eq. GLB_GET ) then
         ivalue = cbmcrimptype(icount)
      elseif ( flag .eq. GLB_SET ) then
         cbmcrimptype(icount) = ivalue
      else
         write(6,*) 'CBMCRIMPTYPE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmctorcount(flag,imolty,iunit,itor,ivalue)
c     ******************************************************************
c     * the number of peaks in a multiple gaussian distribution for    *
c     * each dihedral in the simulation                                *
c     *                                                                *
c     * originally written 05-18-2006 by M.G. Martin                   *
c     * last modified 05-18-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,itor,ivalue
c     --- local variables
      integer imtemp,iutemp,ittemp
c     --- storage variable
      integer cbmctorcount
      dimension cbmctorcount(NTMAX,NUMAX,MAXTOR)
      save cbmctorcount

      if ( flag .eq. GLB_GET ) then
         ivalue = cbmctorcount(imolty,iunit,itor)
      elseif ( flag .eq. GLB_SET ) then
         cbmctorcount(imolty,iunit,itor) = ivalue
      elseif ( flag .eq. GLB_INCR ) then
         cbmctorcount(imolty,iunit,itor) 
     &        = cbmctorcount(imolty,iunit,itor) + ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ittemp = 1,itor
                  cbmctorcount(imtemp,iutemp,ittemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'CBMCTORCOUNT unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmctorinfo(flag,imolty,iunit,itor,ipeak,iparam
     &     ,dvalue)
c     ******************************************************************
c     * information used to generate torsion distribution trials for   *
c     * certain CBMC options                                           *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,itor,ipeak,iparam
      double precision dvalue
c     --- local variables
      integer imtemp,iutemp,ittemp,iptemp,ipartemp
c     --- storage variable
      double precision cbmctorinfo
      dimension cbmctorinfo(NTMAX,NUMAX,MAXTOR
     &,CBMC_MAXTORPEAKS,CBMCMAXTOR)
      save cbmctorinfo

      if ( flag .eq. GLB_GET ) then
         dvalue = cbmctorinfo(imolty,iunit,itor,ipeak,iparam)
      elseif ( flag .eq. GLB_SET ) then
         cbmctorinfo(imolty,iunit,itor,ipeak,iparam) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ittemp = 1,itor
                  do iptemp = 1,ipeak
                     do ipartemp = 1,iparam
                        cbmctorinfo(imtemp,iutemp,ittemp,iptemp
     &                       ,ipartemp) = dvalue
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'CBMCTORINFO unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbphi(flag,ibond,dvalue)
c     ******************************************************************
c     * phi angle values used during CBMC growth moves                 *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modified 10-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibond
      double precision dvalue
c     --- storage variable
      double precision cbphi
      dimension cbphi(NNBOND)
      save cbphi

      if ( flag .eq. GLB_GET ) then
         dvalue = cbphi(ibond)
      elseif ( flag .eq. GLB_SET ) then
         cbphi(ibond) = dvalue
      else
         write(6,*) 'CBPHI unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_charge_assignment(flag,imolty,cvalue)
c     ******************************************************************
c     * charge assignment string for each type of molecule             *
c     *                                                                *
c     * originally written 12-23-2007 by M.G. Martin                   *
c     * last modified 12-23-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty
      character*(*) cvalue
c     --- local variables
      integer imtemp
c     --- storage variable
      character*30 charge_assignment
      dimension charge_assignment(NTMAX)
      save charge_assignment

      if ( flag .eq. GLB_GET ) then
         cvalue = charge_assignment(imolty)
      elseif ( flag .eq. GLB_SET ) then
         charge_assignment(imolty) = cvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            charge_assignment(imtemp) = cvalue
         enddo
      else
         write(6,*) 'CHARGE_ASSIGNMENT unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_chempot(flag,imolty,dvalue)
c     ******************************************************************
c     * the chemical potential specified in a grand canonical ensemble *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer itemp
c     --- storage variable
      double precision chempot
      dimension chempot(NTMAX)
      save chempot

      if ( flag .eq. GLB_GET ) then
         dvalue = chempot(imolty)
      elseif ( flag .eq. GLB_SET ) then
         chempot(imolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do itemp = 1,imolty
            chempot(itemp) = dvalue
         enddo
      else
         write(6,*) 'CHEMPOT unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_chempotperstep(flag,imolty,ivalue)
c     ******************************************************************
c     * the number of additional chemical potential measurements to    *
c     * perform at the end of each Monte Carlo cycle                   *
c     *                                                                *
c     * originally written 12-22-2007 by M.G. Martin                   *
c     * last modified 10-27-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,ivalue
c     --- local variables
      integer imtemp
c     --- storage variable
      integer chempotperstep
      dimension chempotperstep(NTMAX)
      save chempotperstep

      if ( flag .eq. GLB_GET ) then
         ivalue = chempotperstep(imolty)
      elseif ( flag .eq. GLB_SET ) then
         chempotperstep(imolty) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            chempotperstep(imtemp) = ivalue
         enddo
      else
         write(6,*) 'CHEMPOTPERSTEP unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_classical_mixrule(flag,cvalue)
c     ******************************************************************
c     * the classical mixing rule expressed as a test string           *
c     *                                                                *
c     * originally written 04-15-2006 by M.G. Martin                   *
c     * last modified 04-15-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*30 classical_mixrule
      save classical_mixrule

      if ( flag .eq. GLB_GET ) then
         cvalue = classical_mixrule
      elseif ( flag .eq. GLB_SET ) then
         classical_mixrule = cvalue
      else
         write(6,*) 'CLASSICAL_MIXRULE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_classical_potential(flag,cvalue)
c     ******************************************************************
c     * the classical potential expressed as a text string             *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-26-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*30 classical_potential
      save classical_potential

      if ( flag .eq. GLB_GET ) then
         cvalue = classical_potential
      elseif ( flag .eq. GLB_SET ) then
         classical_potential = cvalue
c        --- also set the classical potential number
         call twh_set_classicalpotnum(cvalue)
      else
         write(6,*) 'CLASSICAL_POTENTIAL unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_com(flag,cordflag,ichain,xxx,yyy,zzz)
c     ******************************************************************
c     * get or sets the coordinates of the center of mass for each     *
c     * chain in the system                                            *
c     *                                                                *
c     * originally written 11-02-2007 by M.G. Martin                   *
c     * last modified 06-24-2009 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "functions.h"
c     --- variables passed to/from this subroutine
      integer flag,cordflag
      integer ichain
      double precision xxx,yyy,zzz
c     --- local variables
      integer index
      double precision vector,backup
      dimension vector(3),backup(3)

      index = 3*(ichain-1)
      if ( flag .eq. GLB_GET ) then
         if ( cordflag .eq. CRD_REAL ) then
            call twh_comfield(GLB_GET_TRIPLE,index,vector)
         elseif ( cordflag .eq. CRD_TEMP ) then
            call twh_comtempfield(GLB_GET_TRIPLE,index,vector)
         else
            goto 666
         endif
         xxx = vector(1)
         yyy = vector(2)
         zzz = vector(3)
         return
      elseif ( flag .eq. GLB_SET ) then
         vector(1) = xxx
         vector(2) = yyy
         vector(3) = zzz
         if ( cordflag .eq. CRD_REAL ) then
            call twh_comfield(GLB_SET_TRIPLE,index,vector)
            return
         elseif ( cordflag .eq. CRD_TEMP ) then
            call twh_comtempfield(GLB_SET_TRIPLE,index,vector)
            return
         endif
      elseif ( flag .eq. GLB_INCR ) then
         if ( cordflag .eq. CRD_REAL ) then
            vector(1) = xxx
            vector(2) = yyy
            vector(3) = zzz
            call twh_comfield(GLB_INCR_TRIPLE,index,vector)
            return
         endif
      elseif ( flag .eq. GLB_STORE ) then
         if ( cordflag .eq. CRD_R_TO_T ) then
c           --- store the real coordinates in the temp array
            call twh_comfield(GLB_GET_TRIPLE,index,vector)
            call twh_comtempfield(GLB_SET_TRIPLE,index,vector)
            return
         elseif ( cordflag .eq. CRD_T_TO_R ) then
c           --- store the temp coordinates in the real array
            call twh_comtempfield(GLB_GET_TRIPLE,index,vector)
            call twh_comfield(GLB_SET_TRIPLE,index,vector)
            return
         endif
      elseif ( flag .eq. GLB_EXCHANGE ) then
c        --- exchange the positions in the two arrays
         if ( (cordflag .eq. CRD_R_TO_T) .or.
     &        (cordflag .eq. CRD_T_TO_R) ) then
c           --- swap the coordinates
            call twh_comfield(GLB_GET_TRIPLE,index,vector)
            call twh_comtempfield(GLB_GET_TRIPLE,index,backup)
            call twh_comfield(GLB_SET_TRIPLE,index,backup)
            call twh_comtempfield(GLB_SET_TRIPLE,index,vector)
            return
         endif
      endif

c     --- if we got here there was an error
 666  write(6,*) 'COM: invalid combination of flag and cordflag'
      write(6,*) 'flag:',flag
      write(6,*) 'cordflag:',cordflag
      stop

      return
      end

      subroutine twh_constant_cuberoottwo(flag,dvalue)
c     ******************************************************************
c     * the constant the cube root of two computed to machine precision*
c     * and stored for easy reuse                                      *
c     *                                                                *
c     * originally written 10-30-2007 by M.G. Martin                   *
c     * last modified 10-30-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision cuberoottwo
      save cuberoottwo

      if ( flag .eq. GLB_GET ) then
         dvalue = cuberoottwo
      elseif ( flag .eq. GLB_INIT ) then
         cuberoottwo = 2.0d0**(1.0d0 / 3.0d0)
      else
         write(6,*) 'CONSTANT_CUBEROOTTWO: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_constant_halfpi(flag,dvalue)
c     ******************************************************************
c     * the constant 0.5* Pi computed to machine precision and stored  *
c     * for easy reuse                                                 *
c     *                                                                *
c     * originally written 10-30-2007 by M.G. Martin                   *
c     * last modified 10-30-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision halfpi
      save halfpi

      if ( flag .eq. GLB_GET ) then
         dvalue = halfpi
      elseif ( flag .eq. GLB_INIT ) then
         halfpi = 2.0d0*datan(1.0d0)
      else
         write(6,*) 'CONSTANT_HALF_PI: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_constant_invrootpi(flag,dvalue)
c     ******************************************************************
c     * the constant 1.0 / sqrt(Pi) computed to machine precision and  *
c     * stored for use in the Ewald self interaction                   *
c     *                                                                *
c     * originally written 11-02-2007 by M.G. Martin                   *
c     * last modified 11-02-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision invrootpi
      save invrootpi

      if ( flag .eq. GLB_GET ) then
         dvalue = invrootpi
      elseif ( flag .eq. GLB_INIT ) then
         invrootpi = 1.0d0/dsqrt(4.0d0*datan(1.0d0))
      else
         write(6,*) 'CONSTANT_INVROOTPI: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_constant_onethird(flag,dvalue)
c     ******************************************************************
c     * the constant 1.0/3.0 saved for easy reuse                      *
c     *                                                                *
c     * originally written 10-30-2007 by M.G. Martin                   *
c     * last modified 10-30-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision onethird
      save onethird

      if ( flag .eq. GLB_GET ) then
         dvalue = onethird
      elseif ( flag .eq. GLB_INIT ) then
         onethird = 1.0d0/3.0d0
      else
         write(6,*) 'CONSTANT_ONETHIRD: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_constant_pi(flag,dvalue)
c     ******************************************************************
c     * the constant Pi computed to machine precision and stored for   *
c     * easy reuse                                                     *
c     *                                                                *
c     * originally written 04-17-2006 by M.G. Martin                   *
c     * last modified 04-17-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pi
      save pi

      if ( flag .eq. GLB_INIT ) then
         pi = 4.0d0*datan(1.0d0)
      elseif ( flag .eq. GLB_GET ) then
         dvalue = pi
      else
         write(6,*) 'CONSTANT_PI: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_constant_radtodeg(flag,dvalue)
c     ******************************************************************
c     * a double containing the value to convert radians to degrees    *
c     *                                                                *
c     * originally written 10-30-2007 by M.G. Martin                   *
c     * last modified 10-30-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision radtodeg
      save radtodeg

      if ( flag .eq. GLB_GET ) then
         dvalue = radtodeg
      elseif ( flag .eq. GLB_INIT ) then
         radtodeg = 45.0d0/datan(1.0d0)
      else
         write(6,*) 'CONSTANT_RADTODEG: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_constant_twopi(flag,dvalue)
c     ******************************************************************
c     * a double containing the value of 2*pi                          *
c     *                                                                *
c     * originally written 10-30-2007 by M.G. Martin                   *
c     * last modified 10-30-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision twopi
      save twopi

      if ( flag .eq. GLB_GET ) then
         dvalue = twopi
      elseif ( flag .eq. GLB_INIT ) then
         twopi = 8.0d0*datan(1.0d0)
      else
         write(6,*) 'CONSTANT_TWOPI: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_constant_underflow(flag,dvalue)
c     ******************************************************************
c     * a double containing the value below which we are in danger of  *
c     * getting an underflow error                                     *
c     *                                                                *
c     * originally written 10-30-2007 by M.G. Martin                   *
c     * last modified 10-30-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision underflow
      save underflow

      if ( flag .eq. GLB_GET ) then
         dvalue = underflow
      elseif ( flag .eq. GLB_INIT ) then
         underflow = dlog(SMALLEST)
      else
         write(6,*) 'CONSTANT_UNDERFLOW: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_constant_qqfact(flag,dvalue)
c     ******************************************************************
c     * a double containing the value of the charge prefactor          *
c     *                                                                *
c     * originally written 10-30-2007 by M.G. Martin                   *
c     * last modified 10-30-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision qqfact
      save qqfact

      if ( flag .eq. GLB_GET ) then
         dvalue = qqfact
      elseif ( flag .eq. GLB_INIT ) then
         qqfact = CON_EETOKCAL*CON_KCALTOK
      elseif (flag .eq. GLB_SET ) then
         qqfact = dvalue
      else
         write(6,*) 'CONSTANT_QQFACT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_controlstyle(flag,cvalue)
c     ******************************************************************
c     * a character string containing the control style                *
c     *                                                                *
c     * originally written 08-10-2008 by M.G. Martin                   *
c     * last modified 08-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      character*(*) cvalue
c     --- information variables
      character*20 controlstyle
      save controlstyle

      if ( flag .eq. GLB_GET ) then
         cvalue = controlstyle
      elseif (flag .eq. GLB_SET ) then
         controlstyle = cvalue
      else
         write(6,*) 'CONTROLSTYLE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_coordinates(flag,cordflag,ichain,iunit,xxx,yyy,zzz)
c     ******************************************************************
c     * get or sets the coordinates of a single atom in the            *
c     * coordstorage array                                             *
c     *                                                                *
c     * originally written 08-04-2003 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETATOMNUM
#include "functions.h"
c     --- variables passed to/from this subroutine
      integer flag,cordflag
      integer ichain,iunit
      double precision xxx,yyy,zzz
c     --- local variables
      integer index
      double precision coordvector
      dimension coordvector(3)
c     --- storage variables
      double precision coordsingle
      dimension coordsingle(NUMAX,3)
      save coordsingle

      if ( flag .eq. GLB_SET ) then
         if ( cordflag .eq. CRD_REAL ) then
            index = twh_getatomnum(ichain,iunit)
            coordvector(1) = xxx
            coordvector(2) = yyy
            coordvector(3) = zzz
            call twh_coordstorage(GLB_SET_TRIPLE,index,coordvector)
            return
         elseif ( cordflag .eq. CRD_TEMP ) then
            index = twh_getatomnum(ichain,iunit)
            coordvector(1) = xxx
            coordvector(2) = yyy
            coordvector(3) = zzz
            call twh_coordtemp(GLB_SET_TRIPLE,index,coordvector)
            return
         elseif ( cordflag .eq. CRD_SINGLE ) then
            coordsingle(iunit,1) = xxx
            coordsingle(iunit,2) = yyy
            coordsingle(iunit,3) = zzz
            return
         endif
      elseif ( flag .eq. GLB_GET ) then
         if ( cordflag .eq. CRD_REAL ) then
            index = twh_getatomnum(ichain,iunit)
            call twh_coordstorage(GLB_GET_TRIPLE,index,coordvector)
            xxx = coordvector(1)
            yyy = coordvector(2)
            zzz = coordvector(3)
            return
         elseif ( cordflag .eq. CRD_TEMP ) then
            index = twh_getatomnum(ichain,iunit)
            call twh_coordtemp(GLB_GET_TRIPLE,index,coordvector)
            xxx = coordvector(1)
            yyy = coordvector(2)
            zzz = coordvector(3)
            return
         elseif ( cordflag .eq. CRD_SINGLE ) then
            xxx = coordsingle(iunit,1)
            yyy = coordsingle(iunit,2)
            zzz = coordsingle(iunit,3)
            return
         endif
      elseif (flag .eq. GLB_FIX) then
         if ( cordflag .eq. CRD_TEMP ) then
            index =  iunit
            call twh_coordtemp(GLB_SET,index,xxx)
            return
         endif
      elseif (flag .eq. GLB_CHK) then
         index = iunit 
         if ( cordflag .eq. CRD_TEMP ) then
            call twh_coordtemp(GLB_GET,index,xxx)
            yyy = 0.0d0
            zzz = 0.0d0
            return
         endif
      elseif ( flag .eq. GLB_UST ) then
         index = 3*iunit-2
         if ( cordflag .eq. CRD_TEMP ) then
            coordvector(1) = xxx
            coordvector(2) = yyy
            coordvector(3) = zzz
            call twh_coordtemp(GLB_SET_TRIPLE,index,coordvector)
            return
         endif
      elseif ( flag .eq. GLB_UGT ) then
         index = 3*iunit-2
         if ( cordflag .eq. CRD_TEMP ) then
            call twh_coordtemp(GLB_GET_TRIPLE,index,coordvector)
            xxx = coordvector(1)
            yyy = coordvector(2)
            zzz = coordvector(3)
            return
         endif
      endif

c     --- if we made it this far then there was an error
      write(6,*) 'COORDINATE: invalid flag and cordflag combination'
      write(6,*) 'flag:',flag
      write(6,*) 'cordflag:',cordflag
      stop

      return
      end

      subroutine twh_corp(flag,imolty,jmolty,dvalue)
c     ******************************************************************
c     * correction to the pressure for tail corrections for each pair  *
c     * of molecule types                                              *
c     *                                                                *
c     * originally written 10-03-2008 by M.G. Martin                   *
c     * last modified 10-03-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,jmolty
      double precision dvalue
c     --- storage variable
      double precision corp
      dimension corp(NTMAX,NTMAX)
      save corp 

      if ( flag .eq. GLB_GET ) then
         dvalue = corp(imolty,jmolty)
      elseif ( flag .eq. GLB_SET ) then
         corp(imolty,jmolty) = dvalue
      elseif ( flag .eq. GLB_INCR ) then
         corp(imolty,jmolty) = corp(imolty,jmolty) + dvalue
      else
         write(6,*) 'CORP unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_coru(flag,imolty,jmolty,dvalue)
c     ******************************************************************
c     * correction to the chemical potential for tail corrections for  *
c     * each pair of molecule types                                    *
c     *                                                                *
c     * originally written 10-03-2008 by M.G. Martin                   *
c     * last modified 10-03-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,jmolty
      double precision dvalue
c     --- storage variable
      double precision coru
      dimension coru(NTMAX,NTMAX)
      save coru

      if ( flag .eq. GLB_GET ) then
         dvalue = coru(imolty,jmolty)
      elseif ( flag .eq. GLB_SET ) then
         coru(imolty,jmolty) = dvalue
      elseif ( flag .eq. GLB_INCR ) then
         coru(imolty,jmolty) = coru(imolty,jmolty) + dvalue
      else
         write(6,*) 'CORU unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_coulombstyle(flag,cvalue)
c     ******************************************************************
c     * the style of coulombic interactions                            *
c     *                                                                *
c     * originally written 11-05-2007 by M.G. Martin                   *
c     * last modified 11-05-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*20 coulombstyle
      save coulombstyle

      if ( flag .eq. GLB_GET ) then
         cvalue = coulombstyle
      elseif ( flag .eq. GLB_SET ) then
         coulombstyle = cvalue
      else
         write(6,*) 'COULOMBSTYLE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmc_formulation(flag,cvalue)
c     ******************************************************************
c     * form of the CBMC algorithm used in this simulation             *
c     *                                                                *
c     * originally written 01-22-2009 by M.G. Martin                   *
c     * last modified 01-22-2009 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- local variables
      integer icode
c     --- storage variable
      character*60 cbmc_formulation
      save cbmc_formulation

      if ( flag .eq. GLB_GET ) then
         cvalue = cbmc_formulation
      elseif ( flag .eq. GLB_SET ) then
         cbmc_formulation = cvalue
         if ( cbmc_formulation .eq. 
     &        'Martin and Siepmann 1999 + Martin and Thompson 2004' 
     &        ) then
            icode = CBMC_FORM_MS1999_MT2004
         elseif ( cbmc_formulation .eq. 'Martin and Frischknecht 2006' 
     &           ) then
            icode = CBMC_FORM_MF2006
         else
            write(6,*) 'CBMC_FORMULATION: UNKNOWN CBMC_FORMULATION'
            write(6,*) cbmc_formulation
            stop
         endif
         call twh_cbmc_form_code(GLB_SET,icode)
      else
         write(6,*) 'CBMC_FORMULATION unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cbmc_form_code(flag,ivalue)
c     ******************************************************************
c     * integer code for the form of the CBMC algorithm used in this   *
c     * simulation                                                     *
c     *                                                                *
c     * originally written 01-22-2009 by M.G. Martin                   *
c     * last modified 01-22-2009 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer cbmc_form_code
      save cbmc_form_code

      if ( flag .eq. GLB_GET ) then
         ivalue = cbmc_form_code
      elseif ( flag .eq. GLB_SET ) then
         cbmc_form_code = ivalue
      else
         write(6,*) 'CBMC_FORM_CODE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cp_number(flag,ivalue)
c     ******************************************************************
c     * the classical potential expressed as an integer                *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-26-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer cp_number
      save cp_number

      if ( flag .eq. GLB_GET ) then
         ivalue = cp_number
      elseif ( flag .eq. GLB_SET ) then
         cp_number = ivalue
      else
         write(6,*) 'CP_NUMBER unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_crvalidcount(flag,imolty,ivalue)
c     ******************************************************************
c     * the total number of valid targets for a concerted rotation move*
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,ivalue
c     --- storage variable
      integer crvalidcount
      dimension crvalidcount(NTMAX)
      save crvalidcount

      if ( flag .eq. GLB_GET ) then
         ivalue = crvalidcount(imolty)
      elseif ( flag .eq. GLB_SET ) then
         crvalidcount(imolty) = ivalue
      else
         write(6,*) 'CRVALIDCOUNT unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_crvalidlist(flag,imolty,icount,ivalue)
c     ******************************************************************
c     * the list of valid targets for a concerted rotation move        *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,icount,ivalue
c     --- storage variable
      integer crvalidlist
      dimension crvalidlist(NTMAX,NUMAX)
      save crvalidlist

      if ( flag .eq. GLB_GET ) then
         ivalue = crvalidlist(imolty,icount)
      elseif ( flag .eq. GLB_SET ) then
         crvalidlist(imolty,icount) = ivalue
      else
         write(6,*) 'CRVALIDLIST unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ctname(flag,itype,cvalue)
c     ******************************************************************
c     * the name to use when matching regular torsion types for each   *
c     * type of nonbonded interaction in the system                    *
c     *                                                                *
c     * originally written 04-12-2006 by M.G. Martin                   *
c     * last modified 04-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype
      character*(*) cvalue
c     --- local variables
      integer ittemp
c     --- storage variable
      character*(FFNAMELEN) ctname
      dimension ctname(NNTYPE)
      save ctname

      if ( flag .eq. GLB_GET ) then
         cvalue = ctname(itype)
      elseif ( flag .eq. GLB_SET ) then
         ctname(itype) = cvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            ctname(ittemp) = cvalue
         enddo
      else
         write(6,*) 'CTNAME unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_cubelet(flag,xxx,yyy,zzz)
c     ******************************************************************
c     * the number of energy biasing cublets in the 3 dimensions       *
c     *                                                                *
c     * originally written 11-04-2007 by M.G. Martin                   *
c     * last modified 11-04-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,xxx,yyy,zzz
c     --- storage variable
      integer cubelet
      dimension cubelet(3)
      save cubelet

      if ( flag .eq. GLB_GET ) then
         xxx = cubelet(1)
         yyy = cubelet(2)
         zzz = cubelet(3)
      elseif ( flag .eq. GLB_SET ) then
         cubelet(1) = xxx
         cubelet(2) = yyy
         cubelet(3) = zzz
      else
         write(6,*) 'CUBELET unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_debroglie(flag,imolty,dvalue)
c     ******************************************************************
c     * the deBroglie thermal wavelenth for each molecule type         *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty
      double precision dvalue
c     --- storage variable
      double precision debroglie
      dimension debroglie(NTMAX)
      save debroglie

      if ( flag .eq. GLB_GET ) then
         dvalue = debroglie(imolty)
      elseif ( flag .eq. GLB_SET ) then
         debroglie(imolty) = dvalue
      else
         write(6,*) 'DEBROGLIE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_dihedral_ideal_fraction(flag,dvalue)
c     ******************************************************************
c     * fraction of the generation algorithm that uses the ideal       *
c     * distribution for dihedrals using certain CBMC algorithms       *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      double precision dvalue
c     --- storage variable
      double precision dihedral_ideal_fraction
      save dihedral_ideal_fraction

      if ( flag .eq. GLB_GET ) then
         dvalue = dihedral_ideal_fraction
      elseif ( flag .eq. GLB_SET ) then
         dihedral_ideal_fraction = dvalue
      else
         write(6,*) 'DIHEDRAL_IDEAL_FRACTION unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_dihedral_peak_weight_style(flag,cvalue)
c     ******************************************************************
c     * the peak generation style to use when selecting from multiple  *
c     * peaks in a dihedral distribution with CBMC                     *
c     *                                                                *
c     * originally written 05-18-2006 by M.G. Martin                   *
c     * last modified 05-18-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*30 dihedral_peak_weight_style
      save dihedral_peak_weight_style

      if ( flag .eq. GLB_GET ) then
         cvalue = dihedral_peak_weight_style
      elseif ( flag .eq. GLB_SET ) then
         dihedral_peak_weight_style = cvalue
      else
         write(6,*) 'DIHEDRAL_PEAK_WEIGHT_STYLE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_dihedral_peak_weights(flag,imolty,iunit,itor,ipeak
     &     ,dvalue)
c     ******************************************************************
c     * the probability density for each peak in a multiple gaussian   *
c     * selection for every dihedral in the simulation                 *
c     *                                                                *
c     * originally written 05-18-2006 by M.G. Martin                   *
c     * last modified 05-18-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,itor,ipeak
      double precision dvalue
c     --- local variables
      integer imtemp,iutemp,ittemp,iptemp
c     --- storage variable
      double precision dihedral_peak_weights
      dimension dihedral_peak_weights(NTMAX,NUMAX,MAXTOR
     &     ,CBMC_MAXTORPEAKS)
      save dihedral_peak_weights

      if ( flag .eq. GLB_GET ) then
         dvalue = dihedral_peak_weights(imolty,iunit,itor,ipeak)
      elseif ( flag .eq. GLB_SET ) then
         dihedral_peak_weights(imolty,iunit,itor,ipeak) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ittemp = 1,itor
                  do iptemp = 1,ipeak
                     dihedral_peak_weights(imtemp,iutemp,ittemp,iptemp)
     &                    = dvalue
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'DIHEDRAL_PEAK_WEIGHTS unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_dihedral_sdev_multiplier(flag,dvalue)
c     ******************************************************************
c     * multiplier to rescale the computed standard devations in the   *
c     * dihedral distributions for certain CBMC moves                  *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      double precision dvalue
c     --- storage variable
      double precision dihedral_sdev_multiplier
      save dihedral_sdev_multiplier

      if ( flag .eq. GLB_GET ) then
         dvalue = dihedral_sdev_multiplier
      elseif ( flag .eq. GLB_SET ) then
         dihedral_sdev_multiplier = dvalue
      else
         write(6,*) 'DIHEDRAL_SDEV_MULTIPLIER unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_distprev(flag,ibond,dvalue)
c     ******************************************************************
c     * distance to the prev atoms during a CBMC regrowth move         *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modified 10-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibond
      double precision dvalue
c     --- storage variable
      double precision distprev
      dimension distprev(NNBOND)
      save distprev

      if ( flag .eq. GLB_GET ) then
         dvalue = distprev(ibond)
      elseif ( flag .eq. GLB_SET ) then
         distprev(ibond) = dvalue
      else
         write(6,*) 'DISTPREV unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_donor_acceptor(flag,itype,cvalue)
c     ******************************************************************
c     * the donor/acceptor flag for each nonbonded type.  This is only *
c     * used with a subset of the potentials                           *
c     *                                                                *
c     * originally written 12-08-2006 by M.G. Martin                   *
c     * last modified 12-08-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype
      character*(*) cvalue
c     --- local variables
      integer ittemp
c     --- storage variable
      character*1 donor_acceptor
      dimension donor_acceptor(NNTYPE)
      save donor_acceptor

      if ( flag .eq. GLB_GET ) then
         cvalue = donor_acceptor(itype)
      elseif ( flag .eq. GLB_SET ) then
         donor_acceptor(itype) = cvalue
      elseif ( flag .eq. GLB_INIT ) then 
         do ittemp = 1,itype
            donor_acceptor(ittemp) = cvalue
         enddo
      else
         write(6,*) 'DONOR_ACCEPTOR unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_du_dlambda(flag,ibox,count,dvalue)
c     ******************************************************************
c     * change in energy versus lambda for scaled potentials           *
c     *                                                                *
c     * originally written 12-20-2007 by M.G. Martin                   *
c     * last modified 12-20-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox,count
      double precision dvalue
c     --- local variables
c     --- storage variable
      double precision du_dlambda
      dimension du_dlambda(MAXBOX,2)
      save du_dlambda

      if ( flag .eq. GLB_GET ) then
         dvalue = du_dlambda(ibox,count)
      elseif ( flag .eq. GLB_SET ) then
         du_dlambda(ibox,count) = dvalue
      elseif ( flag .eq. GLB_INCR ) then
         du_dlambda(ibox,count) = du_dlambda(ibox,count) + dvalue
      else
         write(6,*) 'DU_DLAMBDA unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_eam_dens(flag,ientry,itype,jtype,itable,dvalue)
c     ******************************************************************
c     * coefficients for the EAM densities                             *
c     *                                                                *
c     * originally written 10-02-2008 by M.G. Martin                   *
c     * last modified 10-02-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ientry,itype,jtype,itable
      double precision dvalue
c     --- local variables
      integer itemp,jtemp,ttemp,etemp
c     --- storage variable
      double precision eam_dens
      dimension eam_dens(3,MAXTABTYPE,MAXTABTYPE,MAXTABLE)
      save eam_dens

      if ( flag .eq. GLB_GET ) then
         dvalue = eam_dens(ientry,itype,jtype,itable)
      elseif ( flag .eq. GLB_SET ) then
         eam_dens(ientry,itype,jtype,itable) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do etemp = 1,ientry
            do itemp = 1,itype
               do jtemp = 1,jtype
                  do ttemp = 1,itable
                     eam_dens(etemp,itemp,jtemp,ttemp) = dvalue
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'EAM_DENS unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_eam_dens_num(flag,itype,jtype,ivalue)
c     ******************************************************************
c     * integer code for the style of the EAM density potential        *
c     *                                                                *
c     * originally written 07-21-2013 by M.G. Martin                   *
c     * last modified 07-21-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,jtype,ivalue
c     --- storage variable
      integer eam_dens_num
      dimension eam_dens_num(MAXTABTYPE,MAXTABTYPE)
      save eam_dens_num

      if ( flag .eq. GLB_GET ) then
         ivalue = eam_dens_num(itype,jtype)
      elseif ( flag .eq. GLB_SET ) then
         eam_dens_num(itype,jtype) = ivalue
      else
         write(6,*) 'EAM_DENS_NUM unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_eam_dens_style(flag,itype,jtype,cvalue)
c     ******************************************************************
c     * style of the EAM density potential                             *
c     *                                                                *
c     * originally written 10-03-2008 by M.G. Martin                   *
c     * last modified 10-03-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,jtype
      character*(*) cvalue
c     --- storage variable
      character*20 eam_dens_style
      dimension eam_dens_style(MAXTABTYPE,MAXTABTYPE)
      save eam_dens_style

      if ( flag .eq. GLB_GET ) then
         cvalue = eam_dens_style(itype,jtype)
      elseif ( flag .eq. GLB_SET ) then
         eam_dens_style(itype,jtype) = cvalue
c        --- set this eam_dens_num as well
         call twh_set_eamdensnum(cvalue,itype,jtype)
      else
         write(6,*) 'EAM_DENS_STYLE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_eam_embed(flag,ientry,itype,itable,dvalue)
c     ******************************************************************
c     * coefficients for the EAM embedding energies                    *
c     *                                                                *
c     * originally written 10-02-2008 by M.G. Martin                   *
c     * last modified 10-02-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,itable,ientry
      double precision dvalue
c     --- local variables
      integer etemp,itemp,ttemp
c     --- storage variable
      double precision eam_embed
      dimension eam_embed(3,MAXTABTYPE,MAXTABLE)
      save eam_embed

      if ( flag .eq. GLB_GET ) then
         dvalue = eam_embed(ientry,itype,itable)
      elseif ( flag .eq. GLB_SET ) then
         eam_embed(ientry,itype,itable) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do etemp = 1,ientry
            do itemp = 1,itype
               do ttemp = 1,itable
                  eam_embed(etemp,itemp,ttemp) = dvalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'EAM_EMBED unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_eam_embed_num(flag,itype,ivalue)
c     ******************************************************************
c     * integer code for the EAM embedding style                       *
c     *                                                                *
c     * originally written 07-21-2013 by M.G. Martin                   *
c     * last modified 07-21-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,ivalue
c     --- storage variable
      integer eam_embed_num
      dimension eam_embed_num(MAXTABTYPE)
      save eam_embed_num

      if ( flag .eq. GLB_GET ) then
         ivalue = eam_embed_num(itype)
      elseif ( flag .eq. GLB_SET ) then
         eam_embed_num(itype) = ivalue
      else
         write(6,*) 'EAM_EMBED_NUM unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_eam_embed_style(flag,itype,cvalue)
c     ******************************************************************
c     * style of the EAM embedding style                               *
c     *                                                                *
c     * originally written 10-03-2008 by M.G. Martin                   *
c     * last modified 07-22-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype
      character*(*) cvalue
c     --- storage variable
      character*20 eam_embed_style
      dimension eam_embed_style(MAXTABTYPE)
      save eam_embed_style

      if ( flag .eq. GLB_GET ) then
         cvalue = eam_embed_style(itype)
      elseif ( flag .eq. GLB_SET ) then
         eam_embed_style(itype) = cvalue
         call twh_set_eamembednum(cvalue,itype)
      else
         write(6,*) 'EAM_EMBED_STYLE unknown flag:',flag
         stop
      endif

      return
      end


      subroutine twh_eam_ndens(flag,itype,jtype,ivalue)
c     ******************************************************************
c     * number of tabular density terms in the eam style of potentials *
c     *                                                                *
c     * originally written 09-25-2008 by M.G. Martin                   *
c     * last modified 09-25-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,jtype,ivalue
c     --- local variables
      integer itemp,jtemp
c     --- storage variable
      integer eam_ndens
      dimension eam_ndens(MAXTABTYPE,MAXTABTYPE)
      save eam_ndens

      if ( flag .eq. GLB_GET ) then
         ivalue = eam_ndens(itype,jtype)
      elseif ( flag .eq. GLB_SET ) then
         eam_ndens(itype,jtype) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do itemp = 1,itype
            do jtemp = 1,jtype
               eam_ndens(itemp,jtemp) = ivalue
            enddo
         enddo
      else
         write(6,*) 'EAM_NDENS unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_eam_nembed(flag,itype,ivalue)
c     ******************************************************************
c     * number of embedding function entries                           *
c     *                                                                *
c     * originally written 09-25-2008 by M.G. Martin                   *
c     * last modified 09-25-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer itemp
c     --- storage variable
      integer eam_nembed
      dimension eam_nembed(MAXTABTYPE)
      save eam_nembed

      if ( flag .eq. GLB_GET ) then
         ivalue = eam_nembed(itype)
      elseif ( flag .eq. GLB_SET ) then
         eam_nembed(itype) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do itemp = 1,itype
            eam_nembed(itemp) = ivalue
         enddo
      else
         write(6,*) 'EAM_NEMBED unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_eam_pair_number(flag,ipair,ivalue)
c     ******************************************************************
c     * the eam pair potential expressed as an integer                 *
c     *                                                                *
c     * originally written 07-20-2013 by M.G. Martin                   *
c     * last modified 07-20-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ipair,ivalue
c     --- storage variable
      integer eam_pair_number
      dimension eam_pair_number((MAXTABTYPE*MAXTABTYPE+MAXTABTYPE)/2)
      save eam_pair_number

      if ( flag .eq. GLB_GET ) then
         ivalue = eam_pair_number(ipair)
      elseif ( flag .eq. GLB_SET ) then
         eam_pair_number(ipair) = ivalue
      else
         write(6,*) 'EAM_PAIR_NUMBER unknown flag:',flag
         stop
      endif

      return
      end
      subroutine twh_eam_pair_style(flag,ipair,cvalue)
c     ******************************************************************
c     * style for each pair of atoms in the EAM potential              *
c     *                                                                *
c     * originally written 10-03-2008 by M.G. Martin                   *
c     * last modified 07-21-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ipair
      character*(*) cvalue
c     --- storage variable
      character*20 eam_pair_style
      dimension eam_pair_style((MAXTABTYPE*MAXTABTYPE+MAXTABTYPE)/2)
      save eam_pair_style
c     --- local variables
      integer icount

      if ( flag .eq. GLB_GET ) then
         cvalue = eam_pair_style(ipair)
      elseif ( flag .eq. GLB_SET ) then
         eam_pair_style(ipair) = cvalue
c        --- also set eam_pair_num
         call twh_set_eampairnum(cvalue,ipair)
      elseif ( flag .eq. GLB_INIT ) then
         do icount = 1,ipair
            eam_pair_style(icount) = cvalue
c           --- also set eam_pair_num
            call twh_set_eampairnum(cvalue,icount)
         enddo
      else
         write(6,*) 'EAM_PAIR_STYLE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_eam_rho(flag,icrd,ichain,dvalue)
c     ******************************************************************
c     * current embedded atom density for each molecule in the         *
c     * simulation                                                     *
c     *                                                                *
c     * originally written 10-02-2008 by M.G. Martin                   *
c     * last modified 12-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,icrd,ichain
      double precision dvalue
c     --- local variables
      
      if ( icrd .eq. CRD_REAL ) then
         call twh_eam_rho_real(flag,ichain,dvalue)
      elseif ( icrd .eq. CRD_TEMP ) then
         call twh_eam_rho_temp(flag,ichain,dvalue)
      elseif ( icrd .eq. CRD_R_TO_T ) then
         call twh_eam_rho_real(GLB_GET,ichain,dvalue)
         call twh_eam_rho_temp(GLB_SET,ichain,dvalue)
      elseif ( icrd .eq. CRD_T_TO_R ) then
         call twh_eam_rho_temp(GLB_GET,ichain,dvalue)
         call twh_eam_rho_real(GLB_SET,ichain,dvalue)
      else
         write(6,*) 'EAM_RHO unknown flag:',flag
         write(6,*) 'coord:',icrd
         stop
      endif

      return
      end


      subroutine twh_electrostatic_form(flag,cvalue)
c     ******************************************************************
c     * the character string electrostatic form                        *
c     *                                                                *
c     * originally written 12-22-2006 by M.G. Martin                   *
c     * last modified 12-22-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*50 electrostatic_form
      save electrostatic_form

      if ( flag .eq. GLB_GET ) then
         cvalue = electrostatic_form
      elseif ( flag .eq. GLB_SET ) then
         electrostatic_form = cvalue
      else
         write(6,*) 'ELECTROSTATIC_FORM unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_elementname(flag,itype,cvalue)
c     ******************************************************************
c     * the two letter element name for each atom type in the          *
c     * simulation                                                     *
c     *                                                                *
c     * originally written 04-12-2006 by M.G. Martin                   *
c     * last modified 04-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype
      character*(*) cvalue
c     --- local variables
      integer ittemp
c     --- storage variable
      character*2 elementname
      dimension elementname(NNTYPE)
      save elementname

      if ( flag .eq. GLB_GET ) then
         cvalue = elementname(itype)
      elseif ( flag .eq. GLB_SET ) then
         elementname(itype) = cvalue
      elseif ( flag .eq. GLB_INIT ) then 
         do ittemp = 1,itype
            elementname(ittemp) = cvalue
         enddo
      else
         write(6,*) 'ELEMENTNAME unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ensemble(flag,ivalue)
c     ******************************************************************
c     * integer code for the ensemble                                  *
c     *                                                                *
c     * originally written 10-06-2008 by M.G. Martin                   *
c     * last modified 10-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer ensemble
      save ensemble

      if ( flag .eq. GLB_GET ) then
         ivalue = ensemble
      elseif ( flag .eq. GLB_SET ) then
         ensemble = ivalue
      else
         write(6,*) 'ENSEMBLE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_express(flag,dvalue)
c     ******************************************************************
c     * the external pressure in the simulation                        *
c     *                                                                *
c     * originally written 06-12-2006 by M.G. Martin                   *
c     * last modified 06-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      double precision dvalue
c     --- storage variable
      double precision express
      save express

      if ( flag .eq. GLB_GET ) then
         dvalue = express
      elseif ( flag .eq. GLB_SET ) then
         express = dvalue
      else
         write(6,*) 'EXPRESS unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_external_code(flag,cvalue)
c     ******************************************************************
c     * the name of the external code used to compute energies         *
c     *                                                                *
c     * originally written 01-12-2007 by M.G. Martin                   *
c     * last modified 01-12-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*20 external_code
      save external_code

      if ( flag .eq. GLB_GET ) then
         cvalue = external_code
      elseif ( flag .eq. GLB_SET ) then
         external_code = cvalue
      else
         write(6,*) 'EXTERNAL_CODE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ewald_calp(flag,ibox,dvalue)
c     ******************************************************************
c     * calp = kalp / boxlen                                           *
c     *                                                                *
c     * originally written 11-05-2007 by M.G. Martin                   *
c     * last modified 11-05-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox
      double precision dvalue
c     --- local variable
      integer ibtemp
c     --- storage variable
      double precision ewald_calp
      dimension ewald_calp(0:MAXBOX)
      save ewald_calp

      if ( flag .eq. GLB_GET ) then
         dvalue = ewald_calp(ibox)
      elseif ( flag .eq. GLB_SET ) then
         ewald_calp(ibox) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 0,ibox
            ewald_calp(ibtemp) = dvalue
         enddo
      else
         write(6,*) 'EWALD_CALP unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ewald_kalp(flag,dvalue)
c     ******************************************************************
c     * kalp controls the real space Ewald sum cutoff                  *
c     *                                                                *
c     * originally written 11-05-2007 by M.G. Martin                   *
c     * last modified 11-05-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      double precision dvalue
c     --- storage variable
      double precision ewald_kalp
      save ewald_kalp

      if ( flag .eq. GLB_GET ) then
         dvalue = ewald_kalp
      elseif ( flag .eq. GLB_SET ) then
         ewald_kalp = dvalue
      else
         write(6,*) 'EWALD_KALP unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ewald_kvector(flag,cordflag,ivector,ibox
     &     ,xvalue,yvalue,zvalue)
c     ******************************************************************
c     * the k-vectors in each box for the Ewald sum                    *
c     *                                                                *
c     * originally written 11-05-2007 by M.G. Martin                   *
c     * last modified 11-05-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,cordflag,ivector,ibox
      double precision xvalue,yvalue,zvalue
c     --- local variables
      integer ivtemp
c     --- storage variables
      double precision kx
      dimension kx(VECTORMAX,MAXBOX)
      save kx
      double precision ky
      dimension ky(VECTORMAX,MAXBOX)
      save ky
      double precision kz
      dimension kz(VECTORMAX,MAXBOX) 
      save kz
      double precision kxtemp
      dimension kxtemp(VECTORMAX,MAXBOX)
      save kxtemp
      double precision kytemp
      dimension kytemp(VECTORMAX,MAXBOX)
      save kytemp
      double precision kztemp
      dimension kztemp(VECTORMAX,MAXBOX) 
      save kztemp

      if ( flag .eq. GLB_GET ) then
         if ( cordflag .eq. CRD_REAL ) then
            xvalue = kx(ivector,ibox)
            yvalue = ky(ivector,ibox)
            zvalue = kz(ivector,ibox)
            return
         elseif ( cordflag .eq. CRD_TEMP ) then
            xvalue = kxtemp(ivector,ibox)
            yvalue = kytemp(ivector,ibox)
            zvalue = kztemp(ivector,ibox)
            return
         endif
      elseif ( flag .eq. GLB_SET ) then
         if ( cordflag .eq. CRD_REAL ) then
            kx(ivector,ibox) = xvalue
            ky(ivector,ibox) = yvalue
            kz(ivector,ibox) = zvalue
            return
         elseif ( cordflag .eq. CRD_TEMP ) then
            kxtemp(ivector,ibox) = xvalue
            kytemp(ivector,ibox) = yvalue
            kztemp(ivector,ibox) = zvalue
            return
         endif
      elseif ( flag .eq. GLB_STORE ) then
         if ( cordflag .eq. CRD_R_TO_T ) then
            do ivtemp = 1,ivector
               kxtemp(ivtemp,ibox) = kx(ivtemp,ibox)
               kytemp(ivtemp,ibox) = ky(ivtemp,ibox)
               kztemp(ivtemp,ibox) = kz(ivtemp,ibox)
            enddo
            return
         elseif ( cordflag .eq. CRD_T_TO_R ) then
            do ivtemp = 1,ivector
               kx(ivtemp,ibox) = kxtemp(ivtemp,ibox)
               ky(ivtemp,ibox) = kytemp(ivtemp,ibox)
               kz(ivtemp,ibox) = kztemp(ivtemp,ibox)
            enddo
            return
         endif
      endif

c     --- error if we got this far without returning
      write(6,*) 'EWALD_KVECTOR unknown flag and cordflag combo:'
     &     ,flag,cordflag
      stop

      return
      end

      subroutine twh_ewald_numvect(flag,ibox,ivalue)
c     ******************************************************************
c     * the number of vectors in each box for the Ewald sum            *
c     * numvect, the total number of reciprocal vectors                *
c     *                                                                *
c     * originally written 11-05-2007 by M.G. Martin                   *
c     * last modified 11-05-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox,ivalue
c     --- local variables
      integer ibtemp
c     --- storage variable
      integer ewaldnumvect
      dimension ewaldnumvect(MAXBOX)
      save ewaldnumvect

      if ( flag .eq. GLB_GET ) then
         ivalue = ewaldnumvect(ibox)
      elseif ( flag .eq. GLB_SET ) then
         ewaldnumvect(ibox) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 1,ibox
            ewaldnumvect(ibtemp) = ivalue
         enddo
      else
         write(6,*) 'EWALD_NUMVECT unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ewald_prec(flag,dvalue)
c     ******************************************************************
c     * the desired precision of the Ewald sum                         *
c     *                                                                *
c     * originally written 11-05-2007 by M.G. Martin                   *
c     * last modified 11-05-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      double precision dvalue
c     --- storage variable
      double precision ewald_prec
      save ewald_prec

      if ( flag .eq. GLB_GET ) then
         dvalue = ewald_prec
      elseif ( flag .eq. GLB_SET ) then
         ewald_prec = dvalue
      else
         write(6,*) 'EWALD_PREC unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ewald_prefact(flag,cordflag,ivector,ibox,dvalue)
c     ******************************************************************
c     * the prefactor for each vectors in each box for the Ewald sum   *
c     *                                                                *
c     * originally written 11-05-2007 by M.G. Martin                   *
c     * last modified 11-05-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,cordflag,ivector,ibox
      double precision dvalue
c     --- local variables
      integer ivtemp
c     --- storage variable
      double precision ewaldprefact
      dimension ewaldprefact(VECTORMAX,MAXBOX)
      save ewaldprefact
      double precision ewaldprefacttemp
      dimension ewaldprefacttemp(VECTORMAX,MAXBOX)
      save ewaldprefacttemp

      if ( flag .eq. GLB_GET ) then
         if ( cordflag .eq. CRD_REAL ) then
            dvalue = ewaldprefact(ivector,ibox)
            return
         elseif ( cordflag .eq. CRD_TEMP ) then
            dvalue = ewaldprefacttemp(ivector,ibox)
            return
         endif
      elseif ( flag .eq. GLB_SET ) then
         if ( cordflag .eq. CRD_REAL ) then
            ewaldprefact(ivector,ibox) = dvalue
            return
         elseif ( cordflag .eq. CRD_TEMP ) then
            ewaldprefacttemp(ivector,ibox) = dvalue
            return
         endif
      elseif ( flag .eq. GLB_STORE ) then
         if ( cordflag .eq. CRD_R_TO_T ) then
            do ivtemp = 1,ivector
               ewaldprefacttemp(ivtemp,ibox) = ewaldprefact(ivtemp,ibox)
            enddo
            return
         elseif ( cordflag .eq. CRD_T_TO_R ) then
            do ivtemp = 1,ivector
               ewaldprefact(ivtemp,ibox) = ewaldprefacttemp(ivtemp,ibox)
            enddo
            return
         endif
      endif

c     --- error if we got this far without returning
      write(6,*) 'EWALD_PREFACT unknown flag and cordflag combo:'
     &     ,flag,cordflag
      stop

      return
      end

      subroutine twh_ewald_ssumi(flag,cordflag,ivector,ibox,dvalue)
c     ******************************************************************
c     * partial summations used to compute the Ewald sum               *
c     *                                                                *
c     * originally written 11-05-2007 by M.G. Martin                   *
c     * last modified 11-05-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,cordflag,ivector,ibox
      double precision dvalue
c     --- local variables
      integer ivtemp
c     --- storage variable
      double precision ssumi
      dimension ssumi(VECTORMAX,MAXBOX)
      double precision ssuminew
      dimension ssuminew(VECTORMAX,MAXBOX)
      double precision ssumiold
      dimension ssumiold(VECTORMAX,MAXBOX)
      save ssumiold

      if ( flag .eq. GLB_GET ) then
         if ( cordflag .eq. CRD_REAL ) then
            dvalue = ssumi(ivector,ibox)
            return
         elseif ( cordflag .eq. CRD_NEW ) then
            dvalue = ssuminew(ivector,ibox)
            return
         elseif ( cordflag .eq. CRD_OLD ) then
            dvalue = ssumiold(ivector,ibox)
            return
         endif
      elseif ( flag .eq. GLB_SET ) then
         if ( cordflag .eq. CRD_REAL ) then
            ssumi(ivector,ibox) = dvalue
            return
         elseif ( cordflag .eq. CRD_NEW ) then
            ssuminew(ivector,ibox) = dvalue
            return
         elseif ( cordflag .eq. CRD_OLD ) then
            ssumiold(ivector,ibox) = dvalue
            return
         endif
      elseif ( flag .eq. GLB_STORE ) then
         if ( cordflag .eq. CRD_R_TO_N ) then
            do ivtemp = 1,ivector
               ssuminew(ivtemp,ibox) = ssumi(ivtemp,ibox)
            enddo
            return
         elseif ( cordflag .eq. CRD_N_TO_R ) then
            do ivtemp = 1,ivector
               ssumi(ivtemp,ibox) = ssuminew(ivtemp,ibox)
            enddo
            return
         elseif ( cordflag .eq. CRD_R_TO_O ) then
            do ivtemp = 1,ivector
               ssumiold(ivtemp,ibox) = ssumi(ivtemp,ibox)
            enddo
            return
         elseif ( cordflag .eq. CRD_O_TO_R ) then
            do ivtemp = 1,ivector
               ssumi(ivtemp,ibox) = ssumiold(ivtemp,ibox)
            enddo
            return
         endif
      endif

c     --- error if we got this far without returning
      write(6,*) 'EWALD_SSUMI unknown flag and cordflag combo:'
     &     ,flag,cordflag
      stop

      return
      end

      subroutine twh_ewald_ssumr(flag,cordflag,ivector,ibox,dvalue)
c     ******************************************************************
c     * partial summations used to compute the Ewald sum               *
c     *                                                                *
c     * originally written 11-05-2007 by M.G. Martin                   *
c     * last modified 11-05-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,cordflag,ivector,ibox
      double precision dvalue
c     --- local variables
      integer ivtemp
c     --- storage variable
      double precision ssumr
      dimension ssumr(VECTORMAX,MAXBOX)
      double precision ssumrnew
      dimension ssumrnew(VECTORMAX,MAXBOX)
      double precision ssumrold
      dimension ssumrold(VECTORMAX,MAXBOX)
      save ssumrold

      if ( flag .eq. GLB_GET ) then
         if ( cordflag .eq. CRD_REAL ) then
            dvalue = ssumr(ivector,ibox)
            return
         elseif ( cordflag .eq. CRD_NEW ) then
            dvalue = ssumrnew(ivector,ibox)
            return
         elseif ( cordflag .eq. CRD_OLD ) then
            dvalue = ssumrold(ivector,ibox)
            return
         endif
      elseif ( flag .eq. GLB_SET ) then
         if ( cordflag .eq. CRD_REAL ) then
            ssumr(ivector,ibox) = dvalue
            return
         elseif ( cordflag .eq. CRD_NEW ) then
            ssumrnew(ivector,ibox) = dvalue
            return
         elseif ( cordflag .eq. CRD_OLD ) then
            ssumrold(ivector,ibox) = dvalue
            return
         endif
      elseif ( flag .eq. GLB_STORE ) then
         if ( cordflag .eq. CRD_R_TO_N ) then
            do ivtemp = 1,ivector
               ssumrnew(ivtemp,ibox) = ssumr(ivtemp,ibox)
            enddo
            return
         elseif ( cordflag .eq. CRD_N_TO_R ) then
            do ivtemp = 1,ivector
               ssumr(ivtemp,ibox) = ssumrnew(ivtemp,ibox)
            enddo
            return
         elseif ( cordflag .eq. CRD_R_TO_O ) then
            do ivtemp = 1,ivector
               ssumrold(ivtemp,ibox) = ssumr(ivtemp,ibox)
            enddo
            return
         elseif ( cordflag .eq. CRD_O_TO_R ) then
            do ivtemp = 1,ivector
               ssumr(ivtemp,ibox) = ssumrold(ivtemp,ibox)
            enddo
            return
         endif
      endif

c     --- error if we got this far without returning
      write(6,*) 'EWALD_SSUMR unknown flag and cordflag combo:'
     &     ,flag,cordflag
      stop

      return
      end

      subroutine twh_excludelist(flag,imolty,iunit,icount,ivalue)
c     ******************************************************************
c     * contains the arrays that are used for the inclusion table      *
c     * these should only need to be used by the subroutine and        *
c     * function in the file linclude.F                                *
c     *                                                                *
c     * originally written 12-03-2007 by M.G. Martin                   *
c     * last modified 12-03-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,iunit,icount,ivalue
c     --- information variables
      integer excludelist
      dimension excludelist(NTMAX,NUMAX,MAXTOR+MAXBEND+NNBOND)
      save excludelist

      if ( flag .eq. GLB_GET ) then
         ivalue = excludelist(imolty,iunit,icount)
      elseif ( flag .eq. GLB_SET ) then
         excludelist(imolty,iunit,icount) = ivalue
      else
         write(6,*) 'EXCLUDELIST: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_ff_filename(flag,ifile,cvalue)
c     ******************************************************************
c     * string names for the forcefield files                          *
c     *                                                                *
c     * originally written 12-22-2007 by M.G. Martin                   *
c     * last modified 12-22-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ifile
      character*(*) cvalue
c     --- storage variable
      character*(MAXDIRLENGTH) ff_filename
      dimension ff_filename(10)
      save ff_filename

      if ( flag .eq. GLB_GET ) then
         cvalue = ff_filename(ifile)
      elseif ( flag .eq. GLB_SET ) then
         ff_filename(ifile) = cvalue
      else
         write(6,*) 'FF_FILENAME: invalid get/set flag',flag
         stop
      endif

      return
      end

      subroutine twh_ffnumber(flag,ivalue)
c     ******************************************************************
c     * number of forcefield files to read in this simulation          *
c     *                                                                *
c     * originally written 12-22-2007 by M.G. Martin                   *
c     * last modified 12-22-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- storage variable
      integer ffnumber
      save ffnumber

      if ( flag .eq. GLB_GET ) then
         ivalue = ffnumber
      elseif ( flag .eq. GLB_SET ) then
         ffnumber = ivalue
      else
         write(6,*) 'FFNUMBER: invalid get/set flag',flag
         stop
      endif

      return
      end

      subroutine twh_fieldchar(flag,ifield,icount,cvalue)
c     ******************************************************************
c     * character string parameters used to define the fields          *
c     *                                                                *
c     * originally written 11-09-2007 by M.G. Martin                   *
c     * last modified 11-09-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ifield,icount
      character*(*) cvalue
c     --- storage variable
      integer maxfieldchar
      parameter(maxfieldchar=3)
      character*11 fieldchar
      dimension fieldchar(FLDMAX,maxfieldchar)
      save fieldchar

      if ( flag .eq. GLB_GET ) then
         cvalue = fieldchar(ifield,icount)
      elseif ( flag .eq. GLB_SET ) then
         fieldchar(ifield,icount) = cvalue
      else
         write(6,*) 'FIELDCHAR: invalid get/set flag',flag
         stop
      endif

      return
      end

      subroutine twh_fieldcoordinates(flag,ichain,iunit,xxx,yyy,zzz)
c     ******************************************************************
c     * get or sets the initial coordinates of a single atom in the    *
c     * coordinit array                                                *
c     *                                                                *
c     * originally written 08-04-2003 by M.G. Martin                   *
c     * last modified 12-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETATOMNUM
#include "functions.h"
c     --- variables passed to/from this subroutine
      integer flag
      integer ichain,iunit
      double precision xxx,yyy,zzz
c     --- local variables
      integer index
      double precision fieldvector
      dimension fieldvector(3)

      index = twh_getatomnum(ichain,iunit)

      if ( flag .eq. GLB_SET ) then
         fieldvector(1) = xxx
         fieldvector(2) = yyy
         fieldvector(3) = zzz
         call twh_coordfield(GLB_SET_TRIPLE,index,fieldvector)
      elseif ( flag .eq. GLB_GET ) then
         call twh_coordfield(GLB_GET_TRIPLE,index,fieldvector)
         xxx = fieldvector(1)
         yyy = fieldvector(2)
         zzz = fieldvector(3)
      else
         write(6,*) 'FIELDCOORDS: invalid get/set flag',flag
         stop
      endif

      return
      end

      subroutine twh_fielddp(flag,ifield,icount,dvalue)
c     ******************************************************************
c     * double precision parameters used to define the fields          *
c     *                                                                *
c     * originally written 11-08-2007 by M.G. Martin                   *
c     * last modified 10-27-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ifield,icount
      double precision dvalue
c     --- local variables
      integer iftemp,ictemp
c     --- storage variable
      double precision fielddp
      dimension fielddp(FLDMAX,MAXFIELDDP)
      save fielddp

      if ( flag .eq. GLB_GET ) then
         dvalue = fielddp(ifield,icount)
      elseif ( flag .eq. GLB_SET ) then
         fielddp(ifield,icount) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do iftemp = 1,ifield
            do ictemp = 1,icount
               fielddp(iftemp,ictemp) = dvalue
            enddo
         enddo
      else
         write(6,*) 'FIELDDP: invalid get/set flag',flag
         stop
      endif

      return
      end

      subroutine twh_fieldint(flag,ifield,icount,ivalue)
c     ******************************************************************
c     * integer parameters used to define the fields                   *
c     *                                                                *
c     * originally written 11-08-2007 by M.G. Martin                   *
c     * last modified 10-27-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ifield,icount,ivalue
c     --- local variables
      integer iftemp,ictemp
c     --- storage variable
      integer fieldint
      dimension fieldint(FLDMAX,MAXFIELDINT)
      save fieldint

      if ( flag .eq. GLB_GET ) then
         ivalue = fieldint(ifield,icount)
      elseif ( flag .eq. GLB_SET ) then
         fieldint(ifield,icount) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do iftemp = 1,ifield
            do ictemp = 1,icount
               fieldint(iftemp,ictemp) = ivalue
            enddo
         enddo
      else
         write(6,*) 'FIELDINT: invalid get/set flag',flag
         stop
      endif

      return
      end

      subroutine twh_fieldlogical(flag,ifield,lvalue)
c     ******************************************************************
c     * logical parameters used to define the fields                   *
c     *                                                                *
c     * originally written 12-01-2007 by M.G. Martin                   *
c     * last modified 12-02-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ifield
      logical lvalue
c     --- storage variable
      logical fieldlogical
      dimension fieldlogical(FLDMAX)
      save fieldlogical

      if ( flag .eq. GLB_GET ) then
         lvalue = fieldlogical(ifield)
      elseif ( flag .eq. GLB_SET ) then
         fieldlogical(ifield) = lvalue
      else
         write(6,*) 'FIELDLOGICAL: invalid get/set flag',flag
         stop
      endif

      return
      end

      subroutine twh_fieldnbint(flag,ifield,itype,icount,ivalue)
c     ******************************************************************
c     * integer parameters used to define the fields                   *
c     *                                                                *
c     * originally written 11-08-2007 by M.G. Martin                   *
c     * last modified 11-08-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ifield,icount,itype,ivalue
c     --- local variables
      integer iftemp,ittemp,ictemp
c     --- storage variable
      integer maxfieldnbint
      parameter(maxfieldnbint=2)
      integer fieldnbint
      dimension fieldnbint(FLDMAX,NNTYPE,maxfieldnbint)
      save fieldnbint

      if ( flag .eq. GLB_GET ) then
         ivalue = fieldnbint(ifield,itype,icount)
      elseif ( flag .eq. GLB_SET ) then
         fieldnbint(ifield,itype,icount) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do iftemp = 1,ifield
            do ittemp = 1,itype
               do ictemp = 1,maxfieldnbint
                  fieldnbint(iftemp,ittemp,ictemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'FIELDNBINT: invalid get/set flag',flag
         stop
      endif

      return
      end

      subroutine twh_fieldnbparam(flag,ifield,itype,icount,dvalue)
c     ******************************************************************
c     * nonbonded parameters used in the fields                        *
c     *                                                                *
c     * originally written 11-08-2007 by M.G. Martin                   *
c     * last modified 12-14-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ifield,icount,itype
      double precision dvalue
c     --- local variables
      integer iftemp,ittemp,ictemp
c     --- storage variable
      double precision fieldnbparam
      dimension fieldnbparam(FLDMAX,NNTYPE,MAXFIELDNBPARAM)
      save fieldnbparam

      if ( flag .eq. GLB_GET ) then
         dvalue = fieldnbparam(ifield,itype,icount)
      elseif ( flag .eq. GLB_SET ) then
         fieldnbparam(ifield,itype,icount) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do iftemp = 1,ifield
            do ittemp = 1,itype
               do ictemp = 1,icount
                  fieldnbparam(iftemp,ittemp,ictemp) = dvalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'FIELDNBPARAM: invalid get/set flag',flag
         stop
      endif

      return
      end

      subroutine twh_fieldtype(flag,ifield,ivalue)
c     ******************************************************************
c     * integer code type for each field                               *
c     *                                                                *
c     * originally written 11-08-2007 by M.G. Martin                   *
c     * last modified 11-08-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ifield,ivalue
c     --- storage variable
      integer fieldtype
      dimension fieldtype(FLDMAX)
      save fieldtype

      if ( flag .eq. GLB_GET ) then
         ivalue = fieldtype(ifield)
      elseif ( flag .eq. GLB_SET ) then
         fieldtype(ifield) = ivalue
      else
         write(6,*) 'FIELDTYPE: invalid get/set flag',flag
         stop
      endif

      return
      end


      subroutine twh_fieldtypestring(flag,ifield,cvalue)
c     ******************************************************************
c     * character string code type for each field                      *
c     *                                                                *
c     * originally written 11-08-2007 by M.G. Martin                   *
c     * last modified 08-09-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ifield
      character*(*) cvalue
c     --- local variable
      integer itype
c     --- storage variable
      character*20 fieldtypestring
      dimension fieldtypestring(FLDMAX)
      save fieldtypestring


      if ( flag .eq. GLB_GET ) then
         cvalue = fieldtypestring(ifield)
      elseif ( flag .eq. GLB_SET ) then
         fieldtypestring(ifield) = cvalue
c        --- also need to set the integer code
         if ( cvalue .eq. 'Harmonic Attractor' ) then
            itype = FLD_HARMONIC_ATTRACTOR
         elseif ( cvalue .eq. 'Hard Wall' ) then
            itype = FLD_HARD_WALL
         elseif ( cvalue .eq. 'Hard Ball' ) then
            itype = FLD_HARD_BALL
         elseif ( cvalue .eq. 'LJ 9-3 Wall' ) then
            itype = FLD_LJ_9_3_WALL
         elseif ( cvalue .eq. 'Hooper Umbrella' ) then
            itype = FLD_HOOPER_UMBRELLA
         elseif ( cvalue .eq. 'Steele Wall' ) then
            itype = FLD_STEELE_WALL
         else
            write(6,*) 'FIELDTYPESTRING: unknown fieldtype'
            write(6,*) cvalue
            stop
         endif
         call twh_fieldtype(GLB_SET,ifield,itype)
      else
         write(6,*) 'FIELDTYPE: invalid get/set flag',flag
         stop
      endif

      return
      end

      subroutine twh_file_flct(flag,cvalue)
c     ******************************************************************
c     * character string for creating grand canonical output files     *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*30 file_flct
      save file_flct

      if ( flag .eq. GLB_GET ) then
         cvalue = file_flct
      elseif ( flag .eq. GLB_SET ) then
         file_flct = cvalue
      else
         write(6,*) 'FILE_FLCT: invalid get/set flag',flag
         stop
      endif

      return
      end

      subroutine twh_file_his(flag,cvalue)
c     ******************************************************************
c     * character string for creating histogram rewieghting output     *
c     * files                                                          *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*30 file_his
      save file_his

      if ( flag .eq. GLB_GET ) then
         cvalue = file_his
      elseif ( flag .eq. GLB_SET ) then
         file_his = cvalue
      else
         write(6,*) 'FILE_HIS: invalid get/set flag',flag
         stop
      endif

      return
      end

      subroutine twh_file_ndis(flag,imolty,cvalue)
c     ******************************************************************
c     * character string for creating grand canonical output files     *
c     *                                                                *
c     * originally written 10-09-2008 by M.G. Martin                   *
c     * last modified 10-09-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      character*(*) cvalue
c     --- storage variable
      character*30 file_ndis
      dimension file_ndis(NTMAX)
      save file_ndis

      if ( flag .eq. GLB_GET ) then
         cvalue = file_ndis(imolty)
      elseif ( flag .eq. GLB_SET ) then
         file_ndis(imolty) = cvalue
      else
         write(6,*) 'FILE_NDIS: invalid get/set flag',flag
         stop
      endif

      return
      end

      subroutine twh_foreign_energy(flag,ibox,iforeign,dvalue)
c     ******************************************************************
c     * foreign energies for scaled potential                          *
c     *                                                                *
c     * originally written 12-20-2007 by M.G. Martin                   *
c     * last modified 12-20-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ibox,iforeign
      double precision dvalue
c     --- local variables
      integer ibtemp,iftemp
c     --- storage variable
      double precision foreign_energy
      dimension foreign_energy(MAXBOX,MAX_FOREIGN_LAMBDA)
      save foreign_energy

      if ( flag .eq. GLB_GET ) then
         dvalue = foreign_energy(ibox,iforeign)
      elseif ( flag .eq. GLB_SET ) then
         foreign_energy(ibox,iforeign) = dvalue
      elseif ( flag .eq. GLB_INCR ) then
         foreign_energy(ibox,iforeign) = foreign_energy(ibox,iforeign)
     &        + dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 1,ibox
            do iftemp = 1,iforeign
               foreign_energy(ibtemp,iftemp)
     &              = foreign_energy(ibtemp,iftemp) + dvalue
            enddo
         enddo
      else
         write(6,*) 'FOREIGN_ENERGY: invalid flag',flag
         stop
      endif

      return
      end

      subroutine twh_foreign_lambda(flag,flambda,count,dvalue)
c     ******************************************************************
c     * foreign energy lambda values                                   *
c     *                                                                *
c     * originally written 12-22-2007 by M.G. Martin                   *
c     * last modified 12-22-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,flambda,count
      double precision dvalue
#ifndef USE_DYNAMIC_MEMORY
c     --- local variables
      integer ilambda
c     --- storage variable
      double precision foreign_lambda
      dimension foreign_lambda(MAX_FOREIGN_LAMBDA,2)
      save foreign_lambda

      if ( flag .eq. GLB_GET ) then
         dvalue = foreign_lambda(flambda,count)
      elseif ( flag .eq. GLB_SET ) then
         foreign_lambda(flambda,count) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ilambda = 1,flambda
            foreign_lambda(ilambda,count) = dvalue
         enddo
      else
         write(6,*) 'FOREIGN_LAMBDA: invalid flag',flag
         stop
      endif
#endif

      return
      end

      function twh_getatomnum(ichain,iunit)
c     ******************************************************************
c     * returns the atom position in the coordinate arrays             *
c     *                                                                *
c     * originally written 03-25-2002 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETATOMNUM
#include "functions.h"
c     --- variables passed to the function
      integer ichain,iunit,globalpos

      call twh_globalpos(GLB_GET,ichain,globalpos)
      twh_getatomnum = globalpos + 3*(iunit-1)
      return
      end

      subroutine twh_growbond(flag,iunit,ibond,jbond,ivalue)
c     ******************************************************************
c     * growbond data structure used to control the logic of CBMC moves*
c     *                                                                *
c     * originally written 10-04-2008 by M.G. Martin                   *
c     * last modified 10-04-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,iunit,ibond,jbond,ivalue
c     --- storage variable
      integer growbond
      dimension growbond(NUMAX,NNBOND,NNBOND)
      save growbond

      if ( flag .eq. GLB_GET ) then
         ivalue = growbond(iunit,ibond,jbond)
      elseif ( flag .eq. GLB_SET ) then
         growbond(iunit,ibond,jbond) = ivalue
      else
         write(6,*) 'GROWBOND unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_growcount(flag,imolty,iunit,ibond,ivalue)
c     ******************************************************************
c     * the number of atoms that would be regrown when starting from   *
c     * atom iunit on molecule type imolty when ibond is the growprev  *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibond,ivalue
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer growcount
      dimension growcount(NTMAX,NUMAX,NNBOND)
      save growcount

      if ( flag .eq. GLB_GET ) then
         ivalue = growcount(imolty,iunit,ibond)
      elseif ( flag .eq. GLB_SET ) then
         growcount(imolty,iunit,ibond) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibond
                  growcount(imtemp,iutemp,ibtemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'GROWCOUNT unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_growlist(flag,iunit,ibond,ivalue)
c     ******************************************************************
c     * the growlist value for each step of the CBMC growth process    *
c     *                                                                *
c     * originally written 10-03-2008 by M.G. Martin                   *
c     * last modified 10-03-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,iunit,ibond,ivalue
c     --- storage variable
      integer growlist
      dimension growlist(NUMAX,NNBOND)
      save growlist

      if ( flag .eq. GLB_GET ) then
         ivalue = growlist(iunit,ibond)
      elseif ( flag .eq. GLB_SET ) then
         growlist(iunit,ibond) = ivalue
      else
         write(6,*) 'GROWLIST unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_grownbon(flag,iunit,ibond,ivalue)
c     ******************************************************************
c     * the grownbon value for each step of the CBMC growth process    *
c     *                                                                *
c     * originally written 10-04-2008 by M.G. Martin                   *
c     * last modified 10-04-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,iunit,ibond,ivalue
c     --- local variables
      integer iutemp,ibtemp
c     --- storage variable
      integer grownbon
      dimension grownbon(NUMAX,NNBOND)
      save grownbon

      if ( flag .eq. GLB_GET ) then
         ivalue = grownbon(iunit,ibond)
      elseif ( flag .eq. GLB_SET ) then
         grownbon(iunit,ibond) = ivalue
      elseif ( flag .eq. GLB_INCR ) then
         grownbon(iunit,ibond) = grownbon(iunit,ibond) + ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do iutemp = 1,iunit
            do ibtemp = 1,ibond
               grownbon(iutemp,ibtemp) = ivalue
            enddo
         enddo
      else
         write(6,*) 'GROWNBON unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_growvalidlist(flag,imolty,icount,ivalue)
c     ******************************************************************
c     * a list of the valid atoms to originate a CBMC growth from for  *
c     * each molecule type                                             *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,icount,ivalue
c     --- local variables
      integer imtemp,ictemp
c     --- storage variable
      integer growvalidlist
      dimension growvalidlist(NTMAX,NUMAX)
      save growvalidlist

      if ( flag .eq. GLB_GET ) then
         ivalue = growvalidlist(imolty,icount)
      elseif ( flag .eq. GLB_SET ) then
         growvalidlist(imolty,icount) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do ictemp = 1,icount
               growvalidlist(imtemp,ictemp) = ivalue
            enddo
         enddo
      else
         write(6,*) 'GROWVALIDLIST unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_growxtra(flag,iunit,ibond,ivalue)
c     ******************************************************************
c     * a list of the extra atoms (due to ring closing) for CBMC       *
c     * regrowth moves                                                 *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modified 10-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,iunit,ibond,ivalue
c     --- storage variable
      integer growxtra
      dimension growxtra(NUMAX,NNBOND+1)
      save growxtra

      if ( flag .eq. GLB_GET ) then
         ivalue = growxtra(iunit,ibond)
      elseif ( flag .eq. GLB_SET ) then
         growxtra(iunit,ibond) = ivalue
      elseif ( flag .eq. GLB_INCR ) then
         growxtra(iunit,ibond) = growxtra(iunit,ibond) + ivalue
      else
         write(6,*) 'GROWXTRA unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_halfbox(flag,ibox,idim,dvalue)
c     ******************************************************************
c     * half of the box length in each dimension.  Useful for mimage   *
c     *                                                                *
c     * originally written 09-25-2008 by M.G. Martin                   *
c     * last modified 09-25-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox,idim
      double precision dvalue
c     --- storage variable
      double precision halfbox
      dimension halfbox(MAXBOX,3)
      save halfbox

      if ( flag .eq. GLB_GET ) then
         dvalue = halfbox(ibox,idim)
      elseif ( flag .eq. GLB_SET ) then
         halfbox(ibox,idim) = dvalue
      else
         write(6,*) 'HALFBOX unknown flag:',flag
         stop
      endif
      
      return
      end

      subroutine twh_hbondcoeff(flag,itype,icoeff,dvalue)
c     ******************************************************************
c     * the coefficients for special hydrogen bonds                    *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,icoeff
      double precision dvalue
c     --- local variables
      integer ittemp,ictemp
c     --- storage variable
      double precision hbondcoeff
      dimension hbondcoeff(THBONDMAX,MINHBONDCOEFF:MAXHBONDCOEFF)
      save hbondcoeff

      if ( flag .eq. GLB_GET ) then
         dvalue = hbondcoeff(itype,icoeff)
      elseif ( flag .eq. GLB_SET ) then
         hbondcoeff(itype,icoeff) = dvalue
      elseif ( flag .eq. GLB_INIT ) then 
         do ittemp = 1,itype
            do ictemp = MINHBONDCOEFF,icoeff
               hbondcoeff(ittemp,ictemp) = dvalue
            enddo
         enddo
      else
         write(6,*) 'HBONDCOEFF unknown flag:',flag
         stop
      endif
      
      return
      end

      subroutine twh_hbondnames(flag,itype,index,cvalue)
c     ******************************************************************
c     * the two atom names involved in special hydrogen bonds          *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,index
      character*(*) cvalue
c     --- local variables
      integer ittemp,iitemp
c     --- storage variable
      character*(FFNAMELEN) hbondnames
      dimension hbondnames(THBONDMAX,MAXHBONDNAME)
      save hbondnames

      if ( flag .eq. GLB_GET ) then
         cvalue = hbondnames(itype,index)
      elseif ( flag .eq. GLB_SET ) then
         hbondnames(itype,index) = cvalue
      elseif ( flag .eq. GLB_INIT ) then 
         do ittemp = 1,itype
            do iitemp = 1,index
               hbondnames(ittemp,iitemp) = cvalue
            enddo
         enddo
      else
         write(6,*) 'HBONDNAMES unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_helix_angle(flag,imolty,dvalue)
c     ******************************************************************
c     * the angle of incline to use when building a molecule as a helix*
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty
      double precision dvalue
c     --- storage variable
      double precision helix_angle
      dimension helix_angle(NTMAX)
      save helix_angle

      if ( flag .eq. GLB_GET ) then
         dvalue = helix_angle(imolty)
      elseif ( flag .eq. GLB_SET ) then
         helix_angle(imolty) = dvalue
      else
         write(6,*) 'HELIX_ANGLE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_helix_conlen(flag,imolty,dvalue)
c     ******************************************************************
c     * the connection length to use when building a molecule as a     *
c     * helix                                                          *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty
      double precision dvalue
c     --- storage variable
      double precision helix_conlen
      dimension helix_conlen(NTMAX)
      save helix_conlen

      if ( flag .eq. GLB_GET ) then
         dvalue = helix_conlen(imolty)
      elseif ( flag .eq. GLB_SET ) then
         helix_conlen(imolty) = dvalue
      else
         write(6,*) 'HELIX_CONLEN unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_helix_keyname(flag,imolty,cvalue)
c     ******************************************************************
c     * the name to use for matching atoms when building a helix       *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty
      character*(*) cvalue
c     --- storage variable
      character*10 helix_keyname
      dimension helix_keyname(NTMAX)
      save helix_keyname

      if ( flag .eq. GLB_GET ) then
         cvalue = helix_keyname(imolty)
      elseif ( flag .eq. GLB_SET ) then
         helix_keyname(imolty) = cvalue
      else
         write(6,*) 'HELIX_KEYNAME unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_helix_keytype(flag,imolty,cvalue)
c     ******************************************************************
c     * the type of key to use when building a molecule as a helix     *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty
      character*(*) cvalue
c     --- storage variable
      character*10 helix_keytype
      dimension helix_keytype(NTMAX)
      save helix_keytype

      if ( flag .eq. GLB_GET ) then
         cvalue = helix_keytype(imolty)
      elseif ( flag .eq. GLB_SET ) then
         helix_keytype(imolty) = cvalue
      else
         write(6,*) 'HELIX_KEYTYPE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_helix_phase(flag,imolty,dvalue)
c     ******************************************************************
c     * the phase adjustment to use when building a molecule as a helix*
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty
      double precision dvalue
c     --- storage variable
      double precision helix_phase
      dimension helix_phase(NTMAX)
      save helix_phase

      if ( flag .eq. GLB_GET ) then
         dvalue = helix_phase(imolty)
      elseif ( flag .eq. GLB_SET ) then
         helix_phase(imolty) = dvalue
      else
         write(6,*) 'HELIX_PHASE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_helix_radius(flag,imolty,dvalue)
c     ******************************************************************
c     * the radius to use when building a molecule as a helix          *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty
      double precision dvalue
c     --- storage variable
      double precision helix_radius
      dimension helix_radius(NTMAX)
      save helix_radius

      if ( flag .eq. GLB_GET ) then
         dvalue = helix_radius(imolty)
      elseif ( flag .eq. GLB_SET ) then
         helix_radius(imolty) = dvalue
      else
         write(6,*) 'HELIX_RADIUS unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_hinverse(flag,ibox,idim,jdim,dvalue)
c     ******************************************************************
c     * the inverse of the matrix that defines the box dimensions      *
c     *                                                                *
c     * originally written 09-24-2008 by M.G. Martin                   *
c     * last modified 09-24-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox,idim,jdim
      double precision dvalue
c     --- local variables
      integer ibtemp,iitemp,ijtemp
c     --- storage variable
      double precision hinverse
      dimension hinverse(MAXBOX,3,3)
      save hinverse

      if ( flag .eq. GLB_GET ) then
         dvalue = hinverse(ibox,idim,jdim)
      elseif ( flag .eq. GLB_SET ) then
         hinverse(ibox,idim,jdim) = dvalue
      elseif ( flag .eq. GLB_SCALE ) then
         hinverse(ibox,idim,jdim) = hinverse(ibox,idim,jdim) * dvalue
      elseif ( flag .eq. GLB_INIT ) then 
         do ibtemp = 1,ibox
            do iitemp = 1,idim
               do ijtemp = 1,jdim
                  hinverse(ibtemp,iitemp,ijtemp) = dvalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'HINVERSE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_hist_label(flag,ivalue)
c     ******************************************************************
c     * integer code for the histogram label                           *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer hist_label
      save hist_label

      if ( flag .eq. GLB_GET ) then
         ivalue = hist_label
      elseif ( flag .eq. GLB_SET ) then
         hist_label = ivalue
      else
         write(6,*) 'HIST_LABEL unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_hist_nequil(flag,ivalue)
c     ******************************************************************
c     * number of equilibration cycles to discard for a histogram      *
c     * rewieghting run                                                *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer hist_nequil
      save hist_nequil

      if ( flag .eq. GLB_GET ) then
         ivalue = hist_nequil
      elseif ( flag .eq. GLB_SET ) then
         hist_nequil = ivalue
      else
         write(6,*) 'HIST_NEQUIL unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_hist_suffix(flag,cvalue)
c     ******************************************************************
c     * single character that goes at the end of the histogram output  *
c     * files                                                          *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*1 hist_suffix
      save hist_suffix

      if ( flag .eq. GLB_GET ) then
         cvalue = hist_suffix
      elseif ( flag .eq. GLB_SET ) then
         hist_suffix = cvalue
      else
         write(6,*) 'HIST_SUFFIX unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_histcalcfreq(flag,ivalue)
c     ******************************************************************
c     * frequency for calculating the histogram updates                *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer histcalcfreq
      save histcalcfreq

      if ( flag .eq. GLB_GET ) then
         ivalue = histcalcfreq
      elseif ( flag .eq. GLB_SET ) then
         histcalcfreq = ivalue
      else
         write(6,*) 'HISTCALCFREQ unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_histdumpfreq(flag,ivalue)
c     ******************************************************************
c     * frequency for outputting the histogram updates                 *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer histdumpfreq
      save histdumpfreq

      if ( flag .eq. GLB_GET ) then
         ivalue = histdumpfreq
      elseif ( flag .eq. GLB_SET ) then
         histdumpfreq = ivalue
      else
         write(6,*) 'HISTDUMPFREQ unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_hmatrix(flag,ibox,idim,jdim,dvalue)
c     ******************************************************************
c     * the matrix that defines the box dimensions                     *
c     *                                                                *
c     * originally written 09-24-2008 by M.G. Martin                   *
c     * last modified 09-24-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox,idim,jdim
      double precision dvalue
c     --- local variables
      integer ibtemp,iitemp,ijtemp
c     --- storage variable
      double precision hmatrix
      dimension hmatrix(MAXBOX,3,3)
      save hmatrix

      if ( flag .eq. GLB_GET ) then
         dvalue = hmatrix(ibox,idim,jdim)
      elseif ( flag .eq. GLB_SET ) then
         hmatrix(ibox,idim,jdim) = dvalue
      elseif ( flag .eq. GLB_SCALE ) then
         hmatrix(ibox,idim,jdim) = hmatrix(ibox,idim,jdim) * dvalue
      elseif ( flag .eq. GLB_INIT ) then 
         do ibtemp = 1,ibox
            do iitemp = 1,idim
               do ijtemp = 1,jdim
                  hmatrix(ibtemp,iitemp,ijtemp) = dvalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'HMATRIX unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_hyfrom(flag,cvalue)
c     ******************************************************************
c     * the hybridization on the from atom used for certain old        *
c     * CBMC regrowth strategies                                       *
c     *                                                                *
c     * originally written 10-03-2008 by M.G. Martin                   *
c     * last modified 10-03-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*5 hyfrom
      save hyfrom

      if ( flag .eq. GLB_GET ) then
         cvalue = hyfrom
      elseif ( flag .eq. GLB_SET ) then
         hyfrom = cvalue
      else
         write(6,*) 'HYFROM unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_hyprev(flag,cvalue)
c     ******************************************************************
c     * the hybridization on the prev atom used for certain old        *
c     * CBMC regrowth strategies                                       *
c     *                                                                *
c     * originally written 10-03-2008 by M.G. Martin                   *
c     * last modified 10-03-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*5 hyprev
      save hyprev

      if ( flag .eq. GLB_GET ) then
         cvalue = hyprev
      elseif ( flag .eq. GLB_SET ) then
         hyprev = cvalue
      else
         write(6,*) 'HYPREV unknown flag:',flag
         stop
      endif

      return
      end


      subroutine twh_iataa1(flag,imolty,iunit,iaa,ivalue)
c     ******************************************************************
c     * the single angle types for the angle-angle partners from a     *
c     * particular unit and type of molecule                           *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,iaa,ivalue
c     --- local variables
      integer imtemp,iutemp,iatemp
c     --- storage variable
      integer iataa1
      dimension iataa1(NTMAX,NUMAX,MAXAA)
      save iataa1

      if ( flag .eq. GLB_GET ) then
         ivalue = iataa1(imolty,iunit,iaa)
      elseif ( flag .eq. GLB_SET ) then
         iataa1(imolty,iunit,iaa) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do iatemp = 1,iaa
                  iataa1(imtemp,iutemp,iatemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IATAA1 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_iataa2(flag,imolty,iunit,iaa,ivalue)
c     ******************************************************************
c     * the single angle types for the angle-angle partners from a     *
c     * particular unit and type of molecule                           *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,iaa,ivalue
c     --- local variables
      integer imtemp,iutemp,iatemp
c     --- storage variable
      integer iataa2
      dimension iataa2(NTMAX,NUMAX,MAXAA)
      save iataa2

      if ( flag .eq. GLB_GET ) then
         ivalue = iataa2(imolty,iunit,iaa)
      elseif ( flag .eq. GLB_SET ) then
         iataa2(imolty,iunit,iaa) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do iatemp = 1,iaa
                  iataa2(imtemp,iutemp,iatemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IATAA2 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_iattor1(flag,imolty,iunit,itor,ivalue)
c     ******************************************************************
c     * the angle types for the regular torsion, used for some class 2 *
c     * force fields and automatically set by Towhee                   *
c     *                                                                *
c     * originally written 04-03-2006 by M.G. Martin                   *
c     * last modified 04-03-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,itor,ivalue
c     --- local variables
      integer imtemp,iutemp,ittemp
c     --- storage variable
      integer iattor1
      dimension iattor1(NTMAX,NUMAX,MAXTOR)
      save iattor1

      if ( flag .eq. GLB_GET ) then
         ivalue = iattor1(imolty,iunit,itor)
      elseif ( flag .eq. GLB_SET ) then
         iattor1(imolty,iunit,itor) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ittemp = 1,itor
                  iattor1(imtemp,iutemp,ittemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IATTOR1 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_iattor2(flag,imolty,iunit,itor,ivalue)
c     ******************************************************************
c     * the angle types for the regular torsion, used for some class 2 *
c     * force fields and automatically set by Towhee                   *
c     *                                                                *
c     * originally written 04-03-2006 by M.G. Martin                   *
c     * last modified 04-03-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,itor,ivalue
c     --- local variables
      integer imtemp,iutemp,ittemp
c     --- storage variable
      integer iattor2
      dimension iattor2(NTMAX,NUMAX,MAXTOR)
      save iattor2

      if ( flag .eq. GLB_GET ) then
         ivalue = iattor2(imolty,iunit,itor)
      elseif ( flag .eq. GLB_SET ) then
         iattor2(imolty,iunit,itor) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ittemp = 1,itor
                  iattor2(imtemp,iutemp,ittemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IATTOR2 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ibtben1(flag,imolty,iunit,ibend,ivalue)
c     ******************************************************************
c     * bond types in regular angle terms for each molecule type and   *
c     * unit                                                           *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibend,ivalue
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer ibtben1
      dimension ibtben1(NTMAX,NUMAX,MAXBEND)
      save ibtben1

      if ( flag .eq. GLB_GET ) then
         ivalue = ibtben1(imolty,iunit,ibend)
      elseif ( flag .eq. GLB_SET ) then
         ibtben1(imolty,iunit,ibend) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibend
                  ibtben1(imtemp,iutemp,ibtemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IBTBEN1 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ibtben2(flag,imolty,iunit,ibend,ivalue)
c     ******************************************************************
c     * bond types in regular angle terms for each molecule type and   *
c     * unit                                                           *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibend,ivalue
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer ibtben2
      dimension ibtben2(NTMAX,NUMAX,MAXBEND)
      save ibtben2

      if ( flag .eq. GLB_GET ) then
         ivalue = ibtben2(imolty,iunit,ibend)
      elseif ( flag .eq. GLB_SET ) then
         ibtben2(imolty,iunit,ibend) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibend
                  ibtben2(imtemp,iutemp,ibtemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IBTBEN2 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ibtbnbn1(flag,imolty,iunit,ibend,ivalue)
c     ******************************************************************
c     * types of the individual bonds used in the bond-bond terms for  *
c     * each particular unit and type of molecule                      *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibend,ivalue
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer ibtbnbn1
      dimension ibtbnbn1(NTMAX,NUMAX,MAXBEND)
      save ibtbnbn1

      if ( flag .eq. GLB_GET ) then
         ivalue = ibtbnbn1(imolty,iunit,ibend)
      elseif ( flag .eq. GLB_SET ) then
         ibtbnbn1(imolty,iunit,ibend) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibend
                  ibtbnbn1(imtemp,iutemp,ibtemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IBTBNBN1 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ibtbnbn2(flag,imolty,iunit,ibend,ivalue)
c     ******************************************************************
c     * types of the individual bonds used in the bond-bond terms for  *
c     * each particular unit and type of molecule                      *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibend,ivalue
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer ibtbnbn2
      dimension ibtbnbn2(NTMAX,NUMAX,MAXBEND)
      save ibtbnbn2

      if ( flag .eq. GLB_GET ) then
         ivalue = ibtbnbn2(imolty,iunit,ibend)
      elseif ( flag .eq. GLB_SET ) then
         ibtbnbn2(imolty,iunit,ibend) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibend
                  ibtbnbn2(imtemp,iutemp,ibtemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IBTBNBN2 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ibttor1(flag,imolty,iunit,itor,ivalue)
c     ******************************************************************
c     * bond types in a regular torsion.  These are used for some of   *
c     * the class 2 torsional force fields and are set automatically   *
c     * by Towhee                                                      *
c     *                                                                *
c     * originally written 04-03-2006 by M.G. Martin                   *
c     * last modified 04-03-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,itor,ivalue
c     --- local variables
      integer imtemp,iutemp,ittemp
c     --- storage variable
      integer ibttor1
      dimension ibttor1(NTMAX,NUMAX,MAXTOR)
      save ibttor1

      if ( flag .eq. GLB_GET ) then
         ivalue = ibttor1(imolty,iunit,itor)
      elseif ( flag .eq. GLB_SET ) then
         ibttor1(imolty,iunit,itor) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ittemp = 1,itor
                  ibttor1(imtemp,iutemp,ittemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IBTTOR1 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ibttor2(flag,imolty,iunit,itor,ivalue)
c     ******************************************************************
c     * bond types in a regular torsion.  These are used for some of   *
c     * the class 2 torsional force fields and are set automatically   *
c     * by Towhee                                                      *
c     *                                                                *
c     * originally written 04-03-2006 by M.G. Martin                   *
c     * last modified 04-03-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,itor,ivalue
c     --- local variables
      integer imtemp,iutemp,ittemp
c     --- storage variable
      integer ibttor2
      dimension ibttor2(NTMAX,NUMAX,MAXTOR)
      save ibttor2

      if ( flag .eq. GLB_GET ) then
         ivalue = ibttor2(imolty,iunit,itor)
      elseif ( flag .eq. GLB_SET ) then
         ibttor2(imolty,iunit,itor) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ittemp = 1,itor
                  ibttor2(imtemp,iutemp,ittemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IBTTOR2 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ibttor3(flag,imolty,iunit,itor,ivalue)
c     ******************************************************************
c     * bond types in a regular torsion.  These are used for some of   *
c     * the class 2 torsional force fields and are set automatically   *
c     * by Towhee                                                      *
c     *                                                                *
c     * originally written 04-03-2006 by M.G. Martin                   *
c     * last modified 04-03-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,itor,ivalue
c     --- local variables
      integer imtemp,iutemp,ittemp
c     --- storage variable
      integer ibttor3
      dimension ibttor3(NTMAX,NUMAX,MAXTOR)
      save ibttor3

      if ( flag .eq. GLB_GET ) then
         ivalue = ibttor3(imolty,iunit,itor)
      elseif ( flag .eq. GLB_SET ) then
         ibttor3(imolty,iunit,itor) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ittemp = 1,itor
                  ibttor3(imtemp,iutemp,ittemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IBTTOR3 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ijaa0(flag,imolty,iunit,iaa,ivalue)
c     ******************************************************************
c     * the angle-angle partners for each angle from a particular unit *
c     * and type of molecule                                           *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,iaa,ivalue
c     --- local variables
      integer imtemp,iutemp,iatemp
c     --- storage variable
      integer ijaa0
      dimension ijaa0(NTMAX,NUMAX,MAXAA)
      save ijaa0

      if ( flag .eq. GLB_GET ) then
         ivalue = ijaa0(imolty,iunit,iaa)
      elseif ( flag .eq. GLB_SET ) then
         ijaa0(imolty,iunit,iaa) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do iatemp = 1,iaa
                  ijaa0(imtemp,iutemp,iatemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IJAA0 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ijaa1(flag,imolty,iunit,iaa,ivalue)
c     ******************************************************************
c     * the angle-angle partners for each angle from a particular unit *
c     * and type of molecule                                           *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,iaa,ivalue
c     --- local variables
      integer imtemp,iutemp,iatemp
c     --- storage variable
      integer ijaa1
      dimension ijaa1(NTMAX,NUMAX,MAXAA)
      save ijaa1

      if ( flag .eq. GLB_GET ) then
         ivalue = ijaa1(imolty,iunit,iaa)
      elseif ( flag .eq. GLB_SET ) then
         ijaa1(imolty,iunit,iaa) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do iatemp = 1,iaa
                  ijaa1(imtemp,iutemp,iatemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IJAA1 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ijaa2(flag,imolty,iunit,iaa,ivalue)
c     ******************************************************************
c     * the angle-angle partners for each angle from a particular unit *
c     * and type of molecule                                           *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,iaa,ivalue
c     --- local variables
      integer imtemp,iutemp,iatemp
c     --- storage variable
      integer ijaa2
      dimension ijaa2(NTMAX,NUMAX,MAXAA)
      save ijaa2

      if ( flag .eq. GLB_GET ) then
         ivalue = ijaa2(imolty,iunit,iaa)
      elseif ( flag .eq. GLB_SET ) then
         ijaa2(imolty,iunit,iaa) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do iatemp = 1,iaa
                  ijaa2(imtemp,iutemp,iatemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IJAA2 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ijben2(flag,imolty,iunit,ibend,ivalue)
c     ******************************************************************
c     * the bending partners for each angle from a particular unit     *
c     * and type of molecule                                           *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibend,ivalue
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer ijben2
      dimension ijben2(NTMAX,NUMAX,MAXBEND)
      save ijben2

      if ( flag .eq. GLB_GET ) then
         ivalue = ijben2(imolty,iunit,ibend)
      elseif ( flag .eq. GLB_SET ) then
         ijben2(imolty,iunit,ibend) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibend
                  ijben2(imtemp,iutemp,ibtemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IJBEN2 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ijben3(flag,imolty,iunit,ibend,ivalue)
c     ******************************************************************
c     * the bending partners for each angle from a particular unit     *
c     * and type of molecule                                           *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibend,ivalue
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer ijben3
      dimension ijben3(NTMAX,NUMAX,MAXBEND)
      save ijben3

      if ( flag .eq. GLB_GET ) then
         ivalue = ijben3(imolty,iunit,ibend)
      elseif ( flag .eq. GLB_SET ) then
         ijben3(imolty,iunit,ibend) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibend
                  ijben3(imtemp,iutemp,ibtemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IJBEN3 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ijbnbn1(flag,imolty,iunit,ibend,ivalue)
c     ******************************************************************
c     * the partners for each 1-2-3 bond-bond from a particular unit   *
c     * and type of molecule                                           *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibend,ivalue
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer ijbnbn1
      dimension ijbnbn1(NTMAX,NUMAX,MAXBEND)
      save ijbnbn1

      if ( flag .eq. GLB_GET ) then
         ivalue = ijbnbn1(imolty,iunit,ibend)
      elseif ( flag .eq. GLB_SET ) then
         ijbnbn1(imolty,iunit,ibend) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibend
                  ijbnbn1(imtemp,iutemp,ibtemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IJBNBN1 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ijbnbn2(flag,imolty,iunit,ibend,ivalue)
c     ******************************************************************
c     * the partners for each 1-2-3 bond-bond from a particular unit   *
c     * and type of molecule                                           *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibend,ivalue
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer ijbnbn2
      dimension ijbnbn2(NTMAX,NUMAX,MAXBEND)
      save ijbnbn2

      if ( flag .eq. GLB_GET ) then
         ivalue = ijbnbn2(imolty,iunit,ibend)
      elseif ( flag .eq. GLB_SET ) then
         ijbnbn2(imolty,iunit,ibend) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibend
                  ijbnbn2(imtemp,iutemp,ibtemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IJBNBN2 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ijimprop2(flag,imolty,iunit,iimp,ivalue)
c     ******************************************************************
c     * the improper partners for each improper from a particular unit *
c     * and type of molecule                                           *
c     *                                                                *
c     * originally written 03-24-2006 by M.G. Martin                   *
c     * last modified 03-24-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,iimp,ivalue
c     --- local variables
      integer imtemp,iutemp,iitemp
c     --- storage variable
      integer ijimprop2
      dimension ijimprop2(NTMAX,NUMAX,MAXIMPROP)
      save ijimprop2

      if ( flag .eq. GLB_GET ) then
         ivalue = ijimprop2(imolty,iunit,iimp)
      elseif ( flag .eq. GLB_SET ) then
         ijimprop2(imolty,iunit,iimp) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do iitemp = 1,iimp
                  ijimprop2(imtemp,iutemp,iitemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IJIMPROP2 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ijimprop3(flag,imolty,iunit,iimp,ivalue)
c     ******************************************************************
c     * the improper partners for each improper from a particular unit *
c     * and type of molecule                                           *
c     *                                                                *
c     * originally written 03-24-2006 by M.G. Martin                   *
c     * last modified 03-24-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,iimp,ivalue
c     --- local variables
      integer imtemp,iutemp,iitemp
c     --- storage variable
      integer ijimprop3
      dimension ijimprop3(NTMAX,NUMAX,MAXIMPROP)
      save ijimprop3

      if ( flag .eq. GLB_GET ) then
         ivalue = ijimprop3(imolty,iunit,iimp)
      elseif ( flag .eq. GLB_SET ) then
         ijimprop3(imolty,iunit,iimp) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do iitemp = 1,iimp
                  ijimprop3(imtemp,iutemp,iitemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IJIMPROP3 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ijimprop4(flag,imolty,iunit,iimp,ivalue)
c     ******************************************************************
c     * the improper partners for each improper from a particular unit *
c     * and type of molecule                                           *
c     *                                                                *
c     * originally written 03-24-2006 by M.G. Martin                   *
c     * last modified 03-24-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,iimp,ivalue
c     --- local variables
      integer imtemp,iutemp,iitemp
c     --- storage variable
      integer ijimprop4
      dimension ijimprop4(NTMAX,NUMAX,MAXIMPROP)
      save ijimprop4

      if ( flag .eq. GLB_GET ) then
         ivalue = ijimprop4(imolty,iunit,iimp)
      elseif ( flag .eq. GLB_SET ) then
         ijimprop4(imolty,iunit,iimp) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do iitemp = 1,iimp
                  ijimprop4(imtemp,iutemp,iitemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IJIMPROP4 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ijof5(flag,imolty,iunit,iof,ivalue)
c     ******************************************************************
c     * the partner for the one-five interaction for this unit and     *
c     * molecule type                                                  *
c     *                                                                *
c     * originally written 04-03-2006 by M.G. Martin                   *
c     * last modified 04-03-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,iof,ivalue
c     --- local variables
      integer imtemp,iutemp,iotemp
c     --- storage variable
      integer ijof5
      dimension ijof5(NTMAX,NUMAX,MAXOF)
      save ijof5

      if ( flag .eq. GLB_GET ) then
         ivalue = ijof5(imolty,iunit,iof)
      elseif ( flag .eq. GLB_SET ) then
         ijof5(imolty,iunit,iof) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do iotemp = 1,iof
                  ijof5(imtemp,iutemp,iotemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IJOF5 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ijtor2(flag,imolty,iunit,itor,ivalue)
c     ******************************************************************
c     * the partners for each regular torsion from a particular unit  *
c     * and type of molecule                                           *
c     *                                                                *
c     * originally written 03-26-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,itor,ivalue
c     --- local variables
      integer imtemp,iutemp,ittemp
c     --- storage variable
      integer ijtor2
      dimension ijtor2(NTMAX,NUMAX,MAXTOR)
      save ijtor2

      if ( flag .eq. GLB_GET ) then
         ivalue = ijtor2(imolty,iunit,itor)
      elseif ( flag .eq. GLB_SET ) then
         ijtor2(imolty,iunit,itor) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ittemp = 1,itor
                  ijtor2(imtemp,iutemp,ittemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IJTOR2 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ijtor3(flag,imolty,iunit,itor,ivalue)
c     ******************************************************************
c     * the partners for each regular torsion from a particular unit  *
c     * and type of molecule                                           *
c     *                                                                *
c     * originally written 03-26-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,itor,ivalue
c     --- local variables
      integer imtemp,iutemp,ittemp
c     --- storage variable
      integer ijtor3
      dimension ijtor3(NTMAX,NUMAX,MAXTOR)
      save ijtor3

      if ( flag .eq. GLB_GET ) then
         ivalue = ijtor3(imolty,iunit,itor)
      elseif ( flag .eq. GLB_SET ) then
         ijtor3(imolty,iunit,itor) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ittemp = 1,itor
                  ijtor3(imtemp,iutemp,ittemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IJTOR3 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ijtor4(flag,imolty,iunit,itor,ivalue)
c     ******************************************************************
c     * the partners for each regular torsion from a particular unit  *
c     * and type of molecule                                           *
c     *                                                                *
c     * originally written 03-26-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,itor,ivalue
c     --- local variables
      integer imtemp,iutemp,ittemp
c     --- storage variable
      integer ijtor4
      dimension ijtor4(NTMAX,NUMAX,MAXTOR)
      save ijtor4

      if ( flag .eq. GLB_GET ) then
         ivalue = ijtor4(imolty,iunit,itor)
      elseif ( flag .eq. GLB_SET ) then
         ijtor4(imolty,iunit,itor) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ittemp = 1,itor
                  ijtor4(imtemp,iutemp,ittemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IJTOR4 unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ijvib(flag,imolty,iunit,ibond,ivalue)
c     ******************************************************************
c     * the bond partner of each bond emanating from a particular unit *
c     * and type of molecule                                           *
c     *                                                                *
c     * originally written 03-24-2006 by M.G. Martin                   *
c     * last modified 03-24-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibond,ivalue
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer ijvib
      dimension ijvib(NTMAX,NUMAX,NNBOND)
      save ijvib

      if ( flag .eq. GLB_GET ) then
         ivalue = ijvib(imolty,iunit,ibond)
      elseif ( flag .eq. GLB_SET ) then
         ijvib(imolty,iunit,ibond) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibond
                  ijvib(imtemp,iutemp,ibtemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IJVIB unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_impbead(flag,iimp,iatom,ivalue)
c     ******************************************************************
c     * the atoms in the impropers that are used during a CBMC regrowth*
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modified 10-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,iimp,iatom,ivalue
c     --- storage variable
      integer impbead
      dimension impbead(2*MAXIMPROP,4)
      save impbead

      if ( flag .eq. GLB_GET ) then
         ivalue = impbead(iimp,iatom)
      elseif ( flag .eq. GLB_SET ) then
         impbead(iimp,iatom) = ivalue
      else
         write(6,*) 'IMPBEAD unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_impcoeff(flag,itype,icoeff,dvalue)
c     ******************************************************************
c     * the coefficients for each type of improper torsion in the      *
c     * simulation                                                     *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,icoeff
      double precision dvalue
c     --- local variables
      integer ittemp,ictemp
c     --- storage variable
      double precision impcoeff
      dimension impcoeff(TIMPMAX,MINIMPCOEFF:MAXIMPCOEFF)
      save impcoeff

      if ( flag .eq. GLB_GET ) then
         dvalue = impcoeff(itype,icoeff)
      elseif ( flag .eq. GLB_SET ) then
         impcoeff(itype,icoeff) = dvalue
      elseif ( flag .eq. GLB_INIT ) then 
         do ittemp = 1,itype
            do ictemp = MINIMPCOEFF,icoeff
               impcoeff(ittemp,ictemp) = dvalue
            enddo
         enddo
      else
         write(6,*) 'IMPCOEFF unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_impff(flag,itype,cvalue)
c     ******************************************************************
c     * the force field name of each type of improper                  *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype
      character*(*) cvalue
c     --- local variables
      integer ittemp
c     --- storage variable
      character*(FFNAMELEN) impff
      dimension impff(TIMPMAX)
      save impff

      if ( flag .eq. GLB_GET ) then
         cvalue = impff(itype)
      elseif ( flag .eq. GLB_SET ) then
         impff(itype) = cvalue
      elseif ( flag .eq. GLB_INIT ) then 
         do ittemp = 1,itype
            impff(ittemp) = cvalue
         enddo
      else
         write(6,*) 'IMPFF unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_impform(flag,itype,ivalue)
c     ******************************************************************
c     * the form of each type of improper in the simulation            *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer ittemp
c     --- storage variable
      integer impform
      dimension impform(TIMPMAX)
      save impform

      if ( flag .eq. GLB_GET ) then
         ivalue = impform(itype)
      elseif ( flag .eq. GLB_SET ) then
         impform(itype) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do ittemp = 1,itype
            impform(ittemp) = ivalue
         enddo
      else
         write(6,*) 'IMPFORM unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_impnames(flag,itype,isame,index,cvalue)
c     ******************************************************************
c     * the name of each atom in an improper torsion for each type of  *
c     * improper torsion in the system                                 *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,isame,index
      character*(*) cvalue
c     --- local variables
      integer ittemp,istemp,iitemp
c     --- storage variable
      character*10 impnames
      dimension impnames(TIMPMAX,MAXSAMEIMP,MAXIMPNAME)
      save impnames

      if ( flag .eq. GLB_GET ) then
         cvalue = impnames(itype,isame,index)
      elseif ( flag .eq. GLB_SET ) then
         impnames(itype,isame,index) = cvalue
      elseif ( flag .eq. GLB_INIT ) then 
         do ittemp = 1,itype
            do istemp = 1,isame
               do iitemp = 1,index
                  impnames(ittemp,istemp,iitemp) = cvalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'IMPNAMES unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_impropernum(flag,iimp,ivalue)
c     ******************************************************************
c     * improper type for use during CBMC growth                       *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modified 10-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,iimp,ivalue
c     --- storage variable
      integer impropernum
      dimension impropernum(2*MAXIMPROP)
      save impropernum

      if ( flag .eq. GLB_GET ) then
         ivalue = impropernum(iimp)
      elseif ( flag .eq. GLB_SET ) then
         impropernum(iimp) = ivalue
      else
         write(6,*) 'IMPROPERNUM unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_impstyle(flag,itype,ivalue)
c     ******************************************************************
c     * the style of each type of improper in the simulation           *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer ittemp
c     --- storage variable
      integer impstyle
      dimension impstyle(TIMPMAX)
      save impstyle

      if ( flag .eq. GLB_GET ) then
         ivalue = impstyle(itype)
      elseif ( flag .eq. GLB_SET ) then
         impstyle(itype) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do ittemp = 1,itype
            impstyle(ittemp) = ivalue
         enddo
      else
         write(6,*) 'IMPSTYLE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_inaa(flag,imolty,iunit,ivalue)
c     ******************************************************************
c     * the number of angle-angle terms emanating from a particular    *
c     * unit and type of molecule                                      *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ivalue
c     --- local variables
      integer imtemp,iutemp
c     --- storage variable
      integer inaa
      dimension inaa(NTMAX,NUMAX)
      save inaa

      if ( flag .eq. GLB_GET ) then
         ivalue = inaa(imolty,iunit)
      elseif ( flag .eq. GLB_SET ) then
         inaa(imolty,iunit) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               inaa(imtemp,iutemp) = ivalue
            enddo
         enddo
      else
         write(6,*) 'INAA unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_inben(flag,imolty,iunit,ivalue)
c     ******************************************************************
c     * the number of angles emanating from a particular unit and type *
c     * of molecule                                                    *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ivalue
c     --- local variables
      integer imtemp,iutemp
c     --- storage variable
      integer inben
      dimension inben(NTMAX,NUMAX)
      save inben

      if ( flag .eq. GLB_GET ) then
         ivalue = inben(imolty,iunit)
      elseif ( flag .eq. GLB_SET ) then
         inben(imolty,iunit) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               inben(imtemp,iutemp) = ivalue
            enddo
         enddo
      elseif ( flag .eq. GLB_INCR ) then
         inben(imolty,iunit) = inben(imolty,iunit) + ivalue
      else
         write(6,*) 'INBEN unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_inbnbn(flag,imolty,iunit,ivalue)
c     ******************************************************************
c     * the number of 1-2-3 bond-bond terms stored in a similar manner *
c     * as regular angles                                              *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ivalue
c     --- local variables
      integer imtemp,iutemp
c     --- storage variable
      integer inbnbn
      dimension inbnbn(NTMAX,NUMAX)
      save inbnbn

      if ( flag .eq. GLB_GET ) then
         ivalue = inbnbn(imolty,iunit)
      elseif ( flag .eq. GLB_SET ) then
         inbnbn(imolty,iunit) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               inbnbn(imtemp,iutemp) = ivalue
            enddo
         enddo
      elseif ( flag .eq. GLB_INCR ) then
         inbnbn(imolty,iunit) = inbnbn(imolty,iunit) + ivalue
      else
         write(6,*) 'INBNBN unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_inimprop(flag,imolty,iunit,ivalue)
c     ******************************************************************
c     * the number of improper partners for each improper from a       *
c     * particular unit and type of molecule                           *
c     *                                                                *
c     * originally written 03-25-2006 by M.G. Martin                   *
c     * last modified 03-25-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ivalue
c     --- local variables
      integer imtemp,iutemp
c     --- storage variable
      integer inimprop
      dimension inimprop(NTMAX,NUMAX)
      save inimprop

      if ( flag .eq. GLB_GET ) then
         ivalue = inimprop(imolty,iunit)
      elseif ( flag .eq. GLB_SET ) then
         inimprop(imolty,iunit) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               inimprop(imtemp,iutemp) = ivalue
            enddo
         enddo
      elseif ( flag .eq. GLB_DECR ) then
         inimprop(imolty,iunit) = inimprop(imolty,iunit) - ivalue
      elseif ( flag .eq. GLB_INCR ) then
         inimprop(imolty,iunit) = inimprop(imolty,iunit) + ivalue
      else
         write(6,*) 'INIMPROP unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_initboxtype(flag,cvalue)
c     ******************************************************************
c     * the method to use for creating the initial molecule structure  *
c     *                                                                *
c     * originally written 12-09-2007 by M.G. Martin                   *
c     * last modified 12-09-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*20 initboxtype
      save initboxtype

      if ( flag .eq. GLB_GET ) then
         cvalue = initboxtype
      elseif ( flag .eq. GLB_SET ) then
         initboxtype = cvalue
      else
         write(6,*) 'INITBOXTYPE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_initlattice(flag,ibox,imolty,cvalue)
c     ******************************************************************
c     * the type of lattice to use for assigning the inital molecule   *
c     * positions                                                      *
c     *                                                                *
c     * originally written 12-09-2007 by M.G. Martin                   *
c     * last modified 12-09-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox,imolty
      character*(*) cvalue
c     --- local variables
      integer ibtemp,imtemp
c     --- storage variable
      character*20 initlattice
      dimension initlattice(MAXBOX,NTMAX)
      save initlattice

      if ( flag .eq. GLB_GET ) then
         cvalue = initlattice(ibox,imolty)
      elseif ( flag .eq. GLB_SET ) then
         initlattice(ibox,imolty) = cvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 1,ibox
            do imtemp = 1,imolty
               initlattice(ibtemp,imtemp) = cvalue
            enddo
         enddo
      else
         write(6,*) 'INITLATTICE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_initmol(flag,ibox,imolty,ivalue)
c     ******************************************************************
c     * the initial number of molecules of each type in each box       *
c     *                                                                *
c     * originally written 12-09-2007 by M.G. Martin                   *
c     * last modified 12-09-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox,imolty,ivalue
c     --- local variables
      integer ibtemp,imtemp
c     --- storage variable
      integer initmol
      dimension initmol(MAXBOX,NTMAX)
      save initmol

      if ( flag .eq. GLB_GET ) then
         ivalue = initmol(ibox,imolty)
      elseif ( flag .eq. GLB_SET ) then
         initmol(ibox,imolty) = ivalue
      elseif ( flag .eq. GLB_INCR ) then
         initmol(ibox,imolty) = initmol(ibox,imolty) + ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 1,ibox
            do imtemp = 1,imolty
               initmol(ibtemp,imtemp) = ivalue
            enddo
         enddo
      else
         write(6,*) 'INITMOL unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_initstyle(flag,ibox,imolty,cvalue)
c     ******************************************************************
c     * the style to use for generating the initial molecule template  *
c     * for each molecule type in each simulation box                  *
c     *                                                                *
c     * originally written 12-09-2007 by M.G. Martin                   *
c     * last modified 12-09-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox,imolty
      character*(*) cvalue
c     --- local variables
      integer ibtemp,imtemp
c     --- storage variable
      character*20 initstyle
      dimension initstyle(MAXBOX,NTMAX)
      save initstyle

      if ( flag .eq. GLB_GET ) then
         cvalue = initstyle(ibox,imolty)
      elseif ( flag .eq. GLB_SET ) then
         initstyle(ibox,imolty) = cvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 1,ibox
            do imtemp = 1,imolty
               initstyle(ibtemp,imtemp) = cvalue
            enddo
         enddo
      else
         write(6,*) 'INITSTYLE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_inixyz(flag,ibox,ix,iy,iz)
c     ******************************************************************
c     * initial number of molecules in each box in the x,y, and z      *
c     * directions for use when setting up on a lattice                *
c     *                                                                *
c     * originally written 12-09-2007 by M.G. Martin                   *
c     * last modified 12-09-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ibox,ix,iy,iz
c     --- storage variable
      integer inix, iniy, iniz
      dimension inix(MAXBOX),iniy(MAXBOX),iniz(MAXBOX)
      save inix,iniy,iniz

      if ( flag .eq. GLB_GET ) then
         ix = inix(ibox)
         iy = iniy(ibox)
         iz = iniz(ibox)
      elseif ( flag .eq. GLB_SET ) then
         inix(ibox) = ix
         iniy(ibox) = iy
         iniz(ibox) = iz
      else
         write(6,*) 'INIXYZ unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_inof(flag,imolty,iunit,ivalue)
c     ******************************************************************
c     * the number of special one-five interactions for each unit and  *
c     * molecule type.  This is automatically determined by the code   *
c     *                                                                *
c     * originally written 04-03-2006 by M.G. Martin                   *
c     * last modified 04-03-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ivalue
c     --- local variables
      integer imtemp,iutemp
c     --- storage variable
      integer inof
      dimension inof(NTMAX,NUMAX)
      save inof

      if ( flag .eq. GLB_GET ) then
         ivalue = inof(imolty,iunit)
      elseif ( flag .eq. GLB_SET ) then
         inof(imolty,iunit) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               inof(imtemp,iutemp) = ivalue
            enddo
         enddo
      elseif ( flag .eq. GLB_INCR ) then
         inof(imolty,iunit) = inof(imolty,iunit) + ivalue
      else
         write(6,*) 'INOF unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_inring(flag,imolty,iunit,ivalue)
c     ******************************************************************
c     * the ring number for each atom in each molecule type, this is   *
c     * set to 0 if the atom is not in a ring                          *
c     *                                                                *
c     * originally written 12-15-2007 by M.G. Martin                   *
c     * last modified 12-15-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ivalue
c     --- local variables
      integer imtemp,iutemp
c     --- storage variable
      integer inring
      dimension inring(NTMAX,NUMAX)
      save inring

      if ( flag .eq. GLB_GET ) then
         ivalue = inring(imolty,iunit)
      elseif ( flag .eq. GLB_SET ) then
         inring(imolty,iunit) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               inring(imtemp,iutemp) = ivalue
            enddo
         enddo
      else
         write(6,*) 'INRING unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_interpolatestyle(flag,cvalue)
c     ******************************************************************
c     * the style used to interpolate tabular potentials               *
c     *                                                                *
c     * originally written 10-03-2008 by M.G. Martin                   *
c     * last modified 10-03-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*20 interpolatestyle
      save interpolatestyle

      if ( flag .eq. GLB_GET ) then
         cvalue = interpolatestyle
      elseif ( flag .eq. GLB_SET ) then
         interpolatestyle = cvalue
      else
         write(6,*) 'INTERPOLATESTYLE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_intor(flag,imolty,iunit,ivalue)
c     ******************************************************************
c     * the number of partners for each regular torsion from a         *
c     * particular unit and type of molecule                           *
c     *                                                                *
c     * originally written 03-26-2006 by M.G. Martin                   *
c     * last modified 03-26-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ivalue
c     --- local variables
      integer imtemp,iutemp
c     --- storage variable
      integer intor
      dimension intor(NTMAX,NUMAX)
      save intor

      if ( flag .eq. GLB_GET ) then
         ivalue = intor(imolty,iunit)
      elseif ( flag .eq. GLB_SET ) then
         intor(imolty,iunit) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               intor(imtemp,iutemp) = ivalue
            enddo
         enddo
      elseif ( flag .eq. GLB_DECR ) then
         intor(imolty,iunit) = intor(imolty,iunit) - ivalue
      elseif ( flag .eq. GLB_INCR ) then
         intor(imolty,iunit) = intor(imolty,iunit) + ivalue
      else
         write(6,*) 'INTOR unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_invib(flag,imolty,iunit,ivalue)
c     ******************************************************************
c     * the number of bonds emanating from a particular unit and       *
c     * type of molecule                                               *
c     *                                                                *
c     * originally written 03-24-2006 by M.G. Martin                   *
c     * last modified 03-24-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ivalue
c     --- local variables
      integer imtemp,iutemp
c     --- storage variable
      integer invib
      dimension invib(NTMAX,NUMAX)
      save invib

      if ( flag .eq. GLB_GET ) then
         ivalue = invib(imolty,iunit)
      elseif ( flag .eq. GLB_SET ) then
         invib(imolty,iunit) = ivalue
      elseif ( flag .eq. GLB_INCR ) then
         invib(imolty,iunit) = invib(imolty,iunit) + ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               invib(imtemp,iutemp) = ivalue
            enddo
         enddo
      else
         write(6,*) 'INVIB unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_io_directory(flag,cvalue)
c     ******************************************************************
c     * the name of the directory to use for input and output          *
c     *                                                                *
c     * originally written 04-26-2006 by M.G. Martin                   *
c     * last modified 04-26-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage variable
      character*(MAXDIRLENGTH) io_directory
      save io_directory

      if ( flag .eq. GLB_GET ) then
         cvalue = io_directory
      elseif ( flag .eq. GLB_SET ) then
         io_directory = cvalue
      else
         write(6,*) 'IO_DIRECTORY unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_isolvtype(flag,ivalue)
c     ******************************************************************
c     * the integer code for the solvation type                        *
c     *                                                                *
c     * originally written 01-12-2007 by M.G. Martin                   *
c     * last modified 01-12-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer isolvtype
      save isolvtype

      if ( flag .eq. GLB_SET ) then
         isolvtype = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = isolvtype
      else
         write(6,*) 'ISOLVTYPE: unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_itaa(flag,imolty,iunit,iaa,ivalue)
c     ******************************************************************
c     * the type of each angle-angle emanating from a particular unit  *
c     * and type of molecule                                           *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,iaa,ivalue
c     --- local variables
      integer imtemp,iutemp,iatemp
c     --- storage variable
      integer itaa
      dimension itaa(NTMAX,NUMAX,MAXAA)
      save itaa

      if ( flag .eq. GLB_GET ) then
         ivalue = itaa(imolty,iunit,iaa)
      elseif ( flag .eq. GLB_SET ) then
         itaa(imolty,iunit,iaa) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do iatemp = 1,iaa
                  itaa(imtemp,iutemp,iatemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'ITAA unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_itben(flag,imolty,iunit,ibend,ivalue)
c     ******************************************************************
c     * the type of each angle emanating from a particular unit and    *
c     * type of molecule                                               *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibend,ivalue
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer itben
      dimension itben(NTMAX,NUMAX,MAXBEND)
      save itben

      if ( flag .eq. GLB_GET ) then
         ivalue = itben(imolty,iunit,ibend)
      elseif ( flag .eq. GLB_SET ) then
         itben(imolty,iunit,ibend) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibend
                  itben(imtemp,iutemp,ibtemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'ITBEN unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_itbnbn(flag,imolty,iunit,ibend,ivalue)
c     ******************************************************************
c     * the type of each bond-bond emanating from a particular unit and*
c     * type of molecule                                               *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-31-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibend,ivalue
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer itbnbn
      dimension itbnbn(NTMAX,NUMAX,MAXBEND)
      save itbnbn

      if ( flag .eq. GLB_GET ) then
         ivalue = itbnbn(imolty,iunit,ibend)
      elseif ( flag .eq. GLB_SET ) then
         itbnbn(imolty,iunit,ibend) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibend
                  itbnbn(imtemp,iutemp,ibtemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'ITBNBN unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_itimprop(flag,imolty,iunit,iimp,ivalue)
c     ******************************************************************
c     * the type of each improper emanating from a particular unit and *
c     * type of molecule                                               *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,iimp,ivalue
c     --- local variables
      integer imtemp,iutemp,iitemp
c     --- storage variable
      integer itimprop
      dimension itimprop(NTMAX,NUMAX,MAXIMPROP)
      save itimprop

      if ( flag .eq. GLB_GET ) then
         ivalue = itimprop(imolty,iunit,iimp)
      elseif ( flag .eq. GLB_SET ) then
         itimprop(imolty,iunit,iimp) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do iitemp = 1,iimp
                  itimprop(imtemp,iutemp,iitemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'ITIMPROP unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_itscale(flag,imolty,iunit,itor,dvalue)
c     ******************************************************************
c     * the scaling factor to apply to the 1-4 interactions in a       *
c     * dihedral.  This is computed automatically by the code by       *
c     * counting the number of times a 1-4 interaction can occur       *
c     * (normally 1, but different for cyclic systems) and applying any*
c     * force-field dependent scaling factor                           *
c     *                                                                *
c     * originally written 04-10-2006 by M.G. Martin                   *
c     * last modified 04-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,itor
      double precision dvalue
c     --- local variables
      integer imtemp,iutemp,ittemp
c     --- storage variable
      double precision itscale
      dimension itscale(NTMAX,NUMAX,MAXTOR)
      save itscale

      if ( flag .eq. GLB_GET ) then
         dvalue = itscale(imolty,iunit,itor)
      elseif ( flag .eq. GLB_SET ) then
         itscale(imolty,iunit,itor) = dvalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ittemp = 1,itor
                  itscale(imtemp,iutemp,ittemp) = dvalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'ITSCALE unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_itof(flag,imolty,iunit,iof,ivalue)
c     ******************************************************************
c     * the type of each one-five emanating from a particular unit and *
c     * type of molecule                                               *
c     *                                                                *
c     * originally written 04-03-2006 by M.G. Martin                   *
c     * last modified 04-03-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,iof,ivalue
c     --- local variables
      integer imtemp,iutemp,iotemp
c     --- storage variable
      integer itof
      dimension itof(NTMAX,NUMAX,MAXOF)
      save itof

      if ( flag .eq. GLB_GET ) then
         ivalue = itof(imolty,iunit,iof)
      elseif ( flag .eq. GLB_SET ) then
         itof(imolty,iunit,iof) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do iotemp = 1,iof
                  itof(imtemp,iutemp,iotemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'ITOF unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_ittor(flag,imolty,iunit,itor,ivalue)
c     ******************************************************************
c     * the type of each torsion emanating from a particular unit and  *
c     * type of molecule                                               *
c     *                                                                *
c     * originally written 03-26-2006 by M.G. Martin                   *
c     * last modified 03-26-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,itor,ivalue
c     --- local variables
      integer imtemp,iutemp,ittemp
c     --- storage variable
      integer ittor
      dimension ittor(NTMAX,NUMAX,MAXTOR)
      save ittor

      if ( flag .eq. GLB_GET ) then
         ivalue = ittor(imolty,iunit,itor)
      elseif ( flag .eq. GLB_SET ) then
         ittor(imolty,iunit,itor) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ittemp = 1,itor
                  ittor(imtemp,iutemp,ittemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'ITTOR unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_itvib(flag,imolty,iunit,ibond,ivalue)
c     ******************************************************************
c     * the type of each bond emanating from a particular unit and     *
c     * type of molecule                                               *
c     *                                                                *
c     * originally written 03-24-2006 by M.G. Martin                   *
c     * last modified 03-24-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibond,ivalue
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer itvib
      dimension itvib(NTMAX,NUMAX,NNBOND)
      save itvib

      if ( flag .eq. GLB_GET ) then
         ivalue = itvib(imolty,iunit,ibond)
      elseif ( flag .eq. GLB_SET ) then
         itvib(imolty,iunit,ibond) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibond
                  itvib(imtemp,iutemp,ibtemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'ITVIB unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_kim_model(flag,cvalue)
c     ******************************************************************
c     * character string for the kim model                             *
c     *                                                                *
c     * originally written 12-02-2011 by M.G. Martin                   *
c     * last modified 12-02-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      character*(*) cvalue
c     --- information variables
      character*(50) kim_model
      save kim_model

      if ( flag .eq. GLB_SET ) then
         kim_model = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = kim_model
      else
         write(6,*) 'KIM_MODEL: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lbackbone(flag,imolty,iunit,lvalue)
c     ******************************************************************
c     * logical that states whether an atom is along a protein backbone*
c     *                                                                *
c     * originally written 10-03-2008 by M.G. Martin                   *
c     * last modified 10-03-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,iunit
      logical lvalue
c     --- information variables
      logical lbackbone
      dimension lbackbone(NTMAX,NUMAX)
      save lbackbone

      if ( flag .eq. GLB_SET ) then
         lbackbone(imolty,iunit) = lvalue
      elseif ( flag .eq. GLB_GET ) then
         lvalue = lbackbone(imolty,iunit)
      else
         write(6,*) 'LBACKBONE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lbonang(flag,ibtype,lvalue)
c     ******************************************************************
c     * bond-angle interaction logical for each angle type             *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ibtype
      logical lvalue
c     --- local variables
      integer itemp
c     --- information variables
      logical lbonang
      dimension lbonang(TBENMAX)
      save lbonang

      if ( flag .eq. GLB_SET ) then
         lbonang(ibtype) = lvalue
      elseif ( flag .eq. GLB_GET ) then
         lvalue = lbonang(ibtype)
      elseif ( flag .eq. GLB_INIT ) then
         do itemp = 1,ibtype
            lbonang(itemp) = lvalue
         enddo
      else
         write(6,*) 'LBONANG: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lbondgraph(flag,imolty,lvalue)
c     ******************************************************************
c     * logical that is true if the data structure for this molecule   *
c     * type is a fully bonded graph, and therefore a valid molecule   *
c     * false if it it not fully bonded and just a collection of atoms *
c     * this can occur for pretend molecules that are the porous       *
c     * material in an adsorption study                                *
c     *                                                                *
c     * originally written 07-15-2013 by M.G. Martin                   *
c     * last modified 07-15-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      logical lvalue
c     --- information variables
      logical lbondgraph
      dimension lbondgraph(NTMAX)
      save lbondgraph

      if ( flag .eq. GLB_GET ) then
         lvalue = lbondgraph(imolty)
      elseif ( flag .eq. GLB_SET ) then
         lbondgraph(imolty) = lvalue
      else
         write(6,*) 'LBONDGRAPH: invalid flag:',flag
         stop
      endif

      return
      end



      subroutine twh_lbonbon(flag,ibtype,lvalue)
c     ******************************************************************
c     * bond-bond interaction logical for each angle type              *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ibtype
      logical lvalue
c     --- local variables
      integer itemp
c     --- information variables
      logical lbonbon
      dimension lbonbon(TBENMAX)
      save lbonbon

      if ( flag .eq. GLB_SET ) then
         lbonbon(ibtype) = lvalue
      elseif ( flag .eq. GLB_GET ) then
         lvalue = lbonbon(ibtype)
      elseif ( flag .eq. GLB_INIT ) then
         do itemp = 1,ibtype
            lbonbon(itemp) = lvalue
         enddo
      else
         write(6,*) 'LBONBON: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lcao_atomtypes(flag,ivalue)
c     ******************************************************************
c     * number of LCAO atomtypes                                       *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer lcao_atomtypes
      save lcao_atomtypes

      if ( flag .eq. GLB_GET ) then
         ivalue = lcao_atomtypes
      elseif ( flag .eq. GLB_SET ) then
         lcao_atomtypes = ivalue
      else
         write(6,*) 'LCAO_ATOMTYPES: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lcao_filename(flag,itype,cvalue)
c     ******************************************************************
c     * character string used for files with the LCAO external         *
c     * potential                                                      *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype
      character*(*) cvalue
c     --- information variables
      character*50 lcao_filename
      dimension lcao_filename(NNTYPE)
      save lcao_filename

      if ( flag .eq. GLB_GET ) then
         cvalue = lcao_filename(itype)
      elseif ( flag .eq. GLB_SET ) then
         lcao_filename(itype) = cvalue
      else
         write(6,*) 'LCAO_FILENAME: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lcao_functional(flag,cvalue)
c     ******************************************************************
c     * character string used with the LCAO external potential         *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      character*(*) cvalue
c     --- information variables
      character*20 lcao_functional
      save lcao_functional

      if ( flag .eq. GLB_GET ) then
         cvalue = lcao_functional
      elseif ( flag .eq. GLB_SET ) then
         lcao_functional = cvalue
      else
         write(6,*) 'LCAO_FUNCTIONAL: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lcao_gridmultiplier(flag,dvalue)
c     ******************************************************************
c     * LCAO gridmultiplier variable for external potential            *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision lcao_gridmultiplier
      save lcao_gridmultiplier

      if ( flag .eq. GLB_GET ) then
         dvalue = lcao_gridmultiplier
      elseif ( flag .eq. GLB_SET ) then
         lcao_gridmultiplier = dvalue
      else
         write(6,*) 'LCAO_GRIDMULTIPLIER: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lcao_kgridproduct(flag,dvalue)
c     ******************************************************************
c     * LCAO kgridproduct variable for external potential              *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision lcao_kgridproduct
      save lcao_kgridproduct

      if ( flag .eq. GLB_GET ) then
         dvalue = lcao_kgridproduct
      elseif ( flag .eq. GLB_SET ) then
         lcao_kgridproduct = dvalue
      else
         write(6,*) 'LCAO_KGRIDPRODUCT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lcoulomb(flag,lvalue)
c     ******************************************************************
c     * logical for the use of Coulombic pair potential                *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modified 10-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical lcoulomb
      save lcoulomb

      if ( flag .eq. GLB_GET ) then
         lvalue = lcoulomb
      elseif ( flag .eq. GLB_SET ) then
         lcoulomb = lvalue
      else
         write(6,*) 'LCOULOMB: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lcreatemap(flag,lvalue)
c     ******************************************************************
c     * energy map system logical                                      *
c     *                                                                *
c     * originally written 11-02-2007 by M.G. Martin                   *
c     * last modified 11-02-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical lcreatemap
      save lcreatemap

      if ( flag .eq. GLB_GET ) then
         lvalue = lcreatemap
      elseif ( flag .eq. GLB_SET ) then
         lcreatemap = lvalue
      else
         write(6,*) 'LCREATEMAP: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_ldualcutoff(flag,lvalue)
c     ******************************************************************
c     * logical specifying whether to use the dual cutoff algorithm    *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical ldualcutoff
      save ldualcutoff

      if ( flag .eq. GLB_GET ) then
         lvalue = ldualcutoff
      elseif ( flag .eq. GLB_SET ) then
         ldualcutoff = lvalue
      else
         write(6,*) 'LDUALCUTOFF: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lexsched(flag,iunit,lvalue)
c     ******************************************************************
c     * scheduling logical for use during CBMC growths                 *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 12-04-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,iunit
      logical lvalue
c     --- local variables
      integer ivalue

c     --- convert lvalue to integer code
      if ( flag .eq. GLB_SET .or. flag .eq. GLB_INIT ) then
         if ( lvalue ) then
            ivalue = LOG_TRUE
         else
            ivalue = LOG_FALSE
         endif
      endif

      call twh_logical_exsched(flag,iunit,ivalue)

      if ( flag .eq. GLB_GET ) then
c        --- convert ivalue to lvalue
         if ( ivalue .eq. LOG_TRUE ) then
            lvalue = .true.
         else
            lvalue = .false.
         endif
      endif

      return
      end

      subroutine twh_lewald(flag,lvalue)
c     ******************************************************************
c     * logical for whether to use the ewald sum for Coulombic terms   *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modified 10-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical lewald
      save lewald

      if ( flag .eq. GLB_GET ) then
         lvalue = lewald
      elseif ( flag .eq. GLB_SET ) then
         lewald = lvalue
      else
         write(6,*) 'LEWALD: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lexist(flag,iunit,lvalue)
c     ******************************************************************
c     * temporary variable for keeping track of which units have been  *
c     * used for some purpose: in CBMC this is the beads that already  *
c     * exist                                                          *
c     *                                                                *
c     * originally written 12-14-2007 by M.G. Martin                   *
c     * last modified 12-04-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,iunit
      logical lvalue
c     --- local variables
      integer ivalue

c     --- convert lvalue to integer code
      if ( flag .eq. GLB_SET .or. flag .eq. GLB_INIT ) then
         if ( lvalue ) then
            ivalue = LOG_TRUE
         else
            ivalue = LOG_FALSE
         endif
      endif

      call twh_logical_exist(flag,iunit,ivalue)

      if ( flag .eq. GLB_GET ) then
c        --- convert ivalue to lvalue
         if ( ivalue .eq. LOG_TRUE ) then
            lvalue = .true.
         else
            lvalue = .false.
         endif
      endif

      return
      end

      subroutine twh_linit(flag,lvalue)
c     ******************************************************************
c     * system initialization logical                                  *
c     *                                                                *
c     * originally written 12-09-2007 by M.G. Martin                   *
c     * last modified 12-09-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical linit
      save linit

      if ( flag .eq. GLB_GET ) then
         lvalue = linit
      elseif ( flag .eq. GLB_SET ) then
         linit = lvalue
      else
         write(6,*) 'LINIT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lmoveme(flag,ichain,lvalue)
c     ******************************************************************
c     * logical for moving atoms in the atomshift family of moves      *
c     *                                                                *
c     * originally written 06-24-2008 by M.G. Martin                   *
c     * last modified 06-24-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ichain
      logical lvalue
c     --- local variables
      integer ivalue

c     --- convert lvalue to integer code
      if ( flag .eq. GLB_SET .or. flag .eq. GLB_INIT ) then
         if ( lvalue ) then
            ivalue = LOG_TRUE
         else
            ivalue = LOG_FALSE
         endif
      endif

      call twh_logical_moveme(flag,ichain,ivalue)

      if ( flag .eq. GLB_GET ) then
c        --- convert ivalue to lvalue
         if ( ivalue .eq. LOG_TRUE ) then
            lvalue = .true.
         else
            lvalue = .false.
         endif
      endif

      return
      end

      subroutine twh_loftor(flag,ittype,lvalue)
c     ******************************************************************
c     * one-four nonbond interaction logical for each torsion type     *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ittype
      logical lvalue
c     --- local variables
      integer itemp
c     --- information variables
      logical loftor
      dimension loftor(TTORMAX)
      save loftor

      if ( flag .eq. GLB_SET ) then
         loftor(ittype) = lvalue
      elseif ( flag .eq. GLB_GET ) then
         lvalue = loftor(ittype)
      elseif ( flag .eq. GLB_INIT ) then
         do itemp = 1,ittype
            loftor(itemp) = lvalue
         enddo
      else
         write(6,*) 'LOFTOR: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_loutdft(flag,lvalue)
c     ******************************************************************
c     * logical that controls the output in Tramonto DFT format        *
c     *                                                                *
c     * originally written 08-10-2008 by M.G. Martin                   *
c     * last modified 08-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical loutdft
      save loutdft

      if ( flag .eq. GLB_SET ) then
         loutdft = lvalue
      elseif ( flag .eq. GLB_GET ) then
         lvalue = loutdft
      else
         write(6,*) 'LOUTDFT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_loutdlpoly(flag,lvalue)
c     ******************************************************************
c     * logical that controls the output in DLPOLY format              *
c     *                                                                *
c     * originally written 08-10-2008 by M.G. Martin                   *
c     * last modified 08-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical loutdlpoly
      save loutdlpoly

      if ( flag .eq. GLB_SET ) then
         loutdlpoly = lvalue
      elseif ( flag .eq. GLB_GET ) then
         lvalue = loutdlpoly
      else
         write(6,*) 'LOUTDLPOLY: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_louthist(flag,lvalue)
c     ******************************************************************
c     * logical that controls the output of histogram information      *
c     *                                                                *
c     * originally written 08-10-2008 by M.G. Martin                   *
c     * last modified 08-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical louthist
      save louthist

      if ( flag .eq. GLB_SET ) then
         louthist = lvalue
      elseif ( flag .eq. GLB_GET ) then
         lvalue = louthist
      else
         write(6,*) 'LOUTHIST: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_loutlammps(flag,lvalue)
c     ******************************************************************
c     * logical that controls the output of LAMMPS format information  *
c     *                                                                *
c     * originally written 08-10-2008 by M.G. Martin                   *
c     * last modified 08-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical loutlammps
      save loutlammps

      if ( flag .eq. GLB_SET ) then
         loutlammps = lvalue
      elseif ( flag .eq. GLB_GET ) then
         lvalue = loutlammps
      else
         write(6,*) 'LOUTLAMMPS: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_loutput(flag,lvalue)
c     ******************************************************************
c     * logical as to whether to output information to std output      *
c     *                                                                *
c     * originally written 10-06-2008 by M.G. Martin                   *
c     * last modified 10-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical loutput
      save loutput

      if ( flag .eq. GLB_GET ) then
         lvalue = loutput
      elseif ( flag .eq. GLB_SET ) then
         loutput = lvalue
      else
         write(6,*) 'LOUTPUT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lperiodic(flag,ichain,lvalue)
c     ******************************************************************
c     * logical as to whether a molecule currently wraps through the   *
c     * periodic boundaries in a manner that makes it infinite         *
c     *                                                                *
c     * originally written 12-15-2007 by M.G. Martin                   *
c     * last modified 12-15-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ichain
      logical lvalue
c     --- local variables
      integer ivalue

c     --- convert lvalue to integer code
      if ( flag .eq. GLB_SET .or. flag .eq. GLB_INIT ) then
         if ( lvalue ) then
            ivalue = LOG_TRUE
         else
            ivalue = LOG_FALSE
         endif
      endif

      call twh_logical_periodic(flag,ichain,ivalue)

      if ( flag .eq. GLB_GET ) then
c        --- convert ivalue to lvalue
         if ( ivalue .eq. LOG_TRUE ) then
            lvalue = .true.
         else
            lvalue = .false.
         endif
      endif

      return
      end

      subroutine twh_lreadwrapper(flag,lvalue)
c     ******************************************************************
c     * lreadwrapper determines whether the shared read data is used   *
c     * instead of parsed data.  Parameter passed into towheemainloop()*
c     *                                                                *
c     * originally written 06-25-2009 by M.G. Martin                   *
c     * last modified 06-25-2009 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical lreadwrapper
      save lreadwrapper

      if ( flag .eq. GLB_GET ) then
         lvalue = lreadwrapper
      elseif ( flag .eq. GLB_SET ) then
         lreadwrapper = lvalue
      else
         write(6,*) 'LREADWRAPPER: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lrectbox(flag,ibox,lvalue)
c     ******************************************************************
c     * logical for whether each simulation box is rectilinear         *
c     *                                                                *
c     * originally written 09-25-2008 by M.G. Martin                   *
c     * last modified 09-25-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ibox
      logical lvalue
c     --- information variables
      logical lrectbox
      dimension lrectbox(MAXBOX)
      save lrectbox

      if ( flag .eq. GLB_GET ) then
         lvalue = lrectbox(ibox)
      elseif ( flag .eq. GLB_SET ) then
         lrectbox(ibox) = lvalue
      else
         write(6,*) 'LRECTBOX: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lscaleatom(flag,itype,lvalue)
c     ******************************************************************
c     * logical for each type of atom related to scaling energies for  *
c     * the lambda scaling options                                     *
c     *                                                                *
c     * originally written 12-19-2007 by M.G. Martin                   *
c     * last modified 12-18-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype
      logical lvalue
c     --- local variables
      integer ittemp
c     --- information variables
      logical lscaleatom
      dimension lscaleatom(NNTYPE)
      save lscaleatom  

      if ( flag .eq. GLB_GET ) then
         lvalue = lscaleatom(itype)
      elseif ( flag .eq. GLB_SET ) then
         lscaleatom(itype) = lvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            lscaleatom(ittemp) = lvalue
         enddo
      else
         write(6,*) 'LSCALEATOM: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lshift(flag,lvalue)
c     ******************************************************************
c     * logical for whether the vdw potential is shifted to be zero at *
c     * the cutoff                                                     *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modified 10-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical lshift
      save lshift

      if ( flag .eq. GLB_GET ) then
         lvalue = lshift
      elseif ( flag .eq. GLB_SET ) then
         lshift = lvalue
      else
         write(6,*) 'LSHIFT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_ltailc(flag,lvalue)
c     ******************************************************************
c     * logical for whether to use long range tail corrections for the *
c     * vdw potential                                                  *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modified 10-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical ltailc
      save ltailc

      if ( flag .eq. GLB_GET ) then
         lvalue = ltailc
      elseif ( flag .eq. GLB_SET ) then
         ltailc = lvalue
      else
         write(6,*) 'LTAILC: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_luse_c(flag,imolty,iunit,lvalue)
c     ******************************************************************
c     * logical for whether a particular atom has coulombic            *
c     * interactions                                                   *
c     *                                                                *
c     * originally written 09-25-2008 by M.G. Martin                   *
c     * last modified 09-25-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,iunit
      logical lvalue
c     --- local variables
      integer imtemp,iutemp
c     --- information variables
      logical luse_c
      dimension luse_c(NTMAX,NUMAX)
      save luse_c

      if ( flag .eq. GLB_GET ) then
         lvalue = luse_c(imolty,iunit)
      elseif ( flag .eq. GLB_SET ) then
         luse_c(imolty,iunit) = lvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               luse_c(imtemp,iutemp) = lvalue
            enddo
         enddo
      else
         write(6,*) 'LUSE_C: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_luse_v(flag,imolty,iunit,lvalue)
c     ******************************************************************
c     * logical for whether a particular atom has van der Waals        *
c     * interactions                                                   *
c     *                                                                *
c     * originally written 09-25-2008 by M.G. Martin                   *
c     * last modified 09-25-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,iunit
      logical lvalue
c     --- local variables
      integer imtemp,iutemp
c     --- information variables
      logical luse_v
      dimension luse_v(NTMAX,NUMAX)
      save luse_v

      if ( flag .eq. GLB_GET ) then
         lvalue = luse_v(imolty,iunit)
      elseif ( flag .eq. GLB_SET ) then
         luse_v(imolty,iunit) = lvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               luse_v(imtemp,iutemp) = lvalue
            enddo
         enddo
      else
         write(6,*) 'LUSE_V: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_lusecom(flag,lvalue)
c     ******************************************************************
c     * logical for whether we are using the additional center-of-mass *
c     * molecule cutoff as a cost-saving method                        *
c     *                                                                *
c     * originally written 05-03-2008 by M.G. Martin                   *
c     * last modified 05-03-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical lusecom
      save lusecom

      if ( flag .eq. GLB_GET ) then
         lvalue = lusecom
      elseif ( flag .eq. GLB_SET ) then
         lusecom = lvalue
      else
         write(6,*) 'LUSECOM: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_macro_variable(flag,ivalue)
c     ******************************************************************
c     * TMMC macrovariable                                             *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer macro_variable
      save macro_variable

      if ( flag .eq. GLB_GET ) then
         ivalue = macro_variable
      elseif ( flag .eq. GLB_SET ) then
         macro_variable = ivalue
      elseif ( flag .eq. GLB_INCR ) then
         macro_variable = macro_variable + ivalue
      else
         write(6,*) 'MACRO_VARIABLE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_mass(flag,itype,dvalue)
c     ******************************************************************
c     * the mass for each type of atom, or pseudo-atom, in the system  *
c     *                                                                *
c     * originally written 04-12-2006 by M.G. Martin                   *
c     * last modified 04-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype
      double precision dvalue
c     --- local variables
      integer ittemp
c     --- information variables
      double precision mass
      dimension mass(NNTYPE)
      save mass

      if ( flag .eq. GLB_SET ) then
         mass(itype) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = mass(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            mass(ittemp) = dvalue
         enddo
      else
         write(6,*) 'MASS: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_masstotal(flag,imolty,dvalue)
c     ******************************************************************
c     * the total mass for molecule type                               *
c     *                                                                *
c     * originally written 10-02-2008 by M.G. Martin                   *
c     * last modified 10-02-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- information variables
      double precision masstotal
      dimension masstotal(NTMAX)
      save masstotal

      if ( flag .eq. GLB_SET ) then
         masstotal(imolty) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = masstotal(imolty)
      else
         write(6,*) 'MASSTOTAL: invalid global flag',flag
         stop
      endif

      return
      end


      subroutine twh_mapmolty(flag,ivalue)
c     ******************************************************************
c     * the molecule type to use when constructing the energy bias map *
c     *                                                                *
c     * originally written 11-04-2007 by M.G. Martin                   *
c     * last modified 11-04-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer mapmolty
      save mapmolty

      if ( flag .eq. GLB_GET ) then
         ivalue = mapmolty
      elseif ( flag .eq. GLB_SET ) then
         mapmolty = ivalue
      else
         write(6,*) 'MAPMOLTY: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_max_bond_length(flag,dvalue)
c     ******************************************************************
c     * the maximum allowed bond length for all bonds in the system    *
c     *                                                                *
c     * originally written 10-26-2008 by M.G. Martin                   *
c     * last modified 10-26-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision max_bond_length
      save max_bond_length

      if ( flag .eq. GLB_GET ) then
         dvalue = max_bond_length
      elseif ( flag .eq. GLB_SET ) then
         max_bond_length = dvalue
      else
         write(6,*) 'MAX_BOND_LENGTH: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_minboxdim(flag,dvalue)
c     ******************************************************************
c     * the minimum allowed box dimension                              *
c     *                                                                *
c     * originally written 10-06-2008 by M.G. Martin                   *
c     * last modified 10-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision minboxdim
      save minboxdim

      if ( flag .eq. GLB_GET ) then
         dvalue = minboxdim
      elseif ( flag .eq. GLB_SET ) then
         minboxdim = dvalue
      else
         write(6,*) 'MINBOXDIM: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_mixrule_adjust_key(flag,iadjust,cvalue)
c     ******************************************************************
c     * the mixrule adjust key type for each adjust term               *
c     *                                                                *
c     * originally written 04-27-2006 by M.G. Martin                   *
c     * last modified 04-27-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,iadjust
      character*(*) cvalue
c     --- local variables
      integer iatemp
c     --- information variables
      character*(FFNAMELEN) mixrule_adjust_key
      dimension mixrule_adjust_key(MAXMIXADJUST)
      save mixrule_adjust_key

      if ( flag .eq. GLB_SET ) then
         mixrule_adjust_key(iadjust) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = mixrule_adjust_key(iadjust)
      elseif ( flag .eq. GLB_INIT ) then
         do iatemp = 1,iadjust
            mixrule_adjust_key(iatemp) = cvalue
         enddo
      else
         write(6,*) 'MIXRULE_ADJUST_KEY: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_mixrule_adjust_keynames(flag,iadjust,ikey,cvalue)
c     ******************************************************************
c     * the mixrule adjust key type for each adjust term               *
c     *                                                                *
c     * originally written 04-27-2006 by M.G. Martin                   *
c     * last modified 04-27-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,iadjust,ikey
      character*(*) cvalue
c     --- local variables
      integer iatemp,iktemp
c     --- information variables
      character*(FFNAMELEN) mixrule_adjust_keynames
      dimension mixrule_adjust_keynames(MAXMIXADJUST,2)
      save mixrule_adjust_keynames

      if ( flag .eq. GLB_SET ) then
         mixrule_adjust_keynames(iadjust,ikey) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = mixrule_adjust_keynames(iadjust,ikey)
      elseif ( flag .eq. GLB_INIT ) then
         do iatemp = 1,iadjust
            do iktemp = 1,ikey
               mixrule_adjust_keynames(iatemp,iktemp) = cvalue
            enddo
         enddo
      else
         write(6,*) 'MIXRULE_ADJUST_KEYNAMES: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_mixrule_adjust_total(flag,ivalue)
c     ******************************************************************
c     * the total number of mixrule adjust terms                       *
c     *                                                                *
c     * originally written 04-27-2006 by M.G. Martin                   *
c     * last modified 04-21-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer mixrule_adjust_total
      save mixrule_adjust_total

      if ( flag .eq. GLB_SET ) then
c        --- safety check - make sure this is less than the maximum
         if ( ivalue .gt. MAXMIXADJUST ) then
            write(6,*) 'MIXRULE_ADJUST_TOTAL: attempted to set an'
     &           ,' invalid value'
            write(6,*) 'requested value:',ivalue
            write(6,*) 'MAXMIXADJUST (see preproc.h):',MAXMIXADJUST
            stop
         endif
         mixrule_adjust_total = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = mixrule_adjust_total
      else
         write(6,*) 'MIXRULE_ADJUST_TOTAL: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_mixrule_adjustments(flag,iadjust,icoeff,dvalue)
c     ******************************************************************
c     * the mixrule adjustments each adjust term                       *
c     *                                                                *
c     * originally written 04-27-2006 by M.G. Martin                   *
c     * last modified 04-27-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,iadjust,icoeff
      double precision dvalue
c     --- local variables
      integer iatemp,ictemp
c     --- information variables
      double precision mixrule_adjust_value
      dimension mixrule_adjust_value(MAXMIXADJUST,MAXNBCOEFF)
      save mixrule_adjust_value

      if ( flag .eq. GLB_SET ) then
         mixrule_adjust_value(iadjust,icoeff) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = mixrule_adjust_value(iadjust,icoeff)
      elseif ( flag .eq. GLB_INIT ) then
         do iatemp = 1,iadjust
            do ictemp = 1,icoeff
               mixrule_adjust_value(iatemp,ictemp) = dvalue
            enddo
         enddo
      else
         write(6,*) 'MIXRULE_ADJUSTMENTS: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_moviefreq(flag,ivalue)
c     ******************************************************************
c     * frequency for outputing movie files during a simulation        *
c     *                                                                *
c     * originally written 12-09-2007 by M.G. Martin                   *
c     * last modified 12-07-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer moviefreq
      save moviefreq

      if ( flag .eq. GLB_GET ) then
         ivalue = moviefreq
      elseif ( flag .eq. GLB_SET ) then
         moviefreq = ivalue
      else
         write(6,*) 'MOVIEFREQ: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_n_tmmc_max(flag,ivalue)
c     ******************************************************************
c     * maximum number of molecules in a TMMC simulation               *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modiied 10-13-2008 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer n_tmmc_max
      save n_tmmc_max

      if ( flag .eq. GLB_GET ) then
         ivalue = n_tmmc_max
      elseif ( flag .eq. GLB_SET ) then
         n_tmmc_max = ivalue
      else
         write(6,*) 'N_TMMC_MAX: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_n_tmmc_min(flag,ivalue)
c     ******************************************************************
c     * min number of molecules in a TMMC simulation                   *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modiied 10-13-2008 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer n_tmmc_min
      save n_tmmc_min

      if ( flag .eq. GLB_GET ) then
         ivalue = n_tmmc_min
      elseif ( flag .eq. GLB_SET ) then
         n_tmmc_min = ivalue
      else
         write(6,*) 'N_TMMC_MIN: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_naabead(flag,ivalue)
c     ******************************************************************
c     * angle-angle terms in a CBMC growth step                        *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modiied 10-05-2008 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer naabead
      save naabead

      if ( flag .eq. GLB_GET ) then
         ivalue = naabead
      elseif ( flag .eq. GLB_SET ) then
         naabead = ivalue
      elseif ( flag .eq. GLB_INCR ) then
         naabead = naabead + ivalue
      else
         write(6,*) 'NAABEAD: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_naasame(flag,itype,ivalue)
c     ******************************************************************
c     * the number of quartets for each angle-angle type in the        *
c     * simulation                                                     *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modiied 04-13-2006 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer ittemp
c     --- storage variable
      integer naasame
      dimension naasame(TAAMAX)
      save naasame

      if ( flag .eq. GLB_SET ) then
         naasame(itype) = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = naasame(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            naasame(ittemp) = ivalue
         enddo
      else
         write(6,*) 'NAASAME: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_native_param(flag,count,dvalue)
c     ******************************************************************
c     * native params in scaled potential                              *
c     *                                                                *
c     * originally written 12-20-2007 by M.G. Martin                   *
c     * last modiied 12-20-2007 by M.G. Martin                         *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,count
      double precision dvalue
c     --- storage variable
      double precision native_param
      dimension native_param(3)
      save native_param

      if ( flag .eq. GLB_GET ) then
         dvalue = native_param(count)
      elseif ( flag .eq. GLB_SET ) then
         native_param(count) = dvalue
      else
         write(6,*) 'NATIVE_PARAM: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nb_to_hbtype(flag,itype,ivalue)
c     ******************************************************************
c     * the conversion of nonbond to hydrogen bond type for certain    *
c     * potentials                                                     *
c     *                                                                *
c     * originally written 09-25-2008 by M.G. Martin                   *
c     * last modified 09-25-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer ittemp
c     --- information variables
      integer nb_to_hbtype
      dimension nb_to_hbtype(CROSSTYPEMAX)
      save nb_to_hbtype


      if ( flag .eq. GLB_SET ) then
         nb_to_hbtype(itype) = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = nb_to_hbtype(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            nb_to_hbtype(ittemp) = ivalue
         enddo
      else
         write(6,*) 'NB_TO_HBTYPE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nbcoeff(flag,icoeff,itype,dvalue)
c     ******************************************************************
c     * the nonbonded coefficients for the van der waals               *
c     * (non-coulombic) portion of the force field                     *
c     *                                                                *
c     * originally written 04-11-2006 by M.G. Martin                   *
c     * last modified 04-11-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,icoeff
      double precision dvalue
c     --- local variables
      integer ittemp,ictemp
c     --- information variables
      double precision nbcoeff
      dimension nbcoeff(MINNBCOEFF:MAXNBCOEFF,CROSSTYPEMAX)
      save nbcoeff

      if ( flag .eq. GLB_SET ) then
         nbcoeff(icoeff,itype) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = nbcoeff(icoeff,itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ictemp = MINNBCOEFF,icoeff
            do ittemp = 1,itype
               nbcoeff(ictemp,ittemp) = dvalue
            enddo
         enddo
      else
         write(6,*) 'NBCOEFF: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nbff(flag,itype,cvalue)
c     ******************************************************************
c     * the name of the forcefield associated with each nonbond type   *
c     *                                                                *
c     * originally written 04-11-2006 by M.G. Martin                   *
c     * last modified 04-11-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype
      character*(*) cvalue
c     --- local variables
      integer ittemp
c     --- information variables
      character*(FFNAMELEN) nbff
      dimension nbff(NNTYPE)
      save nbff

      if ( flag .eq. GLB_SET ) then
         nbff(itype) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = nbff(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            nbff(ittemp) = cvalue
         enddo
      else
         write(6,*) 'NBFF: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nbname(flag,itype,cvalue)
c     ******************************************************************
c     * the name for determining the nonbonded match for this atom type*
c     *                                                                *
c     * originally written 04-12-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype
      character*(*) cvalue
c     --- local variables
      integer ittemp
c     --- information variables
      character*(FFNAMELEN) nbname
      dimension nbname(NNTYPE)
      save nbname

      if ( flag .eq. GLB_SET ) then
         nbname(itype) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = nbname(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            nbname(ittemp) = cvalue
         enddo
      else
         write(6,*) 'NBNAME: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nboxpair(flag,ivalue)
c     ******************************************************************
c     * the number of pairs of simulation boxes in the simulation      *
c     *                                                                *
c     * originally written 10-30-2007 by M.G. Martin                   *
c     * last modified 10-30-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer nboxpair
      save nboxpair

      if ( flag .eq. GLB_GET ) then
         ivalue = nboxpair
      elseif ( flag .eq. GLB_SET ) then
         nboxpair = ivalue
      else
         write(6,*) 'NBOXPAIR: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nbsame(flag,itype,ivalue)
c     ******************************************************************
c     * the number of triplets of bending names associated with this   *
c     * type of bending angle                                          *
c     *                                                                *
c     * originally written 04-05-2006 by M.G. Martin                   *
c     * last modified 04-05-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer ittemp
c     --- information variables
      integer nbsame
      dimension nbsame(TBENMAX)
      save nbsame

      if ( flag .eq. GLB_SET ) then
         nbsame(itype) = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = nbsame(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            nbsame(ittemp) = ivalue
         enddo
      else
         write(6,*) 'NBSAME: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nch_bend_a(flag,imolty,ivalue)
c     ******************************************************************
c     * number of trials for the generation of the bending angle A     *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,ivalue
c     --- local variables
      integer imtemp
c     --- information variables
      integer nch_bend_a
      dimension nch_bend_a(NTMAX)
      save nch_bend_a

      if ( flag .eq. GLB_GET ) then
         ivalue = nch_bend_a(imolty)
      elseif ( flag .eq. GLB_SET ) then
         nch_bend_a(imolty) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            nch_bend_a(imtemp) = ivalue
         enddo
      else
         write(6,*) 'NCH_BEND_A: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nch_bend_b(flag,imolty,ivalue)
c     ******************************************************************
c     * number of trials for the generation of the bending angle B     *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,ivalue
c     --- local variables
      integer imtemp
c     --- information variables
      integer nch_bend_b
      dimension nch_bend_b(NTMAX)
      save nch_bend_b

      if ( flag .eq. GLB_GET ) then
         ivalue = nch_bend_b(imolty)
      elseif ( flag .eq. GLB_SET ) then
         nch_bend_b(imolty) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            nch_bend_b(imtemp) = ivalue
         enddo
      else
         write(6,*) 'NCH_BEND_B: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nch_nb(flag,imolty,ivalue)
c     ******************************************************************
c     * number of trials for the growth of every atom except for the   *
c     * first atom in CBMC moves                                       *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,ivalue
c     --- local variables
      integer imtemp
c     --- information variables
      integer nch_nb
      dimension nch_nb(NTMAX)
      save nch_nb

      if ( flag .eq. GLB_GET ) then
         ivalue = nch_nb(imolty)
      elseif ( flag .eq. GLB_SET ) then
         nch_nb(imolty) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            nch_nb(imtemp) = ivalue
         enddo
      else
         write(6,*) 'NCH_NB: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nch_nb_one(flag,imolty,ivalue)
c     ******************************************************************
c     * number of trials for the insertion of the first atom in the    *
c     * CBMC moves                                                     *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,ivalue
c     --- local variables
      integer imtemp
c     --- information variables
      integer nch_nb_one
      dimension nch_nb_one(NTMAX)
      save nch_nb_one

      if ( flag .eq. GLB_GET ) then
         ivalue = nch_nb_one(imolty)
      elseif ( flag .eq. GLB_SET ) then
         nch_nb_one(imolty) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            nch_nb_one(imtemp) = ivalue
         enddo
      else
         write(6,*) 'NCH_NB_ONE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nchain(flag,ivalue)
c     ******************************************************************
c     * number of chains in the system                                 *
c     *                                                                *
c     * originally written 10-30-2007 by M.G. Martin                   *
c     * last modified 10-30-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer nchain
      save nchain

      if ( flag .eq. GLB_GET ) then
         ivalue = nchain
      elseif ( flag .eq. GLB_SET ) then
         nchain = ivalue
      else
         write(6,*) 'NCHAIN: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nchbox(flag,ibox,ivalue)
c     ******************************************************************
c     * total number of chains in each box                             *
c     *                                                                *
c     * originally written 11-01-2007 by M.G. Martin                   *
c     * last modified 11-01-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ibox,ivalue
c     --- local variables
      integer ibtemp
c     --- information variables
      integer nchbox
      dimension nchbox(0:MAXBOX)
      save nchbox

      if ( flag .eq. GLB_GET ) then
         ivalue = nchbox(ibox)
      elseif ( flag .eq. GLB_SET ) then
         nchbox(ibox) = ivalue
      elseif ( flag .eq. GLB_INCR ) then
         nchbox(ibox) = nchbox(ibox) + ivalue
      elseif ( flag .eq. GLB_DECR ) then
         nchbox(ibox) = nchbox(ibox) - ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 0,ibox
            nchbox(ibtemp) = ivalue
         enddo
      else
         write(6,*) 'NCHBOX: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nch_pre_nb(flag,imolty,ivalue)
c     ******************************************************************
c     * number of sections for the step that comes after torsions but  *
c     * before nonbonded.  The "pre-nonbond" step of CBMC moves        *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,ivalue
c     --- local variables
      integer imtemp
c     --- information variables
      integer nch_pre_nb
      dimension nch_pre_nb(NTMAX)
      save nch_pre_nb

      if ( flag .eq. GLB_GET ) then
         ivalue = nch_pre_nb(imolty)
      elseif ( flag .eq. GLB_SET ) then
         nch_pre_nb(imolty) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            nch_pre_nb(imtemp) = ivalue
         enddo
      else
         write(6,*) 'NCH_PRE_NB: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nch_tor(flag,itype,ivalue)
c     ******************************************************************
c     * number of torsion selection trials for each type of molecule   *
c     *                                                                *
c     * originally written 11-07-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer ittemp
c     --- information variables
      integer nch_tor
      dimension nch_tor(NTMAX)
      save nch_tor

      if ( flag .eq. GLB_GET ) then
         ivalue = nch_tor(itype)
      elseif ( flag .eq. GLB_SET ) then
         nch_tor(itype) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            nch_tor(ittemp) = ivalue
         enddo
      else
         write(6,*) 'NCH_TOR: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nch_tor_connect(flag,imolty,ivalue)
c     ******************************************************************
c     * number of torsion selection trials when connecting back with a *
c     * fixed-endpoint move for each type of molecule                  *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,ivalue
c     --- local variables
      integer imtemp
c     --- information variables
      integer nch_tor_connect
      dimension nch_tor_connect(NTMAX)
      save nch_tor_connect

      if ( flag .eq. GLB_GET ) then
         ivalue = nch_tor_connect(imolty)
      elseif ( flag .eq. GLB_SET ) then
         nch_tor_connect(imolty) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            nch_tor_connect(imtemp) = ivalue
         enddo
      else
         write(6,*) 'NCH_TOR_CONNECT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nch_vib(flag,imolty,ivalue)
c     ******************************************************************
c     * number of vibrational selection trials for each type of        *
c     * molecule                                                       *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,ivalue
c     --- local variables
      integer imtemp
c     --- information variables
      integer nch_vib
      dimension nch_vib(NTMAX)
      save nch_vib

      if ( flag .eq. GLB_GET ) then
         ivalue = nch_vib(imolty)
      elseif ( flag .eq. GLB_SET ) then
         nch_vib(imolty) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            nch_vib(imtemp) = ivalue
         enddo
      else
         write(6,*) 'NCH_VIB: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_ncmt(flag,ibox,itype,ivalue)
c     ******************************************************************
c     * number of each type of molecule in each simulation box         *
c     *                                                                *
c     * originally written 11-01-2007 by M.G. Martin                   *
c     * last modified 11-01-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ibox,itype,ivalue
c     --- local variables
      integer ibtemp,ittemp
c     --- information variables
      integer ncmt
      dimension ncmt(0:MAXBOX,NTMAX)
      save ncmt

      if ( flag .eq. GLB_GET ) then
         ivalue = ncmt(ibox,itype)
      elseif ( flag .eq. GLB_SET ) then
         ncmt(ibox,itype) = ivalue
      elseif ( flag .eq. GLB_INCR ) then
         ncmt(ibox,itype) = ncmt(ibox,itype) + ivalue
      elseif ( flag .eq. GLB_DECR ) then
         ncmt(ibox,itype) = ncmt(ibox,itype) - ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 0,ibox
            do ittemp = 1,itype
               ncmt(ibtemp,ittemp) = ivalue
            enddo
         enddo
      else
         write(6,*) 'NCMT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nener(flag,ivalue)
c     ******************************************************************
c     * the number of energy terms tracked in the simulation           *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer nener
      save nener

      if ( flag .eq. GLB_GET ) then
         ivalue = nener
      elseif ( flag .eq. GLB_SET ) then
         nener = ivalue
      else
         write(6,*) 'NENER: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nexclude(flag,imolty,iunit,ivalue)
c     ******************************************************************
c     * contains the arrays that are used for the inclusion table      *
c     * these should only need to be used by the subroutine and        *
c     * function in the file linclude.F                                *
c     *                                                                *
c     * originally written 12-03-2007 by M.G. Martin                   *
c     * last modified 12-03-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,iunit,ivalue
c     --- information variables
      integer nexclude
      dimension nexclude(NTMAX,NUMAX)
      save nexclude

      if ( flag .eq. GLB_GET ) then
         ivalue = nexclude(imolty,iunit)
      elseif ( flag .eq. GLB_SET ) then
         nexclude(imolty,iunit) = ivalue
      else
         write(6,*) 'NEXCLUDE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nfield(flag,ivalue)
c     ******************************************************************
c     * the total number of fields in the simulation                   *
c     *                                                                *
c     * originally written 11-08-2007 by M.G. Martin                   *
c     * last modified 11-08-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer nfield
      save nfield

      if ( flag .eq. GLB_GET ) then
         ivalue = nfield
      elseif ( flag .eq. GLB_SET ) then
         nfield = ivalue
      elseif ( flag .eq. GLB_INCR ) then
         nfield = nfield + ivalue
      else
         write(6,*) 'NFIELD: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nimpbead(flag,ivalue)
c     ******************************************************************
c     * number of improper interaction during a CBMC regrowth          *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modified 10-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- storage variable
      integer nimpbead
      save nimpbead

      if ( flag .eq. GLB_GET ) then
         ivalue = nimpbead
      elseif ( flag .eq. GLB_SET ) then
         nimpbead = ivalue
      elseif ( flag .eq. GLB_INCR ) then
         nimpbead = nimpbead + 1
      else
         write(6,*) 'NIMPBEAD unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_nimpsame(flag,itype,ivalue)
c     ******************************************************************
c     * the number of quartets of improper names associated with this  *
c     * type of improper torsion                                       *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer ittemp
c     --- information variables
      integer nimpsame
      dimension nimpsame(TIMPMAX)
      save nimpsame

      if ( flag .eq. GLB_SET ) then
         nimpsame(itype) = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = nimpsame(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            nimpsame(ittemp) = ivalue
         enddo
      else
         write(6,*) 'NIMPSAME: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nmaxcbmc(flag,imolty,ivalue)
c     ******************************************************************
c     * maximum number of atoms to regrow in a normal CBMC regrowth    *
c     * move (not an insertion where we obviously regrow all atoms)    *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,ivalue
c     --- local variables
      integer imtemp
c     --- information variables
      integer nmaxcbmc
      dimension nmaxcbmc(NTMAX)
      save nmaxcbmc

      if ( flag .eq. GLB_GET ) then
         ivalue = nmaxcbmc(imolty)
      elseif ( flag .eq. GLB_SET ) then
         nmaxcbmc(imolty) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            nmaxcbmc(imtemp) = ivalue
         enddo
      else
         write(6,*) 'NMAXCBMC: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nmolectyp(flag,itype,ivalue)
c     ******************************************************************
c     * the number of each molecule types in the entire system         *
c     *                                                                *
c     * originally written 11-01-2007 by M.G. Martin                   *
c     * last modified 11-01-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer ittemp
c     --- information variables
      integer nmolectyp
      dimension nmolectyp(NTMAX)
      save nmolectyp

      if ( flag .eq. GLB_GET ) then
         ivalue = nmolectyp(itype)
      elseif ( flag .eq. GLB_SET ) then
         nmolectyp(itype) = ivalue
      elseif ( flag .eq. GLB_INCR ) then
         nmolectyp(itype) = nmolectyp(itype) + ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            nmolectyp(ittemp) = ivalue
         enddo
      else
         write(6,*) 'NMOLECTYP: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nmolty(flag,ivalue)
c     ******************************************************************
c     * the number of molecule types                                   *
c     *                                                                *
c     * originally written 10-30-2007 by M.G. Martin                   *
c     * last modified 10-30-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer nmolty
      save nmolty

      if ( flag .eq. GLB_GET ) then
         ivalue = nmolty
      elseif ( flag .eq. GLB_SET ) then
         nmolty = ivalue
      else
         write(6,*) 'NMOLTY: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nstep(flag,ivalue)
c     ******************************************************************
c     * the number of Monte Carlo moves steps to perform               *
c     *                                                                *
c     * originally written 12-09-2007 by M.G. Martin                   *
c     * last modified 12-09-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer nstep
      save nstep

      if ( flag .eq. GLB_GET ) then
         ivalue = nstep
      elseif ( flag .eq. GLB_SET ) then
         nstep = ivalue
      else
         write(6,*) 'NSTEP: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_ntsame(flag,itype,ivalue)
c     ******************************************************************
c     * the number of quartests torsion names associated with this     *
c     * type of torsion angle                                          *
c     *                                                                *
c     * originally written 04-10-2006 by M.G. Martin                   *
c     * last modified 04-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer ittemp
c     --- information variables
      integer ntsame
      dimension ntsame(TTORMAX)
      save ntsame

      if ( flag .eq. GLB_SET ) then
         ntsame(itype) = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = ntsame(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            ntsame(ittemp) = ivalue
         enddo
      else
         write(6,*) 'NTSAME: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_ntorloop(flag,itype,ivalue)
c     ******************************************************************
c     * number of similar loops in certain torsion styles              *
c     *                                                                *
c     * originally written 04-05-2006 by M.G. Martin                   *
c     * last modified 04-05-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer ittemp
c     --- information variables
      integer ntorloop
      dimension ntorloop(TTORMAX)
      save ntorloop

      if ( flag .eq. GLB_SET ) then
         ntorloop(itype) = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = ntorloop(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            ntorloop(itype) = ivalue
         enddo
      else
         write(6,*) 'NTORLOOP: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_ntype(flag,itype,iunit,ivalue)
c     ******************************************************************
c     * the nonbonded atom type for each unit in each molecule type    *
c     *                                                                *
c     * originally written 11-01-2007 by M.G. Martin                   *
c     * last modified 11-01-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,iunit,ivalue
c     --- local variables
      integer ittemp,iutemp
c     --- information variables
      integer ntype
      dimension ntype(NTMAX,NUMAX)
      save ntype

      if ( flag .eq. GLB_GET ) then
         ivalue = ntype(itype,iunit)
      elseif ( flag .eq. GLB_SET ) then
         ntype(itype,iunit) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            do iutemp = 1,iunit
               ntype(ittemp,iutemp) = ivalue
            enddo
         enddo
      else
         write(6,*) 'NTYPE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_numboxes(flag,ivalue)
c     ******************************************************************
c     * number of simulation boxes                                     *
c     *                                                                *
c     * originally written 10-06-2008 by M.G. Martin                   *
c     * last modified 10-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer numboxes
      save numboxes

      if ( flag .eq. GLB_GET ) then
         ivalue = numboxes
      elseif ( flag .eq. GLB_SET ) then
         numboxes = ivalue
      else
         write(6,*) 'NUMBOXES: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nunit(flag,itype,ivalue)
c     ******************************************************************
c     * number of units (atoms or pseudo-atoms) in each molecule type  *
c     *                                                                *
c     * originally written 10-31-2007 by M.G. Martin                   *
c     * last modified 10-31-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer ittemp
c     --- information variables
      integer nunit
      dimension nunit(NTMAX)
      save nunit

      if ( itype .gt. NTMAX ) then
         write(6,*) 'flag,itype,ivalue',flag,itype,ivalue
      endif

      if ( flag .eq. GLB_GET ) then
         ivalue = nunit(itype)
      elseif ( flag .eq. GLB_SET ) then
         nunit(itype) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            nunit(itype) = ivalue
         enddo
      else
         write(6,*) 'NUNIT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nvmove(flag,ivalue)
c     ******************************************************************
c     * number of possible boxes or pairs of boxes (depending on the   *
c     * ensemble) for the volume moves                                 *
c     *                                                                *
c     * originally written 10-06-2008 by M.G. Martin                   *
c     * last modified 10-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer nvmove
      save nvmove

      if ( flag .eq. GLB_GET ) then
         ivalue = nvmove
      elseif ( flag .eq. GLB_SET ) then
         nvmove = ivalue
      else
         write(6,*) 'NVMOVE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nvsame(flag,itype,ivalue)
c     ******************************************************************
c     * the number of pairs of vibration names associated with this    *
c     * type of vibration                                              *
c     *                                                                *
c     * originally written 03-30-2006 by M.G. Martin                   *
c     * last modified 03-30-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer ittemp
c     --- information variables
      integer nvsame
      dimension nvsame(TVIBMAX)
      save nvsame

      if ( flag .eq. GLB_SET ) then
         nvsame(itype) = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = nvsame(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            nvsame(ittemp) = ivalue
         enddo
      else
         write(6,*) 'NVSAME: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nxangle(flag,ibond,ivalue)
c     ******************************************************************
c     * number of extra angle terms used in CBMC moves                 *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modified 10-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ibond,ivalue
c     --- information variables
      integer nxangle
      dimension nxangle(NNBOND)
      save nxangle

      if ( flag .eq. GLB_GET ) then
         ivalue = nxangle(ibond)
      elseif ( flag .eq. GLB_SET ) then
         nxangle(ibond) = ivalue
      elseif ( flag .eq. GLB_INCR ) then
         nxangle(ibond) = nxangle(ibond) + ivalue
      else
         write(6,*) 'NXANGLE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_nxbias(flag,ibond,ivalue)
c     ******************************************************************
c     * number of extra bias terms used in CBMC moves                  *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modified 10-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ibond,ivalue
c     --- information variables
      integer nxbias
      dimension nxbias(NNBOND)
      save nxbias

      if ( flag .eq. GLB_GET ) then
         ivalue = nxbias(ibond)
      elseif ( flag .eq. GLB_SET ) then
         nxbias(ibond) = ivalue
      elseif ( flag .eq. GLB_INCR ) then
         nxbias(ibond) = nxbias(ibond) + ivalue
      else
         write(6,*) 'NXBIAS: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_ofcoeff(flag,itype,icoeff,dvalue)
c     ******************************************************************
c     * the coefficients for each type of special 1-5 interaction      *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,icoeff
      double precision dvalue
c     --- local variables
      integer ittemp,ictemp
c     --- information variables
      double precision ofcoeff
      dimension ofcoeff(MAXOF,MINOFCOEFF:MAXOFCOEFF)
      save ofcoeff

      if ( flag .eq. GLB_SET ) then
         ofcoeff(itype,icoeff) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = ofcoeff(itype,icoeff)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            do ictemp = MINOFCOEFF,icoeff
               ofcoeff(ittemp,ictemp) = dvalue
            enddo
         enddo
      else
         write(6,*) 'OFCOEFF: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_offf(flag,itype,cvalue)
c     ******************************************************************
c     * the force field for each of the special one-five interation    *
c     * types                                                          *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype
      character*(*) cvalue
c     --- local variables
      integer ittemp
c     --- information variables
      character*(FFNAMELEN) offf
      dimension offf(TOFMAX)
      save offf

      if ( flag .eq. GLB_SET ) then
         offf(itype) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = offf(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            offf(ittemp) = cvalue
         enddo
      else
         write(6,*) 'OFFF: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_ofnames(flag,itype,index,cvalue)
c     ******************************************************************
c     * the atom names for each of the special one-five interation     *
c     * types                                                          *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,index
      character*(*) cvalue
c     --- local variables
      integer ittemp,iitemp
c     --- information variables
      character*(FFNAMELEN) ofnames
      dimension ofnames(TOFMAX,MAXOFNAME)
      save ofnames

      if ( flag .eq. GLB_SET ) then
         ofnames(itype,index) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = ofnames(itype,index)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            do iitemp = 1,index
               ofnames(ittemp,iitemp) = cvalue
            enddo
         enddo
      else
         write(6,*) 'OFNAMES: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_ofstyle(flag,itype,ivalue)
c     ******************************************************************
c     * the style of the special one-five interation for each of type  *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer ittemp
c     --- information variables
      integer ofstyle
      dimension ofstyle(TOFMAX)
      save ofstyle

      if ( flag .eq. GLB_SET ) then
         ofstyle(itype) = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = ofstyle(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            ofstyle(ittemp) = ivalue
         enddo
      else
         write(6,*) 'OFSTYLE: invalid global flag',flag
         stop
      endif

      return
      end

      function twh_oldgetatomnum(ichain,iunit)
c     ******************************************************************
c     * returns the atom position in the coordinate arrays             *
c     *                                                                *
c     * originally written 03-25-2002 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_OLDGETATOMNUM
#include "functions.h"
c     --- variables passed to the function
      integer ichain,iunit,globalpos

      call twh_globalpos(GLB_GET,ichain,globalpos)
      twh_oldgetatomnum = (globalpos-1)/3 + iunit
      return
      end

      subroutine twh_orderben(flag,imolty,iunit,ibend,ivalue)
c     ******************************************************************
c     * the order (forward or reverse) for each angle in each type of  *
c     * molecule                                                       *
c     *                                                                *
c     * originally written 12-11-2006 by M.G. Martin                   *
c     * last modified 12-11-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty,iunit,ibend,ivalue
c     --- local variables
      integer imtemp,iutemp,ibtemp
c     --- storage variable
      integer orderben
      dimension orderben(NTMAX,NUMAX,MAXBEND)
      save orderben

      if ( flag .eq. GLB_GET ) then
         ivalue = orderben(imolty,iunit,ibend)
      elseif ( flag .eq. GLB_SET ) then
         orderben(imolty,iunit,ibend) = ivalue
      elseif ( flag .eq. GLB_INIT ) then 
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibend
                  orderben(imtemp,iutemp,ibtemp) = ivalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'ORDERBEN unknown flag:',flag
         stop
      endif

      return
      end

      subroutine twh_p_acc(flag,dvalue)
c     ******************************************************************
c     * TMMC variable p_acc                                            *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision p_acc
      save p_acc

      if ( flag .eq. GLB_GET ) then
         dvalue = p_acc
      elseif ( flag .eq. GLB_SET ) then
         p_acc = dvalue
      else
         write(6,*) 'P_ACC: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pairbox(flag,pbox,index,ivalue)
c     ******************************************************************
c     * the list of box pairs in the simulation                        *
c     *                                                                *
c     * originally written 10-31-2007 by M.G. Martin                   *
c     * last modified 10-31-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,pbox,index,ivalue
c     --- information variables
      integer pairbox
      dimension pairbox(MAXPBOX,2)
      save pairbox

      if ( flag .eq. GLB_GET ) then
         ivalue = pairbox(pbox,index)
      elseif ( flag .eq. GLB_SET ) then
         pairbox(pbox,index) = ivalue
      else
         write(6,*) 'PAIRBOX: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_parallelstyle(flag,ivalue)
c     ******************************************************************
c     * the parallel style for use in this simulation                  *
c     *                                                                *
c     * originally written 09-24-2008 by M.G. Martin                   *
c     * last modified 09-24-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer parallelstyle
      save parallelstyle

      if ( flag .eq. GLB_GET ) then
         ivalue = parallelstyle
      elseif ( flag .eq. GLB_SET ) then
         parallelstyle = ivalue
      else
         write(6,*) 'PARALELLSTYLE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pcomp(flag,istress,cvalue)
c     ******************************************************************
c     * names for the components of the stress tensor (pressure)       *
c     *                                                                *
c     * originally written 10-11-2008 by M.G. Martin                   *
c     * last modified 10-11-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,istress
      character*(*) cvalue
c     --- information variables
      character*13 pcomp
      dimension pcomp(0:7) 
      save pcomp

      if ( flag .eq. GLB_GET ) then
         cvalue = pcomp(istress)
      elseif ( flag .eq. GLB_SET ) then
         pcomp(istress) = cvalue
      else
         write(6,*) 'PCOMP: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pdb_output_freq(flag,ivalue)
c     ******************************************************************
c     * frequency for pdb output during a simulation run               *
c     *                                                                *
c     * originally written 12-22-2007 by M.G. Martin                   *
c     * last modified 12-22-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer pdb_output_freq
      save pdb_output_freq

      if ( flag .eq. GLB_GET ) then
         ivalue = pdb_output_freq
      elseif ( flag .eq. GLB_SET ) then
         pdb_output_freq = ivalue
      else
         write(6,*) 'PDB_OUTPUT_FREQ: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pdbname(flag,imolty,iunit,cvalue)
c     ******************************************************************
c     * pdb name for each atom in each molecule type use for outputting*
c     * pdb files                                                      *
c     *                                                                *
c     * originally written 10-02-2008 by M.G. Martin                   *
c     * last modified 10-02-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,iunit
      character*(*) cvalue
c     --- information variables
      character*4 pdbname
      dimension pdbname(NTMAX,NUMAX)
      save pdbname

      if ( flag .eq. GLB_GET ) then
         cvalue = pdbname(imolty,iunit)
      elseif ( flag .eq. GLB_SET ) then
         pdbname(imolty,iunit) = cvalue
      else
         write(6,*) 'PDBNAME: invalid global flag',flag
         stop
      endif

      return
      end


      subroutine twh_planewidth(flag,dvalue)
c     ******************************************************************
c     * the width of the plane for the planeshift move                 *
c     *                                                                *
c     * originally written 08-11-2008 by M.G. Martin                   *
c     * last modified 08-11-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision planewidth
      save planewidth

      if ( flag .eq. GLB_SET ) then
         planewidth = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = planewidth
      else
         write(6,*) 'PLANEWIDTH: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pm1boxcbswap(flag,dvalue)
c     ******************************************************************
c     * the probability sum for the one-box configurational-bias swap  *
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modified 12-24-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pm1boxcbswap
      save pm1boxcbswap

      if ( flag .eq. GLB_SET ) then
         pm1boxcbswap = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = pm1boxcbswap
      else
         write(6,*) 'PM1BOXCBSWAP: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pm1cbswmt(flag,imolty,dvalue)
c     ******************************************************************
c     * the probability sum for selecting the molecule type of the     *
c     * one-box configurational-bias swap move                         *
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modified 12-24-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer imtemp
c     --- information variables
c     --- variables for the cb intra-box swap move
      double precision pm1cbswmt
      dimension pm1cbswmt(NTMAX)
      save pm1cbswmt

      if ( flag .eq. GLB_GET ) then
         dvalue = pm1cbswmt(imolty)
      elseif ( flag .eq. GLB_SET ) then
         pm1cbswmt(imolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            pm1cbswmt(imtemp) = dvalue
         enddo
      else
         write(6,*) 'PM1CBSWMT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pm1boxcomswitch(flag,dvalue)
c     ******************************************************************
c     * the probability sum for the one-box center-of-mass molecule    *
c     * switch move                                                    *
c     *                                                                *
c     * originally written 07-28-2006 by M.G. Martin                   *
c     * last modified 07-28-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pm1boxcomswitch
      save pm1boxcomswitch

      if ( flag .eq. GLB_SET ) then
         pm1boxcomswitch = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = pm1boxcomswitch
      else
         write(6,*) 'PM1BOXCOMSWITCH: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pm1comswbox(flag,ibox,dvalue)
c     ******************************************************************
c     * the probability sum for selecting the box of the one-box       *
c     * center-of-mass molecule switch move                            *
c     *                                                                *
c     * originally written 07-28-2006 by M.G. Martin                   *
c     * last modified 07-28-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ibox
      double precision dvalue
c     --- local variables
      integer ibtemp
c     --- information variables
      double precision pm1comswbox
      dimension pm1comswbox(MAXBOX)
      save pm1comswbox

      if ( flag .eq. GLB_SET ) then
         pm1comswbox(ibox) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = pm1comswbox(ibox)
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 1,ibox
            pm1comswbox(ibtemp) = 0.0d0
         enddo
      else
         write(6,*) 'PM1COMSWBOX: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pm1comswpair(flag,ipair,dvalue)
c     ******************************************************************
c     * the probability sum for selecting the pair for the one-box     *
c     * center-of-mass molecule switch move                            *
c     *                                                                *
c     * originally written 07-28-2006 by M.G. Martin                   *
c     * last modified 06-14-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ipair
      double precision dvalue
c     --- local variables
      integer iptemp
c     --- information variables
      double precision pm1comswpair
      dimension pm1comswpair(NTPAIRMAX)
      save pm1comswpair

      if ( flag .eq. GLB_GET ) then
         dvalue = pm1comswpair(ipair)
      elseif ( flag .eq. GLB_SET ) then
         pm1comswpair(ipair) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do iptemp = 1,ipair
            pm1comswpair(iptemp) = 0.0d0
         enddo
      else
         write(6,*) 'PM1COMSWPAIR: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pm2boxcbswap(flag,dvalue)
c     ******************************************************************
c     * the probability sum for selecting the 2-box                    *
c     * configurational-bias molecule transfer move                    *
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modified 12-24-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pm2boxcbswap
      save pm2boxcbswap

      if ( flag .eq. GLB_GET ) then
         dvalue = pm2boxcbswap
      elseif ( flag .eq. GLB_SET ) then
         pm2boxcbswap = dvalue
      else
         write(6,*) 'PM2BOXCBSWAP: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pm2cbswmt(flag,imolty,dvalue)
c     ******************************************************************
c     * the probability sum for selecting each molecule type for the   *
c     * 2-box configurational-bias molecule transfer move              *
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modified 12-24-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer imtemp
c     --- information variables
      double precision pm2cbswmt
      dimension pm2cbswmt(NTMAX)
      save pm2cbswmt

      if ( flag .eq. GLB_GET ) then
         dvalue = pm2cbswmt(imolty)
      elseif ( flag .eq. GLB_SET ) then
         pm2cbswmt(imolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            pm2cbswmt(imtemp) = dvalue
         enddo
      else
         write(6,*) 'PM2CBSWMT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pm2cbswpr(flag,iboxpair,dvalue)
c     ******************************************************************
c     * the probability sum for selecting each pair of boxes for the   *
c     * 2-box configurational-bias molecule transfer move              *
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modified 12-24-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,iboxpair
      double precision dvalue
c     --- local variables
      integer ibtemp
c     --- information variables
      double precision pm2cbswpr
      dimension pm2cbswpr(MAXPBOX)
      save pm2cbswpr

      if ( flag .eq. GLB_GET ) then
         dvalue = pm2cbswpr(iboxpair)
      elseif ( flag .eq. GLB_SET ) then
         pm2cbswpr(iboxpair) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 1,iboxpair
            pm2cbswpr(ibtemp) = dvalue
         enddo
      else
         write(6,*) 'PM2CBSWPR: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pm2boxrbswap(flag,dvalue)
c     ******************************************************************
c     * the probability sum for selecting the 2-box rotational-bias    *
c     * molecule transfer move                                         *
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modified 12-24-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pm2boxrbswap
      save pm2boxrbswap

      if ( flag .eq. GLB_GET ) then
         dvalue = pm2boxrbswap
      elseif ( flag .eq. GLB_SET ) then
         pm2boxrbswap = dvalue
      else
         write(6,*) 'PM2BOXRBSWAP: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pm2rbswmt(flag,imolty,dvalue)
c     ******************************************************************
c     * the probability sum for selecting each molecule type for the   *
c     * 2-box rotational-bias molecule transfer move                   *
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modified 12-24-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer imtemp
c     --- information variables
      double precision pm2rbswmt
      dimension pm2rbswmt(NTMAX)
      save pm2rbswmt

      if ( flag .eq. GLB_GET ) then
         dvalue = pm2rbswmt(imolty)
      elseif ( flag .eq. GLB_SET ) then
         pm2rbswmt(imolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            pm2rbswmt(imtemp) = dvalue
         enddo
      else
         write(6,*) 'PM2RBSWMT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pm2rbswpr(flag,iboxpair,dvalue)
c     ******************************************************************
c     * the probability sum for selecting each pair of boxes for the   *
c     * 2-box rotational-bias molecule transfer move                   *
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modified 12-24-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,iboxpair
      double precision dvalue
c     --- local variables
      integer ibtemp
c     --- information variables
      double precision pm2rbswpr
      dimension pm2rbswpr(MAXPBOX)
      save pm2rbswpr

      if ( flag .eq. GLB_GET ) then
         dvalue = pm2rbswpr(iboxpair)
      elseif ( flag .eq. GLB_SET ) then
         pm2rbswpr(iboxpair) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 1,iboxpair
            pm2rbswpr(ibtemp) = dvalue
         enddo
      else
         write(6,*) 'PM2RBSWpr: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmall(flag,imolty,dvalue)
c     ******************************************************************
c     * the probability of regrowing all of the atoms in a cbregrow    *
c     * move instead of just regrowing in a single direction.  One     *
c     * entry for each molecule type                                   *
c     *                                                                *
c     * originally written 08-10-2008 by M.G. Martin                   *
c     * last modified 08-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer imtemp
c     --- information variables
      double precision pmall
      dimension pmall(NTMAX)
      save pmall

      if ( flag .eq. GLB_GET ) then
         dvalue = pmall(imolty)
      elseif ( flag .eq. GLB_SET ) then
         pmall(imolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            pmall(imtemp) = dvalue
         enddo
      else
         write(6,*) 'PMALL: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmavb1(flag,dvalue)
c     ******************************************************************
c     * the probability sum for the aggregation-volume-bias move type 1*
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modified 12-24-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmavb1
      save pmavb1

      if ( flag .eq. GLB_GET ) then
         dvalue = pmavb1
      elseif ( flag .eq. GLB_SET ) then
         pmavb1 = dvalue
      else
         write(6,*) 'PMAVB1: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmavb1ct(flag,imolty,jmolty,dvalue)
c     ******************************************************************
c     * the probability sum for the aggregation-volume-bias move type 1*
c     * for selecting the second molecule type                         *
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modified 12-24-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,jmolty
      double precision dvalue
c     --- local variables
      integer imtemp,intemp
c     --- information variables
c     --- variables for aggregation-volume-bias 1
      double precision pmavb1ct
      dimension pmavb1ct(NTMAX,NTMAX)
      save pmavb1ct

      if ( flag .eq. GLB_GET ) then
         dvalue = pmavb1ct(imolty,jmolty)
      elseif ( flag .eq. GLB_SET ) then
         pmavb1ct(imolty,jmolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do intemp = 1,jmolty
               pmavb1ct(imtemp,intemp) = dvalue
            enddo
         enddo
      else
         write(6,*) 'PMAVB1CT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmavb1in(flag,dvalue)
c     ******************************************************************
c     * the probability sum for the aggregation-volume-bias move type 1*
c     * selecting the inner shell as the target                        *
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modified 12-24-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmavb1in
      save pmavb1in

      if ( flag .eq. GLB_GET ) then
         dvalue = pmavb1in
      elseif ( flag .eq. GLB_SET ) then
         pmavb1in = dvalue
      else
         write(6,*) 'PMAVB1IN: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmavb1mt(flag,imolty,dvalue)
c     ******************************************************************
c     * the probability sum for the aggregation-volume-bias move type 1*
c     * for selecting the molecule types                               *
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modified 12-24-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer imtemp
c     --- information variables
      double precision pmavb1mt
      dimension pmavb1mt(NTMAX)
      save pmavb1mt

      if ( flag .eq. GLB_GET ) then
         dvalue = pmavb1mt(imolty)
      elseif ( flag .eq. GLB_SET ) then
         pmavb1mt(imolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            pmavb1mt(imtemp) = dvalue
         enddo
      else
         write(6,*) 'PMAVB1MT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmavb2(flag,dvalue)
c     ******************************************************************
c     * the probability sum for the aggregation-volume-bias move type 2*
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modified 12-24-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmavb2
      save pmavb2

      if ( flag .eq. GLB_GET ) then
         dvalue = pmavb2
      elseif ( flag .eq. GLB_SET ) then
         pmavb2 = dvalue
      else
         write(6,*) 'PMAVB2: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmavb2ct(flag,imolty,jmolty,dvalue)
c     ******************************************************************
c     * the probability sum for the aggregation-volume-bias move type 2*
c     * for selecting the second molecule type                         *
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modified 12-24-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,jmolty
      double precision dvalue
c     --- local variables
      integer imtemp,intemp
c     --- information variables
c     --- variables for aggregation-volume-bias 1
      double precision pmavb2ct
      dimension pmavb2ct(NTMAX,NTMAX)
      save pmavb2ct

      if ( flag .eq. GLB_GET ) then
         dvalue = pmavb2ct(imolty,jmolty)
      elseif ( flag .eq. GLB_SET ) then
         pmavb2ct(imolty,jmolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do intemp = 1,jmolty
               pmavb2ct(imtemp,intemp) = dvalue
            enddo
         enddo
      else
         write(6,*) 'PMAVB2CT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmavb2in(flag,dvalue)
c     ******************************************************************
c     * the probability sum for the aggregation-volume-bias move type 2*
c     * selecting the inner shell as the target                        *
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modified 12-24-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmavb2in
      save pmavb2in

      if ( flag .eq. GLB_GET ) then
         dvalue = pmavb2in
      elseif ( flag .eq. GLB_SET ) then
         pmavb2in = dvalue
      else
         write(6,*) 'PMAVB2IN: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmavb2mt(flag,imolty,dvalue)
c     ******************************************************************
c     * the probability sum for the aggregation-volume-bias move type 2*
c     * for selecting the molecule types                               *
c     *                                                                *
c     * originally written 12-24-2007 by M.G. Martin                   *
c     * last modified 12-24-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer imtemp
c     --- information variables
      double precision pmavb2mt
      dimension pmavb2mt(NTMAX)
      save pmavb2mt

      if ( flag .eq. GLB_GET ) then
         dvalue = pmavb2mt(imolty)
      elseif ( flag .eq. GLB_SET ) then
         pmavb2mt(imolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            pmavb2mt(imtemp) = dvalue
         enddo
      else
         write(6,*) 'PMAVB2MT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmavb3(flag,dvalue)
c     ******************************************************************
c     * the probability sum for the aggregation-volume-bias move type 3*
c     *                                                                *
c     * originally written 02-07-2008 by M.G. Martin                   *
c     * last modified 02-07-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmavb3
      save pmavb3

      if ( flag .eq. GLB_GET ) then
         dvalue = pmavb3
      elseif ( flag .eq. GLB_SET ) then
         pmavb3 = dvalue
      else
         write(6,*) 'PMAVB3: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmavb3ct(flag,imolty,jmolty,dvalue)
c     ******************************************************************
c     * the probability sum for the aggregation-volume-bias move type 3*
c     * for selecting the second molecule type                         *
c     *                                                                *
c     * originally written 02-07-2008 by M.G. Martin                   *
c     * last modified 02-07-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,jmolty
      double precision dvalue
c     --- local variables
      integer imtemp,intemp
c     --- information variables
c     --- variables for aggregation-volume-bias 3
      double precision pmavb3ct
      dimension pmavb3ct(NTMAX,NTMAX)
      save pmavb3ct

      if ( flag .eq. GLB_GET ) then
         dvalue = pmavb3ct(imolty,jmolty)
      elseif ( flag .eq. GLB_SET ) then
         pmavb3ct(imolty,jmolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do intemp = 1,jmolty
               pmavb3ct(imtemp,intemp) = dvalue
            enddo
         enddo
      else
         write(6,*) 'PMAVB3CT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmavb3mt(flag,imolty,dvalue)
c     ******************************************************************
c     * the probability sum for the aggregation-volume-bias move type 3*
c     * for selecting the molecule types                               *
c     *                                                                *
c     * originally written 02-07-2008 by M.G. Martin                   *
c     * last modified 02-07-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer imtemp
c     --- information variables
      double precision pmavb3mt
      dimension pmavb3mt(NTMAX)
      save pmavb3mt

      if ( flag .eq. GLB_GET ) then
         dvalue = pmavb3mt(imolty)
      elseif ( flag .eq. GLB_SET ) then
         pmavb3mt(imolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            pmavb3mt(imtemp) = dvalue
         enddo
      else
         write(6,*) 'PMAVB3MT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmback(flag,dvalue)
c     ******************************************************************
c     * the probability sum for the cbmc backbone regrowth move        *
c     *                                                                *
c     * originally written 08-10-2008 by M.G. Martin                   *
c     * last modified 08-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmback
      save pmback

      if ( flag .eq. GLB_SET ) then
         pmback = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = pmback
      else
         write(6,*) 'PMBACK: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmbkmt(flag,imolty,dvalue)
c     ******************************************************************
c     * the probability sum for selecting molecules types for the cbmc *
c     * backbone regrowth move                                         *
c     *                                                                *
c     * originally written 08-10-2008 by M.G. Martin                   *
c     * last modified 08-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer imtemp
c     --- information variables
      double precision pmbkmt
      dimension pmbkmt(NTMAX)
      save pmbkmt

      if ( flag .eq. GLB_SET ) then
         pmbkmt(imolty) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = pmbkmt(imolty)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            pmbkmt(imtemp) = dvalue
         enddo
      else
         write(6,*) 'PMBKMT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmcb(flag,dvalue)
c     ******************************************************************
c     * the probability sum for the cbmc regrowth move                 *
c     *                                                                *
c     * originally written 08-10-2008 by M.G. Martin                   *
c     * last modified 08-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmcb
      save pmcb

      if ( flag .eq. GLB_SET ) then
         pmcb = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = pmcb
      else
         write(6,*) 'PMCB: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmcbmt(flag,imolty,dvalue)
c     ******************************************************************
c     * the probability sum for selection molecule type for the cbmc   *
c     * regrowth move                                                  *
c     *                                                                *
c     * originally written 08-10-2008 by M.G. Martin                   *
c     * last modified 08-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer imtemp
c     --- information variables
      double precision pmcbmt
      dimension pmcbmt(NTMAX)
      save pmcbmt

      if ( flag .eq. GLB_SET ) then
         pmcbmt(imolty) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = pmcbmt(imolty)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            pmcbmt(imtemp) = dvalue
         enddo
      else
         write(6,*) 'PMCBMT: invalid global flag',flag
         stop
      endif

      return
      end



      subroutine twh_pmcbside(flag,dvalue)
c     ******************************************************************
c     * the probability sum for the cbmc side-chain regrowth move      *
c     *                                                                *
c     * originally written 06-12-2006 by M.G. Martin                   *
c     * last modified 06-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmcbside
      save pmcbside

      if ( flag .eq. GLB_SET ) then
         pmcbside = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = pmcbside
      else
         write(6,*) 'PMCBSIDE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmcbsidemt(flag,imolty,dvalue)
c     ******************************************************************
c     * the probability sum for performing the side-chain cbmc move    *
c     * on each molecule type in the simulation                        *
c     *                                                                *
c     * originally written 06-12-2006 by M.G. Martin                   *
c     * last modified 06-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer imtemp
c     --- information variables
      double precision pmcbsidemt
      dimension pmcbsidemt(NTMAX)
      save pmcbsidemt

      if ( flag .eq. GLB_SET ) then
         pmcbsidemt(imolty) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = pmcbsidemt(imolty)
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            pmcbsidemt(imtemp) = dvalue
         enddo
      else
         write(6,*) 'PMCBSIDEMT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmcell(flag,dvalue)
c     ******************************************************************
c     * the probability sum for performing the volume move that changes*
c     * the shape of the unit cell                                     *
c     *                                                                *
c     * originally written 12-23-2007 by M.G. Martin                   *
c     * last modified 12-23-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmcell
      save pmcell

      if ( flag .eq. GLB_GET ) then
         dvalue = pmcell
      elseif ( flag .eq. GLB_SET ) then
         pmcell = dvalue
      else
         write(6,*) 'PMCELL: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmcellpr(flag,iboxpair,dvalue)
c     ******************************************************************
c     * the probability sum for performing the volume cell move on each*
c     * box (or box pair)                                              *
c     *                                                                *
c     * originally written 12-23-2007 by M.G. Martin                   *
c     * last modified 12-23-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,iboxpair
      double precision dvalue
c     --- local variables
      integer ibtemp
c     --- information variables
      double precision pmcellpr
      dimension pmcellpr(MAXPBOX)
      save pmcellpr

      if ( flag .eq. GLB_GET ) then
         dvalue = pmcellpr(iboxpair)
      elseif ( flag .eq. GLB_SET ) then
         pmcellpr(iboxpair) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 1,iboxpair
            pmcellpr(ibtemp) = dvalue
         enddo
      else
         write(6,*) 'PMCELLPR: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmcellpt(flag,iboxpair,dvalue)
c     ******************************************************************
c     * the probability sum for performing the volume cell move        *
c     * on each box of the pair, while performing an isotropic move on *
c     * the other box                                                  *
c     *                                                                *
c     * originally written 12-23-2007 by M.G. Martin                   *
c     * last modified 12-23-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,iboxpair
      double precision dvalue
c     --- local variables
      integer ibtemp
c     --- information variables
      double precision pmcellpt
      dimension pmcellpt(MAXPBOX)
      save pmcellpt

      if ( flag .eq. GLB_GET ) then
         dvalue = pmcellpt(iboxpair)
      elseif ( flag .eq. GLB_SET ) then
         pmcellpt(iboxpair) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 1,iboxpair
            pmcellpt(ibtemp) = dvalue
         enddo
      else
         write(6,*) 'PMCELLPT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmcomposite(flag,dvalue)
c     ******************************************************************
c     * the probability sum for performing the composite translation   *
c     * and rotation move                                              *
c     *                                                                *
c     * originally written 09-24-2008 by M.G. Martin                   *
c     * last modified 09-24-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmcomposite
      save pmcomposite

      if ( flag .eq. GLB_GET ) then
         dvalue = pmcomposite
      elseif ( flag .eq. GLB_SET ) then
         pmcomposite = dvalue
      else
         write(6,*) 'PMCOMPOSITE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmcomt(flag,imolty,dvalue)
c     ******************************************************************
c     * the probability sum for performing the composite translation   *
c     * and rotation move on each molecule type                        *
c     *                                                                *
c     * originally written 09-24-2008 by M.G. Martin                   *
c     * last modified 09-24-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer imtemp
c     --- information variables
      double precision pmcomt
      dimension pmcomt(NTMAX)
      save pmcomt

      if ( flag .eq. GLB_GET ) then
         dvalue = pmcomt(imolty)
      elseif ( flag .eq. GLB_SET ) then
         pmcomt(imolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            pmcomt(imtemp) = dvalue
         enddo
      else
         write(6,*) 'PMCOMT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmconrot(flag,dvalue)
c     ******************************************************************
c     * the probability sum for performing the concerted rotation move *
c     *                                                                *
c     * originally written 08-11-2008 by M.G. Martin                   *
c     * last modified 08-11-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmconrot
      save pmconrot

      if ( flag .eq. GLB_GET ) then
         dvalue = pmconrot
      elseif ( flag .eq. GLB_SET ) then
         pmconrot = dvalue
      else
         write(6,*) 'PMCONROT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmcrback(flag,dvalue)
c     ******************************************************************
c     * the probability sum for performing the concerted rotation      *
c     * polypeptide backbone move                                      *
c     *                                                                *
c     * originally written 08-11-2008 by M.G. Martin                   *
c     * last modified 08-11-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmcrback
      save pmcrback

      if ( flag .eq. GLB_GET ) then
         dvalue = pmcrback
      elseif ( flag .eq. GLB_SET ) then
         pmcrback = dvalue
      else
         write(6,*) 'PMCRBACK: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmcrbmt(flag,imolty,dvalue)
c     ******************************************************************
c     * the probability sum for performing the concerted rotation      *
c     * polypeptide backbone move on each molecule type                *
c     *                                                                *
c     * originally written 08-11-2008 by M.G. Martin                   *
c     * last modified 08-11-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer imtemp
c     --- information variables
      double precision pmcrbmt
      dimension pmcrbmt(NTMAX)
      save pmcrbmt

      if ( flag .eq. GLB_GET ) then
         dvalue = pmcrbmt(imolty)
      elseif ( flag .eq. GLB_SET ) then
         pmcrbmt(imolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            pmcrbmt(imtemp) = dvalue
         enddo
      else
         write(6,*) 'PMCRBMT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmcrmt(flag,imolty,dvalue)
c     ******************************************************************
c     * the probability sum for performing the concerted rotation move *
c     * on each molecule type in the simulation                        *
c     *                                                                *
c     * originally written 08-11-2008 by M.G. Martin                   *
c     * last modified 08-11-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer imtemp
c     --- information variables
      double precision pmcrmt
      dimension pmcrmt(NTMAX)
      save pmcrmt

      if ( flag .eq. GLB_GET ) then
         dvalue = pmcrmt(imolty)
      elseif ( flag .eq. GLB_SET ) then
         pmcrmt(imolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            pmcrmt(imtemp) = dvalue
         enddo
      else
         write(6,*) 'PMCRMT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmpivmt(flag,imolty,dvalue)
c     ******************************************************************
c     * the probability sum for performing the pivot move on each      *
c     * molecule type in the simulation                                *
c     *                                                                *
c     * originally written 08-11-2008 by M.G. Martin                   *
c     * last modified 08-11-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer imtemp
c     --- information variables
      double precision pmpivmt
      dimension pmpivmt(NTMAX)
      save pmpivmt

      if ( flag .eq. GLB_GET ) then
         dvalue = pmpivmt(imolty)
      elseif ( flag .eq. GLB_SET ) then
         pmpivmt(imolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            pmpivmt(imtemp) = dvalue
         enddo
      else
         write(6,*) 'PMPIVMT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmpivot(flag,dvalue)
c     ******************************************************************
c     * the probability sum for performing the pivot move              *
c     *                                                                *
c     * originally written 08-11-2008 by M.G. Martin                   *
c     * last modified 08-11-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmpivot
      save pmpivot

      if ( flag .eq. GLB_GET ) then
         dvalue = pmpivot
      elseif ( flag .eq. GLB_SET ) then
         pmpivot = dvalue
      else
         write(6,*) 'PMPIVOT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmplane(flag,dvalue)
c     ******************************************************************
c     * the probability sum for performing the plane shift move        *
c     *                                                                *
c     * originally written 08-11-2008 by M.G. Martin                   *
c     * last modified 08-11-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmplane
      save pmplane

      if ( flag .eq. GLB_GET ) then
         dvalue = pmplane
      elseif ( flag .eq. GLB_SET ) then
         pmplane = dvalue
      else
         write(6,*) 'PMPLANE: invalid global flag',flag
         stop
      endif

      return
      end


      subroutine twh_pmplanebox(flag,ibox,dvalue)
c     ******************************************************************
c     * the probability sum for performing the plane shift move in each*
c     * simulation box                                                 *
c     *                                                                *
c     * originally written 08-11-2008 by M.G. Martin                   *
c     * last modified 08-11-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ibox
      double precision dvalue
c     --- local variables
      integer ibtemp
c     --- information variables
      double precision pmplanebox
      dimension pmplanebox(MAXBOX)
      save pmplanebox

      if ( flag .eq. GLB_GET ) then
         dvalue = pmplanebox(ibox)
      elseif ( flag .eq. GLB_SET ) then
         pmplanebox(ibox) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 1,ibox
            pmplanebox(ibtemp) = dvalue
         enddo
      else
         write(6,*) 'PMPLANEBOX: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmromt(flag,imolty,dvalue)
c     ******************************************************************
c     * the probability sum for performing the rotation move on each   *
c     * molecule type in the simulation                                *
c     *                                                                *
c     * originally written 09-24-2008 by M.G. Martin                   *
c     * last modified 09-24-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer imtemp
c     --- information variables
      double precision pmromt
      dimension pmromt(NTMAX)
      save pmromt

      if ( flag .eq. GLB_GET ) then
         dvalue = pmromt(imolty)
      elseif ( flag .eq. GLB_SET ) then
         pmromt(imolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            pmromt(imtemp) = dvalue
         enddo
      else
         write(6,*) 'PMROMT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmrotate(flag,dvalue)
c     ******************************************************************
c     * the probability sum for performing the rotations move          *
c     *                                                                *
c     * originally written 09-24-2008 by M.G. Martin                   *
c     * last modified 09-24-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmrotate
      save pmrotate

      if ( flag .eq. GLB_GET ) then
         dvalue = pmrotate
      elseif ( flag .eq. GLB_SET ) then
         pmrotate = dvalue
      else
         write(6,*) 'PMROTATE: invalid global flag',flag
         stop
      endif

      return
      end


      subroutine twh_pmrow(flag,dvalue)
c     ******************************************************************
c     * the probability sum for performing the row shift move          *
c     *                                                                *
c     * originally written 09-23-2008 by M.G. Martin                   *
c     * last modified 09-23-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmrow
      save pmrow

      if ( flag .eq. GLB_GET ) then
         dvalue = pmrow
      elseif ( flag .eq. GLB_SET ) then
         pmrow = dvalue
      else
         write(6,*) 'PMROW: invalid global flag',flag
         stop
      endif

      return
      end


      subroutine twh_pmrowbox(flag,ibox,dvalue)
c     ******************************************************************
c     * the probability sum for performing the row shift move in each  *
c     * simulation box                                                 *
c     *                                                                *
c     * originally written 09-23-2008 by M.G. Martin                   *
c     * last modified 09-23-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ibox
      double precision dvalue
c     --- local variables
      integer ibtemp
c     --- information variables
      double precision pmrowbox
      dimension pmrowbox(MAXBOX)
      save pmrowbox

      if ( flag .eq. GLB_GET ) then
         dvalue = pmrowbox(ibox)
      elseif ( flag .eq. GLB_SET ) then
         pmrowbox(ibox) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 1,ibox
            pmrowbox(ibtemp) = dvalue
         enddo
      else
         write(6,*) 'PMROWBOX: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmtamt(flag,imolty,dvalue)
c     ******************************************************************
c     * the probability sum for single atom translation move           *
c     *                                                                *
c     * originally written 09-23-2008 by M.G. Martin                   *
c     * last modified 09-23-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer imtemp
c     --- information variables
      double precision pmtamt
      dimension pmtamt(NTMAX)
      save pmtamt

      if ( flag .eq. GLB_GET ) then
         dvalue = pmtamt(imolty)
      elseif ( flag .eq. GLB_SET ) then
         pmtamt(imolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            pmtamt(imtemp) = dvalue
         enddo
      else
         write(6,*) 'PMTAMT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmtcmt(flag,imolty,dvalue)
c     ******************************************************************
c     * the probability of performing a COM translation move on each   *
c     * molecule type in the simulation                                *
c     *                                                                *
c     * originally written 09-24-2008 by M.G. Martin                   *
c     * last modified 09-24-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer imtemp
c     --- information variables
      double precision pmtcmt
      dimension pmtcmt(NTMAX)
      save pmtcmt

      if ( flag .eq. GLB_GET ) then
         dvalue = pmtcmt(imolty)
      elseif ( flag .eq. GLB_SET ) then
         pmtcmt(imolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            pmtcmt(imtemp) = dvalue
         enddo
      else
         write(6,*) 'PMTCMT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmtraat(flag,dvalue)
c     ******************************************************************
c     * the probability sum for single atom translation move           *
c     *                                                                *
c     * originally written 09-23-2008 by M.G. Martin                   *
c     * last modified 09-23-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmtraat
      save pmtraat

      if ( flag .eq. GLB_GET ) then
         dvalue = pmtraat
      elseif ( flag .eq. GLB_SET ) then
         pmtraat = dvalue
      else
         write(6,*) 'PMTRAAT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmtracm(flag,dvalue)
c     ******************************************************************
c     * the probability sum for COM translation move                   *
c     *                                                                *
c     * originally written 09-23-2008 by M.G. Martin                   *
c     * last modified 09-23-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmtracm
      save pmtracm

      if ( flag .eq. GLB_GET ) then
         dvalue = pmtracm
      elseif ( flag .eq. GLB_SET ) then
         pmtracm = dvalue
      else
         write(6,*) 'PMTRACM: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmuvtcbmt(flag,imolty,dvalue)
c     ******************************************************************
c     * chance of performing the grand canonical insertion/deletion    *
c     * move on each molecule type in the simulation                   *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty
      double precision dvalue
c     --- local variables
      integer itemp
c     --- information variables
      double precision pmuvtcbmt
      dimension pmuvtcbmt(NTMAX)
      save pmuvtcbmt

      if ( flag .eq. GLB_GET ) then
         dvalue = pmuvtcbmt(imolty)
      elseif ( flag .eq. GLB_SET ) then
         pmuvtcbmt(imolty) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do itemp = 1,imolty
            pmuvtcbmt(itemp) = dvalue
         enddo
      else
         write(6,*) 'PMUVTCBMT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmuvtcbswap(flag,dvalue)
c     ******************************************************************
c     * the probability sum for the grand canonical insertion/deletion *
c     * move                                                           *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmuvtcbswap
      save pmuvtcbswap

      if ( flag .eq. GLB_GET ) then
         dvalue = pmuvtcbswap
      elseif ( flag .eq. GLB_SET ) then
         pmuvtcbswap = dvalue
      else
         write(6,*) 'PMUVTCBSWAP: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmvlpr(flag,iboxpair,dvalue)
c     ******************************************************************
c     * the probability sum for performing the volume move on each box *
c     * or box pair                                                    *
c     *                                                                *
c     * originally written 06-12-2006 by M.G. Martin                   *
c     * last modified 06-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,iboxpair
      double precision dvalue
c     --- local variables
      integer iptemp
c     --- information variables
      double precision pmvlpr
      dimension pmvlpr(MAXPBOX)
      save pmvlpr

      if ( flag .eq. GLB_SET ) then
         pmvlpr(iboxpair) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = pmvlpr(iboxpair)
      elseif ( flag .eq. GLB_INIT ) then
         do iptemp = 1,iboxpair
            pmvlpr(iptemp) = dvalue
         enddo
      else
         write(6,*) 'PMVLPR: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pmvol(flag,dvalue)
c     ******************************************************************
c     * the probability sum for the isotropic volume move              *
c     *                                                                *
c     * originally written 06-12-2006 by M.G. Martin                   *
c     * last modified 06-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision pmvol
      save pmvol

      if ( flag .eq. GLB_SET ) then
         pmvol = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = pmvol
      else
         write(6,*) 'PMVOL: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_polarizability(flag,itype,dvalue)
c     ******************************************************************
c     * the atomic (or molecular) polarizability for each type of atom *
c     * in the simulation.  This is used for some of the mixing rules  *
c     *                                                                *
c     * originally written 04-18-2006 by M.G. Martin                   *
c     * last modified 04-18-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype
      double precision dvalue
c     --- local variables
      integer ittemp
c     --- information variables
      double precision polarizability
      dimension polarizability(NNTYPE)
      save polarizability

      if ( flag .eq. GLB_SET ) then
         polarizability(itype) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = polarizability(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            polarizability(ittemp) = dvalue
         enddo
      else
         write(6,*) 'POLARIZABILITY: invalid global flag',flag
         stop
      endif

      return
      end

#ifdef USEOPENKIM
      subroutine twh_pointerkim(flag,ivalue)
c     ******************************************************************
c     * special integer pointer for use with openkim                   *
c     *                                                                *
c     * originally written 01-10-2012 by M.G. Martin                   *
c     * last modified 01-10-2012 by M.G. Martin                        *
c     ******************************************************************
      use KIMservice
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      integer (kind=kim_intptr) ivalue
c     --- local variables
      integer (kind=kim_intptr) pointerkim
      save pointerkim
      if ( flag .eq. GLB_GET ) then
         ivalue = pointerkim
      elseif ( flag .eq. GLB_SET ) then
         pointerkim = ivalue
      else
         write(6,*) 'POINTKIM: invalid global flag',flag
         stop
      endif

      return
      end
#endif

      subroutine twh_potentialstyle(flag,ivalue)
c     ******************************************************************
c     * the integer the specifies whether we are using an internal or  *
c     * external potentialstyle
c     *                                                                *
c     * originally written 09-25-2008 by M.G. Martin                   *
c     * last modified 09-25-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer potentialstyle
      save potentialstyle

      if ( flag .eq. GLB_GET ) then
         ivalue = potentialstyle
      elseif ( flag .eq. GLB_SET ) then
         potentialstyle = ivalue
      else
         write(6,*) 'POTENTIALSTYLE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pressurefreq(flag,ivalue)
c     ******************************************************************
c     * the frequency for computing the pressure using the virial      *
c     *                                                                *
c     * originally written 12-09-2007 by M.G. Martin                   *
c     * last modified 12-09-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer pressurefreq
      save pressurefreq

      if ( flag .eq. GLB_GET ) then
         ivalue = pressurefreq
      elseif ( flag .eq. GLB_SET ) then
         pressurefreq = ivalue
      else
         write(6,*) 'PRESSUREFREQ: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_pressurestyle(flag,ivalue)
c     ******************************************************************
c     * integer code for the pressure computation style                *
c     *                                                                *
c     * originally written 10-11-2008 by M.G. Martin                   *
c     * last modified 10-11-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer pressurestyle
      save pressurestyle

      if ( flag .eq. GLB_GET ) then
         ivalue = pressurestyle
      elseif ( flag .eq. GLB_SET ) then
         pressurestyle = ivalue
      else
         write(6,*) 'PRESSURESTYLE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_printfreq(flag,ivalue)
c     ******************************************************************
c     * the frequency that runtime information is printed to stdout    *
c     *                                                                *
c     * originally written 12-09-2007 by M.G. Martin                   *
c     * last modified 12-09-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer printfreq
      save printfreq

      if ( flag .eq. GLB_GET ) then
         ivalue = printfreq
      elseif ( flag .eq. GLB_SET ) then
         printfreq = ivalue
      else
         write(6,*) 'PRINTFREQ: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_qbasevalue(flag,itype,dvalue)
c     ******************************************************************
c     * the default charge on each atom before applying the bond       *
c     * increments to shift charge around.  Normally 0.0 for most atoms*
c     * but nonzero for ions                                           *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype
      double precision dvalue
c     --- local variables
      integer ittemp
c     --- information variables
      double precision qbasevalue
      dimension qbasevalue(NNTYPE)
      save qbasevalue

      if ( flag .eq. GLB_SET ) then
         qbasevalue(itype) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = qbasevalue(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            qbasevalue(ittemp) = dvalue
         enddo
      else
         write(6,*) 'QBASEVALUE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_qbiff(flag,itype,cvalue)
c     ******************************************************************
c     * the name of the forcefield for each charge increment type      *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype
      character*(*) cvalue
c     --- local variables
      integer ittemp
c     --- information variables
      character*(FFNAMELEN) qbiff
      dimension qbiff(TBIMAX)
      save qbiff

      if ( flag .eq. GLB_SET ) then
         qbiff(itype) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = qbiff(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            qbiff(ittemp) = cvalue
         enddo
      else
         write(6,*) 'QBIFF: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_qbinames(flag,itype,index,cvalue)
c     ******************************************************************
c     * the name of the atoms for each charge increment type           *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,index
      character*(*) cvalue
c     --- local variables
      integer ittemp,iitemp
c     --- information variables
      character*(FFNAMELEN) qbinames
      dimension qbinames(TBIMAX,2)
      save qbinames

      if ( flag .eq. GLB_SET ) then
         qbinames(itype,index) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = qbinames(itype,index)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            do iitemp = 1,index
               qbinames(ittemp,iitemp) = cvalue
            enddo
         enddo
      else
         write(6,*) 'QBINAMES: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_qbivalue(flag,itype,dvalue)
c     ******************************************************************
c     * the magnitude of the charge increment for each type            *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype
      double precision dvalue
c     --- local variables
      integer ittemp
c     --- information variables
      double precision qbivalue
      dimension qbivalue(TBIMAX)
      save qbivalue

      if ( flag .eq. GLB_SET ) then
         qbivalue(itype) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = qbivalue(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            qbivalue(ittemp) = dvalue
         enddo
      else
         write(6,*) 'QBIVALUE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_qqatom(flag,itype,iunit,dvalue)
c     ******************************************************************
c     * the charge on each atom of each molecule type in the system    *
c     *                                                                *
c     * originally written 11-02-2007 by M.G. Martin                   *
c     * last modified 11-02-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,iunit
      double precision dvalue
c     --- local variables
      integer ittemp,iutemp
c     --- information variables
      double precision qqatom
      dimension qqatom(NTMAX,NUMAX)
      save qqatom

      if ( flag .eq. GLB_GET ) then
         dvalue = qqatom(itype,iunit)
      elseif ( flag .eq. GLB_SET ) then
         qqatom(itype,iunit) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            do iutemp = 1,iunit
               qqatom(ittemp,iutemp) = dvalue
            enddo
         enddo
      else
         write(6,*) 'QQATOM: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_quantum_gsenergy(flag,itype,dvalue)
c     ******************************************************************
c     * the ground state energy of each molecule type                  *
c     *                                                                *
c     * originally written 05-12-2006 by M.G. Martin                   *
c     * last modified 05-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype
      double precision dvalue
c     --- local variables
      integer ittemp
c     --- information variables
      double precision quantum_gsenergy
      dimension quantum_gsenergy(NTMAX)
      save quantum_gsenergy

      if ( flag .eq. GLB_SET ) then
         quantum_gsenergy(itype) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = quantum_gsenergy(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            quantum_gsenergy(ittemp) = dvalue
         enddo
      else
         write(6,*) 'QUANTUM_GSENERGY: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_radial_pressure_delta(flag,dvalue)
c     ******************************************************************
c     * the size of the shell to use when computing the radial pressure*
c     *                                                                *
c     * originally written 10-11-2008 by M.G. Martin                   *
c     * last modified 10-11-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision radial_pressure_delta
      save radial_pressure_delta

      if ( flag .eq. GLB_GET ) then
         dvalue = radial_pressure_delta
      elseif ( flag .eq. GLB_SET ) then
         radial_pressure_delta = dvalue
      else
         write(6,*) 'RADIAL_PRESSURE_DELTA: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_random_allow_restart(flag,lvalue)
c     ******************************************************************
c     * logical that controls whether initialization as a continuation *
c     * of the last run is allowed                                     *
c     *                                                                *
c     * originally written 11-13-2006 by M.G. Martin                   *
c     * last modified 11-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical random_allow_restart
      save random_allow_restart

      if ( flag .eq. GLB_SET ) then
         random_allow_restart = lvalue
      elseif ( flag .eq. GLB_GET ) then
         lvalue = random_allow_restart
      else
         write(6,*) 'RANDOM_ALLOW_RESTART: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_random_luxlevel(flag,ivalue)
c     ******************************************************************
c     * the random number generator luxury level                       *
c     *                                                                *
c     * originally written 11-13-2006 by M.G. Martin                   *
c     * last modified 11-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      integer*4 ivalue
c     --- information variables
      integer*4 random_luxlevel
      save random_luxlevel

      if ( flag .eq. GLB_SET ) then
         random_luxlevel = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = random_luxlevel
      else
         write(6,*) 'RANDOM_LUXLEVEL: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_random_seed(flag,ivalue)
c     ******************************************************************
c     * the random number generator single integer seed                *
c     *                                                                *
c     * originally written 11-13-2006 by M.G. Martin                   *
c     * last modified 08-16-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      integer*4 ivalue
c     --- information variables
      integer*4 random_seed
      save random_seed

      if ( flag .eq. GLB_SET ) then
         random_seed = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = random_seed
      else
         write(6,*) 'RANDOM_SEED: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_rcelect(flag,ibox,dvalue)
c     ******************************************************************
c     * the electrostatic cutoff in each simulation box                *
c     *                                                                *
c     * originally written 05-04-2008 by M.G. Martin                   *
c     * last modified 08-02-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ibox
      double precision dvalue
c     --- local variables
      integer ibtemp
c     --- information variables
      double precision rcelect
      dimension rcelect(MAXBOX)
      save rcelect

      if ( flag .eq. GLB_GET ) then
         dvalue = rcelect(ibox)
      elseif ( flag .eq. GLB_SET ) then
         rcelect(ibox) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ibtemp = 1,ibox
            rcelect(ibtemp) = dvalue
         enddo
      else
         write(6,*) 'RCELECT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_rcut(flag,dvalue)
c     ******************************************************************
c     * vdw cutoff                                                     *
c     *                                                                *
c     * originally written 10-06-2008 by M.G. Martin                   *
c     * last modified 10-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision rcut
      save rcut

      if ( flag .eq. GLB_GET ) then
         dvalue = rcut
      elseif ( flag .eq. GLB_SET ) then
         rcut = dvalue
      else
         write(6,*) 'RCUT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_rcutin(flag,dvalue)
c     ******************************************************************
c     * the inner cutoff used in CBMC moves                            *
c     *                                                                *
c     * originally written 10-06-2008 by M.G. Martin                   *
c     * last modified 10-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision rcutin
      save rcutin

      if ( flag .eq. GLB_GET ) then
         dvalue = rcutin
      elseif ( flag .eq. GLB_SET ) then
         rcutin = dvalue
      else
         write(6,*) 'RCUTIN: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_rcutinsq(flag,dvalue)
c     ******************************************************************
c     * the square of the inner cutoff used in CBMC moves              *
c     *                                                                *
c     * originally written 10-06-2008 by M.G. Martin                   *
c     * last modified 10-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision rcutinsq
      save rcutinsq

      if ( flag .eq. GLB_GET ) then
         dvalue = rcutinsq
      elseif ( flag .eq. GLB_SET ) then
         rcutinsq = dvalue
      else
         write(6,*) 'RCUTINSQ: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_rcutsq(flag,dvalue)
c     ******************************************************************
c     * square of the vdw cutoff                                       *
c     *                                                                *
c     * originally written 10-06-2008 by M.G. Martin                   *
c     * last modified 10-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision rcutsq
      save rcutsq

      if ( flag .eq. GLB_GET ) then
         dvalue = rcutsq
      elseif ( flag .eq. GLB_SET ) then
         rcutsq = dvalue
      else
         write(6,*) 'RCUTSQ: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_restartfreq(flag,ivalue)
c     ******************************************************************
c     * frequency for outputting a full restart file                   *
c     *                                                                *
c     * originally written 12-22-2007 by M.G. Martin                   *
c     * last modified 12-22-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer restartfreq
      save restartfreq

      if ( flag .eq. GLB_GET ) then
         ivalue = restartfreq
      elseif ( flag .eq. GLB_SET ) then
         restartfreq = ivalue
      else
         write(6,*) 'RESTARTFREQ: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_random_code(flag,ivalue)
c     ******************************************************************
c     * integer code for the random number generator                   *
c     *                                                                *
c     * originally written 08-12-2011 by M.G. Martin                   *
c     * last modified 08-15-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer random_code
      save random_code

      if ( flag .eq. GLB_GET ) then
         ivalue = random_code
      elseif ( flag .eq. GLB_SET ) then
         random_code = ivalue
      else
         write(6,*) 'RANDOM_CODE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_rowwidth(flag,dvalue)
c     ******************************************************************
c     * the width of the row for the row shift move                    *
c     *                                                                *
c     * originally written 09-23-2008 by M.G. Martin                   *
c     * last modified 09-23-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision rowwidth
      save rowwidth

      if ( flag .eq. GLB_SET ) then
         rowwidth = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = rowwidth
      else
         write(6,*) 'ROWWIDTH: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_run_name(flag,cvalue)
c     ******************************************************************
c     * TMMC fileprefix for tmmc output files                          *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      character*(*) cvalue
c     --- information variables
      character*50  run_name
      save run_name

      if ( flag .eq. GLB_GET ) then
         cvalue = run_name
      elseif ( flag .eq. GLB_SET ) then
         run_name = cvalue
      else
         write(6,*) 'RUN_NAME: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_runoutput(flag,cvalue)
c     ******************************************************************
c     * the style of output during the simulation run                  *
c     *                                                                *
c     * originally written 11-08-2006 by M.G. Martin                   *
c     * last modified 11-08-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      character*(*) cvalue
c     --- information variables
      character*20 runoutput
      save runoutput

      if ( flag .eq. GLB_SET ) then
         runoutput = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = runoutput
      else
         write(6,*) 'RUNOUTPUT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_rmcell(flag,ipairbox,idim,jdim,dvalue)
c     ******************************************************************
c     * the maximum volume displacement for hmatrix value in each box  *
c     *                                                                *
c     * originally written 12-23-2007 by M.G. Martin                   *
c     * last modified 12-23-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ipairbox,idim,jdim
      double precision dvalue
c     --- local variables
      integer iptemp,iditemp,idjtemp
c     --- information variables
      double precision rmcell
      dimension rmcell(MAXPBOX,3,3)
      save rmcell

      if ( flag .eq. GLB_GET ) then
         dvalue = rmcell(ipairbox,idim,jdim)
      elseif ( flag .eq. GLB_SET ) then
         rmcell(ipairbox,idim,jdim) = dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do iptemp = 1,ipairbox
            do iditemp = 1,idim
               do idjtemp = 1,jdim
                  rmcell(iptemp,iditemp,idjtemp) = dvalue
               enddo
            enddo
         enddo
      elseif ( flag .eq. GLB_SCALE ) then
         rmcell(ipairbox,idim,jdim) = rmcell(ipairbox,idim,jdim)*dvalue
      else
         write(6,*) 'RMCELL: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_rmin(flag,dvalue)
c     ******************************************************************
c     * the inner hard core for each atom                              *
c     *                                                                *
c     * originally written 10-06-2008 by M.G. Martin                   *
c     * last modified 10-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision rmin
      save rmin

      if ( flag .eq. GLB_GET ) then
         dvalue = rmin
      elseif ( flag .eq. GLB_SET ) then
         rmin = dvalue
      else
         write(6,*) 'RMIN: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_rminsq(flag,dvalue)
c     ******************************************************************
c     * the square of the inner hard core for each atom                *
c     *                                                                *
c     * originally written 10-06-2008 by M.G. Martin                   *
c     * last modified 10-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision rminsq
      save rminsq

      if ( flag .eq. GLB_GET ) then
         dvalue = rminsq
      elseif ( flag .eq. GLB_SET ) then
         rminsq = dvalue
      else
         write(6,*) 'RMINSQ: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_scalecut(flag,dvalue)
c     ******************************************************************
c     * value used to help rescale the rosenbluth weights              *
c     *                                                                *
c     * originally written 10-06-2008 by M.G. Martin                   *
c     * last modified 10-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision scalecut
      save scalecut

      if ( flag .eq. GLB_GET ) then
         dvalue = scalecut
      elseif ( flag .eq. GLB_SET ) then
         scalecut = dvalue
      else
         write(6,*) 'SCALECUT: unknown flag',flag
         stop
      endif

      return
      end

      subroutine twh_scaled_npair(flag,ivalue)
c     ******************************************************************
c     * the number of pairs for scaled potentials                      *
c     *                                                                *
c     * originally written 12-22-2007 by M.G. Martin                   *
c     * last modified 12-22-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer scaled_npair 
      save scaled_npair

      if ( flag .eq. GLB_GET ) then
         ivalue = scaled_npair
      elseif ( flag .eq. GLB_SET ) then
         scaled_npair = ivalue
      else
         write(6,*) 'SCALED_NPAIR: unknown flag',flag
         stop
      endif

      return
      end

      subroutine twh_scaled_pair_list(flag,itype,count,cvalue)
c     ******************************************************************
c     * the list of scaled pairs                                       *
c     *                                                                *
c     * originally written 12-22-2007 by M.G. Martin                   *
c     * last modified 12-22-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,count
      character*(*) cvalue
c     --- local variables
      integer ittemp,ictemp
c     --- information variables
      character*10 scaled_pair_list 
      dimension scaled_pair_list(CROSSTYPEMAX,2)
      save scaled_pair_list

      if ( flag .eq. GLB_GET ) then
         cvalue = scaled_pair_list(itype,count)
      elseif ( flag .eq. GLB_SET ) then
         scaled_pair_list(itype,count) = cvalue
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            do ictemp = 1,count
               scaled_pair_list(ittemp,ictemp) = cvalue
            enddo
         enddo
      else
         write(6,*) 'SCALED_PAIR_LIST: unknown flag',flag
         stop
      endif

      return
      end

      subroutine twh_scalehilog(flag,dvalue)
c     ******************************************************************
c     * value used to scale the rosenbluth weight                      *
c     *                                                                *
c     * originally written 10-06-2008 by M.G. Martin                   *
c     * last modified 10-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision scalehilog
      save scalehilog

      if ( flag .eq. GLB_GET ) then
         dvalue = scalehilog
      elseif ( flag .eq. GLB_SET ) then
         scalehilog = dvalue
      else
         write(6,*) 'SCALEHILOG: unknown flag',flag
         stop
      endif

      return
      end

      subroutine twh_scalelolog(flag,dvalue)
c     ******************************************************************
c     * value used to scale the rosenbluth weight                      *
c     *                                                                *
c     * originally written 10-06-2008 by M.G. Martin                   *
c     * last modified 10-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision scalelolog
      save scalelolog

      if ( flag .eq. GLB_GET ) then
         dvalue = scalelolog
      elseif ( flag .eq. GLB_SET ) then
         scalelolog = dvalue
      else
         write(6,*) 'SCALELOLOG: unknown flag',flag
         stop
      endif

      return
      end

      subroutine twh_scaling_style(flag,ivalue)
c     ******************************************************************
c     * the style of scaling for the scaled potential                  *
c     *                                                                *
c     * originally written 12-20-2007 by M.G. Martin                   *
c     * last modified 12-20-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer scaling_style
      save scaling_style

      if ( flag .eq. GLB_GET ) then
         ivalue = scaling_style
      elseif ( flag .eq. GLB_SET ) then
         scaling_style = ivalue
      else
         write(6,*) 'SCALING_STYLE: unknown flag',flag
         stop
      endif

      return
      end

      subroutine twh_scp_eval_derivative(flag,lvalue)
c     ******************************************************************
c     * logical scp_eval_derivative                                    *
c     * controls whether the derivative of energy with respect to      *
c     * lambda should be evaluated.                                    *
c     *                                                                *
c     * originally written 01-23-2008 by M.G. Martin                   *
c     * last modified 01-23-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical scp_eval_derivative 
      save scp_eval_derivative

      if ( flag .eq. GLB_GET ) then
         lvalue = scp_eval_derivative
      elseif ( flag .eq. GLB_SET ) then
         scp_eval_derivative = lvalue
      else
         write(6,*) 'SCP_EVAL_DERIVATIVE: unknown flag',flag
         stop
      endif

      return
      end

      subroutine twh_scp_num_foreign_lambda(flag,ivalue)
c     ******************************************************************
c     * number of foreign lambda terms for scaled potential            *
c     * number of foreign lambdas at which to evaluate energy.         *
c     * Is set to 0 when no foreign energy evaluations                 *
c     *                                                                *
c     * originally written 01-23-2008 by M.G. Martin                   *
c     * last modified 01-23-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer scp_num_foreign_lambda
      save scp_num_foreign_lambda

      if ( flag .eq. GLB_GET ) then
         ivalue = scp_num_foreign_lambda
      elseif ( flag .eq. GLB_SET ) then
         scp_num_foreign_lambda = ivalue
      else
         write(6,*) 'SCP_NUM_FOREIGN_LAMBDA: unknown flag',flag
         stop
      endif

      return
      end

      subroutine twh_scp_scale_coulomb(flag,lvalue)
c     ******************************************************************
c     * logical for whether coulomb terms are scaled                   *
c     *                                                                *
c     * originally written 01-23-2008 by M.G. Martin                   *
c     * last modified 01-23-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical scp_scale_coulomb 
      save scp_scale_coulomb

      if ( flag .eq. GLB_GET ) then
         lvalue = scp_scale_coulomb
      elseif ( flag .eq. GLB_SET ) then
         scp_scale_coulomb = lvalue
      else
         write(6,*) 'SCP_SCALE_COULOMB: unknown flag',flag
         stop
      endif

      return
      end

      subroutine twh_scp_scale_twobody(flag,lvalue)
c     ******************************************************************
c     * logical for whether twobody terms are scaled                   *
c     *                                                                *
c     * originally written 01-23-2008 by M.G. Martin                   *
c     * last modified 01-23-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical scp_scale_twobody
      save scp_scale_twobody

      if ( flag .eq. GLB_GET ) then
         lvalue = scp_scale_twobody
      elseif ( flag .eq. GLB_SET ) then
         scp_scale_twobody = lvalue
      else
         write(6,*) 'SCP_SCALE_TWOBODY: unknown flag',flag
         stop
      endif

      return
      end

      subroutine twh_sdevbena(flag,dvalue)
c     ******************************************************************
c     * standard deviation scaling for the bending A CBMC generation   *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision sdevbena
      save sdevbena

      if ( flag .eq. GLB_GET ) then
         dvalue = sdevbena
      elseif ( flag .eq. GLB_SET ) then
         sdevbena = dvalue
      else
         write(6,*) 'SDEVBENA: unknown flag'
         write(6,*) flag
         stop
      endif

      return
      end

      subroutine twh_sdevbenb(flag,dvalue)
c     ******************************************************************
c     * standard deviation scaling for the bending B CBMC generation   *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision sdevbenb
      save sdevbenb

      if ( flag .eq. GLB_GET ) then
         dvalue = sdevbenb
      elseif ( flag .eq. GLB_SET ) then
         sdevbenb = dvalue
      else
         write(6,*) 'SDEVBENB: unknown flag'
         write(6,*) flag
         stop
      endif

      return
      end

      subroutine twh_sdevtor(flag,dvalue)
c     ******************************************************************
c     * standard deviation scaling for the torsion CBMC generation     *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision sdevtor
      save sdevtor

      if ( flag .eq. GLB_GET ) then
         dvalue = sdevtor
      elseif ( flag .eq. GLB_SET ) then
         sdevtor = dvalue
      else
         write(6,*) 'SDEVTOR: unknown flag'
         write(6,*) flag
         stop
      endif

      return
      end

      subroutine twh_sdevvib(flag,dvalue)
c     ******************************************************************
c     * standard deviation scaling for the vibration CBMC generation   *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision sdevvib
      save sdevvib

      if ( flag .eq. GLB_GET ) then
         dvalue = sdevvib
      elseif ( flag .eq. GLB_SET ) then
         sdevvib = dvalue
      else
         write(6,*) 'SDEVVIB: unknown flag'
         write(6,*) flag
         stop
      endif

      return
      end

      subroutine twh_solvation_style(flag,cvalue)
c     ******************************************************************
c     * the style of solvation potential                               *
c     *                                                                *
c     * originally written 01-12-2007 by M.G. Martin                   *
c     * last modified 01-12-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      character*(*) cvalue
c     --- information variables
      character*20 solvation_style
      save solvation_style

      if ( flag .eq. GLB_SET ) then
         solvation_style = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = solvation_style
      else
         write(6,*) 'SOLVATION_STYLE: unknown flag'
         write(6,*) flag
         stop
      endif

      return
      end

      subroutine twh_solvation_type(flag,cvalue)
c     ******************************************************************
c     * the type of solvation potential                                *
c     *                                                                *
c     * originally written 01-12-2007 by M.G. Martin                   *
c     * last modified 01-12-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      character*(*) cvalue
c     --- information variables
      character*20 solvation_type
      save solvation_type

      if ( flag .eq. GLB_SET ) then
         solvation_type = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = solvation_type
      else
         write(6,*) 'SOLVATION_TYPE: unknown flag'
         write(6,*) flag
         stop
      endif

      return
      end

      subroutine twh_stepstyle(flag,cvalue)
c     ******************************************************************
c     * the style to use for the Monte Carlo move steps                *
c     *                                                                *
c     * originally written 12-09-2007 by M.G. Martin                   *
c     * last modified 12-09-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      character*(*) cvalue
c     --- information variables
      character*8 stepstyle
      save stepstyle

      if ( flag .eq. GLB_GET ) then
         cvalue = stepstyle
      elseif ( flag .eq. GLB_SET ) then
         stepstyle = cvalue
      else
         write(6,*) 'STEPSTYLE: unknown flag'
         write(6,*) flag
         stop
      endif

      return
      end

      subroutine twh_table_npair(flag,itype,ivalue)
c     ******************************************************************
c     * number of pairs in the tabular potentials                      *
c     *                                                                *
c     * originally written 09-25-2008 by M.G. Martin                   *
c     * last modified 09-25-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer ittemp
c     --- information variables
      integer table_npair
      dimension table_npair((MAXTABTYPE*MAXTABTYPE+MAXTABTYPE)/2)
      save table_npair

      if ( flag .eq. GLB_GET ) then
         ivalue = table_npair(itype)
      elseif ( flag .eq. GLB_SET ) then
         table_npair(itype) = ivalue
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            table_npair(ittemp) = ivalue
         enddo
      else
         write(6,*) 'TABLE_NPAIR: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_table_pair(flag,idim,itype,icount,dvalue)
c     ******************************************************************
c     * coefficients for the tabulated pair potentials                 *
c     *                                                                *
c     * originally written 09-25-2008 by M.G. Martin                   *
c     * last modified 09-26-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,idim,itype,icount
      double precision dvalue
c     --- local variables
      integer idtemp,ittemp,ictemp
c     --- information variables
      double precision table_pair
      dimension table_pair(3,((MAXTABTYPE*MAXTABTYPE+MAXTABTYPE)/2)
     &,MAXTABLE)
      save table_pair

      if ( flag .eq. GLB_GET ) then
         dvalue = table_pair(idim,itype,icount)
      elseif ( flag .eq. GLB_SET ) then
         table_pair(idim,itype,icount) = dvalue
      elseif ( flag .eq. GLB_SCALE ) then
         table_pair(idim,itype,icount)
     &        = table_pair(idim,itype,icount)*dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do idtemp = 1,idim
            do ittemp = 1,itype
               do ictemp = 1,icount
                  table_pair(idtemp,ittemp,ictemp) = dvalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'TABLE_PAIR: invalid global flag',flag
         stop
      endif

      return
      end


      subroutine twh_tacell(flag,dvalue)
c     ******************************************************************
c     * the target acceptance rate for cell volume moves               *
c     *                                                                *
c     * originally written 12-23-2007 by M.G. Martin                   *
c     * last modified 12-23-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision tacell
      save tacell

      if ( flag .eq. GLB_GET ) then
         dvalue = tacell
      elseif ( flag .eq. GLB_SET ) then
         tacell = dvalue
      else
         write(6,*) 'TACELL: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_tarot(flag,dvalue)
c     ******************************************************************
c     * the target acceptance rate for rotation move                   *
c     *                                                                *
c     * originally written 09-24-2008 by M.G. Martin                   *
c     * last modified 09-24-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision tarot
      save tarot

      if ( flag .eq. GLB_GET ) then
         dvalue = tarot
      elseif ( flag .eq. GLB_SET ) then
         tarot = dvalue
      else
         write(6,*) 'TAROT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_tatraa(flag,dvalue)
c     ******************************************************************
c     * the target acceptance rate for single-atom translation moves   *
c     *                                                                *
c     * originally written 09-22-2008 by M.G. Martin                   *
c     * last modified 09-22-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision tatraa
      save tatraa

      if ( flag .eq. GLB_GET ) then
         dvalue = tatraa
      elseif ( flag .eq. GLB_SET ) then
         tatraa = dvalue
      else
         write(6,*) 'TATRAA: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_tatrac(flag,dvalue)
c     ******************************************************************
c     * the target acceptance rate for COM translation move            *
c     *                                                                *
c     * originally written 09-24-2008 by M.G. Martin                   *
c     * last modified 09-24-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision tatrac
      save tatrac

      if ( flag .eq. GLB_GET ) then
         dvalue = tatrac
      elseif ( flag .eq. GLB_SET ) then
         tatrac = dvalue
      else
         write(6,*) 'TATRAC: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_tavol(flag,dvalue)
c     ******************************************************************
c     * the target acceptance rate for isotropic volume moves          *
c     *                                                                *
c     * originally written 06-12-2006 by M.G. Martin                   *
c     * last modified 06-12-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision tavol
      save tavol

      if ( flag .eq. GLB_SET ) then
         tavol = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = tavol
      else
         write(6,*) 'TAVOL: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_temperature(flag,dvalue)
c     ******************************************************************
c     * the temperature in the simulation                              *
c     *                                                                *
c     * originally written 11-02-2007 by M.G. Martin                   *
c     * last modified 11-02-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- local variables
      double precision beta
c     --- information variables
      double precision temperature
      save temperature

      if ( flag .eq. GLB_GET ) then
         dvalue = temperature
      elseif ( flag .eq. GLB_SET ) then
         temperature = dvalue
c        --- also reset the inverse temperature
         beta = 1.0d0/temperature
         call twh_beta(GLB_SET,beta)
      else
         write(6,*) 'TEMPERATURE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_three_bond_fe_bias_style(flag,cvalue)
c     ******************************************************************
c     * bias style for the three-bond fixed-endpoint biasing           *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      character*(*) cvalue
c     --- information variables
      character*50 three_bond_fe_bias_style
      save three_bond_fe_bias_style

      if ( flag .eq. GLB_GET ) then
         cvalue = three_bond_fe_bias_style
      elseif ( flag .eq. GLB_SET ) then
         three_bond_fe_bias_style = cvalue
      else
         write(6,*)'THREE_BOND_FE_BIAS_STYLE: invalid global flag'
     &        ,flag
         stop
      endif

      return
      end

      subroutine twh_three_bond_bias_compute_freq(flag,ivalue)
c     ******************************************************************
c     * frequency for computing the current snapshot of three-bond     *
c     * distributions                                                  *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer three_bond_bias_compute_freq
      save three_bond_bias_compute_freq

      if ( flag .eq. GLB_GET ) then
         ivalue = three_bond_bias_compute_freq
      elseif ( flag .eq. GLB_SET ) then
         three_bond_bias_compute_freq = ivalue
      else
         write(6,*)'THREE_BOND_BIAS_COMPUTE_FREQ: invalid global flag'
     &        ,flag
         stop
      endif

      return
      end

      subroutine twh_three_bond_bias_info(flag,imolty,iunit,itor
     &     ,iparam,dvalue)
c     ******************************************************************
c     * information for the three-bond fixed-endpoint biasing          *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,iunit,itor,iparam
      double precision dvalue
c     --- local varaibles
      integer imtemp,iutemp,ittemp,iptemp
c     --- information variables
      double precision three_bond_bias_info
      dimension three_bond_bias_info(NTMAX,NUMAX
     &,MAXTOR,MAXTHREEBOND)
      save three_bond_bias_info

      if ( flag .eq. GLB_GET ) then
         dvalue = three_bond_bias_info(imolty,iunit,itor,iparam)
      elseif ( flag .eq. GLB_SET ) then
         three_bond_bias_info(imolty,iunit,itor,iparam) = dvalue
      elseif ( flag .eq. GLB_INCR ) then
         three_bond_bias_info(imolty,iunit,itor,iparam) =
     &        three_bond_bias_info(imolty,iunit,itor,iparam) + dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ittemp = 1,itor
                  do iptemp = 1,iparam
                     three_bond_bias_info(imtemp,iutemp,ittemp,iptemp)
     &                    = dvalue
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*)'THREE_BOND_BIAS_INFO: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_three_bond_bias_initial_value(flag,cvalue)
c     ******************************************************************
c     * initial value for the three-bond fixed-endpoint biasing        *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      character*(*) cvalue
c     --- information variables
      character*50 three_bond_bias_initial_value
      save three_bond_bias_initial_value

      if ( flag .eq. GLB_GET ) then
         cvalue = three_bond_bias_initial_value
      elseif ( flag .eq. GLB_SET ) then
         three_bond_bias_initial_value = cvalue
      else
         write(6,*)'THREE_BOND_BIAS_INITIAL_VALUE: invalid global flag'
     &        ,flag
         stop
      endif

      return
      end

      subroutine twh_three_bond_bias_old_fraction(flag,dvalue)
c     ******************************************************************
c     * fraction of the old distribution to keep when computing a new  *
c     * distribution for the three-bond biasing function               *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision three_bond_bias_old_fraction
      save three_bond_bias_old_fraction

      if ( flag .eq. GLB_GET ) then
         dvalue = three_bond_bias_old_fraction
      elseif ( flag .eq. GLB_SET ) then
         three_bond_bias_old_fraction = dvalue
      else
         write(6,*)'THREE_BOND_BIAS_OLD_FRACTION: invalid global flag'
     &        ,flag
         stop
      endif

      return
      end

      subroutine twh_three_bond_bias_sdev_multiplier(flag,dvalue)
c     ******************************************************************
c     * multiplier to rescale the compute standard deviation of the    *
c     * three-bond biasing function                                    *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision three_bond_bias_sdev_multiplier
      save three_bond_bias_sdev_multiplier

      if ( flag .eq. GLB_GET ) then
         dvalue = three_bond_bias_sdev_multiplier
      elseif ( flag .eq. GLB_SET ) then
         three_bond_bias_sdev_multiplier = dvalue
      else
         write(6,*)'THREE_BOND_BIAS_SDEV_MULT: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_three_bond_bias_update_freq(flag,ivalue)
c     ******************************************************************
c     * frequency to update the values of the three-bond biasing       *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer three_bond_bias_update_freq
      save three_bond_bias_update_freq

      if ( flag .eq. GLB_GET ) then
         ivalue = three_bond_bias_update_freq
      elseif ( flag .eq. GLB_SET ) then
         three_bond_bias_update_freq = ivalue
      else
         write(6,*)'THREE_BOND_BIAS_UPDATE_FREQ: invalid global flag'
     &        ,flag
         stop
      endif

      return
      end

      subroutine twh_tmmc_flag(flag,lvalue)
c     ******************************************************************
c     * logical to distinguish between gcmc and gctmmc                 *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      logical lvalue
c     --- information variables
      logical tmmc_flag
      save tmmc_flag

      if ( flag .eq. GLB_GET ) then
         lvalue = tmmc_flag
      elseif ( flag .eq. GLB_SET ) then
         tmmc_flag = lvalue
      else
         write(6,*) 'TMMC_FLAG: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_torcoeff(flag,itype,icoeff,dvalue)
c     ******************************************************************
c     * the coefficients for each type of torsion                      *
c     *                                                                *
c     * originally written 04-10-2006 by M.G. Martin                   *
c     * last modified 04-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,icoeff
      double precision dvalue
c     --- local variables
      integer ittemp,ictemp
c     --- information variables
      double precision torcoeff
      dimension torcoeff(TTORMAX,MINTORCOEFF:MAXTORCOEFF)
      save torcoeff

      if ( flag .eq. GLB_SET ) then
         torcoeff(itype,icoeff) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = torcoeff(itype,icoeff)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            do ictemp = MINTORCOEFF,icoeff
               torcoeff(ittemp,ictemp) = dvalue
            enddo
         enddo
      else
         write(6,*) 'TORCOEFF: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_torff(flag,itype,cvalue)
c     ******************************************************************
c     * name of the forcefield associated with this torsion type       *
c     *                                                                *
c     * originally written 04-10-2006 by M.G. Martin                   *
c     * last modified 05-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype
      character*(*) cvalue
c     --- local variables
      integer ittemp
c     --- information variables
      character*(FFNAMELEN) torff
      dimension torff(TTORMAX)
      save torff

c     --- safety check
      if ( itype .gt. TTORMAX ) then
         write(6,*) 'TORFF: itype exceeds TTORMAX'
         write(6,*) 'itype:',itype
         write(6,*) 'TTORMAX:',TTORMAX
         stop
      endif

      if ( flag .eq. GLB_SET ) then
         torff(itype) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = torff(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            torff(ittemp) = cvalue
         enddo
      else
         write(6,*) 'TORSTRING: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_tornames(flag,itype,isame,index,cvalue)
c     ******************************************************************
c     * name of the atoms associated with this torsion type            *
c     *                                                                *
c     * originally written 04-10-2006 by M.G. Martin                   *
c     * last modified 04-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,isame,index
      character*(*) cvalue
c     --- local variables
      integer ittemp,istemp,iitemp
c     --- information variables
      character*(FFNAMELEN) tornames
      dimension tornames(TTORMAX,MAXSAMET,MAXTNAME)
      save tornames

      if ( flag .eq. GLB_SET ) then
         tornames(itype,isame,index) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = tornames(itype,isame,index)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            do istemp = 1,isame
               do iitemp = 1,index
                  tornames(ittemp,istemp,iitemp) = cvalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'TORNAMES: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_torsionorder(flag,itype,cvalue)
c     ******************************************************************
c     * the torsion order for each angle type for matching in assemble *
c     *                                                                *
c     * originally written 12-13-2006 by M.G. Martin                   *
c     * last modified 12-19-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype
      character*(*) cvalue
c     --- local variables
      integer ittemp
c     --- information variables
      character*(15) torsionorder
      dimension torsionorder(TTORMAX)
      save torsionorder

      if ( flag .eq. GLB_SET ) then
         torsionorder(itype) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = torsionorder(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            torsionorder(ittemp) = cvalue
         enddo
      else
         write(6,*) 'TORSIONORDER: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_torstring(flag,itype,index,cvalue)
c     ******************************************************************
c     * character string to identify the correct atoms in dihedral     *
c     * angles for some class 2 forcefields                            *
c     *                                                                *
c     * originally written 04-06-2006 by M.G. Martin                   *
c     * last modified 04-06-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,index
      character*(*) cvalue
c     --- local variables
      integer ittemp,iitemp
c     --- information variables
      character*(FFNAMELEN) torstring
      dimension torstring(TTORMAX,4)
      save torstring

      if ( flag .eq. GLB_SET ) then
         torstring(itype,index) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = torstring(itype,index)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            do iitemp = 1,index
               torstring(ittemp,iitemp) = cvalue
            enddo
         enddo
      else
         write(6,*) 'TORSTRING: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_torstyle(flag,itype,ivalue)
c     ******************************************************************
c     * an integer representation of the torsion style                 *
c     *                                                                *
c     * originally written 04-05-2006 by M.G. Martin                   *
c     * last modified 04-05-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer ittemp
c     --- information variables
      integer torstyle
      dimension torstyle(TTORMAX)
      save torstyle

      if ( flag .eq. GLB_SET ) then
         torstyle(itype) = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = torstyle(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            torstyle(itype) = ivalue
         enddo
      else
         write(6,*) 'TORSTYLE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_trmaxdispfreq(flag,ivalue)
c     ******************************************************************
c     * the frequency for updating the maximum displacements for the   *
c     * translation and rotation moves                                 *
c     *                                                                *
c     * originally written 12-09-2007 by M.G. Martin                   *
c     * last modified 12-09-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer trmaxdispfreq
      save trmaxdispfreq

      if ( flag .eq. GLB_GET ) then
         ivalue = trmaxdispfreq
      elseif ( flag .eq. GLB_SET ) then
         trmaxdispfreq = ivalue
      else
         write(6,*) 'TRMAXDISPFREQ: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_two_bond_bias_compute_frequency(flag,ivalue)
c     ******************************************************************
c     * frequency for computing the new running averages for the       *
c     * two-bond fixed-endpoint biasing initial value                  *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer two_bond_bias_compute_frequency
      save two_bond_bias_compute_frequency

      if ( flag .eq. GLB_GET ) then
         ivalue = two_bond_bias_compute_frequency
      elseif ( flag .eq. GLB_SET ) then
         two_bond_bias_compute_frequency = ivalue
      else
         write(6,*) 'TWO_BOND_BIAS_COMPUTE_FREQ: invalid global flag'
     &        ,flag
         stop
      endif

      return
      end

      subroutine twh_two_bond_bias_initial_value(flag,cvalue)
c     ******************************************************************
c     * two-bond fixed-endpoint biasing initial value                  *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      character*(*) cvalue
c     --- information variables
      character*50 two_bond_bias_initial_value
      save two_bond_bias_initial_value

      if ( flag .eq. GLB_GET ) then
         cvalue = two_bond_bias_initial_value
      elseif ( flag .eq. GLB_SET ) then
         two_bond_bias_initial_value = cvalue
      else
         write(6,*) 'TWO_BOND_BIAS_INITIAL_VALUE: invalid global flag'
     &        ,flag
         stop
      endif

      return
      end

      subroutine twh_two_bond_fe_bias_style(flag,cvalue)
c     ******************************************************************
c     * fixed-endpoint biasing style for atoms connected by two bonds  *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      character*(*) cvalue
c     --- information variables
      character*50 two_bond_fe_bias_style
      save two_bond_fe_bias_style

      if ( flag .eq. GLB_GET ) then
         cvalue = two_bond_fe_bias_style
      elseif ( flag .eq. GLB_SET ) then
         two_bond_fe_bias_style = cvalue
      else
         write(6,*) 'TWO_BOND_FE_BIAS_STYLE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_two_bond_bias_info(flag,imolty,iunit,ibend,iparam
     &     ,dvalue)
c     ******************************************************************
c     * information about the two-bond biasing for every such term in  *
c     * each molecule type                                             *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,imolty,iunit,ibend,iparam
      double precision dvalue
c     --- local variables
      integer imtemp,iutemp,ibtemp,iptemp
c     --- information variables
      double precision two_bond_bias_info
      dimension two_bond_bias_info(NTMAX,NUMAX
     &,MAXBEND,MAXTWOBOND)
      save two_bond_bias_info

      if ( flag .eq. GLB_GET ) then
         dvalue = two_bond_bias_info(imolty,iunit,ibend,iparam)
      elseif ( flag .eq. GLB_SET ) then
         two_bond_bias_info(imolty,iunit,ibend,iparam) = dvalue
      elseif ( flag .eq. GLB_INCR ) then
         two_bond_bias_info(imolty,iunit,ibend,iparam) =
     &        two_bond_bias_info(imolty,iunit,ibend,iparam) + dvalue
      elseif ( flag .eq. GLB_INIT ) then
         do imtemp = 1,imolty
            do iutemp = 1,iunit
               do ibtemp = 1,ibend
                  do iptemp = 1,iparam
                     two_bond_bias_info(imtemp,iutemp,ibtemp,iptemp)
     &                    = dvalue
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'TWO_BOND_BIAS_INFO: invalid global flag'
     &        ,flag
         stop
      endif

      return
      end

      subroutine twh_two_bond_bias_old_fraction(flag,dvalue)
c     ******************************************************************
c     * fraction of the old distribution to retain when computing the  *
c     * new two-bond biasing functions                                 *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision two_bond_bias_old_fraction
      save two_bond_bias_old_fraction

      if ( flag .eq. GLB_GET ) then
         dvalue = two_bond_bias_old_fraction
      elseif ( flag .eq. GLB_SET ) then
         two_bond_bias_old_fraction = dvalue
      else
         write(6,*) 'TWO_BOND_BIAS_OLD_FRACTION: invalid global flag'
     &        ,flag
         stop
      endif

      return
      end

      subroutine twh_two_bond_bias_sdev_multiplier(flag,dvalue)
c     ******************************************************************
c     * multiplier to scale the compute standard deviations to create  *
c     * the distribution to use in the simulation                      *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision two_bond_bias_sdev_multiplier
      save two_bond_bias_sdev_multiplier

      if ( flag .eq. GLB_GET ) then
         dvalue = two_bond_bias_sdev_multiplier
      elseif ( flag .eq. GLB_SET ) then
         two_bond_bias_sdev_multiplier = dvalue
      else
         write(6,*) 'TWO_BOND_BIAS_SDEV_MULTIPLIER: invalid global flag'
     &        ,flag
         stop
      endif

      return
      end

      subroutine twh_two_bond_bias_update_frequency(flag,ivalue)
c     ******************************************************************
c     * update frequency for adjusting the two-bond bias values        *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ivalue
c     --- information variables
      integer two_bond_bias_update_frequency
      save two_bond_bias_update_frequency

      if ( flag .eq. GLB_GET ) then
         ivalue = two_bond_bias_update_frequency
      elseif ( flag .eq. GLB_SET ) then
         two_bond_bias_update_frequency = ivalue
      else
         write(6,*) 'TWO_BOND_BIAS_UPDATE_FREQ: invalid global flag'
     &        ,flag
         stop
      endif

      return
      end

      subroutine twh_two_bond_bias_vibrange(flag,index,dvalue)
c     ******************************************************************
c     * min and max "bond-length" for the two-bond biasing             *
c     *                                                                *
c     * originally written 10-12-2008 by M.G. Martin                   *
c     * last modified 10-12-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,index
      double precision dvalue
c     --- information variables
      double precision two_bond_bias_vibrange
      dimension two_bond_bias_vibrange(2)
      save two_bond_bias_vibrange

      if ( flag .eq. GLB_GET ) then
         dvalue = two_bond_bias_vibrange(index)
      elseif ( flag .eq. GLB_SET ) then
         two_bond_bias_vibrange(index) = dvalue
      else
         write(6,*) 'TWO_BOND_BIAS_VIBRANGE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_txangle(flag,ibond,ibend,ivalue)
c     ******************************************************************
c     * the type of extra angle for the CBMC regrowth moves            *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modified 10-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer flag,ibond,ibend,ivalue
c     --- information variables
      integer txangle
      dimension txangle(NNBOND,MAXBEND)
      save txangle

      if ( flag .eq. GLB_GET ) then
         ivalue = txangle(ibond,ibend)
      elseif ( flag .eq. GLB_SET ) then
         txangle(ibond,ibend) = ivalue
      else
         write(6,*) 'TXANGLE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_uvtfactor(flag,imolty,dvalue)
c     ******************************************************************
c     * factor used in the grand canonical ensemble for acceptance of  *
c     * insertion or deletion moves                                    *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,imolty
      double precision dvalue
c     --- information variables
      double precision uvtfactor
      dimension uvtfactor(NTMAX)
      save uvtfactor

      if ( flag .eq. GLB_GET ) then
         dvalue = uvtfactor(imolty)
      elseif ( flag .eq. GLB_SET ) then
         uvtfactor(imolty) = dvalue
      else
         write(6,*) 'UVTFACTOR: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_vecprev(flag,xvalue,yvalue,zvalue)
c     ******************************************************************
c     * vector pointing to the previous atom used during CBMC moves    *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modified 10-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      double precision xvalue,yvalue,zvalue
c     --- information variables
      double precision xvecprev,yvecprev,zvecprev
      save xvecprev,yvecprev,zvecprev

      if ( flag .eq. GLB_GET ) then
         xvalue = xvecprev
         yvalue = yvecprev
         zvalue = zvecprev
      elseif ( flag .eq. GLB_SET ) then
         xvecprev = xvalue
         yvecprev = yvalue
         zvecprev = zvalue
      else
         write(6,*) 'VECPREV: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_version(flag,cvalue)
c     ******************************************************************
c     * string describing the version number                           *
c     *                                                                *
c     * originally written 08-09-2011 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      character*(*) cvalue
c     --- storage values
      character*7 version
      save version

      if ( flag .eq. GLB_SET ) then
         version = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = version
      else
         write(6,*) 'VERSION: invalid global flag',flag
      endif

      return
      end

      subroutine twh_vequiv(flag,dvalue)
c     ******************************************************************
c     * energy used to rescale the Rosenbluth weights                  *
c     *                                                                *
c     * originally written 10-06-2008 by M.G. Martin                   *
c     * last modified 10-06-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag
      double precision dvalue
c     --- information variables
      double precision vequiv
      save vequiv

      if ( flag .eq. GLB_GET ) then
         dvalue = vequiv
      elseif ( flag .eq. GLB_SET ) then
         vequiv = dvalue
      else
         write(6,*) 'VEQUIV: invalid global flag',flag
         stop
      endif

      return
      end


      subroutine twh_vibcoeff(flag,itype,icoeff,dvalue)
c     ******************************************************************
c     * the coefficients used to compute the vibrational energies      *
c     *                                                                *
c     * originally written 03-30-2006 by M.G. Martin                   *
c     * last modified 03-30-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,icoeff
      double precision dvalue
c     --- local variables
      integer ittemp,ictemp
c     --- information variables
      double precision vibcoeff
      dimension vibcoeff(TVIBMAX,MINVIBCOEFF:MAXVIBCOEFF)
      save vibcoeff

      if ( flag .eq. GLB_SET ) then
         vibcoeff(itype,icoeff) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = vibcoeff(itype,icoeff)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            do ictemp = MINVIBCOEFF,icoeff
               vibcoeff(ittemp,ictemp) = dvalue
            enddo
         enddo
      else
         write(6,*) 'VIBCOEFF: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_vibff(flag,itype,cvalue)
c     ******************************************************************
c     * the name of the forcefield for matching vibrations             *
c     *                                                                *
c     * originally written 03-30-2006 by M.G. Martin                   *
c     * last modified 03-30-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype
      character*(FFNAMELEN) cvalue
c     --- local variables
      integer ittemp
c     --- information variables
      character*(FFNAMELEN) vibff
      dimension vibff(TVIBMAX)
      save vibff

      if ( flag .eq. GLB_SET ) then
         vibff(itype) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = vibff(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            vibff(ittemp) = cvalue
         enddo
      else
         write(6,*) 'VIBFF: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_vibnames(flag,itype,isame,iatom,cvalue)
c     ******************************************************************
c     * the names of the atoms for matching the vibrations             *
c     *                                                                *
c     * originally written 03-30-2006 by M.G. Martin                   *
c     * last modified 03-30-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,isame,iatom
      character*(FFNAMELEN) cvalue
c     --- local variables
      integer ittemp,istemp,intemp
c     --- information variables
      character*(FFNAMELEN) vibnames
      dimension vibnames(TVIBMAX,MAXSAMEV,MAXVNAME)
      save vibnames

      if ( flag .eq. GLB_SET ) then
         vibnames(itype,isame,iatom) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = vibnames(itype,isame,iatom)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            do istemp = 1,MAXSAMEV
               do intemp = 1,MAXVNAME
                  vibnames(ittemp,istemp,intemp) = cvalue
               enddo
            enddo
         enddo
      else
         write(6,*) 'VIBNAMES: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_viborder(flag,itype,cvalue)
c     ******************************************************************
c     * the bond order for each vibration type for matching in assemble*
c     *                                                                *
c     * originally written 12-09-2006 by M.G. Martin                   *
c     * last modified 12-11-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype
      character*(*) cvalue
c     --- local variables
      integer ittemp
c     --- information variables
      character*(FFNAMELEN) viborder
      dimension viborder(TVIBMAX)
      save viborder

      if ( flag .eq. GLB_SET ) then
         viborder(itype) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = viborder(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do ittemp = 1,itype
            viborder(ittemp) = cvalue
         enddo
      else
         write(6,*) 'VIBORDER: invalid global flag',flag
         stop
      endif

      return
      end


      subroutine twh_vibrang(flag,index,dvalue)
c     ******************************************************************
c     * range for bond lengths generated by some CBMC options          *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,index
      double precision dvalue
c     --- local variables
      integer itemp
c     --- information variables
      double precision vibrang
      dimension vibrang(2)
      save vibrang

      if ( flag .eq. GLB_SET ) then
         vibrang(index) = dvalue
      elseif ( flag .eq. GLB_GET ) then
         dvalue = vibrang(index)
      elseif ( flag .eq. GLB_INIT ) then
         do itemp = 1,index
            vibrang(itemp) = dvalue
         enddo
      else
         write(6,*) 'VIBRANG: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_vibration_order(flag,iatom,ivib,cvalue)
c     ******************************************************************
c     * the bond order for each vibration in the system for the current*
c     * molecule that is being assembled                               *
c     *                                                                *
c     * originally written 12-10-2006 by M.G. Martin                   *
c     * last modified 12-11-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,iatom,ivib
      character*(*) cvalue
c     --- local variables
      integer iatemp,ivtemp
c     --- information variables
      character*(FFNAMELEN) vibration_order
      dimension vibration_order(NUMAX,NNBOND)
      save vibration_order

      if ( flag .eq. GLB_SET ) then
         vibration_order(iatom,ivib) = cvalue
      elseif ( flag .eq. GLB_GET ) then
         cvalue = vibration_order(iatom,ivib)
      elseif ( flag .eq. GLB_INIT ) then
         do iatemp = 1,iatom
            do ivtemp = 1,NNBOND
               vibration_order(iatemp,ivtemp) = cvalue
            enddo
         enddo
      else
         write(6,*) 'VIBRATION_ORDER: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_vibstyle(flag,itype,ivalue)
c     ******************************************************************
c     * bond vibration style for each type of vibration                *
c     *                                                                *
c     * originally written 03-23-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,itype,ivalue
c     --- local variables
      integer itemp
c     --- information variables
      integer vibstyle
      dimension vibstyle(TVIBMAX)
      save vibstyle

      if ( flag .eq. GLB_SET ) then
         vibstyle(itype) = ivalue
      elseif ( flag .eq. GLB_GET ) then
         ivalue = vibstyle(itype)
      elseif ( flag .eq. GLB_INIT ) then
         do itemp = 1,itype
            vibstyle(itemp) = ivalue
         enddo
      else
         write(6,*) 'VIBSTYLE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_vname(flag,iener,cvalue)
c     ******************************************************************
c     * name of the various energy terms                               *
c     *                                                                *
c     * originally written 10-10-2008 by M.G. Martin                   *
c     * last modified 10-10-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,iener
      character*(*) cvalue
c     --- information variables
      character*15 vname
      dimension vname(MAXENER)
      save vname

      if ( flag .eq. GLB_GET ) then
         cvalue = vname(iener)
      elseif ( flag .eq. GLB_SET ) then
         vname(iener) = cvalue
      else
         write(6,*) 'VNAME: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_volmaxdispfreq(flag,ivalue)
c     ******************************************************************
c     * frequency for updating the maximum volume displacement freqs   *
c     *                                                                *
c     * originally written 12-09-2007 by M.G. Martin                   *
c     * last modified 12-09-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- information variables
      integer volmaxdispfreq
      save volmaxdispfreq

      if ( flag .eq. GLB_GET ) then
         ivalue = volmaxdispfreq
      elseif ( flag .eq. GLB_SET ) then
         volmaxdispfreq = ivalue
      else
         write(6,*) 'VOLMAXDISPFREQ: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_weight_freq(flag,ivalue)
c     ******************************************************************
c     * TMMC frequency to update biasing function                      *
c     *                                                                *
c     * originally written 10-13-2008 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,ivalue
c     --- information variables
      integer weight_freq
      save weight_freq

      if ( flag .eq. GLB_GET ) then
         ivalue = weight_freq
      elseif ( flag .eq. GLB_SET ) then
         weight_freq = ivalue
      else
         write(6,*) 'WEIGHT_FREQ: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_xangle(flag,idim,ibond,ibend,ivalue)
c     ******************************************************************
c     * extra angle information for CBMC growth moves                  *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modified 10-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,idim,ibond,ibend,ivalue
c     --- information variables
      integer xangle
      dimension xangle(3,NNBOND,MAXBEND)
      save xangle

      if ( flag .eq. GLB_GET ) then
         ivalue = xangle(idim,ibond,ibend)
      elseif ( flag .eq. GLB_SET ) then
         xangle(idim,ibond,ibend) = ivalue
      else
         write(6,*) 'XANGLE: invalid global flag',flag
         stop
      endif

      return
      end

      subroutine twh_xbias(flag,idim,ibond,itor,ivalue)
c     ******************************************************************
c     * extra bias information for CBMC growth moves                   *
c     *                                                                *
c     * originally written 10-05-2008 by M.G. Martin                   *
c     * last modified 10-05-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer flag,idim,ibond,itor,ivalue
c     --- information variables
      integer xbias
      dimension xbias(2,NNBOND,MAXTOR)
      save xbias

      if ( flag .eq. GLB_GET ) then
         ivalue = xbias(idim,ibond,itor)
      elseif ( flag .eq. GLB_SET ) then
         xbias(idim,ibond,itor) = ivalue
      else
         write(6,*) 'XBIAS: invalid global flag',flag
         stop
      endif

      return
      end

