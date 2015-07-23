#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2006-2012 Marcus G. Martin                       *
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
      subroutine twh_eng_external_field(lfailure,iflag,activebox
     &     ,nchainmove,chaininfo,loverlap,v_ext)
c     ******************************************************************
c     * calculates the external field energy for all of the molecules  *
c     * listed in chaininfo                                            *
c     *                                                                *
c     * Arguments set by subroutine:                                   *
c     *  v_ext                                                         *
c     *  lovrlap                                                       *
c     *                                                                *
c     * Error conditions:                                              *
c     *   returns early if loverlap or lfailure                        *
c     *                                                                *
c     * rewritten from engmolec 08-04-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_VFIELD
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer iflag,activebox,nchainmove
      integer chaininfo,nunit
      dimension chaininfo(nchainmove,3)
      logical loverlap
      double precision v_ext
c     --- local variables
      integer iunit,ntii,imolty,ichainmove,testcrd,testchain
      double precision cord
      dimension cord(3)

      v_ext = 0.0d0
c     --- no external fields in the ideal box
      if ( activebox .eq. 0 ) return

      do ichainmove = 1,nchainmove
         testcrd = chaininfo(ichainmove,iflag)
         if ( testcrd .ne. CRD_NONE ) then
            testchain = chaininfo(ichainmove,3)
            call twh_moltyp(GLB_GET,testchain,imolty)
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
               call twh_coordinates(GLB_GET,testcrd,testchain,iunit
     &              ,cord(1),cord(2),cord(3))
               call twh_ntype(GLB_GET,imolty,iunit,ntii)
               v_ext = v_ext + twh_vfield(lfailure,loverlap,activebox
     &              ,ntii,imolty,testchain,iunit,cord)
               if ( lfailure .or. loverlap ) return
            enddo
         endif
      enddo

      return
      end

      subroutine twh_eng_solvation(lfailure,iflag,activebox,nchainmove
     &     ,chaininfo,v_solv)
c     ******************************************************************
c     * calculates the solvation energy for the box with the changes   *
c     * made for the molecules in chaininfo                            *
c     *                                                                *
c     * Arguments set by subroutine:                                   *
c     *  v_solv                                                        *
c     *                                                                *
c     * Error conditions:                                              *
c     *   lfailure                                                     *
c     *                                                                *
c     * rewritten from engmolec 08-04-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_VSASA
#define FUNCTION_VEEFONE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer iflag,activebox,nchainmove
      integer chaininfo
      dimension chaininfo(nchainmove,3)
      double precision v_solv
c     --- local variables
      logical lcompute
      integer ichainmove,testcrd,testchain,isolvtype
      double precision vnt,vsolvold

      v_solv = 0.0d0
      call twh_isolvtype(GLB_GET,isolvtype)

c     --- no solvation energy in the ideal box
      if ( activebox .eq. 0 ) return

      if (isolvtype .eq. SOLV_NONE ) then
c        --- no solvation energy
         return
      elseif (isolvtype .eq. SOLV_SASA ) then
c           --- SASA implicit water model
         do ichainmove = 1,nchainmove
            testcrd = chaininfo(ichainmove,iflag)
            if ( testcrd .ne. CRD_NONE ) then
               testchain = chaininfo(ichainmove,3)
               v_solv = v_solv + twh_vsasa(testchain,activebox
     &              ,testcrd)
            endif
         enddo
      elseif (isolvtype .eq. SOLV_EEF1 ) then
c        --- EEF1 implicit water model
         do ichainmove = 1,nchainmove
            testcrd = chaininfo(ichainmove,iflag)
            if ( testcrd .ne. CRD_NONE ) then
               testchain = chaininfo(ichainmove,3)
               v_solv = v_solv + twh_veefone(testchain,activebox
     &              ,testcrd)
            endif
         enddo
      elseif (isolvtype .eq. SOLV_TRAMONTO ) then
c        --- Tramonto solvation free energy
c        --- determine if any molecules have moved from the previous
c        --- energy calculation
         lcompute = .false.
         ichainmove = 0
         do while ( ichainmove .lt. nchainmove .and. (.not. lcompute) )
            ichainmove = ichainmove + 1
            if ( chaininfo(ichainmove,iflag) .ne. CRD_REAL ) then
               lcompute = .true.
            endif
         enddo
         if ( lcompute ) then
            call twh_fielddft(lfailure,iflag,activebox,nchainmove
     &           ,chaininfo,vnt)
            if ( lfailure ) return
            v_solv = v_solv + vnt
         else
c           --- use the old solvation energy
            call twh_boxvsolvation(GLB_GET,activebox,vsolvold)
            v_solv = v_solv + vsolvold
         endif
      endif   

      return
      end


      subroutine twh_eng_two_body(lfailure,iflag,activebox,nchainmove
     &     ,chaininfo,lpartial,lgetcom,loverlap
     &     ,v_inter,v_elect,l_vinter,l_velect,v_foreign)
c     ******************************************************************
c     * compute the two-body energy in the active box for all molecules*
c     * interacting with those in the chaininfo                        *
c     *                                                                *
c     * rewritten from engmolec 08-05-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_GT
#define FUNCTION_VCOULOMB
#define FUNCTION_VTWOBODY
#include "functions.h"
c     --- debug toggle
#define DEBUG_ENG_TWO_BODY 0
c     --- variables passed to/from the subroutine
      logical lfailure
      integer iflag,activebox,nchainmove
      integer chaininfo,numiunit,numjunit
      dimension chaininfo(nchainmove,3)
      logical lpartial,lgetcom,loverlap
      double precision v_inter,v_elect,l_vinter,l_velect
      double precision v_foreign
      dimension v_foreign(MAX_FOREIGN_LAMBDA)
c     --- local variables
      logical lvalid,lcompute,leval_v,leval_c,lusecom,lcoulomb
      logical luse_ci,luse_cj,luse_vi,luse_vj
c     --- integer scalars
      integer cutstyle,testcrd,testchain,ichainmove,imolty,jmolty,jlist
      integer jchain,itest,jj,ntjj,iunit,ntii,lambda,boxreal,nmolty,ncmt
      integer scp_num_foreign_lambda,testmolty,liststart,listend,jcrd
c     --- double precision scalars
      double precision rxi,ryi,rzi,rxj,ryj,rzj,rxij,ryij,rzij,rij
      double precision xcmi,ycmi,zcmi,rcmi,rijsq,rcm,rcmsq,dvdllj,dvdlc
      double precision venergy,xcmj,ycmj,zcmj,rcmj,rcutinsq
      double precision rcmreal,qqfact,rcelect
c     --- double precision arrays
      double precision vf
      dimension vf(MAX_FOREIGN_LAMBDA)

c     --- no two-body intermolecular energy in box 0
      if ( activebox .eq. 0 ) return

c     --- retrieve constant
      call twh_constant_qqfact(GLB_GET,qqfact)
      call twh_nmolty(GLB_GET,nmolty)
      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
      call twh_lusecom(GLB_GET,lusecom)
      call twh_rcelect(GLB_GET,activebox,rcelect)
      call twh_lcoulomb(GLB_GET,lcoulomb)
      call twh_rcutinsq(GLB_GET,rcutinsq)

#if DEBUG_ENG_TWO_BODY
      write(6,*) 'nchainmove:',nchainmove
#endif

c     --- this subroutine uses cutstyle 0
      cutstyle = 0

      do ichainmove = 1,nchainmove
         testcrd = chaininfo(ichainmove,iflag)
         if ( testcrd .ne. CRD_NONE ) then
c           --- calculate the center of mass and give it a dummy
c           --- chain identifier
            testchain = chaininfo(ichainmove,3)
#if DEBUG_ENG_TWO_BODY
            write(6,*) 'ichainmove:',ichainmove
            write(6,*) 'testcrd:',testcrd
            write(6,*) 'testchain:',testchain
#endif
            call twh_moltyp(GLB_GET,testchain,imolty)
            call twh_nunit(GLB_GET,imolty,numiunit)
            if ( lusecom ) then
               if ( lgetcom ) then
c                 --- store the real coordinate center of mass radius
                  call twh_rcmu(GLB_GET,testchain,rcmreal)
                  call twh_nboxi(GLB_GET,testchain,boxreal)
c                 --- compute the center of mass of this molecule in
c                 --- the activebox
                  call twh_nboxi(GLB_SET,testchain,activebox)
                  call twh_ctrmas(lfailure,testcrd,activebox,testchain
     &                 ,CTR_ENERGY)
                  if ( lfailure ) return
c                 --- retrieve the center of mass of the test molecule
                  call twh_com(GLB_GET,testcrd,testchain,xcmi,ycmi,zcmi)
                  call twh_rcmu(GLB_GET,testchain,rcmi)
c                 --- reset the real value of the radius and box
                  call twh_rcmu(GLB_SET,testchain,rcmreal)
                  call twh_nboxi(GLB_SET,testchain,boxreal)
               else
                  call twh_com(GLB_GET,testcrd,testchain
     &                 ,xcmi,ycmi,zcmi)
                  call twh_rcmu(GLB_GET,testchain,rcmi)
               endif
            endif

c           --- loop over all other chains in target box
            do testmolty = 0,nmolty
               if ( testmolty .eq. 0 ) then
c                 --- compute against higer index chainmove molecules
                  liststart = ichainmove + 1
                  listend = nchainmove
               else
c                 --- compute against the list of molecules in this box
                  jmolty = testmolty
                  liststart = 1
                  call twh_ncmt(GLB_GET,activebox,jmolty,ncmt)
                  listend = ncmt
               endif
               do jlist = liststart,listend
                  if ( testmolty .eq. 0 ) then
c                    --- check the chaininfo
                     jcrd = chaininfo(jlist,iflag)
                     if ( jcrd .eq. CRD_NONE ) then
                        lcompute = .false.
                     else
                        lcompute = .true.
                        jchain = chaininfo(jlist,3)
                        call twh_moltyp(GLB_GET,jchain,jmolty)
                     endif
                  else
c                    --- check the chainlist
                     call twh_chainlist(GLB_GET,jlist,activebox,jmolty
     &                    ,jchain)
c                    --- use real coordinates for this molecule
                     jcrd = CRD_REAL
c                    --- make sure this is not in the chaininfo list
                     itest = 0
                     lvalid = .true.
                     do while ( itest .lt. nchainmove .and. lvalid )
                        itest = itest + 1
                        if ( jchain .eq. chaininfo(itest,3) ) then
                           lvalid = .false.
                        endif
                     enddo
                     lcompute = .false.
                     if ( lvalid ) then
                        lcompute = .true.
                        if ( lusecom ) then
c                          --- use additional com cut-off
c                           --- minimum image the ctrmas pair distance
                           call twh_com(GLB_GET,jcrd,jchain
     &                          ,xcmj,ycmj,zcmj)
                           call twh_mimage (activebox,xcmi,ycmi,zcmi
     &                          ,xcmj,ycmj,zcmj,rxij,ryij,rzij)
                           rijsq = rxij*rxij + ryij*ryij + rzij*rzij
                           call twh_rcmu(GLB_GET,jchain,rcmj)
c                          --- start with the vdw cutoff
                           call twh_rcut(GLB_GET,rcm)
                           if ( lcoulomb ) then
c                             --- use the coulombic cutoff if larger
                              if ( rcelect .gt. rcm ) rcm = rcelect
                           endif
c                          --- add in the com radius components
                           rcm = rcm + rcmi + rcmj
                           rcmsq = rcm*rcm
                           if ( twh_cmp_gt(rijsq, rcmsq) ) then
                              lcompute = .false.
                           endif
                        endif
                     endif
                  endif

                  if ( lcompute ) then
c                    --- loop over all beads jj of chain jchain 
                     call twh_nunit(GLB_GET,jmolty,numjunit)
                     do jj = 1, numjunit
                        call twh_ntype(GLB_GET,jmolty,jj,ntjj)
                        call twh_coordinates(GLB_GET,jcrd,jchain
     &                       ,jj,rxj,ryj,rzj)
                        call twh_luse_c(GLB_GET,jmolty,jj,luse_cj)
                        call twh_luse_v(GLB_GET,jmolty,jj,luse_vj)
c                       --- loop over all beads iunit of testchain
                        do iunit = 1,numiunit
                           call twh_ntype(GLB_GET,imolty,iunit,ntii)
                           call twh_luse_v(GLB_GET,imolty,iunit,luse_vi)
c                          --- set evaulation logicals
                           leval_v = luse_vi .and. luse_vj
                           call twh_luse_c(GLB_GET,imolty,iunit,luse_ci)
                           leval_c = luse_ci .and. luse_cj
                           if (leval_v .or. leval_c) then
                              call twh_coordinates(GLB_GET,testcrd
     &                             ,testchain,iunit,rxi,ryi,rzi)
c                             --- minimum image the pair separation
                              call twh_mimage (activebox,rxi,ryi,rzi
     &                             ,rxj,ryj,rzj,rxij,ryij,rzij )
                              rijsq = rxij*rxij + ryij*ryij 
     &                             + rzij*rzij
c                             --- if partial correction then we do
c                             --- not compute for rijsq < rcutinsq
                              if ( .not. (lpartial 
     &                             .and. rijsq .lt. rcutinsq) ) then
                                 if (leval_v) then
                                    v_inter = v_inter
     &                                   + twh_vtwobody(lfailure
     &                                   ,rijsq,ntii,ntjj,cutstyle
     &                                   ,loverlap,vf,dvdllj)
                                    if ( lfailure ) return
                                    if ( loverlap ) return
                                    l_vinter = l_vinter + dvdllj
                                    do lambda = 1
     &                                   ,scp_num_foreign_lambda
                                       v_foreign(lambda)
     &                                      = v_foreign(lambda)
     &                                      + vf(lambda)
                                    enddo
                                 endif
c                                --- electrostatics
                                 if ( lcoulomb .and. leval_c) then
                                    rij = dsqrt(rijsq)
                                    venergy = twh_vcoulomb(activebox
     &                                   ,rij,imolty,iunit,jmolty,jj
     &                                   ,.false., vf, dvdlc
     &                                   ,loverlap)
                                    if ( loverlap ) return
c                                   --- vcoulomb returns units [EE]
c                                   --- convert foreign energy to
c                                   --- units [K] here;
c                                   --- v_elect converted to [K]
c                                   --- later.
                                    v_elect = v_elect + venergy
c                                   --- derivatives in units of [K]
                                    l_velect = l_velect+qqfact*dvdlc
                                    do lambda = 1
     &                                   ,scp_num_foreign_lambda
                                       v_foreign(lambda) 
     &                                      = v_foreign(lambda)
     &                                      +qqfact*vf(lambda)
                                    enddo
                                 endif
                              endif
                           endif
                        enddo
                     enddo
                  endif
               enddo
            enddo
         endif
      enddo

      return
      end

      subroutine twh_eng_angle_angle(ichain,ibox,icf,v_ben)
c     ******************************************************************
c     * calculates the angle-angle energy                              *
c     *                                                                *
c     * Arguments set by subroutine:                                   *
c     *  v_ben                                                         *
c     *                                                                *
c     * Error conditions:                                              *
c     *   none                                                         *
c     *                                                                *
c     * rewritten from engmolec 08-25-2006 by M.G. Martin              *
c     * last modified 08-30-2012 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJAA0
#define FUNCTION_GET_IJAA1
#define FUNCTION_GET_IJAA2
#define FUNCTION_GET_INAA
#define FUNCTION_VANGANG
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ichain,icf,ibox
      double precision v_ben
c     --- local variables
c     --- integer scalars
      integer iunit,itor,ibead,imolty,nunit
c     --- integer arrays
      integer ttype
      dimension ttype(4)
c     --- double precision arrays
      double precision xtmp,ytmp,ztmp
      dimension xtmp(4),ytmp(4),ztmp(4)
      double precision dumangle
      dimension dumangle(2)

      if ( icf .eq. CRD_NONE ) return
      call twh_moltyp(GLB_GET,ichain,imolty)
      call twh_nunit(GLB_GET,imolty,nunit)
      do iunit = 1,nunit
         ttype(1) = iunit
         do itor = 1,twh_get_inaa(imolty,iunit)
            ttype(2) = twh_get_ijaa0(imolty,iunit,itor)
            ttype(3) = twh_get_ijaa1(imolty,iunit,itor)
            ttype(4) = twh_get_ijaa2(imolty,iunit,itor)
            do ibead = 1,4
               call twh_coordinates(GLB_GET,icf,ichain,ttype(ibead)
     &              ,xtmp(ibead),ytmp(ibead),ztmp(ibead))
            enddo
c           --- compute the angle-angle energy
            v_ben = v_ben + 
     &            twh_vangang(AA_COORD,ibox,imolty,iunit,itor,dumangle
     &            ,xtmp,ytmp,ztmp)
         enddo
      enddo                   

      return
      end

      subroutine twh_eng_bending(lfailure,ichain,ibox,icf,lovrlap,v_ben)
c     ******************************************************************
c     * calculates the bending energy                                  *
c     *                                                                *
c     * Arguments set by subroutine:                                   *
c     *  lovrlap                                                       *
c     *  v_ben                                                         *
c     *                                                                *
c     * Error conditions:                                              *
c     *   returns early if lovrlap is true                             *
c     *                                                                *
c     * rewritten from engmolec 08-25-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJBEN2
#define FUNCTION_GET_IJBEN3
#define FUNCTION_GET_INBEN
#define FUNCTION_VANGLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ichain,icf,ibox
      logical lfailure,lovrlap
      double precision v_ben
c     --- local variables
c     --- integer scalars
      integer iunit,iben,ip1,ip2,imolty,nunit
c     --- double precision scalars
      double precision duma,dumb,dumc
c     --- double precision arrays
      double precision xtmp,ytmp,ztmp
      dimension xtmp(3),ytmp(3),ztmp(3)

      if ( icf .eq. 0 ) return
      call twh_moltyp(GLB_GET,ichain,imolty)
      call twh_nunit(GLB_GET,imolty,nunit)
      do iunit = 2,nunit
         call twh_coordinates(GLB_GET,icf,ichain,iunit
     &        ,xtmp(1),ytmp(1),ztmp(1))
         do iben = 1, twh_get_inben(imolty,iunit)
            ip2 = twh_get_ijben3(imolty,iunit,iben)
            if ( ip2 .lt. iunit ) then
               ip1 = twh_get_ijben2(imolty,iunit,iben)
c              --- fetch the coordinates
               call twh_coordinates(GLB_GET,icf,ichain,ip1
     &              ,xtmp(2),ytmp(2),ztmp(2))
               call twh_coordinates(GLB_GET,icf,ichain,ip2
     &              ,xtmp(3),ytmp(3),ztmp(3))
               v_ben = v_ben + 
     &               twh_vangle(lfailure,lovrlap,ibox,imolty,iunit
     &               ,iben,1,duma,dumb,dumc,xtmp,ytmp,ztmp)
               if ( lfailure ) return
               if ( lovrlap ) return
            endif
         enddo
      enddo

      return
      end

      subroutine twh_eng_bond_bond(ichain,ibox,icf,v_vib)
c     ******************************************************************
c     * calculates the bond bond energy                                *
c     * Arguments set by subroutine:                                   *
c     *  v_vib                                                         *
c     *                                                                *
c     * Error conditions:                                              *
c     *   none                                                         *
c     *                                                                *
c     * rewritten from engmolec 08-25-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_GET_IJBNBN1
#define FUNCTION_GET_IJBNBN2
#define FUNCTION_GET_INBNBN
#define FUNCTION_VBONBON
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ichain,icf,ibox
      double precision v_vib
c     --- local variables
c     --- integer scalars
      integer iunit,iben,ip1,ip2,imolty,nunit
c     --- double precision scalars
      double precision rxuij,ryuij,rzuij
      double precision rxunit,ryunit,rzunit,rxipone,ryipone,rzipone
      double precision rxiptwo,ryiptwo,rziptwo
c     --- double precision arrays
      double precision length
      dimension length(2)

      if ( icf .eq. CRD_NONE ) return
      call twh_moltyp(GLB_GET,ichain,imolty)
      call twh_nunit(GLB_GET,imolty,nunit)
      do iunit = 1,nunit
         do iben = 1,twh_get_inbnbn(imolty,iunit)
            ip1 = twh_get_ijbnbn1(imolty,iunit,iben)
            ip2 = twh_get_ijbnbn2(imolty,iunit,iben)
            call twh_coordinates(GLB_GET,icf,ichain,iunit
     &           ,rxunit,ryunit,rzunit)
            call twh_coordinates(GLB_GET,icf,ichain,ip1
     &           ,rxipone,ryipone,rzipone)
            call twh_coordinates(GLB_GET,icf,ichain,ip2
     &           ,rxiptwo,ryiptwo,rziptwo)
c           --- compute and store first bond length
            call twh_mimage(ibox,rxunit,ryunit,rzunit
     &           ,rxipone,ryipone,rzipone,rxuij,ryuij,rzuij)
            length(1) = twh_distance(rxuij,ryuij,rzuij)
c           --- compute and store second bond length
            call twh_mimage(ibox,rxunit,ryunit,rzunit
     &           ,rxiptwo,ryiptwo,rziptwo,rxuij,ryuij,rzuij)
            length(2) = twh_distance(rxuij,ryuij,rzuij)
c           --- compute energy as a bond vibration term
            v_vib = v_vib + twh_vbonbon(imolty,iunit,iben,length)
            if ( (iunit .eq. 456) .or. (ip1 .eq. 456) .or. 
     &           (ip2 .eq. 456) ) then
            endif
         enddo
      enddo

      return
      end

      subroutine twh_eng_bond_vibrations(lfailure,ichain,ibox,icf
     &     ,lovrlap,v_vib )
c     ******************************************************************
c     * calculates the vibration bond energy                           *
c     *                                                                *
c     * Arguments set by subroutine:                                   *
c     *  lovrlap                                                       *
c     *  v_vib                                                         *
c     *                                                                *
c     * Error conditions:                                              *
c     *   returns early if lovrlap is true                             *
c     *                                                                *
c     * rewritten from engmolec 08-25-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_GET_INVIB
#define FUNCTION_VBOND
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer icf,ichain,ibox
      logical lfailure,lovrlap
      double precision v_vib

c     --- local variables
c     --- integer scalars
      integer iunit,ivib,ip1,it,imolty,nunit
c     --- double precision scalars
      double precision rxuij,ryuij,rzuij,dist
      double precision rxunit,ryunit,rzunit,rxipone,ryipone,rzipone

      if ( icf .eq. CRD_NONE ) return
      call twh_moltyp(GLB_GET,ichain,imolty)
      call twh_nunit(GLB_GET,imolty,nunit)
      do iunit = 2,nunit
         do ivib = 1, twh_get_invib(imolty,iunit)
            call twh_ijvib(GLB_GET,imolty,iunit,ivib,ip1)
            call twh_itvib(GLB_GET,imolty,iunit,ivib,it)
c           --- increment vibrational energy
            if ( ip1 .lt. iunit ) then
c              --- compute distance
c              --- minimum image the bond vectors. this is 
c              --- required in case the molecule network wraps
c              --- completely through the periodic box
               call twh_coordinates(GLB_GET,icf,ichain,iunit
     &              ,rxunit,ryunit,rzunit)
               call twh_coordinates(GLB_GET,icf,ichain,ip1
     &              ,rxipone,ryipone,rzipone)
               call twh_mimage (ibox,rxunit,ryunit,rzunit
     &              ,rxipone,ryipone,rzipone,rxuij,ryuij,rzuij)
               dist = twh_distance(rxuij,ryuij,rzuij)
               v_vib = v_vib + twh_vbond(lfailure,it,dist,imolty,iunit
     &              ,ip1,lovrlap)
               if ( lfailure ) return
               if ( lovrlap ) return
            endif
         enddo
      enddo

      return
      end

      subroutine twh_eng_improper_torsions(ichain,ibox,icf,v_tor)
c     ******************************************************************
c     * calculates the improper torsions energy                        *
c     *                                                                *
c     * Arguments set by subroutine:                                   *
c     *  v_tor                                                         *
c     *                                                                *
c     * Error conditions:                                              *
c     *   none                                                         *
c     *                                                                *
c     * rewritten from engmolec 08-25-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJIMPROP2
#define FUNCTION_GET_IJIMPROP3
#define FUNCTION_GET_IJIMPROP4
#define FUNCTION_GET_INIMPROP
#define FUNCTION_VIMPROPER
#include "functions.h"
      integer ichain,icf,ibox
      double precision v_tor
c     --- local variables
c     --- logical scalars
      logical lbadstereo
c     --- integer scalars
      integer iunit,itor,ibead,imolty,nunit
c     --- integer arrays
      integer ttype
      dimension ttype(4)
c     --- double precision scalars
      double precision duma
c     --- double precision arrays
      double precision xtmp,ytmp,ztmp
      dimension xtmp(4),ytmp(4),ztmp(4)

      if ( icf .eq. CRD_NONE ) return
      call twh_moltyp(GLB_GET,ichain,imolty)
      call twh_nunit(GLB_GET,imolty,nunit)
c     --- these are only stored on a single atom as they are 
c     --- not symmetrically defined
      do iunit = 1,nunit
         ttype(1) = iunit
         do itor = 1,twh_get_inimprop(imolty,iunit)
            ttype(2) = twh_get_ijimprop2(imolty,iunit,itor)
            ttype(3) = twh_get_ijimprop3(imolty,iunit,itor)
            ttype(4) = twh_get_ijimprop4(imolty,iunit,itor)
c           --- fetch the coordinates
            do ibead = 1,4
               call twh_coordinates(GLB_GET,icf,ichain,ttype(ibead)
     &              ,xtmp(ibead),ytmp(ibead),ztmp(ibead))
            enddo
            v_tor = v_tor + 
     &            twh_vimproper(ibox,imolty,iunit,itor
     &            ,xtmp,ytmp,ztmp,duma,lbadstereo)
         enddo
      enddo

      return
      end

      subroutine twh_eng_nonbonded_intra(lfailure,ichain,ibox,icf
     &     ,lovrlap,v_intra,v_elect )
c     ******************************************************************
c     * calculates the nonbonded intramolecular energy                 *
c     *                                                                *
c     * Arguments set by subroutine:                                   *
c     *  lovrlap                                                       *
c     *  v_intra                                                       *
c     *  v_elect                                                       *
c     *                                                                *
c     * Error conditions:                                              *
c     *   returns early if lovrlap is true                             *
c     *                                                                *
c     * rewritten from engmolec 08-25-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DERFC
#define FUNCTION_EWALD_SELF
#define FUNCTION_GET_QQATOM
#define FUNCTION_LINCLUDE
#define FUNCTION_ONEFIVETYPE
#define FUNCTION_VCOULOMB
#define FUNCTION_VINTRA
#define FUNCTION_VONEFIVE
#define FUNCTION_VTWOBODY
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer icf,ibox,ichain
      logical lfailure,lovrlap
      double precision v_intra, v_elect
c     --- local variables
c     --- logical scalars
      logical lcoulomb,lewald
c     --- integer scalars
      integer iunit,ntii,jj,ntjj,ioftype,cutstyle,imolty,nunit
c     --- double precision scalars
      double precision rxui,rzui,ryui,rxuij,ryuij,rzuij,rijsq,rij,dv
      double precision rxuj,ryuj,rzuj,calp,rcelect
c     --- double precision arrays
      double precision vf
      dimension vf(MAX_FOREIGN_LAMBDA)

c     --- global variables accessed (listed for convenience)
c     dvpm

c     --- retrieve constants
      call twh_rcelect(GLB_GET,ibox,rcelect)
      call twh_lcoulomb(GLB_GET,lcoulomb)
      call twh_lewald(GLB_GET,lewald)

      if ( icf .eq. CRD_NONE ) return
c     --- this subroutine uses cutstyle 0
      cutstyle = 0
c     --- get molecule type
      call twh_moltyp(GLB_GET,ichain,imolty)
      call twh_nunit(GLB_GET,imolty,nunit)
      do iunit = 1,nunit
         call twh_ntype(GLB_GET,imolty,iunit,ntii)
         call twh_coordinates(GLB_GET,icf,ichain,iunit,rxui,ryui,rzui)
         do jj = 1,iunit-1
            call twh_ntype(GLB_GET,imolty,jj,ntjj)
            call twh_coordinates(GLB_GET,icf,ichain,jj
     &           ,rxuj,ryuj,rzuj)
            call twh_mimage(ibox,rxui,ryui,rzui
     &           ,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij)
            rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
            if ( twh_linclude(imolty,iunit,jj) ) then
c              --- get one-five type
               ioftype = twh_onefivetype(imolty,iunit,jj)
               if ( ioftype .eq. 0 ) then
c                 --- normal vtwobody
                  v_intra = v_intra + twh_vtwobody(lfailure,rijsq
     &                 ,ntii,ntjj,cutstyle,lovrlap,vf,dv)
                  if ( lfailure ) return
               else
c                 --- special one-five
                  v_intra = v_intra + twh_vonefive(lfailure,rijsq
     &                 ,ioftype,cutstyle,lovrlap)
                  if ( lfailure ) return
               endif
               if ( lovrlap ) return

               if ( lcoulomb ) then
c                 --- electrostatics
                  rij = dsqrt(rijsq)
                  v_elect = v_elect + twh_vcoulomb(ibox,rij,imolty
     &                 ,iunit,imolty,jj,.false.,vf,dv,lovrlap)
                  if ( lovrlap ) return
               endif
            elseif ( lewald .and. ibox .ne. 0 ) then
c              --- compute the ewald intramolecular (self and
c              --- correction) terms for the interactions of 
c              --- the placed atoms with themselves, and with 
c              --- the rest of their own molecule.
               rij = dsqrt(rijsq)
               if ( rij .lt. rcelect ) then
                  call twh_ewald_calp(GLB_GET,ibox,calp)
                  v_elect = v_elect + twh_get_qqatom(imolty,iunit)
     &                 *twh_get_qqatom(imolty,jj)*
     &                 (twh_derfc(calp*rij)-1.0d0)/rij
               endif
            endif
         enddo
         if ( lewald .and. ibox .ne. 0 ) then
            v_elect = v_elect + twh_ewald_self(ibox,imolty,iunit)
         endif
      enddo

      return
      end

      subroutine twh_eng_regular_torsions(lfailure,ichain,ibox,icf
     &     ,lovrlap,v_tor)
c     ******************************************************************
c     * calculates the regular torsions energy                         *
c     * Arguments set by subroutine:                                   *
c     *  lovrlap                                                       *
c     *  v_tor                                                         *
c     *                                                                *
c     * Error conditions:                                              *
c     *   returns early if lovrlap is true                             *
c     *                                                                *
c     * rewritten from engmolec 10-09-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJTOR2
#define FUNCTION_GET_IJTOR3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_INTOR
#define FUNCTION_VTORSION
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ichain,icf,ibox
      logical lfailure,lovrlap
      double precision v_tor
c     --- local variables
c     --- integer scalars
      integer iunit,itor,ibead,imolty,nunit
c     --- local double precision
      double precision duma
c     --- integer arrays
      integer ttype
      dimension ttype(4)
c     --- double precision arrays
      double precision xtmp,ytmp,ztmp
      dimension xtmp(4),ytmp(4),ztmp(4)

      if ( icf .eq. CRD_NONE ) return
      call twh_moltyp(GLB_GET,ichain,imolty)
      call twh_nunit(GLB_GET,imolty,nunit)
      do iunit = 2,nunit
c        --- store the charge on atom j
         ttype(1) = iunit
         do itor = 1, twh_get_intor(imolty,iunit)
            ttype(4) = twh_get_ijtor4(imolty,iunit,itor)
            if ( ttype(4) .lt. iunit ) then
c                    --- store unit numbers
               ttype(2) = twh_get_ijtor2(imolty,iunit,itor)
               ttype(3) = twh_get_ijtor3(imolty,iunit,itor)
c                    --- store atom positions
               do ibead = 1,4
                  call twh_coordinates(GLB_GET,icf,ichain,ttype(ibead)
     &                 ,xtmp(ibead),ytmp(ibead),ztmp(ibead))
               enddo
c                    --- add in the torsional energy
               v_tor = v_tor + 
     &               twh_vtorsion(lfailure,ibox,1,imolty,iunit,itor
     &               ,xtmp,ytmp,ztmp,duma,lovrlap)
               if ( lfailure ) return
               if ( lovrlap ) return
            endif
         enddo
      enddo

      return
      end


      subroutine twh_eng_stillinger_weber(lfailure,activebox,testcord
     &     ,nchainmove,chaininfo,vthree)
c     ******************************************************************
c     * computes the three-body stillinger-weber energy of the         *
c     * chaininfo atoms with each other and with the other existing    *
c     * atoms in the activebox                                         *
c     *                                                                *
c     * originally written 10-10-2006 by M.G. Martin                   *
c     * last modified 10-14-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_VTHREEBODY
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer activebox,testcord,nchainmove
      integer chaininfo
      dimension chaininfo(nchainmove,3)
      double precision vthree
c     --- local variables
      logical lzero,loverlap,lcompute,kcompute,kzero,jcompute
      integer imove,ichain,icord,ntii,ntjj,imolty,jmolty,jlist
      integer istart,iend,jcord,jchain,testmolty,testmove,testchain
      integer kmolty,kstart,kend,klist,kcord,kchain,ntkk,nmolty
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij
      double precision rxuk,ryuk,rzuk,rxujk,ryujk,rzujk
      double precision duma,bij,rij,rjk,bjk,thetac
      double precision rki,rxuki,ryuki,rzuki,bki

      loverlap = .false.
      vthree = 0.0d0
      call twh_nmolty(GLB_GET,nmolty)

      do imove = 1,nchainmove
         icord = chaininfo(imove,testcord)
         if ( icord .ne. CRD_NONE ) then
c           --- get the chain information and coordinates
            ichain = chaininfo(imove,3)
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_ntype(GLB_GET,imolty,1,ntii)
            call twh_coordinates(GLB_GET,icord,ichain,1,rxui,ryui,rzui)
c           --- run through all existing atoms and also chaininfo atoms
c           --- with a higher index
            do jmolty = 0,nmolty
               if ( jmolty .eq. 0 ) then
                  istart = imove + 1
                  iend = nchainmove
               else
                  istart = 1
                  call twh_ncmt(GLB_GET,activebox,jmolty,iend)
               endif
               do jlist = istart,iend
                  lcompute = .true.
                  if ( jmolty .eq. 0 ) then
c                    --- get from the chaininfo
                     jcord = chaininfo(jlist,testcord)
                     if ( testcord .eq. CRD_NONE ) then
                        lcompute = .false.
                     else
c                       --- get from the chaininfo
                        jchain = chaininfo(jlist,3)
                        call twh_moltyp(GLB_GET,jchain,testmolty)
                        call twh_ntype(GLB_GET,testmolty,1,ntjj)
                     endif
                  else
c                    --- get from real molecules in this box
                     call twh_chainlist(GLB_GET,jlist,activebox,jmolty
     &                    ,jchain)
                     call twh_ntype(GLB_GET,jmolty,1,ntjj)
c                    --- make sure this is not a chaininfo
                     do testmove = 1,nchainmove
                        testchain = chaininfo(testmove,3)
                        if ( testchain .eq. jchain ) lcompute = .false.
                     enddo
                     jcord = CRD_REAL
                  endif
                  if ( lcompute ) then
c                    --- get jchain coordinates
                     call twh_coordinates(GLB_GET,jcord,jchain,1
     &                    ,rxuj,ryuj,rzuj)
c                    --- compute the distance
                     call twh_mimage(activebox,rxui,ryui,rzui
     &                    ,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij)
                     rij = twh_distance(rxuij,ryuij,rzuij)
c                    --- compute the bondorder
                     call twh_bondorder(lfailure,ntii,ntjj,rij,bij,duma
     &                    ,lzero)
                     if ( lfailure ) return
                     if ( .not. lzero ) then
c                       --- scan other atoms to make a triplet
                        do kmolty = 0,nmolty
                           if ( kmolty .eq. 0 ) then
                              kstart = imove + 1
                              kend = nchainmove
                           else
                              kstart = 1
                              call twh_ncmt(GLB_GET,activebox,kmolty
     &                             ,kend)
                           endif
                           do klist = kstart,kend
                              kcompute = .true.
                              if ( kmolty .eq. 0 ) then
                                 kchain = chaininfo(klist,3)
                                 kcord = chaininfo(klist,testcord)
                                 call twh_moltyp(GLB_GET,kchain
     &                                ,testmolty)
                                 call twh_ntype(GLB_GET,testmolty,1
     &                                ,ntkk)
                              else
                                 call twh_chainlist(GLB_GET,klist
     &                                ,activebox,kmolty,kchain)
                                 kcord = CRD_REAL
                                 call twh_ntype(GLB_GET,kmolty,1,ntkk)
c                                --- make sure this is not in chaininfo
                                 do testmove = 1,nchainmove
                                    testchain = chaininfo(testmove,3)
                                    if ( testchain .eq. kchain ) then
                                       kcompute = .false.
                                    endif
                                 enddo
                              endif
                              if ( ichain .eq. jchain .or.
     &                             jchain .eq. kchain .or.
     &                             ichain .eq. kchain ) then
                                 kcompute = .false.
                              endif
                              if ( kcompute ) then
c                                --- see if i-j-k is valid by checking
c                                --- the j-k bondorder
                                 call twh_coordinates(GLB_GET,kcord
     &                                ,kchain,1,rxuk,ryuk,rzuk)
                                 call twh_mimage(activebox,rxuk,ryuk
     &                                ,rzuk,rxuj,ryuj,rzuj,rxujk,ryujk
     &                                ,rzujk)
                                 rjk = twh_distance(rxujk,ryujk,rzujk)

                                 call twh_bondorder(lfailure,ntjj,ntkk
     &                                ,rjk,bjk,duma,kzero)
                                 if ( lfailure ) return
                                 if ( .not. kzero ) then
c                                   --- we have a valid i-j-k so compute
c                                   --- the three-body energy
                                    thetac = (rxuij*rxujk + ryuij*ryujk 
     &                                   + rzuij*rzujk)/(rij*rjk)
                                    vthree = vthree + 
     &                                   twh_vthreebody(lfailure,rij,rjk
     &                                   ,thetac,ntii,ntjj,ntkk)
                                    if ( lfailure ) return
                                 endif
c                                --- consider a k-i-j interaction only
c                                --- if k is later in the list than j
                                 jcompute = .true.
                                 if ( kmolty .lt. jmolty ) then
                                    jcompute = .false.
                                 elseif ( kmolty .eq. jmolty .and.
     &                                   klist .lt. jlist ) then
                                    jcompute = .false.
                                 endif
                                 if ( jcompute ) then
c                                   --- see if k-i-j is valid by
c                                   --- checking the k-i bondorder
                                    call twh_mimage(activebox,rxuk,ryuk
     &                                   ,rzuk,rxui,ryui,rzui
     &                                   ,rxuki,ryuki,rzuki)
                                    rki = twh_distance(rxuki,ryuki
     &                                   ,rzuki)
                                    call twh_bondorder(lfailure,ntkk
     &                                   ,ntii,rki,bki,duma,kzero)
                                    if ( lfailure ) return
                                    if ( .not. kzero ) then
c                                      --- we have a valid k-i-j so 
c                                      --- compute the three-body energy
                                       thetac = (-rxuij*rxuki
     &                                      - ryuij*ryuki 
     &                                      - rzuij*rzuki)/(rki*rij)
                                       vthree = vthree + 
     &                                      twh_vthreebody(lfailure,rki
     &                                      ,rij,thetac,ntkk,ntii,ntjj)
                                       if ( lfailure ) return
                                    endif
                                 endif
                              endif
                           enddo
                        enddo
                     endif
                  endif
               enddo
            enddo
         endif
      enddo

      return
      end


      subroutine twh_ewald_sum(activebox,vrecip)
c     ******************************************************************
c     * computes the new and old ewald sum using the currently stored  *
c     * values                                                         *
c     *                                                                *
c     * originally written 10-24-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer activebox
      double precision vrecip
      dimension vrecip(2)
c     --- local variables
c     ---- integer scalars
      integer ic,numvect
c     --- double precision scalars
      double precision prefact,ssumrn,ssumr,ssumi,ssumin

      call twh_ewald_numvect(GLB_GET,activebox,numvect)
      do ic = 1,numvect
         call twh_ewald_prefact(GLB_GET,CRD_REAL,ic,activebox,prefact)
         call twh_ewald_ssumr(GLB_GET,CRD_REAL,ic,activebox,ssumr)
         call twh_ewald_ssumr(GLB_GET,CRD_NEW,ic,activebox,ssumrn)
         call twh_ewald_ssumi(GLB_GET,CRD_REAL,ic,activebox,ssumi)
         call twh_ewald_ssumi(GLB_GET,CRD_NEW,ic,activebox,ssumin)
         vrecip(CRD_NEW) = vrecip(CRD_NEW)
     &        + (ssumrn*ssumrn + ssumin*ssumin)*prefact
         vrecip(CRD_OLD) = vrecip(CRD_OLD)
     &        + (ssumr*ssumr + ssumi*ssumi)*prefact
      enddo
      return
      end

      function twh_ewald_correct(ibox,imolty,iunit,junit,rij)
c     ******************************************************************
c     * computes the ewald intramolecular correction term in ibox for  *
c     * units iunit and junit of molecule type imolty                  *
c     *                                                                *
c     * originally written 10-24-2006 by M.G. Martin                   *
c     * last modified 11-02-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DERFC
#define FUNCTION_EWALD_CORRECT
#include "functions.h"
c     ---- variables passed to/from the function
      integer ibox,imolty,iunit,junit
      double precision rij,chargei,chargej,calp

      call twh_qqatom(GLB_GET,imolty,iunit,chargei)
      call twh_qqatom(GLB_GET,imolty,junit,chargej)
      call twh_ewald_calp(GLB_GET,ibox,calp)
      twh_ewald_correct = chargei*chargej*(twh_derfc(calp*rij)-1.0d0)
     &     /rij

      return
      end

      function twh_ewald_self(ibox,imolty,iunit)
c     ******************************************************************
c     * computes the ewald self interaction energy                     *
c     *                                                                *
c     * originally written 10-24-2006 by M.G. Martin                   *
c     * last modified 11-05-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EWALD_SELF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ibox,imolty,iunit
c     --- local variables
      double precision charge,invrootpi,calp

      call twh_constant_invrootpi(GLB_GET,invrootpi)
      call twh_qqatom(GLB_GET,imolty,iunit,charge)
      call twh_ewald_calp(GLB_GET,ibox,calp)
c     --- self interaction term
      twh_ewald_self = -charge*charge*calp*invrootpi

      return
      end
