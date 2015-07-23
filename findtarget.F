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
      subroutine twh_findtarget(lfailure,swaptype,box,imolty,ichain
     &,iunit,jmolty,jchain,junit,regstart,regend,pointp,num_in,num_out
     &,lfound,xtarget,ytarget,ztarget)
c     ******************************************************************
c     * wrapper around the findtarget routine that determines some     *
c     * constants required for array dimensions and then calls the     *
c     * findtarget_main routine                                        *
c     *                                                                *
c     * originally written 06-25-2009 by M.G. Martin                   *
c     * last modified 06-25-2009 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lfound
      integer swaptype,box,imolty,ichain,iunit,jmolty,jchain,junit
      integer regstart,regend,num_in,num_out
      integer pointp
      double precision xtarget,ytarget,ztarget
      dimension xtarget(2),ytarget(2),ztarget(2)
c     --- local variables
      integer nchain

      call twh_nchain(GLB_GET,nchain)

      call twh_findtarget_main(lfailure,swaptype,box,imolty,nchain
     &     ,ichain,iunit,jmolty,jchain,junit,regstart,regend,pointp
     &     ,num_in,num_out,lfound,xtarget,ytarget,ztarget)

      return
      end

      subroutine twh_findtarget_main(lfailure,swaptype,box,imolty,nchain
     &     ,ichain,iunit,jmolty,jchain,junit,regstart,regend,pointp
     &     ,num_in,num_out,lfound,xtarget,ytarget,ztarget)
c     ******************************************************************
c     * finds a target molecule for all of the swapmoves               *
c     * in these moves you need to find another molecule that is in a  *
c     * certain region of space.  If a suitable molecule is found then *
c     * it is returned, otherwise lfound is set to false and the move  *
c     * will be rejected                                               *
c     *                                                                *
c     * originally written 09-13-2002 by M.G. Martin                   *
c     * last modified 04-21-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_GET_NBOXI
#define FUNCTION_RANDOM
#include "functions.h"
c     --- debug toggle
#define DEBUG_FINDTARGET 0
c     --- variables passed to/from the subroutine
      logical lfailure,lfound
      integer swaptype,box,imolty,nchain,ichain,iunit,jmolty,jchain
      integer junit,regstart,regend,num_in,num_out
      integer pointp
      double precision xtarget,ytarget,ztarget
      dimension xtarget(2),ytarget(2),ztarget(2)
c     --- local variables
      integer ngood,goodchain,itest,jtest,ktest,pointt,kchain
      integer ncmti,ncmtj,numiunit,numjunit
      dimension goodchain(nchain)
      double precision rxij,ryij,rzij,dist,xtemp,ytemp,ztemp,pmtest
      double precision testrad

c     --- initialize variables
      lfound = .false.
      num_in = 0
      num_out = 0
c     --- retrieve constants
      call twh_nunit(GLB_GET,imolty,numiunit)
      call twh_ncmt(GLB_GET,box,imolty,ncmti)

      if ( swaptype .lt. 2 .or. 
     &     swaptype .eq. 5 .or. 
     &     swaptype .eq. 6) then
         pointp = idint( dble(ncmti)*twh_random() ) + 1
         call twh_chainlist(GLB_GET,pointp,box,imolty,ichain)
#if DEBUG_FINDTARGET
      write(6,*) 'box,imolty,ncmt:',box,imolty,ncmti
      write(6,*) 'pointp,ichain:',pointp,ichain
#endif
c        --- safety checks
         call twh_moltyp(GLB_GET,ichain,itest)
         if ( itest .ne. imolty ) then
            write(6,*) 'FINDTARGET: ichain,imolty,moltyp(ichain)'
     &           ,ichain,imolty,itest
            write(6,*) 'swaptype,box,ncmt',swaptype,box,ncmti
            write(6,*) 'pointp',pointp
            lfailure = .true.
            return
         endif
         if ( twh_get_nboxi(ichain) .ne. box ) then
            write(6,*) 'FINDTARGET: ichain,box,nboxi(ichain)'
     &           ,ichain,box,twh_get_nboxi(ichain)
            lfailure = .true.
            return
         endif
c        --- choose which atom to use as the starting unit
         iunit = idint( twh_random()*numiunit ) + 1
c         --- we are done and we found what we need
         lfound = .true.
      elseif ( swaptype .eq. 2 ) then
c        --- select the chain that will move
         pointp = idint( dble(ncmti)*twh_random() ) + 1
         call twh_chainlist(GLB_GET,pointp,box,imolty,ichain)
c        --- safety checks
         call twh_moltyp(GLB_GET,ichain,itest)
         if ( itest .ne. imolty ) then
            write(6,*) 'FINDTARGET: ichain,imolty,moltyp(ichain)'
     &           ,ichain,imolty,itest
            write(6,*) 'swaptype,box,ncmt',swaptype,box,ncmti
            lfailure = .true.
            return
         endif
         if ( twh_get_nboxi(ichain) .ne. box ) then
            write(6,*) 'FINDTARGET: ichain,box,nboxi(ichain)'
     &           ,ichain,box,twh_get_nboxi(ichain)
            lfailure = .true.
            return
         endif
c        --- choose which atom to use as the starting unit
         iunit = idint( twh_random()*numiunit ) + 1

c        --- we need to select a target molecule of type jmolty at 
c        --- random from the molecules of that type in box
         call twh_ncmt(GLB_GET,box,jmolty,ncmtj)
         if ( ncmtj .eq. 0 ) then
c           --- no molecules of the target type, reject move
            return
         elseif ( imolty .eq. jmolty ) then
c           --- we have to make sure not to select the same molecule
c           --- that we are already using for the move
            if ( ncmti .le. 1 ) then
c              --- there is only one molecule of this type in the box
c              --- therefore we cannot perform the move as there is no
c              --- other molecule to use as a reference
               return
            else
c              --- choose a molecule of the appropriate type that is 
c              --- not the same molecule as the target
               do while ( .not. lfound ) 
                  pointp = idint(dble(ncmtj)*twh_random())+1
                  call twh_chainlist(GLB_GET,pointp,box,jmolty,jchain)
                  if ( jchain .ne. ichain ) lfound = .true.
               enddo
            endif
         else
c           --- jmolty is of a different type than imolty, choose one
c           --- at random
            pointp = idint( dble(ncmtj)*twh_random() ) + 1
            call twh_chainlist(GLB_GET,pointp,box,jmolty,jchain)
         endif
c        --- select an atom from jmolty at random.  this will be the
c        --- target atom reference for insertion
         call twh_nunit(GLB_GET,jmolty,numjunit)
         junit = idint(twh_random()*dble(numjunit))+1
c        --- assign coordinates for the target atom
         call twh_coordinates(GLB_GET,CRD_REAL,jchain,junit
     &        ,xtarget(1),ytarget(1),ztarget(1))
c        --- determine if iunit of molecule ichain is in the inner or 
c        --- outer volume determined by junit
         call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &        ,xtemp,ytemp,ztemp)
         call twh_mimage(box,xtarget(1),ytarget(1),ztarget(1)
     &        ,xtemp,ytemp,ztemp,rxij,ryij,rzij)
         dist = twh_distance(rxij,ryij,rzij)
         call twh_avb1rad(GLB_GET,testrad)
         if ( dist .le. testrad ) then
c           --- we are starting in the inner region
            regstart = 0
         else
c           --- we are starting in the outer region
            regstart = 1
         endif
c        --- determine whether the target region is inner or outer
         call twh_pmavb1in(GLB_GET,pmtest)
         if ( twh_random() .lt. pmtest ) then
            regend = 0
         else
            regend = 1
         endif

      elseif ( swaptype .eq. 3 ) then
c        --- this move logic appears to be reversed from the others 
c        --- as the first chain we select is the one that is not moving
c        --- (refered to as the j chain)
c        --- where the other moves first select the chain that does 
c        --- move.  This is due to the different nature of the AVB2
c        --- algorithm

         call twh_ncmt(GLB_GET,box,jmolty,ncmtj)
         if ( ncmtj .eq. 0 ) then
c           --- no molecules of the target type, reject move
            return
         endif

c        --- select the chain that is stationary (jchain)
         pointp = idint( dble(ncmtj)*twh_random() ) + 1
         call twh_chainlist(GLB_GET,pointp,box,jmolty,jchain)
c        --- safety checks
         call twh_moltyp(GLB_GET,jchain,itest)
         if ( itest .ne. jmolty ) then
            write(6,*) 'FINDTARGET: jchain,jmolty,moltyp(jchain)'
     &           ,jchain,jmolty,itest
            lfailure = .true.
            return
         endif
         if ( twh_get_nboxi(jchain) .ne. box ) then
            write(6,*) 'FINDTARGET: jchain,box,nboxi(jchain)'
     &           ,jchain,box,twh_get_nboxi(jchain)
            lfailure = .true.
            return
         endif

c        --- determine whether to do an in->out or an out->in
         call twh_pmavb2in(GLB_GET,pmtest)
         if ( twh_random() .lt. pmtest ) then
c           --- in -> out
            regstart = 0
            regend = 1
         else
c           --- out -> in
            regstart = 1
            regend = 0
         endif
c        --- choose an atom in jmolty to use as the measuring atom
         call twh_nunit(GLB_GET,jmolty,numjunit)
         junit = idint( twh_random()*dble(numjunit) ) + 1
c        --- choose which atom in imolty to use as the starting unit
         iunit = idint( twh_random()*dble(numiunit) ) + 1

c        --- store the coordinates for the jchain reference molecule
         call twh_coordinates(GLB_GET,CRD_REAL,jchain,junit
     &        ,xtarget(1),ytarget(1),ztarget(1))
c        --- make a list of all the molecule of the appropriate type
c        --- that are in the correct region
         ngood = 0
         do itest = 1,ncmti
            call twh_chainlist(GLB_GET,itest,box,imolty,ichain)
            if ( ichain .ne. jchain ) then
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
               call twh_mimage(box,xtarget(1),ytarget(1),ztarget(1)
     &              ,xtemp,ytemp,ztemp,rxij,ryij,rzij)
               dist = twh_distance(rxij,ryij,rzij)
               call twh_avb2rad(GLB_GET,testrad)
               if ( dist .lt. testrad ) then
                  num_in = num_in + 1
                  if ( regstart .eq. 0 ) then
c                    --- this is in the inner region
                     ngood = ngood + 1
                     goodchain(ngood) = ichain
                  endif
               else
                  num_out = num_out + 1
                  if ( regstart .eq. 1 ) then
c                    --- this is in the outer region
                     ngood = ngood + 1
                     goodchain(ngood) = ichain
                  endif
               endif
            endif
         enddo

c        --- if there are no chains to choose from then reject move
         if ( ngood .eq. 0 ) return
c        --- we made it this far so we must have found one
         lfound = .true.
c        --- choose one of the good chains
         pointp = idint( twh_random()*dble(ngood)) + 1
         ichain = goodchain(pointp)

      elseif ( swaptype .eq. 4 ) then
c        --- the first chain we select is the one that is not moving
c        --- (refered to as the j chain)
c        --- retrieve the avb3rad
         call twh_avb3rad(GLB_GET,testrad)
         call twh_ncmt(GLB_GET,box,jmolty,ncmtj)
         if ( ncmtj .lt. 2 ) then
c           --- not enough molecules of the target type, reject move
            return
         endif

c        --- select the chain that is stationary (jchain)
         pointp = idint( dble(ncmtj)*twh_random() ) + 1
         call twh_chainlist(GLB_GET,pointp,box,jmolty,jchain)
c        --- safety checks
         call twh_moltyp(GLB_GET,jchain,itest)
         if ( itest .ne. jmolty ) then
            write(6,*) 'FINDTARGET: jchain,jmolty,moltyp(jchain)'
     &           ,jchain,jmolty,itest
            lfailure = .true.
            return
         endif
         if ( twh_get_nboxi(jchain) .ne. box ) then
            write(6,*) 'FINDTARGET: jchain,box,nboxi(jchain)'
     &           ,jchain,box,twh_get_nboxi(jchain)
            lfailure = .true.
            return
         endif

c        --- the move is always in -> in
         regstart = 0
         regend = 0
c        --- choose an atom in jmolty to use as the measuring atom
         call twh_nunit(GLB_GET,jmolty,numjunit)
         junit = idint( twh_random()*dble(numjunit) ) + 1
c        --- choose which atom in imolty to use as the starting unit
         iunit = idint( twh_random()*dble(numiunit) ) + 1
c        --- store the coordinates for the jchain reference molecule
         call twh_coordinates(GLB_GET,CRD_REAL,jchain,junit
     &        ,xtarget(1),ytarget(1),ztarget(1))
c        --- first we need to find another molecule of type jmolty that
c        --- does not have overlapping inner regions with jchain
         ngood = 0
         do itest = 1,ncmtj
            call twh_chainlist(GLB_GET,itest,box,jmolty,kchain)
            if ( kchain .ne. jchain ) then
               call twh_coordinates(GLB_GET,CRD_REAL,kchain,junit
     &              ,xtemp,ytemp,ztemp)
               call twh_mimage(box,xtarget(1),ytarget(1),ztarget(1)
     &              ,xtemp,ytemp,ztemp,rxij,ryij,rzij)
               dist = twh_distance(rxij,ryij,rzij)
               if ( dist .gt. 2.0d0*testrad ) then
                  ngood = ngood + 1
                  goodchain(ngood) = kchain
               endif
            endif
         enddo
c        --- if there are no chains to choose from then reject move
         if ( ngood .eq. 0 ) return
c        --- choose one of the good chains
         pointt = idint( twh_random()*dble(ngood)) + 1
         kchain = goodchain(pointt)
c        --- set the secondary targets
         call twh_coordinates(GLB_GET,CRD_REAL,kchain,junit
     &        ,xtarget(2),ytarget(2),ztarget(2))
c        --- choose molecule of type imolty that is in the inner region 
c        --- of jchain and count the number in the inner regions of 
c        --- jchain and kchain
         ngood = 0
         do itest = 1,ncmti
            call twh_chainlist(GLB_GET,itest,box,imolty,ichain)
            if ( ichain .ne. jchain .and. ichain .ne. kchain ) then
c              --- count the number in the inner region of jchain
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
               call twh_mimage(box,xtarget(1),ytarget(1),ztarget(1)
     &              ,xtemp,ytemp,ztemp,rxij,ryij,rzij)
               dist = twh_distance(rxij,ryij,rzij)
               if ( dist .lt. testrad ) then
                  num_in = num_in + 1
                  ngood = ngood + 1
                  goodchain(ngood) = ichain
               endif
c              --- count the number in the inner region of kchain
c              --- this is stored in num_out
               call twh_mimage(box,xtarget(2),ytarget(2),ztarget(2)
     &              ,xtemp,ytemp,ztemp,rxij,ryij,rzij)
               dist = twh_distance(rxij,ryij,rzij)
               if ( dist .lt. testrad ) then
                  num_out = num_out + 1
               endif
            endif
         enddo

c        --- if there are no chains to choose from then reject move
         if ( ngood .eq. 0 ) return
c        --- we made it this far so we must have found one
         lfound = .true.
c        --- choose one of the good chains
         pointp = idint( twh_random()*dble(ngood)) + 1
         ichain = goodchain(pointp)

         call twh_moltyp(GLB_GET,ichain,itest)
         call twh_moltyp(GLB_GET,jchain,jtest)
         call twh_moltyp(GLB_GET,kchain,ktest)
         if ( itest .ne. imolty 
     &        .or. jtest .ne. jmolty
     &        .or. ktest .ne. jmolty ) then
            write(6,*) 'FINDTARGET: problem in avb3'
            write(6,*) 'imolty,jmolty',imolty,jmolty
            write(6,*) 'ichain,moltyp',ichain,itest
            write(6,*) 'jchain,moltyp',jchain,jtest
            write(6,*) 'kchain,moltyp',kchain,ktest
            lfailure = .true.
            return
         endif
      endif

      return
      end
