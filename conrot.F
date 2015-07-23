#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2003-2012 Marcus G. Martin                       *
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
      subroutine twh_conrot(lfailure,crtype)
c     ******************************************************************
c     * concerted rotation move algorithm                              *
c     *                                                                *
c     * crtype = 0: concerted rotation move for protein systems which  *
c     *             uses the lbackbone and backbone related arrays     *
c     *             which are not defined for non-protein systems      *
c     *                                                                *
c     * crtype = 1: concerted rotation move for all other systems      *
c     *             for which the backbone arrays are not defined      *
c     *             and it is not so important to just focus on the    *
c     *             backbone                                           *
c     *                                                                *
c     * originally written 05-01-2003 by C.D. Lorenz                   *
c     * last modified 08-30-2012 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "dvpm.h"
#define FUNCTION_DISTANCE
#define FUNCTION_EXPON
#define FUNCTION_GET_IJAA0
#define FUNCTION_GET_IJAA1
#define FUNCTION_GET_IJAA2
#define FUNCTION_GET_IJBEN2
#define FUNCTION_GET_IJBEN3
#define FUNCTION_GET_IJBNBN1
#define FUNCTION_GET_IJBNBN2
#define FUNCTION_GET_IJIMPROP2
#define FUNCTION_GET_IJIMPROP3
#define FUNCTION_GET_IJIMPROP4
#define FUNCTION_GET_IJTOR2
#define FUNCTION_GET_IJTOR3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_INAA
#define FUNCTION_GET_INBEN
#define FUNCTION_GET_INBNBN
#define FUNCTION_GET_INIMPROP
#define FUNCTION_GET_INTOR
#define FUNCTION_GET_INVIB
#define FUNCTION_RANDOM
#define FUNCTION_VANGANG
#define FUNCTION_VANGLE
#define FUNCTION_VBONBON
#define FUNCTION_VIMPROPER
#define FUNCTION_VTORSION
#include "functions.h"
c     --- debug toggle
#define DEBUG_CONROT 0
c     --- variables passed to/from the subroutine
      logical lfailure
      integer crtype
c     --- local variables
c     --- parameters
      integer searchsize
      parameter (searchsize = 11)
      double precision tol_secant
      parameter (tol_secant = 1.0d-15)
      integer max_soln
      parameter (max_soln=16)
      integer max_iter_sturm
      parameter (max_iter_sturm=100)
      integer max_iter_secant
      parameter (max_iter_secant=20)
c     --- logical scalars
      logical lfound,lgetcom,lcompute,laccept,lbadstereo,ldone,lewald
      logical lpartial,lintra,lbonded,lvalid,scp_eval_derivative
c     -- logical arrays
      logical loverlap
      dimension loverlap(2)
c     --- integer scalars
      integer icmol,imolty,iunit,ichain,ibox,units,iu
      integer unit1,unit2,ivcount
      integer unitss,unitss1,nsoln,i,j,isuccess,solnno
      integer atomid,nsoln2,nmolty
      integer ioldnew,icompute,itor,ibead,itest,ivib,iutwo,nmtype
      integer iend,iuone,ibend,iuthree,iaa,inum,ivalid,testcord
      integer nextback,prevback,junit,currentunit
      integer nchainmove,natommove,nworklist,nvib,testunit,branchunit
      integer iworklist,itestlist,ncmt,ipoint,ensemble
c     --- integer arrays
      integer dcrs,crback,torbead
      dimension dcrs(NUMAX),crback(40),torbead(4)
      integer vibpointer,prevunit,visitcount
      dimension vibpointer(NUMAX),prevunit(NUMAX),visitcount(NUMAX) 
      integer chaininfo
      dimension chaininfo(1,3)
      integer atomlist
      dimension atomlist(searchsize,3)
      integer globalcord
      dimension globalcord(2)
      integer workitem
      dimension workitem(NUMAX)
c     --- double precision scalars
      double precision ranmolty,dchain,olddeter,newdeter,dvdihed
      double precision xtemp,ytemp,ztemp,cr_accept,ddum,dnsoln
      double precision rxuij,ryuij,rzuij,duma,dumb,dumc
      double precision delta_gamma,driver_alpha,cosda,dist12,dist23
      double precision x_unitvec,y_unitvec,z_unitvec,fdum,driver_gamma
      double precision new_gamma,newv_x,newv_y,newv_z,drive_dx,drive_dy
      double precision drive_dz,beta
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,twopi,pmtest
      double precision aa13_min_sqr,aa13_max_sqr
c     --- double precision arrays
      double precision initcoords,solncoords
      dimension initcoords(9,3),solncoords(9,3,max_soln)
      double precision dvtot,dvintra,dvinter,dvext,dvelect,dvsolv
      double precision xtmp,ytmp,ztmp,dvtors,dvimpr
      double precision dvvib,blength,dvbend,dvquantum
      double precision dlinter,dlelect
      dimension xtmp(4),ytmp(4),ztmp(4),blength(2)
      double precision tempx,tempy,tempz,rootunit
      dimension tempx(3),tempy(3),tempz(3),rootunit(4,3)
      double precision len0,b_ang0,t_ang0,numvalid
      dimension len0(6),b_ang0(7),t_ang0(2)
      double precision dvf
      dimension dvf(MAX_FOREIGN_LAMBDA)
      double precision delta
      dimension delta(0:3)
      double precision lcxi
      dimension lcxi(3)
      double precision len_aa
      dimension len_aa(3)
      double precision eta
      dimension eta(3)
      double precision ddum2
      dimension ddum2(2)
      integer iboxa
      dimension iboxa(1)

      external twhinitializesturm

#if DEBUG_CONROT
      write(6,*) 'start CONROT',crtype
#endif

c     --- retrieve constants
      call twh_constant_twopi(GLB_GET,twopi)
      call twh_nmolty(GLB_GET,nmolty)
      call twh_lewald(GLB_GET,lewald)

c     --- select a molecule type at random 
      ranmolty  = twh_random()
      do icmol = 1,nmolty
         if (crtype .eq. 0) then
            call twh_pmcrbmt(GLB_GET,icmol,pmtest)
         else if (crtype .eq. 1) then
            call twh_pmcrmt(GLB_GET,icmol,pmtest)
         endif
         if ( ranmolty .lt. pmtest ) then
            imolty = icmol
            ranmolty = 2.0d0
         endif
      enddo

c     --- select a chain of the determined molecule type
      call twh_ensemble(GLB_GET,ensemble)
      if ( ensemble .eq. ENS_muVT ) then
c        --- box 1
         ibox = 1
c        --- make sure this type exists in box 1
         call twh_ncmt(GLB_GET,ibox,imolty,ncmt)
         if ( ncmt .eq. 0 ) return
c        --- select molecule type imolty in box ibox
         ipoint = idint( dble(ncmt)*twh_random() + 1 )
         call twh_chainlist(GLB_GET,ipoint,ibox,imolty,ichain)
      else
         call twh_nmolectyp(GLB_GET,imolty,nmtype)
         dchain = dble(nmtype)
         ipoint = idint((dchain*twh_random()) + 1)
         call twh_parall(GLB_GET,imolty,ipoint,ichain)
      endif
c     --- determine what box this part of the molecule is in
      call twh_nboxi(GLB_GET,ichain,ibox)

c     --- store the number of units in this molecule
      call twh_nunit(GLB_GET,imolty,units)

c     --- store current coordinates of the molecule as the new trial
      do iu = 1, units
         call twh_coordinates(GLB_GET,CRD_REAL,ichain,iu,xtemp,ytemp
     &        ,ztemp)
         call twh_coordinates(GLB_SET,CRD_TEMP,ichain,iu,xtemp,ytemp
     &        ,ztemp)
      enddo
c     --- make sure there is a valid option for this molecule type
      call twh_crvalidcount(GLB_GET,imolty,ivcount)
      if ( ivcount .eq. 0 ) then
         write(6,*) 'CONROT: molecule type with invalid crvalidcount'
         write(6,*) 'there is no way to perform this move on molecule'
     &        ,' type',imolty
         lfailure = .true.
         return
      endif
c     --- select a unit from the chosen molecule that is part of the 
c     --- backbone as defined in findrings
      numvalid = dble(ivcount)
      ivalid = idint(numvalid*twh_random()) + 1
      call twh_crvalidlist(GLB_GET,imolty,ivalid,unit1)
#if DEBUG_CONROT
      write(6,*) 'calling crstart'
      write(6,*) 'crtype:',crtype
      write(6,*) 'crvalidcount:',ivcount
      write(6,*) 'numvalid:',numvalid
      write(6,*) 'ivalid:',ivalid
      write(6,*) 'imolty,unit1:',imolty,unit1
#endif
      call twh_crstart(crtype,imolty,unit1,dcrs,crback,isuccess)
#if DEBUG_CONROT
      write(6,*) 'back from crstart'
      write(6,*) 'crtype,imolty,unit1,isuccess',crtype,imolty,unit1
     &     ,isuccess
#endif
c     --- store id numbers for the last two backbone atoms
      if (isuccess .eq. 1) then
         unitss1 = crback(searchsize-2) 
         unitss = crback(searchsize-1)
c        --- update number of CONROT/CRBACK move trys
         if (crtype .eq. 0) then
            call twh_bncrback(GLB_INCR,imolty,unit1,ibox,1)
         elseif (crtype .eq. 1) then
            call twh_bnconrot(GLB_INCR,imolty,unit1,ibox,1)
         endif   
c        --- determine second vertice of triangle randomly 
         do iunit = 1, units
            dcrs(iunit) = 0
         enddo 
         unit2 = crback(2)
         dcrs(unit1) = 1
         dcrs(unit2) = 2    
         do iunit = 1,units
c           --- initialize the vibpointer
            vibpointer(iunit) = 0
c           --- initialize the prevunit
            prevunit(iunit) = 0
c           --- initialize visitcount
            visitcount(iunit) = 0
         enddo
         currentunit = unit2
         ldone = .false.
         vibpointer(currentunit) = 0
         do while ( .not. ldone)
            vibpointer(currentunit) = vibpointer(currentunit) + 1
            if (vibpointer(currentunit) .gt.
     &           twh_get_invib(imolty,currentunit)) then
               if (prevunit(currentunit) .eq. 0) then
c                 --- completely finished with this current unit
                  ldone = .true.
               else
c                 --- move back to previous atom to see if can find
c                 --- necessary number of chain members
                  currentunit = prevunit(currentunit)
               endif
            else
               call twh_ijvib(GLB_GET,imolty,currentunit
     &              ,vibpointer(currentunit),junit)
               if (junit .ne. crback(1)) then
c                 --- do not want to double back to first atom in 
c                 --- torsion arm
                  if (dcrs(junit) .eq. 0) then
c                    --- have not visited this atom yet
                     dcrs(junit) = dcrs(currentunit) + 1
                     prevunit(junit) = currentunit
                     currentunit = junit
                     vibpointer(currentunit) = 0
                     visitcount(currentunit) = 
     &                    visitcount(currentunit) + 1
                     if (currentunit .eq. crback(searchsize)) then
c                       --- found path to the last atom in torsion arm
                        ldone = .true.
                     endif
                  endif
               endif
            endif
            if (visitcount(currentunit) .gt. NNBOND*NUMAX) then
c              --- we must be stuck in some strange periodic ring
c              --- structure.  Kick out and hope we figure it out from
c              --- the next startunit
               ldone = .true.
            endif
         enddo

         call twhinitializesturm(tol_secant,max_iter_sturm
     &        ,max_iter_secant)

         do i = 1,3
            tempx(i) = 0.0d0
            tempy(i) = 0.0d0
            tempz(i) = 0.0d0
         enddo

         do iunit = 1,searchsize-1
            atomid = crback(iunit)
            if (iunit .ge. 2) then
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,atomid
     &              ,xtemp,ytemp,ztemp)
               initcoords(iunit-1,1) = xtemp
               initcoords(iunit-1,2) = ytemp
               initcoords(iunit-1,3) = ztemp
            endif
            if (iunit .le. 3) then
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,atomid
     &              ,xtemp,ytemp,ztemp)
               tempx(iunit) = xtemp
               tempy(iunit) = ytemp
               tempz(iunit) = ztemp
            endif
         enddo 

         call twh_initialize_loop_closure(lfailure,ibox,initcoords,len0,
     &        b_ang0,t_ang0,aa13_min_sqr,aa13_max_sqr,delta,lcxi,len_aa
     &        ,eta)
         if ( lfailure ) return

c        --- apply driver angle to the second atom in the 
c        --- chain that has been chosen for the concerted 
c        --- rotation move 
c        --- first need to calculate the current angle 
c        --- between the three atoms
         driver_alpha = (twh_vangle(lfailure,loverlap(1),ibox,imolty,1,1
     &        ,-1,cosda,dist12,dist23,tempx,tempy,tempz))
         if ( lfailure ) return
         x_unitvec = (tempx(2) - tempx(1))/dist12
         y_unitvec = (tempy(2) - tempy(1))/dist12
         z_unitvec = (tempz(2) - tempz(1))/dist12
c        --- then need to calculate the gamma angle which 
c        --- which is the rotation around the z axis that
c        --- we are initially at.

         call twh_cone(lfailure,1,x_unitvec,y_unitvec,z_unitvec
     &        ,fdum,fdum,fdum,fdum,fdum)
         if ( lfailure ) return

         call twh_cone(lfailure,3,fdum,fdum,fdum,driver_alpha
     &        ,driver_gamma,x_unitvec,y_unitvec,z_unitvec)
         if ( lfailure ) return

         delta_gamma = twopi*twh_random() 
         new_gamma = driver_gamma + delta_gamma
c        --- finally get the unit vector that defines the new 
c        --- position along the cone at driver_alpha and new_gamma
         call twh_cone(lfailure,2,fdum,fdum,fdum,driver_alpha,new_gamma
     &        ,newv_x,newv_y,newv_z)           
         if ( lfailure ) return
         drive_dx = newv_x*dist23
         drive_dy = newv_y*dist23
         drive_dz = newv_z*dist23
         
         do i = 1,4
            do j = 1,3
               rootunit(i,j) = 0.0d0
            enddo
         enddo

c        --- assign coordinates of the four root units
c        --- in the concerted rotation move (#1,#2,#8,#9)
         rootunit(1,1) = initcoords(1,1)
         rootunit(1,2) = initcoords(1,2)
         rootunit(1,3) = initcoords(1,3)
         rootunit(2,1) = initcoords(1,1)+drive_dx
         rootunit(2,2) = initcoords(1,2)+drive_dy
         rootunit(2,3) = initcoords(1,3)+drive_dz
         rootunit(3,1) = initcoords(8,1)
         rootunit(3,2) = initcoords(8,2)
         rootunit(3,3) = initcoords(8,3)
         rootunit(4,1) = initcoords(9,1)
         rootunit(4,2) = initcoords(9,2)
         rootunit(4,3) = initcoords(9,3)
         nsoln = 0
         call twh_solv_3pep_poly(lfailure,max_soln,ibox,rootunit
     &        ,len0,b_ang0,solncoords,nsoln,aa13_min_sqr,aa13_max_sqr
     &        ,delta,lcxi,len_aa,eta)
         if ( lfailure ) return
#if DEBUG_CONROT
         write(6,*) "(twh_solv_3pep_poly): nsoln: ", nsoln
#endif

         if (nsoln .gt. 0) then               
            if (nsoln .gt. max_soln) then
               write(6,*) "More solutions found than expected max",
     &              max_soln
               write(6,*) "Increase max_soln in conrot.F"
               lfailure = .true.
               return
            endif 
c           --- because the original position is always returned as a 
c           --- solution we substract the first solution and then pick 
c           --- randomly from the remaining solutions
            dnsoln = dble(nsoln)
            solnno = idint((dnsoln*twh_random())+1)
c           --- place backbone atoms
            do iunit = 2,searchsize-1
               atomid = crback(iunit)
               xtemp = solncoords(iunit-1,1,solnno)
               ytemp = solncoords(iunit-1,2,solnno)
               ztemp = solncoords(iunit-1,3,solnno)
               call twh_coordinates(GLB_SET,CRD_TEMP,ichain,atomid,xtemp
     &              ,ytemp,ztemp)
            enddo

            do iunit = 2,searchsize-1
               atomid = crback(iunit)
c              --- rebuild the sidechains for all atoms after the first
               prevback = crback(iunit-1)
               nextback = crback(iunit+1)
c              --- build a list of all atoms that are bonded to atomid
c              --- in a direction other than prevback or nextback
               nworklist = 0
c              --- check the atoms bonded to atomid to see if there are
c              --- any valid ones to prime the worklist
               call twh_invib(GLB_GET,imolty,atomid,nvib)
               do ivib = 1,nvib
                  call twh_ijvib(GLB_GET,imolty,atomid,ivib,testunit)
                  if ( testunit .ne. prevback .and.
     &                 testunit .ne. nextback ) then
                     nworklist = nworklist + 1
                     workitem(nworklist) = testunit
                  endif
               enddo

               iworklist = 0
               do while ( iworklist .lt. nworklist )
                  iworklist = iworklist + 1
                  branchunit = workitem(iworklist)
                  call twh_invib(GLB_GET,imolty,branchunit,nvib)
                  do ivib = 1,nvib
                     call twh_ijvib(GLB_GET,imolty,branchunit,ivib
     &                    ,testunit)
                     lvalid = .true.
                     if ( testunit .eq. prevback .or.
     &                    testunit .eq. nextback ) then
                        lvalid = .false.
                     endif
                     itestlist = 0
                     do while ( lvalid .and. itestlist .lt. nworklist )
                        itestlist = itestlist + 1
                        if ( testunit .eq. workitem(itestlist) ) then
                           lvalid = .false.
                        endif
                     enddo
                     if ( lvalid ) then
c                       --- found a new valid addition
                        nworklist = nworklist + 1
                        workitem(nworklist) = testunit
                     endif
                  enddo
               enddo
               do iworklist = 1,nworklist
                  branchunit = workitem(iworklist)
                  call cr_angle2cart(ichain,branchunit,prevback
     &                 ,atomid,nextback)
               enddo
            enddo

c           --- conduct the reverse concerted rotation which is required
c           --- in order to determine the acceptance criteria

            call twhinitializesturm(tol_secant,max_iter_sturm
     &           ,max_iter_secant)         

            do iunit = 2,searchsize-1
               atomid = crback(iunit)
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,atomid
     &              ,xtemp,ytemp,ztemp)
               initcoords(iunit-1,1) = xtemp
               initcoords(iunit-1,2) = ytemp
               initcoords(iunit-1,3) = ztemp
            enddo

            call twh_initialize_loop_closure(lfailure,ibox,initcoords
     &           ,len0, b_ang0,t_ang0,aa13_min_sqr,aa13_max_sqr,delta
     &           ,lcxi,len_aa,eta)
            if ( lfailure ) return
            rootunit(1,1) = initcoords(1,1)
            rootunit(1,2) = initcoords(1,2)
            rootunit(1,3) = initcoords(1,3)
            rootunit(2,1) = initcoords(2,1)
            rootunit(2,2) = initcoords(2,2)
            rootunit(2,3) = initcoords(2,3)
            rootunit(3,1) = initcoords(8,1)
            rootunit(3,2) = initcoords(8,2)
            rootunit(3,3) = initcoords(8,3)
            rootunit(4,1) = initcoords(9,1)
            rootunit(4,2) = initcoords(9,2)
            rootunit(4,3) = initcoords(9,3)

            nsoln2 = 0
            call twh_solv_3pep_poly(lfailure,max_soln,ibox,rootunit
     &           ,len0,b_ang0,solncoords,nsoln2,aa13_min_sqr
     &           ,aa13_max_sqr,delta,lcxi,len_aa,eta)
            if ( lfailure ) return
c           --- always have initial configuration as a solution
            if (nsoln2 .eq. 0) then
               nsoln2 = 1
            endif
            lgetcom = .true.
c           --- calculate the determinant of the jacobian matrix at the 
c           --- initial and final positions of the backbone atoms so that 
c           --- these determinants may be used in calculating the 
c           --- acceptance ratio
#if DEBUG_CONROT
      write(6,*) 'calling calcjacob'
#endif
            call twh_calcjacob(ichain,crback,olddeter,newdeter)
c           --- calculate the new and old non-bonded inter molecular
c           --- energies
            dvimpr = 0.0d0
            nchainmove = 1
            chaininfo(1,CRD_OLD) = CRD_REAL
            chaininfo(1,CRD_NEW) = CRD_TEMP
            chaininfo(1,3) = ichain
            lpartial = .false.
            lintra = .true.
            lbonded = .false.
            call twh_energy_change(lfailure,ibox,nchainmove,chaininfo
     &           ,loverlap,lpartial,lgetcom,lintra,lbonded
     &           ,dvtot,dvintra,dvinter,dvext,dvelect,dvsolv,dvtors
     &           ,dvbend,dvvib,dvquantum,dvf,dlinter,dlelect)
            if ( lfailure ) return

            do ioldnew = CRD_NEW,CRD_OLD,(CRD_OLD-CRD_NEW)
               if ( loverlap(ioldnew) ) then
                  if (ioldnew .eq. 2) then
                     goto 500
                  else
                     write(6,*) 'CONROT: Overlap in old conformation'
                     lfailure = .true.
                     return
                  endif
               endif
            enddo

            do iunit = 1, units
               dcrs(iunit) = 0
            enddo 
#if DEBUG_CONROT
      write(6,*) 'searchsize:',searchsize
      do inum = 1,searchsize
         write(6,*) 'inum,crback:',inum,crback(inum)
      enddo
#endif
c           --- set backbone atoms to 1
            do inum = 1,searchsize
               currentunit = crback(inum)
               dcrs(currentunit) = 1
            enddo
c           --- set neighbors of backbone atoms to 2, unless also a 
c           --- backbone atom
            do inum = 1,searchsize
               currentunit = crback(inum)
               do ivib = 1,twh_get_invib(imolty,currentunit)
                  call twh_ijvib(GLB_GET,imolty,currentunit,ivib,iutwo)
                  if ( dcrs(iutwo) .eq. 0 ) dcrs(iutwo) = 2
               enddo
            enddo

            natommove = 0
c           --- include atoms that moved on the backbone
            do inum = 1,searchsize
               currentunit = crback(inum)
               natommove = natommove + 1
               atomlist(natommove,CRD_OLD) = CRD_REAL
               atomlist(natommove,CRD_NEW) = CRD_TEMP
               atomlist(natommove,3) = currentunit
            enddo

            globalcord(CRD_OLD) = CRD_REAL
            globalcord(CRD_NEW) = CRD_TEMP
            call twh_ec_atomlist_bond(lfailure,ichain,ibox
     &           ,natommove,atomlist,globalcord,loverlap,dvvib)
            if ( lfailure ) return
            if ( loverlap(CRD_OLD) .or. loverlap(CRD_NEW) ) then
               write(6,*) 'CONROT: Overlap vibration'
               lfailure = .true.
               return
            endif

            do icompute = 1,units
               if ( dcrs(icompute) .ne. 0) then
c                 --- compute bond-bond interactions
                  do iend = 0,twh_get_invib(imolty,icompute)
                     if ( iend .eq. 0 ) then
                        iuone = icompute
                     else
                        call twh_ijvib(GLB_GET,imolty,icompute,iend
     &                       ,iuone)
                     endif
                     do ibend = 1,twh_get_inbnbn(imolty,iuone)
                        iutwo = twh_get_ijbnbn1(imolty,iuone,ibend)
                        iuthree = twh_get_ijbnbn2(imolty,iuone,ibend)
                        if ( icompute .eq. iuone 
     &                       .or. icompute .eq. iutwo
     &                       .or. icompute .eq. iuthree ) then
                           do ioldnew = 1,2
                              if ( ioldnew .eq. CRD_OLD ) then
                                 testcord = CRD_REAL
                              else
                                 testcord = CRD_TEMP
                              endif
c                             --- compute and store first bond length
                              call twh_coordinates(GLB_GET,testcord
     &                             ,ichain,iuone,rxui,ryui,rzui)
                              call twh_coordinates(GLB_GET,testcord
     &                             ,ichain,iutwo,rxuj,ryuj,rzuj)
                              call twh_mimage(ibox,rxui,ryui,rzui
     &                             ,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij)
                              blength(1)=twh_distance(rxuij,ryuij,rzuij)
c                             --- compute and store second bond length
                              call twh_coordinates(GLB_GET,testcord
     &                             ,ichain,iuthree,rxuj,ryuj,rzuj)
                              call twh_mimage(ibox,rxui,ryui,rzui
     &                             ,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij)
                              blength(2)=twh_distance(rxuij,ryuij,rzuij)
                              dvvib = dvvib + dvpm(ioldnew)*
     &                             twh_vbonbon(imolty,iuone,ibend,
     &                             blength)
                           enddo
                        endif
                     enddo
                  enddo

c                 --- compute bond bending interactions
                  iuone = icompute
                  do ibend = 1,twh_get_inben(imolty,iuone)
                     iutwo = twh_get_ijben2(imolty,iuone,ibend)
                     iuthree = twh_get_ijben3(imolty,iuone,ibend)
                     if ( (dcrs(iuthree) .eq. 0) .or.
     &                    (dcrs(iuthree) .ne. 0 .and.iuone .lt. iuthree)
     &                    ) then
c                       --- compute this bond angle
                        do ioldnew = 1,2
                           if ( ioldnew .eq. CRD_OLD ) then
                              testcord = CRD_REAL
                           else
                              testcord = CRD_TEMP
                           endif
                           loverlap(ioldnew) = .false.
c                          --- store coordinates
                           call twh_coordinates(GLB_GET,testcord,ichain
     &                          ,iuone,xtmp(1),ytmp(1),ztmp(1))
                           call twh_coordinates(GLB_GET,testcord,ichain
     &                          ,iutwo,xtmp(2),ytmp(2),ztmp(2))
                           call twh_coordinates(GLB_GET,testcord,ichain
     &                          ,iuthree,xtmp(3),ytmp(3),ztmp(3))
                           dvbend = dvbend + dvpm(ioldnew)*
     &                          twh_vangle(lfailure,loverlap(ioldnew)
     &                          ,ibox,imolty,iuone,ibend,1
     &                          ,duma,dumb,dumc,xtmp,ytmp,ztmp)
                           if ( lfailure ) return
                           if ( loverlap(ioldnew) ) then
                              write(6,*) 'CONROT: Overlap bending'
                              lfailure = .true.
                              return
                           endif
                        enddo
                     endif
                  enddo

c                 --- compute the angle-angle energies
                  do iend = 0,twh_get_invib(imolty,icompute)
                     if ( iend .eq. 0 ) then
                        iuone = icompute
                     else
                        call twh_ijvib(GLB_GET,imolty,icompute,iend
     &                       ,iuone)
                     endif
                     do iaa = 1,twh_get_inaa(imolty,iuone)
                        torbead(1) = iuone
                        torbead(2) = twh_get_ijaa0(imolty,iuone,iaa)
                        torbead(3) = twh_get_ijaa1(imolty,iuone,iaa)
                        torbead(4) = twh_get_ijaa2(imolty,iuone,iaa)
                        lfound = .false.
                        do ibead = 1,4
                           if (torbead(ibead).eq.icompute) 
     &                          lfound = .true.
                        enddo
                        if ( lfound ) then
                           do ioldnew = 1,2
                              if ( ioldnew .eq. CRD_OLD ) then
                                 testcord = CRD_REAL
                              else
                                 testcord = CRD_TEMP
                              endif
c                             --- assign positions
                              do ibead = 1,4
                                 itest = torbead(ibead)
                                 call twh_coordinates(GLB_GET,testcord
     &                                ,ichain,itest,xtmp(ibead)
     &                                ,ytmp(ibead),ztmp(ibead))
                              enddo
c                             --- compute angle-angle energy   
                              dvbend = dvbend + dvpm(ioldnew)*
     &                             twh_vangang(AA_COORD,ibox,imolty
     &                             ,iuone,iaa,ddum2,xtmp,ytmp,ztmp)
                           enddo
                        endif
                     enddo
                  enddo
           
c                 --- compute the normal torsion energies for all of the
c                 --- backbone atoms.  The side chains should all have 
c                 --- the same torsions and bond angles at least as of 
c                 --- now, may change (CDL 5/21/03)
                  do itor = 1,twh_get_intor(imolty,icompute)
c                    --- assign the unit types for the torsion
                     torbead(1) = icompute
                     torbead(2) = twh_get_ijtor2(imolty,icompute,itor)
                     torbead(3) = twh_get_ijtor3(imolty,icompute,itor)
                     torbead(4) = twh_get_ijtor4(imolty,icompute,itor)
c                    --- check to make sure that not double counting any
c                    --- torsion whose other end is unit1 which has 
c                    --- already been calculated               
                     if ( (dcrs(torbead(4)) .eq. 0) .or.
     &                    (dcrs(torbead(4)) .ne. 0 
     &                    .and. torbead(1) .lt. torbead(4))
     &                    ) then
                        lcompute = .true.
                     else
                        lcompute = .false.
                     endif
                     if ( lcompute ) then
                        do ioldnew = 1,2
                           if ( ioldnew .eq. CRD_OLD ) then
                              testcord = CRD_REAL
                           else
                              testcord = CRD_TEMP
                           endif
c                          --- determine atom positions
                           do ibead = 1,4
                              itest = torbead(ibead)
                              call twh_coordinates(GLB_GET,testcord
     &                             ,ichain,itest,xtmp(ibead),ytmp(ibead)
     &                             ,ztmp(ibead))
                           enddo
c                          --- compute torsion energy
                           dvtors = dvtors + dvpm(ioldnew)*
     &                          twh_vtorsion(lfailure,ibox,1,imolty
     &                          ,icompute,itor,xtmp,ytmp,ztmp,ddum
     &                          ,loverlap(ioldnew))
                           if ( lfailure ) return
                           if ( loverlap(ioldnew) ) then
                              if (ioldnew .eq. 2) then
                                 goto 500
                              else
                                 write(6,*) 'CONROT: Overlap torsion'
     &                                ,' old conformation'
                                 lfailure = .true.
                                 return
                              endif
                           endif
                        enddo
                     endif
                  enddo
              
c                 --- compute the improper torsion energies for all of 
c                 --- the backbone atoms affected by the move.  The 
c                 --- sidechains should not be affected as they are 
c                 --- placed at same bond and torsion angles as 
c                 --- previously found just rotated around the axis.
                  torbead(1) = icompute
                  do itor = 1, twh_get_inimprop(imolty,icompute)
                     torbead(2) = twh_get_ijimprop2(imolty,icompute
     &                    ,itor)
                     torbead(3) = twh_get_ijimprop3(imolty,icompute
     &                    ,itor)
                     torbead(4) = twh_get_ijimprop4(imolty,icompute
     &                    ,itor)
                     do ioldnew = 1,2
                        if ( ioldnew .eq. CRD_OLD ) then
                           testcord = CRD_REAL
                        else
                           testcord = CRD_TEMP
                        endif
c                       --- assign atom positions
                        do ibead = 1,4
                           itest = torbead(ibead)
                           call twh_coordinates(GLB_GET,testcord
     &                          ,ichain,itest,xtmp(ibead),ytmp(ibead)
     &                          ,ztmp(ibead))
                        enddo
                        dvimpr = dvimpr + dvpm(ioldnew)*
     &                       twh_vimproper(ibox,imolty,icompute,itor,
     &                       xtmp,ytmp,ztmp,ddum,lbadstereo)
                     enddo
                  enddo
               endif
            enddo
c           --- calculate total energies of the initial and trial states
c           --- (sum of inter, intra, torsional, coulombic)
            dvtot = dvtot + dvtors + dvimpr + dvvib + dvbend
            if (nsoln2 .gt. 0) then
               call twh_beta(GLB_GET,beta)
               cr_accept = ((1.0d0/(dble(nsoln2)))
     &              *twh_expon(-beta*dvtot)
     &              *newdeter)/((1.0d0/(dble(nsoln)))*olddeter)
            else 
               cr_accept = 0.0d0
            endif

            if (cr_accept .gt. 1.0d0) then
               cr_accept = 1.0d0
            endif
c           --- test for acceptance
            laccept = .false.
            if( twh_random() .lt. cr_accept) then
               laccept = .true.
            endif
#if DEBUG_CONROT
            write(6,*) 'laccept:',laccept
            write(6,*) 'dvtot:',dvtot
            write(6,*) 'dvinter:',dvinter
            write(6,*) 'dvintra:',dvintra
            write(6,*) 'dvext:',dvext
            write(6,*) 'dvsolv:',dvsolv
            write(6,*) 'dvelect:',dvelect
            write(6,*) 'dvimpr:',dvimpr
            write(6,*) 'dvtors:',dvtors
            write(6,*) 'dvbend:',dvbend
            write(6,*) 'dvvib:',dvvib
c           --- force acceptance
            laccept = .true.
#endif
            if (laccept) then
c              --- update count of accepted CONCERTED ROTATION moves
               if (crtype .eq. 0) then
                  call twh_bacrback(GLB_INCR,imolty,unit1,ibox,1)
               elseif (crtype .eq. 1) then
                  call twh_baconrot(GLB_INCR,imolty,unit1,ibox,1)
               endif
c              --- update energies
               call twh_boxvclassic(GLB_INCR,ibox,dvtot)
               call twh_boxvinter(GLB_INCR,ibox,dvinter)
               call twh_boxvintra(GLB_INCR,ibox,dvintra)
               call twh_boxvfield(GLB_INCR,ibox,dvext)
               call twh_boxvelect(GLB_INCR,ibox,dvelect)
               dvdihed = dvtors + dvimpr
               call twh_boxvtorsion(GLB_INCR,ibox,dvdihed)
               call twh_boxvsolvation(GLB_INCR,ibox,dvsolv)
               call twh_boxvbond(GLB_INCR,ibox,dvvib)
               call twh_boxvangle(GLB_INCR,ibox,dvbend)

               call twh_scp_eval_derivative(GLB_GET,scp_eval_derivative)
               if (scp_eval_derivative) then
                  call scp_increment_derivative(ibox, LAMBDA_C, 
     &                  dlelect)
                  call scp_increment_derivative(ibox, LAMBDA_LJ, 
     &                  dlinter)
               endif
               call scp_increment_foreign_energies(ibox, dvf)

c              --- assign new atom positions
               do iunit = 1, units
                  call twh_coordinates(GLB_GET,CRD_TEMP,ichain,iunit
     &                 ,xtemp,ytemp,ztemp)
                  call twh_coordinates(GLB_SET,CRD_REAL,ichain,iunit
     &                 ,xtemp,ytemp,ztemp)
               enddo
               if ( lewald ) then
c                 --- update the reciprical space sum
                  call twh_recip(ibox,2)
               endif
c              --- update the center of mass of chain
               call twh_ctrmas(lfailure,CRD_REAL,ibox,ichain,CTR_CONROT)
               if ( lfailure ) return
c              --- update embedding density
               iboxa(1)=ibox
               call twh_reseteamrho(1,iboxa)
            endif
         endif
      else
         write(6,*) 'unit1',unit1
         write(6,*) 'error: no neighbor found for concerted rotation'
         lfailure = .true.
         return
      endif

#if DEBUG_CONROT
      write (6,*) 'End CONROT Move'
#endif

 500  continue
      return
      end

      subroutine twh_crstart(crtype,imolty,startid,dcrs,crback,isuccess)
c     ******************************************************************
c     * this subroutine finds a possible chain of 11 atoms which are   *
c     * required to do the concerted rotation move                     *
c     *                                                                *
c     * originally written 05-01-2003 by C.D. Lorenz                   *
c     * last modified in 10-13-2008 by M.G. Martin                     *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer crtype,imolty,startid,dcrs,crback,isuccess
      dimension dcrs(NUMAX),crback(40)

c     --- local variables
      integer unit1,searchsize,inbackbone,iiunit,numunits

      logical lfound

      lfound = .false.
      isuccess = -1
      unit1 = startid
c     --- find all units down stream of unit2
      searchsize = 11
      call twh_nunit(GLB_GET,imolty,numunits)
      do iiunit = 1,numunits
         dcrs(iiunit) = 0
      enddo
      do inbackbone = 1,40
         crback(inbackbone) = 0
      enddo
c     --- crback will hold the unit id of all of the backbone atoms that
c     --- make the section of the chain that will undergo the concerted
c     --- rotation move
c     ---  crback(1)crback(2),crback(3) represent the first root of the
c     --- move, and crback(searchsize-2),crback(searchsize-1),
c     --- crback(searchsize) represent the second root
      crback(1) = unit1
      dcrs(unit1) = 1
      call twh_crsearch(unit1,imolty,dcrs,searchsize,crback,unit1,
     &     isuccess,crtype)
      return
      end


      subroutine twh_crsearch(checkunit,imolty,dcrs,searchsize,crback
     &     ,unit3,isuccess,crtype)
c     ******************************************************************
c     * this part of the subroutine conducts a recursive depth first   *
c     * search starting with unit2 which is the second molecule        *
c     * chosen as part of the torsion arm.  This search is built       *
c     * such that it will exclude unit2 and the other part of the      *
c     * torsion arm unit1 and just gets everything 'downstream' of     *
c     * unit2.                                                         *
c     * originally written 05-01-2003 by C.D. Lorenz                   *
c     * last modified in 10-13-2008 by M.G. Martin                     *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_INVIB
#include "functions.h"

c     --- variables passed to/from the subroutine
      integer imolty,checkunit,unit3,searchsize
      integer crback,dcrs,isuccess,crtype
      dimension dcrs(NUMAX),crback(40)
c     --- local variables
      logical ldone,lbackbone
      integer currentunit,iunit,itest
      integer vibpointer,prevunit,visitcount,nunit
      dimension vibpointer(NUMAX),prevunit(NUMAX),visitcount(NUMAX)

      currentunit = checkunit
      call twh_nunit(GLB_GET,imolty,nunit)
      do iunit = 1,nunit
c        --- initialize the vibpointer
         vibpointer(iunit) = 0
c        --- initialize the prevunit
         prevunit(iunit) = 0
c        --- initialize visitcount
         visitcount(iunit) = 0
      enddo
      ldone = .false.

      if (crtype .eq. 0) then
c        --- this type of search only considers atoms designated as
c        --- belonging to the backbone
         call twh_lbackbone(GLB_GET,imolty,currentunit,lbackbone)
         if ( lbackbone ) then
            do while (.not. ldone)
               vibpointer(currentunit) = vibpointer(currentunit) + 1
c              --- we are finished exploring this atom, kick out
               if (vibpointer(currentunit) .gt. 
     &            twh_get_invib(imolty,currentunit)) then
                  if (prevunit(currentunit) .eq. 0) then
c                    --- completely finished with this current unit
                     ldone = .true.
                  else
c                    --- move back to previous atom to see if can find 
c                    --- necessary number of chain members
                     currentunit = prevunit(currentunit)
                  endif
               else
c                 --- see if you can continue along chain to get
c                 --- necessary number of atoms for conrot move
                  call twh_ijvib(GLB_GET,imolty,currentunit,
     &                 vibpointer(currentunit),itest)
                  call twh_lbackbone(GLB_GET,imolty,itest,lbackbone)
                  if ( lbackbone ) then
                     if (dcrs(itest) .eq. 0) then  
c                       --- neighbor has not been visited in this 
c                       --- search yet
                        dcrs(itest) = dcrs(currentunit) + 1
                        if (dcrs(itest) .eq. searchsize) then
c                          --- have found the desired number of atoms in
c                          --- chain
                           if (crback(searchsize) .eq. 0) then
                              ldone = .true.
                              isuccess = 1
                           endif
                        else
c                          --- need to continue search along same path
                           prevunit(itest) = currentunit
                           currentunit = itest
                           vibpointer(currentunit) = 0
                           visitcount(currentunit) = 
     &                         visitcount(currentunit) + 1
                        endif
                     endif 
                  endif
               endif
               if (visitcount(currentunit) .gt. NNBOND*NUMAX) then
c                 --- we must be stuck in some strange periodic ring
c                 --- structure.  Kick out and hope we figure it out
c                 --- from the next startunit
                  ldone = .true.
               endif
            enddo
         endif
      elseif (crtype .eq. 1) then
c        --- this type of search considers all atoms in the molecule
         do while (.not. ldone)
            vibpointer(currentunit) = vibpointer(currentunit) + 1
c           --- we are finished exploring this atom, kick out
            if (vibpointer(currentunit) .gt. 
     &         twh_get_invib(imolty,currentunit)) then
               if (prevunit(currentunit) .eq. 0) then
c                 --- completely finished with this current unit
                  ldone = .true.
               else
c                 --- move back to previous atom to see if can find 
c                 --- necessary number of chain members
                  currentunit = prevunit(currentunit)
               endif
            else
c              --- see if you can continue along chain to get
c              --- necessary number of atoms for conrot move
               call twh_ijvib(GLB_GET,imolty,currentunit,
     &              vibpointer(currentunit),itest)
               if (dcrs(itest) .eq. 0) then  
c                 --- neighbor has not been visited in this 
c                 --- search yet
                  dcrs(itest) = dcrs(currentunit) + 1
                  if (dcrs(itest) .eq. searchsize) then
c                    --- have found the desired number of atoms in
c                    --- chain
                     if (crback(searchsize) .eq. 0) then
                        ldone = .true.
                        isuccess = 1
                     endif
                  else
c                    --- need to continue search along same path
                     prevunit(itest) = currentunit
                     currentunit = itest
                     vibpointer(currentunit) = 0
                     visitcount(currentunit) = 
     &                   visitcount(currentunit) + 1
                  endif
               endif 
            endif
            if (visitcount(currentunit) .gt. NNBOND*NUMAX) then
c              --- we must be stuck in some strange periodic ring
c              --- structure.  Kick out and hope we figure it out from
c              --- the next startunit
               ldone = .true.
            endif
         enddo
      endif
 
      if (isuccess .eq. 1) then
c        --- if found a long enough path for move then call trace
c        --- which labels the atoms for the move
         call twh_crtrace(itest,unit3,dcrs,crback,imolty)
      endif

      return
      end



      subroutine twh_crtrace(checkunit,initunit,dcrs,crback,imolty)
c     *****************************************************************
c     * this part of the subroutine labels all of the units in        *
c     * given segment, uses id of first unit found to be part of      *
c     * the segment as the segment label                              *
c     * originally written 05-01-2003 by C.D. Lorenz                  *
c     * last modified in 10-13-2008 by M.G. Martin                    *
c     *****************************************************************
      implicit none
#define FUNCTION_GET_INVIB
#include "functions.h"

c     --- variables passed to and from the subroutine
      integer checkunit,initunit,imolty
      integer dcrs,crback
      dimension dcrs(NUMAX),crback(40)

c     --- local variables
      integer index1,index2,unitnext,itrace,iback,nnunit,crindex

      index1 = dcrs(checkunit)
      index2 = dcrs(initunit)
      crback(index1) = checkunit
      unitnext = checkunit
      
c     --- backtrace from index1 to index2
      do itrace = index1,index2+1,-1
         do iback = 1,twh_get_invib(imolty,unitnext)
            call twh_ijvib(GLB_GET,imolty,unitnext,iback,nnunit)
            if (nnunit .gt. 0) then
c              --- check to make sure neighbor is previous unit in this 
c              --- search
               if (dcrs(nnunit) .eq. (dcrs(unitnext)-1)) then
                  if (dcrs(nnunit) .gt. 0) then
                     crindex = dcrs(nnunit)
                     crback(crindex) = nnunit
                     unitnext = nnunit                  
                  endif
               endif
            endif
         enddo
      enddo 

      return
      end

      subroutine cr_angle2cart(ichain,atom3,atom1,atom0,atom2)
c     ******************************************************************
c     * originally written 05-01-2003 by C.D. Lorenz                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ichain,atom1,atom2,atom3,atom0
c     --- local variables
c     --- integer scalars
      integer i,j,ioldnew,testcord
c     --- double precision scalars
      double precision crossx,crossy,crossz,dist_2,dist_10,dist_20
c     --- double precision arrays
      double precision vec10,vec1,vec20,vec2,vec3,vec30
      dimension vec10(3),vec1(3,2),vec20(3),vec2(3,2),vec3(3,2),vec30(3)
      double precision trans,coord_a3
      dimension trans(3,3),coord_a3(3)
      double precision xtmp,ytmp,ztmp
      dimension xtmp(0:3),ytmp(0:3),ztmp(0:3)

c     --- determine the body frame coordinates at the 
c     --- original position and the trial position
      do ioldnew = CRD_OLD,CRD_NEW
         if ( ioldnew .eq. CRD_OLD ) then
            testcord = CRD_REAL
         else
            testcord = CRD_TEMP
         endif
         call twh_coordinates(GLB_GET,testcord,ichain,atom0
     &        ,xtmp(0),ytmp(0),ztmp(0))
         call twh_coordinates(GLB_GET,testcord,ichain,atom1
     &        ,xtmp(1),ytmp(1),ztmp(1))
         call twh_coordinates(GLB_GET,testcord,ichain,atom2
     &        ,xtmp(2),ytmp(2),ztmp(2))

         vec10(1) = xtmp(1) - xtmp(0)
         vec10(2) = ytmp(1) - ytmp(0)
         vec10(3) = ztmp(1) - ztmp(0)

         dist_10 = twh_distance(vec10(1),vec10(2),vec10(3))
         do i = 1,3
            vec1(i,ioldnew) = vec10(i)/dist_10
         enddo

         vec20(1) = xtmp(2) - xtmp(0)
         vec20(2) = ytmp(2) - ytmp(0)
         vec20(3) = ztmp(2) - ztmp(0)

         dist_20 = twh_distance(vec20(1),vec20(2),vec20(3))

         call twh_crossproduct(vec10(1),vec10(2),vec10(3),vec20(1),
     &          vec20(2),vec20(3),crossx,crossy,crossz)
 
         dist_2 = twh_distance(crossx,crossy,crossz)

         vec2(1,ioldnew) = crossx/dist_2 
         vec2(2,ioldnew) = crossy/dist_2
         vec2(3,ioldnew) = crossz/dist_2

         call twh_crossproduct(vec1(1,ioldnew),vec1(2,ioldnew)
     &        ,vec1(3,ioldnew),vec2(1,ioldnew),vec2(2,ioldnew)
     &        ,vec2(3,ioldnew),vec3(1,ioldnew),vec3(2,ioldnew)
     &        ,vec3(3,ioldnew))
      enddo

c     --- determine the vector which connects side chain atom to local
c     --- origin in original confirmation (local origin taken to be the 
c     --- atom to which the side chain is connected in the backbone, 
c     --- the C_alpha in in a protein)

      call twh_coordinates(GLB_GET,CRD_REAL,ichain,atom0,xtmp(0)
     &     ,ytmp(0),ztmp(0))
      call twh_coordinates(GLB_GET,CRD_REAL,ichain,atom3,xtmp(3)
     &     ,ytmp(3),ztmp(3))
      vec30(1) = xtmp(3) - xtmp(0)
      vec30(2) = ytmp(3) - ytmp(0)
      vec30(3) = ztmp(3) - ztmp(0)
      
c     --- determine the "transformation matrix" from the old local 
c     --- coordinate system and the new local coordinate system
      do i = 1,3
         do j = 1,3
            trans(i,j) = (vec1(i,2)*vec1(j,1))+(vec2(i,2)*vec2(j,1))+
     &              (vec3(i,2)*vec3(j,1))
         enddo         
      enddo
c     --- determine the coordinates of the side chain atom, assuming 
c     --- that atom0 is the origin.
      do i = 1,3
         coord_a3(i) = 0.0d0
         do j = 1,3
            coord_a3(i) = (vec30(j)*trans(i,j))+coord_a3(i)
         enddo
      enddo

c     --- displace the side chain coordinates by the coordinates of 
c     --- atom0 in the trial configuration
      call twh_coordinates(GLB_GET,CRD_TEMP,ichain,atom0,xtmp(0),ytmp(0)
     &     ,ztmp(0))
      xtmp(3) = coord_a3(1) + xtmp(0)
      ytmp(3) = coord_a3(2) + ytmp(0)
      ztmp(3) = coord_a3(3) + ztmp(0)
      call twh_coordinates(GLB_SET,CRD_TEMP,ichain,atom3,xtmp(3),ytmp(3)
     &     ,ztmp(3))
 
      return
      end

      subroutine twh_calcjacob(ichain,crback,olddeter,newdeter)
c     ******************************************************************
c     * This subroutine calculates the determinants of the Jacobian    *
c     * matrix for the new and old configurations and returns them     *
c     * (newdeter and olddeter)                                        *
c     *                                                                *
c     * originally written 05-01-2003 by C.D. Lorenz                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ichain,crback
      dimension crback(40)
      double precision olddeter,newdeter
c     --- local variables
c     --- integer scalars
      integer i,atomi,atomj,j,ioldnew,testcord
c     --- double precision scalars
      double precision xij,yij,zij,distij,xp,yp,zp,deter
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj
c     --- double precision arrays
      double precision uvec,bmat
      dimension uvec(6,3),bmat(5,5)

c     --- initialize bmat
      do i = 1,5
         do j = 1,5
            bmat(i,j) = 0.0d0
         enddo
      enddo    

c     --- first do the determinant for the old (1) and new (2) 
c     --- configurations
      do ioldnew = CRD_OLD,CRD_NEW
         if ( ioldnew .eq. CRD_OLD ) then
            testcord = CRD_REAL
         else
            testcord = CRD_TEMP
         endif
c        --- calculate the unit vectors u3,u4,u5,u6,u7,u8 and store them
c        --- in uvec
         do i = 3,8
            atomi = crback(i)
            call twh_coordinates(GLB_GET,testcord,ichain,atomi,rxui
     &           ,ryui,rzui)
            atomj = crback(i-1)
            call twh_coordinates(GLB_GET,testcord,ichain,atomj,rxuj
     &           ,ryuj,rzuj)
            xij = rxui - rxuj
            yij = ryui - ryuj
            zij = rzui - rzuj
            distij = twh_distance(xij,yij,zij)
            uvec(i-2,1) = xij/distij
            uvec(i-2,2) = yij/distij
            uvec(i-2,3) = zij/distij
         enddo

c        --- calculate the elements of the B matrix
         atomi = crback(7)
         call twh_coordinates(GLB_GET,testcord,ichain,atomi
     &        ,rxui,ryui,rzui)
         do i = 3,6
            atomj = crback(i)
            call twh_coordinates(GLB_GET,testcord,ichain,atomj
     &           ,rxuj,ryuj,rzuj)
c           --- calculate the r7j's needed for the termsof the first 
c           --- three rows
            xij = rxui - rxuj
            yij = ryui - ryuj
            zij = rzui - rzuj
c           --- calculate the dot products used in the first three rows
            call twh_crossproduct(uvec(i-2,1),uvec(i-2,2),uvec(i-2,3),
     &                     xij,yij,zij,xp,yp,zp)      
            bmat(1,i-2) = xp
            bmat(2,i-2) = yp
            bmat(3,i-2) = zp
         enddo 

         do i = 3,7
c           --- calculate the terms of the fourth row
            call twh_crossproduct(uvec(i-2,1),uvec(i-2,2),uvec(i-2,3),
     &                     uvec(6,1),uvec(6,2),uvec(6,3),xp,yp,zp)
            bmat(4,i-2) = xp
            bmat(5,i-2) = yp
         enddo

         deter = 0.0d0
         call calcdeter(bmat,deter)
         if ( ioldnew .eq. CRD_OLD ) then
            olddeter = 1.0d0/dabs(deter)
         elseif ( ioldnew .eq. 2) then
            newdeter = 1.0d0/dabs(deter)
         endif
      enddo

      return
      end

      subroutine calcdeter(bmat,deter)
c     ******************************************************************
c     * This subroutine calculates the determinant of a 5x5 matrix     *
c     *                                                                *
c     * originally written 05-01-2003 by C.D. Lorenz                   *
c     * last modified in 10-13-2008 by M.G. Martin                     *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision deter,bmat
      dimension bmat(5,5)

c     --- local variables
      integer i,j,k,ii,jj,ninverse
      double precision term 


      do i = 1,5
       do j = 1,5
        do k = 1,5
         do ii = 1,5
          do jj = 1,5
           if (j .ne. i) then
            if (k .ne. j .and. k .ne. i) then
             if (ii .ne. k .and. ii .ne. j .and. ii .ne. i) then
              if (jj .ne. ii .and. jj .ne. k .and. jj .ne. j .and.
     &            jj .ne. i) then
                 ninverse = 0
                 if (i .gt. j) then
                    ninverse = ninverse+1
                 endif
                 if (i .gt. k) then
                    ninverse = ninverse+1
                 endif
                 if (i .gt. ii) then
                    ninverse = ninverse+1
                 endif
                 if (i .gt. jj) then
                    ninverse = ninverse+1
                 endif
                 if (j .gt. k) then
                    ninverse = ninverse+1
                 endif
                 if (j .gt. ii) then
                    ninverse = ninverse+1
                 endif
                 if (j .gt. jj) then
                    ninverse = ninverse+1
                 endif
                 if (k .gt. ii) then
                    ninverse = ninverse+1
                 endif
                 if (k .gt. jj) then
                    ninverse = ninverse+1
                 endif
                 if (ii .gt. jj) then
                    ninverse = ninverse+1
                 endif
                 term = (-1.0d0**ninverse)*bmat(1,i)*bmat(2,j)
     &                 *bmat(3,k)*bmat(4,ii)*bmat(5,jj)
                 deter = deter + term 
              endif
             endif
            endif
           endif
          enddo
         enddo
        enddo
       enddo
      enddo
 
      return
      end
