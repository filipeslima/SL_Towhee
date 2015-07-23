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
      subroutine twh_pivot(lfailure)
c     ******************************************************************
c     * makes a pivot, torsional, movement at some randomly chosen     *
c     * bond whose end points are x1 and x2. x1 is the fixed point     *
c     * while x2 and everything downstream from it are rotated some    *
c     * angle between 0 and 2 pi.  the number of successful pivots are *
c     * stored in bapivot(molty,nunit,nvib,nbox).                      *
c     *                                                                *
c     * originally written by M.G. Martin and C.D. Lorenz 10-16-2002   *
c     * last modifed in 08-09-2011 by M.G. Martin                      *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_EXPON
#define FUNCTION_GET_IJVIB
#define FUNCTION_GET_INVIB
#define FUNCTION_RANDOM
#include "functions.h"
#define DEBUG_PIVOT 0
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
c     --- logical scalars
      logical lfound,lgetcom,laccept,lpartial,lintra,lbonded
      logical scp_eval_derivative,lewald
c     --- logical arrays
      logical loverlap
      dimension loverlap(2)
c     --- integer scalars
      integer icmol,imolty,iunit,ichain,activebox,units,iu
      integer unit1,unit1ptr,unit1nno,unit2,junit,bintheta,binsize
      integer nchainmove,natommove,nmolty,ncmt,ipoint,nmtype
      integer inringone,inringtwo,ensemble
c     --- integer arrays
      integer ds,ntotry,ivalid,validlist
      dimension ds(NUMAX)
      dimension validlist(NNBOND)
      integer chaininfo
      dimension chaininfo(1,3)
      integer atomlist
      dimension atomlist(2,3)
      integer globalcord
      dimension globalcord(2)
      integer iboxa
      dimension iboxa(1)
c     --- double precision scalars
      double precision ranmolty,theta,dchain,beta
      double precision distx,disty,distz,vecx,vecy,vecz
      double precision dvdumb,dvdumc,xtemp,ytemp,ztemp
      double precision dist,thetap,burn
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj
      double precision ddvdl_lj,ddvdl_c,onepi,pmtest
c     --- double precision arrays
      double precision dvtot,dvintra,dvinter,dvext,dvelect,dvsolv
      double precision dvtors,dvquantum,dvatomlist,twopi
      double precision dvf
      dimension dvf(MAX_FOREIGN_LAMBDA)

#if DEBUG_PIVOT
      write(6,*) 'start PIVOT'
#endif
c     --- retrieve constants
      call twh_constant_pi(GLB_GET,onepi)
      call twh_constant_twopi(GLB_GET,twopi)
      call twh_nmolty(GLB_GET,nmolty)
      call twh_lewald(GLB_GET,lewald)

c     --- select a molecule type at random 
      ranmolty  = twh_random()
      do icmol = 1,nmolty
         call twh_pmpivmt(GLB_GET,icmol,pmtest)
         if ( ranmolty .lt. pmtest ) then
            imolty = icmol
            ranmolty = 2.0d0
         endif
      enddo

c     --- select a chain of the determined molecule type
      call twh_ensemble(GLB_GET,ensemble)
      if ( ensemble .eq. ENS_muVT ) then
c        --- take from box 1
         activebox = 1
c        --- see if there are any molecules of this type in active box
         call twh_ncmt(GLB_GET,activebox,imolty,ncmt)
         if ( ncmt .eq. 0 ) return
c        --- select a molecule of type imolty in box activebox
         ipoint = idint( dble(ncmt)*twh_random()) + 1
         call twh_chainlist(GLB_GET,ipoint,activebox,imolty,ichain)
      else
         call twh_nmolectyp(GLB_GET,imolty,nmtype)
         dchain = dble(nmtype)
         ipoint = idint((dchain*twh_random()) + 1)
         call twh_parall(GLB_GET,imolty,ipoint,ichain)
      endif

#if DEBUG_PIVOT
      write(6,*) 'ichain:',ichain
#endif

c     --- determine what box this part of the molecule is in
      call twh_nboxi(GLB_GET,ichain,activebox)

c     --- store the number of units in this molecule
      call twh_nunit(GLB_GET,imolty,units)

      do iu = 1, units
         call twh_coordinates(GLB_GET,CRD_REAL,ichain,iu
     &        ,xtemp,ytemp,ztemp)
         call twh_coordinates(GLB_SET,CRD_TEMP,ichain,iu
     &        ,xtemp,ytemp,ztemp)
      enddo

c     --- initialize ds array for use in finding a valid pair of atoms
c     --- for this move.
c     --- NOTE: this array has a different meaning later in this 
c     --- subroutine.  It is reinitialized after we find a valid target
      ntotry = units
      do iunit = 1, units
         ds(iunit) = iunit
      enddo      

c     --- select two units from the chosen molecule, where at most 
c     --- one of the two units is in a ring, as defined in findrings
      lfound = .false.
      do while ( .not. lfound .and. ntotry .gt. 0 )
         unit1ptr = idint(( dble(ntotry)*twh_random())+1)
         unit1 = ds(unit1ptr)
#if DEBUG_PIVOT
         write(6,*) 'ntotry,unit1ptr,unit1:',ntotry,unit1ptr,unit1
#endif
c        --- select bonded neighbor (unit2) from unit1.
c        --- everything downstream from unit2 will be pivoted.
c        --- make a list of all the atoms bonded to unit1 that would
c        --- be valid
         ivalid = 0
c        --- if the first atom has more than 1 bond then check for 
c        --- second atom, otherwise this is not valid
#if DEBUG_PIVOT
         write(6,*) 'unit1,invib:',unit1,twh_get_invib(imolty,unit1)
#endif
         if ( twh_get_invib(imolty,unit1) .gt. 1 ) then
            do unit1nno = 1, twh_get_invib(imolty,unit1)
               unit2 = twh_get_ijvib(imolty,unit1,unit1nno)
#if DEBUG_PIVOT
               write(6,*) 'unit1nno,unit2,invib(unit2)',unit1nno
     &             ,unit2,twh_get_invib(imolty,unit2)
#endif
c              --- make sure this atom is bonded to something else
c              --- as otherwise this is a meaningless pivot
               if ( twh_get_invib(imolty,unit2) .gt. 1 ) then
c                 --- both atoms are potentially valid for pivot
c                 --- make sure that if they are in rings then they are
c                 --- not in the same ring
                  call twh_inring(GLB_GET,imolty,unit1,inringone)
                  call twh_inring(GLB_GET,imolty,unit2,inringtwo)
#if DEBUG_PIVOT
                  write(6,*) 'inring(unit1),inring(unit2):'
     &               ,inringone,inringtwo
#endif
                  if ( inringone .eq. 0 .or.
     &                 inringone .ne. inringtwo ) then 
                     ivalid = ivalid + 1
                     validlist(ivalid) = unit1nno
                  endif
               endif
            enddo
         endif
#if DEBUG_PIVOT
         write(6,*) 'ivalid:',ivalid
#endif
         if ( ivalid .gt. 0 ) then
c           --- select one of the valid targets at random
            unit1nno = idint( dble(ivalid)*twh_random() ) + 1
            unit1nno = validlist(unit1nno)
            unit2 = twh_get_ijvib(imolty,unit1,unit1nno)
            lfound = .true.
         else
c           --- no valid targets for this unit, remove it from the 
c           --- list and try again
            ds(unit1ptr) = ds(ntotry)
            ntotry = ntotry - 1
         endif
      enddo  
      if ( .not. lfound ) then
         write(6,*) 'PIVOT: no two valid neighbors for pivot found'
         write(6,*) 'molecule type ',imolty
         lfailure = .true.
         return
      endif

#if DEBUG_PIVOT
      write(6,*) 'unit1,unit2:',unit1,unit2
#endif

c     --- reinitialize ds array which is now used to keep track of atoms
c     --- downstream of unit2 
      do iunit = 1, units
         ds(iunit) = 0
      enddo
c     --- update number of PIVOT move trys
      call twh_bnpivot(GLB_INCR,imolty,unit1,unit1nno,activebox,1)
c     --- select an angle to pivot the chain
      theta = twopi * twh_random()
      thetap = theta
      if (theta .gt. onepi) then
         thetap = thetap - onepi
      endif
      thetap = thetap * 180d0/onepi
      binsize = int(180.0d0/dble(MAXPIVOTBIN))
      bintheta = int(thetap/binsize) + 1
      call twh_bnpivbin(GLB_INCR,activebox,bintheta,1)

#if DEBUG_PIVOT
      write(6,*) 'theta:',theta
      write(6,*) 'thetap:',thetap
#endif
         
c     --- find all units down stream of unit2
      ds(unit1) = 1
      ds(unit2) = 2
      call moleculeremains(unit2,imolty,ds)
c     --- remove unit1 and unit2 from the list as units to be moved
c     --- (which will be all units with ds > 0
      ds(unit1) = 0
      ds(unit2) = 0
c     --- determine the vector formed by the unit1-unit2 bond
c     --- this vector will be used to define the rotation axis
      call twh_coordinates(GLB_GET,CRD_REAL,ichain,unit1,rxui,ryui,rzui)
      call twh_coordinates(GLB_GET,CRD_REAL,ichain,unit2,rxuj,ryuj,rzuj)
      call twh_mimage(activebox,rxuj,ryuj,rzuj,rxui,ryui,rzui
     &     ,distx,disty,distz)
      dist = twh_distance(distx,disty,distz)
      vecx = distx/dist
      vecy = disty/dist
      vecz = distz/dist

#if DEBUG_PIVOT
      write(6,*) 'dist,vec*:',dist,vecx,vecy,vecz
#endif
      do junit = 1,units 
         if (ds(junit) .gt. 0) then 
c           --- call pivrot in order to get new possible positions
c           --- for atom after it has been rotated 'theta' degrees 
c           --- around the desired unit vector (vecx,vecy,vecz)  
            call twh_pivrot(ichain,junit,vecx,vecy,vecz,theta,unit2)
         endif
      enddo
         
c     --- calculate the new and old non-bonded inter molecular 
c     --- energies 
      nchainmove = 1
      chaininfo(1,CRD_OLD) = CRD_REAL
      chaininfo(1,CRD_NEW) = CRD_TEMP
      chaininfo(1,3) = ichain
      lpartial = .false.
      lintra = .true.
      lbonded = .false.
      lgetcom = .true.
#if DEBUG_PIVOT
      write(6,*) 'PIVOT: calling energy change',lfailure,loverlap
#endif
      call twh_energy_change(lfailure,activebox,nchainmove,chaininfo
     &     ,loverlap,lpartial,lgetcom,lintra,lbonded,dvtot
     &     ,dvintra,dvinter,dvext,dvelect,dvsolv,dvtors,dvdumb,dvdumc
     &     ,dvquantum,dvf,ddvdl_lj,ddvdl_c)
      if ( lfailure ) return
#if DEBUG_PIVOT
      write(6,*) 'PIVOT: back from energy change',loverlap
#endif

c     --- safety checks for hard overlap
      if ( loverlap(2) ) return
      if ( loverlap(1) ) then
         write(6,*) 'Overlap in old conformation in PIVOT'
         lfailure = .true.
         return
      endif

#if DEBUG_PIVOT
      write(6,*) 'twh_pivot() energy terms: dvtot: ',dvtot
      write(6,*) 'dvintra, dvinter, dvext: ', dvintra, dvinter, dvext
      write(6,*) 'idvelect, dvsolv: ', dvelect, dvsolv
      write(6,*) 'dvtors, dvdumb, dvdumc, dvquantum: ', dvtors, dvdumb, 
     &    dvdumc,dvquantum
#endif

c     --- compute normal torsion energies that include unit1 or unit2
      natommove = 2
      atomlist(1,CRD_OLD) = CRD_REAL
      atomlist(1,CRD_NEW) = CRD_TEMP
      atomlist(1,3) = unit1
      atomlist(2,CRD_OLD) = CRD_REAL
      atomlist(2,CRD_NEW) = CRD_TEMP
      atomlist(2,3) = unit2
      globalcord(CRD_OLD) = CRD_REAL
      globalcord(CRD_NEW) = CRD_TEMP
      dvatomlist = 0.0d0
      call twh_ec_atomlist_torsion(lfailure,ichain,activebox,natommove
     &     ,atomlist,globalcord,loverlap,dvatomlist)
      if ( lfailure ) return
c     --- check for hard overlap
      if ( loverlap(2) ) return
      if ( loverlap(1) ) then
         write(6,*) 'PIVOT: Overlap in old conformation'
         lfailure = .true.
         return
      endif
c     --- add this energy to the torsion and classical energies
#if DEBUG_PIVOT
      write(6,*) 'PIVOT: atomlist energy change:',dvatomlist
#endif
      dvtors = dvtors + dvatomlist
      dvtot = dvtot + dvatomlist

c     --- test for acceptance
      laccept = .false.

      if (dvtot .lt. 0.0d0) then
c        --- accept move if the energy goes down
         laccept = .true.
c        --- burn a random number
         burn = twh_random()
      else
         call twh_beta(GLB_GET,beta)
         if(twh_random() .lt. twh_expon(-beta*dvtot)) then
            laccept = .true.
         endif
      endif
#if DEBUG_PIVOT
      write(6,*) 'laccept:',laccept
      write(6,*) 'dvtot:',dvtot
      write(6,*) 'dvinter:',dvinter
      write(6,*) 'dvintra:',dvintra
      write(6,*) 'dvext:',dvext
      write(6,*) 'dvelect:',dvelect
      write(6,*) 'dvtors:',dvtors
      write(6,*) 'dvsolv:',dvsolv
      laccept = .true.
#endif
      
      if (laccept) then
c        --- update count of accepted PIVOT moves
         call twh_bapivot(GLB_INCR,imolty,unit1,unit1nno,activebox,1)
         call twh_bapivbin(GLB_INCR,activebox,bintheta,1)
c        --- update energies
         call twh_boxvclassic(GLB_INCR,activebox,dvtot)
         call twh_boxvinter(GLB_INCR,activebox,dvinter)
         call twh_boxvintra(GLB_INCR,activebox,dvintra)
         call twh_boxvfield(GLB_INCR,activebox,dvext)
         call twh_boxvelect(GLB_INCR,activebox,dvelect)
         call twh_boxvtorsion(GLB_INCR,activebox,dvtors)
         call twh_boxvsolvation(GLB_INCR,activebox,dvsolv)

         call twh_scp_eval_derivative(GLB_GET,scp_eval_derivative)
         if (scp_eval_derivative) then
            call scp_increment_derivative(activebox,LAMBDA_C,ddvdl_c)
            call scp_increment_derivative(activebox,LAMBDA_LJ,ddvdl_lj)
         endif
         call scp_increment_foreign_energies(activebox, dvf)

c        --- assign new atom positions
         do iunit = 1, units
            call twh_coordinates(GLB_GET,CRD_TEMP,ichain,iunit
     &           ,xtemp,ytemp,ztemp)
            call twh_coordinates(GLB_SET,CRD_REAL,ichain,iunit
     &           ,xtemp,ytemp,ztemp)
         enddo
         if ( lewald ) then
c           --- update the reciprical space sum
            call twh_recip(activebox,2)
         endif            
c        --- update the center of mass of chain
         call twh_ctrmas(lfailure,CRD_REAL,activebox,ichain,CTR_PIVOT)
         if ( lfailure ) return
c        --- update embedding density
         iboxa(1) = activebox
         call twh_reseteamrho(1,iboxa)
      endif

#if DEBUG_PIVOT
      write (6,*) 'End PIVOT Move'
#endif

      return
      end

      subroutine moleculeremains(checkunit,imolty,ds)
c     ******************************************************************
c     * rewritten so that it uses a while loop to find every atom that *
c     * is "downstream" of the checkunit                               *
c     *                                                                *
c     * originally written on 10-17-2002 by C.D. Lorenz                *
c     * last modified on 10-13-2008 by M.G. Martin                     *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJVIB
#define FUNCTION_GET_INVIB
#include "functions.h"

c     --- variables passed to/from the subroutine
      integer imolty,checkunit,ds
      dimension ds(NUMAX)

c     --- local variables
      integer icheck,junit,workleft,testunit
      integer worklist
      dimension worklist(NUMAX)

c     --- set the first element of the worklist
      workleft = 1
      worklist(workleft) = checkunit

      do while ( workleft .gt. 0 ) 
c        --- grab the most recent element on the worklist
         testunit = worklist(workleft)
#if DEBUG_PIVOT
         write(6,*) 'workleft,testunit:',workleft,testunit
#endif
c        --- decrease this from the worklist
         worklist(workleft) = 0
         workleft = workleft - 1
c        --- check all of the neighbors of the testunit and add them
c        --- unless they have already been added to the list
         do icheck = 1,twh_get_invib(imolty,testunit)
            junit =  twh_get_ijvib(imolty,testunit,icheck)
            if ( ds(junit) .eq. 0 ) then
c              --- add this to the list
               workleft = workleft + 1
               worklist(workleft) = junit
c              --- flip the integer switch so we know to rotate this one
               ds(junit) = 3
#if DEBUG_PIVOT
               write(6,*) 'found junit,workleft:',junit,workleft
#endif
            endif
         enddo
      enddo

      return
      end

      subroutine twh_pivrot(ichain,junit,vecx,vecy,vecz,theta,unit2)
c     ******************************************************************
c     * this part of the subroutine calculates the rotation matrix     *
c     * which will be used to deterimine the new locations of the units*
c     * of the chain which is being pivoted (matrix was found in an    *
c     * article entitled "The Mathematics of the 3D Rotation Matrix" by*
c     * Diana Gruber (found at http://www.makegames.com/3drotation/)   *
c     * referenced from Graphics Gems, Glassner, Academic Press, 1990. *
c     *                                                                *
c     * originally written on 10-17-2002 by C.D. Lorenz                *
c     * last modified on 10-13-2008 by M.G. Martin                     *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer junit,unit2,ichain
      double precision vecx,vecy,vecz,theta
c     --- local variables
      double precision oldx,oldy,oldz
      double precision c,s,t,r11,r12,r13,r21,r22,r23
      double precision r31,r32,r33
      double precision newx,newy,newz
      double precision rxj,ryj,rzj,rxu,ryu,rzu

      call twh_coordinates(GLB_GET,CRD_REAL,ichain,junit,rxj,ryj,rzj)
      call twh_coordinates(GLB_GET,CRD_REAL,ichain,unit2,rxu,ryu,rzu)

c     --- determine old coordinates of unit
      oldx = rxj - rxu
      oldy = ryj - ryu
      oldz = rzj - rzu
 
c     --- calculate common terms in rotation matrix
      c = cos(theta)
      s = sin(theta)
      t = 1 - cos(theta) 

c     --- calculate the terms in the rotation matrix
      r11 = (t * vecx * vecx) + c
      r12 = (t * vecx * vecy) - (s * vecz)
      r13 = (t * vecx * vecz) + (s * vecy)
      r21 = (t * vecx * vecy) + (s * vecz)
      r22 = (t * vecy * vecy) + c
      r23 = (t * vecy * vecz) - (s * vecx) 
      r31 = (t * vecx * vecz) - (s * vecy)
      r32 = (t * vecy * vecz) + (s * vecx)
      r33 = (t * vecz * vecz) + c

c     --- calculate the new x,y,z coordinates
      newx = (oldx * r11) + (oldy * r12) + (oldz * r13)
      newy = (oldx * r21) + (oldy * r22) + (oldz * r23)
      newz = (oldx * r31) + (oldy * r32) + (oldz * r33)

c     --- assign new x-, y-, z- coordinates
      rxj = newx + rxu
      ryj = newy + ryu
      rzj = newz + rzu
      call twh_coordinates(GLB_SET,CRD_TEMP,ichain,junit,rxj,ryj,rzj)

      return
      end

