#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2006-2011 Marcus G. Martin                       *
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
      subroutine twh_switch(lfailure)
c     ******************************************************************
c     * switches a pair of molecules at their center-of-mass positions *
c     *                                                                *
c     * completely rewritten 06-07-2008 by I.A. Hijazi                 *
c     * last modified 04-21-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_RANDOM
#define FUNCTION_EXPON
#define FUNCTION_MINBOXLENGTH
#include "functions.h"
c     --- debugging variable
#define DEBUG_SWITCH 0
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
c     --- local logical scalars
      logical lfound,laccept,lgetcom,lintra,lbonded,lpartial
      logical scp_eval_derivative,lewald
c     --- local logical arrays
      logical loverlap
      dimension loverlap(2)
c     --- local integer scalars
      integer ipair,npair,imolty,jmolty,ibox,nmolty,potentialstyle
      integer nimtype,njmtype,ipoint,jpoint,numboxes
      integer ichain,jchain,totiunit,totjunit,iunit,junit,nchainmove
c     --- local integer arrays
      integer chaininfo
      dimension chaininfo(2,3)
c     --- local double precision scalars
      double precision rselect,testprob,rix,riy,riz,rjx,rjy,rjz
      double precision xcmi,ycmi,zcmi,xcmj,ycmj,zcmj,dxcm,dycm,dzcm,rcmu
      double precision dvclassical_trial,dvintra_trial
      double precision dvinter_trial,dvext_trial,dvelect_trial
      double precision dvsolv_trial,dvtor_trial,dvbend_trial
      double precision dvvib_trial,dvquantum_trial,ddvdl_lj,ddvdl_c
      double precision dv_trial,beta,burn,dichain,djchain
c     --- local double precision arrays
      double precision dvf
      dimension dvf(MAX_FOREIGN_LAMBDA)
      integer iboxa
      dimension iboxa(1)

#if DEBUG_SWITCH
      write(6,*) 'start SWITCH'
#endif
c     --- retrieve constants
      call twh_potentialstyle(GLB_GET,potentialstyle)
      call twh_lewald(GLB_GET,lewald)
      call twh_numboxes(GLB_GET,numboxes)

c     --- select a box at random
      rselect = twh_random()
      lfound = .false.
      ibox = 0
      do while ( ibox .lt. numboxes .and. (.not. lfound) )
         ibox = ibox + 1
         call twh_pm1comswbox(GLB_GET,ibox,testprob)
         if ( rselect .lt. testprob ) lfound = .true.
      enddo
      if ( .not. lfound ) then
         write(6,*) 'SWITCH: problem finding a simulation box'
         lfailure = .true.
         return
      endif
#if DEBUG_SWITCH
      write(6,*) 'ibox:',ibox
#endif
c     --- select a molecule pair type at random
      rselect  = twh_random()
      call twh_nmolty(GLB_GET,nmolty)
      npair = (nmolty*(nmolty-1))/2
      ipair = 0
      imolty = 0
      jmolty = nmolty
      lfound = .false.
#if DEBUG_SWITCH
      write(6,*) 'nmolty,npair:',nmolty,npair
#endif
      do while ( ipair .lt. npair .and. (.not. lfound) )
         ipair = ipair + 1
         jmolty = jmolty + 1
         if ( jmolty .gt. nmolty ) then
c           --- increment imolty and reset jmolty
            imolty = imolty + 1
            jmolty = imolty + 1
         endif
         call twh_pm1comswpair(GLB_GET,ipair,testprob)
         if ( rselect .lt. testprob ) lfound = .true.
      enddo
      if ( .not. lfound ) then
         write(6,*) 'SWITCH: pair not found'
         lfailure = .true.
         return
      endif
#if DEBUG_SWITCH
      write(6,*) 'ipair:',ipair
      write(6,*) 'imolty,jmolty:',imolty,jmolty
#endif
c     --- select a chain of imolty type
      call twh_ncmt(GLB_GET,ibox,imolty,nimtype)
      if (nimtype .eq. 0) return
      dichain = dble(nimtype)
      ipoint = int( dichain*twh_random() + 1 )
      call twh_chainlist(GLB_GET,ipoint,ibox,imolty,ichain)
c     --- select a chain of jmolty type
      call twh_ncmt(GLB_GET,ibox,jmolty,njmtype)
      if (njmtype .eq. 0) return
      djchain = dble(njmtype)
      jpoint = int( djchain*twh_random() + 1 )
      call twh_chainlist(GLB_GET,jpoint,ibox,jmolty,jchain)
c     --- store number of units of i in totiunit 
      call twh_nunit(GLB_GET,imolty,totiunit)
c     --- store number of units of j in totjunit 
      call twh_nunit(GLB_GET,jmolty,totjunit)

#if DEBUG_SWITCH
      write(6,*) 'nimtype,njmtype',nimtype,njmtype
      write(6,*) 'dichain,djchain',dichain,djchain
      write(6,*) 'ipoint,jpoint',ipoint,jpoint
      write(6,*) 'totiunit,totjunit',totiunit,totjunit
      write(6,*) 'ichain,jchain',ichain,jchain
#endif

c  --- increment switch move attempt counter
      call twh_acnswitch(GLB_INCR,ipair,ibox,1.0d0)
c  --- switch centers of mass temporarilly 
      call twh_com(GLB_GET,CRD_REAL,ichain,xcmi,ycmi,zcmi)
      call twh_com(GLB_GET,CRD_REAL,jchain,xcmj,ycmj,zcmj)
      call twh_com(GLB_SET,CRD_TEMP,ichain,xcmj,ycmj,zcmj)
      call twh_com(GLB_SET,CRD_TEMP,jchain,xcmi,ycmi,zcmi)
      dxcm = xcmj - xcmi
      dycm = ycmj - ycmi
      dzcm = zcmj - zcmi
#if DEBUG_SWITCH
      write(6,*) 'cmi:',xcmi,ycmi,zcmi
      write(6,*) 'cmj:',xcmj,ycmj,zcmj
#endif
c  --- update chain coordinates temporarilly  
      do iunit = 1,totiunit
         call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit,rix,riy,riz)
         rix = rix + dxcm
         riy = riy + dycm
         riz = riz + dzcm
         call twh_coordinates(GLB_SET,CRD_TEMP,ichain,iunit,rix,riy,riz)
      enddo
      do junit = 1,totjunit
         call twh_coordinates(GLB_GET,CRD_REAL,jchain,junit,rjx,rjy,rjz)
         rjx = rjx - dxcm
         rjy = rjy - dycm
         rjz = rjz - dzcm
         call twh_coordinates(GLB_SET,CRD_TEMP,jchain,junit,rjx,rjy,rjz)
      enddo

      call twh_rcmu(GLB_GET,ichain,rcmu)
      if ( rcmu .gt. twh_minboxlength(ibox)*0.25d0 ) then
c        --- molecule is large enough that it might be interacting with
c        --- iself through the periodic images (occurs for some very 
c        --- short range potentials like hard spheres).   Need to 
c        --- compute the intramolecular and bonded interactions as these
c        --- could change in this case
         lintra  = .true.
         lbonded = .true.
      else
c        --- molecule does not interact with itself through the periodic
c        --- images.  don't bother computing intramolecular or bonded
         lintra  = .false.
         lbonded = .true.
      endif

      nchainmove = 2
      chaininfo(1,CRD_OLD) = CRD_REAL
      chaininfo(1,CRD_NEW) = CRD_TEMP
      chaininfo(1,3) = ichain
      chaininfo(2,CRD_OLD) = CRD_REAL
      chaininfo(2,CRD_NEW) = CRD_TEMP
      chaininfo(2,3) = jchain
      lpartial = .false.
      lgetcom  = .false.
      call twh_energy_change(lfailure,ibox,nchainmove,chaininfo,loverlap
     &     ,lpartial,lgetcom,lintra,lbonded,dvclassical_trial
     &     ,dvintra_trial,dvinter_trial,dvext_trial,dvelect_trial
     &     ,dvsolv_trial,dvtor_trial,dvbend_trial
     &     ,dvvib_trial,dvquantum_trial,dvf,ddvdl_lj,ddvdl_c)
      if ( lfailure ) return
#if DEBUG_SWITCH
      write(6,*) 'dvinter_trial:',dvinter_trial
#endif

c     --- check new hard overlap, if overlap move rejected
      if ( loverlap(2) ) return
c     --- check old hard overlap, if overlap then big problems
      if (loverlap(1)) then
         write(6,*) 'SWITCH: disaster loverlap in old conf'
         lfailure = .true.
         return
      endif

      if ( potentialstyle .eq. POT_INTERNAL ) then
         dv_trial = dvclassical_trial
      elseif ( potentialstyle .eq. POT_EXTERNAL ) then
         dv_trial = dvquantum_trial
      else
         write(6,*) 'SWITCH: unknown potentialstyle'
         write(6,*) potentialstyle
         stop
      endif

c     --- check for acceptance
      laccept = .false.
#if DEBUG_SWITCH
      write(6,*) 'dv_trial:',dv_trial
#endif
      if ( dv_trial .le. 0.0d0 ) then
c        --- accept move
         laccept = .true.
c        --- burn a random number
         burn = twh_random()
      else
         call twh_beta(GLB_GET,beta)
         if ( twh_random() .lt. twh_expon(-beta*dv_trial) ) then
            laccept = .true.
         endif
      endif

#if DEBUG_TRANCOM
      write(6,*) 'SWITCH: dv_trial ',dv_trial
      write(6,*) 'SWITCH: laccept:',laccept
#endif
      if ( .not. laccept ) return
#if DEBUG_TRANCOM
      write(6,*) 'SWITCH: accepted chain',ichain
#endif
c     --- update acceptance accumulators
      call twh_acsswitch(GLB_INCR,ipair,ibox,1.0d0)
c     --- update energies
      if ( potentialstyle .eq. POT_INTERNAL ) then
         call twh_boxvclassic(GLB_INCR,ibox,dvclassical_trial)
         call twh_boxvinter(GLB_INCR,ibox,dvinter_trial)
         call twh_boxvintra(GLB_INCR,ibox,dvintra_trial)
         call twh_boxvfield(GLB_INCR,ibox,dvext_trial)
         call twh_boxvelect(GLB_INCR,ibox,dvelect_trial)
         call twh_boxvsolvation(GLB_INCR,ibox,dvsolv_trial)
         call twh_boxvtorsion(GLB_INCR,ibox,dvtor_trial)
         call twh_boxvangle(GLB_INCR,ibox,dvbend_trial)
         call twh_boxvbond(GLB_INCR,ibox,dvvib_trial)
         call twh_scp_eval_derivative(GLB_GET,scp_eval_derivative)
         if (scp_eval_derivative) then
            call scp_increment_derivative(ibox, LAMBDA_C, ddvdl_c)
            call scp_increment_derivative(ibox, LAMBDA_LJ, ddvdl_lj)
         endif
         call scp_increment_foreign_energies(ibox, dvf)
      elseif ( potentialstyle .eq. POT_EXTERNAL ) then
         call twh_boxvquantum(GLB_INCR,ibox,dvquantum_trial)
      else
         write(6,*) 'SWITCH: unknown potentialstyle'
         write(6,*) potentialstyle
         stop
      endif

c     --- store new centers of mass
      call twh_com(GLB_GET,CRD_TEMP,ichain,xcmi,ycmi,zcmi)
      call twh_com(GLB_SET,CRD_REAL,ichain,xcmi,ycmi,zcmi)
      call twh_com(GLB_GET,CRD_TEMP,jchain,xcmj,ycmj,zcmj)
      call twh_com(GLB_SET,CRD_REAL,jchain,xcmj,ycmj,zcmj)

c     --- store new coordinates
      do iunit = 1,totiunit
         call twh_coordinates(GLB_GET,CRD_TEMP,ichain,iunit,rix,riy,riz)
         call twh_coordinates(GLB_SET,CRD_REAL,ichain,iunit,rix,riy,riz)
      enddo
      do junit = 1,totjunit
         call twh_coordinates(GLB_GET,CRD_TEMP,jchain,junit,rjx,rjy,rjz)
         call twh_coordinates(GLB_SET,CRD_REAL,jchain,junit,rjx,rjy,rjz)
      enddo
      
      if ( lewald ) then
c        --- update reciprocal-space sum
         call twh_recip(ibox,2)
      endif

c     --- update embedding density
      iboxa(1)=ibox
      call twh_reseteamrho(1,iboxa)


#if DEBUG_SWITCH
      write(6,*) 'end SWITCH'
#endif
      return
      end
