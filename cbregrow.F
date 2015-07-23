#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2011 Marcus G. Martin                       *
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
      subroutine twh_cbregrow(lfailure,cbmovetype)
c     ******************************************************************
c     * performs a lengthconserving configurational bias move          *
c     * for linear and  branched molecules                             *
c     * cbmovetype: CBR_ORIGINAL regrows from an atom to the endpoint. *
c     *           If this is a cyclic molecule and one of the cyclic   *
c     *           atoms is selected then the entire molecule is regrown*
c     * cbmovetype: CBR_BACKBONE selects an atom on the protein        *
c     *           backbone, finds another backbone atom that is        *
c     *           connected to the first by a torsion, and then regrows*
c     *           the atoms inbetween the two backbone atoms.          *
c     * cbmovetype: CBR_SIDECHAIN selects an alpha carbon and regrows  *
c     *           the sidechain attached to it.  This can include other*
c     *           sidechaines (eg disulfide linkages) but does not     *
c     *           move the main chain                                  *
c     *                                                                *
c     * rewritten from old config and branch subroutines by            *
c     * M.G. Martin 09-19-1997                                         *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_RANDOM
#define FUNCTION_EXPON
#include "functions.h"
c     --- debug toggle
#define DEBUG_CBREGROW 0
c     --- variables passed to/from the subroutine
      logical lfailure
      integer cbmovetype
c     --- local variables
c     --- logical scalars
      logical lterm,laccept,lnew,lcorrect,lpartial,lintra,lbonded
      logical lgetcom,scp_eval_derivative,lewald,lexschedule,ldualcutoff
c     --- logical arrays
      logical loverlap
      dimension loverlap(2)
c     --- integer scalars
      integer iname,ibox,numunits,icbu,islen,imolty,iustart,moveinfo
      integer ic,total,itest,sch_key,nchainmove,testcord,globalcord
      integer nmolty,ncmt,ipoint,nmtype,ensemble
c     --- integer arrays
      integer chaininfo
      dimension chaininfo(1,3)
      integer iboxa
      dimension iboxa(1)
c     --- double precision scalars
      double precision dpscalenew,dpscaleold,beta
      double precision dchain,rchain,wnlog,wolog,wdlog,wratio
      double precision burn,xtemp,ytemp,ztemp,pmtest,scalecut
      double precision dvtot,dvduma,dvinter,dvdumb,dvelect,dvsolv,dv
      double precision dvdumc,dvdumd,dvdume,dvdumf,ddvdl_lj,ddvdl_c
c     --- double precision arrays
      double precision wrosen
      dimension wrosen(2)
      double precision dvf
      dimension dvf(MAX_FOREIGN_LAMBDA)
      double precision vrostot
      dimension vrostot(2)
      double precision vrosvib
      dimension vrosvib(2)
      double precision vrosang
      dimension vrosang(2)
      double precision vrostor
      dimension vrostor(2)
      double precision vrosext
      dimension vrosext(2)
      double precision vrosintra
      dimension vrosintra(2)
      double precision vrosinter
      dimension vrosinter(2)
      double precision vroselect
      dimension vroselect(2)
      double precision vrosewald
      dimension vrosewald(2)

#if DEBUG_CBREGROW
      write(6,*) 'start CBREGROW'
#endif
c     --- retrieve constants
      call twh_lewald(GLB_GET,lewald)
      call twh_scalecut(GLB_GET,scalecut)

c     --- initialize cbmc growth variables
      call twh_resetcbmc(2,wrosen,vrostot,vrosvib,vrosang,vrostor
     &     ,vrosext,vrosintra,vrosinter,vroselect,vrosewald)
c     --- initialize energy changes
      dvtot = 0.0d0
      dvinter = 0.0d0
      dvelect = 0.0d0
      dvsolv = 0.0d0

c     --- select a molecule type at random using the appropriate array
      rchain  = twh_random()
      call twh_nmolty(GLB_GET,nmolty)
      do icbu = 1,nmolty
         if ( cbmovetype .eq. CBR_ORIGINAL ) then
c           --- use the normal cbmc matrix for finding a molecule type
            call twh_pmcbmt(GLB_GET,icbu,pmtest)
         elseif ( cbmovetype .eq. CBR_BACKBONE ) then
c           --- use the backbone cbmc matrix for finding a molecule type
            call twh_pmbkmt(GLB_GET,icbu,pmtest)
         elseif ( cbmovetype .eq. CBR_SIDECHAIN ) then
c           --- use the sidechainge matrix fo finding a molecule type
            call twh_pmcbsidemt(GLB_GET,icbu,pmtest)
         else
            write(6,*) 'CBREGROW: invalid cbmovetype:',cbmovetype
            lfailure = .true.
            return
         endif
         if ( rchain .lt. pmtest ) then
            imolty = icbu
            rchain = 2.0d0
         endif
      enddo

      call twh_ensemble(GLB_GET,ensemble)
      if ( ensemble .eq. ENS_muVT ) then
c        --- grand canonical, take from box 1
         ibox = 1
c        --- make sure some molecules of this type are in box 1
         call twh_ncmt(GLB_GET,ibox,imolty,ncmt)
         if ( ncmt .eq. 0 )  return
c        --- select molecule of type imolty in box ibox
         dchain = twh_random()*dble(ncmt)
         ipoint = idint(dchain) + 1
         call twh_chainlist(GLB_GET,ipoint,ibox,imolty,iname)
      else
         call twh_nmolectyp(GLB_GET,imolty,nmtype)
         dchain = twh_random()*dble(nmtype)
         ipoint = idint( dchain ) + 1
         call twh_parall(GLB_GET,imolty,ipoint,iname)
         call twh_nboxi(GLB_GET,iname,ibox)
      endif

      call twh_moltyp(GLB_GET,iname,itest)
      if ( itest .ne. imolty ) then
         write(6,*) 'CBREGROW: screwup cbregrow'
         lfailure = .true.
         return
      endif

c     --- store number of units in numunits
      call twh_nunit(GLB_GET,imolty,numunits)

c     --- store initial position for the trial chain as the same as the
c     --- current positions
      do icbu = 1, numunits
         call twh_coordinates(GLB_GET,CRD_REAL,iname,icbu
     &        ,xtemp,ytemp,ztemp)
         call twh_coordinates(GLB_SET,CRD_TEMP,iname,icbu
     &        ,xtemp,ytemp,ztemp)
      enddo

      if ( cbmovetype .eq. CBR_ORIGINAL ) then
c        --- regular cb growth
         sch_key = 1
      elseif ( cbmovetype .eq. CBR_BACKBONE ) then
c        --- backbone cb growth
         sch_key = 4
      elseif ( cbmovetype .eq. CBR_SIDECHAIN ) then
c        --- sidechain cb growth
         sch_key = 6
      endif
      call twh_schedule(lfailure,numunits,imolty,islen,iustart,moveinfo
     &     ,0,sch_key)
      if ( lfailure ) return

      total = 0
      if ( cbmovetype .eq. CBR_ORIGINAL ) then
c        --- determine how many beads are being regrown
         do icbu = 1,numunits
            call twh_lexsched(GLB_GET,icbu,lexschedule)
            if ( .not. lexschedule ) total = total + 1
         enddo
         call twh_bncb(GLB_INCR,imolty,total,ibox,1)
      elseif ( cbmovetype .eq. CBR_BACKBONE ) then
c        --- store according to the atom and torsion
         call twh_bnback(GLB_INCR,imolty,iustart,moveinfo,ibox,1)
      elseif ( cbmovetype .eq. CBR_SIDECHAIN ) then
c        --- store according to the amino acid number (moveinfo)
         call twh_bncbside(GLB_INCR,imolty,moveinfo,ibox,1)
      else
         write(6,*) 'CBREGROW: attempted move unknown cbmovetype:'
     &        ,cbmovetype
         lfailure = .true.
         return
      endif
c     --- always use real coordinates for all other molecules
      globalcord = CRD_REAL

c     --- grow new chain conformation
      lnew = .true.
      testcord = CRD_TEMP
      call twh_rosenbluth (lfailure,lnew,lterm,testcord,globalcord,iname
     &     ,imolty,islen,ibox,numunits,dpscalenew,wrosen(CRD_NEW)
     &     ,vrostot(CRD_NEW),vrosvib(CRD_NEW),vrosang(CRD_NEW)
     &     ,vrostor(CRD_NEW),vrosext(CRD_NEW),vrosintra(CRD_NEW)
     &     ,vrosinter(CRD_NEW),vroselect(CRD_NEW),vrosewald(CRD_NEW))
      if ( lfailure ) return
c     --- termination of cbmc attempt due to walk termination
      if ( lterm ) return

c     --- grow old chain conformation
      lnew = .false.
      testcord = CRD_REAL
      call twh_rosenbluth (lfailure,lnew,lterm,testcord,globalcord,iname
     &     ,imolty,islen,ibox,numunits,dpscaleold,wrosen(CRD_OLD)
     &     ,vrostot(CRD_OLD),vrosvib(CRD_OLD),vrosang(CRD_OLD)
     &     ,vrostor(CRD_OLD),vrosext(CRD_OLD),vrosintra(CRD_OLD)
     &     ,vrosinter(CRD_OLD),vroselect(CRD_OLD),vrosewald(CRD_OLD))
      if ( lfailure ) return

      if ( lterm ) then
         write(6,*) 'CBREGROW: old growth had zero rosenbluth weight'
         write(6,*) 'This means the intitial conformation of the'
         write(6,*) 'molecule is extremely high energy'
         write(6,*) 'the move cannot be accepted because it is not'
     &        ,' reversible'
         return
      endif

C     -----------------------------------------------------------------
c     Begin DC-CBMC and Ewald-sum Corrections 

      lcorrect = .false.
      call twh_ldualcutoff(GLB_GET,ldualcutoff)
      if ( ldualcutoff ) lcorrect = .true.
#if USETRAMONTO
      lcorrect = .true.
#endif
      if ( lcorrect ) then
c        calculate the Full rcut vdw energy for grown beads
c        --- vintra and vext are not recalculated in this call as they
c        --- were already computed properly in engatom.
c        --- only the difference of intermolecular energies is computed
c        --- in this energy call
         nchainmove = 1
         chaininfo(1,CRD_OLD) = CRD_REAL
         chaininfo(1,CRD_NEW) = CRD_TEMP
         chaininfo(1,3) = iname
         lpartial = .true.
         lgetcom = .true.
         lintra = .false.
         lbonded = .false.
         call twh_energy_change(lfailure,ibox,nchainmove,chaininfo
     &        ,loverlap,lpartial,lgetcom,lintra,lbonded,dvtot,dvduma
     &        ,dvinter,dvdumb,dvelect,dvsolv,dvdumc,dvdumd,dvdume
     &        ,dvdumf,dvf,ddvdl_lj,ddvdl_c)
         if ( lfailure ) return
c        --- add all of the energy change into the new array
         call twh_beta(GLB_GET,beta)
         wrosen(2) = wrosen(2) *twh_expon(-(beta*dvtot))
         vrostot(2) = vrostot(2) + dvtot
         vrosinter(2) = vrosinter(2) + dvinter
         vroselect(2) = vroselect(2) + dvelect
      endif

c     --- End of DC-CBMC and Ewald-sum Corrections

c     --- we at least regrew well enough that we have a chance to accept
      if ( cbmovetype .eq. CBR_ORIGINAL ) then
c        --- regular cbmc
         call twh_bscb(GLB_INCR,imolty,total,ibox,1)
      elseif ( cbmovetype .eq. CBR_BACKBONE ) then
c        --- backbone torsion regrowth
         call twh_bsback(GLB_INCR,imolty,iustart,moveinfo,ibox,1)
      elseif ( cbmovetype .eq. CBR_SIDECHAIN ) then
c        --- sidechain regrowth
         call twh_bscbside(GLB_INCR,imolty,moveinfo,ibox,1)
      else
         write(6,*) 'CBREGROW: grown move unknown cbmovetype:'
     &        ,cbmovetype
         lfailure = .true.
         return
      endif

c     --- check for acceptance of trial configuration
      laccept = .false.
c     --- add in the effect of the rescaling
      wnlog = dlog10 ( wrosen(2) ) - dpscalenew*scalecut
c     --- it is possible that wrosen(1) is zero, if so then be careful
      if ( wrosen(1) .gt. 0.0d0 ) then
         wolog = dlog10( wrosen(1) ) - dpscaleold*scalecut
         wdlog = wnlog - wolog
c         write(6,*) 'wdlog,wnlog,wolog',wdlog,wnlog,wolog
c         write(6,*) 'wrosen(1),dpscaleold',wrosen(1),dpscaleold
c         write(6,*) 'wrosen(2),dpscalenew',wrosen(2),dpscalenew
      else
c        --- zero or lower rosenbluth weight should not occur
         write(6,*) 'CBREGROW: old Rosenbluth weight is 0.0 or less'
         write(6,*) wrosen(1)
         write(6,*) 'This should never happen'
         lfailure = .true.
         return
      endif

c     --- make sure we avoid overflow errors and also check to see 
c     --- if we even need to generate a random number
      if ( wdlog .gt. 0.0d0 ) then
c        --- we will accept the move as wrosen(2)/wrosen(1) > 1.0
         laccept = .true.
c        --- burn a random number
         burn = twh_random()
      else
         wratio = 10.0**(wdlog)
         if ( twh_random() .lt. wratio ) laccept = .true.
      endif

#if DEBUG_CBREGROW
      write(6,*) 'vrosinter',vrosinter
      write(6,*) 'vrosintra',vrosintra
      write(6,*) 'vroselect',vroselect
      write(6,*) 'vrosewald',vrosewald
      write(6,*) 'vrosext',vrosext
      write(6,*) 'vrosvib',vrosvib
      write(6,*) 'vrosang',vrosang
      write(6,*) 'vrostor',vrostor
      write(6,*) 'delta vib',vrosvib(2)-vrosvib(1)
      write(6,*) 'delta ang',vrosang(2)-vrosang(1)
      write(6,*) 'delta tor',vrostor(2)-vrostor(1)
      write(6,*) 'delta intra',vrosintra(2)-vrosintra(1)
      write(6,*) 'delta inter',vrosinter(2)-vrosinter(1)
      write(6,*) 'delta elect',vroselect(2)-vroselect(1)
      write(6,*) 'delta total',vrostot(2)-vrostot(1)
      write(6,*) 'wrosen',wrosen
      write(6,*) 'wdlog ',wdlog
      write(6,*) 'accept',laccept
c     --- force acceptance to debug
c      if ( .not. loverlap(CRD_NEW) ) laccept = .true.
#endif

      if ( laccept ) then
#if DEBUG_CBREGROW
      if ( cbmovetype .eq. CBR_ORIGINAL ) then
         write(6,*) 'CBREGROW acpt',iname,ibox,vrostot,wdlog
      elseif ( cbmovetype .eq. CBR_BACKBONE ) then
         write(6,*) 'CB-BACKBONE acpt',iname,ibox,vrostot,wdlog
      elseif ( cbmovetype .eq. CBR_SIDECHAIN ) then
         write(6,*) 'CB-SIDECHAIN acpt',iname,ibox,vrostot,wdlog
      endif
#endif
c        --- we can now accept 
         if ( cbmovetype .eq. CBR_ORIGINAL ) then
            call twh_bacb(GLB_INCR,imolty,total,ibox,1)
         elseif ( cbmovetype .eq. CBR_BACKBONE ) then
            call twh_baback(GLB_INCR,imolty,iustart,moveinfo,ibox,1)
         elseif ( cbmovetype .eq. CBR_SIDECHAIN ) then
            call twh_bacbside(GLB_INCR,imolty,moveinfo,ibox,1)
         else
            write(6,*) 'CBREGROW: accepted unknown cbmovetype:'
     &           ,cbmovetype
         endif
         dv = vrostot(2) - vrostot(1)
         call twh_boxvclassic(GLB_INCR,ibox,dv)
         dv = vrosinter(2) - vrosinter(1)
         call twh_boxvinter(GLB_INCR,ibox,dv)
         dv = vrosintra(2) - vrosintra(1)
         call twh_boxvintra(GLB_INCR,ibox,dv)
         dv = vrosvib(2) - vrosvib(1)
         call twh_boxvbond(GLB_INCR,ibox,dv)
         dv = vrostor(2) - vrostor(1)
         call twh_boxvtorsion(GLB_INCR,ibox,dv)
         dv = vrosext(2) - vrosext(1)
         call twh_boxvfield(GLB_INCR,ibox,dv)
         dv = vrosang(2) - vrosang(1)
         call twh_boxvangle(GLB_INCR,ibox,dv)
         dv = vroselect(2) - vroselect(1) + vrosewald(2) - vrosewald(1)
         call twh_boxvelect(GLB_INCR,ibox,dv)
         call twh_boxvsolvation(GLB_INCR,ibox,dvsolv)

         call twh_scp_eval_derivative(GLB_GET,scp_eval_derivative)
         if (scp_eval_derivative) then
            call scp_increment_derivative(ibox, LAMBDA_C, ddvdl_c)
            call scp_increment_derivative(ibox, LAMBDA_LJ, ddvdl_lj)
         endif
         call scp_increment_foreign_energies(ibox, dvf)
   
         do ic = 1, numunits
            call twh_coordinates(GLB_GET,CRD_TEMP,iname,ic
     &           ,xtemp,ytemp,ztemp)
            call twh_coordinates(GLB_SET,CRD_REAL,iname,ic
     &           ,xtemp,ytemp,ztemp)
         enddo
         
         if ( lewald ) then
c           --- update reciprocal-space sum
            call twh_recip(ibox,2)
         endif
         
c        --- update center of mass
         call twh_ctrmas(lfailure,CRD_REAL,ibox,iname,CTR_CBREGROW)
         if ( lfailure ) return
c        --- update embedding density
         iboxa(1) = ibox
         call twh_reseteamrho(1,iboxa)
      endif

#if DEBUG_CBREGROW
      write(6,*) 'end CBREGROW'
#endif
      return
      end
