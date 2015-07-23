#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2014 Marcus G. Martin                       *
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
      subroutine twh_engtotal(lfailure,lovrlap,vclassical,vinter
     &,vthrinter,vtail
     &,vintra,vvib,vbend,vtg,vext,velect,vsolvation,vquantum,vforeign,
     & dl_inter, dl_elect,ibox, lmolecmove, lmovelist, lprint,lbonded)
c     ******************************************************************
c     * calculates the total potential energy for a configuration.     *
c     *                                                                *
c     * The following are set by the subroutine                        *
c     * lovrlap is true if this subroutine finds substantial overlap   *
c     * vclassical is the total classical energy                       *
c     * vinter is the total intermolecular two-body, three-body and    *
c     *        tail correction non-coulombic energy                    *
c     * vtail: just the tail correction terms of the intermolecular    *
c     *        energy                                                  *
c     * vthrinter: just the 3-body terms of the intermolecular energy  *
c     * vintra is the intramolecular (non-coulombic) non-bonded energy *
c     * vvib is the bond vibration energy                              *
c     * vbend is the bond bending energy                               *
c     * vtg is the dihedral angle energy - including improper torsions *
c     * vext is the energy of external fields/walls                    *
c     * velect is the total coulombic energy - execpt 1-4 coulombics   *
c     *        which are included in the torsional potential           *
c     * vsolvation is the solvation energy                             *
c     * vquantum is the quantum mechanical energy                      *
c     * vforeign is total classical energy evaluated at different      *
c     *     values of lambda for 'Scaled Lennard-Jones' potential.     *
c     *     See scaled_potential.F                                     *
c     * dl_inter -- derivative of two body potential with respect      *
c     *     to lambda_lj, calculated if scp_eval_derivative true       *
c     *     (else 0)                                                   *
c     * dl_elect -- derivative of inter-molecule coulomb potential     *
c     *     with respect to lambda_c                                   *
c     *                                                                *
c     * The following are set by the caller                            *
c     * ibox is the box we are computing the energy in                 *
c     * lmolecmove is true if this call is from a Monte Carlo move,    *
c     *            otherwise this is the call at the start or end of   *
c     *            the simulation.  If this is a move then we use the  *
c     *            additional COM cutoff.                              *
c     * lmovelist is true if we are performing a move and only need the*
c     *            interactions for a subset of the molecules          *
c     * lprint is true if we want to print out all of the energies     *
c     * lbonded is true if we are computing bonded interactions        *
c     *            (vibration, bending, torsion etc)                   *
c     *                                                                *
c     * Error conditions:                                              *
c     *   May exit early if overlap found (error may be printed)       *
c     *                                                                *
c     * originally written prior to 1995                               *
c     * last modified 01-09-2012 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_ELEMENTNAME
#define FUNCTION_GET_NBOXI
#define FUNCTION_GET_NCMT
#define FUNCTION_GET_NTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure,lovrlap, lmolecmove,lmovelist,lprint,lbonded
      integer ibox
      double precision vclassical,vinter,vintra,vthrinter,vtail,vvib
      double precision vbend,vtg,velect,vsolvation,vext,vquantum
      double precision vforeign,dl_inter, dl_elect
      dimension vforeign(MAX_FOREIGN_LAMBDA)

c     --- local variables
c     --- logical scalars
      logical lmadelung,ltailc,lcoulomb,lewald
c     --- character scalars
      character*20 external_code
c     --- integer scalars
      integer ichain, nmcount, cutstyle,ispecial,flagon,l
      integer classical_potnum,nchain,nmolty,nchbox
      integer scp_num_foreign_lambda,potentialstyle,numboxes
c     --- double precision scalars
      double precision maddist,sself,correct,vtmp,qqfact
      double precision velectinter,velectintra,vtwointer,vangleangle
      double precision vimpr,vrealelect,vrecipsum,vbondbond
c     --- double precision arrays
      integer chaininfo
      dimension chaininfo(1,3)

c      write(6,*) 'start ENGTOTAL'

c     --- retrieve constants
      call twh_constant_qqfact(GLB_GET,qqfact)
      call twh_potentialstyle(GLB_GET,potentialstyle)
      call twh_ltailc(GLB_GET,ltailc)
      call twh_lcoulomb(GLB_GET,lcoulomb)
      call twh_lewald(GLB_GET,lewald)

      cutstyle = 0
      lovrlap = .false.
      lmadelung = .false.

      vclassical = 0.0d0
      vinter = 0.0d0
      vtwointer = 0.0d0
      vthrinter = 0.0d0
      vintra = 0.0d0
      vtail = 0.0d0
      vtg = 0.0d0
      vimpr = 0.0d0
      vbend = 0.0d0
      vangleangle = 0.0d0
      vvib = 0.0d0
      vbondbond = 0.0d0 
      vext = 0.0d0
      velect = 0.0d0
      sself = 0.0d0
      correct = 0.0d0
      vrealelect = 0.0d0
      vrecipsum = 0.0d0
      vsolvation = 0.0d0
      velectinter = 0.0d0
      velectintra = 0.0d0
      vquantum = 0.0d0
      dl_inter = 0.0d0
      dl_elect = 0.0d0
      do l=1,MAX_FOREIGN_LAMBDA
         vforeign(l) = 0.0d0
      enddo

c     --- determine whether we want to compute the madelung constant
      call twh_nmolty(GLB_GET,nmolty)
      call twh_numboxes(GLB_GET,numboxes)
      if ( numboxes .eq. 1 .and. nmolty .eq. 2 
     &     .and. twh_get_ncmt(ibox,1) .eq. twh_get_ncmt(ibox,nmolty)
     &     .and. twh_get_elementname(twh_get_ntype(1,1)) .eq. 'Na'
     &     .and. twh_get_elementname(twh_get_ntype(nmolty,1)) .eq. 'Cl'
     &     ) then
         lmadelung = .true.
      endif

c     ---  check the molecule count
      call twh_nchain(GLB_GET,nchain)
      nmcount = 0
      do ichain = 1, nchain
         if ( twh_get_nboxi(ichain) .eq. ibox ) then
            nmcount = nmcount + 1
         endif
      enddo

      call twh_nchbox(GLB_GET,ibox,nchbox)
      if ( nmcount .ne. nchbox ) then
         write(6,*) 'ENGTOTAL: nmcount ne nchbox'
         write(6,*) 'ibox nmcount ',ibox,nmcount
         write(6,*) 'nchbox ',nchbox
         write(6,*) 'nchain ',nchain
         lfailure = .true.
         return
      endif
 
      if ( potentialstyle .eq. POT_INTERNAL ) then
         call twh_engtotal_twobody(lfailure,velectinter,vtwointer
     &        ,vforeign,maddist,dl_inter,dl_elect,lovrlap,lmadelung
     &        ,lmolecmove,lmovelist,cutstyle,ibox )
         if ( lfailure ) return
         if ( lovrlap ) then
            return
         endif

         if ( ltailc ) then
            call twh_engtotal_tail(vtail, ibox)
         endif

         call twh_cp_number(GLB_GET,classical_potnum)
         if ( classical_potnum .eq. CP_STILLINGER_WEBER ) then
            call  twh_engtotal_tb_sw(lfailure,vthrinter, ibox)
            if ( lfailure ) return
         elseif ( classical_potnum .eq. CP_EMBEDDED_ATOM_METHOD ) then
            call twh_engtotal_tb_eam(lfailure,vthrinter, lovrlap, ibox)
            if ( lfailure ) return
         endif

         if ( lewald ) then
c     TODO: add foreign energy terms here
            call twh_engtotal_ewald(lfailure,vrecipsum,sself,correct
     &           , ibox )
            if ( lfailure ) return
         endif

         call twh_engtotal_intrachain(lfailure,vintra, velectintra
     &        , lovrlap,cutstyle, ibox, lmolecmove, lmovelist)
         if ( lfailure ) return
         if ( lovrlap ) then
            return
         endif   
    
         if ( lbonded ) then
            call twh_engtotal_vibbendtors(lfailure,vvib,vbend,vtg
     &           ,vangleangle,vimpr,vbondbond,lovrlap,lprint,ibox)
            if ( lfailure ) return
            if ( lovrlap ) then
               return
            endif
         endif 
         call twh_engtotal_ext(lfailure,vext,lovrlap,ibox,lmolecmove
     &        ,lmovelist )
         if ( lfailure ) return
         if ( lovrlap ) return
         
         call twh_engtotal_solv(lfailure, vsolvation, ibox )
         if ( lfailure ) return
      
         if ( lcoulomb ) then
            vrealelect = velectintra + velectinter
            velect = vrealelect + sself + correct + vrecipsum/qqfact
            velect = qqfact*velect 
         endif

         vinter = vtwointer + vthrinter + vtail

         vtmp = vintra + vvib + vbondbond + vbend 
     &        + vangleangle + vtg + vimpr + vext + vsolvation
         vclassical = vinter + velect + vtmp

c        --- foreign energies add terms which do not change with lambda
         call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
         do l = 1, scp_num_foreign_lambda
            vforeign(l)=vforeign(l) + vtmp
         enddo

      endif

c     --- compute external energies
      if ( potentialstyle .eq. POT_EXTERNAL ) then
         call twh_external_code(GLB_GET,external_code)
         if ( external_code .eq. 'lcao' ) then
            ispecial = 0
            flagon = 0
            chaininfo(1,1)=-1
            chaininfo(1,2)=-1
            chaininfo(1,3)=-1
            call twh_quantumenergy(lfailure,ibox,flagon,1,chaininfo
     &           ,vquantum)
            if ( lfailure ) return
         elseif ( external_code .eq. 'KIM' ) then
            write(6,*) 'hey, should compute this'
            lfailure = .true.
            return
         else
            write(6,*) 'ENGTOTAL: unknown external_code'
            write(6,*) external_code
            lfailure = .true.
            return
         endif
      endif
         
      if ( lprint ) then
         call twh_eng_print(nmcount, vvib, vbondbond, vbend, 
     &       vangleangle, vtg, vimpr, vintra, vinter, vtwointer, 
     &       vthrinter, velect, vrealelect, sself, correct, vrecipsum, 
     &       vext, vsolvation, vclassical, vtail, velectinter, 
     &       velectintra,vquantum,vforeign,dl_inter,dl_elect, lmadelung,
     &       maddist, ibox)
      endif

c     --- add the bond-bond terms into the regular bond
      vvib = vvib + vbondbond
c     --- add the angleangle terms into the regular angle
      vbend = vbend + vangleangle
c     --- add the improper terms into the regular torsions
      vtg = vtg + vimpr
c      write(6,*) 'end ENGTOTAL'
      return
      end

      subroutine twh_engtotal_twobody(lfailure,velectinter,vtwointer
     &     ,vforeign,maddist,dl_twointer,dl_elect,lovrlap,lmadelung
     &     ,lmolecmove,lmovelist,cutstyle,ibox)
c     ******************************************************************
c     * Evaluate two body non-bonded interactions                      *
c     * Arguments set by subroutine:                                   *
c     *  velectinter -- the intermolecular component of real space     *
c     *     coulombic potential                                        *
c     *  vtwointer -- the two-body non-bond potential                  *
c     *  vforeign -- coulomb and twobody energies at different lambda  *
c     *  maddist -- Distance for Madelung constant                     *
c     *  dl_twointer -- derivative of two body potential with respect  *
c     *     to lambda_lj, calculated if scp_eval_derivative true       *
c     *  dl_elect -- derivative of inter-molecule coulomb potential    *
c     *     with respect to lambda_lj                                  *
c     *                                                                *
c     * Arguments modified (but not initialized) by subroutine:        *
c     *  lovrlap -- Atomic overlap found?                              *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  lmadelung -- compute madelung constant?                       *
c     *  lmolecmove -- is true if this call is from a Monte Carlo move *
c     *  lmovelist -- true if we are performing a move and only need   *
c     *            the interactions for a subset of the molecules      *
c     *  cutstyle -- cutoff style.  See twh_vtwobody for definition    *
c     *  ibox -- the box we are computing the energy in                *
c     *                                                                *
c     * All energies and derivatives in units of [K], with the         *
c     * exception of velectinter, which is in units of [EE]            *
c     *                                                                *
c     * Error conditions:                                              *
c     *   May halt program                                             *
c     *   May exit early (sets lovrlap=true)                           *
c     *                                                                *
c     * split from original engtotal 08-11-2005 by MAW                 *
c     * last modified 03-03-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NBOXI
#define FUNCTION_VCOULOMB
#define FUNCTION_VTWOBODY
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision velectinter, vtwointer, maddist
      double precision dl_twointer, dl_elect
      double precision vforeign(MAX_FOREIGN_LAMBDA)
      logical lfailure,lmadelung, lmolecmove, lmovelist, lovrlap
      integer cutstyle, ibox

c     --- local variables
c     --- logical scalars
      logical lcompute,leval_c,leval_v,lusecom,lmovemei,lmovemej
      logical luse_ci,luse_cj,luse_vi,luse_vj,lcoulomb
c     --- integer scalars
      integer ichain,imolty,jchain,jmolty,ii,ntii,jj,ntjj,l,nchain
      integer numiunit,numjunit,scp_num_foreign_lambda
c     --- double precision scalars
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij
      double precision xcmi,ycmi,zcmi,rcmi,rcm,rcmsq,rijsq,rij,dvdl
      double precision qqfact,xcmj,ycmj,zcmj,rcmj,rcelect
c     --- vf holds incremental foreign energies from vtwobody, vcoulomb
c     --- calls
c     --- vfelec and vftwo keep running totals, their sum becomes
c     --- vforeign
      double precision vf, vfelec, vftwo
      dimension vf(MAX_FOREIGN_LAMBDA), vfelec(MAX_FOREIGN_LAMBDA), 
     &      vftwo(MAX_FOREIGN_LAMBDA)

c     --- retrieve constants
      call twh_constant_qqfact(GLB_GET,qqfact)
      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
      call twh_lusecom(GLB_GET,lusecom)
      call twh_rcelect(GLB_GET,ibox,rcelect)
      call twh_lcoulomb(GLB_GET,lcoulomb)

      vtwointer = 0.0d0
      velectinter = 0.0d0
      do l=1,MAX_FOREIGN_LAMBDA
         vfelec(l) = 0.0d0
         vftwo(l) = 0.0d0
      enddo
      maddist = 1d10
      dl_twointer = 0.0d0
      dl_elect = 0.0d0
      rcmi = 0.0d0

c     --- two-body nonbonded interactions
c     --- loop over all chains i 
      call twh_nchain(GLB_GET,nchain)
      do ichain = 1, nchain - 1
c        --- check if i is in relevant box
         if ( twh_get_nboxi(ichain) .eq. ibox ) then
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,numiunit)
            if ( lmolecmove ) then
               call twh_com(GLB_GET,CRD_REAL,ichain,xcmi,ycmi,zcmi)
               call twh_rcmu(GLB_GET,ichain,rcmi)
            endif
c           --- loop over all chains j with j>i 
            do jchain = ichain + 1, nchain
c              --- check for simulation box
               if ( twh_get_nboxi(jchain) .eq. ibox ) then
                  call twh_moltyp(GLB_GET,jchain,jmolty)
                  call twh_nunit(GLB_GET,jmolty,numjunit)
                  lcompute = .true.
                  if ( lmovelist ) then
c                    --- see if at least one of these molecules is on
c                    --- the movelist, otherwise skip
                     call twh_lmoveme(GLB_GET,ichain,lmovemei)
                     call twh_lmoveme(GLB_GET,jchain,lmovemej)
                     if ( .not. (lmovemei .or. lmovemej) ) then
c                       --- neither of these molecules moved, do not
c                       --- compute energy
                        lcompute = .false.
                     endif
                  endif
                  if ( lmolecmove .and. lusecom ) then
c                    --- check if ctrmas within rcmsq
c                    --- minimum image the ctrmas pair separations
                     call twh_com(GLB_GET,CRD_REAL,jchain
     &                    ,xcmj,ycmj,zcmj)
                     call twh_mimage ( ibox,xcmi,ycmi,zcmi,xcmj
     &                    ,ycmj,zcmj,rxuij,ryuij,rzuij)
                     rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                     call twh_rcmu(GLB_GET,jchain,rcmj)
c                    --- start with the vdw cutoff
                     call twh_rcut(GLB_GET,rcm)
                     if ( lcoulomb ) then
c                       --- check the coulombic cutoff
                        if ( rcelect .gt. rcm ) rcm = rcelect
                     endif
c                    --- add in the com radii
                     rcm = rcm + rcmi + rcmj
                     rcmsq = rcm*rcm
                     if ( rijsq .gt. rcmsq ) lcompute = .false.
                  endif
                  if ( lcompute ) then
                     do ii = 1,numiunit
                        call twh_ntype(GLB_GET,imolty,ii,ntii)
                        call twh_coordinates(GLB_GET,CRD_REAL,ichain,ii
     &                       ,rxui,ryui,rzui)
                        call twh_luse_c(GLB_GET,imolty,ii,luse_ci)
                        call twh_luse_v(GLB_GET,imolty,ii,luse_vi)
c                       --- loop over all beads jj of chain j 
                        do jj = 1, numjunit
                           call twh_ntype(GLB_GET,jmolty,jj,ntjj)
                           call twh_coordinates(GLB_GET,CRD_REAL,jchain
     &                          ,jj,rxuj,ryuj,rzuj)
                           call twh_luse_c(GLB_GET,jmolty,jj,luse_cj)
                           call twh_luse_v(GLB_GET,jmolty,jj,luse_vj)
c                          --- check if nonzero interaction sites
                           leval_v = luse_vi .and. luse_vj
                           leval_c = luse_ci .and. luse_cj
                           if ( leval_v .or. leval_c ) then
c                          --- minimum image the pair separation
                              call twh_mimage ( ibox,rxui,ryui,rzui,rxuj
     &                          ,ryuj,rzuj,rxuij,ryuij,rzuij )
                              rijsq = (rxuij*rxuij)+(ryuij*ryuij)
     &                          + (rzuij*rzuij)
                              if (leval_v) then
                                 vtwointer = vtwointer 
     &                              + twh_vtwobody(lfailure,rijsq,ntii
     &                                ,ntjj,cutstyle,lovrlap,vf,dvdl)
                                 if ( lfailure ) return
                                 dl_twointer = dl_twointer + dvdl
                                 do l = 1, scp_num_foreign_lambda
                                    vftwo(l)=vftwo(l) + vf(l)
                                 enddo

                                 if ( lovrlap ) then
                                    call twh_engtotal_overlap(lfailure
     &                                   ,lmolecmove,ichain,ii,jchain,jj
     &                                   ,rijsq,rxui,ryui,rzui,rxuj,ryuj
     &                                   ,rzuj)
                                    if ( lfailure ) return
                                 endif
                              endif
c                             --- electrostatics.
c                             --- qqfact converts dl_elect and vforeign
c                             --- to [K]; velectinter stays [EE]
                              if ( lcoulomb .and. leval_c) then
                                 rij = dsqrt(rijsq)
                                 velectinter = velectinter 
     &                                + twh_vcoulomb(ibox,rij,imolty,ii
     &                                ,jmolty,jj,.false.,vf,dvdl
     &                                ,lovrlap)
                                 if ( lovrlap ) then
                                    call twh_engtotal_overlap(lfailure
     &                                   ,lmolecmove,ichain,ii,jchain,jj
     &                                   ,rijsq,rxui,ryui,rzui,rxuj,ryuj
     &                                   ,rzuj)
                                    if ( lfailure ) return
                                 endif
                                 dl_elect = dl_elect + dvdl * qqfact
                                 do l = 1, scp_num_foreign_lambda
                                    vfelec(l)=vfelec(l)+qqfact*vf(l)
                                 enddo
                                 if ( lmadelung ) then
                                   if ( rij .lt. maddist ) maddist = rij
                                 endif
                              endif
                           endif
c                          --- loop over all beads jj of chain j 
                        enddo
                     enddo
c                    --- lcompute
                  endif
               endif
            enddo
         endif
      enddo

      do l=1,MAX_FOREIGN_LAMBDA
         vforeign(l) = vfelec(l) + vftwo(l)
      enddo

      return
      end


      subroutine twh_engtotal_tail(vtail, ibox)
c     ******************************************************************
c     * Calculate tail corrections beyond the cutoff                   *
c     * Arguments set by subroutine:                                   *
c     *  vtail -- just the tail correction terms of the intermolecular *
c     *        energy                                                  *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  ibox -- the box we are computing the energy in                *
c     *                                                                *
c     * Error conditions:                                              *
c     *   None                                                         *
c     *                                                                *
c     * split from original engtotal 08-11-2005 by MAW                 *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision vtail
      integer ibox
c     --- local variables
      integer imolty,jmolty,nmolty,ncmt
      double precision rho,coru,boxvolume

c     --- retrieve constants
      call twh_boxvolume(GLB_GET,ibox,boxvolume)
      vtail = 0.0d0
c     --- add tail corrections beyond the cutoff
      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1, nmolty
         do jmolty = 1, nmolty
            call twh_ncmt(GLB_GET,ibox,jmolty,ncmt)
            rho = dble(ncmt) / boxvolume
            call twh_ncmt(GLB_GET,ibox,imolty,ncmt)
            call twh_coru(GLB_GET,imolty,jmolty,coru)
            vtail = vtail + dble(ncmt)*rho*coru
         enddo
      enddo
      return
      end

      subroutine twh_engtotal_tb_sw(lfailure,vthrinter,ibox)
c     ******************************************************************
c     * calculates the three-body energy contribution                  *
c     * Arguments set by subroutine:                                   *
c     *  vthrinter -- just the 3-body terms of the intermolecular      *
c     *     energy                                                     *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  ibox -- the box we are computing the energy in                *
c     *                                                                *
c     * Error conditions:                                              *
c     *  None                                                          *
c     *                                                                *
c     * split from original engtotal 08-12-2005 by MAW                 *
c     * last modified 11-01-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_GET_NBOXI
#define FUNCTION_VTHREEBODY
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      double precision vthrinter
      integer ibox
c     --- local variables
      integer ichain,imolty,ntii,jchain,jmolty,ntjj,kchain,kmolty,ntkk
      integer nchain
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij
      double precision rxuk,ryuk,rzuk,rij,bij,duma,thetac
      double precision rxujk,ryujk,rzujk,bjk,rjk
      logical lzero

      vthrinter = 0.0d0

c     --- stillinger-weber three body potential
c     --- this only works for chains that are monatomic 
      call twh_nchain(GLB_GET,nchain)
      do ichain = 1,nchain
         if ( twh_get_nboxi(ichain) .eq. ibox ) then
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_ntype(GLB_GET,imolty,1,ntii)
            call twh_coordinates(GLB_GET,CRD_REAL,ichain,1
     &           ,rxui,ryui,rzui)
            do jchain = 1,nchain
               if ( twh_get_nboxi(jchain) .eq. ibox 
     &              .and. ichain .ne. jchain) then
                  call twh_moltyp(GLB_GET,jchain,jmolty)
                  call twh_ntype(GLB_GET,jmolty,1,ntjj)
                  call twh_coordinates(GLB_GET,CRD_REAL,jchain,1
     &                 ,rxuj,ryuj,rzuj)
c                    --- compute i-j bond order
                  call twh_mimage(ibox,rxui,ryui,rzui,rxuj,ryuj,rzuj
     &                 ,rxuij,ryuij,rzuij)
                  rij = twh_distance(rxuij,ryuij,rzuij)
                  call twh_bondorder(lfailure,ntii,ntjj,rij,bij,duma
     &                 ,lzero)
                  if ( lfailure ) return
                  if ( .not. lzero ) then
                     do kchain = 1,nchain
                        if ( twh_get_nboxi(kchain) .eq. ibox 
     &                       .and. ichain .lt. kchain 
     &                       .and. jchain .ne. kchain ) then
                           call twh_moltyp(GLB_GET,kchain,kmolty)
                           call twh_ntype(GLB_GET,kmolty,1,ntkk)
                           call twh_coordinates(GLB_GET,CRD_REAL
     &                          ,kchain,1,rxuk,ryuk,rzuk)
c                             --- compute j-k bondorder
                           call twh_mimage(ibox,rxuk,ryuk,rzuk,rxuj
     &                          ,ryuj,rzuj,rxujk,ryujk,rzujk)
                           rjk = twh_distance(rxujk,ryujk,rzujk)
                           call twh_bondorder(lfailure,ntjj,ntkk,rjk,bjk
     &                          ,duma,lzero)
                           if ( lfailure ) return
                           if ( .not. lzero ) then
c                                --- compute angle
                              thetac = (rxuij*rxujk + ryuij*ryujk 
     &                             + rzuij*rzujk)/(rij*rjk)
c                                --- compute 3-body energy
                              vthrinter = vthrinter + 
     &                             twh_vthreebody(lfailure,rij,rjk
     &                             ,thetac,ntii,ntjj,ntkk)
                              if ( lfailure ) return
                           endif
                        endif
                     enddo
                  endif
               endif
            enddo
         endif
      enddo
      return
      end

      subroutine twh_engtotal_tb_eam(lfailure,vthrinter, lovrlap, ibox)
c     ******************************************************************
c     * calculates the three-body energy contribution for Embedded     *
c     * Atom Method                                                    *
c     *                                                                *
c     * Arguments set by subroutine:                                   *
c     *  vthrinter -- just the 3-body terms of the intermolecular      *
c     *     energy                                                     *
c     *                                                                *
c     * Arguments modified (but not initialized) by subroutine:        *
c     *  lovrlap -- Atomic overlap found?                              *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  ibox -- the box we are computing the energy in                *
c     *                                                                *
c     * Error conditions:                                              *
c     *   May halt program                                             *
c     *                                                                *
c     * split from original engtotal 08-12-2005 by MAW                 *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_GET_NBOXI
#define FUNCTION_VEMBED
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision vthrinter
      logical lfailure,lovrlap
      integer ibox
c     --- local variables
      logical lzero
      integer ichain,imolty,ntii,jchain,jmolty,ntjj,nchain
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij
      double precision rij,bij,bji,density

      vthrinter = 0.0d0

c     --- Embedded atom potential
c     --- reinitialize all the local densities to zero
      call twh_nchain(GLB_GET,nchain)
      do ichain = 1,nchain
         if ( twh_get_nboxi(ichain) .eq. ibox ) then
            call twh_eam_rho(GLB_SET,CRD_REAL,ichain,0.0d0)
         endif
      enddo

c        --- need to compute local densities for each atom and then
c        --- use that information to compute the embedding energy
      do ichain = 1,nchain-1
         if ( twh_get_nboxi(ichain) .eq. ibox ) then
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_ntype(GLB_GET,imolty,1,ntii)
            call twh_coordinates(GLB_GET,CRD_REAL,ichain,1
     &           ,rxui,ryui,rzui)
            do jchain = ichain+1,nchain
               if ( twh_get_nboxi(jchain) .eq. ibox 
     &              .and. ichain .ne. jchain) then
                  call twh_moltyp(GLB_GET,jchain,jmolty)
                  call twh_ntype(GLB_GET,jmolty,1,ntjj)
                  call twh_coordinates(GLB_GET,CRD_REAL,jchain,1
     &                 ,rxuj,ryuj,rzuj)
c                    --- compute i-j bond density
                  call twh_mimage(ibox,rxui,ryui,rzui,rxuj,ryuj,rzuj
     &                 ,rxuij,ryuij,rzuij)
                  rij = twh_distance(rxuij,ryuij,rzuij)
                  call twh_bondorder(lfailure,ntii,ntjj,rij,bij,bji
     &                 ,lzero)
                  if ( lfailure ) return
                  call twh_eam_rho(GLB_INCR,CRD_REAL,ichain,bij)
                  call twh_eam_rho(GLB_INCR,CRD_REAL,jchain,bji)
                  if ( lovrlap ) then
                     write(6,*) 'ENGTOTAL: undershoot in density'
                     lfailure = .true.
                     return
                  endif
               endif
            enddo
         endif
      enddo
c        --- add up the embedding energy
      do ichain = 1,nchain
         if ( twh_get_nboxi(ichain) .eq. ibox ) then
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_ntype(GLB_GET,imolty,1,ntii)
            call twh_eam_rho(GLB_GET,CRD_REAL,ichain,density)
            vthrinter = vthrinter + twh_vembed(lfailure,ntii,density)
            if ( lfailure ) return
            if ( lovrlap ) then
               write(6,*) 'ENGTOTAL: undershoot in embedding'
               write(6,*) 'density',density
               lfailure = .true.
               return
            endif
         endif
      enddo
      return
      end


      subroutine twh_engtotal_ewald(lfailure,vrecipsum,sself,correct
     &     , ibox )
c     ******************************************************************
c     * calculates the self and correction ewald terms.                *
c     * Note that coulomb potential scaling (scp_scale_coulomb) is not *
c     * yet implemented here                                           *
c     * Ewald corrections do not distinguish between inter- and intra- *
c     * molecular terms.  'Correct' term explicitly excludes real-space*
c     * terms                                                          *
c     *                                                                *
c     * Arguments set by subroutine:                                   *
c     *  vrecipsum -- reciprocal space ewald sum in units [K]          *
c     *  sself -- Coulombic 'self' in units [EE]                       *
c     *  correct -- Coulombic 'correction' in units [EE]               *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  ibox -- the box we are computing the energy in                *
c     *                                                                *
c     * Error conditions:                                              *
c     *  none                                                          *
c     *                                                                *
c     * split from original engtotal 08-12-2005 by MAW                 *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_GT
#define FUNCTION_DISTANCE
#define FUNCTION_EWALD_CORRECT
#define FUNCTION_EWALD_SELF
#define FUNCTION_GET_NBOXI
#define FUNCTION_LINCLUDE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ibox
      double precision vrecipsum, sself, correct
c     --- local variables
      integer ichain,imolty,ii,jj,nchain,numiunit
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij
      double precision rij,rcelect

c     --- retrieve constants
      call twh_rcelect(GLB_GET,ibox,rcelect)

      sself = 0.0d0
      correct = 0.0d0

      call twh_recipsum(lfailure,ibox,vrecipsum)
      if ( lfailure ) return
c        --- update self terms and correction terms
      call twh_nchain(GLB_GET,nchain)
      do ichain = 1,nchain
         if ( twh_get_nboxi(ichain) .eq. ibox) then
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,numiunit)
            do ii = 1,numiunit
               sself = sself + twh_ewald_self(ibox,imolty,ii)
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,ii
     &              ,rxui,ryui,rzui)
               do jj = ii+1,numiunit
                  if ( .not. twh_linclude(imolty,ii,jj) ) then
c                       --- compute the correction term assuming
c                       --- these atoms have no coulombic interations
c                       --- with each other.  For those with 1-4 
c                       --- interactions this assumption is fixed in 
c                       --- vtorsion
                     call twh_coordinates(GLB_GET,CRD_REAL,ichain,jj
     &                    ,rxuj,ryuj,rzuj)
                     call twh_mimage(ibox,rxui,ryui,rzui,rxuj,ryuj
     &                       ,rzuj,rxuij,ryuij,rzuij)
                     rij = twh_distance(rxuij,ryuij,rzuij)
                     if ( twh_cmp_gt(rcelect, rij) ) then 
                        correct = correct
     &                       + twh_ewald_correct(ibox,imolty,ii,jj,rij)
                     endif
                  endif
               enddo
            enddo
         endif
      enddo
      
      return
      end


      subroutine twh_engtotal_intrachain(lfailure,vintra, velectintra
     &     ,lovrlap,cutstyle, ibox, lmolecmove, lmovelist)
c     ******************************************************************
c     * calculates the intrachain interactions.                        *
c     * Arguments set by subroutine:                                   *
c     *  vintra -- the intramolecular (non-coulombic) non-bonded energy*
c     *  velectintra -- the intramolecular component of real space     *
c     *     coulombic potential                                        *
c     *                                                                *
c     * Arguments modified (but not initialized) by subroutine:        *
c     *  lovrlap -- Atomic overlap found?                              *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  cutstyle -- cutoff style.  See twh_vtwobody for definition    *
c     *  ibox -- the box we are computing the energy in                *
c     *  lmolecmove -- is true if this call is from a Monte Carlo move *
c     *  lmovelist -- true if we are performing a move and only need   *
c     *            the interactions for a subset of the molecules      *
c     *                                                                *
c     * Error conditions:                                              *
c     *   May exit early (sets lovrlap=true)                           *
c     *                                                                *
c     * split from original engtotal 08-11-2005 by MAW                 *
c     * last modified 04-22-2009 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NBOXI
#define FUNCTION_LINCLUDE
#define FUNCTION_ONEFIVETYPE
#define FUNCTION_VCOULOMB
#define FUNCTION_VONEFIVE
#define FUNCTION_VTWOBODY
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision vintra, velectintra
      logical lfailure,lovrlap,lmolecmove,lmovelist
      integer cutstyle,ibox

c     --- local variables
      logical lcompute,lmoveme,lcoulomb
      integer ichain,imolty,ii,ntii,jj,ntjj,ioftype,nchain,nunit
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij
      double precision rijsq,rij,dvdl,rminsq
c     --- foreign energies ignored for intramolecular evaluations
      double precision vf
      dimension vf(MAX_FOREIGN_LAMBDA)

c     --- retrieve constants
      call twh_lcoulomb(GLB_GET,lcoulomb)

      vintra = 0.0d0 
      velectintra = 0.0d0 

c     --- loop over all chains
      call twh_nchain(GLB_GET,nchain)
      do ichain = 1, nchain
c        --- check if i is in relevant box
         if ( twh_get_nboxi(ichain) .eq. ibox ) then
            if ( lmovelist ) then
               call twh_lmoveme(GLB_GET,ichain,lmoveme)
               if ( lmoveme ) then
                  lcompute = .true.
               else
                  lcompute = .false.
               endif
            else
               lcompute = .true.
            endif
            
            if ( lcompute ) then 
               call twh_moltyp(GLB_GET,ichain,imolty)
               call twh_nunit(GLB_GET,imolty,nunit)
               do ii = 1, nunit-1
                  call twh_ntype(GLB_GET,imolty,ii,ntii)
                  call twh_coordinates(GLB_GET,CRD_REAL,ichain,ii
     &                 ,rxui,ryui,rzui)
                  do jj = ii+1, nunit
                     if ( twh_linclude(imolty,ii,jj) ) then
                        call twh_ntype(GLB_GET,imolty,jj,ntjj)
                        call twh_coordinates(GLB_GET,CRD_REAL,ichain,jj
     &                       ,rxuj,ryuj,rzuj)
                        call twh_mimage(ibox,rxui,ryui,rzui,rxuj,ryuj
     &                       ,rzuj,rxuij,ryuij,rzuij)
                        rijsq = rxuij*rxuij + ryuij*ryuij 
     &                       + rzuij*rzuij
c                       --- get one-five type
                        ioftype = twh_onefivetype(imolty,ii,jj)
                        if ( ioftype .eq. 0 ) then
c                          --- compute nonbond intramolecular energy
                           vintra = vintra + twh_vtwobody(lfailure,rijsq
     &                          ,ntii,ntjj,cutstyle,lovrlap,vf,dvdl)
                           if ( lfailure ) return
                        else
c                          --- use special one-five terms for energy
                           vintra = vintra + twh_vonefive(lfailure,rijsq
     &                          ,ioftype,cutstyle,lovrlap)
                           if ( lfailure ) return
                        endif
                        if ( lovrlap ) then
                           if ( .not. lmolecmove ) then
                              write(6,*) 'overlap intra'
                              call twh_rminsq(GLB_GET,rminsq)
                              write(6,*) 'rijsq rminsq', rijsq,rminsq
                              write(6,*) 'chain,unit', ichain, ii
                              write(6,*) 'i-pos', rxui,ryui,rzui
                              write(6,*) 'jj', jj
                              write(6,*) 'j-pos',rxuj,ryuj,rzuj
                           endif
                           return
                        endif
                        if ( lcoulomb ) then
c                          --- for volume moves we only need the intra
c                          --- coulomb
                           rij = dsqrt(rijsq)
                           velectintra = velectintra + twh_vcoulomb(ibox
     &                          ,rij,imolty,ii,imolty,jj,.false.,vf
     &                          ,dvdl,lovrlap)
                           if ( lovrlap ) then
                              if ( .not. lmolecmove ) then
                                 write(6,*) 'overlap intra coulomb'
                                 call twh_rminsq(GLB_GET,rminsq)
                                 write(6,*) 'rijsq rminsq',rijsq,rminsq
                                 write(6,*) 'chain,unit', ichain, ii
                                 write(6,*) 'i-pos', rxui,ryui,rzui
                                 write(6,*) 'jj', jj
                                 write(6,*) 'j-pos',rxuj,ryuj,rzuj
                              endif
                              return
                           endif
                        endif
                     endif
                  enddo
               enddo
            endif
         endif
      enddo
      return
      end


      subroutine twh_engtotal_vibbendtors(lfailure,vvib,vbend,vtg
     &     ,vangleangle,vimpr,vbondbond,lovrlap,lprint,ibox)
c     ******************************************************************
c     * Calculate vibrational, bending and torsional energies          *
c     * Arguments set by subroutine:                                   *
c     *  vvib -- the bond vibration energy                             *
c     *  vbend -- the bond bending energy                              *
c     *  vtg -- the dihedral angle energy - including improper torsions*
c     *  vangleangle -- angle-angle potential                          *
c     *  vimpr -- the improper torsional energy                        *
c     *  vbondbond -- the IJ-KL bond-bond cross term                   *
c     *                                                                *
c     * Arguments modified (but not initialized) by subroutine:        *
c     *  lovrlap -- Atomic overlap found?                              *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  ibox -- the box we are computing the energy in                *
c     *  lprint -- print out all of the energies?                      *
c     *                                                                *
c     * Error conditions:                                              *
c     *   May exit early (sets lovrlap=true)                           *
c     *                                                                *
c     * split from original engtotal 08-11-2005 by MAW                 *
c     * last modified 07-24-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
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
#define FUNCTION_GET_NBOXI
#define FUNCTION_VANGANG
#define FUNCTION_VANGLE
#define FUNCTION_VBONBON
#define FUNCTION_VBOND
#define FUNCTION_VIMPROPER
#define FUNCTION_VTORSION
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision vvib, vbend, vtg, vangleangle, vimpr, vbondbond
      logical lfailure,lovrlap, lprint
      integer ibox

c     --- local variables
      logical lbadstereo
      integer ichain,imolty,ii,jjvib,ip1,ip2,it,jjben,ttype,ibead,iitor
      integer nchain,nunit
      dimension ttype(4)
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij
      double precision dist,duma,dumb,dumc
      double precision xtmp,ytmp,ztmp,length,angle
      dimension xtmp(4),ytmp(4),ztmp(4),length(4),angle(2)

      vtg = 0.0d0
      vimpr = 0.0d0
      vbend = 0.0d0
      vangleangle = 0.0d0
      vvib = 0.0d0

      call twh_nchain(GLB_GET,nchain)
      do ichain = 1, nchain
         call twh_moltyp(GLB_GET,ichain,imolty)
         call twh_nunit(GLB_GET,imolty,nunit)
c        --- check if i is in relevant box
         if ( twh_get_nboxi(ichain) .eq. ibox ) then
c              --- calculate all bonds vectors and lengths
c              --- calculate all stretching, bending, and torsional 
c              --- potentials that have an end-bead with an index less
c              --- than the current bead
c              --- bond vibrations
            do ii = 2, nunit
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,ii
     &              ,rxui,ryui,rzui)
               do jjvib = 1, twh_get_invib(imolty,ii)
                  call twh_ijvib(GLB_GET,imolty,ii,jjvib,ip1)
                  call twh_itvib(GLB_GET,imolty,ii,jjvib,it)
c                    --- increment vibrational energy
                  if ( ip1 .lt. ii ) then
                     call twh_coordinates(GLB_GET,CRD_REAL,ichain,ip1
     &                    ,rxuj,ryuj,rzuj)
c                       --- compute distance
c                       --- minimum image the bond vectors. this is 
c                       --- required in case the molecule network wraps
c                       --- completely through the periodic box
                     call twh_mimage ( ibox,rxui,ryui,rzui,rxuj,ryuj
     &                       ,rzuj,rxuij,ryuij,rzuij )
                     dist = twh_distance(rxuij,ryuij,rzuij)
                     vvib = vvib + twh_vbond(lfailure,it,dist,imolty
     &                       ,ii,ip1,lovrlap)
                     if ( lfailure ) return
                     if ( lovrlap ) then
                        write(6,*) 'overlap vibration'
                        write(6,*) 'chain,molecule,atom A,atom B'
     &                          ,ichain,imolty,ii,ip1
                        write(6,*) 'bond length:',dist
                        write(6,*) 'bond type:',it
                        return
                     endif
                  endif
               enddo
            enddo

c           --- molecule with bending angles
            do ii = 2, nunit
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,ii
     &                 ,xtmp(1),ytmp(1),ztmp(1))
               do jjben = 1, twh_get_inben(imolty,ii)
                  ip2 = twh_get_ijben3(imolty,ii,jjben)
                  if ( ip2 .lt. ii ) then
                     ip1 = twh_get_ijben2(imolty,ii,jjben)
c                       --- store the coordinates
                     call twh_coordinates(GLB_GET,CRD_REAL,ichain,ip1
     &                       ,xtmp(2),ytmp(2),ztmp(2))
                     call twh_coordinates(GLB_GET,CRD_REAL,ichain,ip2
     &                       ,xtmp(3),ytmp(3),ztmp(3))
                     vbend = vbend + twh_vangle(lfailure,lovrlap,ibox
     &                    ,imolty,ii,jjben,1,duma,dumb,dumc,xtmp,ytmp
     &                    ,ztmp)
                     if ( lfailure ) return
                     if ( lovrlap ) then
                        write(6,*) 'ENGTOTAL: overlap bending'
                        write(6,*) 'imolty,ii,ip1,ip2',imolty,ii
     &                          ,ip1,ip2
                        return
                     endif
                  endif
               enddo
            enddo

c           --- bond bond cross terms
            do ii = 1,nunit
               if ( twh_get_inbnbn(imolty,ii) .ne. 0 ) then
                  call twh_coordinates(GLB_GET,CRD_REAL,ichain,ii
     &                 ,rxui,ryui,rzui)
                  do jjben = 1,twh_get_inbnbn(imolty,ii)
                     ip1 = twh_get_ijbnbn1(imolty,ii,jjben)
                     ip2 = twh_get_ijbnbn2(imolty,ii,jjben)
c                    --- compute and store first bond length
                     call twh_coordinates(GLB_GET,CRD_REAL,ichain,ip1
     &                    ,rxuj,ryuj,rzuj)
                     call twh_mimage ( ibox,rxui,ryui,rzui,rxuj,ryuj
     &                       ,rzuj,rxuij,ryuij,rzuij )
                     length(1) = twh_distance(rxuij,ryuij,rzuij)
c                    --- compute and store second bond length
                     call twh_coordinates(GLB_GET,CRD_REAL,ichain,ip2
     &                    ,rxuj,ryuj,rzuj)
                     call twh_mimage ( ibox,rxui,ryui,rzui,rxuj,ryuj
     &                       ,rzuj,rxuij,ryuij,rzuij )
                     length(2) = twh_distance(rxuij,ryuij,rzuij)
c                    --- compute energy as a bond vibration term
                     vbondbond = vbondbond+twh_vbonbon(imolty,ii
     &                       ,jjben,length)
                  enddo
               endif
            enddo

c           --- regular torsions 
            do ii = 2, nunit
c              --- store the charge on atom ii
               ttype(1) = ii
               do iitor = 1, twh_get_intor(imolty,ii)
                  ttype(4) = twh_get_ijtor4(imolty,ii,iitor)
                  if ( ttype(4) .lt. ii ) then
c                       --- store unit numbers
                     ttype(2) = twh_get_ijtor2(imolty,ii,iitor)
                     ttype(3) = twh_get_ijtor3(imolty,ii,iitor)
c                       --- store atom positions
                     do ibead = 1,4
                        call twh_coordinates(GLB_GET,CRD_REAL,ichain
     &                          ,ttype(ibead)
     &                          ,xtmp(ibead),ytmp(ibead),ztmp(ibead))
                     enddo
c                       --- add in the torsional energy
                     vtg = vtg + twh_vtorsion(lfailure,ibox,1,imolty,ii
     &                    ,iitor,xtmp,ytmp,ztmp,duma,lovrlap)
                     if ( lfailure ) return
                     if ( lovrlap ) then
                        write(6,*) 'overlap torsion'
                        write(6,*) 'imolty,ttype',imolty,ttype
                        return
                     endif
                  endif
               enddo
            enddo

c           --- angle angle cross terms
            do ii = 1,nunit
               ttype(1) = ii
               do iitor = 1,twh_get_inaa(imolty,ii)
                  ttype(2) = twh_get_ijaa0(imolty,ii,iitor)
                  ttype(3) = twh_get_ijaa1(imolty,ii,iitor)
                  ttype(4) = twh_get_ijaa2(imolty,ii,iitor)
                  do ibead = 1,4
                     call twh_coordinates(GLB_GET,CRD_REAL,ichain
     &                    ,ttype(ibead)
     &                    ,xtmp(ibead),ytmp(ibead),ztmp(ibead))
                  enddo
c                    --- compute the angle-angle energy
                  vangleangle = vangleangle 
     &                 + twh_vangang(AA_COORD,ibox,imolty,ii,iitor
     &                 ,angle,xtmp,ytmp,ztmp)
               enddo
            enddo                   

c           --- improper torsions 
c           --- these are only stored on a single atom as they are 
c           --- not symmetrically defined
            do ii = 1,nunit
               ttype(1) = ii
               do iitor = 1,twh_get_inimprop(imolty,ii)
                  ttype(2) = twh_get_ijimprop2(imolty,ii,iitor)
                  ttype(3) = twh_get_ijimprop3(imolty,ii,iitor)
                  ttype(4) = twh_get_ijimprop4(imolty,ii,iitor)
c                    --- store the coordinates
                  do ibead = 1,4
                     call twh_coordinates(GLB_GET,CRD_REAL,ichain
     &                    ,ttype(ibead)
     &                       ,xtmp(ibead),ytmp(ibead),ztmp(ibead))
                  enddo
                  vimpr = vimpr + twh_vimproper(ibox,imolty,ii,iitor
     &                    ,xtmp,ytmp,ztmp,duma,lbadstereo)
                  if ( lbadstereo .and. lprint ) then
                     write(6,*) 'ENGTOTAL: Bad stereochemistry '
                     write(6,*) 'moltyp chain unit imptor'
     &                       ,imolty,ichain,ii,iitor
                  endif
               enddo
            enddo
         endif
      enddo
      return
      end

      subroutine twh_engtotal_ext(lfailure,vext,lovrlap,ibox,lmolecmove
     &     ,lmovelist)
c     ******************************************************************
c     * calculates the interaction energy with external surface        *
c     * Arguments set by subroutine:                                   *
c     *  vext -- the energy of external fields/walls                   *
c     *                                                                *
c     * Arguments modified (but not initialized) by subroutine:        *
c     *  lovrlap -- Atomic overlap found?                              *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  ibox -- the box we are computing the energy in                *
c     *  lmolecmove -- is true if this call is from a Monte Carlo move *
c     *  lmovelist -- true if we are performing a move and only need   *
c     *            the interactions for a subset of the molecules      *
c     *                                                                *
c     * Error conditions:                                              *
c     *   May exit early (sets lovrlap=true)                           *
c     *                                                                *
c     * split from original engtotal 08-11-2005 by MAW                 *
c     * last modified 10-29-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NBOXI
#define FUNCTION_VFIELD
#include "functions.h"

c     --- variables passed to/from the subroutine
      double precision vext
      logical lfailure,lovrlap, lmolecmove, lmovelist
      integer ibox

c     --- local variables
      logical lcompute,lmoveme
      integer ichain,imolty,ii,ntii,nchain,nunit,nfield
      double precision cord
      dimension cord(3)

      vext = 0.0d0

      call twh_nfield(GLB_GET,nfield)
      if ( nfield .ne. 0 ) then
         call twh_nchain(GLB_GET,nchain)
         do ichain = 1, nchain
            if ( twh_get_nboxi(ichain) .eq. ibox ) then
               if ( lmovelist ) then
                  call twh_lmoveme(GLB_GET,ichain,lmoveme)
                  if ( lmoveme ) then
                     lcompute = .true.
                  else
                     lcompute = .false.
                  endif
               else
                  lcompute = .true.
               endif
               
               if ( lcompute ) then
                  call twh_moltyp(GLB_GET,ichain,imolty)
                  call twh_nunit(GLB_GET,imolty,nunit)
                  do ii = 1, nunit
                     call twh_coordinates(GLB_GET,CRD_REAL,ichain,ii
     &                    ,cord(1),cord(2),cord(3))
                     call twh_ntype(GLB_GET,imolty,ii,ntii)
                     vext = vext + twh_vfield(lfailure,lovrlap,ibox
     &                    ,ntii,imolty,ichain,ii,cord)
                     if ( lfailure ) return
                     if (lovrlap) then
                        if ( .not. lmolecmove ) then
                           write(6,*) 'overlap with field'
                           write(6,*) 'chain,unit',ichain,ii
                           write(6,*) 'x,y,z',cord(1),cord(2),cord(3)
                        endif
                        return
                     endif
                  enddo
               endif
            endif
         enddo
      endif
      return
      end


      subroutine twh_engtotal_solv(lfailure,vsolvation,ibox )
c     ******************************************************************
c     * calculates the solvation energy for a configuration.           *
c     * Arguments set by subroutine:                                   *
c     *  vsolvation -- solvation energy term                           *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  ibox -- the box we are computing the energy in                *
c     *                                                                *
c     * Error conditions:                                              *
c     *   none                                                         *
c     *                                                                *
c     * split from original engtotal 08-11-2005 by MAW                 *
c     * last modified 11-01-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NBOXI
#define FUNCTION_VEEFONE
#define FUNCTION_VSASA
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure 
      double precision vsolvation
      integer ibox

c     --- local variables
      integer ichain,imolty,totunit,isolvtype,nchain
      double precision vnt
      integer chaininfo
      dimension chaininfo(1,3)

      vsolvation = 0.0d0

      call twh_isolvtype(GLB_GET,isolvtype)
      if ( isolvtype .eq. SOLV_SASA ) then
         call twh_nchain(GLB_GET,nchain)
         do ichain = 1, nchain  
            if ( twh_get_nboxi(ichain) .eq. ibox ) then
               call twh_moltyp(GLB_GET,ichain,imolty)
               call twh_nunit(GLB_GET,imolty,totunit)
               vsolvation = vsolvation + twh_vsasa(ichain,ibox
     &              ,CRD_REAL)
            endif
         enddo
      elseif( isolvtype .eq. SOLV_EEF1 ) then
         call twh_nchain(GLB_GET,nchain)
         do ichain = 1, nchain
            if ( twh_get_nboxi(ichain) .eq. ibox ) then
               call twh_moltyp(GLB_GET,ichain,imolty)
               call twh_nunit(GLB_GET,imolty,totunit)
               vsolvation = vsolvation + twh_veefone(ichain,ibox
     &              ,CRD_REAL)
            endif
         enddo
      elseif ( isolvtype .eq. SOLV_TRAMONTO ) then
         chaininfo(1,1)=-1
         chaininfo(1,2)=-1
         chaininfo(1,3)=-1

         call twh_fielddft(lfailure,0,ibox,1,chaininfo,vnt)
c         call twh_fielddft(lfailure,0,0,chaininfo,ibox,vnt)
         if ( lfailure ) return
         vsolvation = vsolvation + vnt
      endif
      return
      end
     

      subroutine twh_eng_print(nmcount, vvib, vbondbond, vbend, 
     &     vangleangle, vtg, vimpr, vintra, vinter, vtwointer, 
     &     vthrinter, velect, vrealelect, sself, correct, vrecipsum, 
     &     vext, vsolvation, vclassical, vtail, velectinter, 
     &     velectintra, vquantum, vforeign, dl_inter, dl_elect, 
     &     lmadelung, maddist, ibox)
c     ******************************************************************
c     * print the output                                               *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  nmcount -- number molecules in this box                       *
c     *  vvib -- the bond vibration energy                             *
c     *  vbondbond -- the IJ-KL bond-bond cross term                   *
c     *  vbend -- the bond bending energy                              *
c     *  vangleangle -- angle-angle potential                          *
c     *  vtg -- the dihedral angle energy - including improper torsions*
c     *  vimpr -- the improper torsional energy                        *
c     *  vintra -- the intramolecular (non-coulombic) non-bonded energy*
c     *  vinter -- the total intermolecular two-body, three-body and   *
c     *        tail correction non-coulombic energy                    *
c     *  vtwointer -- the two-body non-bond potential                  *
c     *  vthrinter -- just the 3-body terms of the intermolecular      *
c     *        energy                                                  *
c     *  velect -- the total coulombic energy - except 1-4 coulombics  *
c     *        which are included in the torsional potential           *
c     *  vrealelect -- real space component of coulombic potential     *
c     *  sself -- Coulombic 'self'                                     *
c     *  correct -- Coulombic 'correction'                             *
c     *  vrecipsum -- reciprocal space ewald sum                       *
c     *  vext -- the energy of external fields/walls                   *
c     *  vsolvation -- solvation energy term                           *
c     *  vclassical -- the total classical energy                      *
c     *  vtail -- just the tail correction terms of the intermolecular *
c     *        energy                                                  *
c     *  velectinter -- the intermolecular component of real space     *
c     *     coulombic potential                                        *
c     *  velectintra -- the intramolecular component of real space     *
c     *     coulombic potential                                        *
c     *  vquantum -- quantum mechanical energy                         *
c     *  vforeign -- classical potential at different values of lambda *
c     *     for Scaled Lennard-Jones potential only; see               *
c     *     scaled_potential.F                                         *
c     *  dl_inter -- derivative of two body potential with respect     *
c     *     to lambda_lj, calculated if scp_eval_derivative true       *
c     *  dl_elect -- derivative of inter-molecule coulomb potential    *
c     *     with respect to lambda_lj                                  *
c     *  lmadelung -- compute madelung constant?                       *
c     *  maddist -- distance for Madelung constant                     *
c     *  ibox -- the box we are computing the energy in                *
c     *                                                                *
c     * All energies and derivatives are in units of [K], with the     *
c     * exception of: vrealelect,velectintra,velectinter,sself,correct *
c     * which are in units of [EE].                                    *
c     *                                                                *
c     * Error conditions:                                              *
c     *  None                                                          *
c     *                                                                *
c     * split from original engtotal 08-11-2005 by MAW                 *
c     * last modified 06-10-2009 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_FOREIGN_LAMBDA
#define FUNCTION_GET_NCMT
#include "functions.h"

c     --- variables passed to/from the subroutine
      double precision vvib, vbondbond, vbend, vangleangle, vtg, vimpr  
      double precision vintra, vinter, vtwointer, vthrinter, velect  
      double precision vrealelect, sself, correct, vrecipsum, vext 
      double precision vsolvation, vclassical, maddist, vtail
     &     ,velectinter,velectintra,vquantum,dl_inter,dl_elect
      double precision vforeign(MAX_FOREIGN_LAMBDA)
      logical lmadelung
      integer nmcount, ibox, l
c     --- local variables
c     --- logical scalars
      logical scp_eval_derivative,ltailc,lcoulomb,lewald
c     --- integer scalars
      integer imolty,classical_potnum,nmolty,scp_num_foreign_lambda
      integer potentialstyle
c     --- double precision scalars
      double precision qqfact

c     --- retrieve constants
      call twh_constant_qqfact(GLB_GET,qqfact)
      call twh_nmolty(GLB_GET,nmolty)
      call twh_potentialstyle(GLB_GET,potentialstyle)
      call twh_ltailc(GLB_GET,ltailc)
      call twh_lcoulomb(GLB_GET,lcoulomb)
      call twh_lewald(GLB_GET,lewald)

      write(6,'(a27,1x,i10)'),'Total molecules in this box', nmcount
      do imolty = 1, nmolty
         write(6,8) 'Molecules of type',imolty,': '
     &        ,twh_get_ncmt(ibox,imolty)
 8       format(a17,1x,i3,1x,a2,i10)
      enddo
c     --- output classical energies
      if ( potentialstyle .eq. POT_INTERNAL ) then
         write(6,10) vvib+vbondbond,(vvib+vbondbond)/CON_KCALTOK
 10      format(' total vibration  ',f15.3,' [K] ',f15.5,' [kcal/mol]')
         write(6,12) vvib,vvib/CON_KCALTOK
 12      format('    regular       ',f15.3,' [K] ',f15.5,' [kcal/mol]')
         write(6,14) vbondbond,vbondbond/CON_KCALTOK
 14      format('    bond-bond(1-2)',f15.3,' [K] ',f15.5,' [kcal/mol]')
         write(6,20) vbend+vangleangle,(vbend+vangleangle)/CON_KCALTOK
 20      format(' total angle      ',f15.3,' [K] ',f15.5,' [kcal/mol]')
         write(6,22) vbend,vbend/CON_KCALTOK
 22      format('    regular       ',f15.3,' [K] ',f15.5,' [kcal/mol]')
         write(6,24) vangleangle,vangleangle/CON_KCALTOK
 24      format('    angle-angle   ',f15.3,' [K] ',f15.5,' [kcal/mol]')
         write(6,30) vtg+vimpr,(vtg+vimpr)/CON_KCALTOK
 30      format(' total torsion    ',f15.3,' [K] ',f15.5,' [kcal/mol]')
         write(6,32) vtg,vtg/CON_KCALTOK
 32      format('    regular       ',f15.3,' [K] ',f15.5,' [kcal/mol]')
         write(6,35) vimpr,vimpr/CON_KCALTOK
 35      format('    improper      ',f15.3,' [K] ',f15.5,' [kcal/mol]')
         write(6,40) vintra+vinter,(vintra+vinter)/CON_KCALTOK
c        --- van der waals
 40      format(' total nonbond    ',f15.3,' [K] ',f15.5,' [kcal/mol]')
         write(6,42) vintra,vintra/CON_KCALTOK
 42      format('    intramolecular',f15.3,' [K] ',f15.5,' [kcal/mol]')
         write(6,44) vtwointer,vtwointer/CON_KCALTOK
 44      format('    2-body nonbond',f15.3,' [K] ',f15.5,' [kcal/mol]')
         write(6,45) vthrinter,vthrinter/CON_KCALTOK
 45      format('    3-body nonbond',f15.3,' [K] ',f15.5,' [kcal/mol]')
         if ( ltailc ) write(6,46) vtail,vtail/CON_KCALTOK
 46      format('    tail correct. ',f15.3,' [K] ',f15.5,' [kcal/mol]')
         if ( lcoulomb ) then
            write(6,60) velect,velect/CON_KCALTOK
 60         format(' total coulombic  ',f15.3,' [K] ',f15.5
     &           ,' [kcal/mol]')
            if ( lewald ) then
               write(6,61) vrealelect*qqfact,vrealelect*qqfact
     &              /CON_KCALTOK
 61            format('    real space    ',f15.3,' [K] ',f15.5
     &              ,' [kcal/mol]')
               write(6,68) velectintra*qqfact,velectintra*qqfact
     &              /CON_KCALTOK
 68            format('      intramolec. ',f15.3,' [K] ',f15.5
     &              ,' [kcal/mol]')
               write(6,67) velectinter*qqfact,velectinter*qqfact
     &              /CON_KCALTOK
 67            format('      intermolec. ',f15.3,' [K] ',f15.5
     &              ,' [kcal/mol]')
               write(6,62) sself*qqfact,sself*qqfact/CON_KCALTOK
 62            format('    self          ',f15.3,' [K] ',f15.5
     &              ,' [kcal/mol]')
               write(6,64) correct*qqfact,correct*qqfact/CON_KCALTOK
 64            format('    correction    ',f15.3,' [K] ',f15.5
     &              ,' [kcal/mol]')
               write(6,66) vrecipsum,vrecipsum/CON_KCALTOK
 66            format('    recip sum     ',f15.3,' [K] ',f15.5
     &              ,' [kcal/mol]')
            endif

         endif
         write(6,70) vext,vext/CON_KCALTOK
 70      format(' external field  ',f16.3,' [K] ',f15.5,' [kcal/mol]')
         write(6,80) vsolvation,vsolvation/CON_KCALTOK
 80      format(' solvation       ',f16.3,' [K] ',f15.5,' [kcal/mol]')
         if ( dabs(vclassical) .gt. 1000.0d0 ) then
c           --- use less precision past the decimal point to avoid
c           --- printing noise
            write(6,89) vclassical,vclassical/CON_KCALTOK
         else
c           --- setting precision high for debugging...
            write(6,90) vclassical,vclassical/CON_KCALTOK
         endif
 89      format(' total classical',f23.4,' [K] ',f20.6,' [kcal/mol]')
 90      format(' total classical',f23.8,' [K] ',f20.10,' [kcal/mol]')

         call twh_scp_eval_derivative(GLB_GET,scp_eval_derivative)
         if (scp_eval_derivative) then
            write(6,91) dl_inter,dl_inter/CON_KCALTOK
            write(6,92) dl_elect,dl_elect/CON_KCALTOK
91          format(' dU/d[lambda_lj]',f17.3,' [K] ',f15.5,
     &         ' [kcal/mol]')
92          format(' dU/d[lambda_c] ',f17.3,' [K] ',f15.5,
     &         ' [kcal/mol]')
         endif

         call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
         if (scp_num_foreign_lambda .ne. 0) then
            write(6,'(1x,a,i1,a)')
     &           "Foreign Lambda energy calculations (box ",ibox,")"
            write(6,901) "index","lambda_lj","lambda_c",
     &            "Total Energy [K]", "Total Energy [kcal/mol]"
            do l = 1,scp_num_foreign_lambda
               write(6,902) l,scp_get_foreign_lambda(LAMBDA_LJ,l),
     &            scp_get_foreign_lambda(LAMBDA_C,l)
     &                 ,vforeign(l),vforeign(l)/CON_KCALTOK
            enddo
901         format(a6,a11,a11,a20,a25)
902         format(i6,f11.3,f11.3,f20.4,f20.4)
         endif

      endif
c     --- output quantum energies
      if ( potentialstyle .eq. POT_EXTERNAL ) then
         write(6,100) vquantum,vquantum/CON_KCALTOK
 100     format(' total quantum  ',f17.3,' [K] ',f15.5,' [kcal/mol]')
      endif
      
c        --- hack to print out the madelung constant in kj/mol
c        --- only want to do this for a system of sodium chloride
      if ( lmadelung ) then
         write(6,'(1x,a,f15.10)') 'Madelung Constant [unitless]: '
     &           ,(-velect*maddist)/(qqfact*dble(twh_get_ncmt(ibox,1)))
      endif
      call twh_cp_number(GLB_GET,classical_potnum)
      if ( classical_potnum .eq. CP_EMBEDDED_ATOM_METHOD ) then
c           --- print total energy in eV for EAM potential
         write(6,'(1x,a,1x,g12.6)') 'Total energy [eV]'
     &        ,vclassical/CON_EVTOK
      endif
      return
      end


      subroutine twh_engtotal_overlap(lfailure,lmolecmove,ichain,ii
     &     ,jchain,jj,rijsq,rxui,ryui,rzui,rxuj,ryuj,rzuj)
c     ******************************************************************
c     * handles the overlap error conditions and assigns the proper    *
c     * degree of problem for an overlap that occurs in engtotal       *
c     *                                                                *
c     * split from eng_twobody 07-17-2006 by M.G. Martin               *
c     * last modified 10-20-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETATOMNUM
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure,lmolecmove
      integer ichain,ii,jchain,jj
      double precision rijsq,rxui,ryui,rzui,rxuj,ryuj,rzuj,rminsq

      if ( .not. lmolecmove ) then
         write(6,*) 'intermolecular overlap'
         call twh_rminsq(GLB_GET,rminsq)
         write(6,*) 'rijsq rminsq',rijsq,rminsq
         write(6,*) 'rij:',dsqrt(rijsq),' rmin:',dsqrt(rminsq)
         write(6,*) 'chain, unit',ichain,ii
         write(6,*) 'i-pos', rxui,ryui,rzui
         write(6,*) 'jchain jj',jchain,jj
         write(6,*) 'j-pos',rxuj,ryuj,rzuj
         write(6,*) 'getatompos',twh_getatomnum(ichain,ii)
     &        ,twh_getatomnum(jchain,jj)
         lfailure = .true.
         return
      endif
      return
      end
