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
      subroutine twh_volnpt(lfailure,movetype)
c     ******************************************************************
c     * movetype: 1                                                    *
c     *    makes an isotropic volume change under const. pressure      *
c     *    the maximum change is controlled by rmvol and the           *
c     *    number of successful trial moves is stored in bsvol         *
c     * movetype: 2                                                    *
c     *    changes a single element of the hmatrix which controls the  *
c     *    unit cell size and shape.  Move is performed on a single    *
c     *    simulation box under constant pressure                      *
c     *    this move was first described by Yashonath and Rao in       *
c     *    Mol. Phys. 54, 1, 245-251 (1985).                           *
c     *                                                                *
c     * originally written in 1995                                     *
c     * last modified 11-27-2011 by M.G.Martin                         *
c     ******************************************************************
      implicit none
c     --- debug toggle
#define DEBUG_VOLNPT 0
c     --- functions
#define FUNCTION_EXPON
#define FUNCTION_GET_NBOXI
#define FUNCTION_GET_PMVLPR
#define FUNCTION_GET_RMVOL
#define FUNCTION_MINBOXLENGTH
#define FUNCTION_RANDOM
#include "functions.h"

c     --- variables passed to/from the subroutine
      logical lfailure
      integer movetype

c     --- local variables
c     --- logical scalars
      logical ovrlap,lvol,lfound,laccept,lreject,lbad,lbonded,lperiodic
      logical lcoulomb,lewald
c     --- integer scalars
      integer ichain,iunit,imolty,boxvch,ibox,cpnumber,numboxes
      integer idim,jdim,movevector,ivector,jvector,icount,nchain
      integer nmolty,nunit,nchbox,ncmt,potentialstyle,keycode
      integer nfield,ifield,fieldbox
c     --- double precision scalars
      double precision volo,deltav,sfact,vinterb,vintrab,vvibb,vbendb
      double precision voln,vintern,vtailn,vextn,velectn,vintran,vtorb
      double precision vsolvn,burn,vquantumn,vextb,velectb,vsolvb
      double precision dx,dy,dz,v,dele,rantest,vduma,vqtotal
      double precision xtemp,ytemp,ztemp,pcheck,external_pressure
      double precision rcmumax,qbox,test,vvibn,vbendn,vtorn,value,rcmu
      double precision dvl,dvc,beta,temperature,xcmi,ycmi,zcmi,rmcell
      double precision arbacmo,arbbcmo,arbccmo,dvalue,density
c     --- double precision arrays
      double precision hmold
      dimension hmold(3,3)
      double precision vf
      dimension vf(MAX_FOREIGN_LAMBDA)

#if DEBUG_VOLNPT
      write(6,*) 'start VOLNPT'
#endif

c     --- retrieve constants
      call twh_cp_number(GLB_GET,cpnumber)
      call twh_nchain(GLB_GET,nchain)
      call twh_potentialstyle(GLB_GET,potentialstyle)
      call twh_lcoulomb(GLB_GET,lcoulomb)
      call twh_lewald(GLB_GET,lewald)
      call twh_numboxes(GLB_GET,numboxes)

c     --- Select a box at  random to change the volume of box
      lfound = .false.
      rantest = twh_random()
      do ibox = 1,numboxes
         if ( movetype .eq. 1 ) then
c           --- isotropic volume move, use pmvlpr
            pcheck = twh_get_pmvlpr(ibox)
         elseif ( movetype .eq. 2 ) then
c           --- unit cell perturbation, use pmcellpr
            call twh_pmcellpr(GLB_GET,ibox,pcheck)
         else
            write(6,*) 'VOLNPT: invalid movetype ',movetype
            lfailure = .true.
            return
         endif
         if ( .not. lfound .and. rantest .lt. pcheck ) then
            lfound = .true.
            boxvch = ibox
         endif
      enddo

      if ( movetype .eq. 1 ) then
c        --- increment the number of attempted 3d volume moves
         call twh_bnvol(GLB_INCR,boxvch,1)
      elseif ( movetype .eq. 2 ) then
c        --- determine which of the 9 quantities to alter
         ivector = idint(twh_random()*3.0d0) + 1
         jvector = idint(twh_random()*3.0d0) + 1
         movevector = 3*(ivector-1) + jvector
c        --- increment the number of attemped unit cell changes
         call twh_bncell(GLB_INCR,boxvch,movevector,1.0d0)
      endif

c     --- check for external fields in this box
      call twh_nfield(GLB_GET,nfield)
      do ifield = 1,nfield
         call twh_fieldint(GLB_GET,ifield,1,fieldbox)
         if ( fieldbox .eq. boxvch ) then
c           --- automatic rejection, print warning
            write(6,*) 'Cannot perform a volume change move in a box'
     &           ,' that has an external field'
            write(6,*) 'Automatic rejection in box:',boxvch
            return
         endif
      enddo

c     --- store old box lengths, volume, configuration etc
c     --- in the box of volume change
      do idim = 1,3
         do jdim = 1,3
            call twh_hmatrix(GLB_GET,boxvch,idim,jdim,hmold(idim,jdim))
         enddo
      enddo
      call twh_boxvolume(GLB_GET,boxvch,volo)
c     --- set the quarter of the minimum box length for old box
      qbox = twh_minboxlength(boxvch)*0.25d0

      if ( lewald ) call twh_recip(boxvch,3)
      do ichain = 1, nchain
c        ---- Check if the chain i is in the correct box
         if ( twh_get_nboxi(ichain). eq. boxvch) then
c           --- make sure the molecule is not periodically bonded
            call twh_lperiodic(GLB_GET,ichain,lperiodic)
            if ( lperiodic ) then
               write(6,*) 'VOLNPT: move attemped in box that has a'
     &              ,' periodic bonded molecule'
               write(6,*) 'Molecule type:',imolty
               write(6,*) 'Chain number:',ichain
               write(6,*) 'Box:',ibox
c              --- automatic rejection
               return
            endif
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,nunit)
c           --- store the current "old" COM in cartesian coords
            call twh_com(GLB_STORE,CRD_R_TO_T,ichain,xtemp,ytemp,ztemp)
c           --- store the current "old" COM in arb coords
            call twh_com(GLB_GET,CRD_REAL,ichain,xtemp,ytemp,ztemp)
            call twh_carttoarb(boxvch,xtemp,ytemp,ztemp,arbacmo
     &           ,arbbcmo,arbccmo)
            call twh_arbcmo(GLB_SET,ichain,arbacmo,arbbcmo,arbccmo)
c           --- store all of the old unit coords in cart
            do iunit = 1, nunit
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
               call twh_coordinates(GLB_SET,CRD_TEMP,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
            enddo
            if ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD ) then
c              --- store the current eam_rho in the temp
               call twh_eam_rho(GLB_STORE,CRD_R_TO_T,ichain,density)
            endif
         endif
      enddo

      lreject = .false.
      if ( movetype .eq. 1 ) then
c        --- calculate new volume         
         voln = volo 
     &        + twh_get_rmvol(boxvch)*( 2.0d0*twh_random() - 1.0d0 )
c        --- make sure the new volume is positive
         if ( voln .lt. 0.0d0 ) then
            lreject = .true.
         else
c           --- isotropic volume change
            sfact = (voln/volo)**(1.0d0/3.0d0)
c           --- set new box dimensions
            do idim = 1,3
               do jdim = 1,3
                  call twh_hmatrix(GLB_SCALE,boxvch,idim,jdim,sfact)
               enddo
            enddo
         endif
      elseif ( movetype .eq. 2 ) then
c        --- modify a single vector of the hmatrix
         call twh_rmcell(GLB_GET,boxvch,ivector,jvector,rmcell)
         call twh_hmatrix(GLB_GET,boxvch,ivector,jvector,dvalue)
         dvalue = dvalue + (2.0d0*twh_random()-1.0d0)*rmcell
         call twh_hmatrix(GLB_SET,boxvch,ivector,jvector,dvalue)
      endif
c     --- check the new box dimensions
      call twh_checkhmatrix(boxvch,.false.,lbad)
      if ( lbad ) lreject = .true.

      if ( .not. lreject ) then
c        --- set the new boxvolume and hinverse
         call twh_inverthmatrix(boxvch)
c        --- make sure the volume is still positive
         call twh_boxvolume(GLB_GET,boxvch,voln)
         if ( voln .lt. 0.0d0 ) then
            lreject = .true.
         endif
      endif

      if ( lreject ) then
         write(6,*) 'VOLNPT: move rejected, negative volume or bad '
     &        ,'trial hmatrix'
c        --- restore old hmatrix
         do idim = 1,3
            do jdim = 1,3
               call twh_hmatrix(GLB_SET,boxvch,idim,jdim
     &              ,hmold(idim,jdim))
            enddo
         enddo
c        --- restore hinverse and boxvolume
         call twh_inverthmatrix(boxvch)
         return
      endif

c     --- determine the maximum rcmu of any chain in this box
      rcmumax = 0.0d0
      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_ncmt(GLB_GET,boxvch,imolty,ncmt)
         do icount = 1,ncmt
            call twh_chainlist(GLB_GET,icount,boxvch,imolty,ichain)
            call twh_rcmu(GLB_GET,ichain,rcmu)
            if ( rcmu .gt. rcmumax ) rcmumax = rcmu
         enddo
      enddo
c     --- reset quarter of minimum box if new minboxlength is smaller 
c     --- than old
      test = twh_minboxlength(boxvch)*0.25d0
      if ( test .lt. qbox ) qbox = test
c     --- determine whether the largest rcmu is greater than the 
c     --- smallest boxlength/4
      if ( rcmumax .gt. qbox ) then
c        --- we need to recompute intramolecular and bonded for volume
c        --- move as the molecules could be interacting with themselves
c        --- through the periodic images
         lbonded = .true.
      else
c        --- we do not need to recompute bonded interactions
         lbonded = .false.
      endif

c     --- determine new positions of the molecules
c     --- calculate centre of mass and its displacement
      do ichain = 1, nchain
         call twh_moltyp(GLB_GET,ichain,imolty)
         call twh_nunit(GLB_GET,imolty,nunit)
         if ( twh_get_nboxi(ichain) .eq. boxvch) then
c           --- transform the old COM in arbitrary coordinates into 
c           --- the new COM in cartesian coordinates
            call twh_arbcmo(GLB_GET,ichain,arbacmo,arbbcmo,arbccmo)
            call twh_arbtocart(boxvch,arbacmo,arbbcmo,arbccmo
     &           ,xtemp,ytemp,ztemp)
c           --- compute the shift from the old to the new coords
            call twh_com(GLB_GET,CRD_TEMP,ichain,xcmi,ycmi,zcmi)
            dx = xtemp - xcmi
            dy = ytemp - ycmi
            dz = ztemp - zcmi
c           --- store the new COM
            call twh_com(GLB_INCR,CRD_REAL,ichain,dx,dy,dz)
c           --- adjust all of the unit coordinates
            do iunit = 1, nunit
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
               xtemp = xtemp + dx
               ytemp = ytemp + dy
               ztemp = ztemp + dz
               call twh_coordinates(GLB_SET,CRD_REAL,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
            enddo
         endif
      enddo

      lvol = .true.
c     --- set new coulombic parameters for the new box length
      if ( lcoulomb ) then
         call twh_resetcoulomb(lfailure,boxvch)
         if ( lfailure ) return
      endif

      call twh_engtotal(lfailure,ovrlap,v,vintern,vduma,vtailn,vintran
     &     ,vvibn,vbendn,vtorn,vextn,velectn,vsolvn,vquantumn
     &     ,vf,dvl,dvc,boxvch,lvol,.false.,.false.,lbonded)
      if ( lfailure ) return
c     --- initialize laccept
      laccept = .false.

      if ( ovrlap ) then
c         write(6,*) 'move rejected due to overlap in VOLNPT'
      else
c        -- store the new total box energy
         if ( potentialstyle .eq. POT_INTERNAL ) then
            call twh_boxvinter(GLB_GET,boxvch,vinterb)
            call twh_boxvintra(GLB_GET,boxvch,vintrab)
            call twh_boxvfield(GLB_GET,boxvch,vextb)
            call twh_boxvelect(GLB_GET,boxvch,velectb)
            call twh_boxvsolvation(GLB_GET,boxvch,vsolvb)
            deltav = (vintern - vinterb) + (vextn  - vextb)
     &           + (velectn - velectb) + (vintran - vintrab)
     &           + (vsolvn - vsolvb)
            if ( lbonded ) then
c              --- add in the bonded terms energy change
               call twh_boxvbond(GLB_GET,boxvch,vvibb)
               call twh_boxvangle(GLB_GET,boxvch,vbendb)
               call twh_boxvtorsion(GLB_GET,boxvch,vtorb)
               deltav = deltav + (vtorn - vtorb) + (vbendn - vbendb)
     &              + (vvibn - vvibb)
            endif
         elseif ( potentialstyle .eq. POT_EXTERNAL ) then
            call twh_boxvquantum(GLB_GET,boxvch,vqtotal)
            deltav = vquantumn - vqtotal
         else
            write(6,*) 'VOLNPT: unknown potentialstyle'
            write(6,*) potentialstyle
            lfailure = .true.
            return
         endif

         if ( movetype .eq. 1 ) then
c           --- store dU/dV in the averages
            keycode = AC_DU_DV
            value = deltav/(voln-volo)
            call twh_averages(lfailure,AVG_INCREMENT,keycode,boxvch,0
     &           ,value,1)
            if ( lfailure ) return
         endif
c        --- add in the number density and pressure terms
         call twh_express(GLB_GET,external_pressure)
         call twh_nchbox(GLB_GET,boxvch,nchbox)
         call twh_temperature(GLB_GET,temperature)
         dele = deltav + external_pressure*(voln - volo) 
     &        - (dble(nchbox)*dlog(voln/volo)*temperature )
c        --- acceptance test 
         if ( dele .lt. 0.0d0 ) then
c           --- accept move
            laccept = .true.
c           --- burn a random number
            burn = twh_random()
         else
            call twh_beta(GLB_GET,beta)
            if ( twh_random() .lt. twh_expon(-(beta*dele)) ) then
               laccept = .true.
            endif
         endif
      endif
#if DEBUG_VOLNPT      
      write(6,*) 'dele, laccept',dele,laccept
#endif
      if ( laccept ) then
         if ( movetype .eq. 1 ) then
c           --- increment the number of accepted volume moves
            call twh_bavol(GLB_INCR,boxvch,1)
         elseif ( movetype .eq. 2 ) then
c           --- increment the number of accepted unit cell moves
            call twh_bacell(GLB_INCR,boxvch,movevector,1.0d0)
         endif
c        --- store the new energies
         if ( potentialstyle .eq. POT_INTERNAL ) then
            call twh_boxvclassic(GLB_INCR,boxvch,deltav)
            call twh_boxvinter(GLB_SET,boxvch,vintern)
            call twh_boxvintra(GLB_SET,boxvch,vintran)
            call twh_boxvtail(GLB_SET,boxvch,vtailn)
            call twh_boxvfield(GLB_SET,boxvch,vextn)
            call twh_boxvelect(GLB_SET,boxvch,velectn)
            call twh_boxvsolvation(GLB_SET,boxvch,vsolvn)
            if ( lbonded ) then
               call twh_boxvtorsion(GLB_SET,boxvch,vtorn)
               call twh_boxvangle(GLB_SET,boxvch,vbendn)
               call twh_boxvbond(GLB_SET,boxvch,vvibn)
            endif
            call scp_set_derivative(boxvch, LAMBDA_C, dvc)
            call scp_set_derivative(boxvch, LAMBDA_LJ, dvl)
            call scp_set_foreign_energies(boxvch, vf)

         elseif ( potentialstyle .eq. POT_EXTERNAL ) then
            call twh_boxvquantum(GLB_INCR,boxvch,deltav)
         else
            write(6,*) 'VOLNPT: unknown potentialstyle'
            write(6,*) potentialstyle
            lfailure = .true.
            return
         endif
      else
c        --- move was rejected
c        --- restore old box lengths, energy, configration etc
         do idim = 1,3
            do jdim = 1,3
               call twh_hmatrix(GLB_SET,boxvch,idim,jdim
     &              ,hmold(idim,jdim))
            enddo
         enddo
c        --- restore hinverse and boxvolume
         call twh_inverthmatrix(boxvch)
         if ( lcoulomb ) then
            call twh_resetcoulomb(lfailure,boxvch)
            if ( lfailure ) return
            if ( lewald ) call twh_recip(boxvch,4)
         endif
c        --- restore coordinates
         do ichain = 1, nchain
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,nunit)
            if ( twh_get_nboxi(ichain) .eq. boxvch) then
               call twh_com(GLB_STORE,CRD_T_TO_R,ichain,xcmi,ycmi,zcmi)
               do iunit = 1, nunit
                  call twh_coordinates(GLB_GET,CRD_TEMP,ichain,iunit
     &                 ,xtemp,ytemp,ztemp)
                  call twh_coordinates(GLB_SET,CRD_REAL,ichain,iunit
     &                 ,xtemp,ytemp,ztemp)
               enddo
               if ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD ) then
c                 --- restore the current eam_rho
                  call twh_eam_rho(GLB_STORE,CRD_T_TO_R,ichain,density)
               endif
            endif
         enddo
      endif

c      write(6,*) 'end VOLNPT'

      return
      end
