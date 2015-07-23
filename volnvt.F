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
      subroutine twh_volnvt(lfailure,movetype)
c     ******************************************************************
c     * movetype 1: makes an isotropic volume change over a range of   *
c     *   ln(vol) the maximum change is controled by rmvol and the     *
c     *   number of successful trial moves is stored in bsvol          *
c     * movetype 2: makes an adjustment of a single hmatrix element in *
c     *   one of the boxes, and a corresponding isotropic volume move  *
c     *   in the other box.  the maximum change is controled by rmcell *
c     *   and the number of successfule trial moves is stored in bscel *
c     *                                                                *
c     * originally written prior to 1995                               *
c     * last modified 08-02-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EXPON
#define FUNCTION_GET_NBOXI
#define FUNCTION_GET_PMVLPR
#define FUNCTION_GET_RMVOL
#define FUNCTION_MINBOXLENGTH
#define FUNCTION_RANDOM
#include "functions.h"
c     --- variables that are passed to/from the subroutines
      logical lfailure
      integer movetype
c     --- local variables
c     --- logical scalars
      logical lovrlap,lfound,ladjust,laccept,lreject,lbad,lperiodic
      logical lcoulomb,lewald
c     --- logical arrays
      logical lbonded
      dimension lbonded(2)
c     --- integer scalars
      integer ibox,imolty,imove,ifind,ilabel,ichain,iunit,keycode
      integer idim,jdim,movevector,ivector,jvector,nmolty,nboxpair
      integer boxiso,boxhmat,labiso,labhmat,icount,cpnumber,nchain
      integer ivalue,jvalue,nunit,nchbox,ncmt,potentialstyle
c     --- double precision scalars
      double precision volt,expdv,dx,dy,dz,vtotal,burn,vvibb
      double precision dele,vinter,vtail,vext,velect,vintra,vvib,vbend
      double precision vbendb,vtorb,vqtotal,boxvolume
      double precision vtor,vquantum,delecurrent
      double precision rantest,vsolv,vduma,rcmumax,test,rcmu
      double precision xtemp,ytemp,ztemp,pcheck,value,dvl,dvc
      double precision beta,temperature,xcmi,ycmi,zcmi,pmtest,rmcell
      double precision arbacmo,arbbcmo,arbccmo,dvalue,density
c     --- double precision arrays
      double precision vf(MAX_FOREIGN_LAMBDA)
      double precision hmold
      dimension hmold(2,3,3)
      double precision volo,voln,sfact
      dimension volo(2),voln(2),sfact(2)
      double precision vboxn,vintern,vtailn,vextn,velectn,vintran
      double precision vsolvn,vquantumn
      double precision vvibn,vbendn,vtorn
      dimension vvibn(2),vbendn(2),vtorn(2)
      dimension vboxn(2),vintern(2),vtailn(2),vextn(2),velectn(2)
      dimension vintran(2),vsolvn(2),vquantumn(2)
      double precision qbox
      dimension qbox(2)
      
c      write(6,*) 'start VOLNVT'
c     --- retrieve constants
      call twh_cp_number(GLB_GET,cpnumber)
      call twh_nchain(GLB_GET,nchain)
      call twh_potentialstyle(GLB_GET,potentialstyle)
      call twh_nboxpair(GLB_GET,nboxpair)
      call twh_lcoulomb(GLB_GET,lcoulomb)
      call twh_lewald(GLB_GET,lewald)

      if ( nboxpair .eq. 0 ) then
         write(6,*) 'VOLNVT: cannot NVT volume move with only 1 box'
         lfailure = .true.
         return
      endif

c     --- determine which boxes to perform the volume exchange
      rantest = twh_random()
      lfound = .false.
      do ifind = 1,nboxpair
         if ( movetype .eq. 1 ) then
            pcheck = twh_get_pmvlpr(ifind)
         elseif ( movetype .eq. 2 ) then
            call twh_pmcellpr(GLB_GET,ifind,pcheck)
         else
            write(6,*) 'VOLNVT: invalid movetype ',movetype
            lfailure = .true.
            return
         endif
         if ( .not. lfound .and. rantest .lt. pcheck ) then
            lfound = .true.
            imove = ifind
         endif
      enddo

      if ( movetype .eq. 1 ) then
c        --- add this to the list of volume moves attemped
         call twh_bnvol(GLB_INCR,imove,1)
      elseif ( movetype .eq. 2 ) then
c        --- determine a random vector
         ivector = idint(3.0d0*twh_random()) + 1
         jvector = idint(3.0d0*twh_random()) + 1
         movevector = 3*(ivector-1) + jvector
         call twh_bncell(GLB_INCR,imove,movevector,1.0d0)
      endif

c     --- store old box dimensions, energy, configuration etc
      do ilabel = 1,2
         call twh_pairbox(GLB_GET,imove,ilabel,ibox)
         do idim = 1,3
            do jdim = 1,3
               call twh_hmatrix(GLB_GET,ibox,idim,jdim
     &              ,hmold(ilabel,idim,jdim))
            enddo
         enddo
         call twh_boxvolume(GLB_GET,ibox,boxvolume)
         volo(ilabel) = boxvolume
         qbox(ilabel) = 0.25d0*twh_minboxlength(ibox)
c        --- store old k vectors and reciprocal sum
         if ( lewald ) call twh_recip(ibox,3)
      enddo

      do ichain = 1, nchain
         call twh_pairbox(GLB_GET,imove,1,ivalue)
         call twh_pairbox(GLB_GET,imove,2,jvalue)
         if ( twh_get_nboxi(ichain) .eq. ivalue .or.
     &        twh_get_nboxi(ichain) .eq. jvalue ) then
c           --- make sure this molecule is not bonded through the
c           --- periodic boundary
            call twh_lperiodic(GLB_GET,ichain,lperiodic)
            if ( lperiodic ) then
               write(6,*) 'VOLNVT: attempted move in a box that has'
     &              ,' a molecule bonded through the periodic boundary'
               write(6,*) 'Molecule type:',imolty
               write(6,*) 'Chain number:',ichain
               write(6,*) 'Box number:',twh_get_nboxi(ichain)
c              --- automatic rejection
               return
            endif
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,nunit)
c           --- store the old COM in cartesian coords
            call twh_com(GLB_STORE,CRD_R_TO_T,ichain,xcmi,ycmi,zcmi)
c           --- store the old COM in arb coords
            call twh_com(GLB_GET,CRD_REAL,ichain,xtemp,ytemp,ztemp)
            call twh_nboxi(GLB_GET,ichain,ibox)
            call twh_carttoarb(ibox,xtemp,ytemp,ztemp,arbacmo
     &           ,arbbcmo,arbccmo)
            call twh_arbcmo(GLB_SET,ichain,arbacmo,arbbcmo,arbccmo)
c           --- store the old cartesian coords for each unit
            do iunit = 1, nunit
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
               call twh_coordinates(GLB_SET,CRD_TEMP,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
            enddo
            if ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD ) then
c              --- store the old embedding density
               call twh_eam_rho(GLB_STORE,CRD_R_TO_T,ichain,density)
            endif
         endif
      enddo

c     --- calculate total volume (using label indexes)
      volt = volo(1) + volo(2)

c     --- initialize lreject, this is used to let us know if one of the
c     --- box dimensions, or box volumes has gone somewhere we will 
c     --- automatically reject.  If so then we skip all the energy 
c     --- computataions and just reset the box dimensions
      lreject = .false.
      if ( movetype .eq. 1 ) then
c        --- calculate new volume
         expdv = twh_expon(dlog(volo(1)/volo(2))
     &        + twh_get_rmvol(imove)*(2.0d0*twh_random()-1.0d0))
         voln(1)= expdv*volt/(1+expdv)
         voln(2)= volt-voln(1)
         if ( voln(1) .lt. 0.0d0 .or. voln(2) .lt. 0.0d0 ) then
            write(6,*) 'VOLNVT: a trial volume went negative'
            write(6,*) 'should not happen in the NVT version of volume'
            lfailure = .true.
            return
         endif
         sfact(1) = (voln(1)/volo(1))**(1.0d0/3.0d0)
         sfact(2) = (voln(2)/volo(2))**(1.0d0/3.0d0)

c        --- adjust the box dimensions of the two boxes
         do ilabel = 1,2
            call twh_pairbox(GLB_GET,imove,ilabel,ibox)
            do idim = 1,3
               do jdim = 1,3
                  call twh_hmatrix(GLB_SCALE,ibox,idim,jdim
     &                 ,sfact(ilabel))
               enddo
            enddo
         enddo
c        --- set the new hinverse and boxvolume for the two boxes
         do ilabel = 1,2
            call twh_pairbox(GLB_GET,imove,ilabel,ibox)
            call twh_inverthmatrix(ibox)
         enddo
      elseif ( movetype .eq. 2 ) then
c        --- we simply adjust a single box dimension for one of the 
c        --- boxes, and cause a corresponding volume change 
c        --- isotropically in the other box
c        --- randomly select which box to change first
         call twh_pmcellpt(GLB_GET,imove,pmtest)
         if ( twh_random() .lt. pmtest ) then
            labhmat = 1
            labiso = 2
         else
            labhmat = 2
            labiso = 1
         endif
c        --- get the box numbers
         call twh_pairbox(GLB_GET,imove,labhmat,boxhmat)
         call twh_pairbox(GLB_GET,imove,labiso,boxiso)
c        --- modify a single vector of the hmatrix in boxhmat
         call twh_rmcell(GLB_GET,boxhmat,ivector,jvector,rmcell)
         call twh_hmatrix(GLB_GET,boxhmat,ivector,jvector,dvalue)
         dvalue = dvalue + (2.0d0*twh_random()-1.0d0)*rmcell
         call twh_hmatrix(GLB_SET,boxhmat,ivector,jvector,dvalue)
c        --- compute the new volume for this box
         call twh_inverthmatrix(boxhmat)
c        --- make sure the new volume is still positive
         call twh_boxvolume(GLB_GET,boxhmat,boxvolume)
         if ( boxvolume .le. 0.0d0 ) then
            lreject = .true.
         else
c           --- compute the volume change and transfer it to the other
c           --- box
            voln(labhmat) = boxvolume
            voln(labiso) = volo(labiso) - (voln(labhmat)-volo(labhmat))
c           --- make sure this new volume is still positive
            if ( voln(labiso) .lt. 0.0d0 ) then
               lreject = .true.
            else
c              --- compute the isotropic scaling factor for hmatrix
               sfact(labiso) =(voln(labiso)/volo(labiso))**(1.0d0/3.0d0)
c              --- scale the hmatrix
               do idim = 1,3
                  do jdim = 1,3
                     call twh_hmatrix(GLB_SCALE,boxiso,idim,jdim
     &                    ,sfact(labiso)) 
                  enddo
               enddo
c              --- compute the new hinverse
               call twh_inverthmatrix(boxiso)
            endif
         endif
      endif

      do ilabel = 1,2
         call twh_pairbox(GLB_GET,imove,ilabel,ibox)
c        --- make sure the box dimensions are still valid for a right
c        --- handed coordinate system, and for rcut
         call twh_checkhmatrix(ibox,.false.,lbad)
         if ( lbad ) lreject = .true.
      enddo

      if ( lreject ) then
         if ( movetype .eq. 1 ) then 
            write(6,*) 'VOLNVT: move rejected, trial volume too low'
         else
            write(6,*) 'VOLNVT: cell move rejected, '
     &           ,'invalid trial hmatrix or volume'
         endif
c        --- restore old hmatrix
         do ilabel = 1,2
            call twh_pairbox(GLB_GET,imove,ilabel,ibox)
            do idim = 1,3
               do jdim = 1,3
                  call twh_hmatrix(GLB_SET,ibox,idim,jdim
     &                 ,hmold(ilabel,idim,jdim))
               enddo
            enddo
c           --- reset the hinverse and boxvolume
            call twh_inverthmatrix(ibox)
         enddo
         return
      endif

c     --- determine new positions of the molecules
c     --- calculate centre of mass and its displacement
      do ichain = 1, nchain
         call twh_nboxi(GLB_GET,ichain,ibox)
         ladjust = .false.
         call twh_pairbox(GLB_GET,imove,1,ivalue)
         call twh_pairbox(GLB_GET,imove,2,jvalue)
         if ( ibox .eq. ivalue ) then
            ilabel = 1
            ladjust = .true.
         elseif (ibox .eq. jvalue ) then
            ilabel = 2
            ladjust = .true.
         endif
         if ( ladjust ) then
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,nunit)
c           --- transform the old COM in arbitrary coordinates into
c           --- the new COM in cartesian coordinates
            call twh_arbcmo(GLB_GET,ichain,arbacmo,arbbcmo,arbccmo)
            call twh_arbtocart(ibox,arbacmo,arbbcmo,arbccmo
     &           ,xtemp,ytemp,ztemp)
c           --- compute the COM displacement
            call twh_com(GLB_GET,CRD_TEMP,ichain,xcmi,ycmi,zcmi)
            dx = xtemp - xcmi
            dy = ytemp - ycmi
            dz = ztemp - zcmi
c           --- increment the difference for the new COM 
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

c     --- compute the energy change
      dele = 0.0d0
      do ilabel = 1,2
         call twh_pairbox(GLB_GET,imove,ilabel,ibox)
c        --- set coulombic parameters for this new box length
         if ( lcoulomb ) then
            call twh_resetcoulomb(lfailure,ibox)
            if ( lfailure ) return
         endif
c        --- determine the maximum rcmu of any chain in this box
         rcmumax = 0.0d0
         call twh_nmolty(GLB_GET,nmolty)
         do imolty = 1,nmolty
            call twh_ncmt(GLB_GET,ibox,imolty,ncmt)
            do icount = 1,ncmt
               call twh_chainlist(GLB_GET,icount,ibox,imolty,ichain)
               call twh_rcmu(GLB_GET,ichain,rcmu)
               if ( rcmu .gt. rcmumax ) rcmumax = rcmu
            enddo
         enddo
c        --- reset quarter of minimum box if new minboxlength is smaller
c        --- than old
         test = twh_minboxlength(ibox)*0.25d0
         if ( test .lt. qbox(ilabel) ) qbox(ilabel) = test
c        --- determine whether the largest rcmu is greater than the 
c        --- smallest boxlength/4
         if ( rcmumax .gt. qbox(ilabel) ) then
c           --- we need to recompute intramolecular and bonded
c           --- move as the molecules could be interacting with 
c           --- themselves through the periodic images
            lbonded(ilabel) = .true.
         else
c           --- we do not need to recompute bonded interactions
            lbonded(ilabel) = .false.
         endif
c        --- calculate the new intermolecular energies
c        --- note that the intramolecular nonbonded term could change
c        --- if the molecule is large and changes from interacting 
c        --- through real space to through an image
         call twh_engtotal(lfailure,lovrlap,vtotal,vinter,vduma,vtail
     &        ,vintra,vvib,vbend,vtor,vext,velect,vsolv,vquantum,vf,dvl
     &        ,dvc,ibox,.true.,.false.,.false.,lbonded(ilabel))
         if ( lfailure ) return
         if ( lovrlap ) then
            laccept = .false.
            goto 500
         endif
         if ( potentialstyle .eq. POT_INTERNAL ) then
            vintern(ilabel) = vinter
            vtailn(ilabel)  = vtail
            vextn(ilabel)   = vext  
            velectn(ilabel) = velect
            vintran(ilabel) = vintra
            vsolvn(ilabel)  = vsolv
            if ( lbonded(ilabel) ) then
c              --- all energies were computed
               vtorn(ilabel) = vtor
               vbendn(ilabel) = vbend
               vvibn(ilabel) = vvib
               vboxn(ilabel) = vtotal
            else
c              --- did not recompute bonded energies
               call twh_boxvbond(GLB_GET,ibox,vvibb)
               call twh_boxvangle(GLB_GET,ibox,vbendb)
               call twh_boxvtorsion(GLB_GET,ibox,vtorb)
               vboxn(ilabel)   = vtotal + vtorb + vbendb + vvibb
            endif
            call scp_set_derivative(ilabel, LAMBDA_C, dvc)
            call scp_set_derivative(ilabel, LAMBDA_LJ, dvl)
            call scp_set_foreign_energies(ilabel, vf)
            
c           --- compute current energy change for this box
            call twh_boxvclassic(GLB_GET,ibox,vtotal)
            delecurrent = vboxn(ilabel) - vtotal
         elseif ( potentialstyle .eq. POT_EXTERNAL ) then
            vquantumn(ilabel) = vquantum
c           --- compute current energy change for this box
            call twh_boxvquantum(GLB_GET,ibox,vqtotal)
            delecurrent = vquantumn(ilabel) - vqtotal
         else
            write(6,*) 'VOLNVT: unknown potentialstyle'
            write(6,*) potentialstyle
            lfailure = .true.
            return
         endif
c        --- increment total dele
         dele = dele + delecurrent

         if ( movetype .eq. 1 ) then
c           --- store dU/dV in the averages
            keycode = AC_DU_DV
            value = delecurrent /(voln(ilabel)-volo(ilabel))
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,0
     &           ,value,1)
            if ( lfailure ) return
         endif
      enddo

      call twh_temperature(GLB_GET,temperature)
      if ( movetype .eq. 1 ) then
         do ilabel = 1,2
            call twh_pairbox(GLB_GET,imove,ilabel,ibox)
            call twh_nchbox(GLB_GET,ibox,nchbox)
            dele = dele - (dble(nchbox+1)*dlog(voln(ilabel)
     &           /volo(ilabel))*temperature)
         enddo
      elseif ( movetype .eq. 2 ) then
         do ilabel = 1,2
            call twh_pairbox(GLB_GET,imove,ilabel,ibox)
            call twh_nchbox(GLB_GET,ibox,nchbox)
            dele = dele - (dble(nchbox)*dlog(voln(ilabel)/volo(ilabel))
     &           *temperature )
         enddo
      endif

c     --- acceptance test
      laccept = .false.
      if ( dele .lt. 0.0d0 ) then
c        --- automatic acceptance
         laccept = .true.
c        --- burn a random number
         burn = twh_random()
      else
         call twh_beta(GLB_GET,beta)
         if ( twh_random() .lt. twh_expon(-beta*dele) ) then
            laccept = .true.
         endif
      endif

 500  if ( laccept ) then
c         write(6,*) 'volume move accepted',dele
c        --- move accepted
         if ( movetype .eq. 1 ) then
            call twh_bavol(GLB_INCR,imove,1)
         elseif ( movetype .eq. 2 ) then
            call twh_bacell(GLB_INCR,imove,movevector,1.0d0)
         endif
         do ilabel = 1,2
            call twh_pairbox(GLB_GET,imove,ilabel,ibox)
c           --- update the energies
            if ( potentialstyle .eq. POT_INTERNAL ) then
               call twh_boxvclassic(GLB_SET,ibox,vboxn(ilabel))
               call twh_boxvinter(GLB_SET,ibox,vintern(ilabel))
               call twh_boxvintra(GLB_SET,ibox,vintran(ilabel))
               call twh_boxvtail(GLB_SET,ibox,vtailn(ilabel))
               call twh_boxvfield(GLB_SET,ibox,vextn(ilabel))
               call twh_boxvelect(GLB_SET,ibox,velectn(ilabel))
               call twh_boxvsolvation(GLB_SET,ibox,vsolvn(ilabel))
               if ( lbonded(ilabel) ) then
                  call twh_boxvtorsion(GLB_SET,ibox,vtorn(ilabel))
                  call twh_boxvangle(GLB_SET,ibox,vbendn(ilabel))
                  call twh_boxvbond(GLB_SET,ibox,vvibn(ilabel))
               endif
            elseif ( potentialstyle .eq. POT_EXTERNAL ) then
               call twh_boxvquantum(GLB_SET,ibox,vquantumn(ilabel))
            else
               write(6,*) 'VOLNVT: unknown potentialstyle'
               write(6,*) potentialstyle
               lfailure = .true.
               return
            endif
         enddo
      else
c        --- volume move was rejected: restore old box lengths
         do ilabel = 1,2
            call twh_pairbox(GLB_GET,imove,ilabel,ibox)
            do idim = 1,3
               do jdim = 1,3
                  call twh_hmatrix(GLB_SET,ibox,idim,jdim
     &                 ,hmold(ilabel,idim,jdim))
               enddo
            enddo
c           --- reset old volume and hinverse
            call twh_inverthmatrix(ibox)

            if ( lcoulomb ) then
               call twh_resetcoulomb(lfailure,ibox)
               if ( lfailure ) return
               if ( lewald ) call twh_recip(ibox,4)
            endif
         enddo

c        --- restore old coordinates
         do ichain = 1, nchain
            call twh_pairbox(GLB_GET,imove,1,ivalue)
            call twh_pairbox(GLB_GET,imove,2,jvalue)
            if ( twh_get_nboxi(ichain) .eq. ivalue .or.
     &           twh_get_nboxi(ichain) .eq. jvalue ) then
               call twh_moltyp(GLB_GET,ichain,imolty)
               call twh_nunit(GLB_GET,imolty,nunit)
               call twh_com(GLB_STORE,CRD_T_TO_R,ichain,xcmi,ycmi,zcmi)
               do iunit = 1, nunit
                  call twh_coordinates(GLB_GET,CRD_TEMP,ichain,iunit
     &                 ,xtemp,ytemp,ztemp)
                  call twh_coordinates(GLB_SET,CRD_REAL,ichain,iunit
     &                 ,xtemp,ytemp,ztemp)
               enddo
               if ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD ) then
c                 --- restore the old embedding density
                  call twh_eam_rho(GLB_STORE,CRD_T_TO_R,ichain,density)
               endif
            endif
         enddo
      endif
c      write(6,*) 'end VOLNVT'

      return
      end

