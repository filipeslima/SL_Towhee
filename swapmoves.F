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
      subroutine twh_swapmoves(lfailure,swaptype,inmoltype,inbox)
c     ******************************************************************
c     * removes a molecule from one box and inserts it into the other  *
c     * using CBMC insertion techniques.                               *
c     * for most of the swaptypes the inmoltype and inbox variables are*
c     * not used.  The exception is swaptype 7 where we are measuring  *
c     * the chemical potential of type inmoltype in box inbox          *
c     * valid swaptype codes are listed here                           *
c     * SWAP_CB_INTER                                                  *
c     *             inter-box configurational-bias molecule transfer   *
c     *             move.  takes a molecule out of one box and inserts *
c     *             it into another box using CBMC                     *
c     * SWAP_CB_INTRA                                                  *
c     *             intra-box configurational-bias molecule transfer   *
c     *             move.  takes a molecule out of one box and inserts *
c     *             it back into the same box using CBMC               *
c     * SWAP_AVB_ONE                                                   *
c     *             aggregation volume bias move type 1.  takes a      *
c     *             molecule in one box and inserts it back into the   *
c     *             same box, using a volume reference related to the  *
c     *             distance from another selected atom in that box    *
c     *             see B. Chen; J.I. Siepmann; J. Phys. Chem. B 104,  *
c     *             8725-8734 (2000)                                   *
c     * SWAP_AVB_TWO                                                   *
c     *             aggregation volume bias move type 2.  takes a      *
c     *             molecule, decides whether to do an in->out or an   *
c     *             out-> in move and then moves another molecule      *
c     *             between those two regions                          *
c     *             see B. Chen; J.I. Siepmann; J. Phys. Chem. B 105,  *
c     *             11275-11282                                        *
c     * SWAP_AVB_THREE                                                 *
c     *             aggregation volume bias move type 3.  takes a      *
c     *             molecule, finds another molecule which is more than*
c     *             avb3rad from it and then moves a third molecule    *
c     *             between the inner and outer regions of the first 2 *
c     *             see B. Chen; J.I. Siepmann; J. Phys. Chem. B 105,  *
c     *             11275-11282                                        *
c     * SWAP_CB_GC                                                     *
c     *             inter-box configurational-bias molecule transfer   *
c     *             move used in the grand canonical ensemble.  takes  *
c     *             a molecule out of the ideal gas reservoir box and  *
c     *             inserts it into the system box using CBMC          *
c     * SWAP_RB_INTER                                                  *
c     *             inter-box rotational-bias molecule transfer move.  *
c     *             takes a molecule out of one box and inserts it into*
c     *             another box using the same configuration, but with *
c     *             a random orientation                               *
c     * SWAP_CHEMPOT_BOX                                               *
c     *             not a Monte Carlo move - instead it uses           *
c     *             the machinery here to compute the chemical         *
c     *             potential for molecule type inmoltype in box inbox *
c     * SWAP_CHEMPOT_ISO                                               *
c     *             computing the chemical potential of an isolated    *
c     *             molecule                                           *
c     *                                                                *
c     * rewritten from older subroutines on 09-19-1997 by M.G. Martin  *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_RANDOM
#include "functions.h"
#define DEBUG_SWAPMOVES 0
c     --- variables passed to/from the subroutine
      logical lfailure
      integer swaptype,inmoltype,inbox
c     --- local variables
c     --- logical scalars
      logical lnew,lempty,ldone,lfound,laccept,lgetcom,ltailc,lembed
      logical scp_eval_derivative,lewald,ldualcutoff,tmmc_flag
      logical lintra,lbonded,lgetbox,lpartial,lcompute
c     --- logical arrays
      logical ltermtest
      dimension ltermtest(2)
c     --- integer scalars
      integer imol,idum,imove,mcount,cordstop,keycode
      integer iutry,isteps,jrem,totunits,iunit,pointp,ipoint
      integer imt,jmt,jmolty,junit,imolty,ibox,ncmt,numboxes
      integer regstart,regend,avbtype
      integer nboxmove,icordflag,num_in,num_out,index,ngood
      integer nchain,nmolty,scp_num_foreign_lambda,potentialstyle
      integer activebox,nchainmove,testcord,globalcord,cpnumber,nboxpair
      integer n_tmmc_min,n_tmmc_max,c_element,macro_variable
      integer nch_nb_one
c     --- integer arrays
      integer boxtest,chaintest
      dimension boxtest(2),chaintest(2)
      integer goodbox
      dimension goodbox(MAXBOX)
      integer chaininfo
      dimension chaininfo(1,3)
c     --- double precision scalars
      double precision rmol,vol_in,vol_out,vtest
      double precision wnlog,wolog,wdlog,wratio,rho,arg,burn
      double precision vdum,xtemp,ytemp,ztemp,wmin,coru,scalecut
      double precision onepi,dpmin,testrad,density,pin,pout,pmtest
      double precision ddvdl_lj,ddvdl_c,ranbox
      double precision scalelolog,boxvolume,debroglie,uvtfactor
      double precision d_tmmc_weight,log_e,dvalue,tmmc_one,tmmc_zero
      double precision tmmc_negone
      double precision ectotal,beta
c     --- double precision arrays
      double precision dvctot
      dimension dvctot(2)
      double precision dvinter
      dimension dvinter(2)
      double precision dvelect
      dimension dvelect(2)
      double precision dvintra
      dimension dvintra(2)
      double precision dvvib
      dimension dvvib(2)
      double precision dvben
      dimension dvben(2)
      double precision dvtor
      dimension dvtor(2)
      double precision dvsolv
      dimension dvsolv(2)
      double precision dvext
      dimension dvext(2)
      double precision dvquantum
      dimension dvquantum(2)
      double precision dpscaletest
      dimension dpscaletest(2)
      double precision vmintest,v1tottest,v1exttest,v1inttest
      dimension vmintest(2),v1tottest(2),v1exttest(2),v1inttest(2)
      double precision v1elctest,v1ewdtest,vtailtest
      dimension v1elctest(2),v1ewdtest(2),vtailtest(2)
      double precision xtarget,ytarget,ztarget
      dimension xtarget(2),ytarget(2),ztarget(2)
      double precision wrosen,w1test,waddtest,voltest
      dimension wrosen(2),w1test(2),waddtest(2),voltest(2)
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

      log_e = 0.4342944819032520D0
c     --- retrieve constants
      call twh_constant_pi(GLB_GET,onepi)
      call twh_nmolty(GLB_GET,nmolty)
      call twh_potentialstyle(GLB_GET,potentialstyle)
      call twh_ltailc(GLB_GET,ltailc)
      call twh_numboxes(GLB_GET,numboxes)
      call twh_tmmc_flag(GLB_GET,tmmc_flag)

#if DEBUG_SWAPMOVES
      write(6,*) 'START SWAPMOVES',swaptype
#endif
c     --- retrieve constants
      call twh_lewald(GLB_GET,lewald)
      call twh_scalecut(GLB_GET,scalecut)

c     --- set variables that depend only upon the swaptype
c     --- set up trickery for DUAL, TRAMONTO, or multibody potentials
      lpartial = .false.
      call twh_ldualcutoff(GLB_GET,ldualcutoff)
      if ( ldualcutoff ) lpartial = .true.
#if USETRAMONTO
      lpartial = .true.
#endif
      lembed = .false.
      if ( potentialstyle .eq. POT_INTERNAL ) then
         call twh_cp_number(GLB_GET,cpnumber)
         if ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD .or.
     &        cpnumber .eq. CP_EAM_PAIR_ONLY .or.
     &        cpnumber .eq. CP_STILLINGER_WEBER .or.
     &        cpnumber .eq. CP_SW_PAIR_ONLY ) then
            lpartial = .true.
         endif
         if ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD ) lembed = .true.
      elseif ( potentialstyle .eq. POT_EXTERNAL ) then
         lpartial = .true.
      else
         write(6,*) 'SWAPMOVES: unknown potentialstyle'
         write(6,*) potentialstyle
         lfailure = .true.
         return
      endif
c     --- set lintra and lbonded for rotational-bias
      if ( swaptype .eq. SWAP_RB_INTER ) then
c        --- we need the intramolecular and bonded energies
         lintra = .true.
         lbonded = .true.
      else
c        --- intramolecular taken care of in engatom
         lintra = .false.
         lbonded = .false.
      endif

      lgetcom = .true.
      lempty = .false.
      lfound = .false.
      lgetbox = .false.
      vtest = 0.0d0
      d_tmmc_weight = 0.0d0
      do icordflag = CRD_OLD,CRD_NEW,(CRD_NEW-CRD_OLD)
c        --- initialize boxtest to 0
         boxtest(icordflag) = 0
c        --- initialize scaling factor
         dpscaletest(icordflag) = 0.0d0
c        --- initialize the rosenbluth weight correction term
         waddtest(icordflag) = 1.0d0
c        --- initialize tail correction energy
         vtailtest(icordflag) = 0.0d0
c        --- initialize growth overlap logical
         ltermtest(icordflag) = .false.
      enddo

c     --- select a box
      if ( swaptype .eq. SWAP_RB_INTER
     &     .or. swaptype .eq. SWAP_CB_INTER ) then
c        --- two box swap moves
         nboxmove = 2
         call twh_nboxpair(GLB_GET,nboxpair)
         if ( nboxpair .eq. 0 ) then
            write(6,*) 'SWAPMOVES: tried to swap with only 1 box'
            lfailure = .true.
            return
         endif
c        --- choose a pair of boxes according to the probabilities
         rmol = twh_random()
         do imove = 1,nboxpair
            if ( swaptype .eq. SWAP_RB_INTER ) then
               call twh_pm2rbswpr(GLB_GET,imove,pmtest)
            elseif ( swaptype .eq. SWAP_CB_INTER ) then
               call twh_pm2cbswpr(GLB_GET,imove,pmtest)
            endif
            if ( .not. lfound .and. rmol .lt. pmtest ) then
c              --- assign the two boxes
               lfound = .true.
c              --- randomly choose the target boxes
               if ( twh_random() .lt. 0.5d0 ) then
                  call twh_pairbox(GLB_GET,imove,1,boxtest(2))
                  call twh_pairbox(GLB_GET,imove,2,boxtest(1))
               else
                  call twh_pairbox(GLB_GET,imove,2,boxtest(2))
                  call twh_pairbox(GLB_GET,imove,1,boxtest(1))
               endif     
            endif
         enddo
      elseif ( swaptype .eq. SWAP_CB_INTRA ) then
c        --- one box swap moves
         nboxmove = 1
c        --- one box move, but we don't want to waste time trying to 
c        --- perform this move in a box that does not contain the
c        --- correct molecule type.  We will wait and choose the 
c        --- molecule types first, and then choose the box randomly 
c        --- from those boxes that contain the desired molecule types
         lgetbox = .true.
      elseif ( swaptype .eq. SWAP_AVB_ONE
     &        .or. swaptype .eq. SWAP_AVB_TWO
     &        .or. swaptype .eq. SWAP_AVB_THREE ) then
c        --- one box swap moves
         nboxmove = 1
c        --- choose a single box at random from the available boxes
         rmol = twh_random()
         do imove = 1,numboxes
            if ( .not. lfound .and. 
     &           rmol .lt. dble(imove)/dble(numboxes) ) then
               lfound = .true.
               boxtest(2) = imove
               boxtest(1) = imove
            endif
         enddo
      elseif (swaptype .eq. SWAP_CB_GC) then
c        --- grand canonical move always between box 1 and 0
         nboxmove = 2
         if (twh_random() .lt. 0.5d0) then
            boxtest(2) = 1
            boxtest(1) = 0
c           --- TMMC Modification:
            call twh_c_element(GLB_SET,1)
         else
            boxtest(2) = 0
            boxtest(1) = 1
c           --- TMMC Modification:
            call twh_c_element(GLB_SET,-1)
         end if
         lfound = .true.
      elseif ( swaptype .eq. SWAP_CHEMPOT_BOX ) then
c        --- box was passed into this routine
         boxtest(2) = inbox
         boxtest(1) = 0
         lfound = .true.
c        --- need to set nboxmove to 2 to be consistent with logic 
c        --- despite the fact only 1 box will be interesting
         nboxmove = 2
      elseif ( swaptype .eq. SWAP_CHEMPOT_ISO ) then
c        --- isolated molecule, set boxes to zero
         boxtest(2) = MAXBOX
         boxtest(1) = 0
         lfound = .true.
         nboxmove = 2
c        --- set the boxvolume in this box to the one in the first box
         call twh_boxvolume(GLB_GET,1,boxvolume)
         call twh_boxvolume(GLB_SET,boxtest(2),boxvolume)
      else
         write(6,*) 'SWAPMOVES: swaptype',swaptype
         write(6,*) 'invalid swaptype in swapmoves'
         lfailure = .true.
         return
      endif

      if ( .not. lfound .and. .not. lgetbox ) then
         write(6,*) 'SWAPMOVES: boxes not found'
         lfailure = .true.
         return
      endif

c     --- initialize cbmc growth variables
      call twh_resetcbmc(2,wrosen,vrostot,vrosvib,vrosang,vrostor
     &     ,vrosext,vrosintra,vrosinter,vroselect,vrosewald)

c     --- initialize some local parameters
      do ibox = 1,nboxmove
c        --- initialize energy changes
         dvctot(ibox) = 0.0d0
         dvinter(ibox) = 0.0d0
         dvelect(ibox) = 0.0d0
         dvintra(ibox) = 0.0d0
         dvext(ibox) = 0.0d0
         dvvib(ibox) = 0.0d0
         dvben(ibox) = 0.0d0
         dvtor(ibox) = 0.0d0
         dvsolv(ibox) = 0.0d0
         dvquantum(ibox) = 0.0d0
      enddo

c     --- compute the volumes of the boxes
      if ( boxtest(1) .ne. 0 ) then
         call twh_boxvolume(GLB_GET,boxtest(1),boxvolume)
         voltest(1) = boxvolume
      else
         voltest(1) = 0.0d0
      endif
      if ( boxtest(2) .ne. 0 ) then
         call twh_boxvolume(GLB_GET,boxtest(2),boxvolume)
         voltest(2) = boxvolume
      else
         voltest(2) = 0.0d0
      endif
      
      if ( swaptype .eq. SWAP_CHEMPOT_BOX
     &     .or. swaptype .eq. SWAP_CHEMPOT_ISO ) then
c        --- use the inmoltype
         imolty = inmoltype
      else
c        --- randomly determine the molecule type according to the 
c        --- fixed  probabilities
         rmol = twh_random()
         ldone = .false.
         do imol = 1,nmolty
            if ( swaptype .eq. SWAP_CB_INTER ) then
c              --- interbox cb swap uses pm2cbswmt
               call twh_pm2cbswmt(GLB_GET,imol,pmtest)
            elseif ( swaptype .eq. SWAP_CB_INTRA ) then
c              --- intrabox cb swap uses pm2cbswmt
               call twh_pm1cbswmt(GLB_GET,imol,pmtest)
            elseif ( swaptype .eq. SWAP_AVB_ONE ) then
c              --- aggregation volume bias move 1 uses pmavb1mt
               call twh_pmavb1mt(GLB_GET,imol,pmtest)
c              --- store the avb1rad as testrad
               call twh_avb1rad(GLB_GET,testrad)
            elseif ( swaptype .eq. SWAP_AVB_TWO ) then
c              --- avb2 uses pmavb2mt
               call twh_pmavb2mt(GLB_GET,imol,pmtest)
c              --- store the avb2rad as testrad
               call twh_avb2rad(GLB_GET,testrad)
            elseif ( swaptype .eq. SWAP_AVB_THREE ) then
c              --- avb3 uses pmavb3mt
               call twh_pmavb3mt(GLB_GET,imol,pmtest)
c              --- store the avb3rad as testrad
               call twh_avb3rad(GLB_GET,testrad)
            elseif ( swaptype .eq. SWAP_CB_GC) then
c              --- uvt ensemble cbmc move uses pmuvtcbmt
               call twh_pmuvtcbmt(GLB_GET,imol,pmtest)
            elseif ( swaptype .eq. SWAP_RB_INTER ) then
c              --- interbox simple swap uses pm2rbswmt
               call twh_pm2rbswmt(GLB_GET,imol,pmtest)
            else
               write(6,*) 'SWAPMOVES: invalid swaptype for pmtest'
               write(6,*) swaptype
               lfailure = .true.
               return
            endif
            if ( rmol .lt. pmtest ) then
               if ( .not. ldone ) then
                  imolty = imol
                  ldone = .true.
               endif
            endif
         enddo
         if ( .not. ldone ) then
            write(6,*) 'SWAPMOVES: problem selecting primary molecule'
            write(6,*) 'swaptype ',swaptype
            lfailure = .true.
            return
         endif
      endif

c     --- select the secondary molecule type for the AVB MC moves
      if ( swaptype .eq. SWAP_AVB_ONE
     &     .or. swaptype .eq. SWAP_AVB_TWO 
     &     .or. swaptype .eq. SWAP_AVB_THREE) then
         rmol = twh_random()
         ldone = .false.
         do imol = 1,nmolty
            if ( swaptype .eq. SWAP_AVB_ONE ) then
c              --- choose jmolty according to pmavb1ct
               call twh_pmavb1ct(GLB_GET,imolty,imol,pmtest)
            elseif ( swaptype .eq. SWAP_AVB_TWO ) then
c              --- choose jmolty according to pmavb2ct
               call twh_pmavb2ct(GLB_GET,imolty,imol,pmtest)
            elseif ( swaptype .eq. SWAP_AVB_THREE ) then
c              --- choose jmolty according to pmavb3ct
               call twh_pmavb3ct(GLB_GET,imolty,imol,pmtest)
            endif
            if ( rmol .lt. pmtest ) then
               if ( .not. ldone ) then
                  jmolty = imol
                  ldone = .true.
               endif
            endif
         enddo
         if ( .not. ldone ) then
            write(6,*) 'SWAPMOVES: problem selecting secondary molecule'
            write(6,*) 'swaptype ',swaptype
            lfailure = .true.
            return
         endif
      endif

c     --- store number of units in totunits 
      call twh_nunit(GLB_GET,imolty,totunits)

      if ( lgetbox ) then
c        --- for this type of move we have not yet selected a box
c        --- determine valid boxes
         ngood = 0
         do ibox = 1,numboxes
            call twh_ncmt(GLB_GET,ibox,imolty,ncmt)
            if ( ncmt .ne. 0 ) then
               ngood = ngood + 1
               goodbox(ngood) = ibox
            endif
         enddo
         if ( ngood .eq. 0 ) then
c           --- strangely there are no valid boxes, just select one 
c           --- at random and a failure will trigger down below
            boxtest(1) = idint(twh_random()*dble(numboxes))+1
         else
c           --- choose an eligible box
            ranbox = twh_random()*dble(ngood)
            boxtest(1) = goodbox(idint(ranbox)+1)
         endif
c        --- single box move, set boxtest(2) to boxtest(1)
         boxtest(2) = boxtest(1)
      endif

c     --- initialize some values
      lfound = .true.
      regstart = -1
      regend = 0
      if ( swaptype .eq. SWAP_CHEMPOT_BOX
     &     .or. swaptype .eq. SWAP_CHEMPOT_ISO ) then
c        --- not bothering with boxtest(1), just set lempty
         lempty = .true.
      else
c        --- select a chain in boxtest(1) at random 
         call twh_ncmt(GLB_GET,boxtest(1),imolty,ncmt)
         if ( ncmt .eq. 0 ) then
            lempty = .true.
c           --- warn users if this is uvt and the empty box is the 
c           --- ideal gas
            if ( boxtest(1) .eq. 0 ) then
            IF (.NOT.tmmc_flag) THEN
               write(6,*) 'WARNING: no molecule remain in the ideal gas'
     &              ,' for molecule type ',imolty
               write(6,*) 'uvt is not valid if it runs out of molecules'
               write(6,*) 'simply add more of this type of molecule and'
     &              ,' continue the simulation'
            END IF
            endif
         else
            call twh_findtarget(lfailure,swaptype,boxtest(1),imolty
     &           ,chaintest(CRD_OLD),iutry,jmolty,jrem,junit,regstart
     &           ,regend,pointp,num_in,num_out,lfound,xtarget,ytarget
     &           ,ztarget)
            if ( lfailure ) return
         endif
      endif
      if ( lempty ) then
c        --- set iutry at random
         iutry = idint( twh_random()*totunits ) + 1
c        --- set the chain number to a warning flag of 0
         chaintest(CRD_OLD) = 0
      endif

      if ( swaptype .eq. SWAP_CB_INTER
     &     .or. swaptype .eq. SWAP_CB_INTRA
     &     .or. swaptype .eq. SWAP_CB_GC ) then
c        --- increment the number of attempted cb swap moves
         call twh_bncbswap(GLB_INCR,imolty,boxtest(2),boxtest(1),1)
      elseif ( swaptype .eq. SWAP_AVB_ONE ) then
c        --- increment the number of attempted avb1 moves
         call twh_bnavb1(GLB_INCR,imolty,jmolty,boxtest(2),-1,1)
c        --- find the target molecule number and unit
c        --- determine the type of the avb1 move
         avbtype = 2*regstart + regend
c        --- increment this type of avb1 move
         call twh_bnavb1(GLB_INCR,imolty,jmolty,boxtest(2),avbtype,1)
      elseif ( swaptype .eq. SWAP_AVB_TWO ) then
c        --- increment the number of attempted avb2 moves
         call twh_bnavb2(GLB_INCR,imolty,jmolty,boxtest(2),-1,1)
c        --- determine the type of the avb2 move
         avbtype = 2*regstart + regend
c        --- increment this type of avb2 move
         call twh_bnavb2(GLB_INCR,imolty,jmolty,boxtest(2),avbtype,1)
      elseif ( swaptype .eq. SWAP_AVB_THREE ) then
c        --- increment the number of attempted avb3 moves
         call twh_bnavb3(GLB_INCR,imolty,jmolty,boxtest(2),1)
      elseif ( swaptype .eq. SWAP_RB_INTER ) then
c        --- increment the number of attempted rb moves
         call twh_bnrbswap(GLB_INCR,imolty,boxtest(2),boxtest(1),1)
      endif
      
      if ( .not. lfound ) return
      if ( lempty .and. swaptype .ne. SWAP_CB_INTER
     &     .and. swaptype .ne. SWAP_CHEMPOT_BOX
     &     .and. swaptype .ne. SWAP_CHEMPOT_ISO ) return

c     --- assign a molecule number for chaintest(CRD_NEW)
      call twh_nchain(GLB_GET,nchain)
      if ( swaptype .eq. SWAP_RB_INTER ) then
         chaintest(CRD_NEW) = nchain + 1
      elseif ( swaptype .eq. SWAP_CB_INTER
     &        .or. swaptype .eq. SWAP_CB_GC 
     &        .or. swaptype .eq. SWAP_CHEMPOT_BOX ) then
         chaintest(CRD_NEW) = nchain + 1
c        --- increment the count of chemical potential computations
c        --- only if we are inserting into a non-ideal box
c        --- this is a bit of trickery as if we automatically reject 
c        --- then we want to count this as an attempt with 0 added 
c        --- to the average
         if ( boxtest(2) .ne. 0 ) then
            arg = 0.0d0
            keycode = AC_NVT_INSERTION_CHEMPOT
            call twh_averages(lfailure,AVG_INCREMENT,keycode,boxtest(2)
     &           ,imolty,arg,1)
            if ( lfailure ) return
            keycode = AC_NPT_INSERTION_CHEMPOT
            call twh_averages(lfailure,AVG_INCREMENT,keycode,boxtest(2)
     &           ,imolty,arg,1)
            if ( lfailure ) return
            keycode = AC_GIBBS_TOT_CHEMPOT
            call twh_averages(lfailure,AVG_INCREMENT,keycode,boxtest(2)
     &           ,imolty,arg,1)
            if ( lfailure ) return
         endif
      elseif ( swaptype .eq. SWAP_CHEMPOT_ISO ) then
c        --- computing isolated molecule chemical potential
         chaintest(CRD_NEW) = nchain + 1
         arg = 0.0d0
         keycode = AC_ISOLATION_CHEMPOT
         call twh_averages(lfailure,AVG_INCREMENT,keycode,boxtest(2)
     &        ,imolty,arg,1)
         if ( lfailure ) return
      elseif ( swaptype .eq. SWAP_CB_INTRA
     &        .or. swaptype .eq. SWAP_AVB_ONE 
     &        .or. swaptype .eq. SWAP_AVB_TWO
     &        .or. swaptype .eq. SWAP_AVB_THREE ) then
c        --- single box move, same chain for new and old
         chaintest(CRD_NEW) = chaintest(CRD_OLD)
      endif

c     --- assign the old coordinates
      if ( .not. lempty ) then
         do iunit = 1,totunits
            call twh_coordinates(GLB_GET,CRD_REAL,chaintest(CRD_OLD)
     &           ,iunit,xtemp,ytemp,ztemp)
            call twh_coordinates(GLB_SET,CRD_TEMP,chaintest(CRD_OLD)
     &           ,iunit,xtemp,ytemp,ztemp)
         enddo
      endif

#if DEBUG_SWAPMOVES
      write(6,*) 'chaintest:',chaintest
#endif

c     --- set some information for the new chain
      call twh_moltyp(GLB_SET,chaintest(CRD_NEW),imolty)
      call twh_nboxi(GLB_SET,chaintest(CRD_NEW),boxtest(2))

c     --- set some common values for new and old
      call twh_nch_nb_one(GLB_GET,imolty,nch_nb_one)

      if ( lempty ) then
c        --- only need to do new if lempty
         cordstop = CRD_NEW
      else
c        --- do old and new
         cordstop = CRD_OLD
      endif

c     --- always use the real coordinates for all other molecules
      globalcord = CRD_REAL

      do icordflag = CRD_NEW,cordstop,-1
c        --- initialize lexist
         call twh_lexist(GLB_INIT,totunits,.false.)
         
c        --- determine testcord and lnew
         if ( icordflag .eq. CRD_NEW ) then 
c           --- use the new coordinates to compute the energy
            lnew = .true.
            testcord = CRD_TEMP
         else
c           --- Compute weights for the molecule to be removed
            testcord = CRD_REAL
c           --- select a position of the first/starting unit at random
c           --- and calculate the boltzmann weight for the chain to be 
c           --- removed
c           --- set the flags for the old move
            if ( boxtest(1) .eq. 0 ) then
c              --- special case for uvt where we are removing a molecule
c              --- from the ideal gas resevoir.  There are no coordinate
c              --- for this molecules so treat as a new growth but still
c              --- use the typical "old" coordinate storage
               lnew = .true.
               testcord = CRD_SINGLE
            else
c              --- everything else, this is the old growth retracing the
c              --- old coordinates
               lnew = .false.
            endif
         endif

         call twh_initial_atom_insertion(lfailure,lnew,ltermtest
     &        ,nch_nb_one,swaptype,imolty,totunits,iutry,globalcord
     &        ,icordflag,testcord,potentialstyle,regend
     &        ,chaintest,boxtest
     &        ,testrad
     &        ,xtarget,ytarget,ztarget
     &        ,vmintest,w1test,v1tottest,v1exttest
     &        ,v1inttest,v1elctest,v1ewdtest)
         if ( lfailure .or. ltermtest(icordflag) ) return

         if ( icordflag .eq. 2 ) then
c           --- set up the growth pattern for the rest of the molecule
c           --- only needs to be done once as both growths must follow 
c           --- the identical schedule
            call twh_schedule(lfailure,totunits,imolty,isteps,iutry
     &           ,idum,0,2)
            if ( lfailure ) return
         endif
#if DEBUG_SWAPMOVES
      write(6,*) 'swapmoves finished schedule imolty,iutry',imolty,iutry
#endif

c        --- perform the molecule growth for the inserted molecule
         if ( swaptype .ne. SWAP_RB_INTER ) then
c           --- grow the rest of the molecule using configurational-bias
            call twh_rosenbluth(lfailure,lnew,ltermtest(icordflag)
     &           ,testcord,globalcord,chaintest(icordflag),imolty,isteps
     &           ,boxtest(icordflag),totunits,dpscaletest(icordflag)
     &           ,wrosen(icordflag),vrostot(icordflag)
     &           ,vrosvib(icordflag),vrosang(icordflag)
     &           ,vrostor(icordflag),vrosext(icordflag)
     &           ,vrosintra(icordflag),vrosinter(icordflag)
     &           ,vroselect(icordflag),vrosewald(icordflag))
            if ( lfailure ) return
         endif

#if DEBUG_SWAPMOVES
      write(6,*) 'swapmoves back from rosenbluth lnew, wrosen:',lnew
     &        ,wrosen(icordflag)
      write(6,*) 'vrostot:',vrostot(icordflag)
      write(6,*) 'vrosvib:',vrosvib(icordflag)
      write(6,*) 'vrosang:',vrosang(icordflag)
      write(6,*) 'vrostor:',vrostor(icordflag)
      write(6,*) 'vrosext:',vrosext(icordflag)
      write(6,*) 'vrosintra:',vrosintra(icordflag)
      write(6,*) 'vrosinter:',vrosinter(icordflag)
      write(6,*) 'vroselect:',vroselect(icordflag)
      write(6,*) 'vrosewald:',vrosewald(icordflag)
#endif

c        --- termination of cbmc attempt due to walk termination
         if ( ltermtest(icordflag) ) then
            if ( icordflag .eq. 2 ) then
               return
            else
c              --- move rejected, kick out and go down below to the
c              --- warning messages
               goto 666
            endif
         endif
      enddo

      if ( lpartial .or. lintra .or. lbonded ) then
c        --- Begin DC-CBMC Corrections
c        --- calculate the true site-site energy
c        --- compute intermolecular energy for interactions between 
c        --- rcutin and rcut
c        --- get the correction energies

         do ibox = 1,nboxmove
            activebox = boxtest(ibox)
c           --- most moves use one chain as the default
            nchainmove = 1
c           --- set the chain number
            chaininfo(1,3) = chaintest(ibox)
            lcompute = .true.
            if ( swaptype .eq. SWAP_CB_INTER 
     &           .or. swaptype .eq. SWAP_CB_GC
     &           .or. swaptype .eq. SWAP_RB_INTER ) then
c              --- two box move
               if ( ibox .eq. 1 ) then
c                 --- this is the deletion box
                  if ( lempty ) then
c                    --- nothing to delete, we are inserting just to
c                    --- measure the chemical potential
                     lcompute = .false.
                  else
c                    --- the old slot includes the molecule and the new
c                    --- does not
                     chaininfo(1,CRD_OLD) = CRD_REAL
                     chaininfo(1,CRD_NEW) = CRD_NONE
                  endif
               else
c                 --- this is the insertion box so the old slot has
c                 --- nothing and the new slot has the molecule
                  chaininfo(1,CRD_OLD) = CRD_NONE
                  chaininfo(1,CRD_NEW) = CRD_TEMP
               endif
            elseif ( swaptype .eq. SWAP_CB_INTRA
     &              .or. swaptype .eq. SWAP_AVB_ONE
     &              .or. swaptype .eq. SWAP_AVB_TWO
     &              .or. swaptype .eq. SWAP_AVB_THREE ) then
c              --- this is a variation of the single box reinsertion
               chaininfo(1,CRD_OLD) = CRD_REAL
               chaininfo(1,CRD_NEW) = CRD_TEMP
            elseif ( swaptype .eq. SWAP_CHEMPOT_BOX
     &              .or. swaptype .eq. SWAP_CHEMPOT_ISO ) then
c              --- just measuring the chemical potential, there is
c              --- nothing to do in the old box
               if ( ibox .eq. 1 ) then
                  nchainmove = 0
               else
                  chaininfo(1,CRD_NEW) = CRD_TEMP
                  chaininfo(1,CRD_OLD) = CRD_NONE
               endif
            endif
            if ( lcompute ) then
               call twh_energy_change(lfailure,activebox,nchainmove
     &              ,chaininfo,ltermtest,lpartial,lgetcom,lintra,lbonded
     &              ,dvctot(ibox),dvintra(ibox),dvinter(ibox)
     &              ,dvext(ibox),dvelect(ibox),dvsolv(ibox),dvtor(ibox)
     &              ,dvben(ibox),dvvib(ibox),dvquantum(ibox)
     &              ,dvf,ddvdl_lj,ddvdl_c)
               if ( lfailure ) return
#if DEBUG_SWAPMOVES
               write(6,*) 'energy change corrections'
               write(6,*) 'ibox,activebox:',ibox,activebox
               write(6,*) 'dvctot:',dvctot(ibox)
               write(6,*) 'dvintra:',dvintra(ibox)
               write(6,*) 'dvinter:',dvinter(ibox)
               write(6,*) 'dvext:',dvext(ibox)
               write(6,*) 'dvelect:',dvelect(ibox)
               write(6,*) 'dvsolv:',dvsolv(ibox)
               write(6,*) 'dvtor:',dvtor(ibox)
               write(6,*) 'dvben:',dvben(ibox)
               write(6,*) 'dvvb:',dvvib(ibox)
#endif
c              --- energy lambda derivatives not implemented here; not 
c              --- sure it even makes sense.  If get here, print warning
c              --- to user that results likely incorrect.  Same goes for
c              --- foreign energies
               call twh_scp_eval_derivative(GLB_GET,scp_eval_derivative)
               if (scp_eval_derivative) then
                  write(6,*) "WARNING: dU/d[lambda] calculations not ",
     &                 "implemented for swapmoves."
                  lfailure = .true.
                  return
               endif
               call twh_scp_num_foreign_lambda(GLB_GET
     &              ,scp_num_foreign_lambda)
               if (scp_num_foreign_lambda .ne. 0) then
                  write(6,*) "WARNING: foreign energy calculations not "
     &                 ,"implemented for swapmoves."
                  lfailure = .true.
                  return
               endif
c              --- safety check
               do icordflag = 2,cordstop,-1
                  if ( ltermtest(icordflag) ) then
                     write(6,*) 'chaintest(CRD_NEW)',chaintest(CRD_NEW)
     &                    ,'chaintest(CRD_OLD)',chaintest(CRD_OLD)
                     write(6,*) 'SWAPMOVES: strange screwup in DC-CBMC'
                     lfailure = .true.
                     return
                  endif
               enddo
            endif
         enddo

         do ibox = 1,nboxmove
c           --- need to be tricky here to get the chemical potential
c           --- correct.  only want the insertion energy change in the 
c           --- "new" portion
            if ( potentialstyle .eq. POT_INTERNAL ) then
               vtest = dvctot(ibox)
            elseif ( potentialstyle .eq. POT_EXTERNAL ) then
               vtest = dvquantum(ibox)
            endif
c           --- For the old box the energy change sign needs to be 
c           --- reversed
            if ( ibox .eq. 1 ) vtest = -vtest
c           --- safely convert this into a rosenbluth weight
            call twh_getweight(lfailure,0,vtest,wmin,dpmin
     &           ,ltermtest(ibox))
            if ( lfailure ) return

            waddtest(ibox) = waddtest(ibox)*wmin
            dpscaletest(ibox) = dpscaletest(ibox) + dpmin
         enddo
c        --- End DC-CBMC Corrections
      endif

      do icordflag = 2,cordstop,-1
         if ( ltailc .and. boxtest(icordflag) .ne. 0 ) then
c           --- Begin Tail corrections
            if ( swaptype .eq. SWAP_RB_INTER
     &           .or. swaptype .eq. SWAP_CB_INTER
     &           .or. swaptype .eq. SWAP_CB_GC
     &           .or. swaptype .eq. SWAP_CHEMPOT_BOX ) then
               if ( icordflag .eq. 2 ) then
c                 --- we have added a molecule to this box
                  mcount = 1
               else
c                 --- we have subtracted a molecule from this box
                  mcount = -1
               endif
               do imt = 1, nmolty
                  do jmt = 1, nmolty
                     call twh_ncmt(GLB_GET,boxtest(icordflag),jmt,ncmt)
                     if ( jmt .eq. imolty ) then
                        rho = dble(ncmt + mcount) / voltest(icordflag)
                     else
                        rho = dble(ncmt) / voltest(icordflag)
                     endif
                     call twh_ncmt(GLB_GET,boxtest(icordflag),imt,ncmt)
                     call twh_coru(GLB_GET,imt,jmt,coru)
                     if ( imt .eq. imolty ) then
                        vtailtest(icordflag) = vtailtest(icordflag) + 
     &                       dble(ncmt + mcount)*rho*coru
                     else
                        vtailtest(icordflag) = vtailtest(icordflag) + 
     &                       dble(ncmt)*rho*coru
                     endif
                  enddo
               enddo
c              --- compute change in tail energy
               call twh_boxvtail(GLB_GET,boxtest(icordflag),vtest)
               if ( icordflag .eq. 2 ) then
                  vtailtest(icordflag) = vtailtest(icordflag) - vtest
               else
                  vtailtest(icordflag) = vtest - vtailtest(icordflag)
               endif
               call twh_getweight(lfailure,0,vtailtest(icordflag)
     &              ,wmin,dpmin,ltermtest(icordflag))
               if ( lfailure ) return
               waddtest(icordflag) = waddtest(icordflag)*wmin
               dpscaletest(icordflag) = dpscaletest(icordflag) + dpmin
               vrostot(icordflag) = vrostot(icordflag) 
     &              + vtailtest(icordflag)
               vrosinter(icordflag) = vrosinter(icordflag) 
     &              + vtailtest(icordflag)
            endif
c           --- End Tail corrections
         endif

c         --- add contributions of the first bead(s)
         vrostot(icordflag) = vrostot(icordflag) + v1tottest(icordflag)
         if ( potentialstyle .eq. POT_INTERNAL ) then
            vrosinter(icordflag) = vrosinter(icordflag) 
     &           + v1inttest(icordflag)
            vrosext(icordflag) = vrosext(icordflag) 
     &           + v1exttest(icordflag)
            vroselect(icordflag) = vroselect(icordflag) 
     &           + v1elctest(icordflag)
            vrosewald(icordflag) = vrosewald(icordflag) 
     &           + v1ewdtest(icordflag)
         endif
c        --- compute wrosen
         wrosen(icordflag) = w1test(icordflag)*waddtest(icordflag)
     &        *wrosen(icordflag)
c        --- account for vmintest
         call twh_getweight(lfailure,0,vmintest(icordflag),wmin,dpmin
     &        ,ltermtest(icordflag))
         if ( lfailure ) return
         wrosen(icordflag) = wrosen(icordflag)*wmin
         dpscaletest(icordflag) = dpscaletest(icordflag) + dpmin
c        --- rescale wrosen
         call twh_getweight(lfailure,1,vdum,wrosen(icordflag),dpmin
     &        ,ltermtest(icordflag))
         if ( lfailure ) return
         dpscaletest(icordflag) = dpscaletest(icordflag) + dpmin

         if ( icordflag .eq. 2 .and. ltermtest(icordflag) ) then
c           --- new weight is too low - reject move immediately
            return
         endif
      enddo

      if ( ( swaptype .eq. SWAP_CB_INTER
     &     .or. swaptype .eq. SWAP_CB_GC
     &     .or. swaptype .eq. SWAP_CHEMPOT_BOX
     &     .or. swaptype .eq. SWAP_CHEMPOT_ISO )
     &     .and. boxtest(2) .ne. 0 ) then
c        --- store the insertion weight for chemical potential
c        --- for insertion into non-ideal boxes
         arg = wrosen(2)
c        --- account for dpscaletest(2)
         if ( dpscaletest(2) .ge. -6.0d0 ) then
            call twh_scalelolog(GLB_GET,scalelolog)
            arg = arg * scalelolog**(dpscaletest(2))
         else
            write(6,*) 'SWAPMOVES: chemical potential has extremely'
     &           ,' large'
            write(6,*) 'argument in swapmoves - very strange'
     &           ,wrosen(2),dpscaletest(2)
         endif
c        --- store rosenbluth weight in the nvt chem potential
c        --- but don't add one to the count as we already did this up
c        --- above
         if ( swaptype .eq. SWAP_CHEMPOT_ISO ) then
c           --- we were computing the isolated molecule chempot
            keycode = AC_ISOLATION_CHEMPOT
            call twh_averages(lfailure,AVG_INCREMENT,keycode,boxtest(2)
     &           ,imolty,arg,0)
            if ( lfailure ) return
         else
            keycode = AC_NVT_INSERTION_CHEMPOT
            call twh_averages(lfailure,AVG_INCREMENT,keycode,boxtest(2)
     &           ,imolty,arg,0)
            if ( lfailure ) return
c           --- npt insertion chemical potential also includes volume
            keycode = AC_NPT_INSERTION_CHEMPOT
            arg = arg * voltest(2)
            call twh_averages(lfailure,AVG_INCREMENT,keycode,boxtest(2)
     &           ,imolty,arg,0)
            if ( lfailure ) return
c           --- gibbs total chem potential includes number of molecules
c           --- and the thermal debroglie wavelength
            call twh_debroglie(GLB_GET,imolty,debroglie)
            arg = arg /(debroglie**3)
            call twh_ncmt(GLB_GET,boxtest(2),imolty,ncmt)
            arg = arg / dble(ncmt+1)
            keycode = AC_GIBBS_TOT_CHEMPOT
            call twh_averages(lfailure,AVG_INCREMENT,keycode,boxtest(2)
     &           ,imolty,arg,0)
            if ( lfailure ) return
         endif
      endif
            
      if ( swaptype .eq. SWAP_CB_INTER
     &     .or. swaptype .eq. SWAP_CB_INTRA
     &     .or. swaptype .eq. SWAP_CB_GC
     &     .or. swaptype .eq. SWAP_CHEMPOT_BOX
     &     .or. swaptype .eq. SWAP_CHEMPOT_ISO ) then
c        --- the config-bias growth was successful (no overlaps)
         call twh_bscbswap(GLB_INCR,imolty,boxtest(2),boxtest(1),1)
      endif
c     --- if lempty then we were only performing this move to 
c     --- compute the chemical potential.  
c     --- We can now exit this move
      if ( lempty ) return

c     --- final stages to determine acceptance of the move
      if ( swaptype .eq. SWAP_RB_INTER
     &     .or. swaptype .eq. SWAP_CB_INTER ) then
c        --- add in the N/V terms
         call twh_ncmt(GLB_GET,boxtest(1),imolty,ncmt)
         wrosen(2) = wrosen(2)*voltest(2)*dble(ncmt) 
         call twh_ncmt(GLB_GET,boxtest(2),imolty,ncmt)
         wrosen(1) = wrosen(1)*voltest(1)*dble(ncmt+1)
      elseif ( swaptype .eq. SWAP_AVB_ONE ) then
c        --- add in the aggregation volume bias 1 terms
         call twh_pmavb1in(GLB_GET,pin)
         pout = 1.0d0-pin
         voltest(2) = onepi*(4.0/3.0)*(testrad**3)
         call twh_boxvolume(GLB_GET,boxtest(2),boxvolume)
         voltest(1) = boxvolume - voltest(2)
         if ( regstart .eq. 0 ) then
c           --- old conformation was in the inner region
            wrosen(1) = wrosen(1)*voltest(2)/pin
         else
c           --- old conformation was in the outer region
            wrosen(1) = wrosen(1)*voltest(1)/pout
         endif
         if ( regend .eq. 0 ) then
c           --- new conformation is in the innner region
            wrosen(2) = wrosen(2)*voltest(2)/pin
         else
c           --- new conformation is in the outer region
            wrosen(2) = wrosen(2)*voltest(1)/pout
         endif
      elseif ( swaptype .eq. SWAP_AVB_TWO ) then
c        --- add in the volume bias 2 terms
         vol_in = onepi*(4.0/3.0)*(testrad**3)
         call twh_boxvolume(GLB_GET,boxtest(2),boxvolume)
         vol_out = boxvolume - vol_in
         call twh_pmavb2in(GLB_GET,pin)
         pout = 1.0d0-pin
         if ( regstart .eq. 0 ) then
c           --- new conformation is in the outer region
            wrosen(2) = pin*vol_out*dble(num_in)*wrosen(2)
c           --- old conformation was in the inner region
            wrosen(1) = pout*vol_in*dble(num_out+1)*wrosen(1)
         else
c           --- new conformation is in the innner region
            wrosen(2) = pout*vol_in*dble(num_out)*wrosen(2)
c           --- old conformation was in the outer region
            wrosen(1) = pin*vol_out*dble(num_in+1)*wrosen(1)
         endif
      elseif ( swaptype .eq. SWAP_AVB_THREE ) then
c        --- add in the volume bias 3 terms
c        --- new conformation is in the k region, which is stored 
c        --- in num_outer.  we divide out the number density so 
c        --- end up multiplying by the j number
         wrosen(2) = wrosen(2)*dble(num_in)
c        --- old conformation was in the j region, which is stored 
c        --- in num_in
         wrosen(1) = wrosen(1)*dble(num_out+1)
      elseif ( swaptype .eq. SWAP_CB_GC ) then
         if (boxtest(2) .eq. 1) then
c           --- molecule added to box 1
            call twh_uvtfactor(GLB_GET,imolty,uvtfactor)
            wrosen(2) = wrosen(2) * voltest(2)*uvtfactor
c           --- rosenbluth weight is 1 in ideal gas box
            call twh_ncmt(GLB_GET,boxtest(2),imolty,ncmt)
            wrosen(1) = dble(ncmt+1)
            IF (tmmc_flag) THEN
c              --- Calculate change in tmmc_weight for biasing.
               call twh_macro_variable(GLB_GET,macro_variable)
               call twh_tmmc_weight(GLB_GET,macro_variable+1,tmmc_one)
               call twh_tmmc_weight(GLB_GET,macro_variable,tmmc_zero)
               d_tmmc_weight = tmmc_one - tmmc_zero
            END IF
         else
c           --- rosenbluth weight is 1 in ideal gas box
            call twh_ncmt(GLB_GET,boxtest(1),imolty,ncmt)
            wrosen(2) = dble(ncmt)
c           --- molecule removed from box 1
            call twh_uvtfactor(GLB_GET,imolty,uvtfactor)
            wrosen(1) = wrosen(1) * voltest(1)*uvtfactor
            if (tmmc_flag) then
c              --- Calculate change in tmmc_weight for biasing.
               call twh_macro_variable(GLB_GET,macro_variable)
               call twh_tmmc_weight(GLB_GET,macro_variable-1
     &              ,tmmc_negone)
               call twh_tmmc_weight(GLB_GET,macro_variable,tmmc_zero)
               d_tmmc_weight = tmmc_negone - tmmc_zero
            endif
         endif
      endif
#if DEBUG_SWAPMOVES
      write(6,*) 'wrosen,dpscaletest new',wrosen(2),dpscaletest(2)
      write(6,*) 'wrosen dpscaletest old',wrosen(1),dpscaletest(1)
#endif

 666  if ( ltermtest(1) ) then
c        --- this really should not happen
         write(6,*) 'SWAPMOVES: rosenbluth old walk terminated'
         write(6,*) 'This invalidates the move and results in',
     &        ' automatic rejection'
         write(6,*) 'If this is anything other than equilibration',
     &        ' you have serious problems with the simulation'
c        -- reject the move
         return
      else
c        --- compute log of weight
         wnlog = dlog10(wrosen(2)) - dpscaletest(2)*scalecut
c        --- compute old log and difference in logs
         wolog = dlog10(wrosen(1)) - dpscaletest(1)*scalecut
         wdlog = wnlog - wolog
      endif
      IF ( ( swaptype .eq. SWAP_CB_GC ).and.tmmc_flag) THEN
c        --- TMMC Modification:
c        --- 1.) calc. unbiased acceptance prob.
c        --- 2.) calc. biased acceptance prob.
c        --- 3.) log_e factor to put everything in base 10
         IF (wdlog .gt. 0.0D0 ) THEN
            call twh_p_acc(GLB_SET,1.0D0)
         ELSE
            dvalue = 10.0D0**(wdlog)
            call twh_p_acc(GLB_SET,dvalue)
         END IF
         wdlog = wdlog - ( d_tmmc_weight * log_e )
      END IF
c
c     --- initialize laccept
      laccept = .false.

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
      IF (tmmc_flag) THEN
c        --- Accept/reject move if it attempts to take the
c        --- the system out of the specified macrostate space.
         call twh_n_tmmc_max(GLB_GET,n_tmmc_max)
         call twh_c_element(GLB_GET,c_element)
         call twh_macro_variable(GLB_GET,macro_variable)
         IF ( (macro_variable+c_element) .gt. n_tmmc_max ) 
     &        laccept = .FALSE.
         call twh_n_tmmc_min(GLB_GET,n_tmmc_min)
         IF ( (macro_variable+c_element) .lt. n_tmmc_min ) 
     &        laccept = .FALSE.
      END IF

#if DEBUG_SWAPMOVES
      write(6,*) 'laccept:',laccept
      write(6,*) 'v1inttest:',v1inttest
      write(6,*) 'vrostot:',vrostot
      write(6,*) 'vrosinter:',vrosinter
      write(6,*) 'vtailtest:',vtailtest
      write(6,*) 'vrosintra:',vrosintra
      write(6,*) 'vrosvib:',vrosvib
      write(6,*) 'vrostor:',vrostor
      write(6,*) 'vrosang:',vrosang
      write(6,*) 'vrosext:',vrosext
#endif

      wratio = 10.0**(wdlog)
      ectotal = vrostot(2)-vrostot(1)+dvctot(1)
      call twh_beta(GLB_GET,beta)

      if ( laccept ) then
#if DEBUG_SWAPMOVES
      write(6,*) 'SWAPMOVES MOVE ACCEPTED'
      write(6,*) 'chaintest',chaintest
      write(6,*) 'boxtest',boxtest
#endif
         if ( swaptype .eq. SWAP_CB_INTER
     &     .or. swaptype .eq. SWAP_CB_INTRA
     &     .or. swaptype .eq. SWAP_CB_GC ) then
            call twh_bacbswap(GLB_INCR,imolty,boxtest(2),boxtest(1),1)
         elseif ( swaptype .eq. SWAP_AVB_ONE ) then
            call twh_baavb1(GLB_INCR,imolty,jmolty,boxtest(2),avbtype,1)
         elseif ( swaptype .eq. SWAP_AVB_TWO ) then
            call twh_baavb2(GLB_INCR,imolty,jmolty,boxtest(2),avbtype,1)
         elseif ( swaptype .eq. SWAP_AVB_THREE ) then
            call twh_baavb3(GLB_INCR,imolty,jmolty,boxtest(2),1)
         elseif ( swaptype .eq. SWAP_RB_INTER ) then
            call twh_barbswap(GLB_INCR,imolty,boxtest(2),boxtest(1),1)
         endif

c        --- update energies from the rosenbluth portion
         if ( boxtest(1) .ne. 0 ) then
            index = boxtest(1)
            if ( potentialstyle .eq. POT_INTERNAL ) then
               vtest =  -vrostot(1)
               call twh_boxvclassic(GLB_INCR,index,vtest)
               vtest = - vrosinter(1)
               call twh_boxvinter(GLB_INCR,index,vtest)
               vtest =  - vtailtest(1)
               call twh_boxvtail(GLB_INCR,index,vtest)
               vtest = - vrosintra(1)
               call twh_boxvintra(GLB_INCR,index,vtest)
               vtest = - vrosvib(1)
               call twh_boxvbond(GLB_INCR,index,vtest)
               vtest = - vrostor(1)
               call twh_boxvtorsion(GLB_INCR,index,vtest)
               vtest =  - vrosext(1)
               call twh_boxvfield(GLB_INCR,index,vtest)
               vtest = - vrosang(1)
               call twh_boxvangle(GLB_INCR,index,vtest)
               vtest = - vroselect(1) -vrosewald(1)
               call twh_boxvelect(GLB_INCR,index,vtest)
            elseif ( potentialstyle .eq. POT_EXTERNAL ) then
               vtest = - vrostot(1)
               call twh_boxvquantum(GLB_INCR,index,vtest)
            endif
         endif

         if ( boxtest(2) .ne. 0 ) then
            index = boxtest(2)
            if ( potentialstyle .eq. POT_INTERNAL ) then
               vtest = vrostot(2)
               call twh_boxvclassic(GLB_INCR,index,vtest)
               vtest = vrosinter(2)
               call twh_boxvinter(GLB_INCR,index,vtest)
               vtest = vtailtest(2)
               call twh_boxvtail(GLB_INCR,index,vtest)
               vtest = vrosintra(2)
               call twh_boxvintra(GLB_INCR,index,vtest)
               vtest = vrosvib(2)
               call twh_boxvbond(GLB_INCR,index,vtest)
               vtest = vrostor(2)
               call twh_boxvtorsion(GLB_INCR,index,vtest)
               vtest = vrosext(2)
               call twh_boxvfield(GLB_INCR,index,vtest)
               vtest = vrosang(2)
               call twh_boxvangle(GLB_INCR,index,vtest)
               vtest = vroselect(2) + vrosewald(2)
               call twh_boxvelect(GLB_INCR,index,vtest)
            elseif ( potentialstyle .eq. POT_EXTERNAL ) then
               vtest = vrostot(2)
               call twh_boxvquantum(GLB_INCR,index,vtest)
            endif
         endif
c        --- update energies from the engmolec correction
         do ibox = 1,nboxmove
            index = boxtest(ibox)
            if ( index .ne. 0 ) then
               if ( potentialstyle .eq. POT_INTERNAL ) then
                  call twh_boxvclassic(GLB_INCR,index,dvctot(ibox))
                  call twh_boxvinter(GLB_INCR,index,dvinter(ibox))
                  call twh_boxvintra(GLB_INCR,index,dvintra(ibox))
                  call twh_boxvelect(GLB_INCR,index,dvelect(ibox))
                  call twh_boxvsolvation(GLB_INCR,index,dvsolv(ibox))
                  call twh_boxvfield(GLB_INCR,ibox,dvext(ibox))
                  call twh_boxvtorsion(GLB_INCR,index,dvtor(ibox))
                  call twh_boxvangle(GLB_INCR,index,dvben(ibox))
                  call twh_boxvbond(GLB_INCR,index,dvvib(ibox))
               elseif ( potentialstyle .eq. POT_EXTERNAL ) then
                  call twh_boxvquantum(GLB_INCR,index,dvquantum(ibox))
               endif
            endif
         enddo

c        --- update book keeping
         if ( swaptype .eq. SWAP_RB_INTER
     &        .or. swaptype .eq. SWAP_CB_INTER
     &        .or. swaptype .eq. SWAP_CB_GC ) then
            call twh_nboxi(GLB_SET,chaintest(CRD_OLD),boxtest(2))
c           --- set the pointer for the inserted chain in the new box 
c           --- to its original value in the deletion box
            call twh_ncmt(GLB_GET,boxtest(2),imolty,ncmt)
            ipoint = ncmt+1
            call twh_chainlist(GLB_SET,ipoint,boxtest(2),imolty
     &           ,chaintest(CRD_OLD))
c           --- switch out the pointer for the deleted molecule with the
c           --- last pointer in the list and then reset that final 
c           --- pointer to the null value of 0
            call twh_ncmt(GLB_GET,boxtest(1),imolty,ncmt)
            call twh_chainlist(GLB_GET,ncmt,boxtest(1),imolty,ipoint)
            call twh_chainlist(GLB_SET,pointp,boxtest(1),imolty,ipoint)
            call twh_chainlist(GLB_SET,ncmt,boxtest(1),imolty,0)
c           --- update molecule in box counters
            call twh_nchbox(GLB_INCR,boxtest(2),1)
            call twh_nchbox(GLB_DECR,boxtest(1),1)
            call twh_ncmt(GLB_INCR,boxtest(2),imolty,1)
            call twh_ncmt(GLB_DECR,boxtest(1),imolty,1)
         endif

c        --- set the new coordinates for the accepted move
c        --- unless we inserted into the ideal gas box
         if ( boxtest(2) .ne. 0 ) then
            do iunit = 1,totunits
               call twh_coordinates(GLB_GET,CRD_TEMP,chaintest(CRD_NEW)
     &              ,iunit,xtemp,ytemp,ztemp)
#if DEBUG_SWAPMOVES
               write(6,*) 'accepted coord iunit:',iunit,xtemp,ytemp
     &              ,ztemp
#endif
               call twh_coordinates(GLB_SET,CRD_REAL,chaintest(CRD_OLD)
     &              ,iunit,xtemp,ytemp,ztemp)
            enddo
         endif
         
         if ( lewald ) then
c           --- update reciprocal-space sum
            if ( swaptype .eq. SWAP_RB_INTER
     &           .or. swaptype .eq. SWAP_CB_INTER
     &           .or. swaptype .eq. SWAP_CB_GC ) then
               call twh_recip(boxtest(2),2)
               call twh_recip(boxtest(1),2)
            else
               call twh_recip(boxtest(2),2)
            endif
         endif

c        --- update center of mass
         if ( boxtest(2) .ne. 0 ) then
            call twh_ctrmas(lfailure,CRD_REAL,boxtest(2)
     &           ,chaintest(CRD_OLD),CTR_SWAP)
            if ( lfailure ) return
         endif

         if ( lembed ) then
c           --- put the insertion density into the removed chain number
            call twh_eam_rho(GLB_GET,CRD_TEMP,chaintest(CRD_NEW)
     &           ,density)
            call twh_eam_rho(GLB_SET,CRD_TEMP,chaintest(CRD_OLD)
     &           ,density)
c           --- update embedding density
            call twh_reseteamrho(nboxmove,boxtest)
         endif
      endif
 
#if DEBUG_SWAPMOVES
      write(6,*) 'end SWAPMOVES'
#endif
      return
      end

      subroutine twh_initial_atom_insertion(lfailure,lnew
     &     ,ltermtest
     &     ,nch_nb_one,swaptype,imolty,totunits,iutry,globalcord
     &     ,icordflag,testcord,potentialstyle,regend
     &     ,chaintest,boxtest
     &     ,testrad
     &     ,xtarget,ytarget,ztarget
     &     ,vmintest,w1test,v1tottest,v1exttest,v1inttest
     &     ,v1elctest,v1ewdtest)
c     ******************************************************************
c     * performs the initial atom insertion for swap moves             *
c     *                                                                *
c     * split from swapmoves.F 12-07-2008 by M.G. Martin               *
c     * last modified 11-16-2010 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EXPON
#define FUNCTION_RANDOM
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure,lnew
      logical ltermtest
      dimension ltermtest(2)
      integer nch_nb_one,swaptype,imolty,totunits,iutry,icordflag
      integer testcord,potentialstyle,globalcord,regend
      integer chaintest,boxtest
      dimension chaintest(2),boxtest(2)
      double precision testrad
      double precision xtarget,ytarget,ztarget
      dimension xtarget(2),ytarget(2),ztarget(2)
      double precision vmintest,w1test,v1tottest,v1exttest
      dimension vmintest(2),w1test(2),v1tottest(2),v1exttest(2)
      double precision v1inttest,v1elctest,v1ewdtest
      dimension v1inttest(2),v1elctest(2),v1ewdtest(2)
c     --- local variables
c     --- logical scalars
      logical lgenerate,lfound
c     --- logical arrays
      logical lnbterm
      dimension lnbterm(nch_nb_one)
c     --- integer scalars
      integer iunit,icbu,nb_one_gencode,engastyle,ntogrow,itemp,ip
      integer isafemax,isafety,iwalk,ibox,idum
c     --- double precision scalars
      double precision xtemp,ytemp,ztemp,rdist,xcmnew,ycmnew,zcmnew
      double precision xcmi,ycmi,zcmi,xdiff,ydiff,zdiff,twopi,onepi
      double precision vtest,beta,bsum,dist,nbbias,rbf
      double precision rxij,ryij,rzij,onethird
c     --- double precision arrays
      double precision bfactor
      dimension bfactor(nch_nb_one)
      double precision vctry
      dimension vctry(nch_nb_one)
      double precision vtrext
      dimension vtrext(nch_nb_one)
      double precision vtrintra
      dimension vtrintra(nch_nb_one)
      double precision vtrinter
      dimension vtrinter(nch_nb_one)
      double precision vtrelect
      dimension vtrelect(nch_nb_one)
      double precision vtrewald
      dimension vtrewald(nch_nb_one)
      double precision vqtry
      dimension vqtry(nch_nb_one)
      double precision rxp,ryp,rzp
      dimension rxp(NUMAX,nch_nb_one)
      dimension ryp(NUMAX,nch_nb_one)
      dimension rzp(NUMAX,nch_nb_one)

c     --- retrieve constants
      call twh_constant_onethird(GLB_GET,onethird)
      call twh_constant_pi(GLB_GET,onepi)
      call twh_constant_twopi(GLB_GET,twopi)
      vtest = 0.0d0

c     --- start of initial atoms insertion
      if ( .not. lnew ) then
c        --- set the first trial coordinates to the old coords
         if ( swaptype .eq. SWAP_RB_INTER ) then
c           --- inserting multiple atoms
            do iunit = 1,totunits
               call twh_coordinates(GLB_GET,CRD_REAL
     &              ,chaintest(CRD_OLD),iunit,xtemp,ytemp,ztemp)
               rxp(iunit,1) = xtemp
               ryp(iunit,1) = ytemp
               rzp(iunit,1) = ztemp
            enddo
         else
c           --- inserting just a single atom
            call twh_coordinates(GLB_GET,CRD_REAL
     &           ,chaintest(CRD_OLD),iutry,xtemp,ytemp,ztemp)
            rxp(1,1) = xtemp
            ryp(1,1) = ytemp
            rzp(1,1) = ztemp
         endif
      endif

c     --- assign positions for all trial choices for all of the 
c     --- atoms inserted in this first step
      do icbu = 1,nch_nb_one
c        --- set lnbterm logical to false for engatom call
         lnbterm(icbu) = .false.
c        --- set a logical that is true if this is the first bead 
c        --- of the old configuration
         if ( .not. lnew .and. icbu .eq. 1 ) then
c            --- don't generate new coordinates
            lgenerate = .false.
         else
c           --- generate new coordinates
            lgenerate = .true.
         endif
c        --- get the nb_one_generation based on the molecule and box
         if ( boxtest(icordflag) .eq. 0 ) then
c           --- always uniform in box 0
            nb_one_gencode = DIST_UNIFORM
         else
            call twh_cbmc_nb_one_gencode(GLB_GET
     &           ,boxtest(icordflag),imolty,nb_one_gencode)
         endif
c        --- initialize the bfactor to 1.0
         bfactor(icbu) = 1.0d0
         if ( swaptype .eq. SWAP_RB_INTER ) then
c           --- rotational-bias move.  this algorithm inserts the 
c           --- entire molecule at once, with a random rotational 
c           --- displacement about the x, y, and z axes through the
c           --- center of mass
            engastyle = 3
            ntogrow = totunits
c           --- set glist
            do iunit = 1,totunits
               itemp = iunit
               call twh_glist(GLB_SET,iunit,itemp)
            enddo
c           --- use the COM maximum distance of the old chain (as 
c           --- the new chain has the same conformation)
            call twh_rcmu(GLB_GET,chaintest(CRD_OLD),rdist)
            if ( lgenerate ) then
c              --- select a position for the center of mass uniformly
c              --- from the entire insertion box
               call twh_uniformbox(boxtest(icordflag),xcmnew,ycmnew
     &              ,zcmnew)
c               --- translate the old conformation so that it has this
c               --- new center of mass
               call twh_com(GLB_GET,CRD_REAL,chaintest(CRD_OLD)
     &              ,xcmi,ycmi,zcmi)
               xdiff = xcmnew - xcmi
               ydiff = ycmnew - ycmi
               zdiff = zcmnew - zcmi
               do iunit = 1,totunits
                  if ( icordflag .eq. CRD_NEW ) then
                     call twh_coordinates(GLB_GET,CRD_REAL
     &                    ,chaintest(CRD_OLD),iunit
     &                    ,xtemp,ytemp,ztemp)
                     xtemp = xtemp + xdiff
                     ytemp = ytemp + ydiff
                     ztemp = ztemp + zdiff
                     call twh_coordinates(GLB_SET,CRD_SINGLE
     &                    ,chaintest(CRD_NEW),iunit
     &                    ,xtemp,ytemp,ztemp)
                  else
                     call twh_coordinates(GLB_GET,CRD_REAL
     &                    ,chaintest(CRD_OLD),iunit
     &                    ,xtemp,ytemp,ztemp)
                     xtemp = xtemp + xdiff
                     ytemp = ytemp + ydiff
                     ztemp = ztemp + zdiff
                     call twh_coordinates(GLB_SET,CRD_SINGLE
     &                    ,chaintest(CRD_OLD),iunit
     &                    ,xtemp,ytemp,ztemp)
                  endif
               enddo
c              --- compute random rotational displacements on -Pi,Pi
               xdiff = twopi*twh_random()-onepi
               ydiff = twopi*twh_random()-onepi
               zdiff = twopi*twh_random()-onepi
c              --- rotate the molecule
               call twh_rotationmatrix(CRD_SINGLE,chaintest(CRD_OLD)
     &              ,xcmnew,ycmnew,zcmnew,xdiff,ydiff,zdiff)
            endif
c           --- transfer the coordinates to r*p arrays
            do iunit = 1,totunits
               if ( lgenerate ) then
c                 --- use our newly generated coordinates
                  call twh_coordinates(GLB_GET,CRD_SINGLE
     &                 ,chaintest(icordflag),iunit,xtemp,ytemp,ztemp)
               else
c                 --- use the real coordinates
                  call twh_coordinates(GLB_GET,CRD_REAL
     &                 ,chaintest(icordflag),iunit,xtemp,ytemp,ztemp)
               endif
               rxp(iunit,icbu) = xtemp
               ryp(iunit,icbu) = ytemp
               rzp(iunit,icbu) = ztemp
            enddo
            
         elseif ( swaptype .eq. SWAP_CB_INTER
     &           .or. swaptype .eq. SWAP_CB_INTRA 
     &           .or. swaptype .eq. SWAP_CB_GC
     &           .or. swaptype .eq. SWAP_CHEMPOT_BOX
     &           .or. swaptype .eq. SWAP_CHEMPOT_ISO ) then
c           --- there is only one bead to be grown in this step so 
c           --- select a random position for it and assign to rxp
            ntogrow = 1
            engastyle = 1
            rdist = 0.0d0
c           --- put the atom number into the first position of glist
            call twh_glist(GLB_SET,1,iutry)
            if ( nb_one_gencode .eq. DIST_UNIFORM ) then
c              --- select a position uniformly from the entire box
               if ( lgenerate ) then
                  call twh_uniformbox(boxtest(icordflag),xtemp,ytemp
     &                 ,ztemp)
                  rxp(1,icbu) = xtemp
                  ryp(1,icbu) = ytemp
                  rzp(1,icbu) = ztemp
               endif
c              --- uniform sampling does not require a correction
               bfactor(icbu) = 1.0d0
c              --- set the 
            elseif ( nb_one_gencode .eq. DIST_ENERGY_BIAS ) then
               if ( .not. lgenerate ) then
c                 --- need to send the first atom postions
                  xtemp = rxp(1,icbu)
                  ytemp = ryp(1,icbu)
                  ztemp = rzp(1,icbu)
               endif
               ibox = boxtest(icordflag)
               call twh_energybiasbox(lfailure,lgenerate,ibox
     &              ,xtemp,ytemp,ztemp,nbbias)
               if ( lfailure ) return
               bfactor(icbu) = nbbias
               if ( lgenerate ) then
c                 --- store the new positions
                  rxp(1,icbu) = xtemp
                  ryp(1,icbu) = ytemp
                  rzp(1,icbu) = ztemp
               endif
            else
               write(6,*) 'INITIAL_ATOM_INSERTION:'
     &              ,' unknown cbmc_nb_one_gencode'
               write(6,*) nb_one_gencode
               write(6,*) 'box:',boxtest(icordflag)
               write(6,*) 'molecule type:',imolty
               lfailure = .true.
               return
            endif
         elseif ( swaptype .eq. SWAP_AVB_ONE
     &           .or. swaptype .eq. SWAP_AVB_TWO ) then
c           --- only 1 atom to be grown in the first step
            ntogrow = 1
            engastyle = 1
            rdist = 0.0d0
c           --- put the atom number into the first position of glist
            call twh_glist(GLB_SET,1,iutry)
            if ( lgenerate ) then
c              --- select a position uniformly from the inner/outer
c              --- volume biasing region
               if ( regend .eq. 0 ) then
c                 --- inner region
c                 --- select a distance less than avb1rad according 
c                 --- to an r^2 probability distribution
                  dist = (twh_random()**onethird)
                  if ( swaptype .eq. SWAP_AVB_ONE
     &                 .or. swaptype .eq. SWAP_AVB_TWO ) then
                     dist = dist*testrad
                  endif
c                 --- select a direction on the unit sphere
                  call twh_sphere(rxij,ryij,rzij)
c                 --- compute the trial positions
                  rxij = xtarget(1) + dist*rxij
                  ryij = ytarget(1) + dist*ryij
                  rzij = ztarget(1) + dist*rzij
c                 --- make sure this site is inside the central box
                  call twh_putcartinbox(boxtest(icordflag)
     &                 ,rxij,ryij,rzij)
c                 --- store the trial position
                  rxp(1,icbu) = rxij
                  ryp(1,icbu) = ryij
                  rzp(1,icbu) = rzij
               else
c                 --- outer region
c                 --- select uniformly from the box, but reject any 
c                 --- that are too close to the target
                  lfound = .false.
                  isafety = 0
                  isafemax = 10000
                  do while ( .not. lfound 
     &                 .and. isafety .lt. isafemax)
                     isafety = isafety + 1
                     call twh_uniformbox(boxtest(icordflag)
     &                    ,rxp(1,icbu),ryp(1,icbu),rzp(1,icbu))
c                    --- compute distance to the target
                     call twh_mimage(boxtest(icordflag),rxp(1,icbu)
     &                    ,ryp(1,icbu),rzp(1,icbu),xtarget(1)
     &                    ,ytarget(1),ztarget(1),rxij,ryij,rzij)
                     dist = dsqrt(rxij*rxij + ryij*ryij + rzij*rzij)
                     if ( swaptype .eq. SWAP_AVB_ONE
     &                    .or. swaptype .eq. SWAP_AVB_TWO ) then
                        if ( dist .gt. testrad ) lfound = .true.
                     endif
                  enddo
                  if ( isafety .eq. isafemax ) then
                     write(6,*) 'INITIAL_ATOM_INSERTION:'
     &                    ,' too many tries to get position'
                     write(6,*) 'isafety,isafemax,imolty'
     &                    ,isafety,isafemax,imolty
                     if ( swaptype .eq. SWAP_AVB_ONE ) then
                        write(6,*) 'avb1rad',testrad
                     endif
                     if ( swaptype .eq. SWAP_AVB_TWO ) then
                        write(6,*) 'avb2rad',testrad
                     endif
                     write(6,*) 'move rejected, swaptype ',swaptype
                     return
                  endif
               endif
            endif
         elseif ( swaptype .eq. SWAP_AVB_THREE ) then
c           --- put the atom number into the first position of glist
            call twh_glist(GLB_SET,1,iutry)
            ntogrow = 1
            engastyle = 1
            rdist = 0.0d0
            if ( lgenerate ) then
c              --- select a position uniformly from the inner biasing
c              --- region of kchain (target 2)
c              --- select a distance less than avb3rad according to 
c              --- an r^2 probability distribution
               dist = (twh_random()**onethird)*testrad
c              --- select a direction on the unit sphere
               call twh_sphere(rxij,ryij,rzij)
c              --- compute the trial positions
               rxij = xtarget(icordflag) + dist*rxij
               ryij = ytarget(icordflag) + dist*ryij
               rzij = ztarget(icordflag) + dist*rzij
c              --- make sure site is inside the central image box
               call twh_putcartinbox(boxtest(icordflag),rxij,ryij
     &              ,rzij)
c              --- store the trial position
               rxp(1,icbu) = rxij
               ryp(1,icbu) = ryij
               rzp(1,icbu) = rzij
            endif
         else
            write(6,*) 'INITIAL_ATOM_INSERTION: unknown swaptype'
     &           ,swaptype
            lfailure = .true.
            return
         endif
      enddo

c     ---  insert the first atoms
      call twh_engatom(lfailure,testcord,globalcord,engastyle
     &     ,ltermtest(icordflag),chaintest(icordflag),imolty
     &     ,boxtest(icordflag),nch_nb_one,idum,ntogrow,rdist,lnbterm
     &     ,NUMAX,nch_nb_one,rxp,ryp,rzp,vctry,vtrext
     &     ,vtrintra,vtrinter,vtrelect,vtrewald,vqtry)
      if ( lfailure ) return

      if ( ltermtest(icordflag) ) then
         if ( icordflag .eq. CRD_NEW ) then
c           --- bail out if all sites are hard overlap
            return
         else
c           --- should not happen for old conformation
            write(6,*) 'INITIAL_ATOM_INSERTION:'
     &           ,' disaster-overlap for 1st old bead'
            lfailure = .true.
            return
         endif
      endif

c     --- factor out a minimum energy to keep bfactors in line
      lfound = .false.
      do ip = 1, nch_nb_one
         if ( potentialstyle .eq. POT_INTERNAL ) then
            vtest = vctry(ip)
         elseif ( potentialstyle .eq. POT_EXTERNAL ) then
            vtest = vqtry(ip)
         else
            write(6,*) 'INITIAL_ATOM_INSERTION: unknown potentialstyle'
            write(6,*) potentialstyle
            lfailure = .true.
            return
         endif
         if ( .not. lnbterm(ip) ) then
            if ( .not. lfound ) then
c              --- this is our first energy, set is as minimum
               vmintest(icordflag) = vtest
               lfound = .true.
            else
               if ( vtest .lt. vmintest(icordflag) ) 
     &              vmintest(icordflag) = vtest
            endif
         endif
      enddo

c     --- compute the rosenbluth weight for part 1 of the move
      w1test(icordflag) = 0.0d0
      do ip = 1, nch_nb_one
         if ( potentialstyle .eq. POT_INTERNAL ) then
            vtest = vctry(ip)
         elseif ( potentialstyle .eq. POT_EXTERNAL ) then
            vtest = vqtry(ip)
         endif
         if ( .not. lnbterm(ip) ) then
            call twh_beta(GLB_GET,beta)
            bfactor(ip) = bfactor(ip)*twh_expon(-beta*(vtest
     &           -vmintest(icordflag)))
         else
            bfactor(ip) = 0.0d0
         endif
         w1test(icordflag) = w1test(icordflag) + bfactor(ip)
      enddo

      if ( lnew ) then
c        --- select one choice at random with a bias toward favorable
c        --- energies.  Note this can happen for either icordflag as 
c        --- there is some trickery in the grand canonical ensemble
         if ( nch_nb_one .gt. 1 ) then
            rbf = w1test(icordflag)*twh_random()
            bsum = 0.0d0 
#if DEBUG_SWAPMOVES
      write(6,*) 'w1test:',w1test
#endif
            do ip = 1, nch_nb_one
               if ( .not. lnbterm(ip) ) then
                  bsum = bsum + bfactor(ip)
                  if ( rbf .lt. bsum ) then
c                    --- select this position
                     iwalk = ip
                     goto 180
                  endif
               endif
            enddo
c           --- error check
            write(6,*) 'big time screwup in INITIAL_ATOM_INSERTION'
            write(6,*) 'w1test(2)',w1test(2)
            lfailure = .true.
            return
c           --- pseudo while loop success kickout
 180        continue
         else
            iwalk = 1
         endif
      else
c        --- take the old conformation
         iwalk = 1
      endif

c     --- save the energies for the selected walk
      if ( potentialstyle .eq. POT_INTERNAL ) then
         v1tottest(icordflag) = vctry(iwalk)  
         v1exttest(icordflag) = vtrext(iwalk)
         v1inttest(icordflag) = vtrinter(iwalk)
         v1elctest(icordflag) = vtrelect(iwalk)
         v1ewdtest(icordflag) = vtrewald(iwalk)
      elseif ( potentialstyle .eq. POT_EXTERNAL ) then
         v1tottest(icordflag) = vqtry(iwalk)
      endif
#if DEBUG_SWAPMOVES
      write(6,*) 'after initial atom insertion'
      write(6,*) 'icordflag,iwalk:',icordflag,iwalk
      write(6,*) 'vctry:',vctry(iwalk)
      write(6,*) 'vtrext:',vtrext(iwalk)
      write(6,*) 'vtrinter:',vtrinter(iwalk)
      write(6,*) 'vtrelect:',vtrelect(iwalk)
      write(6,*) 'vtrewald:',vtrewald(iwalk)
#endif

c     --- divide out the number of trials
      w1test(icordflag) = w1test(icordflag) / dble(nch_nb_one)

      if ( lnew ) then
c        --- transfer the accepted coordinates into the coord arrays
         if ( swaptype .eq. SWAP_RB_INTER ) then
c           --- multiple atoms
            do iunit = 1,totunits
               xtemp = rxp(iunit,iwalk)
               ytemp = ryp(iunit,iwalk)
               ztemp = rzp(iunit,iwalk)
               call twh_coordinates(GLB_SET,testcord
     &              ,chaintest(CRD_NEW),iunit,xtemp,ytemp,ztemp)
            enddo
         else
c           --- single atom to assign
            xtemp = rxp(1,iwalk)
            ytemp = ryp(1,iwalk)
            ztemp = rzp(1,iwalk)
            call twh_coordinates(GLB_SET,testcord,chaintest(CRD_NEW)
     &           ,iutry,xtemp,ytemp,ztemp)
         endif
      endif
c     --- end of initial atom insertions
      return
      end
