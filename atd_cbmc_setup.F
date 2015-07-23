#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
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
      subroutine twh_atd_cbmc_setup(lfailure)
c     ******************************************************************
c     * sets up all of the arbitrary trial distribution information for*
c     * use in the configurational-bias routines                       *
c     *                                                                *
c     * split from checkstruc 05-18-2006 by M.G. Martin                *
c     * last modified 07-24-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- debug logical
#define DEBUG_ATD_CBMC_SETUP 0
c     --- function declarations
#if DEBUG_ATD_CBMC_SETUP
#define FUNCTION_GAUSSPROB
#endif
#define FUNCTION_EXPON
#define FUNCTION_GET_ANGLESTYLE
#define FUNCTION_GET_BENCOEFF
#define FUNCTION_GET_IJBEN2
#define FUNCTION_GET_IJBEN3
#define FUNCTION_GET_IJTOR2
#define FUNCTION_GET_IJTOR3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_IJVIB
#define FUNCTION_GET_INBEN
#define FUNCTION_GET_INVIB
#define FUNCTION_GET_ITBEN
#define FUNCTION_GET_INTOR
#define FUNCTION_LAWOFCOSINE
#define FUNCTION_VANGLE
#define FUNCTION_VBOND
#define FUNCTION_VTORSION
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
      character*20 atest
      character*25 distname
      character*30 peakstyle,cbmc_analysis,cbmc_fit_strategy
      character*50 twobbias_style,twobbias_init,threebbias_init
      character*50 threebbias_style
c     --- local logical scalars
      logical lovrlap,lfound,lrigid,linit,lthreemem
c     --- local logical arrays
      logical lyescbmc
      dimension lyescbmc(NTMAX)
c     --- local integer parameters
      integer maxextrema
      parameter (maxextrema = 2*(CBMC_MAXTORPEAKS+CBMC_MAXBENDPEAKS))
c     --- local integer scalars
      integer imolty,iunit,ivib,junit,itest,itype,ibend,kunit,itor
      integer nextrema,nmaxima,iextrema,lunit,ipeak,npeak
      integer imax,junitb,kunitb,icrossbend,jbend,kbend,ibondmax
      integer bij_type,bjk_type,iscan_bij,iscan_bjk
      integer type,nmolty,nunit,ibox,cbmc_bondgen_code,cbmc_bendgen_code
      integer cbmc_dihedgen_code,bonddist,vibtype,ivalue,ncompute
      integer icompute,numboxes,nummol,invib
      integer count,bendtype,bendstyle,vibstyle
      integer angleadist
#if DEBUG_ATD_CBMC_SETUP
      integer style_code,imgm
#endif
c     --- local integer arrays
      integer computelist
      dimension computelist(NTMAX)
c     --- local dp scalars
      double precision rmax,sum_x,sum_xx,length,theta,truenorm,arbnorm
      double precision energy,weight,sum_w,distone,disttwo,angcos
      double precision phione,phitwo,thetaone,thetatwo,lolimit,hilimit
      double precision xub,yub,zub,gamma,phic,xxx,yyy,zzz
      double precision bij_equil,bjk_equil,costheta
      double precision bij_energy,bjk_energy,aijk_energy,bik_distance
      double precision rik,riksq,onepi,dvalue,twopi,beta
      double precision twobvibone,twobvibtwo,bond_sdev_multiplier
      double precision benda_mult,bendb_mult,dihed_mult
      double precision viblo,vibhi,onethird,mean,std_dev
      double precision radtodeg,ptest
      double precision lovalue,hivalue,dpar
#if DEBUG_ATD_CBMC_SETUP
      double precision ideal_frac
      double precision phidisp,pdf_bias,lodim,hidim
#endif
c     --- local dp arrays
      double precision xbacord,ybacord,zbacord
      dimension xbacord(3),ybacord(3),zbacord(3)
      double precision xtcord,ytcord,ztcord
      dimension xtcord(4),ytcord(4),ztcord(4)
      double precision extrema
      dimension extrema(maxextrema,2)
      double precision prevweight,prevloc,prevenergy
      dimension prevweight(2),prevloc(2),prevenergy(2)
#if DEBUG_ATD_CBMC_SETUP
      double precision delta_dist
      dimension delta_dist(0:CBMCMAXTOR)
#endif

c     --- retrieve constants
      call twh_constant_pi(GLB_GET,onepi)
      call twh_constant_twopi(GLB_GET,twopi)
      call twh_constant_onethird(GLB_GET,onethird)
      call twh_nmolty(GLB_GET,nmolty)
      call twh_beta(GLB_GET,beta)
      call twh_constant_radtodeg(GLB_GET,radtodeg)
      call twh_cbmc_analysis(GLB_GET,cbmc_analysis)
      call twh_cbmc_fit_strategy(GLB_GET,cbmc_fit_strategy)

c     --- open the file for atd cbmc information output
      call twh_openfile(lfailure,44,'towhee_cbmc','DEFAULT','UNKNOWN'
     &     ,'DEFAULT')
      if ( lfailure ) return

c     --- determine which molecules need cbmc parameters
      do imolty = 1,nmolty
         lyescbmc(imolty) = .false.
      enddo

c     --- check for a cbmc init
      call twh_linit(GLB_GET,linit)
      if ( linit ) then
         call twh_numboxes(GLB_GET,numboxes)
         do imolty = 1,nmolty
            do ibox = 1,numboxes
               call twh_initmol(GLB_GET,ibox,imolty,nummol)
               if ( nummol .gt. 0 ) then
                  call twh_initstyle(GLB_GET,ibox,imolty,atest)
                  if ( atest .eq. 'full cbmc' 
     &                 .or. atest .eq. 'helix cbmc' ) then
                     lyescbmc(imolty) = .true.
                  endif
               endif
            enddo
         enddo
      endif

      do imolty = 1,nmolty
         call twh_chempotperstep(GLB_GET,imolty,ivalue)
         if ( ivalue .gt. 0 ) lyescbmc(imolty ) = .true.
      enddo
c     --- 2 box cb swap
      call twh_pm2boxcbswap(GLB_GET,dvalue)
      if ( dvalue .gt. SMALLEST ) then
         do imolty = 1,nmolty
            call twh_pm2cbswmt(GLB_GET,imolty,ptest)
            if ( ptest .gt. SMALLEST ) lyescbmc(imolty) = .true.
         enddo
      endif
c     --- grand canonical
      call twh_pmuvtcbswap(GLB_GET,dvalue)
      if ( dvalue .gt. SMALLEST ) then
         do imolty = 1,nmolty
            call twh_pmuvtcbmt(GLB_GET,imolty,ptest)
            if ( ptest .gt. SMALLEST ) lyescbmc(imolty) = .true.
         enddo
      endif
c     --- 1 box cb swap
      call twh_pm1boxcbswap(GLB_GET,dvalue)
      if ( dvalue .gt. SMALLEST ) then
         do imolty = 1,nmolty
            call twh_pm1cbswmt(GLB_GET,imolty,ptest)
            if ( ptest .gt. SMALLEST ) lyescbmc(imolty) = .true.
         enddo
      endif
c     --- cb regrowth
      call twh_pmcb(GLB_GET,dvalue)
      if ( dvalue .gt. SMALLEST ) then
         do imolty = 1,nmolty
            call twh_pmcbmt(GLB_GET,imolty,ptest)
            if ( ptest .gt. SMALLEST ) lyescbmc(imolty) = .true.
         enddo
      endif
c     --- backbone cb regrowth
      call twh_pmback(GLB_GET,dvalue)
      if ( dvalue .gt. SMALLEST ) then
         do imolty = 1,nmolty
            call twh_pmbkmt(GLB_GET,imolty,ptest)
            if ( ptest .gt. SMALLEST ) lyescbmc(imolty) = .true.
         enddo
      endif
c     --- sidechain cb regrowth
      call twh_pmcbside(GLB_GET,dvalue)
      if ( dvalue .gt. SMALLEST ) then
         do imolty = 1,nmolty
            call twh_pmcbsidemt(GLB_GET,imolty,ptest)
            if ( ptest .gt. SMALLEST ) lyescbmc(imolty) = .true.
         enddo
      endif
c     --- avb1 regrowth
      call twh_pmavb1(GLB_GET,dvalue)
      if ( dvalue .gt. SMALLEST ) then
         do imolty = 1,nmolty
            call twh_pmavb1mt(GLB_GET,imolty,ptest)
            if ( ptest .gt. SMALLEST ) lyescbmc(imolty) = .true.
         enddo
      endif
c     --- avb2 regrowth
      call twh_pmavb2(GLB_GET,dvalue)
      if ( dvalue .gt. SMALLEST ) then
         do imolty = 1,nmolty
            call twh_pmavb2mt(GLB_GET,imolty,ptest)
            if ( ptest .gt. SMALLEST ) lyescbmc(imolty) = .true.
         enddo
      endif
c     --- avb3 regrowth
      call twh_pmavb3(GLB_GET,dvalue)
      if ( dvalue .gt. SMALLEST ) then
         do imolty = 1,nmolty
            call twh_pmavb3mt(GLB_GET,imolty,ptest)
            if ( ptest .gt. SMALLEST ) lyescbmc(imolty) = .true.
         enddo
      endif

      ncompute = 0
      do imolty = 1,nmolty
         if ( lyescbmc(imolty) ) then
            ncompute = ncompute + 1
            computelist(ncompute) = imolty
         endif
      enddo

c     --- set ibox to zero so we use the full coulombic energies instead
c     --- of one scaled by Ewald parameters
      ibox = 0
c     --- initialize
      do ipeak = 1,2
         prevweight(ipeak) = 0.0d0
         prevloc(ipeak) = 0.0d0
         prevenergy(ipeak) = 0.0d0
      enddo

      write(6,'(a,1x,a)')
     &     'Determining Arbitrary Trial Distribution Parameters:'
     &     ,'see towhee_cbmc for full details'

c     --- set up information needed to generation trials in a cbmc move
      write(44,'(a53)')
     &     'Determining CBMC Arbitrary Trial Distributions: Bonds'
      call twh_cbmc_bondgen_code(GLB_GET,cbmc_bondgen_code)
      call twh_max_bond_length(GLB_GET,rmax)
      do icompute = 1,ncompute
         imolty = computelist(icompute)
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            call twh_invib(GLB_GET,imolty,iunit,invib)
            do ivib = 1,invib
               junit = twh_get_ijvib(imolty,iunit,ivib)
               call twh_itvib(GLB_GET,imolty,iunit,ivib,vibtype)
               call twh_vibstyle(GLB_GET,vibtype,vibstyle)
c              --- check the style to see if is one of the special cases
               if ( vibstyle .eq. BOND_STYLE_FIXED ) then
c                 --- fixed bond length, everything else infinite
                  bonddist = DIST_DELTA
c                 --- store the one valid bond length
                  call twh_vibcoeff(GLB_GET,vibtype,0,length)
                  call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib,3
     &                 ,length)
c                 --- only one value, normalization is 1.0 for both
                  truenorm = 1.0d0
                  arbnorm = 1.0d0
               elseif ( vibstyle .eq. BOND_STYLE_NONE
     &                 .or. cbmc_bondgen_code .eq. CBMC_IDEAL) then
c                 --- with no energy terms this is and ideal r squared
c                 --- distribution
                  bonddist = DIST_R_SQ
                  lolimit = 0.0d0
                  hilimit = rmax
c                 --- transform for use generating potentials
                  hivalue = hilimit**3 - lolimit**3
                  lovalue = lolimit**3
c                 --- store these values for use generating
                  call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                 ,3,lovalue)
                  call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                 ,4,hivalue)
c                 --- renorm to the full continuous distribution
                  truenorm = 3.0d0/rmax**3
                  arbnorm = 3.0d0/hivalue
               elseif ( vibstyle .eq. BOND_STYLE_SQUARE_WELL ) then
c                 --- infinite square well with defined upper and lower
c                 --- bounds defined by the potential
c                 --- this is the bounded r square distribution
                  bonddist = DIST_R_SQ
                  call twh_vibcoeff(GLB_GET,vibtype,0,length)
                  call twh_vibcoeff(GLB_GET,vibtype,1,lolimit)
                  call twh_vibcoeff(GLB_GET,vibtype,2,hilimit)
c                 --- transform for use generating potentials
                  hivalue = hilimit**3 - lolimit**3
                  lovalue = lolimit**3
c                 --- store these values for use generating
                  call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                 ,3,lovalue)
                  call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                 ,4,hivalue)
c                 --- renorm to the full continuous distribution
                  truenorm = 3.0d0/rmax**3
                  arbnorm = 3.0d0/hivalue
               else
c                 --- not a special case based on the potential
                  if ( cbmc_bondgen_code 
     &                 .eq. CBMC_R_SQ_WITH_BOUNDS ) then
                     bonddist = DIST_R_SQ
c                    --- reset vibrang values for use in geometry
                     call twh_vibcoeff(GLB_GET,vibtype,0,length)
                     call twh_vibrang(GLB_GET,1,viblo)
                     call twh_vibrang(GLB_GET,2,vibhi)
                     lolimit = viblo*length
                     hilimit = vibhi*length
c                    --- transform for use generating potentials
                     hivalue = hilimit**3 - lolimit**3
                     lovalue = lolimit**3
c                    --- store these values for use generating
                     call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                    ,3,lovalue)
                     call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                    ,4,hivalue)
c                    --- even though our generation is bound, we use
c                    --- the ideal distribution to normalize
                     truenorm = 3.0d0/rmax**3
                     arbnorm = 3.0d0/hivalue
                  elseif ( cbmc_bondgen_code 
     &                    .eq. CBMC_GLOBAL_GAUSSIAN ) then
                     bonddist = DIST_GAUSSIAN
c                    --- uses a gaussian with parameters based on the
c                    --- equilibrium bond length and the default
c                    --- standard deviation
c                    --- make sure it is a vibstyle that has a valid
c                    --- equilibrium bond length
                     if ( vibstyle .eq. BOND_STYLE_NONBOND 
     &                    .or. vibstyle .eq. BOND_STYLE_FENE ) then
                        write(6,*) 'Global Gaussian does not make'
     &                       ,' sense with potentials that do not have'
     &                       ,' a clearly defined equilibrium distance'
                        write(6,*) 'imolty,iunit,vibtype,vibstyle:'
     &                       ,imolty,iunit,vibtype,vibstyle
                        lfailure = .true.
                        return
                     endif
                     call twh_vibcoeff(GLB_GET,vibtype,0,mean)
                     call twh_sdevvib(GLB_GET,std_dev)
                     call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                    ,3,mean)
                     call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                    ,4,std_dev)
c                    --- normalization is determined in the gaussian
c                    --- selection process, true is the usual continuous
                     truenorm = 3.0d0/rmax**3
                     arbnorm = 1.0d0
                  elseif ( cbmc_bondgen_code 
     &                    .eq. CBMC_AUTOFIT_GAUSSIAN ) then
                     bonddist = DIST_GAUSSIAN
c                    --- determine gaussian parameters
                     imax = 1000
                     sum_w = 0.0d0
                     sum_x = 0.0d0
                     sum_xx = 0.0d0
                     do itest = 0,imax
                        length = dble(itest)*rmax/dble(imax)
                        lovrlap = .false.
                        energy = twh_vbond(lfailure,vibtype,length
     &                       ,imolty,iunit,junit,lovrlap)
                        if ( lfailure ) return
                        if ( lovrlap ) then
                           weight = 0.0d0
                        else
                           weight = twh_expon(-beta*energy)
     &                          *length*length
                        endif
                        sum_w = sum_w + weight
                        sum_x = sum_x + weight*length
                        sum_xx = sum_xx + weight*length*length
                     enddo
c                    --- average both quantities
                     sum_x = sum_x / sum_w
                     sum_xx = sum_xx / sum_w
c                    --- assign the mean
                     call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                    ,3,sum_x)
c                    --- assign the standard deviation
                     call twh_bond_sdev_multiplier(GLB_GET
     &                    ,bond_sdev_multiplier)
                     dvalue = bond_sdev_multiplier
     &                    *dsqrt(sum_xx - sum_x*sum_x)
                     call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &                    ,4,dvalue)
c                    --- normalization is ideal based only on rmax
                     truenorm = 3.0d0/rmax**3
                     arbnorm = 1.0d0
                  else
                     write(6,*) 'ATD_CBMC_SETUP: unknown cbmc_bondgen'
                     write(6,*) cbmc_bondgen_code
                     lfailure = .true.
                     return
                  endif
               endif
c              --- store the bonddist
               call twh_cbbonddist(GLB_SET,imolty,iunit,ivib,bonddist)
c              --- cheack bonddist and get its name
               if ( bonddist .eq. DIST_GAUSSIAN ) then
                  distname = 'Gaussian'
               elseif ( bonddist .eq. DIST_DELTA ) then
                  distname = 'Delta'
               elseif ( bonddist .eq. DIST_R_SQ ) then
                  distname = 'r squared with bounds'
               else
                  write(6,*) 'ATD_CBMC_SETUP: unknown bonddist',bonddist
                  lfailure = .true.
                  return
               endif
c              --- output trial generation information
               write(44,'(a7,1x,i5,1x,a6,1x,2(i5,1x),1x,a4,1x,a25)')
     &              'Moltyp:',imolty,'Units:',iunit,junit
     &              ,'ATD:',distname
c              --- store the bondnorms
               call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &              ,1,truenorm)
               call twh_cbmcbondinfo(GLB_SET,imolty,iunit,ivib
     &              ,2,arbnorm)
               if ( lfailure ) return
c              --- output trial generation details
               if ( bonddist .eq. DIST_GAUSSIAN ) then
                  call twh_cbmcbondinfo(GLB_GET,imolty,iunit,ivib
     &                 ,3,mean)
                  call twh_cbmcbondinfo(GLB_GET,imolty,iunit,ivib
     &                 ,4,std_dev)
                  write(44,'(6x,a5,1x,f10.5,1x,a8,1x,f10.7)') 
     &                 ,'Mean:',mean,'Std Dev:',std_dev
               elseif ( bonddist .eq. DIST_DELTA ) then
                  call twh_cbmcbondinfo(GLB_GET,imolty,iunit,ivib,3
     &                 ,length)
                  write(44,'(6x,a7,1x,f10.7)') 
     &                 ,'Length:',length
               elseif ( bonddist .eq. DIST_R_SQ ) then
                  call twh_cbmcbondinfo(GLB_GET,imolty,iunit,ivib
     &                 ,3,lovalue)
                  call twh_cbmcbondinfo(GLB_GET,imolty,iunit,ivib
     &                 ,4,hivalue)
                  write(44,'(6x,a10,1x,f10.7,1x,a11,1x,f10.7)') 
     &                 ,'Low bound:',lovalue**(1.0d0/3.0d0)
     &                 ,'High bound:',(hivalue+lovalue)**(1.0d0/3.0d0)
               else
                  write(6,*) 'ATD_CBMC_SETUP: unknown bonddist',bonddist
                  lfailure = .true.
                  return
               endif
c              --- compute the expected acceptance rate
               if ( cbmc_analysis .eq. 'full' ) then
                  call twh_cbbond_par(lfailure,imolty,iunit,ivib,dpar)
                  write(44,'(6x,a25,1x,f7.3,a1)')
     &                 'Expected Acceptance Rate:',100.0d0*dpar,'%'
               endif
            enddo
         enddo
      enddo

      call twh_cbmc_bendgen_code(GLB_GET,cbmc_bendgen_code)
      call twh_bend_a_sdev_multiplier(GLB_GET,benda_mult)
      call twh_bend_b_sdev_multiplier(GLB_GET,bendb_mult)
      write(44,'(a)') 
     &     'Determining CBMC Arbitrary Trial Distributions: Angle A'
      do icompute = 1,ncompute
         imolty = computelist(icompute)
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            do ibend = 1,twh_get_inben(imolty,iunit)
               call twh_ijben2(GLB_GET,imolty,iunit,ibend,junit)
               call twh_ijben3(GLB_GET,imolty,iunit,ibend,kunit)
               call twh_itben(GLB_GET,imolty,iunit,ibend,bendtype)
               bendstyle = twh_get_anglestyle(bendtype)
c              --- first, check for the special case of 3-membered ring
               call twh_invib(GLB_GET,imolty,iunit,invib)
               lthreemem = .false.
               do ivib = 1,invib
                  call twh_ijvib(GLB_GET,imolty,iunit,ivib,itest)
                  if ( kunit .eq. itest ) then
c                    --- this angle is also a bond
                     lthreemem = .true.
                     call twh_itvib(GLB_GET,imolty,iunit,ivib,vibtype)
                     call twh_vibstyle(GLB_GET,vibtype,vibstyle)
                  endif
               enddo
               if ( lthreemem .and. 
     &              vibstyle .eq. BOND_STYLE_FIXED ) then
c                  --- the rigid bond of the 3 membered ring determines
c                  --- the angle
                  angleadist = DIST_TRIANGLE
c                 --- one peak
                  ipeak = 1
                  call twh_cbmcbendcount(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak)
c                 --- store this bond length
                  call twh_vibcoeff(GLB_GET,vibtype,0,length)
                  call twh_cbmcbendinfo(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak,1,length)
               elseif ( bendstyle .eq. ANG_STYLE_RIGID ) then
c                 --- rigid, delta distribution, single valid angle
                  angleadist = DIST_DELTA
c                 --- one peak
                  ipeak = 1
                  call twh_cbmcbendcount(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak)
c                 --- store the one valid angle
                  call twh_bencoeff(GLB_GET,bendtype,0,mean)
                  call twh_cbmcbendinfo(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak,1,mean)
               elseif ( bendstyle .eq. ANG_STYLE_SQUARE_WELL ) then
c                 --- this is the bounded Sin distribution
                  angleadist = DIST_BOUNDED_SINE
c                 --- one peak
                  ipeak = 1
                  call twh_cbmcbendcount(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak)
c                 --- store the distance lo and hi bounds, will combine
c                 --- those with the actual bond lengths to generate
c                 --- the angle
                  call twh_bencoeff(GLB_GET,bendtype,1,lovalue)
                  call twh_cbmcbendinfo(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak,1,lovalue)
                  call twh_bencoeff(GLB_GET,bendtype,2,hivalue)
                  call twh_cbmcbendinfo(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak,2,hivalue)
               elseif ( cbmc_bendgen_code .eq. CBMC_IDEAL .or.
     &                 bendstyle .eq. ANG_STYLE_NONE ) then
c                 --- generating on the ideal (sin) distribution
                  angleadist = DIST_SINE
c                 --- no peaks, with no data points
                  call twh_cbmcbendcount(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,0)
               elseif ( cbmc_bendgen_code .eq. CBMC_GLOBAL_GAUSSIAN
     &                 ) then
c                  --- gaussian, but uses the equilibrium bond angle
c                  --- and a single global standard deviation
                  angleadist = DIST_GAUSSIAN
c                 --- one peak
                  ipeak = 1
                  call twh_cbmcbendcount(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak)
c                 --- check that this is an angle style that has
c                 --- an equilbrium angle
                  if ( bendstyle .eq. ANG_STYLE_NONBOND ) then
                     write(6,*) 'ATD_CBMC_SETUP: cannot combine'
     &                    ,' Global Gaussian with nonbond style angles'
                     lfailure = .true.
                     return
                  endif
c                 --- assign the mean
                  call twh_bencoeff(GLB_GET,bendtype,0,mean)
                  call twh_cbmcbendinfo(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak,1,mean)
c                 --- assign the standard deviation
                  call twh_sdevbena(GLB_GET,std_dev)
                  call twh_cbmcbendinfo(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak,2,std_dev)
c                 --- sample on the full (0, Pi)
                  call twh_cbmcbendinfo(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak,3,0.0d0)
                  call twh_cbmcbendinfo(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak,4,onepi)

               elseif ( cbmc_bendgen_code .eq. CBMC_AUTOFIT_GAUSSIAN
     &                 .or. cbmc_bendgen_code 
     &              .eq. CBMC_IDEAL_PLUS_AUTOFIT_GAUSSIAN
     &              ) then
                  angleadist = DIST_GAUSSIAN
c                 --- fetch the appropriate bond lengths for the fit
                  call twh_cbmc_fitbond(lfailure,imolty,iunit,junit
     &                 ,distone)
                  if ( lfailure ) return
                  call twh_cbmc_fitbond(lfailure,imolty,junit,kunit
     &                 ,disttwo)
                  if ( lfailure ) return
                  imax = 1000
                  sum_w = 0.0d0
                  sum_x = 0.0d0
                  sum_xx = 0.0d0
                  do itest = 0,imax
                     theta = dble(itest)*onepi/dble(imax)
                     angcos = dcos(theta)
                     lovrlap = .false.
                     energy = 0.0d0
                     if ( lthreemem ) then
c                       --- compute 3-mem ring bond energy
                        length = twh_lawofcosine(distone,disttwo,angcos)
                        energy = energy + twh_vbond(lfailure,vibtype
     &                       ,length,imolty,iunit,kunit,lovrlap)
                     else
                        energy = energy + twh_vangle(lfailure,lovrlap,1
     &                       ,imolty,iunit,ibend,0,angcos,distone
     &                       ,disttwo,xbacord,ybacord,zbacord)
                     endif
                     if ( lfailure ) return
                     if ( lovrlap ) then
                        weight = 0.0d0
                     else
                        weight = twh_expon(-beta*energy)*dsin(theta)
                     endif
                     sum_w = sum_w + weight
                     sum_x = sum_x + weight*theta
                     sum_xx = sum_xx + weight*theta*theta
                  enddo
c                 --- average both quantities
                  sum_x = sum_x / sum_w
                  sum_xx = sum_xx / sum_w
c                 --- one peak
                  ipeak = 1
                  call twh_cbmcbendcount(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak)
c                 --- assign the mean
                  call twh_cbmcbendinfo(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak,1,sum_x)
c                 --- assign the standard deviation
                  dvalue = benda_mult*dsqrt(sum_xx - sum_x*sum_x)
                  call twh_cbmcbendinfo(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak,2,dvalue)
c                 --- assign limits
                  call twh_cbmcbendinfo(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak,3,0.0d0)
                  call twh_cbmcbendinfo(GLB_SET,imolty,iunit,ibend
     &                 ,ibend,ipeak,4,onepi)
               else
                  write(6,*) 'ATD_CBMC_SETUP: no match in Angle A'
                  lfailure = .true.
                  return
               endif
c              --- set the angleadist
               call twh_cbmcbenddist(GLB_SET,imolty,iunit,ibend,ibend
     &              ,angleadist)
c              --- output based upon the angleadist
               if ( angleadist .eq. DIST_TRIANGLE ) then
                  distname = 'Triangle'
               elseif ( angleadist .eq. DIST_DELTA ) then
                  distname = 'Delta'
               elseif ( angleadist .eq. DIST_SINE ) then
                  distname = 'Ideal Sine'
               elseif ( angleadist .eq. DIST_BOUNDED_SINE ) then
                  distname = 'Bounded Sine'
               elseif ( angleadist .eq. DIST_GAUSSIAN ) then
                  distname = 'Gaussian'
               else
                  write(6,*) 'ATD_CBMC_SETUP: unknown Angle A dist'
                  write(6,*) 'angle A Dist code:',angleadist
                  lfailure = .true.
                  return
               endif
c              --- output trial generation information
               call twh_cbmcbendcount(GLB_GET,imolty,iunit,ibend,ibend
     &              ,npeak)
               write(44,'(a7,1x,i5,1x,a6,1x,3(i5,1x),a6,1x,i3
     &              ,1x,a4,1x,a25)')
     &              'Moltyp:',imolty,'Units:',iunit,junit,kunit
     &              ,'Peaks:',npeak,'ATD:',distname
               if ( angleadist .eq. DIST_TRIANGLE ) then
                  do ipeak = 1,npeak
                     call twh_cbmcbendinfo(GLB_GET,imolty,iunit,ibend
     &                    ,ibend,ipeak,1,length)
                     write(44,'(6x,a5,1x,i3,1x,a20,1x,f9.4)')
     &                    'Peak:',ipeak,'Triangle length (A):'
     &                    ,length
                  enddo
               elseif ( angleadist .eq. DIST_DELTA ) then
                  do ipeak = 1,npeak
                     call twh_cbmcbendinfo(GLB_GET,imolty,iunit,ibend
     &                    ,ibend,ipeak,1,mean)
                     write(44,'(6x,a5,1x,i3,1x,a20,1x,f9.4)')
     &                    'Peak:',ipeak,'Allowed Angle (deg):'
     &                    ,mean*radtodeg
                  enddo
               elseif ( angleadist .eq. DIST_GAUSSIAN ) then
                  write(44,'(6x,a21,2(1x,f8.4))')
     &                 'Fitting Bond Lengths:',distone,disttwo
                  do ipeak = 1,npeak
                     call twh_cbmcbendinfo(GLB_GET,imolty,iunit,ibend
     &                    ,ibend,ipeak,1,mean)
                     call twh_cbmcbendinfo(GLB_GET,imolty,iunit,ibend
     &                    ,ibend,ipeak,2,std_dev)
                     call twh_cbmcbendinfo(GLB_GET,imolty,iunit,ibend
     &                    ,ibend,ipeak,3,lovalue)
                     call twh_cbmcbendinfo(GLB_GET,imolty,iunit,ibend
     &                    ,ibend,ipeak,4,hivalue)
                     write(44,'(6x,a5,1x,i3,1x,a17,1x,f9.4,1x
     &                    ,a14,1x,f9.4,1x,a16,1x,f9.4,1x,a17,1x,f9.4)')
     &                    'Peak:',ipeak,'Angle Mean (deg):'
     &                    ,mean*radtodeg,'Std Dev (deg):'
     &                    ,std_dev*radtodeg
     &                    ,'Low Bound (deg):',lovalue*radtodeg
     &                    ,'High Bound (deg):',hivalue*radtodeg
                  enddo
               elseif ( angleadist .eq. DIST_BOUNDED_SINE ) then
                  do ipeak = 1,npeak
                     call twh_cbmcbendinfo(GLB_GET,imolty,iunit,ibend
     &                    ,ibend,ipeak,1,lovalue)
                     call twh_cbmcbendinfo(GLB_GET,imolty,iunit,ibend
     &                    ,ibend,ipeak,2,hivalue)
                     write(44,'(6x,a5,1x,i3,1x,a14,1x,f9.4
     &                    ,1x,a15,1x,f9.4)')
     &                    'Peak:',ipeak
     &                    ,'Low Bound (A):',lovalue
     &                    ,'High Bound (A):',hivalue
                  enddo
               endif

            enddo
         enddo
      enddo


      if ( cbmc_bendgen_code .eq. CBMC_AUTOFIT_GAUSSIAN .or.
     &     cbmc_bendgen_code .eq. CBMC_IDEAL_PLUS_AUTOFIT_GAUSSIAN
     &     ) then
         write(44,'(a46)') 
     &        'Determining Autofit Gaussian bend B parameters'
         do icompute = 1,ncompute
            imolty = computelist(icompute)
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
               do ibend = 1,twh_get_inben(imolty,iunit)
                  junit = twh_get_ijben2(imolty,iunit,ibend)
                  kunit = twh_get_ijben3(imolty,iunit,ibend)
c                 --- fetch the appropriate bond lengths
                  call twh_cbmc_fitbond(lfailure,imolty,junit,kunit
     &                 ,distone)
                  if ( lfailure ) return
                  itype = twh_get_itben(imolty,iunit,ibend)
                  thetaone = twh_get_bencoeff(itype,0)
                  phione = 0.0d0
                  do jbend = 1,twh_get_inben(imolty,iunit)
                     junitb = twh_get_ijben2(imolty,iunit,jbend)
                     kunitb = twh_get_ijben3(imolty,iunit,jbend)
                     if ( ibend .ne. jbend .and.
     &                    junit .eq. junitb ) then
                        call twh_cbmc_fitbond(lfailure,imolty,junit
     &                       ,kunitb,disttwo)
                        if ( lfailure ) return
                        itype = twh_get_itben(imolty,iunit,jbend)
                        thetatwo = twh_get_bencoeff(itype,0)
c                       --- get the values for the combo-angle
                        icrossbend = 0
                        do kbend = 1,twh_get_inben(imolty,kunit)
                           if ( twh_get_ijben2(imolty,kunit,kbend) 
     &                          .eq. junit
     &                          .and. twh_get_ijben3(imolty,kunit,kbend)
     &                          .eq. kunitb ) then
                              icrossbend = kbend
                           endif
                        enddo
                        if ( icrossbend .eq. 0 ) then
                           write(6,*)
     &                          'ATD_CBMC_SETUP: problem finding bend'
                           lfailure = .true.
                           return
                        endif
c                       --- initialize extrema
                        nextrema = 0
c                       --- check to see if we are sampling a rigid bend
                        type = twh_get_itben(imolty,kunit,icrossbend)
                        if ( twh_get_anglestyle(type)
     &                       .eq. ANG_STYLE_RIGID ) then
c                          --- rigid angle, handled in CBMC move
                           lrigid = .true.
                        else
                           lrigid = .false.
                        endif

                        if ( .not. lrigid ) then
                           imax = 1000
                           do itest = -1,imax
                              phitwo = dble(itest)*twopi/dble(imax)
     &                             -onepi
                              lovrlap = .false.
                              call twh_coneangle(thetaone,phione
     &                             ,thetatwo,phitwo,angcos)
                              energy = twh_vangle(lfailure,lovrlap,1
     &                             ,imolty,kunit,icrossbend,0,angcos
     &                             ,distone,disttwo,xbacord,ybacord
     &                             ,zbacord)
                              if ( lfailure ) return
                              if ( lovrlap ) then
                                 weight = 0.0d0
                              else
                                 weight = twh_expon(-beta*energy)
                              endif
                              if ( itest .gt. 0 ) then
c                                --- see if the previous location was an
c                                --- extrema
                                 lfound = .false.
                                 if ( (prevweight(1) .gt. weight .and.
     &                                prevweight(1) .gt. prevweight(2))
     &                                .or.
     &                                (prevweight(1) .lt. weight .and.
     &                                prevweight(1) .lt. prevweight(2))
     &                                ) then
                                    lfound = .true.
                                 elseif ( prevweight(1) .eq. 0.0d0 .and.
     &                                   prevenergy(2) .lt.prevenergy(1)
     &                                   .and. energy .lt. prevenergy(1)
     &                                   ) then
                                    lfound = .true.
                                 endif
                                 if ( lfound ) then
                                    nextrema = nextrema + 1
                                    if ( nextrema .gt. maxextrema ) then
                                       write(6,*)
     &                                      'ATD_CBMC_SETUP:maxextrema'
     &                                      ,' exceeded'
                                       lfailure = .true.
                                       return
                                    endif
                                    extrema(nextrema,1) = prevloc(1)
                                    extrema(nextrema,2) = prevweight(1)
                                 endif
                              endif
c                             --- shift previous weights and locations
                              prevweight(2) = prevweight(1)
                              prevweight(1) = weight
                              prevloc(2) = prevloc(1)
                              prevloc(1) = phitwo
c                             --- shift previous energies
                              prevenergy(2) = prevenergy(1)
                              prevenergy(1) = energy
                           enddo
                           nmaxima = nextrema/2
                           if ( nmaxima .eq. 0 ) then
                              write(6,*)
     &                             'ATD_CBMC_SETUP: insufficient '
     &                             ,'extrema found for bends'
                              write(6,*) 'imolty,iunit,ibend',imolty
     &                             ,iunit,ibend
                              write(6,*) 'nextrema:',nextrema
                              do iextrema = 1,nextrema
                                 write(6,*) 'location, weight'
     &                                ,extrema(iextrema,1)
     &                                ,extrema(iextrema,2)
                              enddo
                              lfailure = .true.
                              return
                           endif
                           if ( nmaxima .gt. CBMC_MAXBENDPEAKS )then
                              write(6,*) 'ATD_CBMC_SETUP: nmaxima >'
     &                             ,' CBMC_MAXBENDPEAKS'
                              write(6,*) 'nmaxima',nmaxima
                              write(6,*) 'CBMC_MAXBENDPEAKS'
     &                             ,CBMC_MAXBENDPEAKS
                              lfailure = .true.
                              return
                           endif
                           iextrema = 0
                           do ipeak = 1,nmaxima
                              iextrema = iextrema + 1
                              lfound = .false.
                              do while ( .not. lfound ) 
                                 if ( iextrema .gt. 1 ) then
c                                   --- compare with previous extrema to
c                                   --- see if this is a maxima
                                    if ( extrema(iextrema,2) .gt. 
     &                                   extrema(iextrema-1,2) ) then
                                       lfound = .true.
                                    endif
                                 else
c                                   --- compare with next extrema to see
c                                   --- if this is a maxima
                                    if ( extrema(iextrema,2) .gt.
     &                                   extrema(iextrema+1,2) ) then
                                       lfound = .true.
                                    endif
                                 endif
                                 if ( .not. lfound ) then
                                    iextrema = iextrema + 1
                                 endif
                              enddo
c                             --- get the lo and hi bounding values
                              if ( iextrema .eq. 1 ) then
                                 lolimit = extrema(nextrema,1) - twopi
                                 hilimit = extrema(iextrema+1,1)
                              elseif ( iextrema .eq. nextrema ) then
                                 lolimit = extrema(nextrema-1,1)
                                 hilimit = extrema(1,1) + twopi
                              else
                                 lolimit = extrema(iextrema-1,1)
                                 hilimit = extrema(iextrema+1,1)
                              endif
                              sum_w = 0.0d0
                              sum_x = 0.0d0
                              sum_xx = 0.0d0
                              do itest = 0,imax
                                 phitwo = lolimit + dble(itest)
     &                                *(hilimit-lolimit)/dble(imax)
                                 lovrlap = .false.
                                 call twh_coneangle(thetaone,phione
     &                                ,thetatwo,phitwo,angcos)
                                 energy = twh_vangle(lfailure,lovrlap,1
     &                                ,imolty,kunit,icrossbend,0,angcos
     &                                ,distone,disttwo,xbacord,ybacord
     &                                ,zbacord)
                                 if ( lfailure ) return
                                 if ( lovrlap ) then
                                    weight = 0.0d0
                                 else
                                    weight = twh_expon(-beta*energy)
                                 endif
                                 sum_w = sum_w + weight
                                 sum_x = sum_x + weight*phitwo
                                 sum_xx = sum_xx + weight*phitwo*phitwo
                              enddo
c                             --- average both quantities
                              sum_x = sum_x / sum_w
                              sum_xx = sum_xx / sum_w
c                             --- increment count
                              call twh_cbmcbendcount(GLB_INCR,imolty
     &                             ,iunit,ibend,jbend,1)
c                             --- assign the mean
                              call twh_cbmcbendinfo(GLB_SET,imolty
     &                             ,iunit,ibend,jbend,ipeak,1,sum_x)
c                             --- assign the standard deviation
                              dvalue = bendb_mult
     &                             *dsqrt(sum_xx - sum_x*sum_x)
                              call twh_cbmcbendinfo(GLB_SET,imolty
     &                             ,iunit,ibend,jbend,ipeak,2,dvalue)
c                             --- assign the limits
                              call twh_cbmcbendinfo(GLB_SET,imolty
     &                             ,iunit,ibend,jbend,ipeak,3,lolimit)
                              call twh_cbmcbendinfo(GLB_SET,imolty
     &                             ,iunit,ibend,jbend,ipeak,4,hilimit)
                              call twh_cbmcbendcount(GLB_GET,imolty
     &                             ,iunit,ibend,jbend,count)
                              write(44,'(a7,1x,i5,1x,a6,1x,4(i5,1x)
     &                             ,a6,1x,i1
     &                             ,1x,a21,1x,2(f8.4,1x)
     &                             ,a11,1x,f7.2,1x,a8,1x,f7.2
     &                             ,1x,a10,1x,f7.2,1x,a10,1x,f7.2)') 
     &                             'Moltyp:',imolty,'Units:',iunit,junit
     &                             ,kunit,kunitb,'Count:',count
     &                             ,'Fitting Bond Lengths:'
     &                             ,distone,disttwo
     &                             ,'Angle Mean:',sum_x*radtodeg
     &                             ,'Std Dev:',dvalue*radtodeg
     &                             ,'Low Bound:',lolimit*radtodeg
     &                             ,'High Bound:',hilimit*radtodeg
                           enddo
                        endif
                     endif
                  enddo
               enddo
            enddo
         enddo
      endif

      call twh_cbmc_dihedgen_code(GLB_GET,cbmc_dihedgen_code)
      if ( cbmc_dihedgen_code .eq. CBMC_AUTOFIT_GAUSSIAN .or.
     &     cbmc_dihedgen_code .eq. CBMC_IDEAL_PLUS_AUTOFIT_GAUSSIAN
     &     ) then
c        --- get the peak weighting style
         call twh_dihedral_peak_weight_style(GLB_GET,peakstyle)
         call twh_dihedral_sdev_multiplier(GLB_GET,dihed_mult)
         write(44,'(a48)') 
     &        'Determining Autofit Gaussian dihedral parameters'
         do icompute = 1,ncompute
            imolty = computelist(icompute)
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
               do itor = 1,twh_get_intor(imolty,iunit)
                  junit = twh_get_ijtor2(imolty,iunit,itor)
                  kunit = twh_get_ijtor3(imolty,iunit,itor)
                  lunit = twh_get_ijtor4(imolty,iunit,itor)
c                 --- set the iunit position
                  xtcord(1) = 0.0d0
                  ytcord(1) = 0.0d0
                  ztcord(1) = 0.0d0
c                 --- use the iunit-junit bond length to get the
c                 --- junit coordinates
                  do ivib = 1,twh_get_invib(imolty,iunit)
                     if ( twh_get_ijvib(imolty,iunit,ivib) .eq. junit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,iunit,ivib,itype)
                        call twh_vibcoeff(GLB_GET,itype,0,length)
                     endif
                  enddo
                  xtcord(2) = length
                  ytcord(2) = 0.0d0
                  ztcord(2) = 0.0d0
c                 --- use equilibrium junit-kunit length and the
c                 --- equilibrium iunit-junit-kunit angle to get
c                 --- the kunit coordinates
                  do ivib = 1,twh_get_invib(imolty,junit)
                     if ( twh_get_ijvib(imolty,junit,ivib) .eq. kunit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,junit,ivib,itype)
                        call twh_vibcoeff(GLB_GET,itype,0,length)
                     endif
                  enddo
                  do ibend = 1,twh_get_inben(imolty,iunit)
                     if ( twh_get_ijben2(imolty,iunit,ibend) .eq. junit
     &                    .and.
     &                    twh_get_ijben3(imolty,iunit,ibend) .eq. kunit 
     &                    ) then
                        itype = twh_get_itben(imolty,iunit,ibend)
                        theta = twh_get_bencoeff(itype,0)
                     endif
                  enddo
c                 --- setup the cone
                  xxx = 1.0d0
                  yyy = 0.0d0
                  zzz = 0.0d0
                  gamma = 0.0d0
                  call twh_cone(lfailure,1,xxx,yyy,zzz,theta,gamma
     &                 ,xub,yub,zub)
                  if ( lfailure ) return
c                 --- set the angle
                  call twh_cone(lfailure,2,xxx,yyy,zzz,theta,gamma
     &                 ,xub,yub,zub)
                  if ( lfailure ) return
                  xtcord(3) = xtcord(2) + length*xub
                  ytcord(3) = ytcord(2) + length*yub
                  ztcord(3) = ztcord(2) + length*zub
c                 --- reset the cone to get ready for the next atom
                  xxx = xub
                  yyy = yub
                  zzz = zub
                  call twh_cone(lfailure,1,xxx,yyy,zzz,theta,gamma
     &                 ,xub,yub,zub)
                  if ( lfailure ) return
                  do ivib = 1,twh_get_invib(imolty,kunit)
                     if ( twh_get_ijvib(imolty,kunit,ivib) .eq. lunit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,kunit,ivib,itype)
                        call twh_vibcoeff(GLB_GET,itype,0,length)
                     endif
                  enddo
                  do ibend = 1,twh_get_inben(imolty,junit)
                     if ( twh_get_ijben2(imolty,junit,ibend) .eq. kunit
     &                    .and.
     &                    twh_get_ijben3(imolty,junit,ibend) .eq. lunit
     &                    ) then
                        itype = twh_get_itben(imolty,junit,ibend)
                        theta = twh_get_bencoeff(itype,0)
                     endif
                  enddo
c                 --- initialize extrema
                  nextrema = 0
                  imax = 1000
c                 --- initialize prevweights
                  prevweight(1) = 0.0d0
                  prevweight(2) = 0.0d0
                  do itest = -1,imax
                     gamma = dble(itest)*twopi/dble(imax) - onepi
                     call twh_cone(lfailure,2,xxx,yyy,zzz,theta,gamma
     &                    ,xub,yub,zub)
                     if ( lfailure ) return
                     xtcord(4) = xtcord(3) + length*xub
                     ytcord(4) = ytcord(3) + length*yub
                     ztcord(4) = ztcord(3) + length*zub
                     lovrlap = .false.
                     energy = twh_vtorsion(lfailure,ibox,1,imolty,iunit
     &                    ,itor,xtcord,ytcord,ztcord,phic,lovrlap)
                     if ( lfailure ) return
                     if ( lovrlap ) then
                        weight = 0.0d0
                     else
                        weight = twh_expon(-beta * energy)
                     endif
                     if ( itest .gt. 0 ) then
c                       --- see if the previous location was an extrema
                        if ( (prevweight(1) .gt. weight .and.
     &                       prevweight(1) .gt. prevweight(2))
     &                       .or.
     &                       (prevweight(1) .lt. weight .and.
     &                       prevweight(1) .lt. prevweight(2))
     &                       ) then
                           nextrema = nextrema + 1
                           if ( nextrema .gt. maxextrema ) then
                              write(6,*) 'ATD_CBMC_SETUP:maxextrema'
     &                             ,' exceeded'
                              lfailure = .true.
                              return
                           endif
                           extrema(nextrema,1) = prevloc(1)
                           extrema(nextrema,2) = prevweight(1)
#if DEBUG_ATD_CBMC_SETUP
                           write(6,*) 'extrema count,loc,weight:'
     &                          ,nextrema,prevloc(1),prevweight(1)
#endif
                        endif
                     endif
c                    --- shift previous weights and locations
                     prevweight(2) = prevweight(1)
                     prevweight(1) = weight
                     prevloc(2) = prevloc(1)
                     prevloc(1) = gamma
                  enddo
                  nmaxima = nextrema/2
                  if ( nmaxima .eq. 0 ) then
                     call twh_cbmctorcount(GLB_GET,imolty,iunit,itor
     &                    ,ivalue)
                     write(44,'(a7,1x,i5,1x,a6,1x,4(i5,1x),a6,1x
     &                    ,i1,1x,a5,1x,f7.2,1x,a8,1x,f7.2
     &                    ,1x,a10,1x,f7.2,1x,a10,1x,f7.2)') 
     &                    'Moltyp:',imolty,'Units:',iunit,junit
     &                    ,kunit,lunit,'Count:',ivalue
                  endif
#if DEBUG_ATD_CBMC_SETUP
                  write(6,*) 'nmaxima:',nmaxima
                  write(6,*) 'nextrema:',nextrema
#endif
                  if ( nmaxima .gt. CBMC_MAXTORPEAKS )then
                     write(6,*) 'ATD_CBMC_SETUP: nmaxima exceeds'
     &                    ,' CBMC_MAXTORPEAKS'
                     write(6,*) 'nmaxima',nmaxima
                     write(6,*) 'CBMC_MAXTORPEAKS',CBMC_MAXTORPEAKS
                     lfailure = .true.
                     return
                  endif
                  iextrema = 0
                  do ipeak = 1,nmaxima
                     iextrema = iextrema + 1
                     lfound = .false.
                     do while ( .not. lfound ) 
                        if ( iextrema .gt. 1 ) then
c                          --- compare with previous extrema to
c                          --- see if this is a maxima
                           if ( extrema(iextrema,2) .gt. 
     &                          extrema(iextrema-1,2) ) then
                              lfound = .true.
                           endif
                        elseif ( iextrema .lt. nextrema ) then
c                          --- compare with next extrema to see if
c                          --- this is a maxima
                           if ( extrema(iextrema,2) .gt.
     &                          extrema(iextrema+1,2) ) then
                              lfound = .true.
                           endif
                        endif
                        if ( .not. lfound ) then
                           iextrema = iextrema + 1
                        endif
                     enddo
c                    --- get the lo and hi bounding values
                     if ( iextrema .eq. 1 ) then
                        lolimit = extrema(nextrema,1) - twopi
                        hilimit = extrema(iextrema+1,1)
                     elseif ( iextrema .eq. nextrema ) then
                        lolimit = extrema(nextrema-1,1)
                        hilimit = extrema(1,1) + twopi
                     else
                        lolimit = extrema(iextrema-1,1)
                        hilimit = extrema(iextrema+1,1)
                     endif
                     sum_w = 0.0d0
                     sum_x = 0.0d0
                     sum_xx = 0.0d0
                     do itest = 0,imax
                        gamma = lolimit + dble(itest)
     &                       *(hilimit-lolimit)/dble(imax)
                        call twh_cone(lfailure,2,xxx,yyy,zzz,theta,gamma
     &                       ,xub,yub,zub)
                        if ( lfailure ) return
                        xtcord(4) = xtcord(3) + length*xub
                        ytcord(4) = ytcord(3) + length*yub
                        ztcord(4) = ztcord(3) + length*zub
                        lovrlap = .false.
                        energy = twh_vtorsion(lfailure,ibox,1,imolty
     &                       ,iunit,itor,xtcord,ytcord,ztcord,phic
     &                       ,lovrlap)
                        if ( lfailure ) return
                        if ( lovrlap ) then
                           weight = 0.0d0
                        else
                           weight = twh_expon(-beta*energy)
                        endif
                        sum_w = sum_w + weight
                        sum_x = sum_x + weight*gamma
                        sum_xx = sum_xx + weight*gamma*gamma
                     enddo
c                    --- average both quantities
                     sum_x = sum_x / sum_w
                     sum_xx = sum_xx / sum_w
c                    --- increment count
                     call twh_cbmctorcount(GLB_INCR,imolty,iunit,itor,1)
c                    --- assign the mean
                     call twh_cbmctorinfo(GLB_SET,imolty,iunit,itor
     &                    ,ipeak,1,sum_x)
c                    --- assign the standard deviation
                     dvalue = dihed_mult*dsqrt(sum_xx - sum_x*sum_x)
                     call twh_cbmctorinfo(GLB_SET,imolty,iunit,itor
     &                    ,ipeak,2,dvalue)
c                     --- assign the limits
                     call twh_cbmctorinfo(GLB_SET,imolty,iunit,itor
     &                    ,ipeak,3,lolimit)
                     call twh_cbmctorinfo(GLB_SET,imolty,iunit,itor
     &                    ,ipeak,4,hilimit)
c                    --- assign the raw weight
                     dvalue = sum_w / (hilimit - lolimit)
                     call twh_dihedral_peak_weights(GLB_SET,imolty
     &                    ,iunit,itor,ipeak,dvalue)

                     call twh_cbmctorcount(GLB_GET,imolty,iunit,itor
     &                    ,ivalue)
                     write(44,'(a7,1x,i5,1x,a6,1x,4(i5,1x),a6,1x
     &,i1,1x,a5,1x,f7.2,1x,a8,1x,f7.2
     &,1x,a10,1x,f7.2,1x,a10,1x,f7.2)') 
     &                    'Moltyp:',imolty,'Units:',iunit,junit
     &                    ,kunit,lunit,'Count:',ivalue
     &                    ,'Mean:',sum_x*radtodeg
     &                    ,'Std Dev:',dvalue*radtodeg
     &                    ,'Low Bound:',lolimit*radtodeg
     &                    ,'High Bound:',hilimit*radtodeg
                  enddo
c                 --- reset the peak weights
                  call twh_cbmctorcount(GLB_GET,imolty,iunit,itor
     &                 ,npeak)
                  if ( peakstyle .eq. 'uniform' ) then
c                    --- weight all peaks evenly
                     do ipeak = 1,npeak
                        dvalue = 1.0d0/dble(npeak)
                        call twh_dihedral_peak_weights(GLB_SET,imolty
     &                       ,iunit,itor,ipeak,dvalue)
                     enddo
                  elseif ( peakstyle .eq. 'isolated Boltzmann' ) then
c                    --- compute the total weight
                     sum_w = 0.0d0
                     do ipeak = 1,npeak
                        call twh_dihedral_peak_weights(GLB_GET,imolty
     &                       ,iunit,itor,ipeak,dvalue)
                        sum_w = sum_w + dvalue
                     enddo
c                    --- normalize
                     do ipeak = 1,npeak
                        call twh_dihedral_peak_weights(GLB_GET,imolty
     &                       ,iunit,itor,ipeak,dvalue)
                        dvalue = dvalue / sum_w
                        call twh_dihedral_peak_weights(GLB_SET,imolty
     &                       ,iunit,itor,ipeak,dvalue)
                     enddo
                  else
                     write(6,*) 'ATD_CBMC_SETUP unknown '
     &                    ,'dihedral_peak_weight_style'
                     write(6,*) peakstyle
                     lfailure = .true.
                     return
                  endif
#if DEBUG_ATD_CBMC_SETUP
      call twh_openfile(lfailure,99,'towhee_dist','DEFAULT','UNKNOWN'
     &                 ,'DEFAULT')
      if ( lfailure ) return
      do ipeak = 1,npeak
         call twh_cbmctorinfo(GLB_GET,imolty,iunit,itor,ipeak,3,lodim)
         call twh_cbmctorinfo(GLB_GET,imolty,iunit,itor,ipeak,4,hidim)
         write(6,*) 'set ipeak,lodim,hidim',ipeak,lodim,hidim
         do imgm = 0,100
            phidisp = (dble(imgm)/dble(100))
     &           *(hidim-lodim) + lodim
            write(6,*) 'phidisp:',phidisp
            call twh_dihedral_distribution(lfailure,.true.,.false.
     &           ,style_code
     &           ,npeak,delta_dist,lolimit,hilimit,mean,std_dev
     &           ,phidisp,dvalue,.true.,imolty,iunit,itor)
            if ( lfailure ) return
            write(6,*) 'lolimit,phidisp,hilimit',lolimit,phidisp,hilimit
            call twh_dihedral_ideal_fraction(GLB_GET,ideal_frac)
            pdf_bias = ideal_frac/twopi + (1.0d0-ideal_frac)
     &           *twh_gaussprob(phidisp,mean,std_dev
     &           ,lolimit,hilimit)*dvalue
            if ( phidisp .lt. 0.0d0 ) phidisp = phidisp + twopi
            write(99,*) phidisp*radtodeg
     &           ,pdf_bias*onepi/180.0d0
         enddo
      enddo
      close(99)
#endif
               enddo
            enddo
         enddo
      endif

c     --- two bond (and three bond that used two bond) fixed endpoint 
c     --- biasing setup
      call twh_two_bond_fe_bias_style(GLB_GET,twobbias_style)
      call twh_two_bond_bias_initial_value(GLB_GET
     &     ,twobbias_init)
      call twh_three_bond_fe_bias_style(GLB_GET,threebbias_style)
      if ( twobbias_style .eq. 'autofit gaussian'
     &     .or. threebbias_style .eq.
     &     'autofit gaussian using max and min 2-4 distance'
     &     .or. twobbias_init .eq. 'autofit gaussian' ) then
         write(44,'(a52,1x,a10)') 
     &        'Determining Autofit Gaussian two bond fixed endpoint'
     &        ,'parameters'
c        --- get the vibrange values
         call twh_two_bond_bias_vibrange(GLB_GET,1,twobvibone)
         call twh_two_bond_bias_vibrange(GLB_GET,2,twobvibtwo)
         do icompute = 1,ncompute
            imolty = computelist(icompute)
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
               do ibend = 1,twh_get_inben(imolty,iunit)
                  junit = twh_get_ijben2(imolty,iunit,ibend)
                  kunit = twh_get_ijben3(imolty,iunit,ibend)
c                 --- take the equilibrium bond lengths for the angle
c                 --- computation
                  ivib = 0
                  lfound = .false.
                  bij_equil = 1.0d0
                  do while ( ivib .lt. twh_get_invib(imolty,iunit) .and.
     &                 .not. lfound )
                     ivib = ivib + 1
                     if ( twh_get_ijvib(imolty,iunit,ivib) .eq. junit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,iunit,ivib
     &                       ,bij_type)
                        call twh_vibcoeff(GLB_GET,itype,0,bij_equil)
                        lfound = .true.
                     endif
                  enddo
                  if ( .not. lfound ) then
                     write(6,*) 'ATD_CBMC_SETUP: bij vib not found'
                     lfailure = .true.
                     return
                  endif
                  ivib = 0
                  lfound = .false.
                  bjk_equil = 1.0d0
                  do while ( ivib .lt. twh_get_invib(imolty,junit) .and.
     &                 .not. lfound )
                     ivib = ivib + 1
                     if ( twh_get_ijvib(imolty,junit,ivib) .eq. kunit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,junit,ivib
     &                       ,bjk_type)
                        call twh_vibcoeff(GLB_GET,itype,0,bjk_equil)
                        lfound = .true.
                     endif
                  enddo
                  if ( .not. lfound ) then
                     write(6,*) 'ATD_CBMC_SETUP: bjk vib not found'
                     lfailure = .true.
                     return
                  endif
                  imax = 180
                  ibondmax = 100
                  sum_w = 0.0d0
                  sum_x = 0.0d0
                  sum_xx = 0.0d0
c                 --- scan the angle on 0,Pi
                  do itest = 0,imax
                     theta = dble(itest)*onepi/dble(imax)
                     angcos = dcos(theta)
                     lovrlap = .false.
c                    --- scan bij on two_bond_bias_range
                     do iscan_bij = 0,ibondmax
                        distone = bij_equil * ( twobvibone
     &                       + (dble(iscan_bij)/dble(ibondmax))
     &                       *(twobvibtwo - twobvibone) )
c                       --- compute the bij energy
                        bij_energy = twh_vbond(lfailure,bij_type,distone
     &                       ,imolty,iunit,junit,lovrlap)
                        if ( lfailure ) return
c                       --- scan bjk on two_bond_bias_range
                        do iscan_bjk = 0,ibondmax
                           disttwo = bjk_equil * ( twobvibone
     &                          + (dble(iscan_bjk)/dble(ibondmax))
     &                          *(twobvibtwo - twobvibone) )
c                          --- compute the bjk energy
                           bjk_energy = twh_vbond(lfailure,bjk_type
     &                          ,disttwo,imolty,junit,kunit,lovrlap)
                           if ( lfailure ) return
c                          --- compute the aijk energy
                           aijk_energy = twh_vangle(lfailure,lovrlap,1
     &                          ,imolty,iunit,ibend,0,angcos,distone
     &                          ,disttwo,xbacord,ybacord,zbacord)
                           if ( lfailure ) return
                           if ( lovrlap ) then
                              weight = 0.0d0
                           else
                              energy = bij_energy + bjk_energy 
     &                             + aijk_energy
                              weight = twh_expon(-beta*energy)
     &                             *dsin(theta)*distone*distone
     &                             *disttwo*disttwo
                           endif
c                          --- compute the ik distance
                           costheta = dcos(theta)
                           bik_distance = twh_lawofcosine(distone
     &                          ,disttwo,costheta)
                           sum_w = sum_w + weight
                           sum_x = sum_x + weight*bik_distance
                           sum_xx = sum_xx + weight*bik_distance**2
                        enddo
                     enddo
                  enddo
c                 --- average both quantities
                  sum_x = sum_x / sum_w
                  sum_xx = sum_xx / sum_w
c                 --- assign the mean
                  call twh_two_bond_bias_info(GLB_SET,imolty,iunit
     &                 ,ibend,1,sum_x)
c                 --- assign the standard deviation
                  dvalue = dsqrt(sum_xx - sum_x*sum_x)
                  call twh_two_bond_bias_info(GLB_SET,imolty,iunit
     &                 ,ibend,2,dvalue) 
                  write(44,'(a7,1x,i5,1x,a6,1x,3(i5,1x),a5,1x,f10.5
     &                 ,1x,a8,1x,f10.7)') 
     &                 'Moltyp:',imolty,'Units:',iunit,junit,kunit
     &                 ,'Mean:',sum_x
     &                 ,'Std Dev:',dvalue
               enddo
            enddo
         enddo
      endif

      call twh_three_bond_bias_initial_value(GLB_GET,threebbias_init)
      if ( threebbias_init  .eq. 'autofit gaussian' ) then
         write(44,'(a50)') 
     &        'Determining Autofit Gaussian Three Bond parameters'
         do icompute = 1,ncompute
            imolty = computelist(icompute)
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
               do itor = 1,twh_get_intor(imolty,iunit)
                  junit = twh_get_ijtor2(imolty,iunit,itor)
                  kunit = twh_get_ijtor3(imolty,iunit,itor)
                  lunit = twh_get_ijtor4(imolty,iunit,itor)
c                 --- set the iunit position
                  xtcord(1) = 0.0d0
                  ytcord(1) = 0.0d0
                  ztcord(1) = 0.0d0
c                 --- use the iunit-junit bond length to get the
c                 --- junit coordinates
                  do ivib = 1,twh_get_invib(imolty,iunit)
                     if ( twh_get_ijvib(imolty,iunit,ivib) .eq. junit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,iunit,ivib,itype)
                        call twh_vibcoeff(GLB_GET,itype,0,length)
                     endif
                  enddo
                  xtcord(2) = length
                  ytcord(2) = 0.0d0
                  ztcord(2) = 0.0d0
c                 --- use equilibrium junit-kunit length and the
c                 --- equilibrium iunit-junit-kunit angle to get
c                 --- the kunit coordinates
                  do ivib = 1,twh_get_invib(imolty,junit)
                     if ( twh_get_ijvib(imolty,junit,ivib) .eq. kunit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,junit,ivib,itype)
                        call twh_vibcoeff(GLB_GET,itype,0,length)
                     endif
                  enddo
                  do ibend = 1,twh_get_inben(imolty,iunit)
                     if ( twh_get_ijben2(imolty,iunit,ibend) .eq. junit
     &                    .and.
     &                    twh_get_ijben3(imolty,iunit,ibend) .eq. kunit
     &                    ) then
                        itype = twh_get_itben(imolty,iunit,ibend)
                        theta = twh_get_bencoeff(itype,0)
                     endif
                  enddo
c                 --- setup the cone
                  xxx = 1.0d0
                  yyy = 0.0d0
                  zzz = 0.0d0
                  gamma = 0.0d0
                  call twh_cone(lfailure,1,xxx,yyy,zzz,theta,gamma
     &                 ,xub,yub,zub)
                  if ( lfailure ) return
c                 --- set the angle
                  call twh_cone(lfailure,2,xxx,yyy,zzz,theta,gamma
     &                 ,xub,yub,zub)
                  if ( lfailure ) return
                  xtcord(3) = xtcord(2) + length*xub
                  ytcord(3) = ytcord(2) + length*yub
                  ztcord(3) = ztcord(2) + length*zub
c                 --- reset the cone to get ready for the next atom
                  xxx = xub
                  yyy = yub
                  zzz = zub
                  call twh_cone(lfailure,1,xxx,yyy,zzz,theta,gamma
     &                 ,xub,yub,zub)
                  if ( lfailure ) return
                  do ivib = 1,twh_get_invib(imolty,kunit)
                     if ( twh_get_ijvib(imolty,kunit,ivib) .eq. lunit 
     &                    ) then
                        call twh_itvib(GLB_GET,imolty,kunit,ivib,itype)
                        call twh_vibcoeff(GLB_GET,itype,0,length)
                     endif
                  enddo
                  do ibend = 1,twh_get_inben(imolty,junit)
                     if ( twh_get_ijben2(imolty,junit,ibend) .eq. kunit
     &                    .and.
     &                    twh_get_ijben3(imolty,junit,ibend) .eq. lunit
     &                    ) then
                        itype = twh_get_itben(imolty,junit,ibend)
                        theta = twh_get_bencoeff(itype,0)
                     endif
                  enddo
                  imax = 1000
                  sum_w = 0.0d0
                  sum_x = 0.0d0
                  sum_xx = 0.0d0
                  lolimit = -onepi
                  hilimit = onepi
                  do itest = 0,imax
                     gamma = lolimit + dble(itest)
     &                    *(hilimit-lolimit)/dble(imax)
                     call twh_cone(lfailure,2,xxx,yyy,zzz,theta,gamma
     &                    ,xub,yub,zub)
                     if ( lfailure ) return
                     xtcord(4) = xtcord(3) + length*xub
                     ytcord(4) = ytcord(3) + length*yub
                     ztcord(4) = ztcord(3) + length*zub
c                    --- compute the 1-4 distance
                     riksq = (xtcord(4)-xtcord(1))**2 
     &                    + (ytcord(4)-ytcord(1))**2
     &                    + (ztcord(4)-ztcord(1))**2
                     rik = dsqrt(riksq)
                     lovrlap = .false.
                     energy = twh_vtorsion(lfailure,ibox,1,imolty,iunit
     &                    ,itor,xtcord,ytcord,ztcord,phic,lovrlap)
                     if ( lfailure ) return
                     if ( lovrlap ) then
                        weight = 0.0d0
                     else
                        weight = twh_expon(-beta*energy)
                     endif
                     sum_w = sum_w + weight
                     sum_x = sum_x + weight*rik
                     sum_xx = sum_xx + weight*riksq
                  enddo
c                 --- average both quantities
                  sum_x = sum_x / sum_w
                  sum_xx = sum_xx / sum_w
c                 --- assign the mean
                  call twh_three_bond_bias_info(GLB_SET,imolty,iunit
     &                 ,itor,1,sum_x)
c                 --- assign the standard deviation
                  dvalue = dsqrt(sum_xx - sum_x*sum_x)
                  call twh_three_bond_bias_info(GLB_SET,imolty,iunit
     &                 ,itor,2,dvalue)
                  write(44,'(a7,1x,i5,1x,a6,1x,4(i5,1x)
     &,a5,1x,f7.2,1x,a8,1x,f7.2)') 
     &                 'Moltyp:',imolty,'Units:',iunit,junit
     &                 ,kunit,lunit,'Mean:',sum_x,'Std Dev:',dvalue
               enddo
            enddo
         enddo
      endif

c     --- close the towhee_cbmc file
      close(44)

      return
      end

      subroutine twh_cbbond_par(lfailure,imolty,iunit,ivib,dpar)
c     ******************************************************************
c     * scans the bond distribution and makes a prediction about the   *
c     * expected acceptance rate                                       *
c     *                                                                *
c     * originally written 07-16-2014 by M.G. Martin                   *
c     * last modified 07-16-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EXPON
#define FUNCTION_GAUSSPROB
#define FUNCTION_VBOND
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer imolty,iunit,ivib
      double precision dpar
c     --- local variables
      logical lovrlap,lscan
      integer junit,bonddist,imax,jmax,itest,jtest,vibtype
      double precision rmax,beta
      double precision ptrue_i,ptrue_j,pgen_i,pgen_j,paccept,sumpaccept
      double precision energy_i,energy_j,length_i,length_j
      double precision sumbias,pbias,truenorm,arbnorm
      double precision hilimit,lolimit,hivalue,lovalue
      double precision mean,std_dev

c     --- retrieve constants
      call twh_beta(GLB_GET,beta)
      call twh_max_bond_length(GLB_GET,rmax)

      call twh_ijvib(GLB_GET,imolty,iunit,ivib,junit)
      call twh_itvib(GLB_GET,imolty,iunit,ivib,vibtype)
      call twh_cbmcbondinfo(GLB_GET,imolty,iunit,ivib
     &     ,1,truenorm)
      call twh_cbmcbondinfo(GLB_GET,imolty,iunit,ivib
     &     ,2,arbnorm)
      lscan = .false.
c     --- get the bonddist
      call twh_cbbonddist(GLB_GET,imolty,iunit,ivib,bonddist)
      if ( bonddist .eq. DIST_DELTA ) then
c        --- one trial, always right, never fails
         dpar = 1.0d0
      elseif ( bonddist .eq. DIST_R_SQ ) then
         lscan = .true.
c        --- scans a range depending on lolimit and hilimit
         call twh_cbmcbondinfo(GLB_GET,imolty,iunit,ivib,3
     &        ,lovalue)
         call twh_cbmcbondinfo(GLB_GET,imolty,iunit,ivib,4
     &        ,hivalue)
         lolimit = lovalue**(1.0d0/3.0d0)
         hilimit = (hivalue+lovalue)**(1.0d0/3.0d0)
      elseif ( bonddist .eq. DIST_GAUSSIAN ) then
         lscan = .true.
c        --- scan the full range
         lolimit = 0.0d0
         hilimit = rmax
c        --- fetch the mean and std dev
         call twh_cbmcbondinfo(GLB_GET,imolty,iunit,ivib,3,mean)
         call twh_cbmcbondinfo(GLB_GET,imolty,iunit,ivib,4,std_dev)
      else
c        --- unknown bonddist
         write(6,*) 'CBBOND_PAR unknown bonddist',bonddist
         lfailure = .true.
         return
      endif

      if ( lscan ) then
         imax = 1000
         jmax = 1000
c        --- sample the acceptance rate
         sumpaccept = 0.0d0
         sumbias = 0.0d0
         do itest = 0,imax
            length_i = lolimit
     &           + (hilimit-lolimit)*dble(itest)/dble(imax)
            lovrlap = .false.
            energy_i = twh_vbond(lfailure,vibtype,length_i
     &           ,imolty,iunit,junit,lovrlap)
            if ( lfailure ) return
            if ( lovrlap ) then
               ptrue_i = 0.0d0
            else
               ptrue_i = twh_expon(-beta*energy_i)
     &              *length_i*length_i*truenorm
            endif
            if ( bonddist .eq. DIST_R_SQ ) then
               pgen_i = length_i*length_i*arbnorm
            elseif ( bonddist .eq. DIST_GAUSSIAN ) then
               pgen_i = twh_gaussprob(length_i,mean,std_dev,lolimit
     &              ,hilimit)
            else
               write(6,*) 'CBBOND_PAR unknown bonddist:',bonddist
               lfailure = .true.
               return
            endif
            do jtest = 0,jmax
               length_j = lolimit
     &              + (hilimit-lolimit)*dble(jtest)/dble(jmax)
               lovrlap = .false.
               energy_j = twh_vbond(lfailure,vibtype,length_j
     &              ,imolty,iunit,junit,lovrlap)
               if ( lfailure ) return
               if ( lovrlap ) then
                  ptrue_j = 0.0d0
               else
                  ptrue_j = twh_expon(-beta*energy_j)
     &                 *length_j*length_j*truenorm
               endif
               if ( bonddist .eq. DIST_R_SQ ) then
                  pgen_j = length_j*length_j*arbnorm
               elseif ( bonddist .eq. DIST_GAUSSIAN ) then
                  pgen_j = twh_gaussprob(length_j,mean,std_dev,lolimit
     &                 ,hilimit)
               else
                  write(6,*) 'CBBOND_PAR unknonwn bonddist:',bonddist
                  lfailure = .true.
                  return
               endif
               if ( pgen_i*ptrue_j .gt. 0.0d0 ) then
                  paccept = (ptrue_i*pgen_j)/(pgen_i*ptrue_j)
               else
                  paccept = 1.0d0
               endif
               if ( paccept .gt. 1.0d0 ) paccept = 1.0d0
               pbias = pgen_i*ptrue_j
               sumbias = sumbias + pbias
               sumpaccept = sumpaccept + pbias*paccept
            enddo
         enddo
c        --- average the sums
         dpar = sumpaccept/sumbias
      endif

      return
      end

      subroutine twh_cbmc_fitbond(lfailure,imolty,iunit,junit,length)
c     ******************************************************************
c     * returns the bond length to use when fitting higher order terms *
c     *                                                                *
c     * originally written 07-21-2014 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer imolty,iunit,junit
      double precision length
c     --- local variables
      logical lfound
      character*30 cbmc_fit_strategy
      integer ivib,nvib,jtest,vibtype,vibstyle
      double precision distlo,disthi

      call twh_cbmc_fit_strategy(GLB_GET,cbmc_fit_strategy)

c     --- find the bond
      lfound = .false.
      call twh_invib(GLB_GET,imolty,iunit,nvib)
      do ivib = 1,nvib
         call twh_ijvib(GLB_GET,imolty,iunit,ivib,jtest)
         if ( jtest .eq. junit ) then
            call twh_itvib(GLB_GET,imolty,iunit,ivib,vibtype)
            lfound = .true.
         endif
      enddo
      if ( .not. lfound ) then
         write(6,*) 'CBMC_FITBOND: bond not found'
         write(6,*) 'imolty,iunit,junit:',imolty,iunit,junit
         lfailure = .true.
         return
      endif

c     --- need to determine which bond lengths to use when fitting
c     --- the angles
      if ( cbmc_fit_strategy .eq. 'equilibrium' ) then
c        --- for most cases, use the equilbrium bond lengths stored
c        --- in vibcoeff(0)
         call twh_vibstyle(GLB_GET,vibtype,vibstyle)
         if ( vibstyle .eq. BOND_STYLE_NONBOND .or.
     &        vibstyle .eq. BOND_STYLE_NONE ) then
c           --- there is no equilibrium bond length, use the max
            call twh_max_bond_length(GLB_GET,length)
         elseif ( vibstyle .eq. BOND_STYLE_SQUARE_WELL ) then
c           --- use the average of the min and max distance
            call twh_vibcoeff(GLB_GET,vibtype,1,distlo)
            call twh_vibcoeff(GLB_GET,vibtype,2,disthi)
            length = 0.5d0*(distlo+disthi)
         elseif ( vibstyle .eq. BOND_STYLE_FENE ) then
c           --- almost certainly this is only used in LJium potentials
c           --- so just use 1.0 as a crude "equilbrium" guess
            length = 1.0d0
         else
            call twh_vibcoeff(GLB_GET,vibtype,0,length)
         endif
      else
         write(6,*) 'CBMC_FITBOND: unknown cbmc_fit_strategy'
         write(6,*) cbmc_fit_strategy
         lfailure = .true.
         return
      endif

      return
      end
