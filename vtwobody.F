#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2000-2013 Marcus G. Martin                       *
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
      function twh_vtwobody(lfailure,rijsq,ntii,ntjj,cutstyle,lovrlap,
     &   vforeign, dvtwobody)
c     ******************************************************************
c     * Computes the nonbond two-body energy for a given type and      *
c     * distance                                                       *
c     * cutstyle has several different values which do different things*
c     * -1: setting up the shifted potential, only done at begining    *
c     *  0: use the normal vdw nonbonded cutoff of rcutsq              *
c     *  1: use the inner cuttoff, unless rcut < rcutin                *
c     *  3: compute 1-4 interaction using rcutin and special 1-4 params*
c     * dvtwobody returns dU/d[lambda_lj], if scp_eval_derivative true,*
c     * 0 otherwise.                                                   *
c     * vforeign returns the potential at 'foreign' lambda values,     *
c     * corresponding to a different hamiltonian at this configuration *
c     * if we are computing foreign lambdas.  It is set to 0.0         *
c     * otherwise.  (see scaled_potential.F)                           *
c     * value returned is in units of [K]                              *
c     *                                                                *
c     * originally written as vnonbond 2-18-2000 by M.G. Martin        *
c     * last update 07-21-2013 by M.G. Martin                          *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_GT
#define FUNCTION_CMP_LT
#define FUNCTION_EXPON
#define FUNCTION_EXTRACTPAIR
#define FUNCTION_GET_HBONDCOEFF
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_NBCOEFF
#define FUNCTION_GET_NATIVE_LAMBDA
#define FUNCTION_GET_FOREIGN_LAMBDA
#define FUNCTION_SCALE_ATOMS
#define FUNCTION_VTWOBODY
#include "functions.h"
c     --- debugging flag
#define DEBUG_VTWOBODY 0
c     --- variables passed to/from the function
      logical lfailure,lovrlap
      integer ntii,ntjj,cutstyle,classical_pot_num
      double precision rijsq, dvtwobody, vforeign
      dimension vforeign(MAX_FOREIGN_LAMBDA)
c     --- local variables
c     --- logical scalars
      logical lonefour,lfound,lcontinue,scp_scale_twobody
      logical scp_eval_derivative,lshift
c     --- integer scalars
      integer ntij,icoeff,ipower,hbtype,l,scp_num_foreign_lambda
      integer tabnpair,aindex,eampairnum
c     --- double precision scalars
      double precision srsix,srthree,srnine,rij,rovers,rmaxsq,adelr
      double precision rinvthree,rinvsix,rinvnine,rinvtwelve,rinvten
      double precision sigoverr,x,acklandphi,scalefactor,wca_cutsq
      double precision rijseven,aaa,bbb,ccc,ddd,onethird,cbrt_two
      double precision tpair,rminsq,rvalue,rcurrent,rdiff
      double precision avalue,powrn
c     --- functions
      double precision vtb_lennard_jones,vtb_scaled_lj

#if DEBUG_VTWOBODY
      write(6,*) 'starting VTWOBODY'
#endif

c     --- retrieve constants
      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
      call twh_rminsq(GLB_GET,rminsq)

c     --- initialize 
      twh_vtwobody = 0.0d0
      dvtwobody = 0.0d0
      lonefour = .false.
c     initialize foreign energies
      do l = 1, scp_num_foreign_lambda
         vforeign(l) = 0.0d0
      enddo

c     --- determine the cutoff
      call twh_getnbinfo(lfailure,cutstyle,rijsq,rmaxsq,lonefour)
      if ( lfailure ) return

c     --- generic check of distance vs. rmaxsq
      if ( twh_cmp_gt(rijsq, rmaxsq) ) return

c     --- get the classical potential number
      call twh_cp_number(GLB_GET,classical_pot_num)
c     -- compute the nonbond energy according to the classical potential
      if (classical_pot_num .eq. CP_LENNARD_JONES .or.
     &     classical_pot_num .eq. CP_12_6_PLUS_SOLVATION .or.
     &     classical_pot_num .eq. CP_UFF_12_6 .or.
     &     classical_pot_num .eq. CP_DRUKKER_H_BOND ) then
c        --- lennard jones 12-6
c        --- note that these nbcoeffs have been modified and 
c        --- no longer correspond to their original meaning in 
c        --- ffnonbond or towhee_ff
         if ( twh_cmp_lt(rijsq, rminsq) ) then
            lovrlap = .true.
         else
            ntij = twh_getnbtype(ntii,ntjj)
            twh_vtwobody = vtb_lennard_jones(ntij, rijsq, lonefour)
         endif

      elseif ( classical_pot_num .eq. CP_SCALED_LENNARD_JONES ) then
c        --- Same as lennard jones 12-6, but with lambda scaling
         if ( twh_cmp_lt(rijsq, rminsq) ) then
            lovrlap = .true.
         else
            ntij = twh_getnbtype(ntii,ntjj)
            call twh_scp_scale_twobody(GLB_GET,scp_scale_twobody)
            if ( lonefour .or. (.not. scp_scale_twobody)) then
c              --- no scaling for 1-4, so foreign energies same
               twh_vtwobody = vtb_lennard_jones(ntij, rijsq, lonefour)
               do l = 1, scp_num_foreign_lambda
                  vforeign(l) = twh_vtwobody
               enddo
            else
c              --- if sigma (aka nbcoeff(13)) = 0, vtwobody = 0        
               if (twh_get_nbcoeff(13, ntij) .eq. 0.0d0) then
                  twh_vtwobody = 0.0d0
               else
                  if (scp_scale_atoms(ntii, ntjj) ) then 
                     twh_vtwobody = vtb_scaled_lj(ntij, rijsq, 
     &                  vforeign, dvtwobody)
                  else
                     twh_vtwobody = vtb_lennard_jones(ntij, rijsq, 
     &                  lonefour)
c                    with no scaling, foreign lambdas have same value
                     do l = 1, scp_num_foreign_lambda
                        vforeign(l) = twh_vtwobody
                     enddo
                  endif
               endif
            endif
         endif

      elseif ( classical_pot_num .eq. CP_WEEKS_CHANDLER_ANDERSON ) then
c         --- Repulsive only form of lennard jones 12-6
         if ( twh_cmp_lt(rijsq, rminsq) ) then
            lovrlap = .true.
         else
            ntij = twh_getnbtype(ntii,ntjj)
            if ( lonefour ) then
c              --- no scaling for 1-4
               srsix = twh_get_nbcoeff(15,ntij) / (rijsq**3)
               twh_vtwobody = srsix*(srsix-1.0d0)
     &              *twh_get_nbcoeff(16,ntij) 
            else
c              --- if sigma (aka twh_get_nbcoeff(13)) = 0, vtwobody = 0
               if (twh_get_nbcoeff(13, ntij) .eq. 0.0d0) then
                  twh_vtwobody = 0.0d0
               else
                  if (scp_scale_atoms(ntii, ntjj)) then
                     call twh_constant_onethird(GLB_GET,onethird)
                     call twh_constant_cuberoottwo(GLB_GET,cbrt_two)
                     wca_cutsq = cbrt_two
     &                    *(twh_get_nbcoeff(13,ntij)**onethird)
                     if (twh_cmp_lt(rijsq, wca_cutsq)) then
                        srsix = twh_get_nbcoeff(13,ntij) / (rijsq**3)
                        twh_vtwobody = srsix*(srsix-1.0d0)*
     &                     twh_get_nbcoeff(14,ntij)
     &                     +(twh_get_nbcoeff(14,ntij) / 4.0d0)
                        call twh_lshift(GLB_GET,lshift)
                        if ( lshift ) then
                           twh_vtwobody = twh_vtwobody
     &                          - twh_get_nbcoeff(0,ntij)
                        endif
                     else 
                        twh_vtwobody = 0.0d0
c                    --- 0 if rijsq .ge. wca_cutsq
                     endif
                  else
c                    --- if not lscale then normal LJ interaction
                     srsix = twh_get_nbcoeff(13,ntij) / (rijsq**3)
                     twh_vtwobody = srsix*(srsix-1.0d0)
     &                    *twh_get_nbcoeff(14,ntij)
                     call twh_lshift(GLB_GET,lshift)
                     if ( lshift ) then
                        twh_vtwobody = twh_vtwobody
     &                       -twh_get_nbcoeff(0,ntij)
                     endif
                  endif
               endif
            endif
         endif
      elseif ( classical_pot_num .eq. CP_12_6_PLUS_12_10_H_BOND ) then
         rinvsix = 0.0d0
c        --- 12-6 same as above, plus additional term
         if (twh_cmp_lt(rijsq, rminsq) ) then
            lovrlap = .true.
         else
            ntij = twh_getnbtype(ntii,ntjj)
            rinvsix = 1.0d0 / ( rijsq**3)
            if ( lonefour ) then
               srsix = twh_get_nbcoeff(15,ntij)*rinvsix
               twh_vtwobody = srsix*(srsix-1.0d0)
     &              *twh_get_nbcoeff(16,ntij) 
            else
               srsix = twh_get_nbcoeff(13,ntij)*rinvsix
               twh_vtwobody = srsix*(srsix-1.0d0)
     &              *twh_get_nbcoeff(14,ntij)
               call twh_lshift(GLB_GET,lshift)
               if ( lshift ) then
                  twh_vtwobody = twh_vtwobody - twh_get_nbcoeff(0,ntij)
               endif
            endif
         endif
c        --- compute the 12-10 hbond if needed
         call twh_nb_to_hbtype(GLB_GET,ntij,hbtype)
         if ( hbtype .ne. 0 ) then
            rinvtwelve = rinvsix*rinvsix
            rinvten = rinvtwelve*rijsq
            if ( lonefour ) then
c              --- use the scaling factor
               scalefactor = twh_get_hbondcoeff(hbtype,3)
            else
c              --- no scaling
               scalefactor = 1.0d0
            endif
            twh_vtwobody = twh_vtwobody + scalefactor*(
     &           twh_get_hbondcoeff(hbtype,1)*rinvtwelve
     &           - twh_get_hbondcoeff(hbtype,2)*rinvten )
         endif
      elseif( classical_pot_num .eq. CP_9_6 ) then
c        --- lennard jones 9-6 
         if (twh_cmp_lt(rijsq, rminsq) ) then
            lovrlap = .true.
         else
            ntij = twh_getnbtype(ntii,ntjj)
            rij = dsqrt(rijsq)
            srthree = twh_get_nbcoeff(1,ntij) / (rij*rijsq)
            srsix = srthree*srthree
            srnine = srsix*srthree
            twh_vtwobody = twh_get_nbcoeff(2,ntij)*((2.0d0)*srnine 
     &           - (3.0d0)*srsix )
            call twh_lshift(GLB_GET,lshift)
            if ( lshift ) then
               twh_vtwobody = twh_vtwobody - twh_get_nbcoeff(0,ntij)
            endif
         endif
      elseif ( classical_pot_num .eq. CP_12_9_6 ) then
c        --- combination 12-9-6 potential
         if (twh_cmp_lt(rijsq, rminsq) ) then
            lovrlap = .true.
         else
            ntij = twh_getnbtype(ntii,ntjj)
            rij = dsqrt(rijsq)
            rinvthree = 1.0d0/(rijsq*rij)
            rinvsix = rinvthree*rinvthree
            rinvnine = rinvsix*rinvthree
            rinvtwelve = rinvsix*rinvsix
            twh_vtwobody = twh_get_nbcoeff(1,ntij)*rinvtwelve 
     &           + twh_get_nbcoeff(2,ntij)*rinvnine
     &           + twh_get_nbcoeff(3,ntij)*rinvsix
            call twh_lshift(GLB_GET,lshift)
            if ( lshift ) then
               twh_vtwobody = twh_vtwobody - twh_get_nbcoeff(0,ntij)
            endif
         endif
      elseif ( classical_pot_num .eq. CP_DOUBLE_EXPONENTIAL ) then
c        --- double exponential
         ntij = twh_getnbtype(ntii,ntjj)
         if ( twh_cmp_lt(rijsq, rminsq) ) then
            lovrlap = .true.
         else
            aaa = twh_get_nbcoeff(1,ntij)
            bbb = twh_get_nbcoeff(2,ntij)
            ccc = twh_get_nbcoeff(3,ntij)
            ddd = twh_get_nbcoeff(4,ntij)
            rij = dsqrt(rijsq)
            twh_vtwobody = aaa * twh_expon(-bbb*rij)
     &           - ccc * twh_expon(-ddd*rij)
            call twh_lshift(GLB_GET,lshift)
            if ( lshift ) then
               twh_vtwobody = twh_vtwobody - twh_get_nbcoeff(0,ntij)
            endif
         endif
      elseif ( classical_pot_num .eq. CP_HARD_2580_MULTISTEP .or.
     &        classical_pot_num .eq. CP_REPULSIVE_2580_MULTISTEP ) then
c         --- the so-called 2,5,8,0 potential of Richard Elliot
c         --- with a hard core
         ntij = twh_getnbtype(ntii,ntjj)
         if (twh_cmp_lt(rijsq, twh_get_nbcoeff(7,ntij)) ) then
            if ( classical_pot_num .eq. CP_HARD_2580_MULTISTEP ) then
c              --- hard overlap
               lovrlap = .true.
            else
c              --- repulsive overlap
               twh_vtwobody = 1.0d10*(twh_get_nbcoeff(7,ntij)-rijsq)
            endif
         elseif (twh_cmp_lt(rijsq, twh_get_nbcoeff(11,ntij)) ) then
c           --- we are in one of the bins, figure out which one
            if (twh_cmp_lt(rijsq, twh_get_nbcoeff(8,ntij)) ) then
c              --- the 1.2 well
               twh_vtwobody = twh_get_nbcoeff(2,ntij)
            elseif (twh_cmp_lt(rijsq, twh_get_nbcoeff(9,ntij)) ) then
c              --- the 1.5 well
               twh_vtwobody = twh_get_nbcoeff(3,ntij)
            elseif (twh_cmp_lt(rijsq, twh_get_nbcoeff(10,ntij)) ) then
c              --- the 1.8 well
               twh_vtwobody = twh_get_nbcoeff(4,ntij)
            else
c              --- we must be in the final 2.0 well
               twh_vtwobody = twh_get_nbcoeff(5,ntij)
            endif
         endif
               
      elseif ( classical_pot_num .eq. CP_EXPONENTIAL_6 ) then
c        --- exponential six
c        --- need to use a special rminsq for this potential as it has
c        --- the annoying feature of going to -infinity as r goes to 
c        --- zero.  Have placed a potential dependent rminsq at the 
c        --- position of maximum energy.
         ntij = twh_getnbtype(ntii,ntjj)
         if (twh_cmp_lt(rijsq, twh_get_nbcoeff(5,ntij)) ) then
            lovrlap = .true.
         else
            rij = dsqrt(rijsq)
            twh_vtwobody = twh_get_nbcoeff(1,ntij)/(rijsq*rijsq*rijsq)
     &           + twh_get_nbcoeff(2,ntij)
     &           *twh_expon(twh_get_nbcoeff(3,ntij)*rij)
            call twh_lshift(GLB_GET,lshift)
            if ( lshift ) then
               twh_vtwobody = twh_vtwobody - twh_get_nbcoeff(0,ntij)
            endif
         endif
      elseif ( classical_pot_num .eq. CP_HARD_SPHERE ) then
c        --- hard sphere
c        --- the square of the diameter is stored in twh_get_nbcoeff(2)
         ntij = twh_getnbtype(ntii,ntjj)
         if (twh_cmp_lt(rijsq, twh_get_nbcoeff(2,ntij)) ) then
            lovrlap = .true.
         endif
      elseif ( classical_pot_num .eq. CP_REPULSIVE_SPHERE ) then
c        --- kind of hard sphere, used to equilibrate an initial 
c        --- hard sphere run, intead of having a hard overlap, it just
c        --- has a high energy
c        --- the square of the diameter is stored in nbcoeff(2)
         ntij = twh_getnbtype(ntii,ntjj)
         if (twh_cmp_lt(rijsq, twh_get_nbcoeff(2,ntij)) ) then
            twh_vtwobody = 1.0d5 + 1.0d5*(twh_get_nbcoeff(2,ntij)-rijsq)
         endif
      elseif ( classical_pot_num .eq. CP_EMBEDDED_ATOM_METHOD
     &        .or. classical_pot_num .eq. CP_EAM_PAIR_ONLY ) then
c        --- embedded atom pair potential portion
         ntij = twh_getnbtype(ntii,ntjj)
         rij = dsqrt(rijsq)
         call twh_eam_pair_number(GLB_GET,ntij,eampairnum)
         if ( eampairnum .eq. EAM_PAIR_TABLE ) then
            twh_vtwobody = twh_extractpair(lfailure,INTERP_NORMAL
     &           ,ntij,rij,lovrlap)
            if ( lfailure ) return
         elseif ( eampairnum .eq. EAM_PAIR_EXPON ) then
            twh_vtwobody = twh_get_nbcoeff(1,ntij)
     &           *twh_expon(twh_get_nbcoeff(2,ntij)*rij)
         elseif ( eampairnum .eq. EAM_PAIR_MORSE ) then
c           --- uses a 3 parameter function that is a generalized
c           --- morse potential
            adelr = twh_get_nbcoeff(2,ntij)
     &           *(rij - twh_get_nbcoeff(3,ntij))
            twh_vtwobody = twh_get_nbcoeff(1,ntij)
     &           *( twh_expon(-2.0d0*adelr) - 2.0d0*twh_expon(-adelr))
         elseif ( eampairnum .eq. EAM_PAIR_ACKLANDTP ) then
c           --- complicated pair potential that has 3 different forms
c           --- depending upon the distances
            if ( rij .eq. 0.0d0 ) then
c              --- infinite energy, set to really big to avoid divide 
c              --- by zero
               twh_vtwobody = 1.0d50
            elseif (twh_cmp_lt(rij, twh_get_nbcoeff(1,ntij)) ) then
c              --- coulombic type repulsion
               x = rij/twh_get_nbcoeff(4,ntij)
               acklandphi = 0.1818d0*twh_expon(-3.2d0*x) 
     &              + 0.5099d0*twh_expon(-0.9423d0*x)
     &              + 0.2802d0*twh_expon(-0.4029d0*x)
     &              + 0.02817d0*twh_expon(-0.2016d0*x)
               twh_vtwobody = acklandphi*twh_get_nbcoeff(3,ntij)/rij
            elseif (twh_cmp_lt(rij, twh_get_nbcoeff(2,ntij)) ) then
c              --- exponential repulsion
               x = twh_get_nbcoeff(5,ntij) 
     &              + twh_get_nbcoeff(6,ntij)*rij
     &              + twh_get_nbcoeff(7,ntij)*rijsq
     &              + twh_get_nbcoeff(8,ntij)*rij*rijsq
               twh_vtwobody = twh_get_nbcoeff(9,ntij)*twh_expon(x)
            else
c              --- cubic sum with heaviside (using a while loop to 
c              --- recreate the heaviside and save some expense)
               call twh_table_npair(GLB_GET,ntij,icoeff)
               if ( icoeff .gt. 0 ) then
                  lcontinue = .true.
               else
                  lcontinue = .false.
               endif
               do while ( lcontinue ) 
                  call twh_table_pair(GLB_GET,1,ntij,icoeff,tpair)
                  x = tpair - rij
c                 --- check that we are still in the right bounds for
c                 --- this potential
                  if ( x .gt. 0.0d0 ) then
c                    --- compute the energy
                     call twh_table_pair(GLB_GET,2,ntij,icoeff,tpair)
                     twh_vtwobody = twh_vtwobody + tpair*x*x*x
c                    --- index for the next round
                     icoeff = icoeff - 1
                     if ( icoeff .eq. 0 ) lcontinue = .false.
                  else
c                    --- out of bounds, this loop is finished
                     lcontinue = .false.
                  endif
               enddo
            endif
         elseif ( eampairnum .eq. EAM_PAIR_ACKLANDPOW ) then
c           --- complicated pair potential that has 3 different forms
c           --- depending upon the distances
            rij = dsqrt(rijsq)
            if ( rij .eq. 0.0d0 ) then
c              --- infinite energy, set to really big to avoid divide 
c              --- by zero
               twh_vtwobody = 1.0d50
            elseif (twh_cmp_lt(rij, twh_get_nbcoeff(1,ntij)) ) then
c              --- coulombic type repulsion
               acklandphi 
     &              = twh_get_nbcoeff(4,ntij)
     &              *twh_expon(twh_get_nbcoeff(5,ntij)*rij) 
     &              + twh_get_nbcoeff(6,ntij)
     &              *twh_expon(twh_get_nbcoeff(7,ntij)*rij)
     &              + twh_get_nbcoeff(8,ntij)
     &              *twh_expon(twh_get_nbcoeff(9,ntij)*rij)
     &              + twh_get_nbcoeff(10,ntij)
     &              *twh_expon(twh_get_nbcoeff(11,ntij)*rij)
               twh_vtwobody = acklandphi*twh_get_nbcoeff(3,ntij)/rij
            elseif (twh_cmp_lt(rij, twh_get_nbcoeff(2,ntij)) ) then
c              --- exponential repulsion
               x = twh_get_nbcoeff(12,ntij) 
     &              + twh_get_nbcoeff(13,ntij)*rij
     &              + twh_get_nbcoeff(14,ntij)*rijsq
     &              + twh_get_nbcoeff(15,ntij)*rij*rijsq
               twh_vtwobody = twh_get_nbcoeff(16,ntij)*twh_expon(x)
            else
c              --- power series with heaviside (using a while loop to
c              --- recreate the heaviside and save some expense)
               call twh_table_npair(GLB_GET,ntij,icoeff)
               if ( icoeff .gt. 0 ) then
                  lcontinue = .true.
               else
                  lcontinue = .false.
               endif
               do while ( lcontinue )
                  call twh_table_pair(GLB_GET,1,ntij,icoeff,tpair)
                  x = tpair - rij
                  if ( x .gt. 0.0d0 ) then
c                    --- compute the energy
                     ipower = icoeff + 3
                     call twh_table_pair(GLB_GET,2,ntij,icoeff,tpair)
                     twh_vtwobody = twh_vtwobody + tpair*(x**ipower)
c                    --- decrement index and check it is still valid
                     icoeff = icoeff - 1
                     if ( icoeff .eq. 0 ) lcontinue = .false.
                  else
c                    --- out of bounds entry, the loop is complete
                     lcontinue = .false.
                  endif
               enddo
            endif
         elseif ( eampairnum .eq. EAM_PAIR_BELASHFE ) then
c           --- Belashchenko Heaviside power series
c           --- doing this in declining distance in an attempt to only
c           --- fetch the minimum required parameters
            if (twh_cmp_gt(rij, twh_get_nbcoeff(6,ntij)) ) then
c              --- beyond the cutoff distance, zero energy
               twh_vtwobody = 0.0d0
            elseif (twh_cmp_lt(rij, twh_get_nbcoeff(1,ntij)) ) then
c              --- below r1, use a different potential
               rvalue = twh_get_nbcoeff(1,ntij)
               twh_vtwobody = twh_get_nbcoeff(7,ntij) 
     &              + twh_get_nbcoeff(8,ntij)*(rvalue-rij)
     &              + twh_get_nbcoeff(9,ntij)*
     &              (twh_expon(twh_get_nbcoeff(10,ntij)*(rvalue-rij))
     &              -1.0d0)
            else
c              --- between r1 and r6, power series
               if ( twh_cmp_lt(rij, twh_get_nbcoeff(2,ntij)) ) then
c                 --- between r1 and r2
                  rvalue = twh_get_nbcoeff(2,ntij)
                  aindex = 11
               elseif ( twh_cmp_lt(rij, twh_get_nbcoeff(3,ntij)) ) then
c                 --- between r2 and r3
                  rvalue = twh_get_nbcoeff(3,ntij)
                  aindex = 20
               elseif ( twh_cmp_lt(rij, twh_get_nbcoeff(4,ntij)) ) then
c                 --- between r3 and r4
                  rvalue = twh_get_nbcoeff(4,ntij)
                  aindex = 29
               elseif ( twh_cmp_lt(rij, twh_get_nbcoeff(5,ntij)) ) then
c                 --- between r4 and r5
                  rvalue = twh_get_nbcoeff(5,ntij)
                  aindex = 38
               else
c                 --- between r5 and r6
                  rvalue = twh_get_nbcoeff(6,ntij)
                  aindex = 47
               endif
c              --- compute the power series
               rdiff = rij - rvalue
               rcurrent = 1.0d0
               do ipower = 0,8
                  call twh_nbcoeff(GLB_GET,aindex+ipower,ntij,avalue)
                  twh_vtwobody = twh_vtwobody + avalue*rcurrent
c                 --- update rcurrent for the next interation
                  rcurrent = rcurrent*rdiff
               enddo
            endif
         elseif ( eampairnum .eq. EAM_PAIR_NONE ) then
c           --- no pair interactions
            twh_vtwobody = 0.0d0
         else
            write(6,*) 'VTWOBODY: non defined eam_pair_num'
            write(6,*) 'ntii,ntjj,ntij',ntii,ntjj,ntij
            write(6,*) 'eam_pair_num:',eampairnum
            lfailure = .true.
            return
         endif
      elseif ( classical_pot_num .eq. CP_EXPONENTIAL_12_6 ) then
c        --- combination of Buckingham (exp-6) and Lennard-Jones
c        --- note that these nbcoeffs have been modified and 
c        --- no longer correspond to their original meaning in 
c        --- ffnonbond or towhee_ff
         if (twh_cmp_lt(rijsq, rminsq )) then
            lovrlap = .true.
         else
            ntij = twh_getnbtype(ntii,ntjj)
            srsix = 1.0d0/(rijsq*rijsq*rijsq)
            rij = dsqrt(rijsq)
            twh_vtwobody = twh_get_nbcoeff(1,ntij)*srsix
     &           + twh_get_nbcoeff(2,ntij)*srsix*srsix
     &           + twh_get_nbcoeff(3,ntij)
     &           *twh_expon(twh_get_nbcoeff(4,ntij)*rij)
            call twh_lshift(GLB_GET,lshift)
            if ( lshift ) then
               twh_vtwobody = twh_vtwobody - twh_get_nbcoeff(0,ntij)
            endif
         endif
      elseif ( classical_pot_num .eq. CP_GORDON_N_6 ) then
c        --- Gordon modified n-6 potential
c        v = 4 e c(n) * [ (sigma/(r-a(n)))^n - (sigma/(r-a(n)))^6 ]
c        nb(1) stores sigma
c        nb(3) stores n
c        nb(4) stores a(n)
c        nb(5) stores 4 e c(n)
         if (twh_cmp_lt(rijsq, rminsq)) then
            lovrlap = .true.
         else
            ntij = twh_getnbtype(ntii,ntjj)
            sigoverr = twh_get_nbcoeff(1,ntij)
     &           /(dsqrt(rijsq) - twh_get_nbcoeff(4,ntij))
            call twh_nbcoeff(GLB_GET,3,ntij,powrn)
            twh_vtwobody = twh_get_nbcoeff(5,ntij)
     &           *( sigoverr**powrn - sigoverr**6 )
            call twh_lshift(GLB_GET,lshift)
            if ( lshift ) then
               twh_vtwobody = twh_vtwobody - twh_get_nbcoeff(0,ntij)
            endif
         endif
      elseif ( classical_pot_num .eq. CP_MULTIWELL ) then
c        --- multiple square wells where the first well is of infinte
c        --- extent
         ntij = twh_getnbtype(ntii,ntjj)
         rij = dsqrt(rijsq)
c        --- check hard core
         call twh_table_pair(GLB_GET,1,ntij,1,tpair)
         if ( rij .le. tpair ) then
c           --- hard overlap
            lovrlap = .true.
         else
c           --- find the first well that is valid
            icoeff = 2
            lfound = .false.
            call twh_table_npair(GLB_GET,ntij,tabnpair)
            do while ( icoeff .le. tabnpair .and. ( .not. lfound ) )
               call twh_table_pair(GLB_GET,1,ntij,icoeff,tpair)
               if ( rij .le. tpair ) then
c                 --- found our well
                  lfound = .true.
c                 --- assign energy
                  call twh_table_pair(GLB_GET,2,ntij,icoeff,tpair)
                  twh_vtwobody = tpair
               endif
c              --- increment icoeff
               icoeff = icoeff + 1
            enddo
         endif
      elseif ( classical_pot_num .eq. CP_REPULSIVE_MULTIWELL ) then
c        --- multiple square wells where the first well is a high and 
c        --- linear repulsive term designed for equilibration of 
c        --- the normal Multiwell potential
         ntij = twh_getnbtype(ntii,ntjj)
         rij = dsqrt(rijsq)
c        --- check repulsive core
         call twh_table_pair(GLB_GET,1,ntij,1,tpair)
         if ( rij .le. tpair ) then
c           --- create a large repulsive term
            call twh_table_pair(GLB_GET,1,ntij,1,tpair)
            twh_vtwobody = 1.0d5 * (tpair-rij)
         else
c           --- find the first well that is valid
            icoeff = 2
            lfound = .false.
            call twh_table_npair(GLB_GET,ntij,tabnpair)
            do while ( icoeff .le. tabnpair .and. ( .not. lfound ) )
               call twh_table_pair(GLB_GET,1,ntij,icoeff,tpair)
               if ( rij .le. tpair ) then
c                 --- found our well
                  lfound = .true.
c                 --- assign energy
                  call twh_table_pair(GLB_GET,2,ntij,icoeff,tpair)
                  twh_vtwobody = tpair
               endif
c              --- increment icoeff
               icoeff = icoeff + 1
            enddo
         endif
      elseif ( classical_pot_num .eq. CP_REPULSIVE_WELL ) then
c        --- uses the same parameters as Square Well, but designed to 
c        --- have a non-infinte, but high and linear, energy in the 
c        --- "hard wall" part of the potential
c        --- Used to equilibrate initial structures of Square Well
c        --- note that these nbcoeffs have been modified and 
c        --- no longer correspond to their original meaning in 
c        --- ffnonbond or towhee_ff
         ntij = twh_getnbtype(ntii,ntjj)
         if ( rijsq .le. twh_get_nbcoeff(1,ntij) ) then
            twh_vtwobody = 1.0d5 + 1.0d5*(twh_get_nbcoeff(1,ntij)-rijsq)
         elseif ( rijsq .le. twh_get_nbcoeff(2,ntij) ) then
            twh_vtwobody = -twh_get_nbcoeff(3,ntij)
         endif
      elseif ( classical_pot_num .eq. CP_SQUARE_WELL ) then
c        --- square well
c        --- note that these nbcoeffs have been modified and 
c        --- no longer correspond to their original meaning in 
c        --- towhee_ff
         ntij = twh_getnbtype(ntii,ntjj)
         if (twh_cmp_lt(rijsq, twh_get_nbcoeff(1,ntij)) ) then
            lovrlap = .true.
         elseif ( rijsq .le. twh_get_nbcoeff(2,ntij) ) then
            twh_vtwobody = -twh_get_nbcoeff(3,ntij)
         endif
      elseif ( classical_pot_num .eq. CP_STILLINGER_WEBER
     &        .or. classical_pot_num .eq. CP_SW_PAIR_ONLY ) then
c        --- stillinger-weber potential
         ntij = twh_getnbtype(ntii,ntjj)
         rovers = dsqrt(rijsq)/twh_get_nbcoeff(4,ntij)
         if (twh_cmp_lt(rovers, twh_get_nbcoeff(6,ntij)) ) then
            twh_vtwobody = twh_get_nbcoeff(1,ntij)
     &           *twh_get_nbcoeff(2,ntij)
     &           *((twh_get_nbcoeff(3,ntij)
     &           *(rovers**(-twh_get_nbcoeff(5,ntij)))-1.0))
     &           *twh_expon( 1.0d0 / (rovers-twh_get_nbcoeff(6,ntij)))
         endif
      elseif ( classical_pot_num .eq. CP_TABULATED_PAIR ) then
c        --- tabulated pair potential
         ntij = twh_getnbtype(ntii,ntjj)
         rij = dsqrt(rijsq)
         twh_vtwobody = twh_extractpair(lfailure,INTERP_NORMAL
     &        ,ntij,rij,lovrlap)
         if ( lfailure ) return
      elseif ( classical_pot_num .eq. CP_BUFFERED_14_7 ) then
c        --- buffered 14-7 potential
c        --- nb(1) through nb(4) contain values used to determine the
c        --- R_IJ and epsilon_IJ values
c        --- nb(5) contains R_IJ and is not used here as the constants
c        --- related to R_IJ are precomputed and stored
c        --- nb(6) contains epsilon_IJ
c        --- nb(7) contains 1.07 R_IJ
c        --- nb(8) contains 0.07 R_IJ
c        --- nb(9) contains 1.12 R_IJ^7
c        --- nb(10) contains 0.12 R_IJ^7
         ntij = twh_getnbtype(ntii,ntjj)
         rij = dsqrt(rijsq)
         rijseven = rijsq*rijsq*rijsq*rij
         twh_vtwobody = twh_get_nbcoeff(6,ntij) *
     &        ( ( twh_get_nbcoeff(7,ntij) 
     &        / (rij + twh_get_nbcoeff(8,ntij)) )**7 )
     &        * ( ( twh_get_nbcoeff(9,ntij) 
     &        / ( rijseven + twh_get_nbcoeff(10,ntij)) ) - 2.0d0 )
      else
         write(6,*)"pot_num =",classical_pot_num
         write(6,*) 'VTWOBODY: non defined classical potential'
         lfailure = .true.
         return
      endif


c        --- do scaling by lambda here.  Any potential can be scaled
c        --- linearly
c        --- by lambda_lj, although Scaled Lennard-Jones potential has a
c        --- special nonlinear scaling which is handled specially in
c        --- vtb_scaled_lj.
c        --- lambda-derivative of potential is also done here for the 
c        --- non-scaled-lj case
      call twh_scp_scale_twobody(GLB_GET,scp_scale_twobody)
      if (scp_scale_twobody) then
         if ( classical_pot_num .ne. CP_SCALED_LENNARD_JONES ) then
c           --- lambda derivative is easy for linear scaling
            call twh_scp_eval_derivative(GLB_GET,scp_eval_derivative)
            if (scp_eval_derivative) dvtwobody = twh_vtwobody
c           --- evaluate energies at foreign lambdas
            do l = 1, scp_num_foreign_lambda
               vforeign(l) = twh_vtwobody * 
     &              scp_get_foreign_lambda(LAMBDA_C, l)
            enddo
c           --- native energy is scaled by natie lambda
            twh_vtwobody = scp_get_native_lambda(LAMBDA_LJ)*twh_vtwobody
         endif
      endif

#if DEBUG_VTWOBODY
      write(6,*) "[",ntii,",",ntjj,"] twh_vtwobody = ", twh_vtwobody,
     &      " (lovrlap: ", lovrlap, ")"
#endif

      return
      end

      subroutine twh_getnbinfo(lfailure,cutstyle,rijsq,rmaxsq,lonefour)
c     ******************************************************************
c     * determines the maximum cutoff squared to use given the cutstyle*
c     * also set the one-four logical                                  *
c     *                                                                *
c     * originally split out of vtwobody 12-15-2004 by M.G. Martin     *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lonefour,ldualcutoff
      integer cutstyle
      double precision rijsq,rmaxsq,rcut,rcutin

c     --- retrieve constants
      call twh_rcut(GLB_GET,rcut)

      if ( cutstyle .eq. -1 ) then
c        --- we are setting up the shifted potential - the distance
c        --- we are passed is exactly equal to rcut.  Set rmaxsq to 
c        --- be larger than this distance to avoid roundoff trouble
         rmaxsq = rijsq + 0.1d0
      elseif ( cutstyle .eq. 0 ) then
c        --- use the normal non-coulombic nonbonded cut-off
         call twh_rcutsq(GLB_GET,rmaxsq)
      elseif ( cutstyle .eq. 1 ) then
         call twh_ldualcutoff(GLB_GET,ldualcutoff)
         if ( ldualcutoff ) then
c           --- use inner cutoff unless it is larger than the rcut
            call twh_rcutin(GLB_GET,rcutin)
            if ( rcutin .lt. rcut ) then
               call twh_rcutinsq(GLB_GET,rmaxsq)
            else
               call twh_rcutsq(GLB_GET,rmaxsq)
            endif
         else
c           --- not using dual cutoff so use rcut
            call twh_rcutsq(GLB_GET,rmaxsq)
         endif
      elseif ( cutstyle .eq. 3 ) then
c        --- one-four interactions, only special for 12-6 LJ
         call twh_rcutsq(GLB_GET,rmaxsq)
         lonefour = .true.
      else
         write(6,*) 'GETNBINFO: cutstyle invalid',cutstyle
         lfailure = .true.
         return
      endif

      return
      end


      double precision function vtb_lennard_jones(ntij, rijsq, lonefour)
c     ******************************************************************
c     * calculates Lennard-Jones twobody potential                     *
c     *                                                                *
c     * originally split out of vtwobody 03-01-2006 by MAW             *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NBCOEFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lonefour
      integer ntij
      double precision rijsq
c     --- local variables
      logical lshift
      double precision srsix

      if ( lonefour ) then
         srsix = twh_get_nbcoeff(15,ntij) / (rijsq**3)
         vtb_lennard_jones = srsix*(srsix-1.0d0)
     &        *twh_get_nbcoeff(16,ntij) 

c    --- Frank Heilmann, BASF AG
c    --- a possible bug? 1-4 interaction not shifted in 
c    --- current release version?         
c    --- --------------------------------------------------------
c    --- email contact with Marcus: it is a feature :)
c    --- --------------------------------------------------------
c    --- I have it shifted here, as I am going to use a 
c    --- 5bead/spring model (freely jointed chain) for hexadecane
         call twh_lshift(GLB_GET,lshift)
         if ( lshift ) then
            vtb_lennard_jones = vtb_lennard_jones 
     &           - twh_get_nbcoeff(0,ntij)
         endif
      else
         srsix = twh_get_nbcoeff(13,ntij) / (rijsq**3)
         vtb_lennard_jones = srsix*(srsix-1.0d0)
     &        *twh_get_nbcoeff(14,ntij)
         call twh_lshift(GLB_GET,lshift)
         if ( lshift ) then
            vtb_lennard_jones = vtb_lennard_jones 
     &           - twh_get_nbcoeff(0,ntij)
         endif
      endif

      return
      end


      double precision function vtb_scaled_lj(ntij, rijsq, 
     &    vforeign, dvtwobody)
c     ******************************************************************
c     * calculates the Scaled Lennard-Jones twobody potential and, if  *
c     * necessary, its lambda derivative                               *
c     * also calculates the foreign energy -- energy at different      *
c     * lambda_lj values                                               *
c     *                                                                *
c     * See, e.g. "Extremely precise free energy calculations..."      *
c     * Shirts, Pitera, Swope, and Pande.  J. Chem. Phys. 119:11(2003) *
c     * pp. 5740--5761.                                                *
c     *  U = B A (A - 1), where                                        *
c     *  B = 4 epsilon lambda_lj^4                                     *
c     *  A = 1/(alpha(1-lambda_lj)^2 + (r_ij / sigma)^6)               *
c     *                                                                *
c     * originally split out of vtwobody 03-01-2006 by MAW             *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_EQ
#define FUNCTION_GET_FOREIGN_LAMBDA
#define FUNCTION_GET_NATIVE_LAMBDA
#define FUNCTION_GET_NBCOEFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ntij
      double precision rijsq, dvtwobody, vforeign
      dimension vforeign(MAX_FOREIGN_LAMBDA)
c     --- local variables
      logical scp_eval_derivative,lshift
      double precision A,B,dA,dB
      double precision rij_sigma, lambda, alpha
      integer l,scp_num_foreign_lambda

c     As defined in setmixterms,
c        * twh_get_nbcoeff(13) and twh_get_nbcoeff(15) are sigma^6 
c        * twh_get_nbcoeff(14) and twh_get_nbcoeff(16) are 4*epsilon

c     --- retrieve constants
      call twh_lshift(GLB_GET,lshift)

      dvtwobody = 0.0d0
c     rij_sigma is r_ij^6 / sigma^6
      rij_sigma = rijsq**3 / twh_get_nbcoeff(13,ntij)
      alpha = scp_get_native_lambda(ALPHA_LJ)

c     Evaluate at foreign lambdas, if necessary
      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
      do l = 1, scp_num_foreign_lambda
         lambda = scp_get_foreign_lambda(LAMBDA_LJ, l)
         B = twh_get_nbcoeff(14, ntij) * lambda**4
         A = 1.0d0 / (alpha * (1.0d0-lambda)**2 + rij_sigma)
         vforeign(l) = A*(A-1.0d0)*B 
         if ( lshift ) then
            vforeign(1) = vforeign(1) - twh_get_nbcoeff(0,ntij)
         endif
      enddo

c     Evaluate at native lambda
      lambda = scp_get_native_lambda(LAMBDA_LJ)
      B = twh_get_nbcoeff(14, ntij) * lambda**4
      A = 1.0d0 / (alpha * (1.0d0-lambda)**2 + rij_sigma)
      vtb_scaled_lj = A*(A-1.0d0)*B
      if ( lshift ) then
         vtb_scaled_lj = vtb_scaled_lj - twh_get_nbcoeff(0,ntij)
      endif

c     --- evaluate dU/dlambda, if necessary.  Derivative is 
c     --- zero in all cases this statement not reached.
      call twh_scp_eval_derivative(GLB_GET,scp_eval_derivative)
      if (scp_eval_derivative) then
         if (twh_cmp_eq(lambda, 0.0d0)) then
               dvtwobody = 0.0d0
         else
            dB = B * 4.0d0 / lambda
            dA = 2.0d0 * alpha * (1.0d0-lambda) * A**2
            dvtwobody = A*dB*(A - 1.0d0) + B*dA*(2*A - 1.0d0)
         endif
      endif

      return
      end

