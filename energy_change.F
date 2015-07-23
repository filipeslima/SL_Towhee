#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2006-2012 Marcus G. Martin                       *
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
      subroutine twh_energy_change(lfailure,activebox,nchainmove
     &     ,chaininfo,loverlap,lpartial,lgetcom,lintra,lbonded
     &     ,dv_classical,dv_intra,dv_inter,dv_ext,dv_elect,dv_solv
     &     ,dv_tor,dv_ben,dv_vib,dv_quantum,dv_foreign
     &     ,dl_inter,dl_elect)
c     ******************************************************************
c     * calculates the change in energy in a single box where an       *
c     * arbitrary number of molecules are listed to have special       *
c     * coordinates (either deleted, or the temporary coords) and the  *
c     * rest use the default CRD_REAL coordinates                      *
c     * Variable summary:                                              *
c     * lfailure: global failure logical (output)                      *
c     * activebox: the box used to compute this energy change          *
c     * nchainmove: the number of chains using special flags.  These   *
c     *             are normally the chains that are being moved       *
c     * chaininfo: two dimensional array that contains the chain number*
c     *            and special flags for all of the chains that are    *
c     *            using special flags                                 *
c     * loverlap: hard overlap logical (output)                        *
c     * lpartial is true then only compute energy for rcutin < r < rcut*
c     * lgetcom is true if we should compute a new center-of-mass for  *
c     *         the molecules.  if false then we just use the com value*
c     *         stored in xcm,ycm,zcm                                  *
c     * lintra is true if we need to compute intramolecular energies   *
c     * lbonded is true if we should compute all of the bonded terms   *
c     * dv_classical: change in total classical energy (output)        *
c     * dv_intra: change in intramolecular nonbonded energy (output)   *
c     * dv_inter: change in intermolecular nonbonded energy (output)   *
c     * dv_elect: change in coulombic energy (output)                  *
c     * dv_solv: change in solvation free energy (output)              *
c     * dv_tor: change in dihedral energy (output)                     *
c     * dv_ben: change in angle energy (output)                        *
c     * dv_vib: change in bond energy (output)                         *
c     * dv_quantum: change in the quantum energy (output)              *
c     * dv_foreign: change in classical energy for different           *
c     *           foreign lambdas (see scaled_potential.F; output)     *
c     * dl_inter: change of d[inter]/d[lambda_lj] (output)             *
c     * dl_elect: change of d[elect]/d[lambda_lj] (output)             *
c     *                                                                *
c     * rewritten from engmolec 08-04-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define DEBUG_ENERGY_CHANGE 0
c     --- variables passed to/from the subroutine
      logical lfailure,lintra,lpartial,lgetcom,lbonded
      logical loverlap
      dimension loverlap(2)
      integer activebox,nchainmove
      integer chaininfo
      dimension chaininfo(nchainmove,3)
      double precision dv_classical,dv_intra,dv_inter,dv_ext,dv_quantum
      double precision dv_elect,dv_solv,dv_vib
      double precision dv_ben,dv_tor,dl_inter,dl_elect
      double precision dv_foreign
      dimension dv_foreign(MAX_FOREIGN_LAMBDA)
c     --- local variables
      logical lcompute,ldualcutoff
      integer ioldnew,ichainmove,lambda,scp_num_foreign_lambda
      integer potentialstyle
      double precision dv_three,dv_tmp
      double precision vquantum
      dimension vquantum(2)


#if DEBUG_ENERGY_CHANGE
      write(6,*) 'start ENERGY_CHANGE'
#endif

c     --- retrieve constants
      call twh_potentialstyle(GLB_GET,potentialstyle)
c     --- initialize variables
      loverlap(1) = .false.
      loverlap(2) = .false.
      dv_classical = 0.0d0
      dv_intra = 0.0d0
      dv_inter = 0.0d0
      dv_ext = 0.0d0
      dv_elect = 0.0d0
      dv_solv = 0.0d0
      dv_tor = 0.0d0
      dv_ben = 0.0d0
      dv_vib = 0.0d0
      dv_quantum = 0.0d0
      dv_three = 0.0d0

      if ( potentialstyle .eq. POT_INTERNAL ) then

c        --- Skip this if computing just the multibody terms that are
c        --- currently missing from engatom
         call twh_ldualcutoff(GLB_GET,ldualcutoff)
         if ( .not.(lpartial .and. .not. ldualcutoff) ) then
#if DEBUG_ENERGY_CHANGE
            write(6,*) 'EC calling one-body-dv:'
#endif
            call twh_ec_one_body_dv(lfailure,activebox,nchainmove
     &           ,chaininfo,lpartial,loverlap,dv_ext,dv_solv)
            if ( lfailure .or. loverlap(CRD_NEW) 
     &           .or. loverlap(CRD_OLD)) return
#if DEBUG_ENERGY_CHANGE
            write(6,*) 'EC calling two-body-dv:'
#endif
            call twh_ec_two_body_dv(lfailure,activebox,nchainmove
     &           ,chaininfo,lpartial,lgetcom,loverlap
     &           ,dv_inter,dv_elect,dl_inter,dl_elect,dv_foreign)
            if ( lfailure .or. loverlap(CRD_NEW) 
     &           .or. loverlap(CRD_OLD)) return
         endif

#if DEBUG_ENERGY_CHANGE
         write(6,*) 'EC calling multi-body-dv:'
#endif
         call twh_ec_multi_body_dv(lfailure,activebox,nchainmove
     &        ,chaininfo,dv_three)
         if ( lfailure ) return

#if DEBUG_ENERGY_CHANGE
         write(6,*) 'EC calling intra-dv:'
#endif
         call twh_ec_intra_dv(lfailure,activebox,nchainmove,chaininfo
     &        ,lintra,lbonded,loverlap,dv_intra,dv_elect,dv_vib
     &        ,dv_ben,dv_tor)
         if ( lfailure ) return
         if (loverlap(CRD_NEW) .or. loverlap(CRD_OLD)) return

#if DEBUG_ENERGY_CHANGE
         write(6,*) '2-body dv_inter:',dv_inter
         write(6,*) '3-body dv_inter:',dv_three
#endif
c        --- add dv_three to dv_inter(bindex)
         dv_inter = dv_inter + dv_three
c        --- note that dv_intra is only computed when
c        --- the flag lintra is true
         dv_tmp = dv_ext + dv_intra + dv_vib + dv_ben + dv_tor + dv_solv
c        --- add parts of energy which do not change with lambda to 
c        --- the foreign energies.
         call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
         do lambda = 1, scp_num_foreign_lambda
            dv_foreign(lambda) = dv_foreign(lambda) + dv_tmp
         enddo

#if DEBUG_ENERGY_CHANGE
         write(6,*) 'EC calling long range elect'
#endif
         call twh_ec_long_range_elect(activebox,nchainmove
     &        ,chaininfo,dv_elect)
c        --- set the total change in classical energy
         dv_classical = dv_inter + dv_elect + dv_tmp

      endif

      if ( potentialstyle .eq. POT_EXTERNAL ) then
c        --- compute the new and old external energies
         do ioldnew = 1,2
            vquantum(ioldnew) = 0.0d0
            lcompute = .false.
c           --- check to see that at least one of the chains is in a new
c           --- location, if not then we can just use the old total
c           --- quantum energy for this portion
            do ichainmove = 1,nchainmove
               if ( chaininfo(ichainmove,ioldnew) .ne. CRD_REAL ) then
                  lcompute = .true.
               endif
            enddo
            if ( lcompute ) then
c              --- we need to compute the total quantum energy for the 
c              --- activebox using the real coordinates for everything 
c              --- except those chains listed in chaininfo
               call twh_quantumenergy(lfailure,activebox,ioldnew
     &              ,nchainmove,chaininfo,vquantum(ioldnew))
               if ( lfailure ) return
            else
c              --- if none of the chains are listed then just use the 
c              --- current energy of that box
               call twh_boxvquantum(GLB_GET,activebox,vquantum(ioldnew))
            endif
         enddo

c        --- compute the quantum energy change
         dv_quantum = vquantum(CRD_NEW) - vquantum(CRD_OLD)

      endif

#if DEBUG_ENERGY_CHANGE
         write(6,*) 'Energy Change'
         write(6,*) 'dv_classical',dv_classical
         write(6,*) 'dv_inter',dv_inter
         write(6,*) 'dv_ext',dv_ext
         write(6,*) 'dv_intra',dv_intra
         write(6,*) 'dv_elect',dv_elect
         write(6,*) 'dv_vib',dv_vib
         write(6,*) 'dv_ben',dv_ben
         write(6,*) 'dv_tor',dv_tor
         write(6,*) 'dv_solv',dv_solv
         write(6,*) 'dv_three',dv_three
         write(6,*) 'dl_inter',dl_inter
         write(6,*) 'dl_elect',dl_elect
#endif

#if DEBUG_ENERGY_CHANGE
      write(6,*) 'end ENERGY_CHANGE'
#endif

      return
      end

      subroutine twh_ec_one_body_dv(lfailure,activebox,nchainmove
     &     ,chaininfo,lpartial,loverlap,dv_ext,dv_solv)
c     ******************************************************************
c     * computes the energy change in box activebox for all one-body   *
c     * interactions involving the molecules listed in chaininfo       *
c     *                                                                *
c     * originally written 08-04-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer activebox,nchainmove
      integer chaininfo
      dimension chaininfo(nchainmove,3)
      logical lpartial,loverlap
      dimension loverlap(2)
      double precision dv_ext,dv_solv
c     --- local variables
      integer ioldnew,nfield
      double precision v_ext,v_solv
      dimension v_ext(2),v_solv(2)

c     --- no energy change in the ideal box
      if ( activebox .eq. 0 ) return

c     --- interaction energy with external field
      call twh_nfield(GLB_GET,nfield)
      if ( .not. lpartial .and. nfield .ne. 0 ) then
c        --- we do not need to compute external field energies if 
c        --- lpartial as this was already done in engatom.  Otherwise
c        --- we need to compute the external field energy (for 
c        --- translate and rotate)
         do ioldnew = 1,2
            call twh_eng_external_field(lfailure,ioldnew,activebox
     &           ,nchainmove,chaininfo,loverlap(ioldnew),v_ext(ioldnew))
            if ( lfailure ) return
            if ( loverlap(ioldnew) ) return
         enddo
         dv_ext = v_ext(CRD_NEW) - v_ext(CRD_OLD)
      endif

c     --- Calculate the solvation energy
      do ioldnew = 1,2
         v_solv(ioldnew) = 0.0d0
         call twh_eng_solvation(lfailure,ioldnew,activebox,nchainmove
     &        ,chaininfo,v_solv(ioldnew))
         if ( lfailure ) return
      enddo
      dv_solv = v_solv(CRD_NEW) - v_solv(CRD_OLD)

      return
      end

      subroutine twh_ec_two_body_dv(lfailure,activebox,nchainmove
     &     ,chaininfo,lpartial,lgetcom,loverlap
     &     ,dv_inter,dv_elect,dl_vinter,dl_velect,dv_foreign)
c     ******************************************************************
c     * compute the energy change for all of the two-body interactions *
c     * involving the molecules in chaininfo                           *
c     *                                                                *
c     * rewritten from engmolec 08-04-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer activebox,nchainmove
      integer chaininfo
      dimension chaininfo(nchainmove,3)
      logical lpartial,lgetcom,loverlap
      dimension loverlap(2)
      double precision dv_inter,dv_elect,dl_vinter,dl_velect
      double precision dv_foreign
      dimension dv_foreign(MAX_FOREIGN_LAMBDA)
c     --- local variables
      integer lambda,ioldnew,scp_num_foreign_lambda
      double precision qqfact
      double precision v_inter,v_elect,l_velect,l_vinter
      dimension v_inter(2),v_elect(2),l_velect(2),l_vinter(2)
      double precision v_foreign
      dimension v_foreign(2,MAX_FOREIGN_LAMBDA)
      double precision vf_temp
      dimension vf_temp(MAX_FOREIGN_LAMBDA)

c     --- no energy changes in the ideal box
      if ( activebox .eq. 0 ) return

c     --- retrieve constants
      call twh_constant_qqfact(GLB_GET,qqfact)
      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)

c     --- initialize
      do ioldnew = CRD_NEW,CRD_OLD,(CRD_OLD - CRD_NEW)
         v_inter(ioldnew) = 0.0d0
         v_elect(ioldnew) = 0.0d0
         l_vinter(ioldnew) = 0.0d0
         l_velect(ioldnew) = 0.0d0
         do lambda = 1,scp_num_foreign_lambda
            v_foreign(ioldnew,lambda) = 0.0d0
         enddo
      enddo

      do ioldnew = CRD_NEW,CRD_OLD,(CRD_OLD - CRD_NEW)
         do lambda = 1,scp_num_foreign_lambda
            vf_temp(lambda) = 0.0d0
         enddo
         call twh_eng_two_body(lfailure,ioldnew,activebox,nchainmove
     &        ,chaininfo,lpartial,lgetcom,loverlap(ioldnew)
     &        ,v_inter(ioldnew),v_elect(ioldnew),l_vinter(ioldnew)
     &        ,l_velect(ioldnew),vf_temp)
         if ( lfailure .or. loverlap(ioldnew) ) return
         do lambda = 1,scp_num_foreign_lambda
            v_foreign(ioldnew,lambda) = vf_temp(lambda)
         enddo
      enddo

#if DEBUG_ENERGY_CHANGE
      write(6,*) 'EC_TWO_BODY_DV v_inter(new):',v_inter(CRD_NEW)
      write(6,*) 'EC_TWO_BODY_DV v_inter(old):',v_inter(CRD_OLD)
#endif
      dv_inter = v_inter(CRD_NEW) - v_inter(CRD_OLD)
      dv_elect = qqfact*(v_elect(CRD_NEW) - v_elect(CRD_OLD))
      dl_vinter = l_vinter(CRD_NEW) - l_vinter(CRD_OLD)
      dl_velect = l_velect(CRD_NEW) - l_velect(CRD_OLD)
      do lambda = 1,scp_num_foreign_lambda
         dv_foreign(lambda) = v_foreign(CRD_NEW,lambda) 
     &        - v_foreign(CRD_OLD,lambda)
      enddo

      return
      end

      subroutine twh_ec_multi_body_dv(lfailure,activebox,nchainmove
     &     ,chaininfo,dv_three)
c     ******************************************************************
c     * computes the change in the multibody energy caused by the      *
c     * molecules in chaininfo                                         *
c     *                                                                *
c     * originally written 08-05-2006 by M.G. Martin                   *
c     * last modified 11-06-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer activebox,nchainmove
      integer chaininfo
      dimension chaininfo(nchainmove,3)
      double precision dv_three
c     --- local variables
      integer classical_pot_num

c     --- no energy change in the ideal box
      if ( activebox .eq. 0 ) return

      call twh_cp_number(GLB_GET,classical_pot_num)
      if ( classical_pot_num .eq. CP_STILLINGER_WEBER ) then
         call twh_ec_stillinger_webber(lfailure,activebox,nchainmove
     &        ,chaininfo,dv_three)
         if ( lfailure ) return
      elseif ( classical_pot_num .eq. CP_EMBEDDED_ATOM_METHOD ) then
         call twh_ec_embedded_atom(lfailure,activebox,nchainmove
     &        ,chaininfo,dv_three)
         if ( lfailure ) return
      endif

      return
      end

      subroutine twh_ec_intra_dv(lfailure,activebox,nchainmove,chaininfo
     &     ,lintra,lbonded,loverlap,dv_intra,dv_elect,dv_vib
     &     ,dv_ben,dv_tor)
c     ******************************************************************
c     * Obtain the intramolecular energy difference for a proposed     *
c     * move.                                                          *
c     *                                                                *
c     * Arguments modified by this subroutine:                         *
c     *  dv_intra                                                      *
c     *  dv_elect                                                      *
c     *  dv_vib                                                        *
c     *  dv_ben                                                        *
c     *  dv_tor                                                        *
c     *  lintra                                                        *
c     * these energy changes are added to the energy changes passed    *
c     * into this subroutine, so initialization needs to be performed  *
c     * before calling this routine.  This enables the addition of     *
c     * energy terms that appear in both the intermolecular and        *
c     * intramolecular portions                                        *
c     *                                                                *
c     * Error conditions:                                              *
c     *   returns early of loverlap is true                            *
c     *                                                                *
c     * rewritten from engmolec 08-25-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer activebox,nchainmove
      integer chaininfo
      dimension chaininfo(nchainmove,3)
      logical lintra,lbonded
      logical loverlap
      dimension loverlap(2)
      double precision dv_intra,dv_elect,dv_vib,dv_ben,dv_tor
c     --- local variables
      integer ichain,imove,iflag,testcrd
      double precision qqfact
      double precision vintra,velect,vvib,vbend,vtor
      dimension vintra(2),velect(2),vvib(2),vbend(2),vtor(2)

#if DEBUG_ENERGY_CHANGE
      write(6,*) 'start TWH_EC_INTRA_DV'
      write(6,*) 'lintra,lbonded:',lintra,lbonded
#endif

c     --- retrieve constant
      call twh_constant_qqfact(GLB_GET,qqfact)

      do imove = 1,nchainmove
         ichain = chaininfo(imove,3)
         do iflag = CRD_NEW,CRD_OLD,-1
c           --- initialize energies
            vintra(iflag) = 0.0d0
            velect(iflag) = 0.0d0
            vvib(iflag) = 0.0d0
            vbend(iflag) = 0.0d0
            vtor(iflag) = 0.0d0
            testcrd = chaininfo(imove,iflag)
            if ( testcrd .ne. CRD_NONE ) then
               if ( lintra ) then
c                 --- compute the nonbonded intramolecular energy
                  call twh_eng_nonbonded_intra(lfailure,ichain,activebox
     &                 ,testcrd,loverlap(iflag),vintra(iflag)
     &                 ,velect(iflag))
                  if ( lfailure ) return
                  if (loverlap(iflag)) return
               endif

               if ( lbonded ) then
c                 --- compute all of the bonded terms for this molecule
c                 --- bond vibrations
                  call twh_eng_bond_vibrations(lfailure,ichain,activebox
     &                 ,testcrd,loverlap(iflag),vvib(iflag))
                  if ( lfailure ) return
                  if (loverlap(iflag)) return

c                 --- molecule with bending angles
                  call twh_eng_bending(lfailure,ichain,activebox,testcrd
     &                 ,loverlap(iflag),vbend(iflag))
                  if ( lfailure ) return
                  if (loverlap(iflag)) return
      
c                 --- bond bond cross terms
                  call twh_eng_bond_bond(ichain,activebox,testcrd
     &                 ,vvib(iflag))

c                 --- regular torsions 
                  call twh_eng_regular_torsions(lfailure,ichain
     &                 ,activebox,testcrd,loverlap(iflag),vtor(iflag))
                  if ( lfailure ) return
                  if (loverlap(iflag)) return

c                 --- angle angle cross terms
                  call twh_eng_angle_angle(ichain,activebox,testcrd
     &                 ,vbend(iflag))

c                 --- improper torsions 
                  call twh_eng_improper_torsions(ichain,activebox
     &                 ,testcrd,vtor(iflag))
               endif
            endif
         enddo

         dv_intra = dv_intra + vintra(CRD_NEW) - vintra(CRD_OLD)
         dv_elect = dv_elect + qqfact*(velect(CRD_NEW)-velect(CRD_OLD))
         dv_vib = dv_vib + vvib(CRD_NEW) - vvib(CRD_OLD)
         dv_ben = dv_ben + vbend(CRD_NEW) - vbend(CRD_OLD)
         dv_tor = dv_tor + vtor(CRD_NEW) - vtor(CRD_OLD)

      enddo

#if DEBUG_ENERGY_CHANGE
      write(6,*) 'finish TWH_EC_INTRA_DV'
#endif
      return
      end

      subroutine twh_ec_long_range_elect(activebox,nchainmove
     &     ,chaininfo,dv_elect)
c     ******************************************************************
c     * computes the change in the long-range electrostatics in the    *
c     * activebox for the molecules listed in chaininfo                *
c     *                                                                *
c     * originally written 10-06-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer activebox,nchainmove
      integer chaininfo
      dimension chaininfo(nchainmove,3)
      double precision dv_elect
c     --- local variables
      logical lewald
      integer ic,imolty,ncount,icrd,iunit,testcrd,imove,ichain,nunit
      integer numvect
      double precision arg,kxfetch,kyfetch,kzfetch,xtemp,ytemp,ztemp
      double precision charge,ssumr,ssumrn,ssumi,ssumin
c     --- double precision arrays
      double precision sumr,sumi,vrecip
      dimension sumr(2),sumi(2),vrecip(2)

c     --- retrieve constants
      call twh_lewald(GLB_GET,lewald)

      if ( lewald ) then
c        --- bail out if this is the ideal chain box
         if ( activebox .eq. 0 ) return
         do icrd = CRD_OLD,CRD_NEW,( CRD_NEW - CRD_OLD )
            vrecip(icrd) = 0.0d0
         enddo

         call twh_ewald_numvect(GLB_GET,activebox,numvect)
         ncount = numvect
         do ic = 1, ncount
c           --- fetch the kx, ky, and kz into local memory
            call twh_ewald_kvector(GLB_GET,CRD_REAL,ic,activebox
     &           ,kxfetch,kyfetch,kzfetch)
            do icrd = CRD_OLD,CRD_NEW,( CRD_NEW - CRD_OLD )
               sumr(icrd) = 0.0d0
               sumi(icrd) = 0.0d0
               do imove = 1,nchainmove
                  testcrd = chaininfo(imove,icrd)
                  if ( testcrd .ne. CRD_NONE ) then
                     ichain = chaininfo(imove,3)
                     call twh_moltyp(GLB_GET,ichain,imolty)
                     call twh_nunit(GLB_GET,imolty,nunit)
                     do iunit = 1,nunit
                        call twh_coordinates(GLB_GET,testcrd,ichain
     &                       ,iunit,xtemp,ytemp,ztemp)
                        arg = kxfetch*xtemp
     &                       + kyfetch*ytemp
     &                       + kzfetch*ztemp
                        call twh_qqatom(GLB_GET,imolty,iunit,charge)
                        sumr(icrd) = sumr(icrd) + charge*dcos(arg)
                        sumi(icrd) = sumi(icrd) + charge*dsin(arg)
                     enddo
                  endif
               enddo
            enddo
            call twh_ewald_ssumr(GLB_GET,CRD_REAL,ic,activebox,ssumr)
            ssumrn = ssumr - sumr(CRD_OLD) + sumr(CRD_NEW)
            call twh_ewald_ssumr(GLB_SET,CRD_NEW,ic,activebox,ssumrn)
            call twh_ewald_ssumi(GLB_GET,CRD_REAL,ic,activebox,ssumi)
            ssumin = ssumi - sumi(CRD_OLD) + sumi(CRD_NEW)
            call twh_ewald_ssumi(GLB_SET,CRD_NEW,ic,activebox,ssumin)
         enddo

         call twh_ewald_sum(activebox,vrecip)

         dv_elect = dv_elect + vrecip(CRD_NEW) - vrecip(CRD_OLD)
      endif

      return
      end

      subroutine twh_ec_embedded_atom(lfailure,activebox,nchainmove
     &        ,chaininfo,dv_three)
c     ******************************************************************
c     * calculates the change in EAM energies                          *
c     *                                                                *
c     * Arguments set by subroutine:                                   *
c     *  dv_three                                                      *
c     *                                                                *
c     * rewritten from engmolec 10-09-2006 by M.G. Martin              *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_VEMBED
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer activebox,nchainmove
      integer chaininfo
      dimension chaininfo(nchainmove,3)
      double precision dv_three
c     --- local variables
c     --- logical scalars
      logical lzero,lovrlap,lcompute
c     --- integer scalars
      integer ntii,jchain,ntjj,imolty,jmolty,kmolty,kchain
      integer jlist,klist,ichain
      integer imove,jmove,icrd,itestcord,jtestcord,nmolty,ncmt
c     --- double precision scalars
      double precision rxui,rzui,ryui,rxuj,rzuj,ryuj,rxuij,ryuij,rzuij
      double precision rij,bij,bji,density
c     --- double precision arrays
      double precision old_eam_rho
      dimension old_eam_rho(nchainmove)

      lovrlap = .false.
      call twh_nmolty(GLB_GET,nmolty)
c     --- initialize the temp densities
      do kmolty = 1,nmolty
         call twh_ncmt(GLB_GET,activebox,kmolty,ncmt)
         do klist = 1,ncmt
            call twh_chainlist(GLB_GET,klist,activebox,kmolty,kchain)
            if ( activebox .eq. 0 ) then
c              --- ideal gas box, all zero densities
               call twh_eam_rho(GLB_SET,CRD_TEMP,kchain,0.0d0)
            else
c              --- start with the current density
               call twh_eam_rho(GLB_STORE,CRD_R_TO_T,kchain,density)
            endif
         enddo
      enddo
c     --- set the temp_eam_rho for the active molecules to 0.0
c     --- set the old_eam_rho for the active molecules to 0.0
      do imove = 1,nchainmove
         old_eam_rho(imove) = 0.0d0
         ichain = chaininfo(imove,3)
         call twh_eam_rho(GLB_SET,CRD_TEMP,ichain,0.0d0)
      enddo

c     --- contributions of the chaininfo list with each other
      do icrd = CRD_OLD,CRD_NEW,(CRD_NEW-CRD_OLD)
         do imove = 1,nchainmove-1
            ichain = chaininfo(imove,3)
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_ntype(GLB_GET,imolty,1,ntii)
            itestcord = chaininfo(imove,icrd)
            if ( itestcord .ne. CRD_NONE ) then
               call twh_coordinates(GLB_GET,itestcord,ichain,1
     &              ,rxui,ryui,rzui)
               do jmove = imove+1,nchainmove
                  jchain = chaininfo(jmove,3)
                  call twh_moltyp(GLB_GET,jchain,jmolty)
                  call twh_ntype(GLB_GET,jmolty,1,ntjj)
                  jtestcord = chaininfo(jmove,icrd)
                  if ( jtestcord .ne. CRD_NONE ) then
c                    --- compute the cross interaction
                     call twh_coordinates(GLB_GET,jtestcord,jchain,1
     &                    ,rxuj,ryuj,rzuj)
c                    --- compute i-j bond distance
                     call twh_mimage(activebox,rxuj,ryuj,rzuj
     &                    ,rxui,ryui,rzui,rxuij,ryuij,rzuij)
                     rij = twh_distance(rxuij,ryuij,rzuij)
c                    --- compute i-j bond density
                     call twh_bondorder(lfailure,ntii,ntjj,rij
     &                    ,bij,bji,lzero)
                     if ( lfailure ) return
c                    --- modify the temp densities appropriately
                     if ( icrd .eq. CRD_NEW ) then
c                       --- new densities so add to the jchain rho
                        call twh_eam_rho(GLB_INCR,CRD_TEMP,jchain,bji)
c                       --- new densities added to the ichain rho
                        call twh_eam_rho(GLB_INCR,CRD_TEMP,ichain,bij)
                     else
c                       --- old densities added to jmove rho
                        old_eam_rho(jmove) = old_eam_rho(jmove) + bji
c                       --- old densities added to imove rho
                        old_eam_rho(imove) = old_eam_rho(imove) + bij
                     endif
                  endif
               enddo
            endif
         enddo
      enddo

c     --- interactions between chaininfo and existing chains
      do jmolty = 1,nmolty
         call twh_ntype(GLB_GET,jmolty,1,ntjj)
         call twh_ncmt(GLB_GET,activebox,jmolty,ncmt)
         do jlist = 1,ncmt
            call twh_chainlist(GLB_GET,jlist,activebox,jmolty,jchain)
            call twh_moltyp(GLB_GET,jchain,jmolty)
            call twh_ntype(GLB_GET,jmolty,1,ntjj)
c           --- make sure this molecule is not in the chaininfo list
            lcompute = .true.
            do imove = 1,nchainmove
               ichain = chaininfo(imove,3)
               if ( ichain .eq. jchain ) lcompute = .false.
            enddo
            if ( lcompute ) then
c              --- get the coordinates for the jchain
               call twh_coordinates(GLB_GET,CRD_REAL,jchain,1
     &              ,rxuj,ryuj,rzuj)
c              --- run through the chainlist and compute the change in 
c              --- the embedding density of jchain
               do icrd = CRD_OLD,CRD_NEW,(CRD_NEW-CRD_OLD)
                  do imove = 1,nchainmove
                     itestcord = chaininfo(imove,icrd)
                     if (  itestcord .ne. CRD_NONE ) then
                        ichain = chaininfo(imove,3)
                        call twh_moltyp(GLB_GET,ichain,imolty)
                        call twh_ntype(GLB_GET,imolty,1,ntii)
                        call twh_coordinates(GLB_GET,itestcord,ichain,1
     &                       ,rxui,ryui,rzui)
c                       --- compute i-j bond distance
                        call twh_mimage(activebox,rxuj,ryuj,rzuj
     &                       ,rxui,ryui,rzui,rxuij,ryuij,rzuij)
                        rij = twh_distance(rxuij,ryuij,rzuij)
c                       --- compute i-j bond density
                        call twh_bondorder(lfailure,ntii,ntjj,rij
     &                       ,bij,bji,lzero)
                        if ( lfailure ) return
c                       --- modify the temp densities appropriately
                        if ( icrd .eq. CRD_NEW ) then
c                          --- new densities so add to the jchain rho
                           call twh_eam_rho(GLB_INCR,CRD_TEMP,jchain
     &                          ,bji)
c                          --- also add to the ichain rho
                           call twh_eam_rho(GLB_INCR,CRD_TEMP,ichain
     &                          ,bij)
                        else
c                          --- old densities so subtract from jchain rho
                           call twh_eam_rho(GLB_DECR,CRD_TEMP,jchain
     &                          ,bji)
c                          --- add to the old eam rho
                           old_eam_rho(imove) = old_eam_rho(imove) + bij
                        endif
                     endif
                  enddo
               enddo
c              --- safety check all of the densities for undershoots
               if ( lovrlap ) then
                  write(6,*) 'EC_EMBEDDED_ATOM: undershoot in density'
                  lfailure = .true.
                  return
               endif
c              --- add the new embedding energy for jchain
               call twh_eam_rho(GLB_GET,CRD_TEMP,jchain,density)
               dv_three = dv_three + twh_vembed(lfailure,ntjj,density)
               if ( lfailure ) return
c              --- subtract the old embedding energy for jchain
               call twh_eam_rho(GLB_GET,CRD_REAL,jchain,density)
               dv_three = dv_three - twh_vembed(lfailure,ntjj,density)
               if ( lfailure ) return
               if ( lovrlap ) then
                  write(6,*) 'EC_EMBEDDED_ATOM: overlap in embedding'
                  lfailure = .true.
                  return
               endif
            endif
         enddo
      enddo

c     --- compute the change in the ebedding energy for chaininfo
      do imove = 1,nchainmove
         ichain = chaininfo(imove,3)
         call twh_moltyp(GLB_GET,ichain,imolty)
         call twh_ntype(GLB_GET,imolty,1,ntii)
c        --- add the new embedding energy
         itestcord = chaininfo(imove,CRD_NEW)
         if ( itestcord .ne. CRD_NONE ) then
            call twh_eam_rho(GLB_GET,CRD_TEMP,ichain,density)
            dv_three = dv_three + twh_vembed(lfailure,ntii,density)
            if ( lfailure ) return
         endif
c        --- subtract the old embedding energy
         itestcord = chaininfo(imove,CRD_OLD)
         if ( itestcord .ne. CRD_NONE ) then
            density = old_eam_rho(imove)
            dv_three = dv_three - twh_vembed(lfailure,ntii,density)
            if ( lfailure ) return
         endif
      enddo

      return
      end

      subroutine twh_ec_stillinger_webber(lfailure,activebox,nchainmove
     &     ,chaininfo,dv_three)
c     ******************************************************************
c     * computes the change in energy of the multibody portion of the  *
c     * stillinger-weber potential                                     *
c     *                                                                *
c     * originally written 10-10-2006 by M.G. Martin                   *
c     * last modified 10-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer activebox,nchainmove
      integer chaininfo
      dimension chaininfo(nchainmove,3)
      double precision dv_three
c     --- local variables
      integer icord
      double precision vthree
      dimension vthree(2)

      do icord = CRD_NEW,CRD_OLD,( CRD_OLD - CRD_NEW )
         call twh_eng_stillinger_weber(lfailure,activebox,icord
     &        ,nchainmove,chaininfo,vthree(icord))
         if ( lfailure ) return
      enddo

      dv_three = vthree(CRD_NEW) - vthree(CRD_OLD)

      return
      end

      subroutine twh_ec_atomlist_bond(lfailure,ichain,activebox
     &     ,natommove,atomlist,globalcord,loverlap,dvatomlist)
c     ******************************************************************
c     * computes the change in the vibration energy for the molecule   *
c     * ichain in the activebox for all terms that involve any of the  *
c     * atoms in the atomlist                                          *
c     *                                                                *
c     * originally written 10-18-2006 by M.G. Martin                   *
c     * last modified 10-25-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "dvpm.h"
#define FUNCTION_DISTANCE
#define FUNCTION_VBOND
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ichain,activebox,natommove
      integer atomlist
      dimension atomlist(natommove,3)
      integer globalcord
      dimension globalcord(2)
      double precision dvatomlist
      logical loverlap
      dimension loverlap(2)
c     --- local variables
c     --- logical scalars
      logical lcompute,lglobal
c     --- integer scalars
      integer icord,ivib,nvib,iatom,imolty,targetatom,jatom,testatom
      integer otheratom,bondtype,atomicord,atomjcord
c     --- double precision scalars
      double precision length,rxui,ryui,rzui,rxuj,ryuj,rzuj
      double precision rxuij,ryuij,rzuij

      call twh_moltyp(GLB_GET,ichain,imolty)
      do iatom = 1,natommove
         targetatom = atomlist(iatom,3)
         call twh_invib(GLB_GET,imolty,targetatom,nvib)
         do ivib = 1,nvib
            call twh_ijvib(GLB_GET,imolty,targetatom,ivib,otheratom)
c           --- avoid double counting by making sure iatom
c           --- is the lowest index of any atoms in the atomlist
            lcompute = .true.
            jatom = 0
            do while ( jatom .lt. iatom-1 .and. lcompute )
               jatom = jatom + 1
               testatom = atomlist(jatom,3)
               if ( otheratom .eq. testatom ) lcompute = .false.
            enddo
            if ( lcompute ) then
c              --- see if the otheratom occurs later in the atomlist
c              --- it could not occur previously in the list due to the
c              --- logic above
               lglobal = .true.
               jatom = iatom
               do while ( lglobal .and. jatom .lt. natommove )
                  jatom = jatom + 1
                  testatom = atomlist(jatom,3)
                  if ( otheratom .eq. testatom ) lglobal = .false.
               enddo
               do icord = CRD_NEW,CRD_OLD,(CRD_OLD - CRD_NEW)
c                 --- determine atom positions for the first atom
c                 --- this atom is in the atomlist at position iatom
                  atomicord = atomlist(iatom,icord)
                  call twh_coordinates(GLB_GET,atomicord,ichain
     &                 ,targetatom,rxui,ryui,rzui)
c                 --- determine atom positions for the second atom
                  if ( lglobal ) then
c                    --- not in atomlist, use global coords
                     atomjcord = globalcord(icord)
                  else
c                    --- at position jatom in the atomlist
                     atomjcord = atomlist(jatom,icord)
                  endif
                  call twh_coordinates(GLB_GET,atomjcord,ichain
     &                 ,otheratom,rxuj,ryuj,rzuj)
c                 --- compute bond distance
                  call twh_mimage(activebox,rxui,ryui,rzui,rxuj,ryuj
     &                 ,rzuj,rxuij,ryuij,rzuij)
                  length = twh_distance(rxuij,ryuij,rzuij)
c                 --- compute bond energy
                  call twh_itvib(GLB_GET,imolty,targetatom,ivib
     &                 ,bondtype)
                  dvatomlist = dvatomlist + dvpm(icord)*
     &                 twh_vbond(lfailure,bondtype,length,imolty
     &                 ,targetatom,otheratom,loverlap(icord))
                  if ( lfailure ) return
                  if ( loverlap(icord) ) return
               enddo
            endif
         enddo
      enddo

      return
      end

      subroutine twh_ec_atomlist_bond_bond(ichain,activebox
     &     ,natommove,atomlist,globalcord,dvatomlist)
c     ******************************************************************
c     * computes the change in the bond-bond energy for the molecule   *
c     * ichain in the activebox for terms that involve any of the atoms*
c     * in the atomlist                                                *
c     *                                                                *
c     * originally written 10-25-2006 by M.G. Martin                   *
c     * last modified 08-30-2012 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "dvpm.h"
#define FUNCTION_DISTANCE
#define FUNCTION_IN_ATOMLIST
#define FUNCTION_GET_ATOMLIST_CORD
#define FUNCTION_VBONBON
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ichain,activebox,natommove
      integer atomlist
      dimension atomlist(natommove,3)
      integer globalcord
      dimension globalcord(2)
      double precision dvatomlist
c     --- local variables
c     --- logical scalars
      logical lcompute
c     --- integer scalars
      integer iatomlist,iatom,jatom,ibnbn,nbnbn,testlist,testunit
      integer imolty,ioldnew,targetcord,targetatom,iatomcord,jatomcord
      integer nvib,ivib,tbbone,tbbtwo,keyatom
      integer keyatomcord,oneatomcord,twoatomcord
c     --- double precision scalars
      double precision rx,ry,rz
c     --- double precision arrays
      double precision xpos,ypos,zpos
      dimension xpos(3),ypos(3),zpos(3)
      double precision length
      dimension length(2)

      call twh_moltyp(GLB_GET,ichain,imolty)

c     --- on the first pass compute all bond-bond interactions that
c     --- have an end on an atomlist atom
c     --- just need to check to make sure we do not double count when
c     --- both ends have an atomlist atom
      do iatomlist = 1,natommove
         targetatom = atomlist(iatomlist,3)
         call twh_inbnbn(GLB_GET,imolty,targetatom,nbnbn)
         do ibnbn = 1,nbnbn
            call twh_ijbnbn1(GLB_GET,imolty,targetatom,ibnbn,iatom)
            call twh_ijbnbn2(GLB_GET,imolty,targetatom,ibnbn,jatom)
c           --- make sure the jatom was not previously a target atom
            lcompute = .true.
            testlist = 0
            do while ( testlist .lt. iatomlist-1 .and. lcompute )
               testlist = testlist + 1
               testunit = atomlist(testlist,3)
               if ( testunit .eq. jatom ) lcompute = .false.
            enddo
            if ( lcompute ) then
               do ioldnew = CRD_NEW,CRD_OLD,(CRD_OLD-CRD_NEW)
c                 --- targetatom is always in the atomlist
                  targetcord = atomlist(iatomlist,ioldnew)
c                 --- determine the coordinates to use for iatom
                  iatomcord = twh_get_atomlist_cord(iatom,ioldnew
     &                 ,natommove,atomlist,globalcord)
c                 --- determine the coordinates to use for jatom
                  jatomcord = twh_get_atomlist_cord(jatom,ioldnew
     &                 ,natommove,atomlist,globalcord)
                  if ( targetcord .ne. CRD_NONE .and.
     &                 iatomcord .ne. CRD_NONE .and.
     &                 jatomcord .ne. CRD_NONE ) then
c                    --- fetch the coordinates for targetatom
                     call twh_coordinates(GLB_GET,targetcord,ichain
     &                    ,targetatom,xpos(1),ypos(1),zpos(1))
c                    --- fetch the coordinates for iatom
                     call twh_coordinates(GLB_GET,iatomcord,ichain,iatom
     &                    ,xpos(2),ypos(2),zpos(2))
c                    --- fetch the coordinates for jatom
                     call twh_coordinates(GLB_GET,jatomcord,ichain,jatom
     &                    ,xpos(3),ypos(3),zpos(3))
c                    --- compute the first distance between targetatom
c                    --- and iatom
                     call twh_mimage(activebox,xpos(1),ypos(1),zpos(1)
     &                    ,xpos(2),ypos(2),zpos(2),rx,ry,rz)
                     length(1) = twh_distance(rx,ry,rz)
c                    --- compute the second distance between iatom and
c                    --- jatom
                     call twh_mimage(activebox,xpos(1),ypos(1),zpos(1)
     &                    ,xpos(3),ypos(3),zpos(3),rx,ry,rz)
                     length(2) = twh_distance(rx,ry,rz)
c                    --- compute bond-bond energy
                     dvatomlist = dvatomlist + dvpm(ioldnew)
     &                    *twh_vbonbon(imolty,targetatom,ibnbn,length)
                  endif
               enddo
            endif
         enddo
      enddo

c     --- on the second pass we want to find bond-bond terms that have
c     --- the atomlist in the bnbn1 or bnbn2 position
c     --- to avoid double counting we must make sure that the keyatom
c     --- is not one of the other atomlist atoms (handled above)
c     --- and if both the bnbn1 and bnbn2 atoms are atomlist atoms we
c     --- only compute if our current targetatom is in the bnbn1 spot
      do iatomlist = 1,natommove
         targetatom = atomlist(iatomlist,3)
c        --- check the neighbors of the targetatom
         call twh_invib(GLB_GET,imolty,targetatom,nvib)
         do ivib = 1,nvib
c           --- get the neighboring atom
            call twh_ijvib(GLB_GET,imolty,targetatom,ivib,keyatom)
c           --- make sure this is not in the atomlist
            if ( .not. twh_in_atomlist(keyatom,natommove,atomlist)
     &           ) then
c              --- determine the number of bond-bond terms coming from
c              --- the keyatom
               call twh_inbnbn(GLB_GET,imolty,keyatom,nbnbn)
               do ibnbn = 1,nbnbn
c                 --- set default to not compute
                  lcompute = .false.
c                 --- check the first atom
                  call twh_ijbnbn1(GLB_GET,imolty,keyatom,ibnbn,tbbone)
c                 --- check the second atom
                  call twh_ijbnbn2(GLB_GET,imolty,keyatom,ibnbn,tbbtwo)
                  if ( tbbone .eq. targetatom ) then
c                    --- we always compute if bnbn1 is targetatom
                     lcompute = .true.
                  elseif ( tbbtwo .eq. targetatom ) then
c                    --- compute so long as tbbone is not in atomlist
                     if ( .not. twh_in_atomlist(tbbone,natommove
     &                    ,atomlist) ) then
                        lcompute = .true.
                     endif
                  endif

                  if ( lcompute ) then
c                    --- valid bond-bond to evalulate
                     do ioldnew = CRD_NEW,CRD_OLD,(CRD_OLD-CRD_NEW)
c                       --- keyatom cord is not in the atomlist
                        keyatomcord = globalcord(ioldnew)
                        if ( tbbone .eq. targetatom ) then
c                          --- oneatom is in the atomlist at iatomlist
                           oneatomcord = atomlist(iatomlist,ioldnew)
                        else
c                          --- only valid compute if targetatom is not
c                          --- in tbbone is non-atomlist atom in tbbone
                           oneatomcord = globalcord(ioldnew)
                        endif
                        if ( tbbtwo .eq. targetatom ) then
c                          --- twoatom is in the atomlist at iatomlist
                           twoatomcord = atomlist(iatomlist,ioldnew)
                        else
c                          --- could be either atomlist or globalcord
                           twoatomcord = twh_get_atomlist_cord(tbbtwo
     &                          ,ioldnew,natommove,atomlist,globalcord)
                        endif
                        if ( keyatomcord .ne. CRD_NONE .and.
     &                       oneatomcord .ne. CRD_NONE .and.
     &                       twoatomcord .ne. CRD_NONE ) then
c                          --- fetch keyatom
                           call twh_coordinates(GLB_GET,keyatomcord
     &                          ,ichain,keyatom,xpos(1),ypos(1),zpos(1))
c                          --- fetch oneatom
                           call twh_coordinates(GLB_GET,oneatomcord
     &                          ,ichain,tbbone,xpos(2),ypos(2),zpos(2))
c                          --- fetch twoatom
                           call twh_coordinates(GLB_GET,twoatomcord
     &                          ,ichain,tbbtwo,xpos(3),ypos(3),zpos(3))
c                          --- compute the distance between keyatom
c                          --- and targetatom
                           call twh_mimage(activebox,xpos(1),ypos(1)
     &                          ,zpos(1),xpos(2),ypos(2),zpos(2)
     &                          ,rx,ry,rz)
                           length(1) = twh_distance(rx,ry,rz)
c                          --- compute the distance between 
c                          --- keyatom and jatom
                           call twh_mimage(activebox,xpos(1),ypos(1)
     &                          ,zpos(1),xpos(3),ypos(3),zpos(3)
     &                          ,rx,ry,rz)
                           length(2) = twh_distance(rx,ry,rz)
c                          --- compute energy change
                           dvatomlist = dvatomlist + dvpm(ioldnew)
     &                          *twh_vbonbon(imolty,keyatom,ibnbn
     &                          ,length)
                        endif
                     enddo
                  endif
               enddo
            endif
         enddo
      enddo

      return
      end

      subroutine twh_ec_atomlist_angle(lfailure,ichain,activebox
     &     ,natommove,atomlist,globalcord,loverlap,dvatomlist)
c     ******************************************************************
c     * computes the change in the bending energy for the molecule     *
c     * ichain in the activebox for terms that involve any of the atoms*
c     * in the atomlist                                                *
c     *                                                                *
c     * originally written 10-25-2006 by M.G. Martin                   *
c     * last modified 10-14-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "dvpm.h"
#define FUNCTION_IN_ATOMLIST
#define FUNCTION_GET_ATOMLIST_CORD
#define FUNCTION_VANGLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ichain,activebox,natommove
      integer atomlist
      dimension atomlist(natommove,3)
      integer globalcord
      dimension globalcord(2)
      double precision dvatomlist
      logical loverlap
      dimension loverlap(2)
c     --- local variables
c     --- logical scalars
      logical lcompute
c     --- integer scalars
      integer iatomlist,iatom,jatom,ibend,nbend,testlist,testunit
      integer imolty,ioldnew,targetcord,targetatom,iatomcord,jatomcord
      integer nvib,ivib,styleflag
c     --- double precision scalars
      double precision angcos,distone,disttwo
c     --- double precision arrays
      double precision xpos,ypos,zpos
      dimension xpos(3),ypos(3),zpos(3)

      call twh_moltyp(GLB_GET,ichain,imolty)
c     --- always use the coordinates to compute the angle energy
      styleflag = 1

c     --- on the first pass compute all angle interactions that
c     --- have an end on an atomlist atom
c     --- just need to check to make sure we do not double count when
c     --- both ends have an atomlist atom
      do iatomlist = 1,natommove
         targetatom = atomlist(iatomlist,3)
         call twh_inben(GLB_GET,imolty,targetatom,nbend)
         do ibend = 1,nbend
            call twh_ijben2(GLB_GET,imolty,targetatom,ibend,iatom)
            call twh_ijben3(GLB_GET,imolty,targetatom,ibend,jatom)
c           --- make sure the jatom was not previously a target atom
            lcompute = .true.
            testlist = 0
            do while ( testlist .lt. iatomlist-1 .and. lcompute )
               testlist = testlist + 1
               testunit = atomlist(testlist,3)
               if ( testunit .eq. jatom ) lcompute = .false.
            enddo
            if ( lcompute ) then
               do ioldnew = CRD_NEW,CRD_OLD,(CRD_OLD-CRD_NEW)
c                 --- targetatom is always in the atomlist
                  targetcord = atomlist(iatomlist,ioldnew)
c                 --- determine the coordinates to use for iatom
                  iatomcord = twh_get_atomlist_cord(iatom,ioldnew
     &                 ,natommove,atomlist,globalcord)
c                 --- determine the coordinates to use for jatom
                  jatomcord = twh_get_atomlist_cord(jatom,ioldnew
     &                 ,natommove,atomlist,globalcord)
                  if ( targetcord .ne. CRD_NONE .and.
     &                 iatomcord .ne. CRD_NONE .and.
     &                 jatomcord .ne. CRD_NONE ) then
c                    --- fetch the coordinates for targetatom
                     call twh_coordinates(GLB_GET,targetcord,ichain
     &                    ,targetatom,xpos(1),ypos(1),zpos(1))
c                    --- fetch the coordinates for iatom
                     call twh_coordinates(GLB_GET,iatomcord,ichain,iatom
     &                    ,xpos(2),ypos(2),zpos(2))
c                    --- fetch the coordinates for jatom
                     call twh_coordinates(GLB_GET,jatomcord,ichain,jatom
     &                    ,xpos(3),ypos(3),zpos(3))
c                    --- compute angle energy
                     dvatomlist = dvatomlist + dvpm(ioldnew)
     &                    *twh_vangle(lfailure,loverlap(ioldnew)
     &                    ,activebox,imolty,targetatom,ibend,styleflag
     &                    ,angcos,distone,disttwo,xpos,ypos,zpos)
                     if ( lfailure ) return
                  endif
               enddo
            endif
         enddo
      enddo

c     --- on the second pass we want to find angle terms that have
c     --- the atomlist in the central position
c     --- to avoid double counting we must make sure that none of the
c     --- other atomlist atoms are involved as those would have been
c     --- handled above
      do iatomlist = 1,natommove
         targetatom = atomlist(iatomlist,3)
c        --- check the neighbors of the targetatom
         call twh_invib(GLB_GET,imolty,targetatom,nvib)
         do ivib = 1,nvib
c           --- get the neighboring atom
            call twh_ijvib(GLB_GET,imolty,targetatom,ivib,iatom)
c           --- make sure this is not in the atomlist
            if ( .not. twh_in_atomlist(iatom,natommove,atomlist) ) then
c              --- determine the number of angle terms coming from
c              --- that atom
               call twh_inben(GLB_GET,imolty,iatom,nbend)
               do ibend = 1,nbend
c                 --- see if the central bond-bond is targetatom
                  call twh_ijben2(GLB_GET,imolty,iatom,ibend,testunit)
                  if ( testunit .eq. targetatom ) then
c                    --- find the jatom
                     call twh_ijben3(GLB_GET,imolty,iatom,ibend,jatom)
c                    --- make sure jatom is not in the atomlist
                     if ( .not. twh_in_atomlist(jatom,natommove
     &                    ,atomlist) 
     &                    .and. iatom .lt. jatom ) then
c                       --- valid angle to evalulate
                        do ioldnew = CRD_NEW,CRD_OLD,(CRD_OLD-CRD_NEW)
c                          --- targetatom cord is in the atomlist
                           targetcord = atomlist(iatomlist,ioldnew)
c                          --- iatom is not in the atomlist
                           iatomcord = globalcord(ioldnew)
c                          --- jatom is not in the atomlist
                           jatomcord = globalcord(ioldnew)
c                          --- fetch iatom
                           call twh_coordinates(GLB_GET,iatomcord,ichain
     &                          ,iatom,xpos(1),ypos(1),zpos(1))
c                          --- fetch targetatom
                           call twh_coordinates(GLB_GET,targetcord
     &                          ,ichain,targetatom,xpos(2),ypos(2)
     &                          ,zpos(2))
c                          --- fetch jatom
                           call twh_coordinates(GLB_GET,jatomcord,ichain
     &                          ,jatom,xpos(3),ypos(3),zpos(3))
c                          --- compute energy change
                           dvatomlist = dvatomlist + dvpm(ioldnew)
     &                          *twh_vangle(lfailure,loverlap(ioldnew)
     &                          ,activebox,imolty,iatom,ibend,styleflag
     &                          ,angcos,distone,disttwo,xpos,ypos,zpos)
                           if ( lfailure ) return
                        enddo
                     endif
                  endif
               enddo
            endif
         enddo
      enddo

      return
      end

      subroutine twh_ec_atomlist_angle_angle(ichain,activebox
     &     ,natommove,atomlist,globalcord,dvatomlist)
c     ******************************************************************
c     * computes the change in the angle-angle energy for the          *
c     * molecule ichain in the activebox for terms that involve any of *
c     * the atoms in the atomlist                                      *
c     *                                                                *
c     * originally written 10-25-2006 by M.G. Martin                   *
c     * last modified 08-30-2012 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "dvpm.h"
#define FUNCTION_IN_ATOMLIST
#define FUNCTION_GET_ATOMLIST_CORD
#define FUNCTION_VANGANG
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ichain,activebox,natommove
      integer atomlist
      dimension atomlist(natommove,3)
      integer globalcord
      dimension globalcord(2)
      double precision dvatomlist
c     --- local variables
c     --- logical scalars
      logical lcompute
c     --- integer scalars
      integer iatomlist,iatom,jatom,katom,iaa,naa,testlist,testunit
      integer imolty,ioldnew,targetcord,targetatom,iatomcord,jatomcord
      integer katomcord,nvib,ivib,latom,latomcord,aastyle
c     --- double precision arrays
      double precision xpos,ypos,zpos
      dimension xpos(4),ypos(4),zpos(4)
      double precision angle
      dimension angle(2)

      call twh_moltyp(GLB_GET,ichain,imolty)
c     --- compute from coordinates
      aastyle = AA_COORD

c     --- on the first pass compute all angle-angle interactions that
c     --- start from an atomlist atom
      do iatomlist = 1,natommove
         targetatom = atomlist(iatomlist,3)
         call twh_inaa(GLB_GET,imolty,targetatom,naa)
         do iaa = 1,naa
            call twh_ijaa0(GLB_GET,imolty,targetatom,iaa,iatom)
            call twh_ijaa1(GLB_GET,imolty,targetatom,iaa,jatom)
            call twh_ijaa2(GLB_GET,imolty,targetatom,iaa,katom)
            do ioldnew = CRD_NEW,CRD_OLD,(CRD_OLD-CRD_NEW)
c              --- targetatom is always in the atomlist
               targetcord = atomlist(iatomlist,ioldnew)
c              --- determine the coordinates to use for iatom
               iatomcord = twh_get_atomlist_cord(iatom,ioldnew
     &              ,natommove,atomlist,globalcord)
c              --- determine the coordinates to use for jatom
               jatomcord = twh_get_atomlist_cord(jatom,ioldnew
     &              ,natommove,atomlist,globalcord)
c              --- determine the coordinates to use for katom
               katomcord = twh_get_atomlist_cord(katom,ioldnew
     &              ,natommove,atomlist,globalcord)

               if ( targetcord .ne. CRD_NONE .and.
     &              iatomcord .ne. CRD_NONE .and.
     &              jatomcord .ne. CRD_NONE .and.
     &              katomcord .ne. CRD_NONE ) then
c                 --- fetch the coordinates for targetatom
                  call twh_coordinates(GLB_GET,targetcord,ichain
     &                 ,targetatom,xpos(1),ypos(1),zpos(1))
c                 --- fetch the coordinates for iatom
                  call twh_coordinates(GLB_GET,iatomcord,ichain,iatom
     &                 ,xpos(2),ypos(2),zpos(2))
c                 --- fetch the coordinates for jatom
                  call twh_coordinates(GLB_GET,jatomcord,ichain,jatom
     &                 ,xpos(3),ypos(3),zpos(3))
c                 --- fetch the coordinates for katom
                  call twh_coordinates(GLB_GET,katomcord,ichain,katom
     &                 ,xpos(4),ypos(4),zpos(4))
c                 --- compute angle-angle energy
                  dvatomlist = dvatomlist + dvpm(ioldnew)
     &                 *twh_vangang(aastyle,activebox,imolty,targetatom
     &                 ,iaa,angle,xpos,ypos,zpos)
               endif
            enddo
         enddo
      enddo

c     --- on the second pass we want to find angle terms that have
c     --- at least one atomlist atom in one of the other three positions
c     --- to avoid double counting we must make sure that none of the
c     --- other atomlist atoms are in the first position and that none
c     --- of the previously scanned atom list atoms are in any other
c     --- position
      do iatomlist = 1,natommove
         targetatom = atomlist(iatomlist,3)
c        --- check the neighbors of the targetatom
         call twh_invib(GLB_GET,imolty,targetatom,nvib)
         do ivib = 1,nvib
c           --- get the neighboring atom
            call twh_ijvib(GLB_GET,imolty,targetatom,ivib,iatom)
c           --- make sure this is not in the atomlist
            if ( .not. twh_in_atomlist(iatom,natommove,atomlist) ) then
c              --- determine the number of angle-angle terms coming from
c              --- this atom
               call twh_inaa(GLB_GET,imolty,iatom,naa)
               do iaa = 1,naa
c                 --- get the unit numbers for the other atoms
                  call twh_ijaa0(GLB_GET,imolty,iatom,iaa,jatom)
                  call twh_ijaa1(GLB_GET,imolty,iatom,iaa,katom)
                  call twh_ijaa2(GLB_GET,imolty,iatom,iaa,latom)
c                 --- see if any of the other atoms are the targetatom
                  lcompute = .false.
                  if ( jatom .eq. targetatom .or.
     &                 katom .eq. targetatom .or.
     &                 latom .eq. targetatom ) then
                     lcompute = .true.
                  endif
c                 --- make sure none of the other atoms occur
c                 --- previously in the atomlist
                  testlist = 0
                  do while ( lcompute .and. testlist .lt. iatomlist-1 )
                     testlist = testlist + 1
                     testunit = atomlist(testlist,3)
                     if ( jatom .eq. testlist .or.
     &                    katom .eq. testlist .or.
     &                    latom .eq. testlist ) then
                        lcompute = .false.
                     endif
                  enddo

                  if ( lcompute ) then
c                    --- valid angle-angle, compute energy
                     do ioldnew = CRD_NEW,CRD_OLD,(CRD_OLD-CRD_NEW)
c                       --- iatom is never in the atomlist
                        iatomcord = globalcord(ioldnew)
c                       --- fetch jatomcord
                        jatomcord = twh_get_atomlist_cord(jatom,ioldnew
     &                       ,natommove,atomlist,globalcord)
c                       --- fetch katomcord
                        katomcord = twh_get_atomlist_cord(katom,ioldnew
     &                       ,natommove,atomlist,globalcord)
c                       --- fetch latomcord
                        latomcord = twh_get_atomlist_cord(latom,ioldnew
     &                       ,natommove,atomlist,globalcord)
c                       --- fetch all of the coordinates
                        call twh_coordinates(GLB_GET,iatomcord,ichain
     &                       ,iatom,xpos(1),ypos(1),zpos(1))
                        call twh_coordinates(GLB_GET,jatomcord,ichain
     &                       ,jatom,xpos(2),ypos(2),zpos(2))
                        call twh_coordinates(GLB_GET,katomcord,ichain
     &                       ,katom,xpos(3),ypos(3),zpos(3))
                        call twh_coordinates(GLB_GET,latomcord,ichain
     &                       ,latom,xpos(4),ypos(4),zpos(4))
c                       --- compute angle-angle energy
                        dvatomlist = dvatomlist + dvpm(ioldnew)
     &                       *twh_vangang(aastyle,activebox,imolty
     &                       ,iatom,iaa,angle,xpos,ypos,zpos)
                     enddo
                  endif
               enddo   
            endif
         enddo
      enddo

      return
      end


      subroutine twh_ec_atomlist_torsion(lfailure,ichain,activebox
     &     ,natommove,atomlist,globalcord,loverlap,dvatomlist)
c     ******************************************************************
c     * computes the change in the regular torsion energy for the      *
c     * molecule ichain in the activebox for all regular torsions that *
c     * involve any of the atoms in the atomlist                       *
c     *                                                                *
c     * originally written 10-13-2006 by M.G. Martin                   *
c     * last modified 10-25-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "dvpm.h"
#define FUNCTION_GET_ATOMLIST_CORD
#define FUNCTION_GET_IJTOR2
#define FUNCTION_GET_IJTOR3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_IJVIB
#define FUNCTION_GET_INTOR
#define FUNCTION_GET_INVIB
#define FUNCTION_VTORSION
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ichain,activebox,natommove
      integer atomlist
      dimension atomlist(natommove,3)
      integer globalcord
      dimension globalcord(2)
      double precision dvatomlist
      logical loverlap
      dimension loverlap(2)
c     --- local variables
c     --- logical scalars
      logical lcompute
c     --- integer scalars
      integer icord,itry,iatom,imolty,targetatom,jatom,testatom
      integer icompute,itor,ibead,testunit,testcord
c     --- integer arrays
      integer torbead
      dimension torbead(4)
c     --- double precision scalars
      double precision ddum
c     --- double precision arrays
      double precision xtemp,ytemp,ztemp
      dimension xtemp(4),ytemp(4),ztemp(4)

      call twh_moltyp(GLB_GET,ichain,imolty)
      do iatom = 1,natommove
         targetatom = atomlist(iatom,3)
         do itry = 0,twh_get_invib(imolty,targetatom)
            if (itry .eq. 0) then
c             --- compute torsions starting at this atom
               icompute = targetatom
            else
c             --- compute torsions for neighbors of this atom
               icompute = twh_get_ijvib(imolty,targetatom,itry)
            endif
            do itor = 1,twh_get_intor(imolty,icompute)
c              --- assign the unit types for the torsion
               torbead(1) = icompute
               torbead(2) = twh_get_ijtor2(imolty,icompute,itor)
               torbead(3) = twh_get_ijtor3(imolty,icompute,itor)
               torbead(4) = twh_get_ijtor4(imolty,icompute,itor)
               lcompute = .true.
c              --- only compute if torbead 1 or 2 is targetatom
               if ( torbead(1) .ne. targetatom 
     &              .and. torbead(2) .ne. targetatom ) then
                  lcompute = .false.
               endif
c              --- avoid double counting by making sure iatom
c              --- is the lowest index of any atoms in the atomlist
               do jatom = 1,iatom-1
                  testatom = atomlist(jatom,3)
                  do ibead = 1,4
                     if ( torbead(ibead) .eq. testatom ) then
                        lcompute = .false.
                     endif
                  enddo
               enddo
               if ( lcompute ) then
                  do icord = CRD_NEW,CRD_OLD,(CRD_OLD - CRD_NEW)
c                     --- determine atom positions
                     do ibead = 1,4
c                       --- see if this atom is in the atomlist
                        testunit = torbead(ibead)
                        testcord = twh_get_atomlist_cord(testunit
     &                       ,icord,natommove,atomlist,globalcord)
                        call twh_coordinates(GLB_GET,testcord,ichain
     &                       ,testunit,xtemp(ibead),ytemp(ibead)
     &                       ,ztemp(ibead))
                     enddo
c                    --- compute torsion energy
                     dvatomlist = dvatomlist + dvpm(icord)*
     &                    twh_vtorsion(lfailure,activebox,1,imolty
     &                    ,icompute,itor
     &                    ,xtemp,ytemp,ztemp,ddum,loverlap(icord))   
                     if ( lfailure ) return
                     if ( loverlap(icord) ) return
                  enddo
               endif
            enddo
         enddo
      enddo

      return
      end


      subroutine twh_ec_atomlist_improper(ichain,activebox
     &     ,natommove,atomlist,globalcord,dvatomlist)
c     ******************************************************************
c     * computes the change in the improper torsion energy for the     *
c     * molecule ichain in the activebox for terms that involve any of *
c     * the atoms in the atomlist                                      *
c     *                                                                *
c     * originally written 10-25-2006 by M.G. Martin                   *
c     * last modified 11-08-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "dvpm.h"
#define FUNCTION_IN_ATOMLIST
#define FUNCTION_GET_ATOMLIST_CORD
#define FUNCTION_VIMPROPER
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ichain,activebox,natommove
      integer atomlist
      dimension atomlist(natommove,3)
      integer globalcord
      dimension globalcord(2)
      double precision dvatomlist
c     --- local variables
c     --- logical scalars
      logical lbadstereo,lcompute
c     --- integer scalars
      integer iatomlist,iatom,jatom,katom,itor,ntor,testlist,testunit
      integer imolty,ioldnew,targetcord,targetatom,iatomcord,jatomcord
      integer katomcord,nvib,ivib,latom,latomcord
c     --- double precision scalars
      double precision phic
c     --- double precision arrays
      double precision xpos,ypos,zpos
      dimension xpos(4),ypos(4),zpos(4)

      call twh_moltyp(GLB_GET,ichain,imolty)

c     --- on the first pass compute all improper interactions that
c     --- start from an atomlist atom
      do iatomlist = 1,natommove
         targetatom = atomlist(iatomlist,3)
         call twh_inimprop(GLB_GET,imolty,targetatom,ntor)
         do itor = 1,ntor
            call twh_ijimprop2(GLB_GET,imolty,targetatom,itor,iatom)
            call twh_ijimprop3(GLB_GET,imolty,targetatom,itor,jatom)
            call twh_ijimprop4(GLB_GET,imolty,targetatom,itor,katom)
            do ioldnew = CRD_NEW,CRD_OLD,(CRD_OLD-CRD_NEW)
c              --- targetatom is always in the atomlist
               targetcord = atomlist(iatomlist,ioldnew)
c              --- determine the coordinates to use for iatom
               iatomcord = twh_get_atomlist_cord(iatom,ioldnew
     &              ,natommove,atomlist,globalcord)
c              --- determine the coordinates to use for jatom
               jatomcord = twh_get_atomlist_cord(jatom,ioldnew
     &              ,natommove,atomlist,globalcord)
c              --- determine the coordinates to use for katom
               katomcord = twh_get_atomlist_cord(katom,ioldnew
     &              ,natommove,atomlist,globalcord)

               if ( targetcord .ne. CRD_NONE .and.
     &              iatomcord .ne. CRD_NONE .and.
     &              jatomcord .ne. CRD_NONE .and.
     &              katomcord .ne. CRD_NONE ) then
c                 --- fetch the coordinates for targetatom
                  call twh_coordinates(GLB_GET,targetcord,ichain
     &                 ,targetatom,xpos(1),ypos(1),zpos(1))
c                 --- fetch the coordinates for iatom
                  call twh_coordinates(GLB_GET,iatomcord,ichain,iatom
     &                 ,xpos(2),ypos(2),zpos(2))
c                 --- fetch the coordinates for jatom
                  call twh_coordinates(GLB_GET,jatomcord,ichain,jatom
     &                 ,xpos(3),ypos(3),zpos(3))
c                 --- fetch the coordinates for katom
                  call twh_coordinates(GLB_GET,katomcord,ichain,katom
     &                 ,xpos(4),ypos(4),zpos(4))
c                 --- compute improper energy
                  dvatomlist = dvatomlist + dvpm(ioldnew)
     &                 *twh_vimproper(activebox,imolty,targetatom
     &                 ,itor,xpos,ypos,zpos,phic,lbadstereo)
               endif
            enddo
         enddo
      enddo

c     --- on the second pass we want to find angle terms that have
c     --- at least one atomlist atom in one of the other three positions
c     --- to avoid double counting we must make sure that none of the
c     --- other atomlist atoms are in the first position and that none
c     --- of the previously scanned atom list atoms are in any other
c     --- position
      do iatomlist = 1,natommove
         targetatom = atomlist(iatomlist,3)
c        --- check the neighbors of the targetatom
         call twh_invib(GLB_GET,imolty,targetatom,nvib)
         do ivib = 1,nvib
c           --- get the neighboring atom
            call twh_ijvib(GLB_GET,imolty,targetatom,ivib,iatom)
c           --- make sure this is not in the atomlist
            if ( .not. twh_in_atomlist(iatom,natommove,atomlist) ) then
c              --- determine the number of improper terms coming from
c              --- this atom
               call twh_inimprop(GLB_GET,imolty,iatom,ntor)
               do itor = 1,ntor
c                 --- get the unit numbers for the other atoms
                  call twh_ijimprop2(GLB_GET,imolty,iatom,itor,jatom)
                  call twh_ijimprop3(GLB_GET,imolty,iatom,itor,katom)
                  call twh_ijimprop4(GLB_GET,imolty,iatom,itor,latom)
c                 --- see if any of the other atoms are the targetatom
                  lcompute = .false.
                  if ( jatom .eq. targetatom .or.
     &                 katom .eq. targetatom .or.
     &                 latom .eq. targetatom ) then
                     lcompute = .true.
                  endif
c                 --- make sure none of the other atoms occur
c                 --- previously in the atomlist
                  testlist = 0
                  do while ( lcompute .and. testlist .lt. iatomlist-1 )
                     testlist = testlist + 1
                     testunit = atomlist(testlist,3)
                     if ( jatom .eq. testlist .or.
     &                    katom .eq. testlist .or.
     &                    latom .eq. testlist ) then
                        lcompute = .false.
                     endif
                  enddo
                  if ( lcompute ) then
c                    --- valid improper, compute improper energy
                     do ioldnew = CRD_NEW,CRD_OLD,(CRD_OLD-CRD_NEW)
c                       --- iatom is never in the atomlist
                        iatomcord = globalcord(ioldnew)
c                       --- fetch jatomcord
                        jatomcord = twh_get_atomlist_cord(jatom,ioldnew
     &                       ,natommove,atomlist,globalcord)
c                       --- fetch katomcord
                        katomcord = twh_get_atomlist_cord(katom,ioldnew
     &                       ,natommove,atomlist,globalcord)
c                       --- fetch latomcord
                        latomcord = twh_get_atomlist_cord(latom,ioldnew
     &                       ,natommove,atomlist,globalcord)
c                       --- fetch all of the coordinates
                        call twh_coordinates(GLB_GET,iatomcord,ichain
     &                       ,iatom,xpos(1),ypos(1),zpos(1))
                        call twh_coordinates(GLB_GET,jatomcord,ichain
     &                       ,jatom,xpos(2),ypos(2),zpos(2))
                        call twh_coordinates(GLB_GET,katomcord,ichain
     &                       ,katom,xpos(3),ypos(3),zpos(3))
                        call twh_coordinates(GLB_GET,latomcord,ichain
     &                       ,latom,xpos(4),ypos(4),zpos(4))
c                       --- compute improper energy
                        dvatomlist = dvatomlist + dvpm(ioldnew)
     &                       *twh_vimproper(activebox,imolty,iatom
     &                       ,itor,xpos,ypos,zpos,phic,lbadstereo)
                     enddo
                  endif
               enddo   
            endif
         enddo
      enddo

      return
      end

      subroutine twh_ec_atomlist_nonbond(lfailure,ichain,activebox
     &     ,natommove,atomlist,globalcord,loverlap,dvquantum
     &     ,dvclassical,dvinter,dvintra,dvelect,dvsolv,dvext
     &     ,dvforeign,ddvdl_vdw,ddvdl_c)
c     ******************************************************************
c     * computes the nonbonded (quantum or classical plus coulomb)     *
c     * energy change for a list of atoms on a single molecule ichain  *
c     * in box activebox                                               *
c     *                                                                *
c     * originally written 10-20-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ichain,activebox,natommove
      integer atomlist
      dimension atomlist(natommove,3)
      integer globalcord
      dimension globalcord(2)
      logical loverlap
      dimension loverlap(2)
      double precision dvquantum,dvclassical,dvinter,dvintra,dvelect
      double precision dvsolv,dvext
      double precision dvforeign
      dimension dvforeign(MAX_FOREIGN_LAMBDA)
      double precision ddvdl_vdw,ddvdl_c
c     --- local variables
      integer lambda,scp_num_foreign_lambda,potentialstyle

c     --- retrive constant
      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
      call twh_potentialstyle(GLB_GET,potentialstyle)

c     --- initialize energy changes
      dvquantum = 0.0d0
      dvclassical = 0.0d0
      dvinter = 0.0d0
      dvintra = 0.0d0
      dvelect = 0.0d0
      dvsolv = 0.0d0
      dvext = 0.0d0
      ddvdl_vdw = 0.0d0
      ddvdl_c = 0.0d0
      do lambda = 1,scp_num_foreign_lambda
         dvforeign(lambda) = 0.0d0
      enddo

      if ( potentialstyle .eq. POT_EXTERNAL ) then
c        --- not yet implemented
         write(6,*) 'EC_ATOMLIST_NONBOND: not yet implemented for'
     &        ,' quantum potentialstyles'
         lfailure = .true.
         return
      endif

      if ( potentialstyle .eq. POT_INTERNAL ) then
c        --- get the one-body energy change
         call twh_ec_atomlist_one_body(lfailure,ichain,activebox
     &        ,natommove,atomlist,globalcord,loverlap,dvsolv,dvext)
         if ( lfailure .or. loverlap(CRD_NEW) .or. loverlap(CRD_OLD)
     &        ) return
c        --- get the two-body energy change
         call twh_ec_atomlist_two_body(lfailure,ichain,activebox
     &        ,natommove,atomlist,globalcord,loverlap,dvinter,dvelect
     &        ,dvforeign,ddvdl_vdw,ddvdl_c)
         if ( lfailure .or. loverlap(CRD_NEW) .or. loverlap(CRD_OLD)
     &        ) return
c        --- get the multi-body energy change
         call twh_ec_atomlist_multi_body(lfailure,ichain,activebox
     &        ,natommove,atomlist,dvinter)
         if ( lfailure ) return
c        --- get the nonbonded intramolecular energy change
         call twh_ec_atomlist_nb_intra(lfailure,ichain,activebox
     &        ,natommove,atomlist,globalcord,loverlap,dvintra,dvelect)
         if ( lfailure .or. loverlap(CRD_NEW) .or. loverlap(CRD_OLD)
     &        ) return
c        --- get the long-range electrostatic energy change
         call twh_ec_atomlist_lr_elect(ichain,activebox
     &        ,natommove,atomlist,dvelect)

         dvclassical = dvinter + dvintra + dvelect + dvsolv + dvext

      endif

      return
      end

      subroutine twh_ec_atomlist_nb_intra(lfailure,ichain,activebox
     &     ,natommove,atomlist,globalcord,loverlap,dvintra,dvelect)
c     ******************************************************************
c     * computes the change in the nonbonded intramolecular energies   *
c     * for the molecule ichain in box activebox.                      *
c     * this is typically broken into a nonbonded vdw intramolecular   *
c     * energy (dvintra) and an intramolecular and self-term electro   *
c     * energy (dvelect)                                               *
c     *                                                                *
c     * originally written 10-23-2006 by M.G. Martin                   *
c     * last modified 10-14-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "dvpm.h"
#define FUNCTION_EWALD_CORRECT
#define FUNCTION_EWALD_SELF
#define FUNCTION_LINCLUDE
#define FUNCTION_ONEFIVETYPE
#define FUNCTION_VCOULOMB
#define FUNCTION_VONEFIVE
#define FUNCTION_VTWOBODY
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ichain,activebox,natommove
      integer atomlist
      dimension atomlist(natommove,3)
      integer globalcord
      dimension globalcord(2)
      logical loverlap
      dimension loverlap(2)
      double precision dvintra,dvelect
c     --- local variables
c     --- logical scalars
      logical ltorsion,lcompute,lvalidintra,luse_ct,luse_ci
      logical luse_vt,luse_vi,lcoulomb,lewald
c     --- integer scalars
      integer testunit,iatomlist,icord,ioldnew,ntii,ntjj,imolty
      integer iunit,gcord,cutstyle,ioftype,nunit
c     --- double precision scalars
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,dvdlv,dvdlc
      double precision rxuij,ryuij,rzuij,rijsq,rij,qqfact,rcelect
c     --- double precision arrays
      double precision vf
      dimension vf(MAX_FOREIGN_LAMBDA)

c     --- set defaults
      ltorsion = .false.
      cutstyle = CUTSTYLE_FULL

c     --- retrieve constants
      call twh_constant_qqfact(GLB_GET,qqfact)
      call twh_moltyp(GLB_GET,ichain,imolty)
      call twh_nunit(GLB_GET,imolty,nunit)
      call twh_rcelect(GLB_GET,activebox,rcelect)
      call twh_lcoulomb(GLB_GET,lcoulomb)
      call twh_lewald(GLB_GET,lewald)

      do testunit = 1,nunit
c        --- make sure this unit is not in the atomlist
         lcompute = .true.
         do iatomlist = 1,natommove
            iunit = atomlist(iatomlist,3)
            if ( testunit .eq. iunit ) lcompute = .false.
         enddo
         if ( lcompute ) then
            call twh_luse_c(GLB_GET,imolty,testunit,luse_ct)
            call twh_luse_v(GLB_GET,imolty,testunit,luse_vt)
            call twh_ntype(GLB_GET,imolty,testunit,ntjj)
            do iatomlist = 1,natommove
               iunit = atomlist(iatomlist,3)
               call twh_luse_c(GLB_GET,imolty,iunit,luse_ci)
               call twh_luse_v(GLB_GET,imolty,iunit,luse_vi)
               lvalidintra = twh_linclude(imolty,testunit,iunit)
               if ( lvalidintra .or. lewald ) then
c                 --- need to compute distance either for use in the
c                 --- direct intramolecular interactions or ewald terms
                  ioftype = twh_onefivetype(imolty,testunit,iunit)
                  call twh_ntype(GLB_GET,imolty,iunit,ntii)
                  do ioldnew = CRD_NEW,CRD_OLD,(CRD_OLD-CRD_NEW)
                     icord = atomlist(iatomlist,ioldnew)
                     if ( icord .ne. CRD_NONE ) then
                        call twh_coordinates(GLB_GET,icord,ichain,iunit
     &                       ,rxui,ryui,rzui)
                        gcord = globalcord(ioldnew)
                        call twh_coordinates(GLB_GET,gcord,ichain
     &                       ,testunit,rxuj,ryuj,rzuj)
                        call twh_mimage(activebox,rxui,ryui,rzui,rxuj
     &                       ,ryuj,rzuj,rxuij,ryuij,rzuij)
                        rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                        if ( lvalidintra ) then
c                          --- compute nonbonded intramolecular energy
                           if ( luse_vt .and. luse_vi ) then
                              if ( ioftype .eq. 0 ) then
                                 dvintra = dvintra 
     &                                + dvpm(ioldnew)
     &                                *twh_vtwobody(lfailure,rijsq,ntii
     &                                ,ntjj,cutstyle,loverlap(ioldnew)
     &                                ,vf,dvdlv)
                                 if ( lfailure ) return
                              else
                                 dvintra = dvintra + dvpm(ioldnew)
     &                                *twh_vonefive(lfailure
     &                                ,rijsq,ioftype,cutstyle
     &                                ,loverlap(ioldnew))
                                 if ( lfailure ) return
                              endif
                              if (loverlap(ioldnew)) return
                           endif
c                          --- compute nonbonded electrostatic energy
                           if ( lcoulomb ) then
                              if ( luse_ct .and. luse_ci ) then
                                 rij = dsqrt(rijsq)
                                 dvelect = dvelect
     &                                + dvpm(ioldnew)*qqfact
     &                                *twh_vcoulomb(activebox,rij,imolty
     &                                ,iunit,imolty,testunit,ltorsion,vf
     &                                ,dvdlc,loverlap(ioldnew))
                                 if (loverlap(ioldnew)) return
                              endif
                           endif
                        elseif ( lewald .and. activebox .ne. 0 ) then
c                          --- compute the ewald intramolecular
c                          --- (self and correction) terms
                           rij = dsqrt(rijsq)
                           if ( rij .lt. rcelect ) then
                              dvelect = dvelect + dvpm(ioldnew)
     &                             *qqfact*twh_ewald_correct(activebox
     &                             ,imolty,iunit,testunit,rij)
                           endif
                        endif
                     endif
                  enddo
               endif
            enddo
         endif
      enddo

c     --- atomlist self interactions
      if ( lewald .and. activebox .ne. 0 ) then
         do iatomlist = 1,natommove
            iunit = atomlist(iatomlist,3)
            do ioldnew = CRD_NEW,CRD_OLD,(CRD_OLD-CRD_NEW)
               icord = atomlist(iatomlist,ioldnew)
               if ( icord .ne. CRD_NONE ) then
                  dvelect = dvelect + dvpm(ioldnew)*qqfact
     &                 *twh_ewald_self(activebox,imolty,iunit)
               endif
            enddo
         enddo
      endif

      return
      end

      subroutine twh_ec_atomlist_one_body(lfailure,ichain,activebox
     &     ,natommove,atomlist,globalcord,loverlap,dvsolv,dvext)
c     ******************************************************************
c     * computes the change in the one-body energy for a list of atoms *
c     * on molecule ichain in box activebox                            *
c     *                                                                *
c     * originally written 10-20-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "dvpm.h"
#define FUNCTION_VFIELD
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ichain,activebox,natommove
      integer atomlist
      dimension atomlist(natommove,3)
      integer globalcord
      dimension globalcord(2)
      logical loverlap
      dimension loverlap(2)
      double precision dvsolv,dvext
c     --- local variables
c     --- logical scalars
      logical ltransfer
c     --- integer scalars
      integer iunit,ntii,imolty,testcord,ioldnew,testkey,testunit
      integer iatommove,cordkey,nchainmove,isolvtype,nunit,nfield
c     --- integer arrays
      integer chaininfo
      dimension chaininfo(1,3)
c     --- double precision scalars
      double precision vsolv
c     --- double precision arrays
      double precision cord
      dimension cord(3)

c     --- interaction energy with external field
      call twh_nfield(GLB_GET,nfield)
      if ( nfield .ne. 0 ) then
         call twh_moltyp(GLB_GET,ichain,imolty)
         do iatommove = 1,natommove
            iunit = atomlist(iatommove,3)
            call twh_ntype(GLB_GET,imolty,iunit,ntii)
            do ioldnew = CRD_NEW,CRD_OLD,(CRD_OLD - CRD_NEW)
               testcord = atomlist(iatommove,ioldnew)
               if ( testcord .ne. CRD_NONE ) then
                  call twh_coordinates(GLB_GET,testcord,ichain,iunit
     &                 ,cord(1),cord(2),cord(3))
                  dvext = dvext + dvpm(ioldnew)*twh_vfield(lfailure
     &                 ,loverlap(ioldnew),activebox,ntii,imolty,ichain
     &                 ,iunit,cord)
                  if ( lfailure .or. loverlap(ioldnew) ) return
               endif
            enddo
         enddo
      endif

c     --- Calculate the solvation energy
c     --- this is currently implemented in a manner that recomputes the
c     --- solvation energy for the entire molecule instead of just the
c     --- atoms that moved.  This is required due to the form of some of
c     --- the solvation potentials.
      call twh_isolvtype(GLB_GET,isolvtype)
      if ( isolvtype .ne. SOLV_NONE ) then
         call twh_moltyp(GLB_GET,ichain,imolty)
         call twh_nunit(GLB_GET,imolty,nunit)
         do ioldnew = CRD_OLD,CRD_NEW
c           --- need to store the entire coordinates of this chain
c           --- this is only going to work if all of the atoms have the
c           --- same coordinate key
            cordkey = atomlist(1,ioldnew)
            do iatommove = 1,natommove
               testkey = atomlist(iatommove,ioldnew)
               if ( testkey .ne. cordkey ) then
                  write(6,*) 'EC_ATOMLIST_ONEBODY: problem with cordkey'
                  write(6,*) 'testkey:',testkey
                  write(6,*) 'cordkey:',cordkey
                  lfailure = .true.
                  return
               endif
            enddo
c           --- see if the globalcord is the same as the cordkey
            if ( globalcord(ioldnew) .ne. cordkey ) then
c              --- this is a problem that can be fixed if cordkey is
c              --- CRD_TEMP and globalcord is CRD_REAL by copying the
c              --- CRD_REAL into the CRD_TEMP for the atoms that are 
c              --- not in the atommove list
               if ( globalcord(ioldnew) .eq. CRD_REAL
     &              .and. cordkey .eq. CRD_TEMP ) then
                  do iunit = 1,nunit
c                    --- make sure this is not in the atomlist
                     ltransfer = .true.
                     do iatommove = 1,natommove
                        testunit = atomlist(iatommove,3)
                        if ( testunit .eq. iunit ) ltransfer = .false.
                     enddo
                     if ( ltransfer ) then
                        call twh_coordinates(GLB_GET,CRD_REAL,ichain
     &                       ,iunit,cord(1),cord(2),cord(3))
                        call twh_coordinates(GLB_SET,CRD_TEMP,ichain
     &                       ,iunit,cord(1),cord(2),cord(3))
                     endif
                  enddo
               else
c                 --- this is a problem
                  write(6,*) 'EC_ATOMLIST_ONE_BODY: incompatible'
     &                 ,' globalcord and cordkey'
                  write(6,*) globalcord(ioldnew)
                  write(6,*) cordkey
                  lfailure = .true.
                  return
               endif
            endif
c           --- compute solvation energy for this molecule
            nchainmove = 1
            chaininfo(1,ioldnew) = cordkey
            chaininfo(1,3) = ichain
            call twh_eng_solvation(lfailure,ioldnew,activebox,nchainmove
     &           ,chaininfo,vsolv)
            if ( lfailure ) return
            dvsolv = dvsolv + dvpm(ioldnew)*vsolv
         enddo
      endif

      return
      end

      subroutine twh_ec_atomlist_two_body(lfailure,ichain,activebox
     &     ,natommove,atomlist,globalcord,loverlap,dvinter,dvelect
     &     ,dvforeign,dvdl_vdw,dvdl_c)
c     ******************************************************************
c     * computes the two-body energy change for a list of atoms on     *
c     * molecule ichain in the activebox                               *
c     *                                                                *
c     * originally written 10-20-2006 by M.G. Martin                   *
c     * last modified 10-14-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_VCOULOMB
#define FUNCTION_VTWOBODY
#include "functions.h"
#include "dvpm.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ichain,activebox,natommove
      integer atomlist
      dimension atomlist(natommove,3)
      integer globalcord
      dimension globalcord(2)
      logical loverlap
      dimension loverlap(2)
      double precision dvinter,dvelect
      double precision dvforeign
      dimension dvforeign(MAX_FOREIGN_LAMBDA)
      double precision dvdl_vdw,dvdl_c
c     --- logical scalars
      logical lsame,lfetch,leval_vdw,leval_coul,leval_v_j,leval_c_j
      logical ltorsion,leval_c_i,leval_v_i
c     --- integer scalars
      integer ioldnew,testmolty,testlist,testchain,iunit,numunit
      integer testtype,testcord,ilist,icord,iatom,cutstyle,imolty,itype
      integer lambda,nmolty,ncmt,scp_num_foreign_lambda
c     --- double precision scalars
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij
      double precision rijsq,dl_vdw,dl_c,rij,qqfact
c     --- double precision arrays
      double precision vfor
      dimension vfor(MAX_FOREIGN_LAMBDA)

c     --- retrieve constants
      call twh_constant_qqfact(GLB_GET,qqfact)
      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)

c     --- use the full potential
      cutstyle = CUTSTYLE_FULL
c     --- not a torsion interaction
      ltorsion = .false.
c     --- set the molecule type for ichain
      call twh_moltyp(GLB_GET,ichain,imolty)
c     --- retrieve nmolty
      call twh_nmolty(GLB_GET,nmolty)
c     --- see if the globalcords are the same
      if ( globalcord(CRD_NEW) .eq. globalcord(CRD_OLD) ) then
         lsame = .true.
      else
         lsame = .false.
      endif

c     --- run through all atoms on other molecules
      do testmolty = 1,nmolty
         call twh_nunit(GLB_GET,testmolty,numunit)
         call twh_ncmt(GLB_GET,activebox,testmolty,ncmt)
         do testlist = 1,ncmt
            call twh_chainlist(GLB_GET,testlist,activebox,testmolty
     &           ,testchain)
            if ( testchain .ne. ichain ) then
               do iunit = 1,numunit
                  call twh_ntype(GLB_GET,testmolty,iunit,testtype)
                  call twh_luse_v(GLB_GET,testmolty,iunit,leval_v_j)
                  call twh_luse_c(GLB_GET,testmolty,iunit,leval_c_j)
                  do ioldnew = CRD_NEW,CRD_OLD,(CRD_OLD-CRD_NEW)
c                    --- fetch the coordinates for the iunit
                     if ( ioldnew .eq. CRD_NEW ) then
                        testcord = globalcord(CRD_NEW)
                        lfetch = .true.
                     else
                        if ( lsame ) then
c                          --- we already fetched these coordinates
                           lfetch = .false.
                        else
                           testcord = globalcord(CRD_OLD)
                           lfetch = .true.
                        endif
                     endif
                     if ( lfetch ) then
                        call twh_coordinates(GLB_GET,testcord
     &                       ,testchain,iunit,rxuj,ryuj,rzuj)
                     endif
c                    --- run through all of the atoms in atomlist
                     do ilist = 1,natommove
                        icord = atomlist(ilist,ioldnew)
                        if ( icord .ne. CRD_NONE ) then
c                          --- set atom type
                           iatom = atomlist(ilist,3)
                           call twh_luse_v(GLB_GET,imolty,iatom
     &                          ,leval_v_i)
                           leval_vdw = leval_v_j .and. leval_v_i
                           call twh_luse_c(GLB_GET,imolty,iatom
     &                          ,leval_c_i)
                           leval_coul = leval_c_j .and. leval_c_i
                           if ( leval_vdw .or. leval_coul ) then
                              call twh_ntype(GLB_GET,imolty,iatom,itype)
c                             --- fetch coordinates for list
                              call twh_coordinates(GLB_GET,icord,ichain
     &                             ,iatom,rxui,ryui,rzui)
c                             --- minimum image the vector
                              call twh_mimage(activebox,rxui,ryui,rzui
     &                             ,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij)
c                             --- compute distance squared for vdw
                              rijsq = rxuij*rxuij + ryuij*ryuij
     &                             + rzuij*rzuij
                              if ( leval_vdw ) then
c                                --- compute change in vdw energy
                                 dvinter = dvinter + dvpm(ioldnew)*
     &                                twh_vtwobody(lfailure,rijsq
     &                                ,itype,testtype,cutstyle
     &                                ,loverlap(ioldnew),vfor,dl_vdw)
                                 if ( lfailure ) return
                                 if ( scp_num_foreign_lambda .ne. 0 
     &                                ) then
                                    dvdl_vdw = dvdl_vdw
     &                                   + dvpm(ioldnew)*dl_vdw
                                    do lambda = 1,scp_num_foreign_lambda
                                       dvforeign(lambda) 
     &                                      = dvforeign(lambda)
     &                                      + dvpm(ioldnew)*vfor(lambda)
                                    enddo
                                 endif
                              endif
                              if ( leval_coul ) then
c                                --- compute change in real coulomb
                                 rij = dsqrt(rijsq)
                                 ltorsion = .false.
                                 dvelect = dvelect 
     &                                + dvpm(ioldnew)*qqfact*
     &                                twh_vcoulomb(activebox,rij
     &                                ,imolty,iatom,testmolty,iunit
     &                                ,ltorsion,vfor,dl_c
     &                                ,loverlap(ioldnew))
                                 if ( scp_num_foreign_lambda .ne. 0
     &                                ) then
                                    dvdl_c = dvdl_c
     &                                   + dvpm(ioldnew)*dl_c
                                    do lambda = 1,scp_num_foreign_lambda
                                       dvforeign(lambda) =
     &                                      dvforeign(lambda)
     &                                      + dvpm(ioldnew)*vfor(lambda)
                                    enddo
                                 endif
                              endif
                           endif
                        endif
                     enddo
                  enddo
               enddo
            endif
         enddo
      enddo

      return
      end

      subroutine twh_ec_atomlist_multi_body(lfailure,ichain,activebox
     &     ,natommove,atomlist,dvinter)
c     ******************************************************************
c     * computes the multi-body energy change in the activebox for the *
c     * atoms in ichain                                                *
c     * currently this is not very different from the whole chain      *
c     * version as the multibody potentials currently implmented all   *
c     * are monatomic                                                  *
c     *                                                                *
c     * originally written 10-20-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ichain,activebox,natommove
      integer atomlist
      dimension atomlist(natommove,3)
      double precision dvinter
c     --- local variables
c     --- integer scalars
      integer nchainmove,classical_pot_num
c     --- integer arrays
      integer chaininfo
      dimension chaininfo(1,3)

      call twh_cp_number(GLB_GET,classical_pot_num)
      if ( classical_pot_num .eq. CP_STILLINGER_WEBER ) then
         if ( natommove .ne. 1 .or. atomlist(1,3) .ne. 1 ) then
            write(6,*) 'EC_ATOMLIST_MULTI_BODY: invalid potential'
     &           ,' and molecule geometry combination'
            lfailure = .true.
            return
         endif
         nchainmove = 1
         chaininfo(1,CRD_NEW) = atomlist(1,CRD_NEW)
         chaininfo(1,CRD_OLD) = atomlist(1,CRD_OLD)
         chaininfo(1,3) = ichain
         call twh_ec_stillinger_webber(lfailure,activebox,nchainmove
     &        ,chaininfo,dvinter)
         if ( lfailure ) return
      elseif ( classical_pot_num .eq. CP_EMBEDDED_ATOM_METHOD ) then
         if ( natommove .ne. 1 .or. atomlist(1,3) .ne. 1 ) then
            write(6,*) 'EC_ATOMLIST_MULTI_BODY: invalid potential'
     &           ,' and molecule geometry combination'
            lfailure = .true.
            return
         endif
         nchainmove = 1
         chaininfo(1,CRD_NEW) = atomlist(1,CRD_NEW)
         chaininfo(1,CRD_OLD) = atomlist(1,CRD_OLD)
         chaininfo(1,3) = ichain
         call twh_ec_embedded_atom(lfailure,activebox,nchainmove
     &        ,chaininfo,dvinter)
         if ( lfailure ) return
      endif

      return
      end

      subroutine twh_ec_atomlist_lr_elect(ichain,activebox,natommove
     &     ,atomlist,dvelect)
c     ******************************************************************
c     * long range electrostatic energy change for the atoms listed in *
c     * atomlist on molecule ichain in box activebox                   *
c     *                                                                *
c     * originally written 10-24-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer ichain,activebox,natommove
      integer atomlist
      dimension atomlist(natommove,3)
      double precision dvelect
c     --- local variables
c     --- logical scalars
      logical lewald
c     --- integer scalars
      integer ic,imolty,ncount,icrd,iunit,testcrd,iatomlist
c     --- double precision scalars
      double precision arg,kxfetch,kyfetch,kzfetch,xtemp,ytemp,ztemp
      double precision charge,ssumr,ssumrn,ssumi,ssumin
c     --- double precision arrays
      double precision sumr,sumi,vrecip
      dimension sumr(2),sumi(2),vrecip(2)

c     --- retrieve constants
      call twh_lewald(GLB_GET,lewald)

      if ( lewald ) then
c        --- bail out if this is the ideal chain box
         if ( activebox .eq. 0 ) return
c        --- initialze vrecip
         do icrd = CRD_OLD,CRD_NEW,( CRD_NEW - CRD_OLD )
            vrecip(icrd) = 0.0d0
         enddo

         call twh_ewald_numvect(GLB_GET,activebox,ncount)
c        --- get molecule type
         call twh_moltyp(GLB_GET,ichain,imolty)
         do ic = 1, ncount
c           --- fetch the kx, ky, and kz
            call twh_ewald_kvector(GLB_GET,CRD_REAL,ic,activebox
     &           ,kxfetch,kyfetch,kzfetch)
            do icrd = CRD_OLD,CRD_NEW,( CRD_NEW - CRD_OLD )
               sumr(icrd) = 0.0d0
               sumi(icrd) = 0.0d0
               do iatomlist = 1,natommove
                  testcrd = atomlist(iatomlist,icrd)
                  if ( testcrd .ne. CRD_NONE ) then
                     iunit = atomlist(iatomlist,3)
                     call twh_coordinates(GLB_GET,testcrd,ichain
     &                    ,iunit,xtemp,ytemp,ztemp)
                     arg = kxfetch*xtemp
     &                    + kyfetch*ytemp
     &                    + kzfetch*ztemp
                     call twh_qqatom(GLB_GET,imolty,iunit,charge)
                     sumr(icrd) = sumr(icrd) + charge*dcos(arg)
                     sumi(icrd) = sumi(icrd) + charge*dsin(arg)
                  endif
               enddo
            enddo
            call twh_ewald_ssumr(GLB_GET,CRD_REAL,ic,activebox,ssumr)
            ssumrn = ssumr - sumr(CRD_OLD) + sumr(CRD_NEW)
            call twh_ewald_ssumr(GLB_SET,CRD_NEW,ic,activebox,ssumrn)
            call twh_ewald_ssumi(GLB_GET,CRD_REAL,ic,activebox,ssumi)
            ssumin = ssumi - sumi(CRD_OLD) + sumi(CRD_NEW)
            call twh_ewald_ssumi(GLB_SET,CRD_NEW,ic,activebox,ssumin)
         enddo

         call twh_ewald_sum(activebox,vrecip)

         dvelect = dvelect + vrecip(CRD_NEW) - vrecip(CRD_OLD)
      endif

      return
      end

      function twh_get_atomlist_cord(testunit,oldnewflag,natommove
     &     ,atomlist,globalcord)
c     ******************************************************************
c     * utility function for the often performed determination of the  *
c     * coordinates to use for a particular atom in combination with   *
c     * the atomlist energy change subroutines                         *
c     *                                                                *
c     * originally written 10-25-2006 by M.G. Martin                   *
c     * last modified 10-25-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_ATOMLIST_CORD
#include "functions.h"
c     --- variables passed to/from the function
      integer testunit,oldnewflag,natommove
      integer atomlist
      dimension atomlist(natommove,3)
      integer globalcord
      dimension globalcord(2)
c     --- local variables
      integer jatom,testatom

      twh_get_atomlist_cord = globalcord(oldnewflag)
      do jatom = 1,natommove
         testatom = atomlist(jatom,3)
         if ( testatom .eq. testunit ) then
            twh_get_atomlist_cord = atomlist(jatom,oldnewflag)
         endif
      enddo

      return
      end

      function twh_in_atomlist(iatom,natommove,atomlist)
c     ******************************************************************
c     * utility function that returns a logical to indicate whether the*
c     * iatom is in the atomlist                                       *
c     *                                                                *
c     * originally written 10-25-2006 by M.G. Martin                   *
c     * last modified 10-25-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_IN_ATOMLIST
#include "functions.h"
c     --- variables passed to/from the function
      integer iatom,natommove
      integer atomlist
      dimension atomlist(natommove,3)
c     --- local variables
      integer iatomlist,testunit

      twh_in_atomlist = .false.
      do iatomlist = 1,natommove
         testunit = atomlist(iatomlist,3)
         if ( testunit .eq. iatom ) then
            twh_in_atomlist = .true.
            return
         endif
      enddo
      return
      end
