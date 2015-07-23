#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2014 Marcus G. Martin                       *
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
      subroutine twh_initialize(lfailure)
c     ******************************************************************
c     * initializes all of the force field parameters                  *
c     * initializes many of the input parameters                       *
c     *                                                                *
c     * originally written 03-01-2002 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
      integer ivalue,ij,ibox,idim,jdim
      integer inix,iniy,iniz
      character*(FFNAMELEN) ctemp
      double precision dvalue,twopi

c     --- default values for random number generation
      call twh_random_code(GLB_SET,RNG_DX_1597_2_7)
      call twh_random_luxlevel(GLB_SET,3)
      call twh_random_seed(GLB_SET,1302002)
      call twh_random_allow_restart(GLB_SET,.true.)
c     --- set the names of the stress tensor components
      call twh_pcomp(GLB_SET,0,'Pressure')
      call twh_pcomp(GLB_SET,1,'S_xx')
      call twh_pcomp(GLB_SET,2,'S_yy')
      call twh_pcomp(GLB_SET,3,'S_zz')
      call twh_pcomp(GLB_SET,4,'S_xy')
      call twh_pcomp(GLB_SET,5,'S_xz')
      call twh_pcomp(GLB_SET,6,'S_yz')
      call twh_pcomp(GLB_SET,7,'P_tail')
c     --- initialite nener
      call twh_nener(GLB_SET,0)
c     --- set up constants and conversion factors
      call twh_constant_pi(GLB_INIT,dvalue)
      call twh_constant_twopi(GLB_INIT,dvalue)
      call twh_constant_halfpi(GLB_INIT,dvalue)
      call twh_constant_invrootpi(GLB_INIT,dvalue)
      call twh_constant_onethird(GLB_INIT,dvalue)
      call twh_constant_underflow(GLB_INIT,dvalue)
      call twh_constant_radtodeg(GLB_INIT,dvalue)
      call twh_constant_cuberoottwo(GLB_INIT,dvalue)
c     --- safety check maxtabtype
      if ( MAXTABTYPE .gt. NNTYPE ) then
         write(6,*) 'MAXTABTYPE must be <= NNTYPE'
         lfailure = .true.
         return
      endif
c     --- initialize io_directory='' to signify writing to the current
c     --- directory.  this variable is modified in mainloop if we are
c     --- writing to a different directory
      call twh_io_directory(GLB_SET,'')

c     --- initialize input arrays
      call twh_invib(GLB_INIT,NTMAX,NUMAX,0)
      call twh_itvib(GLB_INIT,NTMAX,NUMAX,NNBOND,0)
      call twh_ijvib(GLB_INIT,NTMAX,NUMAX,NNBOND,0)
      ctemp = 'null'
      call twh_vibration_order(GLB_INIT,NUMAX,NNBOND,ctemp)
      call twh_ijben2(GLB_INIT,NTMAX,NUMAX,MAXBEND,0)
      call twh_ijben3(GLB_INIT,NTMAX,NUMAX,MAXBEND,0)
      call twh_inben(GLB_INIT,NTMAX,NUMAX,0)
      call twh_itben(GLB_INIT,NTMAX,NUMAX,MAXBEND,0)
      call twh_orderben(GLB_INIT,NTMAX,NUMAX,MAXBEND,ORD_NULL)
      call twh_inimprop(GLB_INIT,NTMAX,NUMAX,0)
      call twh_ijimprop2(GLB_INIT,NTMAX,NUMAX,MAXIMPROP,0)
      call twh_ijimprop3(GLB_INIT,NTMAX,NUMAX,MAXIMPROP,0)
      call twh_ijimprop4(GLB_INIT,NTMAX,NUMAX,MAXIMPROP,0)
      call twh_itimprop(GLB_INIT,NTMAX,NUMAX,MAXIMPROP,0)
      call twh_intor(GLB_INIT,NTMAX,NUMAX,0)
      call twh_ijtor2(GLB_INIT,NTMAX,NUMAX,MAXTOR,0)
      call twh_ijtor3(GLB_INIT,NTMAX,NUMAX,MAXTOR,0)
      call twh_ijtor4(GLB_INIT,NTMAX,NUMAX,MAXTOR,0)
      call twh_ittor(GLB_INIT,NTMAX,NUMAX,MAXTOR,0)
      call twh_inaa(GLB_INIT,NTMAX,NUMAX,0)
      call twh_itaa(GLB_INIT,NTMAX,NUMAX,MAXAA,0)
      call twh_ijaa0(GLB_INIT,NTMAX,NUMAX,MAXAA,0)
      call twh_ijaa1(GLB_INIT,NTMAX,NUMAX,MAXAA,0)
      call twh_ijaa2(GLB_INIT,NTMAX,NUMAX,MAXAA,0)
      call twh_inbnbn(GLB_INIT,NTMAX,NUMAX,0)
      call twh_itbnbn(GLB_INIT,NTMAX,NUMAX,MAXBEND,0)
      call twh_ijbnbn1(GLB_INIT,NTMAX,NUMAX,MAXBEND,0)
      call twh_ijbnbn2(GLB_INIT,NTMAX,NUMAX,MAXBEND,0)
      call twh_ibtbnbn1(GLB_INIT,NTMAX,NUMAX,MAXBEND,0)
      call twh_ibtbnbn2(GLB_INIT,NTMAX,NUMAX,MAXBEND,0)
      call twh_inof(GLB_INIT,NTMAX,NUMAX,0)
      call twh_itof(GLB_INIT,NTMAX,NUMAX,MAXOF,0)
      call twh_ijof5(GLB_INIT,NTMAX,NUMAX,MAXOF,0)
      call twh_ibttor1(GLB_INIT,NTMAX,NUMAX,MAXTOR,0)
      call twh_ibttor2(GLB_INIT,NTMAX,NUMAX,MAXTOR,0)
      call twh_ibttor3(GLB_INIT,NTMAX,NUMAX,MAXTOR,0)
      call twh_iattor1(GLB_INIT,NTMAX,NUMAX,MAXTOR,0)
      call twh_iattor2(GLB_INIT,NTMAX,NUMAX,MAXTOR,0)
      call twh_growcount(GLB_INIT,NTMAX,NUMAX,NNBOND,0)

c     --- initialize all of the structure arrays
      call twh_nunit(GLB_INIT,NTMAX,1)
      call twh_inring(GLB_INIT,NTMAX,NUMAX,0)
      call twh_ntype(GLB_INIT,NTMAX,NUMAX,0)
      call twh_chempot(GLB_INIT,NTMAX,1.0d0)
      call twh_chempotperstep(GLB_INIT,NTMAX,0)
      call twh_nmaxcbmc(GLB_INIT,NTMAX,1)
      call twh_cbmcbondinfo(GLB_INIT,NTMAX,NUMAX,NNBOND
     &     ,CBMCMAXBOND,0.0d0)
      call twh_cbmcbendcount(GLB_INIT,NTMAX,NUMAX,MAXBEND,MAXBEND,0)
      call twh_cbmcbenddist(GLB_INIT,NTMAX,NUMAX,MAXBEND,MAXBEND
     &     ,DIST_NULL)
      call twh_cbmcbendinfo(GLB_INIT,NTMAX,NUMAX,MAXBEND
     &     ,MAXBEND,CBMC_MAXBENDPEAKS,CBMCMAXBEND,0.0d0)
      call twh_ncmt(GLB_INIT,MAXBOX,NTMAX,0)
      call twh_nchbox(GLB_INIT,MAXBOX,0)

c     --- initialize energy bias cube variables
      call twh_lcreatemap(GLB_SET,.false.)
      call twh_cubelet(GLB_SET,0,0,0)
      call twh_mapmolty(GLB_SET,0)

c     --- initialize potential scaling.
      call scp_initialize()

c     --- initialize input parameters
      call twh_nmolty(GLB_SET,0)
      call twh_nmolectyp(GLB_INIT,NTMAX,0)
      call twh_charge_assignment(GLB_INIT,NTMAX,'manual')
      call twh_louthist(GLB_SET,.false.)
      call twh_histcalcfreq(GLB_SET,100000)
      call twh_histdumpfreq(GLB_SET,100000)
      call twh_temperature(GLB_SET,1.0d0)
      call twh_express(GLB_SET,1.0d0)
      call twh_stepstyle(GLB_SET,'cycles')
      call twh_nstep(GLB_SET,0)
      call twh_printfreq(GLB_SET,100)
      call twh_blocksize(GLB_SET,200)
      call twh_moviefreq(GLB_SET,100)
      call twh_backupfreq(GLB_SET,100)
      call twh_runoutput(GLB_SET,'full')
      call twh_pdb_output_freq(GLB_SET,0)
      call twh_loutdft(GLB_SET,.false.)
      call twh_loutlammps(GLB_SET,.false.)
      call twh_loutdlpoly(GLB_SET,.false.)
      call twh_pressurefreq(GLB_SET,10)
      call twh_trmaxdispfreq(GLB_SET,100)
      call twh_volmaxdispfreq(GLB_SET,100)
      call twh_potentialstyle(GLB_SET,POT_INTERNAL)
      call twh_ffnumber(GLB_SET,1)
      call twh_ff_filename(GLB_SET,1,'towhee_ff')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')
      call twh_ldualcutoff(GLB_SET,.true.)
      call twh_max_bond_length(GLB_SET,3.0d0)
      call twh_rcut(GLB_SET,1.0d0)
      call twh_rcutin(GLB_SET,1.0d0)
      call twh_rmin(GLB_SET,0.0d0)
      call twh_rminsq(GLB_SET,0.0d0)
      call twh_electrostatic_form(GLB_SET,'coulomb')
      call twh_coulombstyle(GLB_SET,'ewald_fixed_kmax')
      call twh_lcoulomb(GLB_SET,.true.)
      call twh_lewald(GLB_SET,.true.)
      call twh_ewald_kalp(GLB_SET,5.6d0)
      call twh_ewald_calp(GLB_INIT,MAXBOX,0.0d0)
      call twh_rcelect(GLB_INIT,MAXBOX,1.0d0)
      call twh_constant_qqfact(GLB_INIT,dvalue)
      call twh_ewald_prec(GLB_SET,1.0d-4)
      call twh_ltailc(GLB_SET,.false.)
      call twh_lshift(GLB_SET,.false.)
      call twh_nfield(GLB_SET,0)

c     --- initialize all ljfld parameters
      call twh_fieldint(GLB_INIT,FLDMAX,MAXFIELDINT,0)
      call twh_fielddp(GLB_INIT,FLDMAX,MAXFIELDDP,0.0d0)
      call twh_fieldnbparam(GLB_INIT,FLDMAX,NNTYPE,MAXFIELDNBPARAM
     &     ,0.0d0)

      call twh_solvation_style(GLB_SET,'none')
      call twh_isolvtype(GLB_SET,SOLV_NONE)
      call twh_linit(GLB_SET,.true.)
      call twh_initboxtype(GLB_SET,'dimensions')
      call twh_initstyle(GLB_INIT,MAXBOX,NTMAX,'coords')
      call twh_initlattice(GLB_INIT,MAXBOX,NTMAX,'none')
c     --- initialize hmatrix
      call twh_hmatrix(GLB_INIT,MAXBOX,3,3,0.0d0)
c     --- set diagonal elements to 100.0
      do ibox = 1,MAXBOX
         do idim = 1,3
            jdim = idim
            call twh_hmatrix(GLB_SET,ibox,idim,jdim,100.0d0)
         enddo
      enddo
c     --- initialize and then set hinverse
      call twh_hinverse(GLB_INIT,MAXBOX,3,3,0.0d0)
      do ibox = 1,MAXBOX
         call twh_inverthmatrix(ibox)
      enddo
      ivalue = 100
      do ibox = 1,MAXBOX
         inix = ivalue
         iniy = ivalue
         iniz = ivalue
         call twh_inixyz(GLB_SET,ibox,inix,iniy,iniz)
      enddo
      call twh_ensemble(GLB_SET,ENS_NVT)

c     --- Give initial values to move probabilities.
      call twh_pmvol(GLB_SET,0.0d0)
      call twh_pmvlpr(GLB_INIT,MAXPBOX,1.0d0)
      call twh_tavol(GLB_SET,0.5d0)
      call twh_pmcell(GLB_SET,0.0d0)
      call twh_pmcellpr(GLB_INIT,MAXPBOX,1.0d0)
      call twh_pmcellpt(GLB_INIT,MAXPBOX,0.5d0)
      call twh_rmcell(GLB_INIT,MAXPBOX,3,3,1.0d0)
      call twh_tacell(GLB_SET,0.5d0)
      
      call twh_pm2boxrbswap(GLB_SET,0.0d0)
      call twh_pm2rbswmt(GLB_INIT,NTMAX,1.0d0)
      call twh_pm2rbswpr(GLB_INIT,MAXPBOX,1.0d0)

      call twh_pm2boxcbswap(GLB_SET,0.0d0)
      call twh_pm2cbswmt(GLB_INIT,NTMAX,1.0d0)
      call twh_pm2cbswpr(GLB_INIT,MAXPBOX,1.0d0)

      call twh_pm1boxcbswap(GLB_SET,0.0d0)
      call twh_pm1cbswmt(GLB_INIT,NTMAX,1.0d0)
      
      call twh_pm1boxcomswitch(GLB_SET,0.0d0)
      call twh_pm1comswbox(GLB_INIT,MAXBOX,1.0d0)
      call twh_pm1comswpair(GLB_INIT,NTPAIRMAX,1.0d0)
      
      call twh_pmavb1(GLB_SET,0.0d0)
      call twh_pmavb1in(GLB_SET,0.5d0)
      call twh_pmavb1mt(GLB_INIT,NTMAX,1.0d0)
      call twh_pmavb1ct(GLB_INIT,NTMAX,NTMAX,1.0d0)
      call twh_avb1rad(GLB_SET,0.1d0)
      
      call twh_pmavb2(GLB_SET,0.0d0)
      call twh_pmavb2in(GLB_SET,0.5d0)
      call twh_pmavb2mt(GLB_INIT,NTMAX,1.0d0)
      call twh_pmavb2ct(GLB_INIT,NTMAX,NTMAX,1.0d0)
      call twh_avb2rad(GLB_SET,0.1d0)
      
      call twh_pmavb3(GLB_SET,0.0d0)
      call twh_pmavb3mt(GLB_INIT,NTMAX,1.0d0)
      call twh_pmavb3ct(GLB_INIT,NTMAX,NTMAX,1.0d0)
      call twh_avb3rad(GLB_SET,0.1d0)
      
      call twh_pmuvtcbswap(GLB_SET,0.0d0)
      call twh_pmuvtcbmt(GLB_INIT,NTMAX,1.0d0)
      
      call twh_pmcb(GLB_SET,0.0d0)
      call twh_pmcbmt(GLB_INIT,NTMAX,1.0d0)
      call twh_pmall(GLB_INIT,NTMAX,0.0d0)
      
      call twh_pmback(GLB_SET,0.0d0)
      call twh_pmbkmt(GLB_INIT,NTMAX,1.0d0)

      call twh_pmcbside(GLB_SET,0.0d0)
      call twh_pmcbsidemt(GLB_INIT,NTMAX,1.0d0)

      call twh_pmpivot(GLB_SET,0.0d0)
      call twh_pmpivmt(GLB_INIT,NTMAX,1.0d0)

      call twh_pmconrot(GLB_SET,0.0d0)
      call twh_pmcrmt(GLB_INIT,NTMAX,1.0d0)

      call twh_pmcrback(GLB_SET,0.0d0)
      call twh_pmcrbmt(GLB_INIT,NTMAX,1.0d0)
      
      call twh_pmplane(GLB_SET,0.0d0)
      call twh_pmplanebox(GLB_INIT,MAXBOX,1.0d0)
      call twh_planewidth(GLB_SET,1.0d0)
      
      call twh_pmrow(GLB_SET,0.0d0)
      call twh_pmrowbox(GLB_INIT,MAXBOX,1.0d0)
      call twh_rowwidth(GLB_SET,1.0d0)
      
      call twh_pmtraat(GLB_SET,0.0d0)
      call twh_pmtamt(GLB_INIT,NTMAX,1.0d0)
      call twh_tatraa(GLB_SET,0.5d0)
      
      call twh_pmtracm(GLB_SET,0.0d0)
      call twh_pmtcmt(GLB_INIT,NTMAX,1.0d0)
      call twh_tatrac(GLB_SET,0.5d0)

      call twh_pmrotate(GLB_SET,0.0d0)
      call twh_pmromt(GLB_INIT,NTMAX,1.0d0)
      call twh_tarot(GLB_SET,0.5d0)

      call twh_pmcomposite(GLB_SET,0.0d0)
      call twh_pmcomt(GLB_INIT,NTMAX,1.0d0)

c     --- initialize vibration parameters
      call twh_vibcoeff(GLB_INIT,TVIBMAX,MAXVIBCOEFF,0.0d0)
      call twh_vibstyle(GLB_INIT,TVIBMAX,BOND_STYLE_NULL)
      ctemp = 'init-error'
      call twh_vibff(GLB_INIT,TVIBMAX,ctemp)
      call twh_vibnames(GLB_INIT,TVIBMAX,MAXSAMEV,MAXVNAME,ctemp)
      call twh_nvsame(GLB_INIT,TVIBMAX,0)
c     --- default the viborder to wild so it works unless we are
c     --- deliberately being more sophisticated
      ctemp = 'wild'
      call twh_viborder(GLB_INIT,TVIBMAX,ctemp)
      call twh_angleorder(GLB_INIT,TBENMAX,ctemp)
      call twh_torsionorder(GLB_INIT,TTORMAX,ctemp)
      call twh_biorder(GLB_INIT,TBIMAX,ctemp)

c     --- initialize bond increments
      ctemp = 'init-error'
      call twh_qbiff(GLB_INIT,TBIMAX,ctemp)
      call twh_qbivalue(GLB_INIT,TBIMAX,0.0d0)
      call twh_qbinames(GLB_INIT,TBIMAX,2,ctemp)

      ctemp = 'init-error'
      call twh_lbonang(GLB_INIT,TBENMAX,.false.)
      call twh_lbonbon(GLB_INIT,TBENMAX,.false.)
      call twh_anglestyle(GLB_INIT,TBENMAX,ANG_STYLE_NULL)
      call twh_bencoeff(GLB_INIT,TBENMAX,MAXBENDCOEFF,0.0d0)
      call twh_nbsame(GLB_INIT,TBENMAX,0)
      call twh_bendff(GLB_INIT,TBENMAX,ctemp)
      call twh_bendnames(GLB_INIT,TBENMAX,MAXSAMEB,MAXBNAME,ctemp)
      call twh_aastyle(GLB_INIT,TAAMAX,-1)
      call twh_aacoeff(GLB_INIT,TAAMAX,MAXAACOEFF,0.0d0)

      call twh_loftor(GLB_INIT,TTORMAX,.false.)
      call twh_ntorloop(GLB_INIT,TTORMAX,0)
      call twh_torstyle(GLB_INIT,TTORMAX,0)
      call twh_torstring(GLB_INIT,TTORMAX,4,ctemp)
      call twh_tornames(GLB_INIT,TTORMAX,MAXSAMET,MAXTNAME,ctemp)
      call twh_torff(GLB_INIT,TTORMAX,ctemp)
      call twh_ntsame(GLB_INIT,TTORMAX,0)
      call twh_torcoeff(GLB_INIT,TTORMAX,MAXTORCOEFF,0.0d0)

      call twh_impstyle(GLB_INIT,TIMPMAX,0)
      call twh_impform(GLB_INIT,TIMPMAX,0)
      call twh_impcoeff(GLB_INIT,TIMPMAX,MAXIMPCOEFF,0.0d0)
      call twh_nimpsame(GLB_INIT,TIMPMAX,0)
      ctemp = 'init-error'
      call twh_impnames(GLB_INIT,TIMPMAX,MAXSAMEIMP,MAXIMPNAME,ctemp)
c     --- initialize special one-five 
      call twh_ofstyle(GLB_INIT,TOFMAX,0)
      call twh_ofcoeff(GLB_INIT,TOFMAX,MAXOFCOEFF,0.0d0)
      ctemp = 'init-error'
      call twh_ofnames(GLB_INIT,TOFMAX,MAXOFNAME,ctemp)

c     --- initialize hbond
      ctemp = 'init-error'
      call twh_hbondnames(GLB_INIT,THBONDMAX,MAXHBONDNAME,ctemp)
      call twh_hbondcoeff(GLB_INIT,THBONDMAX,MAXHBONDCOEFF,0.0d0)

c     --- initialize nonbonded parameters for traditional potentials
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      ij = twh_getnbtype(NNTYPE,NNTYPE)
      call twh_nbcoeff(GLB_INIT,MAXNBCOEFF,ij,0.0d0)
      call twh_mass(GLB_INIT,NNTYPE,0.0d0)
      call twh_elementname(GLB_INIT,NNTYPE,'NA')
      call twh_donor_acceptor(GLB_INIT,NNTYPE,'-')
      call twh_bondpatt(GLB_INIT,NNTYPE,'null')
      call twh_quantum_gsenergy(GLB_INIT,NTMAX,0.0d0)
      ctemp = 'init-error'
      call twh_nbff(GLB_INIT,NNTYPE,ctemp)
      call twh_nbname(GLB_INIT,NNTYPE,ctemp)
      call twh_caname(GLB_INIT,NNTYPE,ctemp)
      call twh_cbname(GLB_INIT,NNTYPE,ctemp)
      call twh_ctname(GLB_INIT,NNTYPE,ctemp)

      call twh_qbasevalue(GLB_INIT,NNTYPE,0.0d0)
      call twh_polarizability(GLB_INIT,NNTYPE,1.0d0)
      call twh_nb_to_hbtype(GLB_INIT,CROSSTYPEMAX,0)

c     --- initialize eam parameters
      call twh_classical_potential(GLB_SET,'Embedded Atom Method')
      call twh_table_npair(GLB_INIT
     &     ,((MAXTABTYPE*MAXTABTYPE+MAXTABTYPE)/2),0)
      call twh_eam_ndens(GLB_INIT,MAXTABTYPE,MAXTABTYPE,0)
      call twh_eam_nembed(GLB_INIT,MAXTABTYPE,0)
      call twh_eam_pair_style(GLB_INIT
     &     ,(MAXTABTYPE*MAXTABTYPE+MAXTABTYPE)/2,'none')
      call twh_table_pair(GLB_INIT,3
     &     ,((MAXTABTYPE*MAXTABTYPE+MAXTABTYPE)/2),MAXTABLE,0.0d0)
      call twh_eam_dens(GLB_INIT,3,MAXTABTYPE,MAXTABTYPE,MAXTABLE,0.0d0)
      call twh_eam_embed(GLB_INIT,3,MAXTABTYPE,MAXTABLE,0.0d0)
c     --- reinitialize classical_potential and number
      call twh_classical_potential(GLB_SET,'init-error')
      call twh_cbmc_nb_one_generation(GLB_INIT,MAXBOX,NTMAX,'uniform')
      call twh_cbmc_nb_one_gencode(GLB_INIT,MAXBOX,NTMAX,DIST_UNIFORM)
      call twh_cbmc_bond_generation(GLB_SET,'ideal')
      call twh_cbmc_bondgen_code(GLB_SET,CBMC_IDEAL)
      call twh_bond_sdev_multiplier(GLB_SET,1.0d0)
      call twh_cbmc_bend_generation(GLB_SET,'ideal')
      call twh_cbmc_bendgen_code(GLB_SET,CBMC_IDEAL)
      call twh_bend_a_sdev_multiplier(GLB_SET,1.0d0)
      call twh_bend_b_sdev_multiplier(GLB_SET,1.0d0)
      call twh_cbmc_dihedral_generation(GLB_SET,'ideal')
      call twh_cbmc_dihedgen_code(GLB_SET,CBMC_IDEAL)
      call twh_dihedral_peak_weight_style(GLB_SET,'init-error')
      call twh_cbmctorcount(GLB_INIT,NTMAX,NUMAX,MAXTOR,0)
      call twh_dihedral_peak_weights(GLB_INIT,NTMAX,NUMAX,MAXTOR
     &     ,CBMC_MAXTORPEAKS,1.0d0)

      call twh_nch_tor(GLB_INIT,NTMAX,360) 
      call twh_nch_nb_one(GLB_INIT,NTMAX,10)
      call twh_nch_nb(GLB_INIT,NTMAX,10)
      call twh_nch_tor_connect(GLB_INIT,NTMAX,360)
      call twh_nch_bend_a(GLB_INIT,NTMAX,100)
      call twh_nch_bend_b(GLB_INIT,NTMAX,100)
      call twh_nch_vib(GLB_INIT,NTMAX,1000)
      call twh_nch_pre_nb(GLB_INIT,NTMAX,1)
      call twh_cbmc_formulation(GLB_SET,'Martin and Frischknecht 2006')
      call twh_cbmc_analysis(GLB_SET,'normal')
      call twh_cbmc_setting_style(GLB_SET,'Martin and Frischknecht')
      call twh_cbmc_fit_strategy(GLB_SET,'equilibrium')
      call twh_two_bond_fe_bias_style(GLB_SET,'none')
      call twh_two_bond_bias_vibrange(GLB_SET,1,0.85d0)
      call twh_two_bond_bias_vibrange(GLB_SET,2,1.15d0)
      call twh_two_bond_bias_sdev_multiplier(GLB_SET,1.0d0)
      call twh_two_bond_bias_initial_value(GLB_SET,'none')
      call twh_two_bond_bias_update_frequency(GLB_SET,0)
      call twh_two_bond_bias_compute_frequency(GLB_SET,0)
      call twh_two_bond_bias_old_fraction(GLB_SET,0.5d0)
      call twh_two_bond_bias_info(GLB_INIT,NTMAX,NUMAX
     &     ,MAXBEND,MAXTWOBOND,0.0d0)
      call twh_three_bond_bias_info(GLB_INIT,NTMAX,NUMAX
     &     ,MAXTOR,MAXTHREEBOND,0.0d0)
      call twh_three_bond_fe_bias_style(GLB_SET,'none')
      call twh_three_bond_bias_initial_value(GLB_SET,'none')
      call twh_three_bond_bias_update_freq(GLB_SET,0)
      call twh_three_bond_bias_compute_freq(GLB_SET,0)
      call twh_three_bond_bias_old_fraction(GLB_SET,0.5d0)
      call twh_three_bond_bias_sdev_multiplier(GLB_SET,1.0d0)
      dvalue = 0.85d0**3
      call twh_vibrang(GLB_SET,1,dvalue)
      dvalue = 1.15d0**3 - dvalue
      call twh_vibrang(GLB_SET,2,dvalue)

      call twh_constant_twopi(GLB_GET,twopi)
      dvalue = 10.0d0*twopi/180.0d0
      call twh_sdevtor(GLB_SET,dvalue)
      dvalue = 5.0d0*twopi/180.0d0
      call twh_sdevbena(GLB_SET,dvalue)
      dvalue = 5.0d0*twopi/180.0d0
      call twh_sdevbenb(GLB_SET,dvalue)
      call twh_sdevvib(GLB_SET,0.1d0)

      call twh_lusecom(GLB_SET,.true.)

      call twh_minboxdim(GLB_SET,0.0d0)
c     --- initialize temporary accumulators
      call twh_baavb1(GLB_INIT,NTMAX,NTMAX,MAXBOX,3,0)
      call twh_bnavb1(GLB_INIT,NTMAX,NTMAX,MAXBOX,3,0)
      call twh_baavb2(GLB_INIT,NTMAX,NTMAX,MAXBOX,3,0)
      call twh_bnavb2(GLB_INIT,NTMAX,NTMAX,MAXBOX,3,0)
      call twh_baavb3(GLB_INIT,NTMAX,NTMAX,MAXBOX,0)
      call twh_bnavb3(GLB_INIT,NTMAX,NTMAX,MAXBOX,0)
      call twh_baback(GLB_INIT,NTMAX,NUMAX,MAXTOR,MAXBOX,0)
      call twh_bnback(GLB_INIT,NTMAX,NUMAX,MAXTOR,MAXBOX,0)
      call twh_bsback(GLB_INIT,NTMAX,NUMAX,MAXTOR,MAXBOX,0)
      call twh_bacb(GLB_INIT,NTMAX,NUMAX,MAXBOX,0)
      call twh_bncb(GLB_INIT,NTMAX,NUMAX,MAXBOX,0)
      call twh_bscb(GLB_INIT,NTMAX,NUMAX,MAXBOX,0)
      call twh_bacbside(GLB_INIT,NTMAX,NUMAX,MAXBOX,0)
      call twh_bncbside(GLB_INIT,NTMAX,NUMAX,MAXBOX,0)
      call twh_bscbside(GLB_INIT,NTMAX,NUMAX,MAXBOX,0)
      call twh_bacbswap(GLB_INIT,NTMAX,MAXBOX,MAXBOX,0)
      call twh_bncbswap(GLB_INIT,NTMAX,MAXBOX,MAXBOX,0)
      call twh_bscbswap(GLB_INIT,NTMAX,MAXBOX,MAXBOX,0)
      call twh_baconrot(GLB_INIT,NTMAX,NUMAX,MAXBOX,0)
      call twh_bnconrot(GLB_INIT,NTMAX,NUMAX,MAXBOX,0)
      call twh_bacrback(GLB_INIT,NTMAX,NUMAX,MAXBOX,0)
      call twh_bncrback(GLB_INIT,NTMAX,NUMAX,MAXBOX,0)
      call twh_bapivot(GLB_INIT,NTMAX,NUMAX,NNBOND,MAXBOX,0)
      call twh_bnpivot(GLB_INIT,NTMAX,NUMAX,NNBOND,MAXBOX,0)
      call twh_bnpivbin(GLB_INIT,MAXBOX,MAXPIVOTBIN,0)
      call twh_bapivbin(GLB_INIT,MAXBOX,MAXPIVOTBIN,0)
      call twh_barbswap(GLB_INIT,NTMAX,MAXBOX,MAXBOX,0)
      call twh_bnrbswap(GLB_INIT,NTMAX,MAXBOX,MAXBOX,0)
      call twh_bashift(GLB_INIT,2,MAXBOX,3,0)
      call twh_bnshift(GLB_INIT,2,MAXBOX,3,0)
      call twh_bavol(GLB_INIT,MAXPBOX,0)
      call twh_bnvol(GLB_INIT,MAXPBOX,0)

c     --- initialize use coulomb & VdW logical arrays
      call twh_luse_c(GLB_INIT,NTMAX,NUMAX,.true.)
      call twh_luse_v(GLB_INIT,NTMAX,NUMAX,.true.)

c     --- initialize aminonum
      call twh_aminonum(GLB_INIT,NTMAX,NUMAX,0)

c     --- initialize interal storage in averages
      call twh_averages(lfailure,AVG_INITIALIZE,0,MAXBOX
     &     ,NTMAX,dvalue,ivalue)

      return
      end
