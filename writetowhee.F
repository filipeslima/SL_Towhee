#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2011 Marcus G. Martin                       *
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
      subroutine twh_writetowhee(lfailure,lprotein,imolty)
c     ******************************************************************
c     * outputs information suitable for inpstyle 0 to towhee_altinp   *
c     *                                                                *
c     * rewritten from assemble 03-19-2002 by M.G. Martin              *
c     * last modified 02-05-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_CLASSICAL_MIXRULE
#define FUNCTION_GET_CLASSICAL_POTENTIAL
#define FUNCTION_GET_ELEMENTNAME
#define FUNCTION_GET_IJAA0
#define FUNCTION_GET_IJAA1
#define FUNCTION_GET_IJAA2
#define FUNCTION_GET_IJBEN2
#define FUNCTION_GET_IJBEN3
#define FUNCTION_GET_IJIMPROP2
#define FUNCTION_GET_IJIMPROP3
#define FUNCTION_GET_IJIMPROP4
#define FUNCTION_GET_IJTOR2
#define FUNCTION_GET_IJTOR3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_IJVIB
#define FUNCTION_GET_INAA
#define FUNCTION_GET_INBEN
#define FUNCTION_GET_INIMPROP
#define FUNCTION_GET_INTOR
#define FUNCTION_GET_INVIB
#define FUNCTION_GET_ITAA
#define FUNCTION_GET_ITBEN
#define FUNCTION_GET_ITIMPROP
#define FUNCTION_GET_ITTOR
#define FUNCTION_GET_ITVIB
#define FUNCTION_GET_NATIVE_LAMBDA
#define FUNCTION_GET_NBNAME
#define FUNCTION_GET_NTYPE
#define FUNCTION_GET_PMVLPR
#define FUNCTION_GET_PMVOL
#define FUNCTION_GET_RMVOL
#define FUNCTION_GET_SCALING_STYLE
#define FUNCTION_GET_SCALING_STYLE_STRING
#define FUNCTION_GET_TAVOL
#define FUNCTION_GET_QQATOM
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure,lprotein
      integer imolty
c     --- local variables
      character*1 hist_suffix
      character*3 aminoshort
      character*4 pdbname
      character*8 stepstyle
      character*10 tenchar
      character*11 fchar
      character*20 routput,cvalue,coulombstyle,fieldstring,initboxtype
      character*20 initstyle,initlattice,interpolatestyle
      character*20 isarray
      dimension isarray(NTMAX)
      character*30 potentialstyle_str,ensemble_str
      character*50 eform
      character*(MAXDIRLENGTH) ffname
c     --- logical scalars
      logical lfound,ldoinitstyle,loutdft,louthist,lshift
      logical lwritemixrule,lwritelshift,lwriteltailc,lwritermin
      logical lwritercut,lwritercutin,lwriteinterpolatestyle,lwriterpd
      logical lwritescaling,lvalue,loutdlpoly,loutlammps,ltailc
c     --- character string arrays
      character*10 pair_list
      dimension pair_list(CROSSTYPEMAX,2)
c     --- small integer scalars
      integer*4 random_seed,luxlevel
c     --- integer scalars
      integer iatom,ivib,ibend,itor,iaa,ibox,jmolty,idim,jdim,iorder
      integer kmolty,ifile,ifld,itype,ntii,ientry,ipair,rwfile,nmolty
      integer nboxpair,nunit,kmax,nfield,ivalue,ffnumber,npairs
      integer aminonum,potentialstyle,numboxes,rng_code
      integer npair,inix,iniy,iniz,nvmove,ensemble
c     --- integer arrays
      integer idummy
      dimension idummy(NTMAX)
c     --- double precision scalars
      double precision dvalue,qqfact,temperature,ewald_prec,kalp,rcut
      double precision fieldtwo,fieldthree,fieldfour,fnbparam,pmcb,rmin
      double precision rcutin
c     --- double precision arrays
      double precision dpboxarray
      dimension dpboxarray(MAXPBOX)
      double precision dmarray
      dimension dmarray(NTMAX)
      double precision dpairarray
      dimension dpairarray(NTPAIRMAX)
      double precision dboxarray
      dimension dboxarray(MAXBOX)
      double precision dimtemp
      dimension dimtemp(3)

c      write(6,*) 'starting WRITETOWHEE ',imolty
c     --- retrieve constants
      call twh_potentialstyle(GLB_GET,potentialstyle)
      call twh_numboxes(GLB_GET,numboxes)
      call twh_nvmove(GLB_GET,nvmove)
      call twh_ensemble(GLB_GET,ensemble)

      if ( imolty .eq. 0 ) then
c        --- output all of the towhee_input file except for the 
c        --- molecule information
         write(15,'(a11)') 'inputformat'
         write(15,'(a1,a6,a1)') char(39),'Towhee',char(39)
         call twh_random_code(GLB_GET,rng_code)
         write(15,'(a23)') 'random_number_generator'
         if ( rng_code .eq. RNG_RANLUX ) then
            write(15,'(a1,a6,a1)') char(39),'RANLUX',char(39)
            call twh_random_luxlevel(GLB_GET,luxlevel)
            write(15,'(a15)') 'random_luxlevel'
            write(15,*) luxlevel
         elseif ( rng_code .eq. RNG_DX_1597_2_7 ) then
            write(15,'(a1,a11,a1)') char(39),'DX-1597-2-7',char(39)
         elseif ( rng_code .eq. RNG_KISS99 ) then
            write(15,'(a1,a6,a1)') char(39),'KISS99',char(39)
         elseif ( rng_code .eq. RNG_MRG32K3A ) then
            write(15,'(a1,a8,a1)') char(39),'MRG32k3a',char(39)
         elseif ( rng_code .eq. RNG_DEBUG ) then
            write(15,'(a1,a5,a1)') char(39),'debug',char(39)
         else
            write(6,*) 'WRITETOWHEE: unknown random_code'
            write(6,*) rng_code
            lfailure = .true.
            return
         endif
         call twh_random_seed(GLB_GET,random_seed)
         write(15,'(a11)') 'random_seed'
         write(15,*) random_seed
         call twh_random_allow_restart(GLB_GET,lvalue)
         write(15,'(a20)') 'random_allow_restart'
         write(15,'(L1)') lvalue

         if (ensemble .eq. ENS_NPT) ensemble_str='npt'
         if (ensemble .eq. ENS_NVT) ensemble_str='nvt'
         if (ensemble .eq. ENS_muVT) ensemble_str='uvt'
         write(15,'(a8)') 'ensemble'
         write(15,'(a1,a30,a1)') char(39),ensemble_str,char(39)

         write(15,'(a11)') 'temperature'
         call twh_temperature(GLB_GET,temperature)
         write(15,*) temperature
         if ( ensemble .eq. ENS_NPT ) then
            write(15,'(a8)') 'pressure'
            call twh_express(GLB_GET,dvalue)
            dvalue = dvalue/CON_KPATOSIM
            write(15,*) dvalue
         endif
         write(15,'(a6)') 'nmolty'
         call twh_nmolty(GLB_GET,nmolty)
         write(15,*) nmolty
         write(15,'(a9)') 'nmolectyp'
         do jmolty = 1,nmolty
            call twh_nmolectyp(GLB_GET,jmolty,idummy(jmolty))
         enddo
         write(15,*) (idummy(jmolty),jmolty=1,nmolty)
         if ( ensemble .eq. ENS_muVT ) then
            write(15,'(a7)') 'chempot'
            do jmolty = 1,nmolty
               call twh_chempot(GLB_GET,jmolty,dmarray(jmolty))
            enddo
            write(15,*) (dmarray(jmolty),jmolty=1,nmolty)
         endif
         write(15,'(a8)') 'numboxes'
         write(15,*) numboxes
         write(15,'(a9)') 'stepstyle'
         call twh_stepstyle(GLB_GET,stepstyle)
         write(15,'(a1,a8,a1)') char(39),stepstyle,char(39)
         if ( stepstyle .eq. 'cycles' .or. stepstyle .eq. 'moves' ) then
            write(15,'(a5)') 'nstep'
            call twh_nstep(GLB_GET,ivalue)
            write(15,*) ivalue
         else
            write(6,*) 'WRITETOWHEE: unknown stepstyle'
         endif
         write(15,'(a9)') 'printfreq'
         call twh_printfreq(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a9)') 'blocksize'
         call twh_blocksize(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a9)') 'moviefreq'
         call twh_moviefreq(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a10)') 'backupfreq'
         call twh_backupfreq(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a11)') 'restartfreq'
         call twh_restartfreq(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a9)') 'runoutput'
         call twh_runoutput(GLB_GET,routput)
         write(15,'(a1,a20,a1)') char(39),routput,char(39)
         write(15,'(a15)') 'pdb_output_freq'
         call twh_pdb_output_freq(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a7)') 'loutdft'
         call twh_loutdft(GLB_GET,loutdft)
         write(15,*) loutdft
         write(15,'(a10)') 'loutlammps'
         call twh_loutlammps(GLB_GET,loutlammps)
         write(15,*) loutlammps
         write(15,'(a10)') 'loutdlpoly'
         call twh_loutdlpoly(GLB_GET,loutdlpoly)
         write(15,*) loutdlpoly
         if ( ensemble .eq. ENS_muVT ) then
            write(15,'(a8)') 'louthist'
            call twh_louthist(GLB_GET,louthist)
            write(15,*) louthist
            if ( louthist ) then
               write(15,'(a10)') 'hist_label'
               call twh_hist_label(GLB_GET,ivalue)
               write(15,*) ivalue
               write(15,'(a11)') 'hist_suffix'
               call twh_hist_suffix(GLB_GET,hist_suffix)
               write(15,'(a1,a1,a1)') char(39),hist_suffix,char(39)
               write(15,'(a11)') 'hist_nequil'
               call twh_hist_nequil(GLB_GET,ivalue)
               write(15,*) ivalue
               write(15,'(a12)') 'histcalcfreq'
               call twh_histcalcfreq(GLB_GET,ivalue)
               write(15,*) ivalue
               write(15,'(a12)') 'histdumpfreq'
               call twh_histdumpfreq(GLB_GET,ivalue)
               write(15,*) ivalue
            endif
         endif
         write(15,'(a12)') 'pressurefreq'
         call twh_pressurefreq(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a13)') 'trmaxdispfreq'
         call twh_trmaxdispfreq(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a14)') 'volmaxdispfreq'
         call twh_volmaxdispfreq(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a14)') 'chempotperstep'
         do jmolty = 1,nmolty
            call twh_chempotperstep(GLB_GET,jmolty,idummy(jmolty))
         enddo
         write(15,*) (idummy(jmolty),jmolty=1,nmolty)
         write(15,'(a14)') 'potentialstyle'

         if ( potentialstyle .eq. POT_INTERNAL) then
            potentialstyle_str = 'internal'
         elseif ( potentialstyle .eq. POT_EXTERNAL ) then
            potentialstyle_str = 'external'
         else
            write(6,*) 'Unrecognized potentialstyle in WRITETOWHEE'
            stop
         endif    
         write(15,'(a1,a30,a1)') char(39),potentialstyle_str,char(39)

         if ( potentialstyle .eq. POT_INTERNAL ) then
            write(15,'(a8)') 'ffnumber'
            call twh_ffnumber(GLB_GET,ffnumber)
            write(15,*) ffnumber
            if ( ffnumber .gt. 0 ) then
               write(15,'(a11)') 'ff_filename'
               do ifile = 1,ffnumber
                  call twh_ff_filename(GLB_GET,ifile,ffname)
                  write(15,11) ffname
               enddo
            endif
 11         format(a MAXDIRLENGTH)
            write(15,'(a19)') 'classical_potential'
            write(15,'(a1,a30,a1)') char(39)
     &           ,twh_get_classical_potential(),char(39)

            call twh_setpotentiallogic(lfailure,lwritemixrule
     &           ,lwritelshift,lwriteltailc,lwritermin,lwritercut
     &           ,lwritercutin,lwriteinterpolatestyle,lwriterpd
     &           ,lwritescaling)
            if ( lfailure ) return

            if ( lwritemixrule ) then
               write(15,'(a17)') 'classical_mixrule'
               write(15,'(a1,a30,a1)') char(39)
     &              ,twh_get_classical_mixrule(),char(39)

               if ( lwritescaling ) then
                  write(15,'(a20)') 'cmix_rescaling_style'
                  write(15,'(a1,a30,a1)') char(39),
     &                 scp_get_scaling_style_string() ,char(39)
                  if ( scp_get_scaling_style() .eq. SCP_SC ) then
                     write(15,'(a14)') 'cmix_lambda_lj'
                     write(15,*) scp_get_native_lambda(LAMBDA_LJ)
                     write(15,'(a13)') 'cmix_alpha_lj'
                     write(15,*) scp_get_native_lambda(ALPHA_LJ)
                     write(15,'(a13)') 'cmix_lambda_c'
                     write(15,*) scp_get_native_lambda(LAMBDA_C)
                     write(15,'(a15)') 'cmix_lprintdudl'
                     call twh_scp_eval_derivative(GLB_GET,lvalue)
                     write(15,*) lvalue
                  endif
                  if (scp_get_scaling_style() .eq. SCP_SC .or.
     &                 scp_get_scaling_style() .eq. SCP_WCA ) then
                     call scp_get_scaled_atoms(npair, pair_list)
                     write(15,'(a10)') 'cmix_npair'
                     write(15,*) npair
                     write(15,'(a14)') 'cmix_pair_list'
                     do ipair = 1,npair
                        write(15,'(2(a1,a10,a1,1x))') (char(39)
     &                       ,pair_list(ipair,ientry),char(39)
     &                       ,ientry=1,2)
                     enddo
                  elseif ( scp_get_scaling_style() .eq. SCP_NONE ) then
c                    --- nothing to output
                  else
                     write(6,*) 
     &                    'WRITETOWHEE: unknown cmix_rescaling_style'
                     stop
                  endif
               endif
            endif
            if ( lwriteinterpolatestyle ) then
               write(15,'(a16)') 'interpolatestyle'
               call twh_interpolatestyle(GLB_GET,interpolatestyle)
               write(15,'(a1,a20,a1)') char(39),interpolatestyle
     &              ,char(39)
            endif
            if ( lwritelshift ) then
               write(15,'(a6)') 'lshift'
               call twh_lshift(GLB_GET,lshift)
               write(15,*) lshift
            endif
            if ( lwriteltailc ) then
               write(15,'(a6)') 'ltailc'
               call twh_ltailc(GLB_GET,ltailc)
               write(15,*) ltailc
            endif
            if ( lwritermin ) then
               write(15,'(a4)') 'rmin'
               call twh_rmin(GLB_GET,rmin)
               write(15,*) rmin
            endif
            if ( lwritercut ) then
               write(15,'(a4)') 'rcut'
               call twh_rcut(GLB_GET,rcut)
               write(15,*) rcut
            endif
            if ( lwritercutin ) then
               write(15,'(a6)') 'rcutin'
               call twh_rcutin(GLB_GET,rcutin)
               write(15,*) rcutin
            endif
            if ( lwriterpd ) then
               write(15,'(a21)') 'radial_pressure_delta'
               call twh_radial_pressure_delta(GLB_GET,dvalue)
               write(15,*) dvalue
            endif
            write(15,'(a18)') 'electrostatic_form'
            call twh_electrostatic_form(GLB_GET,eform)
            write(15,'(a1,a50,a1)') char(39),eform,char(39)
            if ( eform .eq. 'none' ) then
c              --- nothing else to write here, but not an error
            elseif ( eform .eq. 'coulomb' ) then
               write(15,'(a12)') 'coulombstyle'
               call twh_coulombstyle(GLB_GET,coulombstyle)
               write(15,'(a1,a20,a1)') char(39),coulombstyle,char(39)
               if ( coulombstyle .eq. 'ewald_fixed_kmax' ) then
                  write(15,'(a4)') 'kalp'
                  call twh_ewald_kalp(GLB_GET,kalp)
                  write(15,*) kalp
                  call twh_ewald_kmax(GLB_GET,1,kmax)
                  write(15,'(a4)') 'kmax'
                  write(15,*) kmax
               elseif ( coulombstyle .eq. 'ewald_fixed_cutoff' ) then
                  write(15,'(a10)') 'ewald_prec'
                  call twh_ewald_prec(GLB_GET,ewald_prec)
                  write(15,*) ewald_prec
                  write(15,'(a7)') 'rcelect'
                  call twh_rcelect(GLB_GET,1,dvalue)
                  write(15,*) dvalue
               elseif ( coulombstyle .eq. 'minimum image' ) then
c                 --- not an error, but nothing to output
               else
                  write(6,*) 'WRITETOWHEE: unknown coulombstyle'
                  stop
               endif

               write(15,'(a7)') 'dielect'
               call twh_constant_qqfact(GLB_GET,qqfact)
               write(15,*) CON_EETOKCAL*CON_KCALTOK/qqfact
            else
               write(6,*) 'WRITETOWHEE: unknown electrostatic_form'
               write(6,*) eform
               lfailure = .true.
               return
            endif

            write(15,'(a15)') 'max_bond_length'
            call twh_max_bond_length(GLB_GET,dvalue)
            write(15,*) dvalue

            write(15,'(a6)') 'nfield'
            call twh_nfield(GLB_GET,nfield)
            write(15,*) nfield
            do ifld = 1,nfield
               write(15,'(a9)') 'fieldtype'
               call twh_fieldtypestring(GLB_GET,ifld,fieldstring)
               write(15,'(a1,a20,a1)') char(39),fieldstring,char(39)
               if ( fieldstring .eq. 'Hard Wall' ) then
                  write(15,'(a6)') 'hrdbox'
                  call twh_fieldint(GLB_GET,ifld,1,ivalue)
                  write(15,*) ivalue
                  write(15,'(a6)') 'hrdxyz'
                  call twh_fieldchar(GLB_GET,ifld,1,fchar)
                  write(15,'(3a1)') char(39),fchar,char(39)
                  write(15,'(a6)') 'hrdcen'
                  call twh_fielddp(GLB_GET,ifld,1,dvalue)
                  write(15,*) dvalue
                  write(15,'(a6)') 'hrdrad'
                  call twh_fielddp(GLB_GET,ifld,2,dvalue)
                  write(15,*) dvalue
                  write(15,'(a19)') 'hrd_repulsion_style'
                  call twh_fieldchar(GLB_GET,ifld,2,fchar)
                  write(15,'(a1,a11,a1)') char(39),fchar,char(39)
                  write(15,'(a15)') 'hrd_energy_type'
                  call twh_fieldchar(GLB_GET,ifld,3,fchar)
                  write(15,'(a1,a11,a1)') char(39),fchar,char(39)
                  if ( fchar .eq. 'finite' ) then
                     write(15,'(a15)') 'hrd_wall_energy'
                     call twh_fielddp(GLB_GET,ifld,4,dvalue)
                     write(15,*) dvalue
                  endif
               elseif ( fieldstring .eq. 'Hard Ball' ) then
                  write(15,'(a6)') 'hrbbox'
                  call twh_fieldint(GLB_GET,ifld,FLD_HRBBOX,ivalue)
                  write(15,*) ivalue
                  write(15,'(a6)') 'hrbpos'
                  call twh_fielddp(GLB_GET,ifld,FLD_HRBPOSX,fieldtwo)
                  call twh_fielddp(GLB_GET,ifld,FLD_HRBPOSY,fieldthree)
                  call twh_fielddp(GLB_GET,ifld,FLD_HRBPOSZ,fieldfour)
                  write(15,*) fieldtwo, fieldthree, fieldfour
                  write(15,'(a6)') 'hrbrad'
                  call twh_fielddp(GLB_GET,ifld,FLD_HRBRAD,dvalue)
                  write(15,*) dvalue
                  write(15,'(a6)') 'hrbfac'
                  call twh_fielddp(GLB_GET,ifld,FLD_HRBFAC,dvalue)
                  write(15,*) dvalue
                  write(15,'(a15)') 'hrb_energy_type'
                  call twh_fieldchar(GLB_GET,ifld,1,fchar)
                  write(15,'(a1,a11,a1)') char(39),fchar,char(39)
                  if ( fchar .eq. 'finite' ) then
                     write(15,'(a15)') 'hrb_energy_type'
                     call twh_fielddp(GLB_GET,ifld,FLD_HRBENERGY,dvalue)
                     write(15,*) dvalue
                  endif
               elseif ( fieldstring .eq. 'Hooper Umbrella' ) then
                  write(15,'(a6)') 'umbbox'
                  call twh_fieldint(GLB_GET,ifld,1,ivalue)
                  write(15,*) ivalue
                  write(15,'(a6)') 'umbxyz'
                  call twh_fieldchar(GLB_GET,ifld,1,fchar)
                  write(15,'(3a1)') char(39),fchar,char(39)
                  write(15,'(a9)') 'umbcenter'
                  call twh_fielddp(GLB_GET,ifld,1,dvalue)
                  write(15,*) dvalue
                  write(15,'(a4)') 'umba'
                  call twh_fielddp(GLB_GET,ifld,2,dvalue)
                  write(15,*) dvalue
               elseif ( fieldstring .eq. 'LJ 9-3 Wall' ) then
                  write(15,'(a6)') 'ljfbox'
                  call twh_fieldint(GLB_GET,ifld,1,ivalue)
                  write(15,*) ivalue
                  write(15,'(a6)') 'ljfxyz'
                  call twh_fieldchar(GLB_GET,ifld,1,fchar)
                  write(15,'(3a1)') char(39),fchar,char(39)
                  write(15,'(a6)') 'ljfcen'
                  call twh_fielddp(GLB_GET,ifld,1,dvalue)
                  write(15,*) dvalue
                  write(15,'(a6)') 'ljfdir'
                  call twh_fieldint(GLB_GET,ifld,2,ivalue)
                  write(15,*) ivalue
                  write(15,'(a6)') 'ljfcut'
                  call twh_fielddp(GLB_GET,ifld,2,dvalue)
                  write(15,*) dvalue
                  write(15,'(a8)') 'ljfshift'
                  call twh_fieldlogical(GLB_GET,ifld,lvalue)
                  write(15,*) lvalue
                  write(15,'(a6)') 'ljfrho'
                  call twh_fielddp(GLB_GET,ifld,3,dvalue)
                  write(15,*) dvalue
                  write(15,'(a9)') 'ljfntypes'
                  call twh_fieldint(GLB_GET,ifld,1,ivalue)
                  write(15,*) ivalue
                  do itype = 1,ivalue
                     call twh_fieldnbint(GLB_GET,ifld,itype,1,ntii)
                     write(15,'(a7)') 'ljfname'
                     write(15,'(a1,a10,a1)') char(39)
     &                    ,twh_get_nbname(ntii),char(39)
                     write(15,'(a6)') 'ljfsig'
                     call twh_fieldnbparam(GLB_GET,ifld,ntii,4,fnbparam)
                     write(15,*) fnbparam
                     write(15,'(a6)') 'ljfeps'
                     call twh_fieldnbparam(GLB_GET,ifld,ntii,5,fnbparam)
                     write(15,*) fnbparam
                  enddo
               elseif ( fieldstring .eq. 'Steele Wall' ) then
                  write(15,'(a10)') 'steele box'
                  call twh_fieldint(GLB_GET,ifld,1,ivalue)
                  write(15,*) ivalue
                  write(15,'(a10)') 'steele xyz'
                  call twh_fieldchar(GLB_GET,ifld,1,fchar)
                  write(15,'(3a1)') char(39),fchar,char(39)
                  write(15,'(a14)') 'steele surface'
                  call twh_fielddp(GLB_GET,ifld,1,dvalue)
                  write(15,*) dvalue
                  write(15,'(a10)') 'steele dir'
                  call twh_fieldint(GLB_GET,ifld,2,ivalue)
                  write(15,*) ivalue
                  write(15,'(a13)') 'steele cutoff'
                  call twh_fielddp(GLB_GET,ifld,2,dvalue)
                  write(15,*) dvalue
                  write(15,'(a12)') 'steele shift'
                  call twh_fieldlogical(GLB_GET,ifld,lvalue)
                  write(15,*) lvalue
                  write(15,'(a12)') 'steele delta'
                  call twh_fielddp(GLB_GET,ifld,3,dvalue)
                  write(15,*) dvalue
                  write(15,'(a12)') 'steele rho_s'
                  call twh_fielddp(GLB_GET,ifld,4,dvalue)
                  write(15,*) dvalue
                  write(15,'(a12)') 'steele ntype'
                  call twh_fieldint(GLB_GET,ifld,1,ivalue)
                  write(15,*) ivalue
                  do itype = 1,ivalue
                     call twh_fieldnbint(GLB_GET,ifld,itype,1,ntii)
                     write(15,'(a11)') 'steele name'
                     write(15,'(a1,a10,a1)') char(39)
     &                    ,twh_get_nbname(ntii),char(39)
                     write(15,'(a8)') 'sigma_sf'
                     call twh_fieldnbparam(GLB_GET,ifld,ntii,2,fnbparam)
                     write(15,*) fnbparam
                     write(15,'(a10)') 'epsilon_sf'
                     call twh_fieldnbparam(GLB_GET,ifld,ntii,5,fnbparam)
                     write(15,*) fnbparam
                  enddo
               elseif ( fieldstring .eq. 'Harmonic Attractor' ) then
                  write(15,'(a6)') 'hafbox'
                  call twh_fieldint(GLB_GET,ifld,1,ivalue)
                  write(15,*) ivalue
                  write(15,'(a4)') 'hafk'
                  call twh_fielddp(GLB_GET,ifld,1,dvalue)
                  write(15,*) dvalue
                  write(15,'(a11)') 'hafentries'
                  call twh_fieldint(GLB_GET,ifld,2,ivalue)
                  write(15,*) ivalue
                  write(15,'(a9)') 'hafrefpos'
                  call twh_fieldchar(GLB_GET,ifld,1,fchar)
                  write(15,'(a1,a7,a1)') char(39),fchar,char(39)
                  if ( fchar .eq. 'Global') then
                     write(15,'(a10)') 'hafglobxyz'
                     call twh_fielddp(GLB_GET,ifld,2,fieldtwo)
                     call twh_fielddp(GLB_GET,ifld,3,fieldthree)
                     call twh_fielddp(GLB_GET,ifld,4,fieldfour)
                     write(15,'(3f10.6)') fieldtwo,fieldthree,fieldfour
                  elseif ( fchar .ne. 'Initial') then
                     write(6,*) 'WRITETOWHEE: unknown hafrefpos'
                     stop
                     return
                  endif
                  write(15,'(a6)') 'hafkey'
                  call twh_fieldchar(GLB_GET,ifld,2,fchar)
                  write(15,'(a1,a7,a1)') char(39),fchar,char(39)
                  call twh_fieldint(GLB_GET,ifld,2,ivalue)
                  do ientry = 1,ivalue
                     write(15,'(a8)') 'hafmolec'
                     call twh_fieldnbint(GLB_GET,ifld,ientry,1,ivalue)
                     write(15,*) ivalue
                     call twh_fieldnbint(GLB_GET,ifld,ientry,2,ntii)
                     if ( fchar .eq. 'Element') then
                        write(15,'(a10)') 'hafelement'
                        write(15,'(a1,a2,a1)') char(39)
     &                       ,twh_get_elementname(ntii),char(39)
                     elseif ( fchar .eq. 'FFtype') then
                        write(15,'(a7)') 'hafname'
                        write(15,'(a1,a10,a1)') char(39)
     &                       ,twh_get_nbname(ntii),char(39)
                     endif
                  enddo
               else
                  write(6,*) 'WRITETOWHEE: unknown fieldtype'
                  write(6,*) fieldstring
                  stop
               endif
            enddo
         endif

c        --- quantum related output
         if ( potentialstyle .eq. POT_EXTERNAL ) then
            call twh_rwexternal(lfailure,WRITE_FLAG,15)
            if ( lfailure ) return
         endif

         write(15,'(a15)') 'solvation_style'
         call twh_solvation_style(GLB_GET,cvalue)
         write(15,'(a1,a20,a1)') char(39),cvalue,char(39)
         if ( cvalue .ne. 'none' ) then
            call twh_solvation_type(GLB_GET,cvalue)
            write(15,'(a14)') 'solvation_type'
            write(15,'(a1,a20,a1)') char(39),cvalue,char(39)
         endif

         write(15,'(a5)') 'linit'
         call twh_linit(GLB_GET,lvalue)
         write(15,*) lvalue

         write(15,'(a11)') 'initboxtype'
         call twh_initboxtype(GLB_GET,initboxtype)
         write(15,'(a1,a20,a1)') char(39),initboxtype,char(39)
         if ( initboxtype .eq. 'dimensions' ) then
            ldoinitstyle = .true.
         elseif ( initboxtype .eq. 'number density' ) then
            ldoinitstyle = .true.
         elseif ( initboxtype .eq. 'unit cell' ) then
            ldoinitstyle = .false.
         else
            write(6,*) 'WRITETOWHEE: unknown initboxtype'
            write(6,*) initboxtype
            stop
            return
         endif

         if ( ldoinitstyle ) then
            write(15,'(a9)') 'initstyle'
            do ibox = 1,numboxes
               do jmolty = 1,nmolty
                  call twh_initstyle(GLB_GET,ibox,jmolty,initstyle)
                  isarray(jmolty) = initstyle
               enddo
               write(15,*) (char(39),isarray(jmolty),char(39)
     &              ,' ',jmolty=1,nmolty)
            enddo
c           --- helix_* variables are only output if the initstyle
c           --- for that molecule type is 'COM helix' or 'helix'
            do jmolty = 1,nmolty
               lfound = .false.
               do ibox = 1,numboxes
                  call twh_initstyle(GLB_GET,ibox,jmolty,initstyle)
                  if ( initstyle .eq. 'COM helix' .or.
     &                 initstyle .eq. 'helix' ) then
                     lfound = .true.
                  endif
               enddo
               if ( lfound ) then
c                 --- output the helix_* varibles
                  write(15,'(a12)') 'helix_moltyp'
                  write(15,*) jmolty
                  write(15,'(a12)') 'helix_radius'
                  call twh_helix_radius(GLB_GET,jmolty,dvalue)
                  write(15,*) dvalue
                  write(15,'(a11)') 'helix_angle'
                  call twh_helix_angle(GLB_GET,jmolty,dvalue)
                  write(15,*) dvalue
                  write(15,'(a13)') 'helix_keytype'
                  call twh_helix_keytype(GLB_GET,jmolty,tenchar)
                  write(15,'(a1,a10,a1)') char(39),tenchar,char(39)
                  write(15,'(a13)') 'helix_keyname'
                  call twh_helix_keyname(GLB_GET,jmolty,tenchar)
                  write(15,'(a1,a10,a1)') char(39),tenchar,char(39)
                  write(15,'(a12)') 'helix_conlen'
                  call twh_helix_conlen(GLB_GET,jmolty,dvalue)
                  write(15,*) dvalue
                  write(15,'(a11)') 'helix_phase'
                  call twh_helix_phase(GLB_GET,jmolty,dvalue)
                  write(15,*) dvalue
               endif
            enddo
            write(15,'(a11)') 'initlattice'
            do ibox = 1,numboxes
               do jmolty = 1,nmolty
                  call twh_initlattice(GLB_GET,ibox,jmolty,initlattice)
                  isarray(jmolty) = initlattice
               enddo
               write(15,22) (char(39),isarray(jmolty),char(39)
     &              ,jmolty=1,nmolty)
            enddo
 22         format(NTMAX(a1,a20,a1,1x))
            write(15,'(a7)') 'initmol'
            do ibox = 1,numboxes
               do jmolty = 1,nmolty
                  call twh_initmol(GLB_GET,ibox,jmolty,idummy(jmolty))
               enddo
               write(15,*) (idummy(jmolty),jmolty=1,nmolty)
            enddo
         endif
         write(15,'(a14)') 'inix iniy iniz'
         do ibox = 1,numboxes
            call twh_inixyz(GLB_GET,ibox,inix,iniy,iniz)
            write(15,*) inix,iniy,iniz
         enddo

         if ( initboxtype .eq. 'dimensions' ) then
            ldoinitstyle = .true.
            write(15,'(a7)') 'hmatrix'
            do ibox = 1,numboxes
               do idim = 1,3
                  do jdim = 1,3
                     call twh_hmatrix(GLB_GET,ibox,idim,jdim
     &                    ,dimtemp(jdim))
                  enddo
                  write(15,*) (dimtemp(jdim),jdim=1,3)
               enddo
            enddo
         elseif ( initboxtype .eq. 'number density' ) then
            ldoinitstyle = .true.
            write(15,'(a18)') 'box_number_density'
            do ibox = 1,numboxes
               call twh_box_number_density(GLB_GET,ibox,dvalue)
               write(15,*) dvalue
            enddo
         endif

         call twh_nboxpair(GLB_GET,nboxpair)
         if ( ( ensemble .eq. ENS_NVT .and. numboxes .gt. 1 ) 
     &        .or. ensemble .eq. ENS_NPT ) then
            write(15,'(a5)') 'pmvol'
            write(15,44) twh_get_pmvol()
            write(15,'(10x,a6)') 'pmvlpr'
            if ( nvmove .eq. 0 ) then
               write(15,100) 1.0d0
            else
               write(15,100) (twh_get_pmvlpr(ibox),ibox=1,nvmove)
            endif
            write(15,'(10x,a5)') 'rmvol'
            write(15,100) twh_get_rmvol(1)
            write(15,'(10x,a5)') 'tavol'
            write(15,100) twh_get_tavol()
            write(15,'(a6)') 'pmcell'
            call twh_pmcell(GLB_GET,dvalue)
            write(15,44) dvalue
            write(15,'(10x,a8)') 'pmcellpr'
            if ( nvmove .eq. 0 ) then
               write(15,100) 1.0d0
            else
               do ibox=1,nvmove
                  call twh_pmcellpr(GLB_GET,ibox,dvalue)
                  dpboxarray(ibox) = dvalue
               enddo
               write(15,100) (dpboxarray(ibox),ibox=1,nvmove)
            endif
            write(15,'(10x,a8)') 'pmcellpt'
            if ( nboxpair .eq. 0 ) then
               write(15,100) 1.0d0
            else
               do ibox = 1,nboxpair
                  call twh_pmcellpt(GLB_GET,ibox,dvalue)
                  write(15,100) dvalue
               enddo
            endif
            write(15,'(10x,a6)') 'rmcell'
            call twh_rmcell(GLB_GET,1,1,1,dvalue)
            write(15,100) dvalue
            write(15,'(10x,a6)') 'tacell'
            call twh_tacell(GLB_GET,dvalue)
            write(15,100) dvalue
         endif

         if ( numboxes .gt. 1 ) then
            write(15,'(a12)') 'pm2boxrbswap'
            call twh_pm2boxrbswap(GLB_GET,dvalue)
            write(15,44) dvalue
            write(15,'(10x,a9)') 'pm2rbswmt'
            do jmolty = 1,nmolty
               call twh_pm2rbswmt(GLB_GET,jmolty,dvalue)
               dmarray(jmolty) = dvalue
            enddo
            write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
            write(15,'(10x,a9)') 'pm2rbswpr'
            if ( nboxpair .eq. 0 ) then
               write(15,100) 1.0d0
            else
               do ibox = 1,nboxpair
                  call twh_pm2rbswpr(GLB_GET,ibox,dvalue)
                  dpboxarray(ibox) = dvalue
               enddo
               write(15,100) (dpboxarray(ibox),ibox=1,nboxpair)
            endif

            write(15,'(a12)') 'pm2boxcbswap'
            call twh_pm2boxcbswap(GLB_GET,dvalue)
            write(15,44) dvalue
            write(15,'(10x,a9)') 'pm2cbswmt'
            do jmolty = 1,nmolty
               call twh_pm2cbswmt(GLB_GET,jmolty,dvalue)
               dmarray(jmolty) = dvalue
            enddo
            write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
            write(15,'(10x,a9)') 'pm2cbswpr'
            if ( nboxpair .eq. 0 ) then
               write(15,100) 1.0d0
            else
               do ibox = 1,nboxpair
                  call twh_pm2cbswpr(GLB_GET,ibox,dvalue)
                  dpboxarray(ibox) = dvalue
               enddo
               write(15,100) (dpboxarray(ibox),ibox=1,nboxpair)
            endif
         endif

         if ( ensemble .eq. ENS_muVT ) then
            write(15,'(a11)') 'pmuvtcbswap'
            call twh_pmuvtcbswap(GLB_GET,dvalue)
            write(15,*) dvalue
            write(15,'(a9)') 'pmuvtcbmt'
            do jmolty = 1,nmolty
               call twh_pmuvtcbmt(GLB_GET,jmolty,dmarray(jmolty))
            enddo
            write(15,*) (dmarray(jmolty),jmolty=1,nmolty)
         endif

         write(15,'(a12)') 'pm1boxcbswap'
         call twh_pm1boxcbswap(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a9)') 'pm1cbswmt'
         do jmolty = 1,nmolty
            call twh_pm1cbswmt(GLB_GET,jmolty,dvalue)
            dmarray(jmolty) = dvalue
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)

         npairs = (nmolty*(nmolty-1))/2
         if ( npairs .ge. 1 ) then
            write(15,'(a15)') 'pm1boxcomswitch'
            call twh_pm1boxcomswitch(GLB_GET,dvalue)
            write(15,44) dvalue
            write(15,'(10x,a11)') 'pm1comswbox'
            do ibox = 1,numboxes
               call twh_pm1comswbox(GLB_GET,ibox,dvalue)
               dpboxarray(ibox) = dvalue
            enddo
            write(15,100) (dpboxarray(ibox),ibox=1,numboxes)
            write(15,'(10x,a12)') 'pm1comswpair'
            do ipair = 1,npairs
               call twh_pm1comswpair(GLB_GET,ipair,dvalue)
               dpairarray(ipair) = dvalue
            enddo
            write(15,100) (dpairarray(ipair),ipair=1,npairs)
         endif

         write(15,'(a6)') 'pmavb1'
         call twh_pmavb1(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a8)') 'pmavb1in'
         call twh_pmavb1in(GLB_GET,dvalue)
         write(15,100) dvalue
         write(15,'(10x,a8)') 'pmavb1mt'
         do jmolty = 1,nmolty
            call twh_pmavb1mt(GLB_GET,jmolty,dvalue)
            dmarray(jmolty) = dvalue
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         write(15,'(10x,a8)') 'pmavb1ct'
         do jmolty = 1,nmolty
            do kmolty = 1,nmolty
               call twh_pmavb1ct(GLB_GET,jmolty,kmolty,dvalue)
               dmarray(kmolty) = dvalue
            enddo
            write(15,100) (dmarray(kmolty),kmolty=1,nmolty)
         enddo
         write(15,'(10x,a7)') 'avb1rad'
         call twh_avb1rad(GLB_GET,dvalue)
         write(15,101) dvalue
         
         write(15,'(a6)') 'pmavb2'
         call twh_pmavb2(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a8)') 'pmavb2in'
         call twh_pmavb2in(GLB_GET,dvalue)
         write(15,100) dvalue
         write(15,'(10x,a8)') 'pmavb2mt'
         do jmolty = 1,nmolty
            call twh_pmavb2mt(GLB_GET,jmolty,dvalue)
            dmarray(jmolty) = dvalue
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         write(15,'(10x,a8)') 'pmavb2ct'
         do jmolty = 1,nmolty
            do kmolty = 1,nmolty
               call twh_pmavb2ct(GLB_GET,jmolty,kmolty,dvalue)
               dmarray(kmolty) = dvalue
            enddo
            write(15,100) (dmarray(kmolty),kmolty=1,nmolty)
         enddo
         write(15,'(10x,a7)') 'avb2rad'
         call twh_avb2rad(GLB_GET,dvalue)
         write(15,101) dvalue
         
         write(15,'(a6)') 'pmavb3'
         call twh_pmavb3(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a8)') 'pmavb3mt'
         do jmolty = 1,nmolty
            call twh_pmavb3mt(GLB_GET,jmolty,dvalue)
            dmarray(jmolty) = dvalue
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         write(15,'(10x,a8)') 'pmavb3ct'
         do jmolty = 1,nmolty
            do kmolty = 1,nmolty
               call twh_pmavb3ct(GLB_GET,jmolty,kmolty,dvalue)
               dmarray(kmolty) = dvalue
            enddo
            write(15,100) (dmarray(kmolty),kmolty=1,nmolty)
         enddo
         write(15,'(10x,a7)') 'avb3rad'
         call twh_avb3rad(GLB_GET,dvalue)
         write(15,101) dvalue
         
         write(15,'(a4)') 'pmcb'
         call twh_pmcb(GLB_GET,pmcb)
         write(15,44) pmcb
         write(15,'(10x,a6)') 'pmcbmt'
         do jmolty = 1,nmolty
            call twh_pmcbmt(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)

         write(15,'(10x,a5)') 'pmall'
         do jmolty = 1,nmolty
            call twh_pmall(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         
         write(15,'(a6)') 'pmback'
         call twh_pmback(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a6)') 'pmbkmt'
         do jmolty = 1,nmolty
            call twh_pmbkmt(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)

         write(15,'(a8)') 'pmcbside'
         call twh_pmcbside(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a10)') 'pmcbsidemt'
         do jmolty = 1,nmolty
            call twh_pmcbsidemt(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         
         write(15,'(a7)') 'pmpivot'
         call twh_pmpivot(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a7)') 'pmpivmt'
         do jmolty = 1,nmolty
            call twh_pmpivmt(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         
         write(15,'(a8)') 'pmconrot'
         call twh_pmconrot(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a6)') 'pmcrmt'
         do jmolty = 1,nmolty
            call twh_pmcrmt(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)

         write(15,'(a8)') 'pmcrback'
         call twh_pmcrback(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a7)') 'pmcrbmt'
         do jmolty = 1,nmolty
            call twh_pmcrbmt(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         
         write(15,'(a7)') 'pmplane'
         call twh_pmplane(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a10)') 'pmplanebox'
         do ibox = 1,numboxes
            call twh_pmplanebox(GLB_GET,ibox,dboxarray(ibox))
         enddo
         write(15,100) (dboxarray(ibox),ibox=1,numboxes)
         write(15,'(10x,a10)') 'planewidth'
         call twh_planewidth(GLB_GET,dvalue)
         write(15,100) dvalue
         
         write(15,'(a5)') 'pmrow'
         call twh_pmrow(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a8)') 'pmrowbox'
         do ibox = 1,numboxes
            call twh_pmrowbox(GLB_GET,ibox,dboxarray(ibox))
         enddo
         write(15,100) (dboxarray(ibox),ibox=1,numboxes)
         write(15,'(10x,a8)') 'rowwidth'
         call twh_rowwidth(GLB_GET,dvalue)
         write(15,100) dvalue
         
         write(15,'(a7)') 'pmtraat'
         call twh_pmtraat(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a6)') 'pmtamt'
         do jmolty = 1,nmolty
            call twh_pmtamt(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         write(15,'(10x,a6)') 'rmtraa'
         call twh_rmtraa(GLB_GET,1,1,dvalue)
         write(15,100) dvalue
         write(15,'(10x,a6)') 'tatraa'
         call twh_tatraa(GLB_GET,dvalue)
         write(15,100) dvalue

         write(15,'(a11)') 'pmcomposite'
         call twh_pmcomposite(GLB_GET,dvalue)
         write(15,100) dvalue
         write(15,'(a6)') 'pmcomt'
         do jmolty = 1,nmolty
            call twh_pmcomt(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         write(15,'(a8)') 'rmcomtra'
         call twh_rmcomtra(GLB_GET,1,1,dvalue)
         write(15,100) dvalue
         write(15,'(a8)') 'rmcomrot'
         call twh_rmcomrot(GLB_GET,1,1,dvalue)
         write(15,100) dvalue

         write(15,'(a7)') 'pmtracm'
         call twh_pmtracm(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a6)') 'pmtcmt'
         do jmolty =1,nmolty
            call twh_pmtcmt(GLB_GET,jmolty,dvalue)
            dmarray(jmolty) = dvalue
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         write(15,'(10x,a6)') 'rmtrac'
         call twh_rmtrac(GLB_GET,1,1,dvalue)
         write(15,100) dvalue
         write(15,'(10x,a6)') 'tatrac'
         call twh_tatrac(GLB_GET,dvalue)
         write(15,100) dvalue
         
         write(15,'(a8)') 'pmrotate'
         call twh_pmrotate(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a6)') 'pmromt'
         do jmolty = 1,nmolty
            call twh_pmromt(GLB_GET,jmolty,dvalue)
            dmarray(jmolty) = dvalue
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         write(15,'(10x,a5)') 'rmrot'
         call twh_rmrot(GLB_GET,1,1,dvalue)
         write(15,100) dvalue
         write(15,'(10x,a5)') 'tarot'
         call twh_tarot(GLB_GET,dvalue)
         write(15,100) dvalue

         rwfile = 15
         call twh_rwcbmc(lfailure,WRITE_FLAG,rwfile)
         if ( lfailure ) then
            write(6,*) 'WRITETOWHEE: problem with cbmc output'
            return
         endif
      else
c        --- output the molecule information for towhee_input
         write(15,'(a11)') 'input_style'
         write(15,'(a8)') 'explicit'
         write(15,'(a5)') 'nunit'
         call twh_nunit(GLB_GET,imolty,nunit)
         write(15,*) nunit
         write(15,'(a8)') 'nmaxcbmc'
         call twh_nmaxcbmc(GLB_GET,imolty,ivalue)
         write(15,*) ivalue
         write(15,'(a9)') 'lpdbnames'
         write(15,*) lprotein
         
         do iatom = 1,nunit
            write(15,'(a17)') 'unit ntype qqatom'
            write(15,*) iatom,twh_get_ntype(imolty,iatom)
     &           ,twh_get_qqatom(imolty,iatom)
            if ( lprotein ) then
c              --- include pdb stuff
               write(15,'(a27)') 'pdbname aminonum aminoshort'
               call twh_aminoshort(GLB_GET,imolty,iatom,aminoshort)
               call twh_aminonum(GLB_GET,imolty,iatom,aminonum)
               call twh_pdbname(GLB_GET,imolty,iatom,pdbname)
               write(15,*) char(39),pdbname,char(39),' '
     &              ,aminonum,' ',char(39),aminoshort,char(39)
            endif
            write(15,'(a9)') 'vibration'
            write(15,*) twh_get_invib(imolty,iatom)
            do ivib = 1,twh_get_invib(imolty,iatom)
               write(15,*) twh_get_ijvib(imolty,iatom,ivib)
     &              ,twh_get_itvib(imolty,iatom,ivib)
            enddo
            write(15,'(a7)') 'bending'
            write(15,*) twh_get_inben(imolty,iatom)
            do ibend = 1,twh_get_inben(imolty,iatom)
               call twh_orderben(GLB_GET,imolty,iatom,ibend,iorder)
               write(15,*) twh_get_ijben2(imolty,iatom,ibend)
     &              ,twh_get_ijben3(imolty,iatom,ibend)
     &              ,twh_get_itben(imolty,iatom,ibend),iorder
            enddo
            write(15,'(a7)') 'torsion'
            write(15,*) twh_get_intor(imolty,iatom)
            do itor = 1,twh_get_intor(imolty,iatom)
               write(15,*) twh_get_ijtor2(imolty,iatom,itor)
     &              ,twh_get_ijtor3(imolty,iatom,itor)
     &              ,twh_get_ijtor4(imolty,iatom,itor)
     &              ,twh_get_ittor(imolty,iatom,itor)
            enddo
            write(15,'(a11)') 'angle-angle'
            write(15,*) twh_get_inaa(imolty,iatom)
            do iaa = 1,twh_get_inaa(imolty,iatom)
               write(15,*) twh_get_ijaa0(imolty,iatom,iaa)
     &              ,twh_get_ijaa1(imolty,iatom,iaa)
     &              ,twh_get_ijaa2(imolty,iatom,iaa)
     &              ,twh_get_itaa(imolty,iatom,iaa)
            enddo
            write(15,'(a16)') 'improper torsion'
            write(15,*) twh_get_inimprop(imolty,iatom)
            do itor = 1,twh_get_inimprop(imolty,iatom)
               write(15,*) twh_get_ijimprop2(imolty,iatom,itor)
     &              ,twh_get_ijimprop3(imolty,iatom,itor)
     &              ,twh_get_ijimprop4(imolty,iatom,itor)
     &              ,twh_get_itimprop(imolty,iatom,itor)
            enddo
         enddo
      endif

 44   format(f6.4)
 100  format(10x,20(f6.4,1x))
 101  format(10x,20(f7.4,1x))
         
c     write(6.*) 'finished WRITETOWHEE'
         
      return
      end
