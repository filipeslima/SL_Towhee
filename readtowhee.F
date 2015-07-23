#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2013 Marcus G. Martin                       *
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
      subroutine twh_readtowhee(lfailure)
c     ******************************************************************
c     * reads in all of the input parameters from towhee_input         *
c     * reads in the previous configuration from towhee_initial or     *
c     * calls initconf to set up the initial conformation              *
c     * performs some simple error checking to find input errors       *
c     *                                                                *
c     * rewritten from old readinput 05-15-2002 by M.G. Martin         *
c     * last modified 07-30-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EXPON
#define FUNCTION_LEN_TRIM
#define FUNCTION_CHECK_LABEL
#define FUNCTION_GET_ELEMENTNAME
#define FUNCTION_GET_GROWVALIDCOUNT
#define FUNCTION_GET_INVIB
#define FUNCTION_GET_MASS
#define FUNCTION_GET_NBOXI
#define FUNCTION_GET_NTYPE
#define FUNCTION_GET_PMVLPR
#define FUNCTION_GET_PMVOL
#define FUNCTION_READ_FLOAT
#define FUNCTION_READ_LABELED_FLOAT
#define FUNCTION_READ_LABELED_INTEGER
#define FUNCTION_READ_LABELED_LOGICAL
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
c     --- logical scalars
      logical lpdbnames,lfound,lvalid,lreadnunit,lreadnmaxcbmc,louthist
      logical lhenry,lgetinitstyle,linit,lusecom,lvalue
      logical lreadpdbnames,lreadforcefield,lreadchargeassign
      logical lreadimpassign,lreadmatchstyle,loutdft,loutdlpoly
      logical loutlammps,lshift,ltailc,lcoulomb,lewald,ldualcutoff
      logical in_c_flag,lfieldcoords,lallow,lbondgraph
c     --- logical arrays
      logical lpdbset
      dimension lpdbset(NTMAX)
c     --- character string scalars
      character*1 hist_suffix
      character*2 ename
      character*3 fname3,fname4
      character*4 pdbname
      character*8 stepstyle
      character*10 tenchar
      character*11 fieldchar
      character*(FFNAMELEN) forcefield
      character*15 idname
      character*20 runoutput,solvation_style,solvation_type,initboxtype
      character*20 initstyle,initlattice,controlstyle,keyname
      character*20 external_code
      character*30 ftemp,filename,ctemp,ensemble_str,rng_str
      character*30 potentialstyle_str,charge_assignment
      character*30 imp_assign,match_style
      character*50 tmpfile,input_style,run_name,in_cfile
c     --- character string arrays
      character*20 satmp
      dimension satmp(NTMAX)
c     --- integer scalars
      integer imolty,jmolty,iunit,ivib,itor,imove,cpnumber,nchain,kmax
      integer ichain,ielement,npairs,ipair,nmolty,ivalue,jvalue
      integer itest,count,idim,jdim,nboxpair,nchbox,nfield
      integer icount,itmp,crtype,iutemp,idum, numunits,potentialstyle
      integer nmcount,ibox,jbox,isuccess,numtypes,isolvtype,nmtype
      integer inix,iniy,iniz,initmol,nstep,printfreq,moviefreq
      integer trmaxdispfreq,blocksize,backupfreq,volmaxdispfreq
      integer pressurefreq,pdb_output_freq,restartfreq,cperstep
      integer backlist,numboxes,nvmove,ensemble,histcalcfreq,nener
      integer histdumpfreq,n_tmmc_min,ii,n_tmmc_max,jj,maxunit,maxchain
      integer maxboxes,natoms,ifield,fieldtype,total,atomcount
      integer maxnumtor,numtor,numblock,rng_code
      integer*4 luxvalue,iseed
c     --- integer arrays
      integer idummy
      dimension idummy(NTMAX)
      integer dcrs
      dimension dcrs(NUMAX)
      integer crback
      dimension crback(40)
      integer iatmp
      dimension iatmp(NTMAX+MAXBOX)
c     --- double precision scalars
      double precision pressure,pmax,pmmin,volume,length,localmass
      double precision dvalue,pmcb,pmtest,vequiv,scalecut,scalehilog
      double precision moveprob,eefreftemp,temperature,beta,kalp,calp
      double precision box_number_density,pmcell,masstotal,scalelolog
      double precision minboxdim,chempot,uvtfactor,pmtype
      double precision cm_ii_one,cm_ii_zero,cm_ii_negone
      double precision cm_jj_one,cm_jj_zero,cm_jj_negone
      double precision xtemp,ytemp,ztemp
c     --- double precision arrays
      double precision fatmp
      dimension fatmp(3)
      double precision dtarray
      dimension dtarray(MAXPBOX)
      double precision dmarray
      dimension dmarray(NTMAX)
      double precision dmaltarray
      dimension dmaltarray(NTMAX)
      double precision dboxarray
      dimension dboxarray(MAXBOX)
      double precision prob_pair_array
      dimension prob_pair_array((NTMAX*(NTMAX-1))/2)

c     --- set subroutine idname
      idname = 'readtowhee'

c     --- optional read in the random number generator
      if ( twh_check_label(lfailure,4,'random_number_generator') ) then
         call twh_read_labeled_string(lfailure,4
     &        ,'random_number_generator',.true.,idname
     &        ,rng_str)
         if ( rng_str .eq. 'RANLUX' ) then
            rng_code = RNG_RANLUX

c           --- optional get the random number generator luxury level
            if ( twh_check_label(lfailure,4,'random_luxlevel') ) then
               luxvalue = twh_read_labeled_integer(lfailure,4
     &              ,'random_luxlevel',.true.,idname)
               if ( lfailure ) return
c              --- error check
               if ( luxvalue .lt. 0 .or. luxvalue .gt. 4 ) then
                  write(6,*) 'READINPUT: random_luxlevel must be'
     &                 ,' between 0 and 4 inclusive'
                  lfailure = .true.
                  return
               endif
               call twh_random_luxlevel(GLB_SET,luxvalue)
            endif
         elseif ( rng_str .eq. 'DX-1597-2-7' ) then
            rng_code = RNG_DX_1597_2_7
         elseif ( rng_str .eq. 'KISS99' ) then
            rng_code = RNG_KISS99
         elseif ( rng_str .eq. 'MRG32k3a' ) then
            rng_code = RNG_MRG32K3A
         elseif ( rng_str .eq. 'debug' ) then
            rng_code = RNG_DEBUG
         else
            write(6,*) 'READTOWHEE: unknown random_number_generator'
            write(6,'(a)') rng_str
            lfailure = .true.
            return
         endif
c        --- store the random number generator code
         call twh_random_code(GLB_SET,rng_code)
      endif

c     --- optional read in the random number seed
      if (twh_check_label(lfailure,4,'random_seed')) then
         call twh_get_random_seed(lfailure,4, iseed, .false., idname)
         if ( lfailure ) return
         if ( iseed .lt. 0 ) then
            write(6,*) 'READINPUT: random_seed must be postitive'
            lfailure = .true.
            return
         endif
c        --- set the global random_seed
         call twh_random_seed(GLB_SET,iseed)
      endif

c      --- optional read random_allow_restart
      if ( twh_check_label(lfailure,4,'random_allow_restart') ) then
         lallow = twh_read_labeled_logical(lfailure,4
     &        ,'random_allow_restart',.true.,idname)
         if ( lfailure ) return
         call twh_random_allow_restart(GLB_SET,lallow)
      else
c        --- use the default value
         call twh_random_allow_restart(GLB_GET,lallow)
      endif
c     --- see if we are forcing initialiation here
      if ( .not. lallow ) then
c        --- initialize random number generator now
         call twh_init_random_single(lfailure)
         if ( lfailure ) return
      endif

c     --- read in ensemble string
      call twh_read_labeled_string(lfailure,4,'ensemble',.true.,idname
     &     ,ensemble_str)
      if ( lfailure ) return
c     --- error check ensemble
      if (ensemble_str .eq. 'nvt') then
         ensemble = ENS_NVT
      elseif (ensemble_str .eq. 'npt') then
         ensemble = ENS_NPT
      elseif (ensemble_str .eq. 'uvt') then
         ensemble = ENS_muVT
      elseif ( ensemble_str .eq. 'pseudo-bubble' ) then
         ensemble = ENS_BUBBLE
      else
         write(6,*) 'READTOWHEE: invalid ensemble: ', ensemble_str
         lfailure = .true.
         return
      endif
c     --- set the ensemble
      call twh_ensemble(GLB_SET,ensemble)

c     --- read in temperature 
      temperature = twh_read_labeled_float(lfailure,4,'temperature'
     &     ,.true., idname)
      if ( lfailure ) return
c     --- Note: setting temperature automatically sets beta as well
      call twh_temperature(GLB_SET,temperature)
c     --- error check temperature
      if ( temperature .le. 0.0d0 ) then
         write(6,*) 'READTOWHEE: temperature must be positive'
         lfailure = .true.
         return
      endif

      if ( ensemble .eq. ENS_NPT ) then
c        --- read in external pressure
         pressure = twh_read_labeled_float(lfailure,4,'pressure',
     &     .true., idname)
         if ( lfailure ) return
c        --- error check pressure
         if ( pressure .le. 0.0d0 ) then
            write(6,*) 'READTOWHEE: pressure must be positive'
            lfailure = .true.
            return
         endif
c        --- convert pressure from kPa to simulation units
         dvalue = pressure*CON_KPATOSIM
         call twh_express(GLB_SET,dvalue)
      endif

c     --- tmmc_flag = .true. if doing GC-TMMC
      call twh_tmmc_flag(GLB_SET,.false.)
      if ( ensemble .eq. ENS_muVT ) then
         if (twh_check_label(lfailure,4,'tmmc_flag')) then
            lvalue = twh_read_labeled_logical(lfailure,4, 'tmmc_flag'
     &           ,.true., idname)
            if ( lfailure ) return
            call twh_tmmc_flag(GLB_SET,lvalue)
         endif
         call twh_tmmc_flag(GLB_GET,lvalue)
         if ((ensemble.ne.ENS_muVT).and. lvalue ) then
            write(6,*) 'READTOWHEE: can only do tmmc in uVT Ens.'
            lfailure = .true.
            return
         endif
      endif

c     --- read in number of molecule types
      nmolty = twh_read_labeled_integer(lfailure,4,'nmolty',.true.,
     &    idname)
      if ( lfailure ) return
      call twh_nmolty(GLB_SET,nmolty)
c     --- error check nmolty
      if ( nmolty .lt. 1 ) then
         write(6,*) 'READTOWHEE: nmolty must be at least 1'
         lfailure = .true.
         return
      endif
      if ( nmolty .gt. NTMAX ) then
         write(6,*) 'READTOWHEE: nmolty exceeds ntmax'
         write(6,*) 'set a higher value of ntmax in preproc.h'
         lfailure = .true.
         return
      endif
c     --- allocate arrays based on nmolty
      call twh_allocate_maxmolty(nmolty)

c     --- Error check number of molecule types for currently
c     --- implemented TMMC variant (for now, pure fluids only).
      call twh_tmmc_flag(GLB_GET,lvalue)
      if ( lvalue ) then
         if (nmolty.GT.1) then
            write(6,*) 'READTOWHEE: tmmc only works for '
            write(6,*) 'pure fluids now '
            lfailure = .true.
            return
         endif
      endif

c     --- read in number of each type of molecule
      call twh_read_labeled_integer_array(lfailure,4, 'nmolectyp',
     &    nmolty, .true., idname, idummy)
      if ( lfailure ) return

c     --- set the total number of molecules
      nchain = 0
      do imolty = 1, nmolty
c        --- store the nmolectyp variable
         call twh_nmolectyp(GLB_SET,imolty,idummy(imolty))
c        --- increment nchain
         nchain = nchain + idummy(imolty)
      enddo
      call twh_nchain(GLB_SET,nchain)
c     --- allocate arrays based on maxchain
      maxchain = nchain + 1
      call twh_allocate_maxchain(maxchain)
c     --- allocate based on maxchain and nmolty
      call twh_allocate_maxmolty_maxchain(nmolty,maxchain)

c     --- error check nmolectyp
      lhenry = .false.
      do imolty = 1,nmolty
         call twh_nmolectyp(GLB_GET,imolty,nmtype)
         if ( nmtype .lt. 0 ) then
            write(6,*) 
     &           'READTOWHEE: all values of moltyp must be at least 0'
            lfailure = .true.
            return
         endif
         if ( nmtype .eq. 0 ) lhenry = .true.
      enddo

      if ( ensemble .eq. ENS_muVT ) then
c        --- read in the chemical potential for each molecule
         call twh_read_labeled_float_array(lfailure,4, 'chempot',
     &      nmolty, .true., idname, dmarray)
         if ( lfailure ) return
c        --- set the values
         do jmolty = 1,nmolty
            dvalue = dmarray(jmolty)
            call twh_chempot(GLB_SET,jmolty,dvalue)
         enddo
      endif

c     --- read in the number of simulation boxes
      numboxes = twh_read_labeled_integer(lfailure,4,'numboxes',
     &    .true.,idname)
      if ( lfailure ) return
      call twh_numboxes(GLB_SET,numboxes)
      maxboxes = numboxes
c     --- error check
      if ( numboxes .gt. MAXBOX ) then
         write(6,*) 'READTOWHEE: numboxes is greater than maxbox'
         write(6,*) 'set a higher value of maxbox in preproc.h'
         lfailure = .true.
         return
      elseif ( numboxes .lt. 1 ) then
         write(6,*) 'READTOWHEE: must have at least 1 simulation box'
         lfailure = .true.
         return
      endif
      if ( lhenry ) then
c        --- numboxes must be less than maxboxes when we are computing
c        --- henry's law coefficients
         maxboxes = numboxes + 1
      endif
c     --- allocate arrays based on numboxes
      call twh_allocate_numboxes(maxboxes)
c     --- allocate arrays based on nmolty and numboxes
      call twh_allocate_maxmolty_numboxes(nmolty,maxboxes)
c     --- allocate arrays based on nmolty, numboxes, and nchain
      call twh_allocate_maxmolty_numboxes_nchain(nmolty,maxboxes,nchain)

c     --- read the step style
      call twh_read_labeled_string(lfailure,4,'stepstyle',.true.,idname
     &     ,stepstyle)
      if ( lfailure ) return
      call twh_stepstyle(GLB_SET,stepstyle)
c     --- error check
      if ( ensemble .eq. ENS_muVT .and. stepstyle .ne. 'moves' ) then
         write(6,*) 'READTOWHEE: stepstyle must be moves when using the'
         write(6,*) 'grand canonical ensemble. cycles is ill defined'
         lfailure = .true. 
         return
      endif
         
      if ( stepstyle .eq. 'cycles' .or. stepstyle .eq. 'moves' ) then
c        --- read in number of MC steps
         call twh_get_nstep(lfailure,4, nstep, .true., idname)
         if ( lfailure ) return
         call twh_nstep(GLB_SET,nstep)
c        --- error check
         if ( nstep .lt. 0 ) then
            write(6,*) 'READTOWHEE: must have at least 0 MC steps'
            lfailure = .true.
            return
         endif
      else
         write(6,*) 'READTOWHEE: unknown stepstyle'
         lfailure = .true.
         return
      endif

c     --- determine the controlstyle
      if ( twh_check_label(lfailure,4,'controlstyle') ) then
         call twh_read_labeled_string(lfailure,4,'controlstyle'
     &        ,.true.,idname,controlstyle)
         if ( lfailure ) return
      else
c        --- if not specified then default to manual
         controlstyle = 'manual'
      endif
      call twh_controlstyle(GLB_SET,controlstyle)

      if ( controlstyle .eq. 'production' ) then
c        --- setting designed for production work
         printfreq = nstep/10
         call twh_printfreq(GLB_SET,printfreq)
         blocksize = nstep/5
         call twh_blocksize(GLB_SET,blocksize)
         moviefreq = nstep/20
         call twh_moviefreq(GLB_SET,moviefreq)
         backupfreq = nstep/5
         call twh_backupfreq(GLB_SET,backupfreq)
         runoutput = 'full'
         call twh_runoutput(GLB_SET,runoutput)
         pdb_output_freq = nstep/5
         call twh_pdb_output_freq(GLB_SET,pdb_output_freq)
         loutdft = .false.
         call twh_loutdft(GLB_SET,loutdft)
         loutlammps = .false.
         call twh_loutlammps(GLB_SET,loutlammps)
         loutdlpoly = .false.
         call twh_loutdlpoly(GLB_SET,loutdlpoly)
         louthist = .false.
         call twh_louthist(GLB_SET,louthist)
         pressurefreq = 10
         call twh_pressurefreq(GLB_SET,pressurefreq)
         trmaxdispfreq = nstep/10
         call twh_trmaxdispfreq(GLB_SET,trmaxdispfreq)
         volmaxdispfreq = nstep/10
         call twh_volmaxdispfreq(GLB_SET,volmaxdispfreq)

      elseif ( controlstyle .eq. 'equilibration' ) then
c        --- setting designed to rapidly equilibrate the settings
         printfreq = nstep/10
         call twh_printfreq(GLB_SET,printfreq)
         blocksize = nstep/5
         call twh_blocksize(GLB_SET,blocksize)
         moviefreq = 0
         call twh_moviefreq(GLB_SET,moviefreq)
         backupfreq = nstep/5
         call twh_backupfreq(GLB_SET,backupfreq)
         runoutput = 'full'
         call twh_runoutput(GLB_SET,runoutput)
         pdb_output_freq = nstep/5
         call twh_pdb_output_freq(GLB_SET,pdb_output_freq)
         loutdft = .false.
         call twh_loutdft(GLB_SET,loutdft)
         loutlammps = .false.
         call twh_loutlammps(GLB_SET,loutlammps)
         loutdlpoly = .false.
         call twh_loutdlpoly(GLB_SET,loutdlpoly)
         louthist = .false.
         call twh_louthist(GLB_SET,louthist)
         pressurefreq = 10
         call twh_pressurefreq(GLB_SET,pressurefreq)
         trmaxdispfreq = nstep/100
         if ( trmaxdispfreq .eq. 0 ) trmaxdispfreq = 1
         call twh_trmaxdispfreq(GLB_SET,trmaxdispfreq)
         volmaxdispfreq = nstep/100
         if ( volmaxdispfreq .eq. 0 ) volmaxdispfreq = 1 
         call twh_volmaxdispfreq(GLB_SET,volmaxdispfreq)

      elseif ( controlstyle .eq. 'manual' ) then
c        --- read steps between printing output to stdout
         if (twh_check_label(lfailure,4,'printfreq')) then
            printfreq = twh_read_labeled_integer(lfailure,4,'printfreq'
     &           ,.true.,idname)
            if ( lfailure ) return
         else
            printfreq = 0
         endif
         call twh_printfreq(GLB_SET,printfreq)
c        --- error check (if = 0, no output)
         if ( printfreq .lt. 0 ) then
            write(6,*) 'READTOWHEE: printfreq must positive'
            lfailure = .true.
            return
         endif
c        --- weight_freq and c_matrix_freq are specified in input file
c        --- run_name is also specified, otherwise it is set to noname
c        --- n_tmmc_min is optionally specified here (for breaking
c        ---- macrostate space into windows).
         call twh_weight_freq(GLB_SET,0)
         call twh_c_matrix_freq(GLB_SET,0)
         call twh_tmmc_flag(GLB_GET,lvalue)
         if ( lvalue ) then
c           --- allocate TMMC variables
            maxchain = nchain+1
            call twh_allocate_tmmc(maxchain)
c           --- initialize transition matrix Monte Carlo variables
            call twh_n_tmmc_min(GLB_SET,0)
            call twh_n_tmmc_max(GLB_SET,nchain)
c           --- read in optional tmmc variabls
            call twh_run_name(GLB_SET,'noname')
            if ( twh_check_label(lfailure,4,'n_tmmc_min') ) then
               n_tmmc_min = twh_read_labeled_integer(lfailure,4
     &              ,'n_tmmc_min',.false.,idname)
               if ( lfailure ) return
               call twh_n_tmmc_min(GLB_SET,n_tmmc_min)
            endif
            call twh_n_tmmc_min(GLB_GET,n_tmmc_min)
            call twh_n_tmmc_max(GLB_GET,n_tmmc_max)
            write(6,*) 'N_TMMC_MAX: ', n_tmmc_max
            write(6,*) 'N_TMMC_MIN: ', n_tmmc_min
            if ( twh_check_label(lfailure,4,'weight_freq') ) then
               ivalue = twh_read_labeled_integer(lfailure,4
     &              ,'weight_freq',.true.,idname)
               if ( lfailure ) return
               call twh_weight_freq(GLB_SET,ivalue)
            endif
            if ( twh_check_label(lfailure,4,'c_matrix_freq') ) then
               ivalue = twh_read_labeled_integer(lfailure,4
     &              ,'c_matrix_freq',.true.,idname)
               if ( lfailure ) return
               call twh_c_matrix_freq(GLB_SET,ivalue)
            endif
            if ( twh_check_label(lfailure,4,'run_name') ) then
               call twh_read_labeled_string(lfailure,4,'run_name',
     &              .true.,idname,run_name)
               if ( lfailure ) return
               call twh_run_name(GLB_SET,run_name)
            endif
            in_c_flag = .false.
            if (twh_check_label(lfailure,4,'in_c_flag')) then
               in_c_flag = twh_read_labeled_logical(lfailure,4
     &              , 'in_c_flag',.true., idname)
               if ( lfailure ) return
            endif
            if ( in_c_flag ) then
               if ( twh_check_label(lfailure,4,'in_cfile') ) then
                  call twh_read_labeled_string(lfailure,4,'in_cfile',
     &                 .true., idname,in_cfile)
                  if ( lfailure ) return
                  tmpfile  = in_cfile(1:twh_len_trim(in_cfile))
                  call twh_openfile(lfailure,51,tmpfile,'DEFAULT'
     &                 ,'UNKNOWN','DEFAULT')
                  if ( lfailure ) return
                  do icount = n_tmmc_min, n_tmmc_max
                     read(51,*) itmp,cm_ii_negone,cm_ii_zero,cm_ii_one
                     call twh_c_matrix(GLB_SET,icount,-1,cm_ii_negone)
                     call twh_c_matrix(GLB_SET,icount,0,cm_ii_zero)
                     call twh_c_matrix(GLB_SET,icount,1,cm_ii_one)
                  enddo
                  CLOSE(51)
               endif
            endif
c           --- stuff that was in initialize and I think never used
            call twh_n_tmmc_min(GLB_GET,n_tmmc_min)
            if ( in_c_flag ) then
               do itor = 0, n_tmmc_max
                  call twh_v_semigrand(GLB_SET,itor,0.0D0)
               enddo
               call twh_tmmc_weight(GLB_SET,n_tmmc_min,0.0D0)
               do ii = n_tmmc_min + 1, n_tmmc_max
                  jj = ii-1
                  call twh_c_matrix(GLB_GET,jj,1,cm_jj_one)
                  call twh_c_matrix(GLB_GET,ii,-1,cm_ii_negone)
                  if ( cm_jj_one .gt. 0.0 .and. 
     &                 cm_ii_negone .gt. 0.0 ) then
                     call twh_c_matrix(GLB_GET,jj,0,cm_jj_zero)
                     call twh_c_matrix(GLB_GET,jj,-1,cm_jj_negone)
                     call twh_c_matrix(GLB_GET,ii,1,cm_ii_one)
                     call twh_c_matrix(GLB_GET,ii,0,cm_ii_zero)
                     call twh_tmmc_weight(GLB_GET,jj,dvalue)
                     dvalue = dvalue + dlog( cm_jj_one / ( cm_jj_one
     &                    + cm_jj_zero + cm_jj_negone ) ) -
     &                    dlog( cm_jj_negone / ( cm_ii_one
     &                    + cm_ii_zero + cm_ii_negone ) )
                     call twh_tmmc_weight(GLB_SET,ii,dvalue)
                  else
                     call twh_tmmc_weight(GLB_GET,jj,dvalue)
                     call twh_tmmc_weight(GLB_SET,ii,dvalue)
                  endif
               enddo
            else
               do itor = n_tmmc_min, n_tmmc_max
                  call twh_c_matrix(GLB_SET,itor,-1,0.0D0)
                  call twh_c_matrix(GLB_SET,itor,0,0.0D0)
                  call twh_c_matrix(GLB_SET,itor,+1,0.0D0)
                  call twh_tmmc_weight(GLB_SET,itor,0.0D0)
                  call twh_v_semigrand(GLB_SET,itor,0.0D0)
               enddo
            endif
c           --- end of stuff that I think was never used
         endif

c        --- read in size of blocks for block averages
         if (twh_check_label(lfailure,4,'blocksize')) then
            blocksize = twh_read_labeled_integer(lfailure,4,'blocksize'
     &           ,.true.,idname)
            if ( lfailure ) return
c           --- error check (if = 0, no output)
            if ( blocksize .lt. 0 ) then
               write(6,*) 'READTOWHEE: blocksize must be postive'
               lfailure = .true.
               return
            endif
         else
            blocksize = 0
         endif
         call twh_blocksize(GLB_SET,blocksize)

c        --- read in steps between movie file output
         if (twh_check_label(lfailure,4,'moviefreq')) then
            moviefreq = twh_read_labeled_integer(lfailure,4,'moviefreq'
     &           ,.true.,idname)
            if ( lfailure ) return
         else
            moviefreq = 0
         endif
         call twh_moviefreq(GLB_SET,moviefreq)
c        --- error check (if = 0 => no movies written)
         if ( moviefreq .lt. 0 ) then
            write(6,*) 'READTOWHEE: moviefreq must be positive'
            lfailure = .true.
            return
         endif

c        --- read in steps between backup of configuration
c        --- this is an optional label
         if (twh_check_label(lfailure,4,'backupfreq')) then
            backupfreq = twh_read_labeled_integer(lfailure,4
     &           ,'backupfreq',.true.,idname)
            if ( lfailure ) return
         else
c           --- default value if not specified
            backupfreq = 0
         endif
         call twh_backupfreq(GLB_SET,backupfreq)
c        --- error check (if = 0, no output)
         if ( backupfreq .lt. 0 ) then
            write(6,*) 'READTOWHEE: backupfreq must be positive'
            lfailure = .true.
            return
         endif

c         --- read in steps between restart save of configuration
         if (twh_check_label(lfailure,4,'restartfreq')) then
            restartfreq = twh_read_labeled_integer(lfailure,4
     &           ,'restartfreq',.true.,idname)
            if ( lfailure ) return
         else
            restartfreq = 0
         endif
         call twh_restartfreq(GLB_SET,restartfreq)
c        --- error check (if = 0, no output)
         if ( backupfreq .lt. 0 ) then
            write(6,*) 'READTOWHEE: backupfreq must be positive'
            lfailure = .true.
            return
         endif

c        --- read in the runoutput variable
         call twh_read_labeled_string(lfailure,4,'runoutput',
     &        .true.,idname,runoutput)
         if ( lfailure ) return
c        --- store the runoutput variable
         call twh_runoutput(GLB_SET,runoutput)
c        --- error check runoutput
         if ( runoutput .eq. 'full' ) then
            write(6,'(a41)') 'Full output of updates and block averages'
         elseif ( runoutput .eq. 'blocks' ) then
            write(6,'(a29)') 'Output running block averages'
         elseif ( runoutput .eq. 'updates' ) then
            write(6,'(a35)') 'Output maximum displacement updates'
         elseif ( runoutput .eq. 'none' ) then
            write(6,'(a38)') 'No output of updates or block averages'
         else
            write(6,*) 'READTOWHEE: invalid runoutput'
            lfailure = .true.
            return
         endif

c        --- read pdb output style logical
         if (twh_check_label(lfailure,4,'pdb_output_freq')) then
            pdb_output_freq = twh_read_labeled_integer(lfailure,4
     &           ,'pdb_output_freq',.true.,idname)
            if ( lfailure ) return
            call twh_pdb_output_freq(GLB_SET,pdb_output_freq)
            if ( pdb_output_freq .lt. 0) then
               write(6,*) 'READTOWHEE: pdb_output_freq must be positive'
               write(6,*) 'to disable pdb output set this variable to 0'
               lfailure = .true.
               return
            endif
         else
            call twh_pdb_output_freq(GLB_SET,0)
         endif

c        --- optional read in dft output style logical
         if (twh_check_label(lfailure,4,'loutdft') ) then
            loutdft = twh_read_labeled_logical(lfailure,4,'loutdft',
     &           .true.,idname)
            if ( lfailure ) return
            call twh_loutdft(GLB_SET,loutdft)
         endif

c        --- optional read in lammps output style logical
         if ( twh_check_label(lfailure,4,'loutlammps') ) then
            loutlammps = twh_read_labeled_logical(lfailure,4
     &           ,'loutlammps', .true.,idname)
            if ( lfailure ) return
            call twh_loutlammps(GLB_SET,loutlammps)
         endif
c        --- optional read in dlpoly output style logical
         if (twh_check_label(lfailure,4,'loutdlpoly')) then
            loutdlpoly = twh_read_labeled_logical(lfailure,4
     &           ,'loutdlpoly', .true.,idname)
            if ( lfailure ) return
            call twh_loutdlpoly(GLB_SET,loutdlpoly)
         endif

         if ( ensemble .eq. ENS_muVT 
     &        .and. twh_check_label(lfailure,4,'louthist') ) then
c           --- read logical for computation of histogram
            louthist = twh_read_labeled_logical(lfailure,4,'louthist',
     &           .true.,idname)
            if ( lfailure ) return
            call twh_louthist(GLB_SET,louthist)
            if (louthist) then
c              --- setup files for histogram reweighting
               ivalue = twh_read_labeled_integer(lfailure,4
     &              ,'hist_label',.true.,idname)
               if ( lfailure ) return
               call twh_hist_label(GLB_SET,ivalue)
               write(ftemp,*) ivalue
               call twh_read_labeled_string(lfailure,4,'hist_suffix',
     &              .true.,idname,hist_suffix)
               if ( lfailure ) return
               call twh_hist_suffix(GLB_SET,hist_suffix)
               read(ftemp,*) fname3
               filename =
     &              'towhee_his'
     &              //fname3(1:twh_len_trim(fname3))
     &              //hist_suffix//'.dat'
               call twh_file_his(GLB_SET,filename)
               filename = 
     &              'towhee_nfl'//fname3(1:twh_len_trim(fname3))
     &              //hist_suffix//'.dat'
               call twh_file_flct(GLB_SET,filename)
               ivalue = twh_read_labeled_integer(lfailure,4
     &              ,'hist_nequil',.true.,idname)
               if ( lfailure ) return
               call twh_hist_nequil(GLB_SET,ivalue)
               call twh_file_flct(GLB_GET,filename)
               call twh_openfile(lfailure,50,filename,'DEFAULT'
     &              ,'UNKNOWN','DEFAULT')
               if ( lfailure ) return
               close(50)
c              --- build filenames for the P(N) distributions
               do imolty = 1,nmolty
                  write(ftemp,*) imolty
                  read(ftemp,*) fname4
c                 --- compute the length of this character string
                  filename = 
     &                 'towhee_n'//fname4(1:twh_len_trim(fname4))
     &                 //'dis'//fname3(1:twh_len_trim(fname3))
     &                 //hist_suffix//'.dat'
                  call twh_file_ndis(GLB_SET,imolty,filename)
               enddo
            endif
            if ( louthist ) then
c              --- read frequency which the histogram information is 
c              --- calculated
               histcalcfreq = twh_read_labeled_integer(lfailure,4
     &              ,'histcalcfreq',.true.,idname)
               if ( lfailure ) return
               call twh_histcalcfreq(GLB_SET,histcalcfreq)
c              --- error check histcalcfreq
               if ( histcalcfreq .le. 0) then
                  write(6,*) 'READTOWHEE: histcalcfreq must be positive'
                  lfailure = .true.
                  return
               endif 

c              --- read frequency which the histogram information is 
c              --- output to a file
               histdumpfreq = twh_read_labeled_integer(lfailure,4
     &              ,'histdumpfreq',.true.,idname)
               if ( lfailure ) return
               call twh_histdumpfreq(GLB_SET,histdumpfreq)
c              --- error check histdumpfreq
               if (histdumpfreq .le. 0) then
                  write(6,*) 'READTOWHEE: histdumpfreq must be positive'
                  lfailure = .true.
                  return
               endif
               call twh_tmmc_flag(GLB_GET,lvalue)
               if ( .NOT. lvalue ) then
c                 --- error check combination of histdumpfreq and
c                 --- histcalcfreq
                  if (((histdumpfreq/histcalcfreq)+1) .gt. NDUMPHIST)
     &                 then
                     write(6,*) 'READTOWHEE: histdumpfreq > NDUMPHIST'
                     write(6,*) 'either increase NDUMPHIST in preproc.h'
                     write(6,*) 'or decrease histdumpfreq in'
     &                    ,' towhee_input.'
                     lfailure = .true.
                     return
                  endif
               endif
            endif
         endif

c        --- read in number of steps between pressure calculation
         if (twh_check_label(lfailure,4,'pressurefreq')) then
            pressurefreq = twh_read_labeled_integer(lfailure,4
     &           ,'pressurefreq',.true.,idname)
            if ( lfailure ) return
c           --- error check
            if ( pressurefreq .lt. 0 ) then
               write(6,*) 'READTOWHEE: pressurefreq must be'
     &              ,' non-negative'
               lfailure = .true.
               return
            endif
         else
            pressurefreq = 0
         endif
         call twh_pressurefreq(GLB_SET,pressurefreq)

c        --- read in number of steps between update of translational and
c        --- rotational maximum displacements
         if (twh_check_label(lfailure,4,'trmaxdispfreq')) then
            trmaxdispfreq = twh_read_labeled_integer(lfailure,4
     &           ,'trmaxdispfreq',.true.,idname)
            if ( lfailure ) return
c           --- error check (if zero, no updates)
            if ( trmaxdispfreq .lt. 0 ) then
               write(6,*) 'READTOWHEE: trmaxdispfreq must be positive'
               lfailure = .true.
               return
            endif
         else
c           --- default to a value of zero to turn off updates
            trmaxdispfreq = 0
         endif
         call twh_trmaxdispfreq(GLB_SET,trmaxdispfreq)

c        --- read in number of steps between update of volume move
c        --- maximum displacements
         if (twh_check_label(lfailure,4,'volmaxdispfreq')) then
            volmaxdispfreq = twh_read_labeled_integer(lfailure,4
     &           ,'volmaxdispfreq',.true.,idname)
            if ( lfailure ) return
c           --- error check (if zero, no updates)
            if ( volmaxdispfreq .lt. 0 ) then
               write(6,*) 'READTOWHEE: volmaxdispfreq must be positive'
               lfailure = .true.
               return
            endif
         else
c           --- default to a value of zero to turn off updates
            volmaxdispfreq = 0
         endif
         call twh_volmaxdispfreq(GLB_SET,volmaxdispfreq)
      else
c        --- unknown controlstyle
         write(6,*) 'READTOWHEE: unknown controlstyle'
         write(6,*) controlstyle
         lfailure = .true.
         return
      endif

c     --- optional read in the number of chemical potential
c     --- measurements to take per step for all of the molecule types
      if ( twh_check_label(lfailure,4,'chempotperstep') ) then
         call twh_read_labeled_integer_array(lfailure,4
     &        , 'chempotperstep',nmolty, .true., idname, idummy)
         if ( lfailure ) return
c        --- set
         do imolty = 1,nmolty
            call twh_chempotperstep(GLB_SET,imolty,idummy(imolty))
         enddo
      endif

c     --- error check
      do imolty = 1,nmolty
         call twh_chempotperstep(GLB_GET,imolty,cperstep)
         if ( cperstep .lt. 0 ) then
            write(6,*) 'READTOWHEE: chempotperstep must be nonnegative'
            lfailure = .true.
            return
         endif

         call twh_nmolectyp(GLB_GET,imolty,nmtype)
         if ( cperstep .eq. 0 .and. nmtype .eq. 0 ) then
            write(6,*) 'READTOWHEE: No reason to include a molecule if'
     &           ,' nmolectyp and chempotperstep are both zero'
            write(6,*) 'molecule type ',imolty
            lfailure = .true.
            return
         endif
         if ( cperstep .ne. 0 .and. numboxes .eq. MAXBOX ) then
            write(6,*) 'READTOWHEE: Can only use the extra chemical'
     &           ,' potential measurement when'
            write(6,*) 'the number of boxes is less than MAXBOX'
            write(6,*) 'numboxes:',numboxes
            write(6,*) 'maxbox:',MAXBOX
            lfailure = .true.
            return
         endif
      enddo

c     --- set numblock and allocate based on the blocksize
      call twh_blocksize(GLB_GET,blocksize)
      if (blocksize .ne. 0) then
         numblock = nstep/blocksize
c        --- check for partial block
         if ( mod(nstep,blocksize) .ne. 0 ) then
            numblock = numblock + 1
         endif
         call twh_allocate_blocks(MAXBAPROP,numblock)
      endif

c     --- set up some safety features, these are now hardwired after
c     --- running some tests to find good values
      scalecut = 50.0d0
      call twh_scalecut(GLB_SET,scalecut)
      scalelolog = 10.0d0**(-scalecut)
      call twh_scalelolog(GLB_SET,scalelolog)
      scalehilog = 10.0d0**(scalecut)
      call twh_scalehilog(GLB_SET,scalehilog)
c     --- compute vequiv
      vequiv = -dlog(scalelolog)*temperature
      call twh_vequiv(GLB_SET,vequiv)

c     --- read in the potential style
      call twh_read_labeled_string(lfailure,4,'potentialstyle',.true.
     &     ,idname,potentialstyle_str)
      if ( lfailure ) return
      if (potentialstyle_str .eq. 'internal') then
         call twh_potentialstyle(GLB_SET,POT_INTERNAL)
      elseif (potentialstyle_str .eq. 'external') then
         call twh_potentialstyle(GLB_SET,POT_EXTERNAL)
      else
         write(6,*) 'unknown potentialstyle: ', potentialstyle_str
         lfailure = .true.
         return
      endif

      call twh_potentialstyle(GLB_GET,potentialstyle)
c     --- read the classical potential if required
      if ( potentialstyle .eq. POT_INTERNAL ) then
         call twh_readclassical(lfailure,numtypes)
         if ( lfailure ) return
      endif

c     --- read the external potential if required
      if ( potentialstyle .eq. POT_EXTERNAL ) then
         call twh_rwexternal(lfailure,READ_FLAG,4)
         if ( lfailure ) return
      endif

      if ( twh_check_label(lfailure,4,'solvation_style') ) then
c        --- read in the solvation_style
         call twh_read_labeled_string(lfailure,4,'solvation_style'
     &        ,.true.,idname,solvation_style)
         if ( lfailure ) return
         call twh_solvation_style(GLB_SET,solvation_style)
         if ( solvation_style .eq. 'internal' .or.
     &        solvation_style .eq. 'external' ) then
c           --- get the solvation_type
            call twh_read_labeled_string(lfailure,4,'solvation_type'
     &           ,.true.,idname,solvation_type)
            if ( lfailure ) return
            call twh_solvation_type(GLB_SET,solvation_type)
         endif
      else
c        --- retrieve the default value
         call twh_solvation_style(GLB_GET,solvation_style)
      endif
c     --- convert solvation_style to an integer code
      call twh_set_isolvtype(lfailure,solvation_style,solvation_type)
      if ( lfailure ) return

c     --- get the integer code for the solvation type
      call twh_isolvtype(GLB_GET,isolvtype)
      if ( isolvtype .eq. SOLV_NONE ) then
         write(6,*) 'No solvation model used'
      elseif ( isolvtype .eq. SOLV_SASA ) then
c        --- solvent-assessible surface area implicit water model
c        --- that has parameters for CHARMM 19 atom types 
c        --- Ferrara, Aposolakis, Caflisch Proteins 46, 24-33 (2002)
         write(6,*) 'Solvent assessible surface area (SASA) implicit '
     &        ,'water applied'
      elseif ( isolvtype .eq. SOLV_EEF1 ) then
c        --- effective energy function (EEF1) implicit water model 
c        --- that has parameters for CHARMM 19 atom types 
c        --- Lazardis and Karplus, Proteins 35,133-152 (1999)
         write(6,*) 'Effective energy function (EEF1) implicit water '
     &        ,'applied'
         eefreftemp = 298.15d0
         call twh_temperature(GLB_GET,temperature)
         if (temperature .ne. eefreftemp) then
            call twh_thermocorrect(temperature,eefreftemp)
         endif
      elseif (isolvtype .eq. SOLV_TRAMONTO ) then
c        --- classical density functional theory using TRAMANTO 
         write(6,*) 'Tramonto solvation'
#ifndef USETRAMONTO
         write(6,*) 'READCLASSICAL: you must use the --enable-tramonto'
         write(6,*) 'configuration variable to use Tramonto solvation'
         lfailure = .true.
         return
#endif
      endif    

c     --- error check use of ewald sum with implicit solvent
      if ( isolvtype .ne. SOLV_NONE ) then
         call twh_ewald_kmax(GLB_GET,1,kmax)
         call twh_ewald_kalp(GLB_GET,kalp)
         if ((kalp .ne. 0.0d0) .or. (kmax .ne. 0)) then
            write(6,*) 'READCLASSICAL: When using implicit water models'
     &           ,' both kalp and kmax must be set to zero'
            lfailure = .true.  
            if ( lfailure ) return
         endif
      endif

c     --- read initial setup information.  linit value may come from
c     --- file or from wrapper data passed from C wrapper
      call twh_get_linit(lfailure,4, linit, .true., idname)
      if ( lfailure ) return
      call twh_linit(GLB_SET,linit)

c     --- read in the initboxtype
      call twh_read_labeled_string(lfailure,4,'initboxtype',.true.
     &     ,idname,initboxtype)
      if ( lfailure ) return
      call twh_initboxtype(GLB_SET,initboxtype)
      
c     --- initialize lgetinitstyle
      lgetinitstyle = .false.
      if ( initboxtype .eq. 'dimensions' ) then
c        --- we will also need all the initstyle information
         lgetinitstyle = .true.
      elseif ( initboxtype .eq. 'number density' ) then
c        --- we will also need all the initstyle information
         lgetinitstyle = .true.
      elseif ( initboxtype .eq. 'unit cell' ) then
c        --- we do not need the initstyle information as we are
c        --- building from unit cells
         lgetinitstyle = .false.
c        --- set all of the initstyle and initlattice
         do ibox = 1,numboxes
            do imolty = 1,nmolty
               call twh_initstyle(GLB_SET,ibox,imolty,'unit cell')
               call twh_initlattice(GLB_SET,ibox,imolty,'none')
            enddo
         enddo
      else
c        --- unknown initboxtype
         write(6,*) 'READTOWHEE: unknown initboxtype'
         write(6,*) initboxtype
         lfailure = .true.
         return
      endif

      if ( lgetinitstyle ) then
         call twh_read_label(lfailure,4,'initstyle',.false.,idname)
         if ( lfailure ) return
         do ibox = 1,numboxes
            write(6,'(a14,i5)') 'initstyle Box:',ibox
            call twh_read_string_array(lfailure,4,'initstyle',nmolty
     &           ,.true.,satmp)
            if ( lfailure ) return
            do imolty=1,nmolty
               initstyle = satmp(imolty)
               call twh_initstyle(GLB_SET,ibox,imolty,initstyle)
            enddo
         enddo

c        --- read in additional information for each molecule that has 
c        --- initstyle of 'helix cbmc' in either of the boxes
         do imolty = 1,nmolty
            lfound = .false.
            do ibox = 1,numboxes
               call twh_initstyle(GLB_GET,ibox,imolty,initstyle)
               if ( initstyle .eq. 'helix cbmc' ) then
                  lfound = .true.
               endif
            enddo
            if ( lfound ) then
               itest = twh_read_labeled_integer(lfailure,4
     &              ,'helix_moltyp',.true.,idname)
               if ( lfailure ) return
               if ( itest .ne. imolty ) then
                  write(6,*) 'Mismatch of helix_moltyp'
                  write(6,*) 'Looking for type ',imolty
                  write(6,*) 'Found type ',itest
                  lfailure = .true.
                  return
               endif
               dvalue = twh_read_labeled_float(lfailure,4
     &              ,'helix_radius',.true. ,idname)
               if ( lfailure ) return
               call twh_helix_radius(GLB_SET,imolty,dvalue)
               dvalue = twh_read_labeled_float(lfailure,4
     &              ,'helix_angle',.true. ,idname)
               if ( lfailure ) return
               call twh_helix_angle(GLB_SET,imolty,dvalue)
               call twh_read_labeled_string(lfailure,4,'helix_keytype'
     &              ,.true.,idname,tenchar)
               if ( lfailure ) return
               call twh_helix_keytype(GLB_SET,imolty,tenchar)

               call twh_read_labeled_string(lfailure,4,'helix_keyname'
     &              ,.true.,idname,tenchar)
               if ( lfailure ) return
               call twh_helix_keyname(GLB_SET,imolty,tenchar)

               dvalue = twh_read_labeled_float(lfailure,4
     &              ,'helix_conlen',.true. ,idname)
               if ( lfailure ) return
               call twh_helix_conlen(GLB_SET,imolty,dvalue)

               dvalue = twh_read_labeled_float(lfailure,4
     &              ,'helix_phase',.true. ,idname)
               if ( lfailure ) return
               call twh_helix_phase(GLB_SET,imolty,dvalue)
            endif
         enddo

         call twh_read_label(lfailure,4,'initlattice',.false.,idname)
         if ( lfailure ) return
         do ibox = 1,numboxes
            call twh_read_string_array(lfailure,4,'initlattice',nmolty
     &           ,.false.,satmp)
            if ( lfailure ) return
            write(6,66) ibox,(satmp(imolty),imolty=1,nmolty)
 66         format('Box:',1x,i2,1x,'initlattice:',1x,NTMAX(a20,1x))
            do imolty=1,nmolty
               initlattice = satmp(imolty)
               call twh_initlattice(GLB_SET,ibox,imolty,initlattice)
            enddo
         enddo

c        --- read in initial number of each type of molecule in each box
         call twh_read_label(lfailure,4,'initmol',.false.,idname)
         if ( lfailure ) return
         do ibox = 1,numboxes
c           --- read from the file
            call twh_read_integer_array(lfailure,4,'initmol',nmolty
     &           ,.false.,iatmp)
            write(6,6) ibox,'initmol:'
     &           ,(iatmp(imolty),imolty=1,nmolty)
 6          format('Box:',1x,i2,1x,a8,NTMAX(i10,1x))
            do imolty=1,nmolty
               initmol = iatmp(imolty)
               call twh_initmol(GLB_SET,ibox,imolty,initmol)
               if ( initmol .lt. 0 ) then
                  write(6,*) 'READTOWHEE: initmol must be non-negative'
                  lfailure = .true.
                  return
               endif
            enddo
         enddo
      
c        --- error check initmol against nmolectyp
         do imolty = 1,nmolty
            itest = 0
            do ibox = 1,numboxes
               call twh_initmol(GLB_GET,ibox,imolty,initmol)
               itest = itest + initmol
            enddo
            call twh_nmolectyp(GLB_GET,imolty,nmtype)
            if ( ( ensemble .ne. ENS_muVT .and. 
     &           itest .ne. nmtype )
     &           .or. (ensemble .eq. ENS_muVT 
     &           .and. itest .gt. nmtype) ) then
               write(6,*) 'READTOWHEE: moltyp does not match initmol'
               write(6,*) 'imolty ',imolty
               write(6,*) 'nmolectyp',nmtype
               lfailure = .true.
               return
            endif
         enddo
      endif
      call twh_tmmc_flag(GLB_GET,lvalue)
      if ( lvalue ) then
c        --- make sure initial number of molecules
c        --- is in the range (n_tmmc_min, n_tmmc_max).
         if ( itest .lt. n_tmmc_min ) then
            write(6,*) 'init num mol. LT n_tmmc_min'
            write(6,*) 'init num mol. ', itest
            write(6,*) 'n_tmmc_min    ', n_tmmc_min
            lfailure = .true.
            return
         endif
         if ( itest .gt. n_tmmc_max ) then
            write(6,*) 'init num mol. GT n_tmmc_max'
            write(6,*) 'init num mol. ', itest
            write(6,*) 'n_tmmc_max    ', n_tmmc_max
            lfailure = .true.
            return
         endif
      endif

c     --- read initial number of molecules in each direction
      call twh_read_label(lfailure,4,'inix iniy iniz',.false.,idname)
      if ( lfailure ) return
      do ibox = 1,numboxes
         call twh_read_integer_array(lfailure,4,'inix iniy iniz'
     &        ,3,.false.,iatmp)
         if ( lfailure ) return
         inix = iatmp(1)
         iniy = iatmp(2)
         iniz = iatmp(3)
         call twh_inixyz(GLB_SET,ibox,inix,iniy,iniz)
         write(6,'(a4,1x,i2,1x,a15,3(i10,1x))') 
     &        'Box:',ibox,'inix,iniy,iniz:',inix,iniy,iniz
c        --- error check
         if ( inix .lt. 0 ) then
            write(6,*) 'READTOWHEE: inix must be non-negative'
            lfailure = .true.
            return
         endif
         if ( iniy .lt. 0 ) then
            write(6,*) 'READTOWHEE: iniy must be non-negative'
            lfailure = .true.
            return
         endif
         if ( iniz .lt. 0 ) then
            write(6,*) 'READTOWHEE: iniz must be non-negative'
            lfailure = .true.
            return
         endif
      enddo

      if ( initboxtype .eq. 'dimensions' ) then
c        --- read in initial box dimensions
         call twh_read_label(lfailure,4,'hmatrix',.false.,idname)
         if ( lfailure ) return
         do ibox = 1,numboxes
            do idim = 1,3
c              --- read box dimensions from a file
               call twh_read_float_array(lfailure,4,'hmatrix',3
     &              ,.false.,fatmp)
               if ( lfailure ) return
               do jdim=1,3
                  call twh_hmatrix(GLB_SET,ibox,idim,jdim,fatmp(jdim))
               enddo
               write(6,'(a17,i5,1x,i1,3(1x,f10.5))') 
     &              'Box idim hmatrix:',ibox,idim
     &              ,(fatmp(jdim),jdim=1,3)
            enddo
c           --- error check box dimensions
            call twh_checkhmatrix(ibox,.true.,lfailure)
            if ( lfailure ) return
         enddo
      elseif ( initboxtype .eq. 'number density' ) then
c        --- read in initial box number density
         call twh_read_label(lfailure,4,'box_number_density',.false.
     &        ,idname)
         if ( lfailure ) return
         do ibox = 1,numboxes
            box_number_density = twh_read_float(lfailure,4
     &           ,'box_number_density',.true.)
            if ( lfailure ) return
            call twh_box_number_density(GLB_SET,ibox,box_number_density)
c           --- error check
            if ( box_number_density .le. 0.0d0 ) then
               write(6,*) 'READTOWHEE: box_number_density must be'
     &              ,' positive'
               lfailure = .true.
               return
            endif
c           --- transform the initial number density into a cubic box
c           --- count the number of molecules in the box
            count = 0
            do imolty = 1,nmolty
               call twh_initmol(GLB_GET,ibox,imolty,initmol)
               count = count + initmol
            enddo
            volume = 1000.0d0*(dble(count)/box_number_density)
            length = volume**(1.0d0/3.0d0)
            do idim = 1,3
               call twh_hmatrix(GLB_SET,ibox,idim,idim,length)
            enddo
         enddo
      endif

c     --- determine the number of possible box pairs
      nboxpair = numboxes*(numboxes-1)/2
      call twh_nboxpair(GLB_SET,nboxpair)
      if ( ensemble .eq. ENS_NPT ) then
         call twh_nvmove(GLB_SET,numboxes)
      else
         call twh_nvmove(GLB_SET,nboxpair)
      endif
c     --- fetch the nvmove for later use in this routine
      call twh_nvmove(GLB_GET,nvmove)
c     --- allocate arrays based upon nvmove
      call twh_allocate_nvmove(nvmove)

c     --- set up pairbox array (used for volume and swap moves)
      itest = 0
      do ibox = 1,numboxes-1
         do jbox = ibox+1,numboxes
            itest = itest + 1
            call twh_pairbox(GLB_SET,itest,1,ibox)
            call twh_pairbox(GLB_SET,itest,2,jbox)
            write(6,'(a6,i3,1x,a8,i3,1x,i3,1x)') 'itest:',itest
     &           ,'pairbox:',ibox,jbox
         enddo
      enddo
      if ( itest .ne. nboxpair ) then
         write(6,*) 'READTOWHEE: pair box setup wrong'
         lfailure = .true.
         return
      endif

c     --- initialize pmmin
      pmmin = 0.0d0

c     A note about default values: all move probabilities (pm*) are now
c     optional, implying that they will take on the default values if
c     not specified in towhee_input.  The default values are specified
c     in initialize.F, and should give zero probailites to moves not
c     otherwise defined.  When adding new moves, please add them to
c     initialize.F

c     --- isotropic volume move
      if ( ( ensemble .eq. ENS_NVT .and. numboxes .gt. 1 ) 
     &     .or. ensemble .eq. ENS_NPT ) then

c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c        --- read in probabilities for the 3d volume move
         if (twh_check_label(lfailure,4,'pmvol')) then
            dvalue = twh_read_labeled_float(lfailure,4,'pmvol',.true.
     &           ,idname)
            if ( lfailure ) return
            call twh_pmvol(GLB_SET,dvalue)
            if ( twh_get_pmvol() .gt. pmmin ) pmmin = twh_get_pmvol()

            call twh_read_labeled_float_array(lfailure,4,'pmvlpr'
     &           ,nvmove, .true. , idname, dtarray)
            if ( lfailure ) return
            do imove = 1,nvmove
               call twh_pmvlpr(GLB_SET,imove,dtarray(imove))
            enddo
c           --- error check pmvlpr
            pmax = 0.0d0
            do ibox = 1,nvmove
               if ( twh_get_pmvlpr(ibox) .gt. pmax ) then
                  pmax = twh_get_pmvlpr(ibox)
               endif
            enddo
            if ( pmax .lt. 1.0d0 ) then
               write(6,*) 'READTOWHEE: at least one value of pmvlpr'
     &              ,' must be 1.0d0'
               lfailure = .true.
               return
            endif

c           --- read in the initial 3d volume move maximum displacements
            dvalue = twh_read_labeled_float(lfailure,4,'rmvol',.true.
     &           ,idname)
            if ( lfailure ) return
c           --- safety check
            if ( dvalue .le. 0.0d0 ) then
               write(6,*) 'READTOWHEE: rmvol must be positive'
               lfailure = .true.
               return
            endif
c           --- set all to the initial value
            call twh_rmvol(GLB_INIT,nvmove,dvalue)
            
c           --- read in the target acceptance rate for 3d volume
            dvalue = twh_read_labeled_float(lfailure,4,'tavol',.true.
     &           ,idname)
            if ( lfailure ) return
c           --- error check
            if ( dvalue .le. 0.0d0 .or. dvalue .ge. 1.0d0 ) then
               write(6,*) 'READTOWHEE: tavol must be between 0.0 ',
     &               'and 1.0'
               lfailure = .true.
               return
            endif
            call twh_tavol(GLB_SET,dvalue)
         endif
      else
c     user-friendly sanity check.  If pmvol is accidentally in
c     towhee_input, the optional-label logic will skip all other pm*,
c     not what we expect.
c     Catch this condition and complain.
         if (twh_check_label(lfailure,4,'pmvol')) then
            write(6,*) 'READTOWHEE: pmvol is permitted only for NPT; ',
     &         'or, NVT and > 1 numboxes'
            lfailure = .true.
            return
         endif
      endif

c     --- anisotropic volume move
      if ( ( ensemble .eq. ENS_NVT .and. numboxes .gt. 1 ) 
     &     .or. ensemble .eq. ENS_NPT ) then

c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c        --- read in probabilities for the unit cell adjust move
         if (twh_check_label(lfailure,4,'pmcell')) then
            pmcell = twh_read_labeled_float(lfailure,4,'pmcell',.true.
     &           ,idname)
            if ( lfailure ) return
            call twh_pmcell(GLB_SET,pmcell)
            if ( pmcell .gt. pmmin ) pmmin = pmcell

            call twh_read_labeled_float_array(lfailure,4,'pmcellpr',
     &         nvmove, .true., idname, dtarray)
            if ( lfailure ) return
c           --- error check pmcellpr
            pmax = 0.0d0
            do ibox = 1,nvmove
               dvalue = dtarray(ibox)
               call twh_pmcellpr(GLB_SET,ibox,dvalue)
               if ( dvalue .gt. pmax ) pmax = dvalue
            enddo
            if ( pmax .lt. 1.0d0 ) then
               write(6,*) 'READTOWHEE: at least one value of pmcellpr'
     &              ,' must be 1.0d0'
               lfailure = .true.
               return
            endif
            if ( ensemble .eq. ENS_NPT ) then
c              --- ignore pmcellpt
               call twh_skip_lines(lfailure,4,2)
               if ( lfailure ) return
            else
c              --- read in pmcellpt
               call twh_read_label(lfailure,4,'pmcellpt',.false.,idname)
               if ( lfailure ) return
               do ibox = 1,nboxpair
                  dvalue = twh_read_float(lfailure,4,'pmcellpt',.false.)
                  if ( lfailure ) return
                  call twh_pmcellpt(GLB_SET,ibox,dvalue)
                  call twh_pairbox(GLB_GET,ibox,1,ivalue)
                  call twh_pairbox(GLB_GET,ibox,2,jvalue)
                  write(6,9) 'pairbox:',ivalue,jvalue
     &                 ,'pmcellpt:',dvalue
 9                format(a8,i2,1x,i2,1x,a9,f9.6)
                  if ( dvalue .lt. 0.0d0 ) then
                     write(6,*) 'READTOWHEE: pmcellpt must be'
     &                    ,' non-negative'
                     lfailure = .true.
                     return
                  endif
               enddo
            endif

c           --- read in the initial unit cell adjustment max
c           --- displacements
            dvalue = twh_read_labeled_float(lfailure,4,'rmcell',.true.
     &           ,idname)
            if ( lfailure ) return
            if ( dvalue .le. 0.0d0 ) then
               write(6,*) 'READTOWHEE: rmcell must be positive'
               lfailure = .true.
               return
            endif
c           --- set all to the initial value
            call twh_rmcell(GLB_INIT,nvmove,3,3,dvalue)

c           --- read in the target acceptance rate for unit cell move
            dvalue = twh_read_labeled_float(lfailure,4,'tacell',.true.
     &           ,idname)
            if ( lfailure ) return
            call twh_tacell(GLB_SET,dvalue)
c           --- error check
            if ( dvalue .le. 0.0d0 .or. dvalue .ge. 1.0d0 ) then
               write(6,*) 'READTOWHEE: tacell must be between 0.0 ',
     &               'and 1.0'
               lfailure = .true.
               return
            endif
         endif
      endif

c     --- 2 box rotational-bias swap
      if ( numboxes .gt. 1 ) then 

c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c        --- read 2-box rotational-bias swap move
         if (twh_check_label(lfailure,4,'pm2boxrbswap')) then
            dvalue = twh_read_labeled_float(lfailure,4,'pm2boxrbswap'
     &           ,.true. ,idname)
            if ( lfailure ) return
            call twh_pm2boxrbswap(GLB_SET,dvalue)
            if ( dvalue .gt. pmmin ) pmmin = dvalue
            call twh_read_labeled_float_array(lfailure,4,'pm2rbswmt',
     &            nmolty, .true., idname, dmarray)
            if ( lfailure ) return
            do imolty = 1,nmolty
               dvalue = dmarray(imolty)
               call twh_pm2rbswmt(GLB_SET,imolty,dvalue)
            enddo
c           --- error check pm2rbswmt
            pmax = 0.0d0
            do imolty = 1,nmolty
               call twh_pm2rbswmt(GLB_GET,imolty,dvalue)
               if ( dvalue .gt. pmax ) pmax = dvalue
            enddo
            if ( pmax .lt. 1.0d0 ) then
               write(6,*) 'READTOWHEE: at least one value of pm2rbswmt'
     &              ,' must be 1.0d0'
               lfailure = .true.
               return
            endif

            call twh_read_labeled_float_array(lfailure,4,'pm2rbswpr',
     &            nboxpair, .true., idname, dtarray)
            if ( lfailure ) return
            do imove = 1,nboxpair
               dvalue = dtarray(imove)
               call twh_pm2rbswpr(GLB_SET,imove,dvalue)
            enddo
c           --- error check pm2rbswpr
            call twh_pm2rbswpr(GLB_GET,nboxpair,dvalue)
            if ( dvalue .lt. 1.0d0 ) then
               write(6,*) 'READTOWHEE: the final value of pm2rbswpr'
     &              ,' must be 1.0'
               lfailure = .true.
               return
            endif
         endif
      endif

c     --- 2 box configurational-bias swap
      if ( numboxes .gt. 1 ) then

c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c        --- read 2boxcbswap info
         if (twh_check_label(lfailure,4,'pm2boxcbswap')) then
            dvalue = twh_read_labeled_float(lfailure,4,'pm2boxcbswap',
     &         .true., idname)
            if ( lfailure ) return
            call twh_pm2boxcbswap(GLB_SET,dvalue)

            if ( dvalue .gt. pmmin ) pmmin = dvalue
            call twh_read_labeled_float_array(lfailure,4,'pm2cbswmt',
     &            nmolty, .true., idname, dmarray)
            if ( lfailure ) return
c           --- error check and set pm2cbswmt
            pmax = 0.0d0
            do imolty = 1,nmolty
               dvalue = dmarray(imolty)
               call twh_pm2cbswmt(GLB_SET,imolty,dvalue)
               if ( dvalue .gt. pmax ) pmax = dvalue
            enddo
            if ( pmax .lt. 1.0d0 ) then
               write(6,*) 'READTOWHEE: at least one value of pm2cbswmt'
     &              ,' must be 1.0d0'
               lfailure = .true.
               return
            endif

            call twh_read_labeled_float_array(lfailure,4,'pm2cbswpr',
     &            nboxpair, .true., idname, dtarray)
            if ( lfailure ) return
            do imove = 1,nboxpair
               dvalue = dtarray(imove)
               call twh_pm2cbswpr(GLB_SET,imove,dvalue)
            enddo
c           --- error check pm2cbswpr
            call twh_pm2cbswpr(GLB_GET,nboxpair,dvalue)
            if ( dvalue .lt. 1.0d0 ) then
               write(6,*) 'READTOWHEE: the final value of pm2cbswpr'
     &              ,' must be 1.0'
               lfailure = .true.
               return
            endif
         endif
      endif

c     --- grand canonical insertion/deletion
      if ( ensemble .eq. ENS_muVT ) then
c        --- read in grand canonical swap move info
         if (twh_check_label(lfailure,4,'pmuvtcbswap')) then
            dvalue = twh_read_labeled_float(lfailure,4,'pmuvtcbswap'
     &           ,.true., idname)
            if ( lfailure ) return
            call twh_pmuvtcbswap(GLB_SET,dvalue)

            if ( dvalue .gt. pmmin ) pmmin = dvalue
            call twh_read_labeled_float_array(lfailure,4,'pmuvtcbmt'
     &           ,nmolty ,.true.,idname,dmarray)
            if ( lfailure ) return
c           --- set and error check pmuvtcbmt
            pmax = 0.0d0
            do imolty = 1,nmolty
               dvalue = dmarray(imolty)
               call twh_pmuvtcbmt(GLB_SET,imolty,dvalue)
               if ( dvalue .gt. pmax ) pmax = dvalue
            enddo
            if ( pmax .lt. 1.0d0 ) then
               write(6,*) 'READTOWHEE: at least one value of pmuvtcbmt'
     &              ,' must be 1.0d0'
               lfailure = .true.
               return
            endif
         endif
      endif

c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     --- read intra-box cbswap (single box) info
      if (twh_check_label(lfailure,4,'pm1boxcbswap')) then
         dvalue = twh_read_labeled_float(lfailure,4,'pm1boxcbswap'
     &        ,.true.,idname)
         if ( lfailure ) return
         call twh_pm1boxcbswap(GLB_SET,dvalue)
         if ( dvalue .gt. pmmin ) pmmin = dvalue
         call twh_read_labeled_float_array(lfailure,4,'pm1cbswmt',
     &       nmolty, .true., idname,dmarray)
         if ( lfailure ) return

c        --- save and error check pm1cbswmt
         pmax = 0.0d0
         do imolty = 1,nmolty
            dvalue = dmarray(imolty)
            call twh_pm1cbswmt(GLB_SET,imolty,dvalue)
            if ( dvalue .gt. pmax ) pmax = dvalue
         enddo
         if ( pmax .lt. 1.0d0 ) then
            write(6,*) 'READTOWHEE: at least one value of pm1cbswmt'
     &           ,' must be 1.0d0'
            lfailure = .true.
            return
         endif
      endif


c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      npairs = nmolty*(nmolty-1)/2
      if ( npairs .gt. 0 ) then
c        --- read intra-box center-of-mass switch (single box) info
         if (twh_check_label(lfailure,4,'pm1boxcomswitch')) then
            moveprob = twh_read_labeled_float(lfailure,4
     &           ,'pm1boxcomswitch',.true. ,idname)
            if ( lfailure ) return
            call twh_pm1boxcomswitch(GLB_SET,moveprob)
            if ( moveprob .gt. pmmin ) pmmin = moveprob
            
            call twh_read_labeled_float_array(lfailure,4,'pm1comswbox',
     &           numboxes, .true., idname,dboxarray)
            if ( lfailure ) return
c           --- error check box probabilities
            pmax = 0.0d0
            do ibox = 1,numboxes
               if ( dboxarray(ibox) .gt. pmax ) then
                  pmax = dboxarray(ibox)
               endif
            enddo
            if ( pmax .lt. 1.0d0 ) then
               write(6,*) 'READTOWHEE: at least one value of'
     &              ,' pm1comswbox '
     &              ,'must be 1.0d0'
               lfailure = .true.
               return
            endif
c           --- store the values in pm1comswbox
            do ibox = 1,numboxes
               call twh_pm1comswbox(GLB_SET,ibox,dboxarray(ibox))
            enddo
            
            call twh_read_labeled_float_array(lfailure,4,'pm1comswpair',
     &           npairs, .true., idname,prob_pair_array)
            if ( lfailure ) return
c           --- error check box probabilities
            pmax = 0.0d0
            do ipair = 1,npairs
               if ( prob_pair_array(ipair) .gt. pmax ) then
                  pmax = prob_pair_array(ipair)
               endif
            enddo
            if ( pmax .lt. 1.0d0 ) then
               write(6,*) 'READTOWHEE: at least one value of'
     &              ,' pm1comswpair '
     &              ,'must be 1.0d0'
               lfailure = .true.
               return
            endif
c           --- store the values in pm1comswpair
            do ipair = 1,npairs
               call twh_pm1comswpair(GLB_SET,ipair
     &              ,prob_pair_array(ipair))
            enddo
         endif
      endif

c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     --- read aggregation volume bias move (AVB-1) info
      if (twh_check_label(lfailure,4,'pmavb1')) then
         dvalue = twh_read_labeled_float(lfailure,4,'pmavb1',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_pmavb1(GLB_SET,dvalue)
         if ( dvalue .gt. pmmin ) pmmin = dvalue
         dvalue = twh_read_labeled_float(lfailure,4,'pmavb1in',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_pmavb1in(GLB_SET,dvalue)
c        --- error check pmavb1in
         if ( dvalue .le. 0.0d0 .or. dvalue .ge. 1.0d0 ) then
            write(6,*) 'READTOWHEE: pmavb1in must be strictly between'
     &           ,' 0.0 and 1.0'
            lfailure = .true.
            return
         endif

         call twh_read_labeled_float_array(lfailure,4, 'pmavb1mt',
     &        nmolty, .true., idname, dmarray)
         if ( lfailure ) return
c        --- set and error check pmavb1mt
         pmax = 0.0d0
         do imolty = 1,nmolty
            dvalue = dmarray(imolty)
            call twh_pmavb1mt(GLB_SET,imolty,dvalue)
            if ( dvalue .gt. pmax ) pmax = dvalue
         enddo
         if ( pmax .lt. 1.0d0 ) then
            write(6,*) 'READTOWHEE: at least one value of pmavb1mt must'
     &           ,' be 1.0d0'
            lfailure = .true.
            return
         endif

         call twh_read_label(lfailure,4,'pmavb1ct',.false.,idname)
         if ( lfailure ) return
         do imolty = 1,nmolty
c        look at pmavbXct
            call twh_read_float_array(lfailure,4, 'pmavb1ct',nmolty
     &           ,.false.,dmarray)
            if ( lfailure ) return
c           --- error check pmavb1ct
            pmax = 0.0d0
            do jmolty = 1,nmolty
               dvalue = dmarray(jmolty)
               call twh_pmavb1ct(GLB_SET,imolty,jmolty,dvalue)
               if ( dvalue .gt. pmax ) pmax = dvalue
            enddo
            write(6,14) 'moltyp:',imolty
     &           ,'pmavb1ct:',(dmarray(jmolty),jmolty=1,nmolty)
 14         format(a7,i3,1x,a9,NTMAX(f9.6,1x))
            if ( pmax .lt. 1.0d0 ) then
               write(6,*) 'READTOWHEE: for each moltyp at least one'
     &              ,' value of pmavb1ct must be 1.0d0'
               lfailure = .true.
               return
            endif
         enddo

         dvalue = twh_read_labeled_float(lfailure,4,'avb1rad',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_avb1rad(GLB_SET,dvalue)
c        --- error check avb1rad
         if ( dvalue .le. 0.0d0 ) then
            write(6,*) 'READTOWHEE: avb1rad must be greater than 0.0'
            lfailure = .true.
            return
         endif
      endif


c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     --- read aggregation volume bias move (AVB-2) info
      if (twh_check_label(lfailure,4,'pmavb2')) then
         dvalue = twh_read_labeled_float(lfailure,4,'pmavb2',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_pmavb2(GLB_SET,dvalue)
         if ( dvalue .gt. pmmin ) pmmin = dvalue
         dvalue = twh_read_labeled_float(lfailure,4,'pmavb2in',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_pmavb2in(GLB_SET,dvalue)
c        --- error check pmavb2in
         if ( dvalue .le. 0.0d0 .or. dvalue .ge. 1.0d0 ) then
            write(6,*) 'READTOWHEE: pmavb2in must be between 0.0 ',
     &            'and 1.0'
            lfailure = .true.
            return
         endif

         call twh_read_labeled_float_array(lfailure,4, 'pmavb2mt',
     &         nmolty, .true., idname, dmarray)
         if ( lfailure ) return
c        --- set and error check pmavb2mt
         pmax = 0.0d0
         do imolty = 1,nmolty
            dvalue = dmarray(imolty)
            call twh_pmavb2mt(GLB_SET,imolty,dvalue)
            if ( dvalue .gt. pmax ) pmax = dvalue
         enddo
         if ( pmax .lt. 1.0d0 ) then
            write(6,*) 'READTOWHEE: at least one value of pmavb2mt must'
     &           ,' be 1.0d0'
            lfailure = .true.
            return
         endif

         call twh_read_label(lfailure,4,'pmavb2ct',.false.,idname)
         if ( lfailure ) return
         do imolty = 1,nmolty
            call twh_read_float_array(lfailure,4, 'pmavb2ct',nmolty
     &           ,.false.,dmarray)
            if ( lfailure ) return
c           --- set and error check pmavb2ct
            pmax = 0.0d0
            do jmolty = 1,nmolty
               dvalue = dmarray(jmolty)
               call twh_pmavb2ct(GLB_SET,imolty,jmolty,dvalue)
               if ( dvalue .gt. pmax ) pmax = dvalue
            enddo
            write(6,16) 'moltyp:',imolty
     &           ,'pmavb2ct:',(dmarray(jmolty),jmolty=1,nmolty)
 16         format(a7,i3,1x,a9,NTMAX(f9.6,1x))
            if ( pmax .lt. 1.0d0 ) then
               write(6,*) 'READTOWHEE: for each moltyp at least one'
     &              ,' value of pmavb2ct must be 1.0d0'
               lfailure = .true.
               return
            endif
         enddo

         dvalue = twh_read_labeled_float(lfailure,4,'avb2rad',.true.
     &        , idname)
         if ( lfailure ) return
         call twh_avb2rad(GLB_SET,dvalue)
c        --- error check avb2rad
         if ( dvalue .le. 0.0d0 ) then
            write(6,*) 'READTOWHEE: avb2rad must be greater than 0.0'
            lfailure = .true.
            return
         endif
      endif


c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     --- read aggregation volume bias move (AVB-3) info
      if (twh_check_label(lfailure,4,'pmavb3')) then
         dvalue = twh_read_labeled_float(lfailure,4,'pmavb3',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_pmavb3(GLB_SET,dvalue)
         if ( dvalue .gt. pmmin ) pmmin = dvalue

         call twh_read_labeled_float_array(lfailure,4, 'pmavb3mt',
     &         nmolty, .true., idname, dmarray)
         if ( lfailure ) return
c        --- store and error check pmavb3mt
         pmax = 0.0d0
         do imolty = 1,nmolty
            dvalue = dmarray(imolty)
            call twh_pmavb3mt(GLB_SET,imolty,dvalue)
            if ( dvalue .gt. pmax ) pmax = dvalue
         enddo
         if ( pmax .lt. 1.0d0 ) then
            write(6,*) 'READTOWHEE: at least one value of pmavb3mt must'
     &           ,' be 1.0d0'
            lfailure = .true.
            return
         endif

         call twh_read_label(lfailure,4,'pmavb3ct',.false.,idname)
         if ( lfailure ) return
         do imolty = 1,nmolty
            call twh_read_float_array(lfailure,4, 'pmavb3ct',
     &         nmolty, .false., dmarray)
            if ( lfailure ) return
c           --- fill in and error check pmavb3ct
            pmax = 0.0d0
            do jmolty = 1,nmolty
               dvalue =  dmarray(jmolty)
               call twh_pmavb3ct(GLB_SET,imolty,jmolty,dvalue)
               if ( dvalue .gt. pmax ) pmax = dvalue
            enddo
            write(6,18) 'moltyp:',imolty
     &           ,'pmavb3ct:',(dmarray(jmolty),jmolty=1,nmolty)
 18         format(a7,i3,1x,a9,NTMAX(f9.6,1x))
            if ( pmax .lt. 1.0d0 ) then
               write(6,*) 'READTOWHEE: for each moltyp at least one'
     &              ,' value of pmavb3ct must be 1.0d0'
               lfailure = .true.
               return
            endif
         enddo

         dvalue = twh_read_labeled_float(lfailure,4,'avb3rad',.true.
     &        , idname)
         if ( lfailure ) return
         call twh_avb3rad(GLB_SET,dvalue)
c        --- error check avb3rad
         if ( dvalue .le. 0.0d0 ) then
            write(6,*) 'READTOWHEE: avb3rad must be greater than 0.0'
            lfailure = .true.
            return
         endif
      endif


c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     --- read regular cbmc info
      if (twh_check_label(lfailure,4,'pmcb')) then
         pmcb = twh_read_labeled_float(lfailure,4,'pmcb',.true.,idname)
         if ( lfailure ) return
         call twh_pmcb(GLB_SET,pmcb)
         if ( pmcb .gt. pmmin ) pmmin = pmcb

         call twh_read_labeled_float_array(lfailure,4,'pmcbmt',nmolty
     &        ,.true.,idname,dmarray)
         if ( lfailure ) return
c        --- set and error check pmcbmt
         pmax = 0.0d0
         do imolty = 1,nmolty
            pmtest = dmarray(imolty)
            call twh_pmcbmt(GLB_SET,imolty,pmtest)
            if ( pmtest .gt. pmax ) pmax = pmtest
         enddo
         if ( pmax .lt. 1.0d0 ) then
            write(6,*) 'READTOWHEE: at least one value of pmcbmt must'
     &           ,' be 1.0d0'
            lfailure = .true.
            return
         endif

         call twh_read_labeled_float_array(lfailure,4,'pmall',nmolty
     &        ,.true., idname,dmarray)
         if ( lfailure ) return
c        --- set and error check pmall
         do imolty = 1,nmolty
            pmtest = dmarray(imolty)
            call twh_pmall(GLB_SET,imolty,pmtest)
            if ( pmtest .lt. 0.0d0 .or. pmtest .gt. 1.0d0)
     &           then
               write(6,*) 'READTOWHEE: all values of pmall must be'
     &              ,' between 0.0 and 1.0 (inclusive)'
               lfailure = .true.
               return
            endif
         enddo
      endif


c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     --- read backbone regrowth cbmc info
      if (twh_check_label(lfailure,4,'pmback')) then
         pmtest = twh_read_labeled_float(lfailure,4,'pmback',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_pmback(GLB_SET,pmtest)
         if ( pmtest .gt. pmmin ) pmmin = pmtest
         call twh_read_labeled_float_array(lfailure,4,'pmbkmt',nmolty
     &        ,.true.,idname,dmarray)
         if ( lfailure ) return
c        --- set and error check pmbkmt
         pmax = 0.0d0
         do imolty = 1,nmolty
            pmtest = dmarray(imolty)
            call twh_pmbkmt(GLB_SET,imolty,pmtest)
            if ( pmtest .gt. pmax ) pmax = pmtest
         enddo
         if ( pmax .lt. 1.0d0 ) then
            write(6,*) 'READTOWHEE: at least one value of pmbkmt must'
     &           ,' be 1.0d0'
            lfailure = .true.
            return
         endif
      endif

      
c     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     --- read cbmc sidechain regrowth information
      if (twh_check_label(lfailure,4,'pmcbside')) then
         dvalue = twh_read_labeled_float(lfailure,4,'pmcbside',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_pmcbside(GLB_SET,dvalue)
         if ( dvalue .gt. pmmin ) pmmin = dvalue

         call twh_read_labeled_float_array(lfailure,4,'pmcbsidemt'
     &        ,nmolty,.true.,idname,dmarray)
         if ( lfailure ) return
         do imolty = 1,nmolty
            call twh_pmcbsidemt(GLB_SET,imolty,dmarray(imolty))
         enddo
c        --- error check pmcbsidemt
         pmax = 0.0d0
         do imolty = 1,nmolty
            call twh_pmcbsidemt(GLB_GET,imolty,pmtest)
            if ( pmtest .gt. pmax ) pmax = pmtest
         enddo
         if ( pmax .lt. 1.0d0 ) then
            write(6,*) 'READTOWHEE: at least one value of pmcbsidemt'
     &           ,' must be 1.0d0'
            lfailure = .true.
            return
         endif
      endif


c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     --- read pivot move info
      if (twh_check_label(lfailure,4,'pmpivot')) then
         dvalue = twh_read_labeled_float(lfailure,4,'pmpivot',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_pmpivot(GLB_SET,dvalue)
         if ( dvalue .gt. pmmin ) pmmin = dvalue
         call twh_read_labeled_float_array(lfailure,4,'pmpivmt',nmolty
     &        ,.true.,idname,dmarray)
         if ( lfailure ) return
c        --- set and error check pmpivmt
         pmax = 0.0d0
         do imolty = 1,nmolty
            dvalue = dmarray(imolty)
            call twh_pmpivmt(GLB_SET,imolty,dvalue)
            if ( dvalue .gt. pmax ) pmax = dvalue
         enddo
         if ( pmax .lt. 1.0d0 ) then
            write(6,*) 'READTOWHEE: at least one value of pmpivmt must'
     &           ,' be 1.0d0'
            lfailure = .true.
            return
         endif
      endif

c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     --- read concerted rotation move info
      if (twh_check_label(lfailure,4,'pmconrot')) then
         dvalue = twh_read_labeled_float(lfailure,4,'pmconrot',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_pmconrot(GLB_SET,dvalue)
         if ( dvalue .gt. pmmin ) pmmin = dvalue
         call twh_read_labeled_float_array(lfailure,4,'pmcrmt',nmolty
     &        ,.true.,idname,dmarray)
         if ( lfailure ) return
c        --- set and error check pmcrmt
         pmax = 0.0d0
         do imolty = 1,nmolty
            dvalue = dmarray(imolty)
            call twh_pmcrmt(GLB_SET,imolty,dvalue)
            if ( dvalue .gt. pmax ) pmax = dvalue
         enddo
         if ( pmax .lt. 1.0d0 ) then
            write(6,*) 'READTOWHEE: at least one value of pmcrmt must'
     &           ,' be 1.0d0'
            lfailure = .true.
            return
         endif
      endif

c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     --- read backbone concerted rotation move info
      if (twh_check_label(lfailure,4,'pmcrback')) then
         dvalue = twh_read_labeled_float(lfailure,4,'pmcrback',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_pmcrback(GLB_SET,dvalue)
         if ( dvalue .gt. pmmin ) pmmin = dvalue
         call twh_read_labeled_float_array(lfailure,4,'pmcrbmt',nmolty
     &        ,.true.,idname,dmarray)
         if ( lfailure ) return
c        --- set and error check pmcrbmt
         pmax = 0.0d0
         do imolty = 1,nmolty
            dvalue = dmarray(imolty)
            call twh_pmcrbmt(GLB_SET,imolty,dvalue)
            if ( dvalue .gt. pmax ) pmax = dvalue
         enddo
         if ( pmax .lt. 1.0d0 ) then
            write(6,*) 'READTOWHEE: at least one value of pmcrbmt must'
     &           ,' be 1.0d0'
            lfailure = .true.
            return
         endif
      endif

c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     --- read in plane shift info
      if (twh_check_label(lfailure,4,'pmplane')) then
         dvalue = twh_read_labeled_float(lfailure,4,'pmplane',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_pmplane(GLB_SET,dvalue)
         if ( dvalue .gt. pmmin ) pmmin = dvalue
         call twh_read_labeled_float_array(lfailure,4,'pmplanebox'
     &        ,numboxes, .true., idname,dboxarray)
         if ( lfailure ) return
c        --- set and error check pmplanebox
         pmax = 0.0d0
         do ibox = 1,numboxes
            dvalue = dboxarray(ibox)
            call twh_pmplanebox(GLB_SET,ibox,dvalue)
            if ( dvalue .gt. pmax ) pmax = dvalue
         enddo
         if ( pmax .lt. 1.0d0 ) then
            write(6,*) 'READTOWHEE: at least one value of pmplanebox'
     &           ,' must be 1.0d0'
            lfailure = .true.
            return
         endif
         dvalue = twh_read_labeled_float(lfailure,4,'planewidth',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_planewidth(GLB_SET,dvalue)
         if ( dvalue .lt. 0.0d0 ) then
            write(6,*) 'READTOWHEE: planewidth must be greater than',
     &            ' zero'
            lfailure = .true.
            return
         endif
      endif

c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     --- read in row shift info
      if (twh_check_label(lfailure,4,'pmrow')) then
         dvalue = twh_read_labeled_float(lfailure,4,'pmrow',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_pmrow(GLB_SET,dvalue)
         if ( dvalue .gt. pmmin ) pmmin = dvalue
         call twh_read_labeled_float_array(lfailure,4,'pmrowbox'
     &        ,numboxes,.true.,idname,dboxarray)
         if ( lfailure ) return
c        --- set and error check pmrowbox
         pmax = 0.0d0
         do ibox = 1,numboxes
            dvalue = dboxarray(ibox)
            call twh_pmrowbox(GLB_SET,ibox,dvalue)
            if ( dvalue .gt. pmax ) pmax = dvalue
         enddo
         if ( pmax .lt. 1.0d0 ) then
            write(6,*) 'READTOWHEE: at least one value of pmrowbox must'
     &           ,' be 1.0d0'
            lfailure = .true.
            return
         endif

         dvalue = twh_read_labeled_float(lfailure,4,'rowwidth',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_rowwidth(GLB_SET,dvalue)
         if ( dvalue .lt. 0.0d0 ) then
            write(6,*) 'READTOWHEE: rowwidth must be greater than zero'
            lfailure = .true.
            return
         endif
      endif

c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     --- read translation of single atom info
      if (twh_check_label(lfailure,4,'pmtraat')) then
         dvalue = twh_read_labeled_float(lfailure,4,'pmtraat',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_pmtraat(GLB_SET,dvalue)
         if ( dvalue .gt. pmmin ) pmmin = dvalue
         call twh_read_labeled_float_array(lfailure,4,'pmtamt',nmolty
     &        ,.true.,idname,dmarray)
         if ( lfailure ) return
c        --- set and error check pmtamt
         pmax = 0.0d0
         do imolty = 1,nmolty
            dvalue = dmarray(imolty)
            call twh_pmtamt(GLB_SET,imolty,dvalue)
            if ( dvalue .gt. pmax ) pmax = dvalue
         enddo
         if ( pmax .lt. 1.0d0 ) then
            write(6,*) 'READTOWHEE: at least one value of pmtamt must'
     &           ,' be 1.0d0'
            lfailure = .true.
            return
         endif

c        --- read in the initial atom translation maximum displacements
         dvalue = twh_read_labeled_float(lfailure,4,'rmtraa',.true.
     &        ,idname)
         if ( lfailure ) return
c        --- error check
         if ( dvalue .le. 0.0d0 ) then
            write(6,*) 'READTOWHEE: rmtraax must be positive'
            lfailure = .true.
            return
         endif
c        --- set all of the values for the moltyps to the same input
c        --- value.  For a restart this will be overwritten later
         call twh_rmtraa(GLB_INIT,nmolty,numboxes,dvalue)

c        --- read in the target acceptance rate for atom translate
         dvalue = twh_read_labeled_float(lfailure,4,'tatraa',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_tatraa(GLB_SET,dvalue)
c        --- error check
         if ( dvalue .le. 0.0d0 .or. dvalue .ge. 1.0d0 ) then
            write(6,*) 'READTOWHEE: tatraa must be between 0.0 and 1.0'
     &           ,' exclusive'
            lfailure = .true.
            return
         endif
      endif

c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     --- composite move info here
      if (twh_check_label(lfailure,4,'pmcomposite')) then
         dvalue = twh_read_labeled_float(lfailure,4,'pmcomposite'
     &        ,.true.,idname)
         if ( lfailure ) return
         call twh_pmcomposite(GLB_SET,dvalue)
         if ( dvalue .gt. pmmin ) pmmin = dvalue
c        --- error check pmcomposite
         if ( dvalue .lt. 0.0d0 ) then
            write(6,*) 'READTOWHEE: pmcomposite must be non-negative'
            lfailure = .true.
            return
         endif

         call twh_read_labeled_float_array(lfailure,4,'pmcomt',nmolty
     &        ,.true.,idname,dmarray)
         if ( lfailure ) return
c        --- set and error check pmcomt
         pmax = 0.0d0
         do imolty = 1,nmolty
            dvalue = dmarray(imolty)
            call twh_pmcomt(GLB_SET,imolty,dvalue)
            if ( dvalue .gt. pmax ) pmax = dvalue
         enddo
         if ( pmax .lt. 1.0d0 ) then
            write(6,*) 'READTOWHEE: at least one value of pmcomt must'
     &           ,' be 1.0d0'
            lfailure = .true.
            return
         endif

c        --- read in the initial translation maximum displacements
         dvalue = twh_read_labeled_float(lfailure,4,'rmcomtra',.true.
     &        ,idname)
         if ( lfailure ) return
c        --- error check
         if ( dvalue .lt. 0.0d0 ) then
            write(6,*) 'READTOWHEE: rmcomtra must be positive'
            lfailure = .true.
            return
         endif
c        --- set all of the values for the moltyps to the same value
         do ibox = 1,numboxes
            do imolty = 1,nmolty
               call twh_rmcomtra(GLB_SET,imolty,ibox,dvalue)
            enddo
         enddo

c        --- read in the initial rotation maximum displacements
         dvalue = twh_read_labeled_float(lfailure,4,'rmcomrot',.true.
     &        ,idname)
         if ( lfailure ) return
c        --- error check
         if ( dvalue .lt. 0.0d0 ) then
            write(6,*) 'READTOWHEE: rmcomrot must be positive'
            lfailure = .true.
            return
         endif
c        --- set all of the values for the moltyps to the same value
         do ibox = 1,numboxes
            do imolty = 1,nmolty
               call twh_rmcomrot(GLB_SET,imolty,ibox,dvalue)
            enddo
         enddo

      endif

c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     --- read translation of center-of-mass info
      if (twh_check_label(lfailure,4,'pmtracm')) then
         dvalue = twh_read_labeled_float(lfailure,4,'pmtracm',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_pmtracm(GLB_SET,dvalue)
         if ( dvalue .gt. pmmin ) pmmin = dvalue
         call twh_read_labeled_float_array(lfailure,4,'pmtcmt',nmolty
     &        ,.true.,idname,dmarray)
         if ( lfailure ) return
c        --- set and error check pmtcmt
         pmax = 0.0d0
         do imolty = 1,nmolty
            dvalue = dmarray(imolty)
            call twh_pmtcmt(GLB_SET,imolty,dvalue)
            if ( dvalue .gt. pmax ) pmax = dvalue
         enddo
         if ( pmax .lt. 1.0d0 ) then
            write(6,*) 'READTOWHEE: at least one value of pmtcmt must'
     &           ,' be 1.0d0'
            lfailure = .true.
            return
         endif

c        --- read in the initial COM translation maximum displacements
         dvalue = twh_read_labeled_float(lfailure,4,'rmtrac',.true.
     &        ,idname)
         if ( lfailure ) return
c        --- error check
         if ( dvalue .le. 0.0d0 ) then
            write(6,*) 'READTOWHEE: rmtrac must be positive'
            lfailure = .true.
            return
         endif
c        --- set all of the values for the moltyps to the same input 
c        --- value
         call twh_rmtrac(GLB_INIT,nmolty,numboxes,dvalue)

c        --- read in the target acceptance rate for COM translate
         dvalue = twh_read_labeled_float(lfailure,4,'tatrac',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_tatrac(GLB_SET,dvalue)
c        --- error check
         if ( dvalue .le. 0.0d0 .or. dvalue .ge. 1.0d0 ) then
            write(6,*) 'READTOWHEE: tatrac must be between 0.0 and 1.0'
     &           ,' exclusive'
            lfailure = .true.
            return
         endif
      endif

c     --- read rotation info
      if (twh_check_label(lfailure,4,'pmrotate')) then
         dvalue = twh_read_labeled_float(lfailure,4,'pmrotate',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_pmrotate(GLB_SET,dvalue)
         if ( dvalue .gt. pmmin ) pmmin = dvalue
c        --- error check pmrotate
         if ( dvalue .lt. 0.0d0 ) then
            write(6,*) 'READTOWHEE: pmrotate must be non-negative'
            lfailure = .true.
            return
         endif
         call twh_read_labeled_float_array(lfailure,4,'pmromt',nmolty
     &        ,.true.,idname,dmarray)
         if ( lfailure ) return
c        --- set and error check pmromt
         pmax = 0.0d0
         do imolty = 1,nmolty
            dvalue = dmarray(imolty)
            call twh_pmromt(GLB_SET,imolty,dvalue)
            if ( dvalue .gt. pmax ) pmax = dvalue
         enddo
         if ( pmax .lt. 1.0d0 ) then
            write(6,*) 'READTOWHEE: at least one value of pmromt must'
     &           ,' be 1.0d0'
            lfailure = .true.
            return
         endif
c        --- read in the initial rotation maximum displacements
         dvalue = twh_read_labeled_float(lfailure,4,'rmrot',.true.
     &        ,idname)
         if ( lfailure ) return
c        --- error check
         if ( dvalue .le. 0.0d0 ) then
            write(6,*) 'READTOWHEE: rmrot must be positive'
            lfailure = .true.
            return
         endif
c        --- set all of the values for the moltyps to the same value
         do ibox = 1,numboxes
            do imolty = 1,nmolty
               call twh_rmrot(GLB_SET,imolty,ibox,dvalue)
            enddo
         enddo

c        --- read in the target acceptance rate for rotate
         dvalue = twh_read_labeled_float(lfailure,4,'tarot',.true.
     &        ,idname)
         if ( lfailure ) return
         call twh_tarot(GLB_SET,dvalue)
c        --- error check
         if ( dvalue .le. 0.0d0 .or. dvalue .ge. 1.0d0 ) then
            write(6,*) 'READTOWHEE: tarot must be between 0.0 and 1.0'
     &           ,' exclusive'
            lfailure = .true.
            return
         endif
      endif


c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     --- move error checking
c        --- make sure at least one move has a probability >= 1.0d0
         if ( pmmin .lt. 1.0d0 ) then
            write(6,*) 'READTOWHEE: at least one move must have a'
     &           ,' probability of 1.0'
            lfailure = .true.
            return
         endif



c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     --- get all of the variables associated with the details of the 
c     --- configurational-bias monte carlo moves
      call twh_rwcbmc(lfailure,READ_FLAG,4)
      if ( lfailure) return

c     --- output the towhee_altinp header information
      call twh_writetowhee(lfailure,.false.,0)
      if ( lfailure ) return

c     --- read bead potential information
      do imolty = 1, nmolty
c        --- initialize lpdbset to false to show we have not setup 
c        --- pdbname
         lpdbset(imolty) = .false.
c        --- initialize the molecular mass
         call twh_masstotal(GLB_SET,imolty,0.0d0)

c        --- read input style
         call twh_read_labeled_string(lfailure,4,'input_style'
     &        , .true.,idname,input_style)
         if ( lfailure ) return
c        --- error check input_style with potentialstyle
         if ( ( potentialstyle .eq. POT_EXTERNAL
     &        .and. input_style .ne. 'atomic' ) .or.
     &        ( potentialstyle .ne. POT_EXTERNAL
     &        .and. input_style .eq. 'atomic' ) ) then
            write(6,*)
     &           'READTOWHEE: input_style atomic and potentialstyle'
     &           ,' external must be used together'
            lfailure = .true.
            return
         endif

c        --- initialize
         lpdbnames = .false.
         lreadnunit = .false.
         lreadnmaxcbmc = .false.
         lreadpdbnames = .false.
         lreadforcefield = .false.
         lreadchargeassign = .false.
         lreadimpassign = .false.
         lreadmatchstyle = .false.
c        --- initialize the vibration order to the defaults
         call twh_vibration_order(GLB_INIT,NUMAX,NNBOND,'any')

         if ( input_style .eq. 'explicit' ) then
            lreadnunit = .true.
            lreadnmaxcbmc = .true.
         elseif ( input_style .eq. 'polypeptide builder' ) then
            lpdbnames = .true.
            lreadnunit = .true.
            lreadnmaxcbmc = .true.
         elseif ( input_style .eq. 'basic connectivity map' ) then
            lreadnunit = .true.
            lreadnmaxcbmc = .true.
            lreadpdbnames = .true.
            lreadforcefield = .true.
            lreadchargeassign = .true.
         elseif ( input_style .eq. 'nucleic acid builder' ) then
            lpdbnames = .true.
            lreadnunit = .true.
            lreadnmaxcbmc = .true.
         elseif ( input_style .eq. 'nanotube builder' ) then
c           --- all false, but not an error
         elseif ( input_style .eq. 'atomic' ) then
c           --- all false, but not an error
         elseif ( input_style .eq. 'advanced connectivity map' ) then
            lreadnunit = .true.
            lreadnmaxcbmc = .true.
            lreadpdbnames = .true.
            lreadforcefield = .true.
            lreadchargeassign = .true.
            lreadimpassign = .true.
            lreadmatchstyle = .true.
         else
            write(6,*) 'READTOWHEE: unknown input_style'
            write(6,*) input_style
            lfailure = .true.
            return
         endif

         if ( lreadnunit ) then
c           --- read number of units
            numunits = twh_read_labeled_integer(lfailure,4,'nunit'
     &           ,.true.,idname)
            if ( lfailure ) return
            call twh_nunit(GLB_SET,imolty,numunits)
c           --- error check
            if ( numunits .gt. NUMAX ) then
               write(6,*) 'READTOWHEE: nunit exceeds NUMAX'
               write(6,*) 'set numax higher in preproc.h'
               lfailure = .true.
               return
            endif
            if ( numunits .le. 0 ) then 
               write(6,*) 'READTOWHEE: nunit must be positive'
               lfailure = .true.
               return
            endif
            call twh_cp_number(GLB_GET,cpnumber)
            if ( cpnumber .eq. CP_STILLINGER_WEBER .or.
     &           cpnumber .eq. CP_EMBEDDED_ATOM_METHOD ) then
               if ( numunits .gt. 1 ) then
                  write(6,*) 'READTOWHEE: cannot use polyatomic '
     &                 ,'molecules with cplocals'
                  write(6,*) 'That have multibody nonbonded terms'
                  lfailure = .true.
                  return
               endif
            endif
         endif

         if ( lreadnmaxcbmc ) then
c           --- read max units for cbmc
            ivalue = twh_read_labeled_integer(lfailure,4,'nmaxcbmc'
     &           ,.true., idname)
            if ( lfailure ) return
            call twh_nmaxcbmc(GLB_SET,imolty,ivalue)
c           --- error check
            if ( ivalue .le. 0 ) then
               write(6,*) 'READTOWHEE: nmaxcbmc must be positive'
               lfailure = .true.
               return
            endif
            if ( ivalue .gt. NUMAX ) then
               write(6,*) 'READTOWHEE: nmaxcbmc > NUMAX'
               write(6,*) 'set numax higher in preproc.h'
               lfailure = .true.
               return
            endif
         endif

         if ( lreadpdbnames ) then
c           --- read in the pdbnames logical
            lpdbnames = twh_read_labeled_logical(lfailure,4
     &           ,'lpdbnames',.true. ,idname)
            if ( lfailure ) return
         endif

         if ( lreadforcefield ) then
c           --- read in the forcefield
            call twh_read_labeled_string(lfailure,4,'forcefield'
     &           ,.false.,idname,forcefield)
            if ( lfailure ) return
            write(6,'(a13,a10,a12)') '   using the ',forcefield
     &           ,' force field'
         endif

         if ( lreadchargeassign ) then
c           --- read in the charge_assignment
            call twh_read_labeled_string(lfailure,4,'charge_assignment'
     &           ,.true. , idname,charge_assignment)
            if ( lfailure ) return
            call twh_charge_assignment(GLB_SET,imolty,charge_assignment)
         endif

         if ( lreadimpassign ) then
c           --- read in the improper_assignment 
            call twh_read_labeled_string(lfailure,4
     &           ,'improper_assignment',.true.,idname,imp_assign)
            if ( lfailure ) return
         else
c           --- default to manual
            imp_assign = 'manual'
         endif

         if ( lreadmatchstyle  ) then
c           --- read in the match_style
            call twh_read_labeled_string(lfailure,4,'match_style'
     &           ,.true.,idname,match_style)
            if ( lfailure ) return
         else
c           --- default to standard
            match_style = 'standard'
         endif

c        --- cycle through each atom/bead 
         if ( input_style .eq. 'explicit' ) then
c           --- complete declaration of all atom types, vibrations
c           --- , angles, torsions, impropers, angle-angle in terms of 
c           --- the parameter numbers used in ffbond and ffnonbond
c           --- see if we need to read in pdb names
            call twh_read_inpstyle_explicit(lfailure,imolty,lpdbnames)
            if ( lfailure ) return
            if ( lpdbnames ) lpdbset(imolty) = .true.
         elseif ( input_style .eq. 'polypeptide builder' ) then
c           --- input is as peptides in a protein using a two letter
c           --- code where the first letter is the peptide and the 
c           --- second letter is the protonation state (if applicable)
c           --- this is handled by the subroutine buildprot which then
c           --- turns this information into the input parameters needed
c           --- by the code
            call twh_buildprot(lfailure,imolty)
            if ( lfailure ) return
c           --- this will have setup pdbname
            lpdbset(imolty) = .true.
         elseif ( input_style .eq. 'basic connectivity map' .or.
     &           input_style .eq. 'advanced connectivity map' ) then
c           --- input is a collection of atoms with names using a 
c           --- four letter code, and a set of vibrations.
            call twh_buildmolec(lfailure,imolty,lpdbnames,input_style
     &           ,forcefield,imp_assign,match_style)
            if ( lfailure ) return
         elseif ( input_style .eq. 'nucleic acid builder' ) then
c           --- input is as nucleotides in a nucleic acide strand 
c           --- using a single or two letter
c           --- code where the first letter tells if it is a 
c           --- deoxyribonucleotide (d)
c           --- or an oxyribonucleotide (no first letter) and the 
c           --- second letter tells what type of nucleotide it is
            call twh_buildna(lfailure,imolty)
            if ( lfailure ) return
c           --- this will have setup pdbname
            lpdbset(imolty) = .true.
         elseif ( input_style .eq. 'nanotube builder' ) then
c           --- input is the information needed to build a nanotube
            call twh_buildnanotube(lfailure,imolty)
            if ( lfailure ) return
         elseif ( input_style .eq. 'atomic' ) then
c           --- elemental input - only valid for pure quantum approach
c           --- monatomic molecule
            call twh_nunit(GLB_SET,imolty,1)
            call twh_nmaxcbmc(GLB_SET,imolty,1)
            call twh_ntype(GLB_SET,imolty,1,imolty)
c           --- get the element name
            call twh_read_labeled_string(lfailure,4,'element',.true.
     &           ,idname,ename)
            if ( lfailure ) return
            call twh_elementname(GLB_SET,imolty,ename)
c           --- assign mass
            tenchar = 'getmass'
            call twh_getelement(lfailure,tenchar,localmass,ename,
     &           ielement)
            if ( lfailure ) return
            call twh_mass(GLB_SET,imolty,localmass)
c           --- ground state quantum energy per atom
            dvalue = twh_read_labeled_float(lfailure,4
     &           ,'quantum_gsenergy' ,.true.,idname)
            if ( lfailure ) return
            call twh_quantum_gsenergy(GLB_SET,imolty,dvalue)
         else
            write(6,*) 'READTOWHEE: non-valid input style for moltyp'
     &           ,imolty
            write(6,*) input_style
            lfailure = .true.
            return
         endif
c        --- output to the towhee_altinp
         call twh_writetowhee(lfailure,lpdbnames,imolty)
         if ( lfailure ) return
      enddo

c     --- finished reading input from towhee_input
c     --- determine the largest value of nunit
      maxunit = 0
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,numunits)
         if ( numunits .gt. maxunit ) maxunit = numunits
      enddo
c     --- allocate arrays based maxunit
      call twh_allocate_maxunit(maxunit)
c     --- compute the total number of atoms
      natoms = 0
      do imolty = 1,nmolty
         call twh_nmolectyp(GLB_GET,imolty,ivalue)
         call twh_nunit(GLB_GET,imolty,numunits)
         natoms = natoms + ivalue*numunits
      enddo
c     --- allocate arrays for coords
      call twh_allocate_coords(natoms,maxunit)

c     --- determine maximum number of torsions on any atom of any
c     --- molecule in the system
      maxnumtor = 0
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,numunits)
         do iunit = 1,numunits
            call twh_intor(GLB_GET,imolty,iunit,numtor)
            if ( numtor .gt. maxnumtor ) maxnumtor = numtor
         enddo
      enddo
c     --- allocate arrays based on nmolty,maxunit and maxnumtor
      call twh_allocate_maxmolty_maxunit_maxnumtor(nmolty
     &     ,maxunit,maxnumtor)

      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,numunits)
         if ( .not. lpdbset(imolty) ) then
c           --- setup pdbname 
            do iunit = 1,numunits
               pdbname = twh_get_elementname(
     &              twh_get_ntype(imolty,iunit))
               call twh_pdbname(GLB_SET,imolty,iunit,pdbname)
            enddo
         endif
      enddo

c     --- compute the molecular weight of the molecules
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,numunits)
         masstotal = 0.0d0
         do iunit = 1,numunits
            call twh_ntype(GLB_GET,imolty,iunit,iutemp)
            masstotal = masstotal + twh_get_mass(iutemp)
         enddo
         call twh_masstotal(GLB_SET,imolty,masstotal)
      enddo

c     --- read in the initial SAFE parameters if needed
      call twh_rwsafe(lfailure,READ_FLAG)
      if ( lfailure ) return
c     --- check the molecule input to make sure it is consistent
      call twh_checkstruc(lfailure)
      if ( lfailure ) return
      call twh_atd_cbmc_setup(lfailure)
      if ( lfailure ) return

c     --- work through the molecular architecture and see if there 
c     --- are rings
c     --- and if they are simple enough to be sampled in this code
      do imolty = 1,nmolty
         call twh_findrings(imolty)
c        --- error checking for the backbone concerted rotation move 
c        --- meant for proteins
         call twh_pmcrback(GLB_GET,dvalue)
         if ( dvalue .gt. 0.0d0) then
            call twh_pmcrbmt(GLB_GET,imolty,pmtest)
            if ( pmtest .gt. 0.0d0) then
               call twh_backlist(GLB_GET,imolty,0,backlist)
               if ( backlist .eq. 0) then
                   write(6,*) 'READTOWHEE: Using backbone concerted' 
                   write(6,*) 'rotation move (pmcrback) for a system'
                   write(6,*) 'which does not consist of proteins.' 
                   write(6,*) 'Please use regular concerted rotation'
                   write(6,*) 'move (pmconrot).'     
                   lfailure = .true.
                   return
               endif
            endif
         endif
      enddo

      if ( potentialstyle .eq. POT_INTERNAL ) then
         call twh_setclassical(lfailure,numtypes)
         if ( lfailure ) return
      endif
      if ( potentialstyle .eq. POT_EXTERNAL ) then
c        --- not much external stuff to set up so still doing it here
c        --- if this gets larger then will split into a separate routine
c        --- set the minimum box dimension to 1.0 just for use when
c        --- updating the maximum translational displacements
         call twh_minboxdim(GLB_SET,1.0d0)
c        --- set the pressurestyle to thermodynamic
         call twh_pressurestyle(GLB_SET,PRES_THERMO)
         call twh_nener(GLB_GET,nener)
         nener = nener + 1
         call twh_external_code(GLB_GET,external_code)
         if ( external_code .eq. 'lcao' ) then
            keyname = 'Total Quantum'
            call twh_vname(GLB_SET,nener,keyname)
            call twh_nener(GLB_SET,nener)
         elseif ( external_code .eq. 'KIM' ) then
            keyname = 'Total KIM'
            call twh_vname(GLB_SET,nener,keyname)
            call twh_nener(GLB_SET,nener)
         else
            write(6,*) 'READTOWHEE: unknown external_code'
            write(6,*) external_code
            lfailure = .true.
            return
         endif

      endif

c     --- output potentialstyle
      if ( potentialstyle .eq. POT_INTERNAL ) then
         write(6,'(a)') 'Energies exclusively from internal potentials'
      elseif ( potentialstyle .eq. POT_EXTERNAL ) then
         write(6,'(a)') 'Energies exclusively from external potentials'
      endif

c     --- output information about nonbonded
      call twh_writenonbond(6,lfailure)
      if ( lfailure ) return

c     --- fetch the minboxdim
      call twh_minboxdim(GLB_GET,minboxdim)

c     --- re-check box dimensions
      do ibox = 1,numboxes
         call twh_checkhmatrix(ibox,.true.,lfailure)
         if ( lfailure ) then
            write(6,*) 'READTOWHEE: bad box dimensions in box',ibox
            write(6,*) 'all of the dimensions of hmatrix must be '
            write(6,*) 'larger than minboxdim'
            write(6,*) 'minboxdim ',minboxdim
            do idim = 1,3
               do jdim = 1,3
                  call twh_hmatrix(GLB_GET,ibox,idim,jdim,fatmp(jdim))
               enddo
               write(6,*) 'idim hmatrix ',idim,(fatmp(jdim),jdim=1,3)
            enddo
            return
         endif
      enddo

      if ( potentialstyle .eq. POT_INTERNAL ) then
c        --- error check avb1rad
         call twh_avb1rad(GLB_GET,dvalue)
         if ( dvalue .gt. 0.5*minboxdim ) then
            write(6,*) 'minboxdim ',minboxdim
            write(6,*) 'avb1rad must be less than minboxdim/2'
            lfailure = .true.
c           --- failure return is down below after more error messages
         endif
c        --- error check avb2rad
         call twh_avb2rad(GLB_GET,dvalue)
         if ( dvalue .gt. 0.5*minboxdim ) then
            write(6,*) 'minboxdim ',minboxdim
            write(6,*) 'avb2rad must be less than minboxdim/2'
            lfailure = .true.
c           --- failure return is down below after more error messages
         endif
c        --- error check avb3rad
         call twh_avb3rad(GLB_GET,dvalue)
         if ( dvalue .gt. 0.5*minboxdim ) then
            write(6,*) 'minboxdim ',minboxdim
            write(6,*) 'avb3rad must be less than minboxdim/2'
            lfailure = .true.
c           --- failure return is down below after more error messages
         endif
      endif

      if ( lfailure ) then
         call twh_cp_number(GLB_GET,cpnumber)
         if ( cpnumber .eq. CP_HARD_SPHERE .or.
     &        cpnumber .eq. CP_REPULSIVE_SPHERE ) then
            write(6,*) 'minboxdim is equal to twice the largest'
            write(6,*) 'hard sphere diameter'
         elseif ( cpnumber .eq. CP_STILLINGER_WEBER .or.
     &           cpnumber .eq. CP_SW_PAIR_ONLY ) then
            write(6,*) 'minboxdim is four times the largest a parameter'
            write(6,*) 'a is nbcoeff(6)'
         elseif ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD .or.
     &           cpnumber .eq. CP_EAM_PAIR_ONLY ) then
            write(6,*) 'minboxdim is twice the largest distance in '
            write(6,*) 'the lookup tables'
         elseif ( cpnumber .eq. CP_SQUARE_WELL .or.
     &           cpnumber .eq. CP_REPULSIVE_WELL ) then
            write(6,*) 'minboxdim is equal to twice the largest'
            write(6,*) 'attractive well diameter'
         else
            write(6,*) 'minboxdim/2 is equal to rcut'
         endif
         return
      endif

c     --- determine if we should use the additional COM cutoff
c     --- this is not useful for monatomic systems
      lusecom = .false.
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,numunits)
         if ( numunits .gt. 1 ) lusecom = .true.
      enddo
c     --- set lusecom
      call twh_lusecom(GLB_SET,lusecom)

c     --- set up the inclusion table
      call twh_setinclude
 
c     --- Begin outputting information about the simulation
      if ( ensemble .eq. ENS_NPT ) then
         if ( numboxes .eq. 1 ) then
            write(6,*) 'Isobaric-isothermal ensemble'
         else
            write(6,*) 'Isobaric-isothermal Gibbs ensemble'
         endif
      elseif ( ensemble .eq. ENS_NVT ) then
         if ( numboxes .eq. 1 ) then
            write(6,*) 'Canonical ensemble'
         else
            write(6,*) 'Canonical Gibbs ensemble'
         endif
      elseif ( ensemble .eq. ENS_muVT ) then
         write(6,*) 'Grand Canonical ensemble'
      else
         write(6,*) 'READTOWHEE: unknown ensemble'
         lfailure = .true.
         return
      endif

      write(6,*) '3-dimensional periodic box'

      write(6,*) 'Additional Center-of-Mass cutoff'
      call twh_ldualcutoff(GLB_GET,ldualcutoff)
      if ( ldualcutoff ) then
         write(6,*) 'Dual Cutoff Configurational-bias Monte Carlo'
      endif

      write(6,*) 'Coupled-decoupled Configurational-bias MC'

      call twh_nfield(GLB_GET,nfield)
      if ( nfield .ne. 0 ) write(6,*) 'External fields specified'

      call twh_lshift(GLB_GET,lshift)
      call twh_ltailc(GLB_GET,ltailc)
      if ( lshift .and. ltailc ) then
         write(6,*) 'READTOWHEE: cannot have shift and tail'
         lfailure = .true.
         return
      endif
         
      call twh_lcoulomb(GLB_GET,lcoulomb)
      if ( lcoulomb ) then
         write(6,*) 
     &        'Coulombic inter- and intra-molecular interactions'
         call twh_lewald(GLB_GET,lewald)
         if ( lewald ) write(6,*) '    with an Ewald sum '
         write(6,*) '    including the real-space terms up to half '
     &        ,'the shortest box length'
      endif
      write(6,*)
      do imolty = 1,nmolty
         call twh_masstotal(GLB_GET,imolty,masstotal)
         write(6,40) 'Molecular mass for molecule type',imolty
     &        ,'is',masstotal,'g/mol'
 40      format(a32,1x,i5,1x,a2,1x,f10.4,1x,a5)
      enddo

c     --- check to make sure we are not attempting moves on any
c     --- molecule that is not a bonded graph
 67   format('Molecule type',1x,i4,1x
     &     ,'is not a true bonded graph molecule and'
     &     ,1x
     &     ,'has a nonzero move probability of'
     &     ,1x,f12.10,1x,'for the',1x,a) 
      do imolty = 1,nmolty
         call twh_lbondgraph(GLB_GET,imolty,lbondgraph)
         if ( .not. lbondgraph ) then
c           --- check every individual move type to make sure none
c           --- are performed on this molecule type

c           --- 2 box rotational bias
            call twh_pm2boxrbswap(GLB_GET,pmtest)
            call twh_pm2rbswmt(GLB_GET,imolty,pmtype)
            dvalue = pmtest*pmtype
            if ( (dvalue) .gt. SMALLEST ) then
               write(6,67) imolty,dvalue,'2 box rbswap move'
               lfailure = .true.
               return
            endif
c           --- 2 box configurational bias
            call twh_pm2boxcbswap(GLB_GET,pmtest)
            call twh_pm2cbswmt(GLB_GET,imolty,pmtype)
            dvalue = pmtest*pmtype
            if ( (dvalue) .gt. SMALLEST ) then
               write(6,67) imolty,dvalue,'2 box cbswap move'
               lfailure = .true.
               return
            endif
c           --- uvt swap
            call twh_pmuvtcbswap(GLB_GET,pmtest)
            call twh_pmuvtcbmt(GLB_GET,imolty,pmtype)
            dvalue = pmtest*pmtype
            if ( (dvalue) .gt. SMALLEST ) then
               write(6,67) imolty,dvalue,'uvtswap move'
               lfailure = .true.
               return
            endif
c           --- 1box cb swap
            call twh_pm1boxcbswap(GLB_GET,pmtest)
            call twh_pm1cbswmt(GLB_GET,imolty,pmtype)
            dvalue = pmtest*pmtype
            if ( (dvalue) .gt. SMALLEST ) then
               write(6,67) imolty,dvalue,'1 box cbswap move'
               lfailure = .true.
               return
            endif
c           --- cb move
            call twh_pmcb(GLB_GET,pmtest)
            call twh_pmcbmt(GLB_GET,imolty,pmtype)
            dvalue = pmtest*pmtype
            if ( (dvalue) .gt. SMALLEST ) then
               write(6,67) imolty,dvalue,'configurational-bias move'
               lfailure = .true.
               return
            endif
c           --- backbone cb move
            call twh_pmback(GLB_GET,pmtest)
            call twh_pmbkmt(GLB_GET,imolty,pmtype)
            dvalue = pmtest*pmtype
            if ( (dvalue) .gt. SMALLEST ) then
               write(6,67) imolty,dvalue,'backbone cb move'
               lfailure = .true.
               return
            endif
c           --- sidechain cb move
            call twh_pmcbside(GLB_GET,pmtest)
            call twh_pmcbsidemt(GLB_GET,imolty,pmtype)
            dvalue = pmtest*pmtype
            if ( (dvalue) .gt. SMALLEST ) then
               write(6,67) imolty,dvalue,'sidechain cb move'
               lfailure = .true.
               return
            endif
c           --- pivot move
            call twh_pmpivot(GLB_GET,pmtest)
            call twh_pmpivmt(GLB_GET,imolty,pmtype)
            dvalue = pmtest*pmtype
            if ( (dvalue) .gt. SMALLEST ) then
               write(6,67) imolty,dvalue,'pivot move'
               lfailure = .true.
               return
            endif
c           --- concerted rotation move
            call twh_pmconrot(GLB_GET,pmtest)
            call twh_pmcrmt(GLB_GET,imolty,pmtype)
            dvalue = pmtest*pmtype
            if ( (dvalue) .gt. SMALLEST ) then
               write(6,67) imolty,dvalue,'concerted rotation move'
               lfailure = .true.
               return
            endif
c           --- cr backbone move
            call twh_pmcrback(GLB_GET,pmtest)
            call twh_pmcrbmt(GLB_GET,imolty,pmtype)
            dvalue = pmtest*pmtype
            if ( (dvalue) .gt. SMALLEST ) then
               write(6,67) imolty,dvalue,'cr backbone move'
               lfailure = .true.
               return
            endif
c           --- single atom translation move
            call twh_pmtraat(GLB_GET,pmtest)
            call twh_pmtamt(GLB_GET,imolty,pmtype)
            dvalue = pmtest*pmtype
            if ( (dvalue) .gt. SMALLEST ) then
               write(6,67) imolty,dvalue,'single-atom translation move'
               lfailure = .true.
               return
            endif
c           --- composite move
            call twh_pmcomposite(GLB_GET,pmtest)
            call twh_pmcomt(GLB_GET,imolty,pmtype)
            dvalue = pmtest*pmtype
            if ( (dvalue) .gt. SMALLEST ) then
               write(6,67) imolty,dvalue,'composite move'
               lfailure = .true.
               return
            endif
c           --- com translation move
            call twh_pmtracm(GLB_GET,pmtest)
            call twh_pmtcmt(GLB_GET,imolty,pmtype)
            dvalue = pmtest*pmtype
            if ( (dvalue) .gt. SMALLEST ) then
               write(6,67) imolty,dvalue,'com translation move'
               lfailure = .true.
               return
            endif
c           --- rotation move
            call twh_pmrotate(GLB_GET,pmtest)
            call twh_pmromt(GLB_GET,imolty,pmtype)
            dvalue = pmtest*pmtype
            if ( (dvalue) .gt. SMALLEST ) then
               write(6,67) imolty,dvalue,'rotation move'
               lfailure = .true.
               return
            endif

         endif
      enddo

c     --- get or generate the initial configuration of the system
      if ( linit ) then
c        --- generate a new initial configuration
         write(6,*) 'Calling initconf'
         call twh_initconf(lfailure)
         if ( lfailure ) then
            write(6,*) 'READTOWHEE: there was an error in initconf'
            return
         endif
      else
c        --- read the initial conformation from a file
         call twh_rwconf(lfailure,nchain,READ_FLAG,0,atomcount)
         if ( lfailure ) return
         call twh_lcoulomb(GLB_GET,lcoulomb)
         call twh_lewald(GLB_GET,lewald)
         if ( lcoulomb ) then
c           --- set ewald parameters: 
            do ibox = 1, numboxes
               call twh_resetcoulomb(lfailure,ibox)
               if ( lfailure ) return
               if ( lewald ) then
                  call twh_ewald_calp(GLB_GET,ibox,calp)
                  write(6,41) 'Box:',ibox,'Initial calp:',calp
 41               format(a4,1x,i5,1x,a13,1x,f10.5)
                  call twh_ewald_kmax(GLB_GET,ibox,kmax)
                  write(6,'(a4,1x,i5,1x,a13,1x,i5)')
     &                 'Box:',ibox,'Initial kmax:',kmax
               endif
            enddo
         endif
      endif

      if ( lhenry ) then
         ibox = numboxes+1
c        --- disable the ewald sum in MAXBOX
         call twh_ewald_kmax(GLB_SET,ibox,0)
         call twh_ewald_calp(GLB_SET,ibox,0.0d0)
c        --- set extremely large hmatrix for MAXBOX
         do idim = 1,3
            call twh_hmatrix(GLB_SET,ibox,idim,idim,1.0d5)
         enddo
c        --- error check box dimensions
         call twh_checkhmatrix(ibox,.true.,lfailure)
         if ( lfailure ) return
c        --- invert this matrix
         call twh_inverthmatrix(ibox)
      endif

c     --- check to see if we need the fieldcoords
      lfieldcoords = .false.
      call twh_nfield(GLB_GET,nfield)
      do ifield = 1,nfield
         call twh_fieldtype(GLB_GET,ifield,fieldtype)
         if ( fieldtype .eq. FLD_HARMONIC_ATTRACTOR ) then
c           --- check to see if this is an Initial field
            call twh_fieldchar(GLB_GET,ifield,1,fieldchar)
            if ( fieldchar .eq. 'Initial' ) lfieldcoords = .true.
         endif
      enddo
      if ( lfieldcoords ) then
c        --- allocate space to store the initial coords of all atoms
         total = 3*natoms
         dvalue = 0.0d0
         call twh_coordfield(GLB_ALLOC,total,dvalue)
c        --- store the initial coordinates of all of atoms in the system
         call twh_nchain(GLB_GET,nchain)
         do ichain = 1,nchain
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,numunits)
            do iunit = 1,numunits
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
               call twh_fieldcoordinates(GLB_SET,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
            enddo
         enddo
      endif

c     --- book keeping arrays
      do ibox = 0,numboxes
         do imolty = 1,nmolty
            idummy(imolty) = 0
         enddo      
         do ichain = 1, nchain
            if ( twh_get_nboxi(ichain) .eq. ibox ) then
               call twh_moltyp(GLB_GET,ichain,imolty)
               if ( imolty .lt. 1 .or. imolty .gt. NTMAX ) then
                  write(6,*) 'READTOWHEE: error with assigned moltyp'
                  write(6,*) 'ichain:',ichain,' moltyp:',imolty
                  lfailure = .true.
                  return
               endif
               idummy(imolty) = idummy(imolty)+1
               call twh_chainlist(GLB_SET,idummy(imolty),ibox,imolty
     &              ,ichain)
            endif
         enddo
         nmcount = 0
         do imolty = 1,nmolty
            nmcount = nmcount + idummy(imolty)
         enddo
         call twh_nchbox(GLB_GET,ibox,nchbox)
         if ( nmcount .ne. nchbox ) then
            write(6,*) 'READTOWHEE: nmcount(',nmcount,') ne nchbox ('
     &           ,nchbox,') in box:',ibox
            lfailure = .true.
            return
         endif
      enddo

c     --- set idummy counter to 0
      do imolty = 1,nmolty
         idummy(imolty)=0
      enddo

c     --- set up parall
c     --- parall is a map from a given molecule type (imolty) and nth 
c     --- instance of that type to ichain, which identifies a given
c     --- molecule/chain.  Used to get ichain when imolty and 
c     --- N (range: 1..ncmt(ibox,imolty) ) are known, often in the 
c     --- beginning of an MC move.
      do ichain = 1,nchain
         call twh_moltyp(GLB_GET,ichain,imolty)
         idummy(imolty) = idummy(imolty) + 1
         call twh_parall(GLB_SET,imolty,idummy(imolty),ichain)
      enddo

      do ibox = 1,numboxes
         call twh_ctrmas(lfailure,0,ibox,0,CTR_INITIAL)
         if ( lfailure ) return
      enddo

c     --- if energy biasing requested call twh_mapconf
      lfound = .false.
      do ibox = 1,numboxes
         do imolty = 1,nmolty
            call twh_cbmc_nb_one_generation(GLB_GET,ibox,imolty,ctemp)
            if ( ctemp .eq. 'energy bias' ) lfound = .true.
         enddo
      enddo

      if ( lfound ) then
         call twh_mapconf(lfailure)
         if ( lfailure ) then
            write(6,*) 'READTOWHEE: there was an error in mapconf'
            return
         endif
      endif

c     --- initialize movie file 
      if (moviefreq .gt. 0) then
         call twh_writemovie(lfailure,0)
         if ( lfailure ) return
      endif

c     OLD LOCATION OF NONBOND OUTPUT

c     --- write input data to stdout
      write(6,*)
      if ( stepstyle .eq. 'cycles' ) then 
         write(6,'(a20,1x,i20)') 'Number of MC cycles:',nstep
      elseif ( stepstyle .eq. 'moves' ) then
         write(6,'(a19,1x,i20)') 'Number of MC moves:',nstep
      endif
      write(6,'(a20,1x,i20)') 'Number of molecules:',nchain
      call twh_temperature(GLB_GET,temperature)
      write(6,'(a16,1x,f15.5)') 'Temperature [K]:', temperature
      if ( ensemble .eq. ENS_NPT ) then
         write(6,'(a24,1x,f15.5)') 'External pressure [kPa]:',pressure
      endif
      write(6,*)

      call twh_pmconrot(GLB_GET,dvalue)
      if ( dvalue .gt. 0.0d0) then
c        --- set up array with possible conrot starting points
         crtype = 1               
         do imolty = 1,nmolty
            count = 0
            call twh_pmcrmt(GLB_GET,imolty,pmtest)
            if ( pmtest .gt. 0.0d0) then               
               call twh_nunit(GLB_GET,imolty,numunits)
               do iunit = 1,numunits                  
                  call twh_crstart(crtype,imolty,iunit,dcrs,crback
     &                 ,isuccess)
                  if (isuccess .eq. 1) then
                     count = count + 1
                     call twh_crvalidlist(GLB_SET,imolty,count,iunit)
                  endif
               enddo
            endif
            call twh_crvalidcount(GLB_SET,imolty,count)
         enddo
      endif

     
      call twh_pmcrback(GLB_GET,dvalue)
      if ( dvalue .gt. 0.0d0) then 
         crtype = 0        
c        --- set up array with possible conrot starting points
         do imolty = 1,nmolty
            count = 0
            call twh_pmcrbmt(GLB_GET,imolty,pmtest)
            if ( pmtest .gt. 0.0d0) then
               call twh_nunit(GLB_GET,imolty,numunits)
               do iunit = 1,numunits                  
                  call twh_crstart(crtype,imolty,iunit,dcrs,crback
     &                 ,isuccess)
                  if (isuccess .eq. 1) then
                     count = count + 1
                     call twh_crvalidlist(GLB_SET,imolty,count,iunit)
                  endif
               enddo
            endif
            call twh_crvalidcount(GLB_SET,imolty,count)
         enddo
      endif

c     --- set up the growlogic arrays
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,numunits)
         count = 0
         call twh_nmaxcbmc(GLB_GET,imolty,itest)
         if ( itest .lt. numunits ) then
            do iunit = 1,numunits
               lvalid = .false.
               do ivib = 1,twh_get_invib(imolty,iunit)
                  call twh_schedule(lfailure,numunits,imolty,idum,iunit
     &                 ,itor,ivib,0)
                  if ( lfailure ) return
                  call twh_growcount(GLB_SET,imolty,iunit,ivib,itor)
                  if ( itor .le. itest ) lvalid = .true.
               enddo
               if ( lvalid ) then
                  count = count + 1
                  call twh_growvalidlist(GLB_SET,imolty,count,iunit)
               endif
            enddo
         else
c           --- all units are equally valid when not using nmaxcbmc
            do iunit = 1,numunits
               count = count + 1
               call twh_growvalidlist(GLB_SET,imolty,count,iunit)
            enddo
         endif
c        --- set the growvalidcount
         call twh_growvalidcount(GLB_SET,imolty,count)
c        --- safety check to make sure this is positive
         if ( twh_get_growvalidcount(imolty) .eq. 0 ) then
            call twh_nmaxcbmc(GLB_GET,imolty,ivalue)
            write(6,*) 'READTOWHEE: nmaxcbmc of ',ivalue
     &           ,' for molecule type ',imolty
            write(6,*) 'results in no possible growth patterns'
            write(6,*) 'please select a higher value'
            lfailure = .true.
            return
         endif
      enddo

c     --- set the thermal debroglie wavelength and uvtfactor
      call twh_beta(GLB_GET,beta)
      do imolty = 1,nmolty
         call twh_masstotal(GLB_GET,imolty,masstotal)
         dvalue = 17.458d0/( dsqrt(masstotal*temperature ))
         call twh_debroglie(GLB_SET,imolty,dvalue)
         if ( ensemble .eq. ENS_muVT ) then
            call twh_chempot(GLB_GET,imolty,chempot)
            uvtfactor = twh_expon(chempot*beta)/
     &           (dvalue*dvalue*dvalue)
            call twh_uvtfactor(GLB_SET,imolty,uvtfactor)
         endif
      enddo 
      call twh_louthist(GLB_GET,louthist)
      if (ensemble .eq. ENS_muVT .and. louthist) then
         call twh_file_his(GLB_GET,filename)
         call twh_openfile(lfailure,51,filename,'DEFAULT','UNKNOWN'
     &        ,'DEFAULT')
         if ( lfailure ) return
         call twh_temperature(GLB_GET,temperature)
         do jmolty = 1,nmolty
            call twh_masstotal(GLB_GET,jmolty,dmarray(jmolty))
            call twh_chempot(GLB_GET,jmolty,dmaltarray(jmolty))
         enddo
         write(51,'(g12.4,2x,i5,2x,4g15.6)') temperature,
     &        nmolty,(dmaltarray(jmolty),jmolty=1,nmolty),
     &        (dmarray(jmolty),jmolty=1,nmolty)   
         close(51)
      endif
      return
      end



      subroutine twh_read_inpstyle_explicit(lfailure,imolty,lpdbnames)
c     ******************************************************************
c     * This subroutine reads in data for inpstyle=0 case (explicit    *
c     * declaration of all terms).                                     *
c     *                                                                *
c     * Variables:                                                     *
c     * imolty -- the molecule type we are dealing with                *
c     * lpdbnames                                                      *
c     * lpdbset                                                        *
c     * lfailure                                                       *
c     * lpdb, lpdbset, lerror are modified in this subroutine          *
c     *                                                                *
c     * last modified 05-29-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_INAA
#define FUNCTION_GET_INBEN
#define FUNCTION_GET_INIMPROP
#define FUNCTION_GET_INTOR
#define FUNCTION_GET_NTYPE
#define FUNCTION_READ_LABELED_INTEGER
#define FUNCTION_READ_LABELED_LOGICAL
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer imolty
      logical lpdbnames,lfailure

c     --- local variables
c     --- character string scalars
      character*3 aminoshort
      character*4 pdbname
      character*30 idname
c     --- integer scalars
      integer iunit,ibead,ivib,iben,ivalue,itor,iaa,numvib
      integer numunits,aminonum
c     --- integer arrays
      integer iatmp
      dimension iatmp(4)
c     --- double precision scalars
      double precision charge

      idname = 'twh_read_inpstyle_explicit'
      lpdbnames = twh_read_labeled_logical(lfailure,4,'lpdbnames'
     &     ,.false.,idname)
      if ( lfailure ) return

c     --- if lpdbnames then don't reset pdbnames 
      call twh_nunit(GLB_GET,imolty,numunits)
      do iunit = 1, numunits
c        --- read in the atom type
         call twh_read_atom_type(lfailure,4,'unit ntype qqatom',ibead
     &        ,ivalue,charge)
         if ( lfailure ) return
         call twh_ntype(GLB_SET,imolty,iunit,ivalue)
         call twh_qqatom(GLB_SET,imolty,iunit,charge)
         write(6,36) 'bead:',ibead
     &        ,'beadtype:',twh_get_ntype(imolty,iunit)
     &        ,'charge:',charge
36            format(a5,1x,i10,1x,a9,1x,i10,1x,a7,1x,f10.5)

         if ( lpdbnames ) then
c           --- also read in the pdb name
            call twh_read_pdbname(lfailure,4
     &           ,'pdbname aminonum aminoshort',pdbname,aminonum
     &           ,aminoshort)
            if ( lfailure ) return
            call twh_pdbname(GLB_SET,imolty,iunit,pdbname)
            call twh_aminonum(GLB_SET,imolty,iunit,aminonum)
            call twh_aminoshort(GLB_SET,imolty,iunit,aminoshort)
            write(6,37) 'pdbname:',pdbname
     &           ,'aminonum:',aminonum,'aminoshort:',aminoshort
37               format(a8,1x,a4,1x,a9,1x,i10,1x,a11,1x,a3)
         endif

c        --- read in bond vibrations
         numvib = twh_read_labeled_integer(lfailure,4,
     &      'vibration', .false., idname)
         if ( lfailure ) return
         if ( numvib .gt. NNBOND ) then
            write(6,*) 'imolty:',imolty,'iunit:',iunit,'invib:',numvib
            write(6,*) 'READTOWHEE: too many vibrations'
            lfailure = .true.
            return
         endif
         call twh_invib(GLB_SET,imolty,iunit,numvib)
         do ivib = 1,numvib
            call twh_read_integer_array(lfailure,4,'vibration',
     &         2,.false.,iatmp)
            if ( lfailure ) return
            call twh_ijvib(GLB_SET,imolty,iunit,ivib,iatmp(1))
            call twh_itvib(GLB_SET,imolty,iunit,ivib,iatmp(2))
         enddo

c        --- read in bond bending
         iatmp(1) = twh_read_labeled_integer(lfailure,4,'bending'
     &        , .false., idname)
         if ( lfailure ) return
         call twh_inben(GLB_SET,imolty,iunit,iatmp(1))
         if ( twh_get_inben(imolty,iunit) .gt. MAXBEND ) then
            write(6,*) 'READTOWHEE: too many bends'
            lfailure = .true.
            return
         endif
         do iben = 1, twh_get_inben(imolty,iunit)
            call twh_read_integer_array(lfailure,4,'bending',
     &         4,.false.,iatmp)
            if ( lfailure ) return
            call twh_ijben2(GLB_SET,imolty,iunit,iben,iatmp(1))
            call twh_ijben3(GLB_SET,imolty,iunit,iben,iatmp(2))
            if ( iatmp(3) .lt. 1 ) then
               write(6,*) 'READTOWHEE: invalid bending type'
               write(6,*) 'imolty,iunit,iben,ijben2,ijben3,itben'
     &              ,imolty,iunit,iben,iatmp(1),iatmp(2),iatmp(3)
               lfailure = .true.
               return
            endif
            call twh_itben(GLB_SET,imolty,iunit,iben,iatmp(3))
            call twh_orderben(GLB_SET,imolty,iunit,iben,iatmp(4))
         enddo

c        --- read in torsions
         iatmp(1) = twh_read_labeled_integer(lfailure,4,'torsion'
     &        , .false., idname)
         if ( lfailure ) return
         call twh_intor(GLB_SET,imolty,iunit,iatmp(1))
         if ( twh_get_intor(imolty,iunit) .gt. MAXTOR ) then
            write(6,*) 'READTOWHEE: too many torsions'
            lfailure = .true.
            return
         endif
         do itor = 1, twh_get_intor(imolty,iunit)
            call twh_read_integer_array(lfailure,4,'torsion',
     &         4,.false.,iatmp)
            if ( lfailure ) return
            call twh_ijtor2(GLB_SET,imolty,iunit,itor,iatmp(1))
            call twh_ijtor3(GLB_SET,imolty,iunit,itor,iatmp(2))
            call twh_ijtor4(GLB_SET,imolty,iunit,itor,iatmp(3))
            call twh_ittor(GLB_SET,imolty,iunit,itor,iatmp(4))
         enddo

c        --- angle - angle (compass)
c        --- specified from atom J in the following format K I L
c        --- where the angle - angle is K-J-I and K-J-L
         iatmp(1) = twh_read_labeled_integer(lfailure,4,'angle-angle'
     &        , .false., idname)
         if ( lfailure ) return
         call twh_inaa(GLB_SET,imolty,iunit,iatmp(1))
         if ( twh_get_inaa(imolty,iunit) .gt. MAXAA ) then
            write(6,*) 'READTOWHEE: too many angle-angle terms'
            lfailure = .true.
            return
         endif
         do iaa = 1,twh_get_inaa(imolty,iunit)
            call twh_read_integer_array(lfailure,4,'angle-angle',
     &         4,.false.,iatmp)
            if ( lfailure ) return
            call twh_ijaa0(GLB_SET,imolty,iunit,iaa,iatmp(1))
            call twh_ijaa1(GLB_SET,imolty,iunit,iaa,iatmp(2))
            call twh_ijaa2(GLB_SET,imolty,iunit,iaa,iatmp(3))
            call twh_itaa(GLB_SET,imolty,iunit,iaa,iatmp(4))
         enddo

c        --- read in improper torsions
         iatmp(1) = twh_read_labeled_integer(lfailure,4
     &        ,'improper torsion', .false., idname)
         if ( lfailure ) return
         call twh_inimprop(GLB_SET,imolty,iunit,iatmp(1))
         if ( twh_get_inimprop(imolty,iunit) .gt. MAXIMPROP ) then
            write(6,*) 'READTOWHEE: too many improper torsions'
            lfailure = .true.
            return
         endif
         do itor = 1, twh_get_inimprop(imolty,iunit)
            call twh_read_integer_array(lfailure,4,'improper torsion'
     &           ,4,.false.,iatmp)
            if ( lfailure ) return
            call twh_ijimprop2(GLB_SET,imolty,iunit,itor,iatmp(1))
            call twh_ijimprop3(GLB_SET,imolty,iunit,itor,iatmp(2))
            call twh_ijimprop4(GLB_SET,imolty,iunit,itor,iatmp(3))
            call twh_itimprop(GLB_SET,imolty,iunit,itor,iatmp(4))
         enddo
      enddo
      end

      subroutine twh_read_atom_type(lfailure,io_unit,label,ibead
     &     ,ntype,charge)
c     ******************************************************************
c     * A convenience function for reading atom type information       *
c     * Skips all comments.                                            *
c     * Reads two integers and a double precision                      *
c     *                                                                *
c     * Variables:                                                     *
c     * io_unit specifies the file number to read from                 *
c     * label is the label we expect                                   *
c     * ibead, ntype, charge are returned values                       *
c     *                                                                *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(*) label
      logical lfailure
      integer io_unit, ibead, ntype
      double precision charge

c     --- Confirm correct label
      call twh_read_label(lfailure,io_unit,label,.false.
     &     ,'twh_read_atom_type')
      if ( lfailure ) return
      call twh_skip_comments(lfailure,io_unit)
      if ( lfailure ) return
      read(io_unit,*,ERR=100) ibead, ntype, charge
      return

c     --- get here in case of read error
100   write(6, *) "Error reading data for label ",label
      lfailure = .true.
      return
      end

      subroutine twh_read_pdbname(lfailure,io_unit,label,pdbname
     &     ,aminonum,aminoshort)
c     ******************************************************************
c     * A convenience function for reading atom type information       *
c     * Skips all comments.                                            *
c     * Reads two integers and a double precision                      *
c     *                                                                *
c     * Variables:                                                     *
c     * io_unit specifies the file number to read from                 *
c     * label is the label we expect                                   *
c     * pdbname, aminonum, aminoshort are returned values              *
c     *                                                                *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(*) label, pdbname, aminoshort
      logical lfailure
      integer io_unit, aminonum

c     --- Confirm correct label
      call twh_read_label(lfailure,io_unit,label,.false.
     &     ,'twh_read_pdbname')
      if ( lfailure ) return
      call twh_skip_comments(lfailure,io_unit)
      if ( lfailure ) return
      read(io_unit,*,ERR=100) pdbname, aminonum, aminoshort
      return

c     --- get here in case of read error
100   write(6, *) "Error reading data array for label ",label
      lfailure = .true.
      return
      end

      subroutine twh_get_linit(lfailure,io_unit, l_init, print_value, 
     &   idname)
c     ******************************************************************
c     * Returns linit.  This value may come from one of two sources:   *
c     * from the towhee_input file or, if lreadwrapper is set (in a    *
c     * call to towheemainloop) from variables set (by the C wrapper)  *
c     * in the wrapper common block.                                   *
c     *                                                                *
c     * Variables:                                                     *
c     * io_unit specifies the file number to read from                 *
c     * l_init is set by this subroutine.  This local variable will    *
c     *   typically be assigned to the appropriate global variable by  *
c     *   the caller.                                                  *
c     * print_value specifies whether the label and value output to    *
c     * screen.                                                        *
c     * idname is the calling subroutine, printed only in case of error*
c     *                                                                *
c     * originally prior to 2007                                       *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "wrapperF.h"
#define FUNCTION_CHECK_LABEL
#define FUNCTION_READ_LABELED_LOGICAL
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer io_unit
      logical lfailure,l_init,print_value
      character*(*) idname
c     --- local variables
      logical lreadwrapper

      call twh_lreadwrapper(GLB_GET,lreadwrapper)
      if (lreadwrapper) then
c     --- if linit comes from wrapper, then this parameter is allowed,
c     --- but optional, in towhee_input.  See if it exists,  read it
c     --- so that its skipped over.
         if (twh_check_label(lfailure,io_unit,'linit')) then
            l_init = twh_read_labeled_logical(lfailure,io_unit
     &           ,'linit', .false.,idname)
            if ( lfailure ) return
         endif

c     --- actual value comes from wrapper
         l_init = wrap_linit
         if (print_value) then
            write(6,*) 'linit (from wrapper): ',l_init
         endif
      else
         l_init = twh_read_labeled_logical(lfailure,io_unit,'linit'
     &        , print_value,idname)
         if ( lfailure ) return
      endif
      return
      end

      subroutine twh_get_nstep(lfailure,io_unit, num_steps, print_value
     &     , idname)
c     ******************************************************************
c     * Returns num_steps.  This value may come from one of two sources*
c     * from the towhee_input file or,                                 *
c     * if lreadwrapper is set (in a call to towheemainloop) from      *
c     * variables set (by the C wrapper) in the wrapper common block.  *
c     *                                                                *
c     * Variables:                                                     *
c     * io_unit specifies the file number to read from                 *
c     * num_steps is set by this subroutine.  This local variable will *
c     *   typically be assigned to the appropriate global variable by  *
c     *   the caller (nstep).                                          *
c     * print_value specifies whether the label and value output to    *
c     *   the screen.                                                  *
c     * idname is the calling subroutine, printed only in case of error*
c     *                                                                *
c     * originally written prior to 2007                               *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "wrapperF.h"
#define FUNCTION_CHECK_LABEL
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer io_unit,num_steps
      logical lfailure,print_value
      character*(*) idname
c     --- local variables
      logical lreadwrapper

      call twh_lreadwrapper(GLB_GET,lreadwrapper)
      if (lreadwrapper) then
c     --- if num_steps comes from wrapper, then this parameter is allowed,
c     --- but optional, in towhee_input.  See if it exists,  read it
c     --- so that its skipped over.
         if (twh_check_label(lfailure,io_unit,'nstep')) then
            num_steps = twh_read_labeled_integer(lfailure,io_unit
     &           ,'nstep',.false.,idname)
            if ( lfailure ) return
         endif

c     --- actual value comes from wrapper
         num_steps = wrap_nstep
         if (print_value) then
            write(6,*) 'nstep (from wrapper): ',num_steps
         endif
      else
         num_steps = twh_read_labeled_integer(lfailure,io_unit,'nstep'
     &        ,print_value,idname)
         if ( lfailure ) return
      endif
      return
      end

      subroutine twh_get_random_seed(lfailure,io_unit, randseed
     &     , print_value, idname)
c     ******************************************************************
c     * Returns random_seed.  This value may come from one of two      *
c     * sources: from the towhee_input file or,                        *
c     * if lreadwrapper is set (in a call to towheemainloop) from      *
c     * variables set (by the C wrapper) in the wrapper common block.  *
c     *                                                                *
c     * Variables:                                                     *
c     * io_unit specifies the file number to read from randseed is set *
c     * by this subroutine.  This local variable will typically be     *
c     * assigned to the appropriate global variable by the caller.     *
c     * print_value specifies whether the label and value output to    *
c     * screen.                                                        *
c     * idname is the calling subroutine, printed only in case of error*
c     *                                                                *
c     * originally written prior to 2006                               *
c     * last modified 08-16-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "wrapperF.h"
#define FUNCTION_CHECK_LABEL
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer randseed, io_unit
      logical lfailure,print_value
      character*(*) idname
c     --- local variables
      logical lreadwrapper

      call twh_lreadwrapper(GLB_GET,lreadwrapper)
      if (lreadwrapper) then
c     --- if random_seed comes from wrapper this parameter is allowed,
c     --- but optional, in towhee_input.  See if it exists,  read it
c     --- so that its skipped over.
         if (twh_check_label(lfailure,io_unit,'random_seed')) then
            randseed = twh_read_labeled_integer(lfailure,io_unit
     &           ,'random_seed',.false.,idname)
            if ( lfailure ) return
         endif

c     --- actual value comes from wrapper
         randseed = wrap_random_seed
         if (print_value) then
            write(6,*) 'random_seed (from wrapper): ',randseed
         endif
      else
         randseed = twh_read_labeled_integer(lfailure,io_unit
     &        ,'random_seed',print_value, idname)
         if ( lfailure ) return
      endif
      return
      end

      
      subroutine twh_init_random_single(lfailure)
c     ******************************************************************
c     * initializes the random number generator using a single integer *
c     * seed                                                           *
c     *                                                                *
c     * originally written 11-13-2006 by M.G. Martin                   *
c     * last modified 02-05-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_RANDOM
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
      integer*4 iseed
      integer itest,rng_code,luxlevel
      double precision tenmilsum
      double precision rtest
      dimension rtest(5)

c     --- retrieve constants
      call twh_random_code(GLB_GET,rng_code)
      call twh_random_seed(GLB_GET,iseed)
      if ( rng_code .eq. RNG_RANLUX ) then
         call twh_random_luxlevel(GLB_GET,luxlevel)
         write(6,'(a,1x,a,1x,i1,1x,a,1x,i10)')
     &        'Selected random number generator RANLUX with luxury'
     &        ,'level',luxlevel,'and single integer seed',iseed
      elseif ( rng_code .eq. RNG_DEBUG ) then
         write(6,'(a46,1x,a19)')
     &        'Selected random number generator DEBUG that is'
     &        ,'not remotely random'
      elseif ( rng_code .eq. RNG_DX_1597_2_7 ) then
         write(6,'(a,1x,a,1x,i10)')
     &        'Selected random number generator DX-1597-2-7 and single'
     &        ,'integer seed',iseed
      elseif ( rng_code .eq. RNG_KISS99 ) then
         write(6,'(a,1x,a,1x,i10)')
     &        'Selected random number generator KISS99 and single'
     &        ,'integer seed',iseed
      elseif ( rng_code .eq. RNG_MRG32K3A ) then
         write(6,'(a,1x,a,1x,i10)')
     &        'Selected random number generator MRG32k3a and single'
     &        ,'integer seed',iseed
      else
         write(6,*) 'TWH_INIT_RANDOM_SINGLE: tried to initialize'
     &        ,' unknown random number generator with code:',rng_code
         lfailure = .true.
         return
      endif
      call rnginit(rng_code,iseed)
      write(6,'(a,1x,a)')
     &     'Testing random number generator using'
     &     ,'the single integer seed'
c     --- print random numbers
      do itest=1,5
         rtest(itest) = twh_random()
      enddo
      write(6,'(2x,5f10.6)') (rtest(itest),itest=1,5)
c     --- random number test
      tenmilsum = 0.0d0
      do itest = 1,10000000
         tenmilsum = tenmilsum + twh_random()
      enddo
      write(6,'(a,1x,f20.10)')
     &     '10 million RNG sum (approximately 5 million):'
     &     ,tenmilsum

      return
      end

      
      subroutine twh_restart_random(lfailure,rng_code,nrng,irng)
c     ******************************************************************
c     * restarts the random number generator using the full array of   *
c     * integer seeds                                                  *
c     *                                                                *
c     * originally written 08-16-2011 by M.G. Martin                   *
c     * last modified 08-16-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer rng_code,nrng
      integer irng
      dimension irng(nrng)
c     --- local variables
      integer luxlevel

c     --- retrieve constants
      if ( rng_code .eq. RNG_RANLUX ) then
         call twh_random_luxlevel(GLB_GET,luxlevel)
         write(6,'(a,1x,a,1x,i1,1x,a)')
     &        'Restarting random number generator RANLUX with luxury'
     &        ,'level',luxlevel,'from integer seed array'
      elseif ( rng_code .eq. RNG_DX_1597_2_7 ) then
         write(6,'(a,1x,a)')
     &        'Restarting random number generator DX-1597-2-7'
     &        ,'from integer seed array'
      elseif ( rng_code .eq. RNG_KISS99 ) then
         write(6,'(a,1x,a)')
     &        'Restarting random number generator KISS99'
     &        ,'from integer seed array'
      elseif ( rng_code .eq. RNG_MRG32K3A ) then
         write(6,'(a,1x,a)')
     &        'Restarting random number generator MRG32k3a'
     &        ,'from integer seed array'
      else
         write(6,*) 'TWH_RESTART_RANDOM: tried to restart'
     &        ,' unknown random number generator with code:',rng_code
         lfailure = .true.
         return
      endif
      call rngrestart(rng_code,nrng,irng)

      return
      end
