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
      subroutine towheemainloop(outputlocation,outputmode,
     &   int_read_wrapper,pstyle)
c     ******************************************************************
c     * serves only to call the real main loop and to free all of the  *
c     * allocated memory after the simulation is complete              *
c     * rewritten as a shell 12-03-2008 by M.G. Martin                 *
c     * last modified 11-16-2010 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer outputlocation,outputmode,int_read_wrapper,pstyle
c     --- local variables
      logical lreadwrapper

c     --- int_read_wrapper is immediately translated and stored in
c     --- lreadwrapper
c     --- This variable indicates whether certain values 
c     --- (e.g. cmix_lambda_lj) are read from towhee_input, or from
c     --- the wrapper data structure defined in wrapperF/C.h and set
c     --- in the C part of the code.
      if ( int_read_wrapper .eq. LOG_TRUE ) then
         lreadwrapper = .true.
      elseif ( int_read_wrapper .eq. LOG_FALSE ) then
         lreadwrapper = .false.
      else
         write(6,*) 'TOWHEEMAINLOOP: invalid value of int_read_wrapper:'
     &        ,int_read_wrapper
         return
      endif
      call twh_lreadwrapper(GLB_SET,lreadwrapper)
c     --- map the local pstyle to the global parallelstyle
      call twh_parallelstyle(GLB_SET,pstyle)

c     --- start the mainloop
      call  twh_mainloop(outputlocation,outputmode)
      call twh_finalize()

      return
      end

      subroutine twh_mainloop(outputlocation,outputmode)
c     ******************************************************************
c     * calls readinput to get info from towhee_input and starting conf*
c     * then either returns, or calls mainloop_main to control the     *
c     * main calculation                                               *
c     * outputlocation -- used for parallel routines.  The rank of this*
c     *     node, which becomes the index to the directory listing in  *
c     *     towhee_parallel.                                           *
c     * outputmode -- if .ne. 0, loutput = false -> quiet mode         *
c     *                                                                *
c     * split from the original twh_mainloop 06-25-2009 by M.G. Martin *
c     * last modified 03-03-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- debug toggle
#define DEBUG_MAINLOOP 0
c     --- variables passed to/from the subroutine
      integer outputlocation,outputmode
c     --- local variables
      logical lfailure,loutput,lfinish,lredirect_stdout
      integer nchain

#if DEBUG_MAINLOOP
#define DEBUG_WRITE(T) write(6,*)T
#else
#define DEBUG_WRITE(T)
#endif
c     --- initialize failure logical
      lfailure = .false.

c     --- initialize many of the variables and accumulators in towhee
      call twh_initialize(lfailure)
      if ( lfailure ) return

c     --- set output logical
      if ( outputmode .eq. 0 ) then
         loutput = .true.
      else
         loutput = .false.
      endif
      call twh_loutput(GLB_SET,loutput)

c     I/O redirection takes place during parallel runs.  When this 
c     happens,
c     * STDOUT is directed to file instead of screen 
c     (and lredirect_stdout = true)
c     * the 'io_directory' global variable is set, and twh_openfile
c     is asked to prepend this directory to all files it opens.
c     This means that all input (towhee_initial, etc.) and output
c     (towhee_final, etc.) files are written to/read from this directory
c     The exception is that when C wrappers are used (described below),
c     * the location of towhee_input is specified explicitly, and may be
c     different from the io_directory 
c     * wrap_lredirect_stdout controls whether stdout is written to file
c     or screen.  If this variable is true, redirected output gets 
c     written to the file towhee.out in the io_directory.  
c     (for MPI parallel,
c     stdout goes to filename specified in towhee_parallel)
c     Currently, such redirection takes place in two situations:
c     * MPI parallel processing -- configuration takes place via 
c     towhee_parallel file
c     * Replica exchange runs -- configuration takes place mostly in 
c     C wrapper, and data passed here via 'wrapper' data.
c     
      call twh_main_get_IO_redirect_info(lfailure,outputlocation
     &     , lredirect_stdout)
      if ( lfailure ) return
c     --- Output some information about the code
      call twh_loutput(GLB_GET,loutput)
      if ( loutput ) call twh_writeintro()

c     --- read in the input files
      call twh_readinput(lfailure,lfinish)
      if ( lfailure ) return
      if ( lfinish ) then
c        --- we are done, this happens for input styles
c        --- other than Towhee
         if ( lredirect_stdout ) then
c           --- close the general output file
            close(6)
         endif
         return
      endif

c     --- assign all of the averages
      call twh_initaverages(lfailure)
      if ( lfailure ) return
c     --- fetch some useful values
      call twh_nchain(GLB_GET,nchain)

      call twh_mainloop_main(lfailure,lredirect_stdout,nchain)

      return
      end

      subroutine twh_mainloop_main(lfailure,lredirect_stdout,nchain)
c     ******************************************************************
c     * controls the bulk of the Monte Carlo simulation given the      *
c     * information that has now been populated into the data          *
c     * structures (usually via readinput, but potentially by an       *
c     * external program)                                              *
c     * runs over all of the cycles and decides which move to perform  *
c     * writes out final configuration to towhee_final at the end      *
c     *                                                                *
c     * portions originally written prior to 1995                      *
c     * rewritten as twh_mainloop 12-03-2008 by M.G. Martin            *
c     * split from the old twh_mainloop 06-25-2009 by M.G. Martin      *
c     * last modified 08-12-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c#include "wrapperF.h"
#define FUNCTION_GET_IJTOR2
#define FUNCTION_GET_IJTOR3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_IJVIB
#define FUNCTION_GET_INTOR
#define FUNCTION_GET_INVIB
#define FUNCTION_GET_NCMT
#define FUNCTION_GET_RMVOL
#define FUNCTION_LEN_TRIM
#include "functions.h"
c     --- variables passded to/from the subroutine
      logical lfailure,lredirect_stdout
      integer nchain
c     --- energy tolerance threshold
#define ENERGY_THRESH 1.d-2
c     --- local variables
#if DEBUG_MAINLOOP
      logical ldebug_fail,lprint
#endif
c     --- logical scalars
      logical lovrlap,lnone,loutdft,louthist,lreadwrapper
      logical loutdlpoly,loutlammps,tmmc_flag
c     --- logical arrays
      logical lmolmove
      dimension lmolmove(MAXBOX,NTMAX)
c     --- character strings
      character*3 aminoshort
      character*3 target
      dimension target(0:1)
      character*8 stepstyle
      character*14 nconfig_char
      character*30 filename
      character*50 run_name
      character*100 c_filename
c     --- integer scalars
      integer ibox,jbox,nnn,ii,idim,jdim,imolty,jmolty,imove,l
      integer iunit,inercycle,itor,itype,istart,iend,scount,ncount
      integer tot_try,tot_rgn,tot_acp,movevector,ivib,aminonum
      integer ichain,ibin,nentry,nconfig,itmax,hist_nequil
      integer istress,atemp,ntemp,stemp,ensemble
      integer itrial,ientry,incount,nmtype,nstep,printfreq
      integer total_nmols, i_tmp,nmolty,ivalue,jvalue,nunit,ncmt
      integer moviefreq,trmaxdispfreq,blocksize,backupfreq
      integer volmaxdispfreq,pressurefreq,pdboutfreq,restartfreq
      integer scp_num_foreign_lambda,ipair,npair,numboxes,nvmove
      integer histcalcfreq,histdumpfreq,pressurestyle
      integer two_bond_bias_update_frequency,three_bond_bias_update_freq
      integer n_tmmc_min,n_tmmc_max,c_element,c_matrix_freq,jj
      integer macro_variable,frequency,keycode
#if DEBUG_MAINLOOP
      integer cp_number
#endif
c     --- integer arrays
      integer nminp,nmaxp
      dimension nminp(NTMAX)
      dimension nmaxp(NTMAX)
      integer ndist
      dimension ndist(0:nchain,NTMAX)
      integer ncmt_list
      dimension ncmt_list(NDUMPHIST,NTMAX)
      integer chaininfo
      dimension chaininfo(1,3)
c     --- double precision scalars
      double precision v,vinter,vtail,vintra,vvib,vbend,vtg
      double precision vext,velect,vsolvation,vdum,vquantum
      double precision dvinter,dvelect,datemp,dntemp
      double precision pscb1,pscb2,ratvol,ratcell,dvalue
      double precision pspiv2,pscr,pscrb,binsize,pspiv
      double precision vhist,value,rmcell,p_acc
      double precision cm_ii_one,cm_ii_zero,cm_ii_negone
      double precision cm_jj_one,cm_jj_zero,cm_jj_negone
#if DEBUG_MAINLOOP
      double precision vtest,density,dtemp
#endif
c     --- double precision arrays
      double precision currentpressure
      dimension currentpressure(MAXBOX)
      double precision stress
      dimension stress(0:7)
      double precision eng_list
      dimension eng_list(NDUMPHIST)
      double precision speden
      dimension speden(MAXBOX)
      double precision tempdim
      dimension tempdim(3)
      double precision vforeign
      dimension vforeign(MAX_FOREIGN_LAMBDA)


c     --- retrieve constants now that they have been set by the 
c     --- input files
      call twh_numboxes(GLB_GET,numboxes)
      call twh_nvmove(GLB_GET,nvmove)
      call twh_ensemble(GLB_GET,ensemble)
      call twh_two_bond_bias_update_frequency(GLB_GET
     &     ,two_bond_bias_update_frequency)
      call twh_three_bond_bias_update_freq(GLB_GET
     &     ,three_bond_bias_update_freq)
      call twh_tmmc_flag(GLB_GET,tmmc_flag)

c     --- initialize currentpressure
      do ibox = 1,numboxes
         currentpressure(ibox) = 0.0d0
      enddo
c     --- initialize local accumulators for grand-canonical ensemble
      nentry = 0
      nconfig = 0
      do itmax = 1, NTMAX
         nminp(itmax) = nchain + 1
         nmaxp(itmax) = -1
         do ichain = 0,nchain
            ndist(ichain,itmax) = 0
         enddo
      enddo
c     --- initialize stress tensor
      do istress = 0,7
         stress(istress) = 0.0d0
      enddo


c     --- calculate initial energy and check for overlaps
      do ibox = 1, numboxes
         write(6,'(a24,1x,i5)') 'Initial Energies for Box',ibox
         call twh_engtotal(lfailure,lovrlap, v, vinter,vdum,vtail
     &        ,vintra,vvib,vbend,vtg,vext,velect,vsolvation,vquantum
     &        ,vforeign,dvinter,dvelect,ibox, .false.,.false.,.true.
     &        ,.true.)
         if ( lfailure ) return
         if( lovrlap ) then
            write(6,*) 'MAINLOOP: overlap in initial configuration '
            return
         endif
         call twh_boxvclassic(GLB_SET,ibox,v)
         call twh_boxvinter(GLB_SET,ibox,vinter)
         call twh_boxvtail(GLB_SET,ibox,vtail)
         call twh_boxvintra(GLB_SET,ibox,vintra)
         call twh_boxvbond(GLB_SET,ibox,vvib)
         call twh_boxvangle(GLB_SET,ibox,vbend)
         call twh_boxvtorsion(GLB_SET,ibox,vtg)
         call twh_boxvfield(GLB_SET,ibox,vext)
         call twh_boxvelect(GLB_SET,ibox,velect)
         call twh_boxvsolvation(GLB_SET,ibox,vsolvation)
         call twh_boxvquantum(GLB_SET,ibox,vquantum)
         call scp_set_derivative(ibox, LAMBDA_LJ, dvinter)
         call scp_set_derivative(ibox, LAMBDA_C, dvelect)
         call scp_set_foreign_energies(ibox, vforeign)
      enddo

      call twh_pressurefreq(GLB_GET,pressurefreq)
      call twh_pressurestyle(GLB_GET,pressurestyle)
      if ( pressurefreq .ne. 0 ) then
c        --- calculate initial stress tensor
         if ( pressurestyle .eq. PRES_VIRIAL ) then
            do ibox = 1, numboxes
               call twh_stresstensor(lfailure,stress,ibox)
               if ( lfailure ) return
               currentpressure(ibox) = stress(0)
               write(6,'(1x,a30,1x,i2,1x,a1,1x,f14.2)') 
     &              'initial virial pressure in box',ibox,'=',stress(0)
            enddo
         elseif ( pressurestyle .eq. PRES_RADIAL ) then
            do ibox = 1,numboxes
               call twh_radialpressure(lfailure,stress(0),ibox)
               if ( lfailure ) return
               write(6,'(1x,a30,1x,i2,1x,a1,1x,f14.2)') 
     &              'initial radial pressure in box',ibox,'=',stress(0)
            enddo
         endif
      endif

      call twh_pdb_output_freq(GLB_GET,pdboutfreq)
      if ( pdboutfreq .ne. 0 ) then
c        --- output initial snapshot
         nnn = 0
         call twh_writepdb(lfailure,nnn)
         if ( lfailure ) return
      endif

      write(6,*)
      write(6,*) '+++++ start of markov chain +++++'
      write(6,*)

c     write the column headings of runtime output. (in writeruntime.F)
      call twh_writeruntime_header()
c     --- retrieve some useful constants
      call twh_nmolty(GLB_GET,nmolty)

c     --- Begin main Monte Carlo loop over all cycles and all molecules 
      call twh_stepstyle(GLB_GET,stepstyle)
      if ( stepstyle .eq. 'cycles' ) then
c        --- we want to perform nchain (a full cycle) moves per step
         inercycle = nchain
      elseif ( stepstyle .eq. 'moves' ) then
c        --- we only want to do one move per step
         inercycle = 1
      else
         write(6,*) 'MAINLOOP: unknown stepstyle'
         return
      endif

      call twh_nstep(GLB_GET,nstep)
      call twh_printfreq(GLB_GET,printfreq)
      call twh_moviefreq(GLB_GET,moviefreq)
      call twh_trmaxdispfreq(GLB_GET,trmaxdispfreq)
      call twh_blocksize(GLB_GET,blocksize)
      call twh_backupfreq(GLB_GET,backupfreq)
      call twh_restartfreq(GLB_GET,restartfreq)
      call twh_volmaxdispfreq(GLB_GET,volmaxdispfreq)
      do nnn = 1, nstep
         do ii = 1, inercycle 

c           Perform one monte carlo move
            call main_do_MC_move(lfailure)
            if ( lfailure ) return

c           --- accumulate averages
            call twh_accumulate_averages(lfailure,speden)
            if ( lfailure ) return
#if DEBUG_MAINLOOP       
c           --- debug stuff
            ldebug_fail = .false.
            call twh_cp_number(GLB_GET,cp_number)
            if ( cp_number .eq. CP_EMBEDDED_ATOM_METHOD ) then
c              --- store the eam_rho
               do ichain = 1,nchain
                  call twh_eam_rho(GLB_STORE,CRD_R_TO_T,ichain,density)
               enddo
            endif
            do ibox = 1,numboxes
c              --- toggle lprint for engtotal output
               lprint = .true.
c              --- compute the total energies
               call twh_engtotal(lfailure,lovrlap,v,vinter,vdum,vtail
     &              ,vintra
     &              ,vvib,vbend,vtg,vext,velect,vsolvation,vquantum
     &              ,vforeign,dvinter,dvelect,ibox,.false.,.false.
     &              ,lprint,.true.)
               if ( lfailure ) return
               call twh_boxvclassic(GLB_GET,ibox,vtest)
               if ( (v - vtest)**2 .gt. ENERGY_THRESH ) then
                  write(6,*) 'ibox ',ibox
                  write(6,*) 'v, boxvclassic:',v,vtest,v-vtest
                  call twh_boxvinter(GLB_GET,ibox,vtest)
                  write(6,*) 'vinter,boxvinter',vinter,vtest
     &                 ,vinter-vtest
                  call twh_boxvtail(GLB_GET,ibox,vtest)
                  write(6,*) 'vtail,boxvtail',vtail,vtest
     &                 ,vtail-vtest
                  call twh_boxvintra(GLB_GET,ibox,vtest)
                  write(6,*) 'vintra,vintrab',vintra,vtest
     &                 ,vintra-vtest
                  call twh_boxvbond(GLB_GET,ibox,vtest)
                  write(6,*) 'vvib,vvibb',vvib,vtest,vvib-vtest
                  call twh_boxvangle(GLB_GET,ibox,vtest)
                  write(6,*) 'vbend,boxvangle',vbend,vtest
     &                 ,vbend-vtest
                  call twh_boxvtorsion(GLB_GET,ibox,vtest)
                  write(6,*) 'vtg,boxvtorsion',vtg,vtest,vtg-vtest
                  call twh_boxvfield(GLB_GET,ibox,vtest)
                  write(6,*) 'vext,boxvfield',vext,vtest,vext-vtest
                  call twh_boxvelect(GLB_GET,ibox,vtest)
                  write(6,*) 'velect,boxvelect ',velect,vtest
     &                 ,velect-vtest
                  call twh_boxvsolvation(GLB_GET,ibox,vtest)
                  write(6,*) 'vsolvation,boxvsolvation',vsolvation
     &                 ,vtest,vsolvation-vtest
                  ldebug_fail = .true.
               endif
            enddo
c           --- check the eam_rho
            if ( cp_number .eq. CP_EMBEDDED_ATOM_METHOD ) then
               do ichain = 1,nchain
                  call twh_eam_rho(GLB_GET,CRD_REAL,ichain,density)
                  call twh_eam_rho(GLB_GET,CRD_TEMP,ichain,dtemp)
                  if ( ( dtemp - density)**2
     &                 .gt. ENERGY_THRESH ) then
                     write(6,*) 'mismatch eam_rho ichain',ichain
     &                    ,density,dtemp
                     call twh_nboxi(GLB_GET,ichain,jbox)
                     write(6,*) 'nboxi(ichain)',jbox
                     ldebug_fail = .true.
                  endif
               enddo
            endif
            if ( ldebug_fail ) then
               write(6,*) 'debug_mainloop error detected'
               return
            endif
c           --- end of debug stuff
#endif
         enddo
c        --- end of a single full Monte Carlo cycle
c
c        --- TMMC Modification:
c        --- These calculations are not needed when doing
c        --- tmmc, and only slow down the simulation.
c        --- Rg calculation might be useful to include down the road
         if ( .NOT. tmmc_flag ) then
c           --- perform additional chemical potential measurements
            do ibox = 1,numboxes
               do imolty = 1,nmolty
                  call twh_chempotperstep(GLB_GET,imolty,ivalue)
                  do itrial = 1,ivalue
                     call twh_swapmoves(lfailure,SWAP_CHEMPOT_BOX
     &                    ,imolty,ibox)
                     if ( lfailure ) return
                  enddo
               enddo
            enddo
c           --- perform computation of isolated molecule chemical
c           --- potential
            ibox = numboxes+1
            do imolty = 1,nmolty
               call twh_nmolectyp(GLB_GET,imolty,nmtype)
               if ( nmtype .eq. 0 ) then
c                 --- also compute the isolated molecule chempot
                  call twh_chempotperstep(GLB_GET,imolty,ivalue)
                  do itrial = 1,ivalue
                     call twh_swapmoves(lfailure,SWAP_CHEMPOT_ISO
     &                    ,imolty,ibox)
                     if ( lfailure ) return
                  enddo
               endif
            enddo

c           --- accumulate average radius of gyration.  Note, this is 
c           --- only done once per cycle instead of after every move as
c           --- this calculation is slightly more expensive.
            keycode = AC_RADIUS_OF_GYRATION
            do ichain = 1,nchain
               call twh_moltyp(GLB_GET,ichain,imolty)
               call twh_nboxi(GLB_GET,ichain,ibox)
               if ( ibox .ne. 0 ) then
                  call twh_gyration(GLB_GET,ichain,value)
                  incount = 1
                  call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox
     &                 ,imolty,value,incount)
                  if ( lfailure ) return
               endif
            enddo
         endif

c        --- histogram/tmmc data processing
         call twh_louthist(GLB_GET,louthist)
         if ( ensemble .eq. ENS_muVT .and. louthist) then
            call twh_boxvinter(GLB_GET,1,vinter)
            call twh_boxvintra(GLB_GET,1,vintra)
            call twh_boxvbond(GLB_GET,1,vvib)
            call twh_boxvangle(GLB_GET,1,vbend)
            call twh_boxvtorsion(GLB_GET,1,vtg)
            call twh_boxvelect(GLB_GET,1,velect)
            vhist = vinter + velect + vvib + vbend + vintra + vtg
            nconfig = nconfig + 1
c           --- Main tmmc stuff takes place within this IF construct
            IF (tmmc_flag) THEN
c              --- retrieve tmmc variables
               call twh_run_name(GLB_GET,run_name)
               call twh_n_tmmc_min(GLB_GET,n_tmmc_min)
               call twh_n_tmmc_max(GLB_GET,n_tmmc_max)
               total_nmols = 0
               DO i_tmp = 1, nmolty
                  call twh_ncmt(GLB_GET,1,i_tmp,ncmt)
                  total_nmols = total_nmols + ncmt
               END DO
c              --- Update C-matrix
               call twh_c_element(GLB_GET,c_element)
               call twh_p_acc(GLB_GET,p_acc)
               call twh_macro_variable(GLB_GET,macro_variable)
               call twh_c_matrix(GLB_INCR,macro_variable,c_element
     &              ,p_acc)
               dvalue = 1.0D0 - p_acc
               call twh_c_matrix(GLB_INCR,macro_variable,0,dvalue)
c              --- Update semigrand accumulator
c              --- For now, it's the potential energy
               call twh_v_semigrand(GLB_INCR,total_nmols,vhist)
c              --- Update tmmc_weight if necessary
               call twh_weight_freq(GLB_GET,frequency)
               if ( mod(nconfig, frequency) .eq. 0) then
                  call twh_tmmc_weight(GLB_SET,n_tmmc_min,0.0D0)
                  do ii = (n_tmmc_min+1), n_tmmc_max
                     jj = ii - 1
                     call twh_c_matrix(GLB_GET,jj,1,cm_jj_one)
                     call twh_c_matrix(GLB_GET,ii,-1,cm_ii_negone)
                     if ( ( cm_jj_one .gt. 0.0d0 ) .and. 
     &                    ( cm_ii_negone .gt. 0.0d0) ) then
                        call twh_c_matrix(GLB_GET,jj,0,cm_jj_zero)
                        call twh_c_matrix(GLB_GET,jj,-1,cm_jj_negone)
                        call twh_c_matrix(GLB_GET,ii,1,cm_ii_one)
                        call twh_c_matrix(GLB_GET,ii,0,cm_ii_zero)
                        call twh_tmmc_weight(GLB_GET,jj,dvalue)
                        dvalue = dvalue + dlog( cm_jj_one / (cm_jj_one +
     &                       cm_jj_zero + cm_jj_negone ) ) -
     &                       dlog( cm_ii_negone / ( cm_ii_one +
     &                       cm_ii_zero + cm_ii_negone ) )
                        call twh_tmmc_weight(GLB_SET,ii,dvalue)
                     else
                        call twh_tmmc_weight(GLB_GET,jj,dvalue)
                        call twh_tmmc_weight(GLB_SET,ii,dvalue)
                     endif
                  enddo
                  
                  filename = run_name(1:twh_len_trim(run_name))//
     &                 '.tmmc_weights.dat'
                  call twh_openfile(lfailure,51,filename,'DEFAULT'
     &                 ,'UNKNOWN','DEFAULT')
                  if ( lfailure ) return
                  DO ientry = n_tmmc_min, n_tmmc_max
                     call twh_tmmc_weight(GLB_GET,ientry,dvalue)
                     write(51,'(I5,2X,E16.8)') ientry,dvalue
                  END DO
                  close(51)
               end if
c              --- Write out c_matrix and v_semigrand if necessary
               call twh_c_matrix_freq(GLB_GET,c_matrix_freq)
               IF ( MOD(nconfig,c_matrix_freq).eq.0 ) THEN
                  WRITE(nconfig_char,'(I14.14)') (nconfig/c_matrix_freq)
                  c_filename=run_name(1:twh_len_trim(run_name))//
     &                 '.c.'//nconfig_char//'.dat'
                  call twh_openfile(lfailure,51,c_filename,'DEFAULT'
     &                 ,'UNKNOWN','DEFAULT')
                  if ( lfailure ) return
                  do ientry = n_tmmc_min, n_tmmc_max
                     call twh_c_matrix(GLB_GET,ientry,-1,cm_ii_negone)
                     call twh_c_matrix(GLB_GET,ientry,0,cm_ii_zero)
                     call twh_c_matrix(GLB_GET,ientry,1,cm_ii_one)
                     write(51,'(I5,E16.8,2X,E16.8,2X,E16.8)') ientry
     &                    ,cm_ii_negone,cm_ii_zero,cm_ii_one
                  enddo
                  close(51)
            
                  c_filename=run_name(1:twh_len_trim(run_name))//
     &                 '.vsg.'//nconfig_char//'.dat'
                  call twh_openfile(lfailure,51,c_filename,'DEFAULT'
     &                 ,'UNKNOWN','DEFAULT')
                  if ( lfailure ) return
                  DO ientry = n_tmmc_min, n_tmmc_max
                     call twh_v_semigrand(GLB_GET,ientry,dvalue)
                     WRITE(51,'(I5,2X,E16.8)') ientry,dvalue
                  END DO
                  close(51)
               endif
            endif
c           --- Some of this histogram  stuff is not needed when doing
c           --- GC-TMMC, so they have been cut by virtue of the
c           --- tmmc_flag variable
            if ( .NOT. tmmc_flag) then
               if (mod(nconfig, printfreq) .eq. 0) then
                  call twh_file_flct(GLB_GET,filename)
                  call twh_openfile(lfailure,50,filename,'DEFAULT','OLD'
     &                 ,'APPEND')
                  if ( lfailure ) return
                  write(50,* ) nconfig
     &                 ,(twh_get_ncmt(1,imolty), imolty=1,nmolty),vhist
                  close(50)
               endif
               call twh_hist_nequil(GLB_GET,hist_nequil)
               call twh_histcalcfreq(GLB_GET,histcalcfreq)
               if ( nnn .gt. hist_nequil .and. 
     &              mod(nconfig,histcalcfreq).eq. 0 ) then
c                 --- gather up histogram information
                  do imolty = 1,nmolty
                     nminp(imolty) = min(nminp(imolty)
     &                    ,twh_get_ncmt(1,imolty))
                     nmaxp(imolty) = max(nmaxp(imolty)
     &                    ,twh_get_ncmt(1,imolty))
                  enddo
                  nentry = nentry + 1
                  do imolty = 1,nmolty
                     call twh_ncmt(GLB_GET,1,imolty,ncmt)
                     ncmt_list(nentry,imolty) = ncmt
                     ndist(ncmt,imolty) = ndist(ncmt,imolty) + 1
                  enddo
                  eng_list(nentry) = vhist
               endif
               call twh_histdumpfreq(GLB_GET,histdumpfreq)
               if ( mod(nconfig,histdumpfreq) .eq. 0 ) then
c                 --- output histogram information to a file
                  call twh_file_his(GLB_GET,filename)
                  call twh_openfile(lfailure,51,filename,'DEFAULT','OLD'
     &                 ,'APPEND')
                  if ( lfailure ) return
                  do ientry = 1,nentry
                     write(51,*) 
     &                    (ncmt_list(ientry,imolty),imolty=1,nmolty),
     &                    eng_list(ientry)
                  enddo
                  close(51)
                  nentry = 0
                  do imolty = 1,nmolty
                     call twh_file_ndis(GLB_GET,imolty,filename)
                     call twh_openfile(lfailure,52,filename,'DEFAULT'
     &                    ,'UNKNOWN','DEFAULT')
                     if ( lfailure ) return
                     do ientry = nminp(imolty),nmaxp(imolty)
                        write(52,*) ientry,ndist(ientry,imolty)
                     enddo
                     close(52)
                  enddo
               endif
            endif
         endif

c        --- test for output of conformation to pdb file
         if ( pdboutfreq .ne. 0 ) then
            if ( mod(nnn,pdboutfreq) .eq. 0 ) then
               call twh_writepdb(lfailure,nnn)
               if ( lfailure ) return
            endif
         endif
         
c        --- test to see if we need to perform periodic operations 
c        --- test for computation of pressure
         if ( pressurefreq .ne. 0 ) then
            if (  mod(nnn,pressurefreq) .eq. 0 .and.
     &           ( pressurestyle .eq. PRES_VIRIAL
     &           .or. pressurestyle .eq. PRES_RADIAL ) ) then
               call twh_evaluate_pressure(lfailure,speden
     &              ,currentpressure,stress)
               if ( lfailure ) return
            endif
         endif
c        --- see if any periodic operations need to be performed
c        --- note that for many of these, freq=0 (implying, never) is 
c        --- permitted, and the mod operator treats these cases 
c        --- correctly
         if (trmaxdispfreq .ne. 0) then
            if ( mod(nnn,trmaxdispfreq) .eq. 0 ) 
     &           call twh_updatetrmaxdisp()
         endif

         if (volmaxdispfreq .ne. 0) then
            if ( mod(nnn,volmaxdispfreq) .eq. 0 )
     &           call twh_updatevolmaxdisp()
         endif

         if (printfreq .ne. 0) then
            if ( mod(nnn,printfreq) .eq. 0 ) 
     &        call twh_writeruntime(currentpressure,nnn)
         endif

         if (moviefreq .ne. 0) then
            if ( mod(nnn,moviefreq) .eq. 0 ) then
               call twh_writemovie(lfailure,nnn)
               if ( lfailure ) return
            endif
         endif

         if (backupfreq .ne. 0) then
            if ( mod(nnn,backupfreq) .eq. 0 ) then
c           --- write out the restart configurations to SAFETY-file 
               call twh_rwconf(lfailure,nchain,WRITE_FLAG,nnn
     &              ,RWCONF_BACKUP)
               if ( lfailure ) then
                  write(6,*) 'MAINLOOP: problem writing backup'
                  return
               endif
            endif
         endif

         if (restartfreq .ne. 0) then
            if ( mod(nnn,restartfreq) .eq. 0 ) then
c              --- write out the restart configurations to the restart
c              --- file
               call twh_rwconf(lfailure,nchain,WRITE_FLAG,nnn
     &              ,RWCONF_RESTART)
               if ( lfailure ) then
                  write(6,*) 'MAINLOOP: problem writing restart'
                  return
               endif
            endif
         endif

         if (blocksize .ne. 0) then
            if ( mod(nnn,blocksize) .eq. 0 ) call twh_pupdate(lfailure)
            if ( lfailure ) return
         endif

         call twh_updatetwobondbias(nnn)
         call twh_updatethreebondbias(nnn)
      enddo

c     --- end of the loop over the cycles  
      close(51)

      write(6,*)
      write(6,*) '+++++ end of markov chain +++++'
      write(6,*)
 
c     --- output hmatrix information
      write(6,*) 'Final hmatrix (general box dimensions) '
      do ibox = 1,numboxes
         write(6,'(a4,1x,i5)') 'Box:',ibox
         do idim = 1,3
            do jdim = 1,3
               call twh_hmatrix(GLB_GET,ibox,idim,jdim,tempdim(jdim))
            enddo
            write(6,126) idim,(tempdim(jdim),jdim=1,3)
         enddo
      enddo
 126  format(2x,'hmatrix(',i1,',x)',3f14.5)

c     --- determine whether any 3D volume moves were performed
      lnone = .true.
      do imove = 1, nvmove
c        --- add the last batch of volume attempts into the total
         call twh_bnvol(GLB_GET,imove,ntemp)
         call twh_acnvol(GLB_INCR,imove,ntemp)
c        --- add the last batch of volume accepts into the total
         call twh_bavol(GLB_GET,imove,atemp)
         call twh_acsvol(GLB_INCR,imove,atemp)
c        --- check to see if we have any attempted moves
         call twh_acnvol(GLB_GET,imove,ncount)
         if ( ncount .ne. 0 ) lnone = .false.
      enddo

      if ( .not. lnone ) then
c        --- output information on volume move 
         write(6,*)
         write(6,'(a26)') '* 3D Volume Change Moves *'
         do imove = 1, nvmove
            call twh_acnvol(GLB_GET,imove,ncount)
            if ( ncount .ne. 0 ) then
               call twh_acsvol(GLB_GET,imove,scount)
               ratvol = dble(scount) / dble(ncount)
            else
               ratvol = 0.0d0
            endif
            if ( ensemble .eq. ENS_NPT ) then
               write(6,40) imove,ncount,ratvol,twh_get_rmvol(imove)
            else
               call twh_acsvol(GLB_GET,imove,scount)
               call twh_pairbox(GLB_GET,imove,1,ivalue)
               call twh_pairbox(GLB_GET,imove,2,jvalue)
               write(6,41) ivalue,jvalue,ncount,scount,ratvol
     &              ,twh_get_rmvol(imove)
            endif
         enddo
      endif
 40   format(' Box ',i2,' Tries: ',i8,' Acp. Ratio: ',f6.3,
     &     ' Max Disp.: ',e9.3)
 41   format(' Box ',i2,' and ',i2,' Tries: ',i8
     &     ,' Accepted: ',i8,' Acp. Ratio: ',f6.3,' Max Disp.: ',e9.3)

c     --- determine whether any unit cell perturbation moves were tried
      lnone = .true.
      do imove = 1,nvmove
         do movevector = 1,9
c           --- add the last batch of attempted moves into the total
            call twh_bncell(GLB_GET,imove,movevector,dntemp)
            call twh_acncell(GLB_INCR,imove,movevector,dntemp)
c           --- add the last batch of accepted moves into the total
            call twh_bacell(GLB_GET,imove,movevector,datemp)
            call twh_acscell(GLB_INCR,imove,movevector,datemp)
c           --- see if we have any attempts
            call twh_acncell(GLB_GET,imove,movevector,dntemp)
            if ( dntemp .gt. 0.5d0 ) lnone = .false.
         enddo
      enddo

      if ( .not. lnone ) then
c        --- output information on unit cell perturbation move 
         write(6,*)
         write(6,'(a32)') '* Unit Cell Perturbation Moves *'
         do imove = 1, nvmove
            do idim = 1,3
               do jdim = 1,3
                  movevector = 3*(idim-1) + jdim
                  call twh_acncell(GLB_GET,imove,movevector,dntemp)
                  if ( dntemp .gt. 0.5d0 ) then
                     call twh_acscell(GLB_GET,imove,movevector,datemp)
                     ratcell = datemp / dntemp
                  else
                     ratcell = 0.0d0
                  endif
                  if ( ensemble .eq. ENS_NPT ) then
                     call twh_rmcell(GLB_GET,imove,idim,jdim,rmcell)
                     write(6,44) imove,idim,jdim,dntemp,ratcell
     &                    ,rmcell
                  else
                     call twh_pairbox(GLB_GET,imove,1,ivalue)
                     call twh_pairbox(GLB_GET,imove,2,jvalue)
                     call twh_rmcell(GLB_GET,imove,idim,jdim,rmcell)
                     write(6,45) ivalue,jvalue,idim,jdim,dntemp
     &                    ,ratcell,rmcell
                  endif
               enddo
            enddo
         enddo
      endif
 44   format(' Box ',i2,' H-matrix(',i1,',',i1,') Tries: ',f12.0
     &     ,' Acp. Ratio: ',f6.3,' Max Disp.: ',e9.3)
 45   format(' Box',i2,' and ',i2,' H-matrix(',i1,',',i1,') Tries: '
     &     ,f12.0,' Acp. Ratio: ',f6.3,' Max Disp.: ',e9.3)

c     --- determine if any rbswaps were attempted
      lnone = .true.
      do imolty = 1, nmolty
         do ibox = 1, numboxes
            do jbox = 1,numboxes
               call twh_bnrbswap(GLB_GET,imolty,ibox,jbox,ntemp)
               if ( ntemp .ne. 0 ) lnone = .false.
            enddo
         enddo
      enddo

      if ( .not. lnone ) then
c        --- report on rbswap performance
         write(6,*)
         write(6,'(a30)') '* Rotational-Bias SWAP Moves *'
         do imolty = 1, nmolty
            write(6,*) 'Molecule type:',imolty
            do ibox = 1, numboxes
               do jbox = 1,numboxes
                  call twh_bnrbswap(GLB_GET,imolty,ibox,jbox,ntemp)
                  if ( ntemp .ne. 0 ) then
                     call twh_barbswap(GLB_GET,imolty,ibox,jbox,atemp)
                     write(6,50) jbox,ibox,ntemp,atemp
                  endif
               enddo
            enddo
         enddo
 50      format('  From box ',i2,' to box ',i2,' Attempted:',i8
     &        ,' Accepted:',i8)
      endif

c     --- determine if any cbswaps were attempted
      lnone = .true.
      do imolty = 1, nmolty
         do ibox = 1, numboxes
            do jbox = 1,numboxes
               call twh_bncbswap(GLB_GET,imolty,ibox,jbox,ntemp)
               if ( ntemp .ne. 0 ) lnone = .false.
            enddo
         enddo
      enddo

      if ( .not. lnone ) then
c        --- report on cbswap performance
         write(6,*)
         write(6,'(a35)') '* Configurational-Bias SWAP Moves *'
         do imolty = 1, nmolty
            write(6,'(a14,1x,i5)') 'Molecule type:',imolty
            do ibox = 1, numboxes
               do jbox = 1,numboxes
                  call twh_bncbswap(GLB_GET,imolty,ibox,jbox,ntemp)
                  if ( ntemp .ne. 0 ) then
                     call twh_bacbswap(GLB_GET,imolty,ibox,jbox,atemp)
                     call twh_bscbswap(GLB_GET,imolty,ibox,jbox,stemp)
                     write(6,51) jbox,ibox,ntemp,stemp,atemp
                  endif
               enddo
            enddo
         enddo
 51      format('  From box ',i2,' to box ',i2,' Attempted:',i8
     &        ,' Grown:',i8,' Accepted:',i8)
      endif

c     --- switch moves
      npair = (nmolty*(nmolty-1))/2
c     --- determine if any SWITCH moves were attempted
      lnone = .true.
      do ipair = 1,npair
         do ibox = 1,numboxes
            call twh_acnswitch(GLB_GET,ipair,ibox,dntemp)
            if ( dntemp .gt. 0.5d0 ) lnone = .false.
         enddo
      enddo
c     --- output performance of SWITCH moves
      if ( .not. lnone ) then
         write(6,*)
         write(6,'(a25)') '* SWITCH Moves *'
         do ibox = 1,numboxes
            do ipair = 1,npair
               call twh_acnswitch(GLB_GET,ipair,ibox,dntemp)
               if ( dntemp .gt. 0.5d0 ) then
                  call twh_acsswitch(GLB_GET,ipair,ibox,datemp)
                  write(6,'(a16,i3,1x,a4,i2,1x,a9,1x,f15.0
     &                 ,1x,a9,1x,f15.0,1x,a9,1x,f7.3,1x,a1)') 
     &                 'Molecule pair:',ipair,'Box:',ibox
     &                 ,'Attempts:',dntemp
     &                 ,'Accepted:',datemp
     &                 ,'Accepted:',100.0d0*datemp/dntemp,'%'
               endif
            enddo
         enddo
      endif

c     --- determine if any avb1 moves were performed
      lnone = .true.
      do ibox = 1, numboxes
         do imolty = 1, nmolty
            do jmolty = 1,nmolty
               call twh_bnavb1(GLB_GET,imolty,jmolty,ibox,-1,ntemp)
               if ( ntemp .ne. 0 ) lnone = .false.
            enddo
         enddo
      enddo

      if ( .not. lnone ) then
c        --- write some information about cbregrow performance 
         write(6,*)
         write(6,'(a35)') '* Aggregation-Volume-Bias-1 Moves *'
         target(0) = ' in'
         target(1) = 'out'
         do ibox = 1, numboxes
            do imolty = 1, nmolty
               do jmolty = 1,nmolty
                  call twh_bnavb1(GLB_GET,imolty,jmolty,ibox,-1,ntemp)
                  if ( ntemp .ne. 0 ) then
                     write(6,*) 'Box Moved Around From To  Attempted'
     &                    ,'  Accepted  Percent'
                     do istart = 0,1
                        do iend = 0,1
                           itype = 2*istart + iend
                           call twh_bnavb1(GLB_GET,imolty,jmolty,ibox
     &                          ,itype,stemp)
                           if ( stemp .ne. 0 ) then
                              call twh_baavb1(GLB_GET,imolty,jmolty,ibox
     &                             ,itype,atemp)
                              write(6,'(1x,i3,1x,i5,1x,i6,2x,a3,1x,a3
     &                             ,1x,i9,2x,i8,2x,f7.2)') 
     &                             ibox,imolty,jmolty
     &                             ,target(istart),target(iend)
     &                             ,stemp,atemp
     &                             ,1d2*dble(atemp)/dble(stemp)
                           endif
                        enddo
                     enddo
c                    --- output number not found
                     itype = -2
                     call twh_bnavb1(GLB_GET,imolty,jmolty,ibox,itype
     &                    ,stemp)
                     if ( stemp .ne. 0 ) then
                        write(6,'(1x,i3,1x,i5,1x,i6,1x
     &                       ,a23,1x,i9,1x,a6)') 
     &                       ibox,imolty,jmolty
     &                       ,'Eligable pair not found'
     &                       ,stemp
     &                       ,'times.'
                     endif
                  endif
               enddo
            enddo
         enddo
      endif

c     --- determine if any avb2 moves were performed
      lnone = .true.
      do ibox = 1, numboxes
         do imolty = 1, nmolty
            do jmolty = 1,nmolty
               call twh_bnavb2(GLB_GET,imolty,jmolty,ibox,-1,ntemp)
               if ( ntemp .ne. 0 ) lnone = .false.
            enddo
         enddo
      enddo

      if ( .not. lnone ) then
c        --- write some information about cbregrow performance 
         write(6,*)
         write(6,'(a34)') '* Aggregation-Volume-Bias-2 Moves *'
         target(0) = ' in'
         target(1) = 'out'
         do ibox = 1, numboxes
            do imolty = 1, nmolty
               do jmolty = 1,nmolty
                  call twh_bnavb2(GLB_GET,imolty,jmolty,ibox,-1,ntemp)
                  if ( ntemp .ne. 0 ) then
                     write(6,*) 'Box Moved Around From To  Attempted'
     &                    ,'  Accepted  Percent'
                     do istart = 0,1
                        do iend = 0,1
                           itype = 2*istart + iend
                           call twh_bnavb2(GLB_GET,imolty,jmolty,ibox
     &                          ,itype,stemp)
                           if ( stemp .ne. 0 ) then
                              call twh_baavb2(GLB_GET,imolty,jmolty
     &                             ,ibox,itype,atemp)
                              write(6,'(1x,i3,1x,i5,1x,i6,2x,a3,1x,a3
     &                             ,1x,i9,2x,i8,2x,f7.2)') 
     &                             ibox,imolty,jmolty
     &                             ,target(istart),target(iend)
     &                             ,stemp,atemp
     &                             ,1d2*dble(atemp)/dble(stemp)
     &                             
                           endif
                        enddo
                     enddo
c                    --- output number not found
                     itype = -2
                     call twh_bnavb2(GLB_GET,imolty,jmolty,ibox,itype
     &                    ,stemp)
                     if ( stemp .ne. 0 ) then
                        write(6,'(1x,i3,1x,i5,1x,i6,1x
     &                       ,a23,1x,i9,1x,a6)') 
     &                       ibox,imolty,jmolty
     &                       ,'Eligable pair not found'
     &                       ,stemp
     &                       ,'times.'
                     endif
                  endif
               enddo
            enddo
         enddo
      endif

c     --- determine if any avb3 moves were performed
      lnone = .true.
      do ibox = 1, numboxes
         do imolty = 1, nmolty
            do jmolty = 1,nmolty
               call twh_bnavb3(GLB_GET,imolty,jmolty,ibox,ntemp)
               if ( ntemp .ne. 0 ) lnone = .false.
            enddo
         enddo
      enddo

      if ( .not. lnone ) then
c        --- write some information about avb3 performance 
         write(6,*)
         write(6,'(a34)') '* Aggregation-Volume-Bias-3 Moves *'
         do ibox = 1, numboxes
            do imolty = 1, nmolty
               do jmolty = 1,nmolty
                  call twh_bnavb3(GLB_GET,imolty,jmolty,ibox,ntemp)
                  if ( ntemp .ne. 0 ) then
                     call twh_baavb3(GLB_GET,imolty,jmolty,ibox,atemp)
                     write(6,*) 'Box Moved Around Attempted'
     &                    ,'  Accepted  Percent'
                     write(6,'(1x,i3,1x,i5,1x,i5
     &                    ,1x,i9,2x,i8,2x,f7.2)') 
     &                    ibox,imolty,jmolty,ntemp,atemp
     &                    ,1d2*dble(atemp)/dble(ntemp)
                  endif
               enddo
            enddo
         enddo
      endif

c     --- determine if any grand canonical ensemble swap moves were 
c     --- performed
      lnone = .true.
      do imolty = 1,nmolty
         do ibox = 0,numboxes
            do jbox = 0,numboxes
c              --- only grand canonical if one of the boxes is box zero
               if ( ibox*jbox .eq. 0 ) then
                  call twh_bncbswap(GLB_GET,imolty,ibox,jbox,ntemp)
                  if ( ntemp .ne. 0 ) lnone = .false.
               endif
            enddo
         enddo
      enddo

      if ( .not. lnone ) then
c        --- report on uvtcb performance
         write(6,*)
         write(6,'(a39)') '* Grand Canonical Ensemble SWAP Moves *'
         do imolty = 1, nmolty
            write(6,'(a14,1x,i5)') 'Molecule type:',imolty
            do ibox = 0, numboxes
               do jbox = 0,numboxes
c                 --- only grand canonical if one of the boxes is 0
                  if ( ibox*jbox .eq. 0 ) then
                     call twh_bncbswap(GLB_GET,imolty,ibox,jbox,ntemp)
                     if ( ntemp .ne. 0 ) then
                        call twh_bacbswap(GLB_GET,imolty,ibox,jbox
     &                       ,atemp)
                        call twh_bscbswap(GLB_GET,imolty,ibox,jbox
     &                       ,stemp)
                        write(6,51) jbox,ibox,ntemp,stemp,atemp
                     endif
                  endif
               enddo
            enddo
         enddo
      endif

c     --- determine if any regular cbregrow moves were performed
      lnone = .true.
      do ibox = 1, numboxes
         do imolty = 1, nmolty
            call twh_nunit(GLB_GET,imolty,nunit)
            lmolmove(ibox,imolty) = .false.
            do iunit = 1, nunit
               call twh_bncb(GLB_GET,imolty,iunit,ibox,ntemp)
               if ( ntemp .ne. 0 ) then
                  lnone = .false.
                  lmolmove(ibox,imolty) = .true.
               endif
            enddo
         enddo
      enddo

      if ( .not. lnone ) then
c        --- write some information about cbregrow performance 
         write(6,*)
         write(6,'(a39)') '* Configurational-Bias REGROWTH Moves *'
         do ibox = 1, numboxes
            do imolty = 1, nmolty
               if ( lmolmove(ibox,imolty) ) then
                  write(6,'(a15,i5,1x,a4,i5)') 'Molecule type:',imolty
     &                 ,'Box:',ibox
                  write(6,*) '    Length  Attempts  Regrown'
     &                 ,'  Accepted  %Regrown   %Accep.'
                  call twh_nunit(GLB_GET,imolty,nunit)
                  do iunit = 1, nunit
                     call twh_bncb(GLB_GET,imolty,iunit,ibox,ntemp)
                     if ( ntemp .ne. 0 ) then
                        call twh_bacb(GLB_GET,imolty,iunit,ibox,atemp)
                        call twh_bscb(GLB_GET,imolty,iunit,ibox,stemp)
                        pscb1 = dble(stemp)*100.0d0/dble(ntemp)
                        pscb2 = dble(atemp)*100.0d0/dble(ntemp)
                        write(6,'(i9,3i10,2f10.2)') iunit
     &                       ,ntemp,stemp,atemp,pscb1,pscb2
                     endif
                  enddo
               endif
            enddo
         enddo
      endif

c     --- determine if any backbone cbregrow moves were performed
      lnone = .true.
      do ibox = 1, numboxes
         do imolty = 1, nmolty
            lmolmove(ibox,imolty) = .false.
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1, nunit
               do itor = 1,twh_get_intor(imolty,iunit)
                  call twh_bnback(GLB_GET,imolty,iunit,itor,ibox,ntemp)
                  if ( ntemp .ne. 0 ) then
                     lnone = .false.
                     lmolmove(ibox,imolty) = .true.
                  endif
               enddo
            enddo
         enddo
      enddo

      if ( .not. lnone ) then
c        --- write some information about backbone cbregrow performance 
         write(6,*)
         write(6,'(a48)') 
     &        '* Configurational-Bias BACKBONE REGROWTH Moves *'
         do ibox = 1, numboxes
            do imolty = 1, nmolty
               if ( lmolmove(ibox,imolty) ) then
                  write(6,'(a14,i5,1x,a4,i5)') 'Molecule type:',imolty
     &                 ,'Box:',ibox
                  write(6,*) 'Backbone Atoms    Attempts  Regrown'
     &                 ,'  Accepted  %Regrown   %Accep.'
                  tot_try = 0
                  tot_rgn = 0
                  tot_acp = 0
                  call twh_nunit(GLB_GET,imolty,nunit)
                  do iunit = 1, nunit
                     do itor = 1,twh_get_intor(imolty,iunit)
                        call twh_bnback(GLB_GET,imolty,iunit,itor,ibox
     &                       ,ntemp)
                        tot_try = tot_try + ntemp
                        if ( ntemp .ne. 0 ) then
                           call twh_bsback(GLB_GET,imolty,iunit,itor
     &                          ,ibox,stemp)
                           tot_rgn = tot_rgn + stemp 
                           pscb1 = dble(stemp)*1d2/dble(ntemp)
                           call twh_baback(GLB_GET,imolty,iunit,itor
     &                          ,ibox,atemp)
                           tot_acp = tot_acp + atemp 
                           pscb2 = dble(atemp)*1d2/dble(ntemp)
                           write(6,'(4i4,3i10,2f10.2)') iunit
     &                          ,twh_get_ijtor2(imolty,iunit,itor)
     &                          ,twh_get_ijtor3(imolty,iunit,itor)
     &                          ,twh_get_ijtor4(imolty,iunit,itor)
     &                          ,ntemp,stemp,atemp,pscb1,pscb2
                        endif
                     enddo
                  enddo
                  if ( tot_try .ne. 0 ) then
                     pscb1 = 1d2*dble(tot_rgn)/dble(tot_try)
                     pscb2 = 1d2*dble(tot_acp)/dble(tot_try)
                     write(6,'(a16,3i10,2f10.2)') 'totals:'
     &                    ,tot_try,tot_rgn,tot_acp,pscb1,pscb2
                  endif
               endif
            enddo
         enddo
      endif

c     --- determine if any cb sidechain moves were performed
      lnone = .true.
      do ibox = 1, numboxes
         do imolty = 1, nmolty
            lmolmove(ibox,imolty) = .false.
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1, nunit
               call twh_bncbside(GLB_GET,imolty,iunit,ibox,ntemp)
               if ( ntemp .ne. 0 ) then
                  lnone = .false.
                  lmolmove(ibox,imolty) = .true.
               endif
            enddo
         enddo
      enddo

      if ( .not. lnone ) then
c        --- write some information about cb side-chain performance 
         write(6,*)
         write(6,'(a33)') '* CB Side-Chain Regrowth Moves *'
         do ibox = 1, numboxes
            do imolty = 1, nmolty
               tot_try = 0
               tot_rgn = 0
               tot_acp = 0
               if ( lmolmove(ibox,imolty) ) then
                  write(6,'(a15,i5,1x,a4,i5)') 'Molecule type:',imolty
     &                 ,'Box:',ibox
                  write(6,*) '    Peptide  Attempts  Regrown'
     &                 ,'  Accepted  %Regrown   %Accep.'
                  call twh_nunit(GLB_GET,imolty,nunit)
                  do iunit = 1, nunit
                     call twh_bncbside(GLB_GET,imolty,iunit,ibox,ntemp)
                     if ( ntemp .ne. 0 ) then
                        call twh_bacbside(GLB_GET,imolty,iunit,ibox
     &                       ,atemp)
                        call twh_bscbside(GLB_GET,imolty,iunit,ibox
     &                       ,stemp)
                        tot_try = tot_try + ntemp
                        tot_rgn = tot_rgn + stemp
                        tot_acp = tot_acp + atemp
                        pscb1 = dble(stemp)*100.0d0/dble(ntemp)
                        pscb2 = dble(atemp)*100.0d0/dble(ntemp)
                        call twh_aminoshort(GLB_GET,imolty,iunit
     &                       ,aminoshort)
                        call twh_aminonum(GLB_GET,imolty,iunit,aminonum)
                        write(6,'(1x,i6,a1,a3,3i10,2f10.2)') 
     &                       aminonum,'-',aminoshort
     &                       ,ntemp,stemp,atemp,pscb1,pscb2
                     endif
                  enddo
                  if ( tot_try .ne. 0 ) then
                     pscb1 = 1d2*dble(tot_rgn)/dble(tot_try)
                     pscb2 = 1d2*dble(tot_acp)/dble(tot_try)
                     write(6,'(5x,a5,1x,3i10,2f10.2)') 'Total'
     &                    ,tot_try,tot_rgn,tot_acp,pscb1,pscb2
                  endif
               endif
            enddo
         enddo
      endif

c     --- determine if any pivot moves were performed
      lnone = .true.
      do ibox = 1,numboxes
         do imolty = 1, nmolty
            lmolmove(ibox,imolty) = .false.
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1, nunit
               do ivib = 1, twh_get_invib(imolty,iunit)
                  call twh_bnpivot(GLB_GET,imolty,iunit,ivib,ibox,ntemp)
                  if ( ntemp .ne. 0) then
                     lnone = .false.
                     lmolmove(ibox,imolty) = .true.
                  endif
               enddo
            enddo
         enddo
      enddo

      if ( .not. lnone) then
c        --- write some information about the pivot move performance
         write(6,*)
         write(6,'(a25)') '* Pivot Torsional Moves *'
         do ibox = 1, numboxes
            do imolty = 1, nmolty          
               if (lmolmove(ibox,imolty)) then
                  write(6,*) 'Molecule type: ',imolty,' Box: ',ibox
                  write(6,*) 'Torsional Atoms  Attempts  Accepted'
     &                 ,'  % Accepted'
                  tot_try = 0
                  tot_acp = 0
                  call twh_nunit(GLB_GET,imolty,nunit)
                  do iunit = 1, nunit
                     do ivib = 1, twh_get_invib(imolty,iunit)
                        call twh_bnpivot(GLB_GET,imolty,iunit,ivib,ibox
     &                       ,ntemp)
                        tot_try = tot_try + ntemp
                        if ( ntemp .ne. 0) then
                           call twh_bapivot(GLB_GET,imolty,iunit,ivib
     &                          ,ibox,atemp)
                           tot_acp = tot_acp + atemp
                           pspiv = dble(atemp) * 1d2/dble(ntemp)
                           write(6,'(2i4,2i10,f10.2)') iunit,
     &                          twh_get_ijvib(imolty,iunit,ivib),
     &                          ntemp,atemp,pspiv
                        endif
                     enddo
                  enddo
                  if (tot_try .ne. 0) then
                     pspiv = 1d2*dble(tot_acp)/dble(tot_try)
                     write(6, '(a7,1x,2i10,f10.2)') 'totals:',
     &                    tot_try,tot_acp,pspiv
                  endif
               endif
            enddo
         enddo

c        --- write some information about the pivot move performance 
c        --- (breakdown of size of move by increments of 30 degrees 
         write(6,'(a35)') '* Pivot breakdown by size of move *'
         binsize = 180.0d0/dble(MAXPIVOTBIN)
         do ibox = 1, numboxes
            write(6,*) 'Box: ',ibox
            write(6,*) 'Size of Move  Attempts  Accepted'
     &                 ,'  % Accepted'
            tot_try = 0
            tot_acp = 0
            do ibin = 1, MAXPIVOTBIN
               call twh_bnpivbin(GLB_GET,ibox,ibin,ntemp)
               tot_try = tot_try + ntemp
               if ( ntemp .ne. 0) then
                  call twh_bapivbin(GLB_GET,ibox,ibin,atemp)
                  tot_acp = tot_acp + atemp
                  pspiv2 = dble(atemp) * 1d2/dble(ntemp)
                  write(6,'(2f8.2,2i10,f10.2)') ((ibin-1)*binsize),
     &                  (ibin*binsize),ntemp,atemp,pspiv2
               endif
            enddo
            if (tot_try .ne. 0) then
               pspiv2 = 1d2*dble(tot_acp)/dble(tot_try)
               write(6, '(a16,2i10,f10.2)') 'totals:',
     &              tot_try,tot_acp,pspiv
            endif
         enddo
      endif   

c     --- determine if any concerted rotation moves were performed
      lnone = .true.
      do ibox = 1,numboxes
         do imolty = 1, nmolty
            lmolmove(ibox,imolty) = .false.
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1, nunit
               call twh_bnconrot(GLB_GET,imolty,iunit,ibox,ntemp)
               if ( ntemp .ne. 0) then
                  lnone = .false.
                  lmolmove(ibox,imolty) = .true.
               endif
            enddo
         enddo
      enddo

      if ( .not. lnone) then
c        --- write some information about the concerted rotation move 
c        --- performance
         write(6,*)
         write(6,'(a28)') '* Concerted Rotation Moves *'
         do ibox = 1, numboxes
            do imolty = 1, nmolty          
               if (lmolmove(ibox,imolty)) then
                  write(6,*) 'Molecule type: ',imolty,' Box: ',ibox
                  write(6,*) 'Atom  Attempts  Accepted'
     &                 ,'  % Accepted'
                  tot_try = 0
                  tot_acp = 0
                  call twh_nunit(GLB_GET,imolty,nunit)
                  do iunit = 1, nunit
                     call twh_bnconrot(GLB_GET,imolty,iunit,ibox,ntemp)
                     tot_try = tot_try + ntemp
                     if ( ntemp .ne. 0) then
                        call twh_baconrot(GLB_GET,imolty,iunit,ibox
     &                       ,atemp)
                        tot_acp = tot_acp + atemp
                        pscr = 
     &                      dble(atemp) * 1d2/dble(ntemp)
                        write(6,'(i4,2i10,f10.2)') iunit,
     &                          ntemp,atemp,pscr
                     endif
                  enddo
                  if (tot_try .ne. 0) then
                     pscr = 1d2*dble(tot_acp)/dble(tot_try)
                     write(6, '(a16,2i10,f10.2)') 'totals:',
     &                    tot_try,tot_acp,pscr
                  endif
               endif
            enddo
         enddo
      endif   

c     --- determine if any backbone concerted rotation moves were 
c     --- performed
      lnone = .true.
      do ibox = 1,numboxes
         do imolty = 1, nmolty
            lmolmove(ibox,imolty) = .false.
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1, nunit
               call twh_bncrback(GLB_GET,imolty,iunit,ibox,ntemp)
               if ( ntemp .ne. 0 ) then
                  lnone = .false.
                  lmolmove(ibox,imolty) = .true.
               endif
            enddo
         enddo
      enddo

      if ( .not. lnone) then
c        --- write some information about the backbone concerted 
c        --- rotation move performance
         write(6,*)
         write(6,'(a37)') '* Backbone Concerted Rotation Moves *'
         do ibox = 1, numboxes
            do imolty = 1, nmolty          
               if (lmolmove(ibox,imolty)) then
                  write(6,*) 'Molecule type: ',imolty,' Box: ',ibox
                  write(6,*) 'Atom  Attempts  Accepted'
     &                 ,'  % Accepted'
                  tot_try = 0
                  tot_acp = 0
                  call twh_nunit(GLB_GET,imolty,nunit)
                  do iunit = 1, nunit
                     call twh_bncrback(GLB_GET,imolty,iunit,ibox,ntemp)
                     tot_try = tot_try + ntemp
                     if ( ntemp .ne. 0) then
                        call twh_bacrback(GLB_GET,imolty,iunit,ibox
     &                       ,atemp)
                         tot_acp = tot_acp + atemp
                         pscrb = dble(atemp)*1d2/dble(ntemp)
                         write(6,'(i4,2i10,f10.2)') iunit,
     &                          ntemp,atemp,pscrb
                     endif
                  enddo
                  if (tot_try .ne. 0) then
                     pscrb = 1d2*dble(tot_acp)/dble(tot_try)
                     write(6, '(a16,2i10,f10.2)') 'totals:',
     &                    tot_try,tot_acp,pscrb
                  endif
               endif
            enddo
         enddo
      endif   

c     --- determine if any planeshift moves were attempted
      lnone = .true.
      do ibox = 1,numboxes
         do idim = 1,3
            call twh_bnshift(GLB_GET,1,ibox,idim,ntemp)
            if ( ntemp .ne. 0 ) lnone = .false.
         enddo
      enddo
c     --- output performance of planeshift moves
      if ( .not. lnone ) then
         write(6,*)
         write(6,'(a21)') '* Plane Shift Moves *'
         do ibox = 1,numboxes
            do idim = 1,3
               call twh_bnshift(GLB_GET,1,ibox,idim,ntemp)
               call twh_bashift(GLB_GET,1,ibox,idim,atemp)
               write(6,192) ibox,idim,ntemp,atemp
 192           format('Box ',i2,1x,'Perpendicular Dimension '
     &              ,i1,1x,'Attempts '
     &              ,i10,1x,'Accepted ',i10)
            enddo
         enddo
      endif

c     --- determine if any row shift moves were attempted
      lnone = .true.
      do ibox = 1,numboxes
         do idim = 1,3
            call twh_bnshift(GLB_GET,2,ibox,idim,ntemp)
            if ( ntemp .ne. 0 ) lnone = .false.
         enddo
      enddo
c     --- output performance of row shift moves
      if ( .not. lnone ) then
         write(6,*)
         write(6,'(a19)') '* Row Shift Moves *'
         do ibox = 1,numboxes
            do idim = 1,3
               call twh_bnshift(GLB_GET,2,ibox,idim,ntemp)
               call twh_bashift(GLB_GET,2,ibox,idim,atemp)
               write(6,193) ibox,idim,ntemp,atemp
 193           format('Box ',i2,1x,'Row Dimension ',i1,1x,'Attempts '
     &              ,i10,1x,'Accepted ',i10)
            enddo
         enddo
      endif

c     --- add the temporary atom translation move information into the
c     --- total just in case we did not end evenly on an update
      do imolty = 1,nmolty
         do ibox = 1,numboxes
c           --- increment total number of moves attempted
            call twh_bntraa(GLB_GET,imolty,ibox,dntemp)
            call twh_acntraa(GLB_INCR,imolty,ibox,dntemp)
c           --- increment total number of moves accepted
            call twh_batraa(GLB_GET,imolty,ibox,datemp)
            call twh_acstraa(GLB_INCR,imolty,ibox,datemp)
         enddo
      enddo
c     --- determine if any single atom tranlation moves were attempted
      lnone = .true.
      do imolty = 1,nmolty
         do ibox = 1,numboxes
            call twh_acntraa(GLB_GET,imolty,ibox,dntemp)
            if ( dntemp .gt. 0.5d0 ) lnone = .false.
         enddo
      enddo
c     --- output performance of atom translation moves
      if ( .not. lnone ) then
         write(6,*) 
         write(6,'(a26)') '* Atom Translation Moves *'
         do ibox = 1,numboxes
            do imolty = 1,nmolty
               call twh_acntraa(GLB_GET,imolty,ibox,dntemp)
               if ( dntemp .gt. 0.5d0 ) then
                  call twh_acstraa(GLB_GET,imolty,ibox,datemp)
                  write(6,'(a9,i3,1x,a4,i2,1x,a9,1x,f15.0
     &                 ,1x,a9,1x,f15.0,1x,a9,1x,f7.3,1x,a1)') 
     &                 'Molecule:',imolty,'Box:',ibox
     &                 ,'Attempts:',dntemp
     &                 ,'Accepted:',datemp
     &                 ,'Accepted:',100.0d0*datemp/dntemp,'%'
               endif
            enddo
         enddo
      endif
         
c     --- add the temporary COM translation move information into the
c     --- total just in case we did not end evenly on an update
      do imolty = 1,nmolty
         do ibox = 1,numboxes
c           --- increment total attempted moves
            call twh_bntrac(GLB_GET,imolty,ibox,dntemp)
            call twh_acntrac(GLB_INCR,imolty,ibox,dntemp)
c           --- increment total accepted moves
            call twh_batrac(GLB_GET,imolty,ibox,datemp)
            call twh_acstrac(GLB_INCR,imolty,ibox,datemp)
         enddo
      enddo
c     --- determine if any COM tranlation moves were attempted
      lnone = .true.
      do imolty = 1,nmolty
         do ibox = 1,numboxes
            call twh_acntrac(GLB_GET,imolty,ibox,dntemp)
            if ( dntemp .gt. 0.5d0 ) lnone = .false.
         enddo
      enddo
c     --- output performance of COM translation moves
      if ( .not. lnone ) then
         write(6,*)
         write(6,'(a25)') '* COM Translation Moves *'
         do ibox = 1,numboxes
            do imolty = 1,nmolty
               call twh_acntrac(GLB_GET,imolty,ibox,dntemp)
               if ( dntemp .gt. 0.5d0 ) then
                  call twh_acstrac(GLB_GET,imolty,ibox,datemp)
                  write(6,'(a9,i3,1x,a4,i2,1x,a9,1x,f15.0
     &                 ,1x,a9,1x,f15.0,1x,a9,1x,f7.3,1x,a1)') 
     &                 'Molecule:',imolty,'Box:',ibox
     &                 ,'Attempts:',dntemp
     &                 ,'Accepted:',datemp
     &                 ,'Accepted:',100.0d0*datemp/dntemp,'%'
               endif
            enddo
         enddo
      endif
         
c     --- add the temporary rotation move information into the
c     --- total just in case we did not end evenly on an update
      do imolty = 1,nmolty
         do ibox = 1,numboxes
c           --- increment total attempted moves
            call twh_bnrot(GLB_GET,imolty,ibox,dntemp)
            call twh_acnrot(GLB_INCR,imolty,ibox,dntemp)
c           --- increment total accepted moves
            call twh_barot(GLB_GET,imolty,ibox,datemp)
            call twh_acsrot(GLB_INCR,imolty,ibox,datemp)
         enddo
      enddo
c     --- determine if any rotation moves were attempted
      lnone = .true.
      do imolty = 1,nmolty
         do ibox = 1,numboxes
            call twh_acnrot(GLB_GET,imolty,ibox,dntemp)
            if ( dntemp .gt. 0.5d0 ) lnone = .false.
         enddo
      enddo
c     --- output performance of rotation moves
      if ( .not. lnone ) then
         write(6,*)
         write(6,'(a18)') '* Rotation Moves *'
         do ibox = 1,numboxes
            do imolty = 1,nmolty
               call twh_acnrot(GLB_GET,imolty,ibox,dntemp)
               if ( dntemp .gt. 0.5d0 ) then
                  call twh_acsrot(GLB_GET,imolty,ibox,datemp)
                  write(6,'(a9,i3,1x,a4,i2,1x,a9,1x,f15.0
     &                 ,1x,a9,1x,f15.0,1x,a9,1x,f7.3,1x,a1)') 
     &                 'Molecule:',imolty,'Box:',ibox
     &                 ,'Attempts:',dntemp
     &                 ,'Accepted:',datemp
     &                 ,'Accepted:',100.0d0*datemp/dntemp,'%'
               endif
            enddo
         enddo
      endif

      call print_composite_move_summary()
         
c     --- write out the final configuration from the run to towhee_final
      call twh_rwconf(lfailure,nchain,WRITE_FLAG,nstep,RWCONF_FINAL)
      if ( lfailure ) then
         write(6,*) 'problem writing final configuration'
         return
      endif
      call twh_check_energies(lfailure)
      if ( lfailure ) return
      call twh_writeaverages(lfailure)
      if ( lfailure ) return
c     --- write foreign energies to the wrapper structure, to make it
c     --- visible to the enclosing C wrapper, if necessary.  Energies
c     --- in box 1 only.
      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
      call twh_lreadwrapper(GLB_GET,lreadwrapper)
      if (scp_num_foreign_lambda .gt. 0 .and. lreadwrapper) then
         call scp_get_foreign_energies(1, vforeign)
         do l = 1,scp_num_foreign_lambda
            call twh_wrap_foreign_energy(GLB_SET,l,vforeign(l))
         enddo
      endif

c     --- close the movie file if it was opened
      if (moviefreq .ne. 0) then
         if ( moviefreq .le. nstep ) close(10)
      endif

c     --- output the final configuration as a pdb file
      if ( pdboutfreq .ne. 0  ) then
         call twh_writepdb(lfailure,nstep)
         if ( lfailure ) return
      endif
c     --- output the citations information
      call twh_citations(lfailure)
      if ( lfailure ) return
         
c     --- output the final configuration as a tramonto file
      call twh_loutdft(GLB_GET,loutdft)
      if ( loutdft ) then
         do ibox = 1,numboxes
            call twh_writetramonto(lfailure,CRD_REAL,ibox,0,chaininfo)
            if ( lfailure ) then
               write(6,*) 'MAINLOOP: Problem writing final DFT files'
               return
            endif
         enddo
      endif
c     --- output the final configuration as a lammps file
      call twh_loutlammps(GLB_GET,loutlammps)
      if ( loutlammps ) then
         call twh_writelammps(lfailure)
         if ( lfailure ) then
            write(6,*) 'MAINLOOP: problem writing final lammps file'
            return
         endif
      endif
c     --- output the final configuration as a dlpoly file
      call twh_loutdlpoly(GLB_GET,loutdlpoly)
      if ( loutdlpoly ) then
         call twh_writedlpoly(lfailure)
         if ( lfailure ) then
            write(6,*) 'MAINLOOP: problem writing final dlpoly file'
            return
         endif
      endif
c     --- output the final fixed endpoint SAFE parameters
      if ( two_bond_bias_update_frequency .ne. 0 .or.
     &     three_bond_bias_update_freq .ne. 0 ) then
         call twh_rwsafe(lfailure,WRITE_FLAG)
         if ( lfailure ) return
      endif

      if ( lredirect_stdout ) then
c        --- close the general output file
         close(6)
      endif
      return
      end

      subroutine twh_evaluate_pressure(lfailure,speden,currentpressure
     &     ,stress)
c     ******************************************************************
c     * Evaluate and store pressure terms                              *
c     *                                                                *
c     * Variables modified:                                            *
c     *  speden -- specific density for each box                       *
c     *  currentpressure -- eponymous, per box                         *
c     *  stress -- the stress tensor                                   *
c     * Split from mainloop 04-06-2006 by MAW                          *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision speden
      dimension speden(MAXBOX)
      double precision currentpressure
      dimension currentpressure(MAXBOX)
      double precision stress
      dimension stress(0:7)
      logical lfailure
c     --- local variables
c     --- integer scalars
      integer ibox,imolty,ncmt_vap,ncmt_liq,numboxes
      integer istart,iend,istress,nmolty,keycode
      integer incount,ivap,iliq,pressurestyle
c     --- double precision scalars
      double precision value,boxvolume,volliq,volvap
c     --- double precision arrays
      double precision uterm
      dimension uterm(MAXBOX)
      double precision pvterm
      dimension pvterm(MAXBOX)

c     --- retrieve constants
      call twh_numboxes(GLB_GET,numboxes)
      call twh_pressurestyle(GLB_GET,pressurestyle)

      do ibox = 1, numboxes
         if ( pressurestyle .eq. PRES_VIRIAL ) then
            istart = 0
            iend = 7
            call twh_stresstensor(lfailure,stress,ibox)
            if ( lfailure ) return
         elseif ( pressurestyle .eq. PRES_RADIAL ) then
            istart = 0
            iend = 0
            call twh_radialpressure(lfailure,stress(0),ibox)
            if ( lfailure ) return
         else
            write(6,*) 'MAINLOOP: unknown pressurestyle'
            write(6,*) pressurestyle
            lfailure = .true.
            return
         endif
c        --- set the current pressure to this pressure
         currentpressure(ibox) = stress(0)
c        --- accumulate pressure averages
         imolty = 0
         incount = 1
         do istress = istart,iend
            keycode = AC_PRESSURE - istress
            value = stress(istress)
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &               ,value,incount)
            if ( lfailure ) return
         enddo
c        --- transform potential energy into kJ/mol
         call twh_boxvclassic(GLB_GET,ibox,value)
         uterm(ibox) = value*8.314d0/1000.0d0
c        --- pressure*volume
         call twh_boxvolume(GLB_GET,ibox,boxvolume)
         pvterm(ibox) = 6.02d23*currentpressure(ibox)*boxvolume/1.0d30
         keycode = AC_PV
         value = pvterm(ibox)
         call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &            ,value,incount)
         if ( lfailure ) return
c        --- enthalpy
         value = uterm(ibox) + pvterm(ibox)
         keycode = AC_ENTHALPY
         call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &        ,value,incount)
               if ( lfailure ) return
c        --- enthalpy squared
         value = value**2
         keycode = AC_ENTHALPY_SQ
         call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &            ,value,incount)
         if ( lfailure ) return
      enddo
      call twh_nmolty(GLB_GET,nmolty)
      if ( numboxes .eq. 2 .and. nmolty .eq. 1 ) then
         ibox = 0
         imolty = 1
c                 --- determine the vapor and liquid boxes
         if ( speden(1) .gt. speden(2) ) then
            ivap = 2
            iliq = 1
         else
            ivap = 1
            iliq = 2
         endif
         call twh_ncmt(GLB_GET,ivap,imolty,ncmt_vap)
         call twh_ncmt(GLB_GET,iliq,imolty,ncmt_liq)
         if ( ncmt_vap * ncmt_liq .ne. 0 ) then
c                    --- heat of vaporization (Direct)
            keycode = AC_H_VAP_DIRECT
            value = ((uterm(ivap) + pvterm(ivap)) / dble(ncmt_vap)) - 
     &               ( (uterm(iliq) + pvterm(iliq)) / dble(ncmt_liq))
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &               ,value,incount)
            if ( lfailure ) return
c           --- heat of vaporization (vapor p)
            keycode = AC_H_VAP_VAPOR_P
            call twh_boxvolume(GLB_GET,ivap,volvap)
            call twh_boxvolume(GLB_GET,iliq,volliq)
            value = (uterm(ivap)/dble(ncmt_vap))
     &           - (uterm(iliq)/dble(ncmt_liq)) 
     &           + (6.02d23/1.0d30)*currentpressure(ivap)
     &           *( (volvap/dble(ncmt_vap)) - (volliq/dble(ncmt_liq)) )
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,imolty
     &               ,value,incount)
            if ( lfailure ) return
         endif
      endif

      return
      end

      subroutine main_do_MC_move(lfailure)
c     ******************************************************************
c     * Perform one Monte Carlo move.                                  *
c     *                                                                *
c     * Variables:                                                     *
c     *   lfailure -- set to true if problem                           *
c     * last modified 10-14-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_RANDOM
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
      logical tmmc_flag
      double precision rm,pmtest
c     --- the values of these parameters do not matter here.
      integer imolty,ibox,i_tmp,nmolty,ncmt,ensemble

c     --- retrieve constants
      call twh_ensemble(GLB_GET,ensemble)
      call twh_tmmc_flag(GLB_GET,tmmc_flag)
c     --- Initialize tmmc variables
      if (tmmc_flag) then
         call twh_macro_variable(GLB_SET,0)
         call twh_nmolty(GLB_GET,nmolty)
         do i_tmp = 1, nmolty
            call twh_ncmt(GLB_GET,1,i_tmp,ncmt)
            call twh_macro_variable(GLB_INCR,ncmt)
         enddo
         call twh_c_element(GLB_SET,0)
         call twh_p_acc(GLB_SET,0.0d0)
      END IF
c     --- determine which move to perform according to the pm**
      rm = twh_random()
#if DEBUG_MAINLOOP
      write(6,*) 'rm ',rm
#endif

c     --- isotropic volume move
      call twh_pmvol(GLB_GET,pmtest)
      if  ( rm .le. pmtest ) then
         if ( ensemble .eq. ENS_NPT ) then
c           --- isobaric-isothermal volume move
            DEBUG_WRITE('MOVE: isobaric-isothermal volume move')
            call twh_volnpt(lfailure,1)
         else
c           --- canonical Gibbs ensemble volume exchange
            DEBUG_WRITE('MOVE: canonical Gibbs ens vol ex.')
            call twh_volnvt(lfailure,1)
         endif
         return
      endif

c     --- unit cell perturbation move
      call twh_pmcell(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         if ( ensemble .eq. ENS_NPT ) then
            DEBUG_WRITE('MOVE: unit cell perturbation NPT')
            call twh_volnpt(lfailure,2)
         else
            DEBUG_WRITE('MOVE: unit cell perturbation NVT')
            call twh_volnvt(lfailure,2)
         endif
         return
      endif

c     --- inter-box rotational-bias swap move
      call twh_pm2boxrbswap(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         DEBUG_WRITE('MOVE: inter-box rot-bias swap')
c        --- imolty and ibox are dummy variables to twh_swapmoves here.
         call twh_swapmoves(lfailure,SWAP_RB_INTER,imolty,ibox)
         return
      endif

c     --- inter-box configurational-bias swap move
      call twh_pm2boxcbswap(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         DEBUG_WRITE('MOVE: inter-box conf-bias swap')
         call twh_swapmoves(lfailure,SWAP_CB_INTER,imolty,ibox)
         return
      endif

c     --- insertion/removal move used in grand canonical ensemble
      call twh_pmuvtcbswap(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         DEBUG_WRITE('MOVE: insert/remove move')
         call twh_swapmoves(lfailure,SWAP_CB_GC,imolty,ibox)
         return
      endif

c     --- intra-box configurational-bias cbswap move 
      call twh_pm1boxcbswap(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         DEBUG_WRITE('MOVE: intra-box conf-bias cbswap')
         call twh_swapmoves(lfailure,SWAP_CB_INTRA,imolty,ibox)
         return
      endif

c     --- intra-box center-of-mass molecule switch move
      call twh_pm1boxcomswitch(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         DEBUG_WRITE('MOVE: intra-box center-of-mass switch')
         call twh_switch(lfailure)
         return
      endif

c     --- aggregation volume bias move 1
c     --- B. Chen, J.I. Siepmann "A Novel Monte Carlo Algorithm
c     --- for Simulating Strongly Associating Fluids:
c     --- Applications to Water, Hydrogen Fluoride, and Acetic
c     --- Acid", J. Phys. Chem. B 104, 8725-8734 (2000).
      call twh_pmavb1(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         DEBUG_WRITE('MOVE: aggregation volume bias move 1')
         call twh_swapmoves(lfailure,SWAP_AVB_ONE,imolty,ibox)
         return
      endif

c     --- aggregation volume bias move 2
c     --- B. Chen; J.I. Siepmann "Improving the Efficiency 
c     --- of the Aggregation-Volume-Bias Monte Carlo Algorithm"
c     --- J. Phys. Chem. B 105, 11275-11282 (2001).
      call twh_pmavb2(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         DEBUG_WRITE('MOVE: aggregation volume bias move 2')
         call twh_swapmoves(lfailure,SWAP_AVB_TWO,imolty,ibox)
         return
      endif

c     --- aggregation volume bias move 3
c     --- B. Chen; J.I. Siepmann "Improving the Efficiency 
c     --- of the Aggregation-Volume-Bias Monte Carlo Algorithm"
c     --- J. Phys. Chem. B 105, 11275-11282 (2001).
      call twh_pmavb3(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         DEBUG_WRITE('MOVE: aggregation volume bias move 3')
         call twh_swapmoves(lfailure,SWAP_AVB_THREE,imolty,ibox)
         return
      endif

c     --- single molecule configurational bias regrowth move
      call twh_pmcb(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         DEBUG_WRITE('MOVE: one molecule conf bias regrowth')
         call twh_cbregrow(lfailure,CBR_ORIGINAL)
         return
      endif

c     --- single molecule configurational bias regrowth move
c     --- which performs a fixed-endpoint regrowth of interior
c     --- atoms along the protein backbone
      call twh_pmback(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         DEBUG_WRITE('MOVE: one molecule conf bias regrowth backbone')
         call twh_cbregrow(lfailure,CBR_BACKBONE)
         return
      endif

c     --- polypeptide sidechain regrowth move
      call twh_pmcbside(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         DEBUG_WRITE('MOVE: one molecule conf bias regrowth sidechain')
         call twh_cbregrow(lfailure,CBR_SIDECHAIN)
         return
      endif

c     --- pivot move 
      call twh_pmpivot(GLB_GET,pmtest)
      if ( rm .le. pmtest) then
         DEBUG_WRITE('MOVE: pivot move')
         call twh_pivot(lfailure)
         return
      endif

c     --- concerted rotation move 
      call twh_pmconrot(GLB_GET,pmtest)
      if ( rm .le. pmtest) then
         DEBUG_WRITE('MOVE: concerted rotation 1')
         call twh_conrot(lfailure,1)
         return
      endif

c     --- concerted rotation move 
      call twh_pmcrback(GLB_GET,pmtest)
      if ( rm .le. pmtest) then
         DEBUG_WRITE('MOVE: concerted rotation 0')
         call twh_conrot(lfailure,0)
         return
      endif

c     --- plane shift move
      call twh_pmplane(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         DEBUG_WRITE('MOVE: plane shift')
         call twh_atomshift(lfailure,1)
         return
      endif

c     --- row shift move
      call twh_pmrow(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         DEBUG_WRITE('MOVE: row shift')
         call twh_atomshift(lfailure,2)
         return
      endif

c     --- translation of atom positions
      call twh_pmtraat(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         DEBUG_WRITE('MOVE: translation atom positions')
         call twh_tranatom(lfailure)
         return
      endif

c     --- rotation around x,y, or z axis 
      call twh_pmcomposite(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         DEBUG_WRITE('MOVE: composite rotation and translation')
         call twh_composite(lfailure)
         return
      endif

c     --- translation of com in x,y, or z direction 
      call twh_pmtracm(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         DEBUG_WRITE('MOVE: translation of com in x,y,z')
         call twh_trancom(lfailure)
         return
      endif

c     --- rotation around x,y, or z axis 
      call twh_pmrotate(GLB_GET,pmtest)
      if ( rm .le. pmtest ) then
         DEBUG_WRITE('MOVE: rotation about x,y,z')
         call twh_rotate(lfailure)
         return
      endif

c     --- this should not happen as we should have found a move by now
      write(6,*) 'MAINLOOP: no move selected ',rm
      lfailure = .true.
      return
      end


      subroutine print_composite_move_summary()
c     ******************************************************************
c     * Output the composite move summary                              *
c     *                                                                *
c     * split from mainloop.F on 03-25-2006 by MAW                     *
c     * last modified 06-23-2009 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- local variables
      logical lnone
      integer ibox,imolty,nmolty,numboxes
      double precision datemp,dntemp

c     --- retrieve constants
      call twh_numboxes(GLB_GET,numboxes)

c     --- add the temporary composite move information into the
c     --- total just in case we did not end evenly on an update
      call twh_nmolty(GLB_GET,nmolty)
c     --- determine if any composite moves were attempted
      lnone = .true.
      do imolty = 1,nmolty
         do ibox = 1,numboxes
            call twh_acncomp(GLB_GET,imolty,ibox,dntemp)
            if ( dntemp .gt. 0.5d0 ) lnone = .false.
         enddo
      enddo
c     --- output performance of composite moves
      if ( .not. lnone ) then
         write(6,*)
         write(6,*) '* Composite Moves *'
         do ibox = 1,numboxes
            do imolty = 1,nmolty
               call twh_acncomp(GLB_GET,imolty,ibox,dntemp)
               if ( dntemp .gt. 0.5d0 ) then
                  call twh_acscomp(GLB_GET,imolty,ibox,datemp)
                  write(6,'(a9,i3,1x,a4,i2,1x,a9,1x,f15.0
     &                 ,1x,a9,1x,f15.0,1x,a9,1x,f7.3,1x,a1)') 
     &                 'Molecule:',imolty,'Box:',ibox
     &                 ,'Attempts:',dntemp
     &                 ,'Accepted:',datemp
     &                 ,'Accepted:',100.0d0*datemp/dntemp,'%'
               endif
            enddo
         enddo

      endif
      end

      subroutine twh_check_energies(lfailure)
c     ******************************************************************
c     * Check the running total of energies against engtotal           *
c     * calculation, and report if there are any problems.             *
c     *                                                                *
c     * Variables modified:                                            *
c     * Split from mainloop 04-26-2006 by MAW                          *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_FOREIGN_LAMBDA
#define FUNCTION_GET_DERIVATIVE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
c     --- logical scalars
      logical lproblem,lforeign_prob,lovrlap
c     --- integer scalars
      integer ibox,l,scp_num_foreign_lambda,numboxes
c     --- double precision scalars
      double precision v,vinter,vtail,vintra,vvib,vbend,vtg
      double precision vext,velect,vsolvation,vdum,vquantum,vtest
      double precision dvinter,dvelect
c     --- double precision arrays
      double precision vforeign
      dimension vforeign(MAX_FOREIGN_LAMBDA)
      double precision vf
      dimension vf(MAX_FOREIGN_LAMBDA)

c     --- retrieve constants
      call twh_numboxes(GLB_GET,numboxes)

c     --- check if the final value of the potential energy is consistent
c     --- with the running total of the energy changes
      write(6,*)
c     --- initialize lproblem
      lproblem = .false.
      do ibox = 1, numboxes
         write(6,'(a22,1x,i5)') 'Final Energies for Box',ibox
         call twh_engtotal(lfailure,lovrlap,v,vinter,vdum,vtail,vintra
     &        ,vvib,vbend,vtg,vext,velect,vsolvation,vquantum,vforeign
     &        ,dvinter,dvelect,ibox, .false.,.false.,.true.,.true. )
         if ( lfailure ) return
c        ---  need to check each of the types of energy 
         call twh_boxvclassic(GLB_GET,ibox,vtest)
         if ( (v - vtest)**2 .gt. ENERGY_THRESH ) then
            lproblem = .true.
            write(6,*) '### problem with Total Energy ###  box ',ibox
            write(6,*) ' current, running, difference: ',v,vtest
     &           ,v-vtest
         endif
         call twh_boxvinter(GLB_GET,ibox,vtest)
         if ( (vinter - vtest)**2 .gt. ENERGY_THRESH) then
            lproblem = .true.
            write(6,*) '### problem with Intermolecular Energy ###',
     &            '  box ',ibox
            write(6,*) ' current, running: ',vinter,vtest
         endif
         call twh_boxvtail(GLB_GET,ibox,vtest)
         if ( (vtail - vtest)**2 .gt. ENERGY_THRESH) then
            lproblem = .true.
            write(6,*) '### problem with Tail Correction Energy ###',
     &            '  box ',ibox
            write(6,*) ' current, running: ',vtail,vtest
         endif
         call twh_boxvintra(GLB_GET,ibox,vtest)
         if ( (vintra - vtest)**2 .gt. ENERGY_THRESH) then
            lproblem = .true.
            write(6,*) '### problem with Intramolecular Energy ###',
     &            '  box ',ibox
            write(6,*) ' current, running, difference: ',vintra,
     &            vtest,vintra-vtest
         endif
         call twh_boxvbond(GLB_GET,ibox,vtest)
         if ( (vvib - vtest)**2 .gt. ENERGY_THRESH) then
            lproblem = .true.
            write(6,*) '### problem with Bond Vibration Energy ###',
     &            '  box ',ibox
            write(6,*) ' current, running: ',vvib,vtest
         endif
         call twh_boxvangle(GLB_GET,ibox,vtest)
         if ( (vbend - vtest)**2 .gt. ENERGY_THRESH) then
            lproblem = .true.
            write(6,*) '### problem with Bond Bending Energy ###',
     &            '  box ',ibox
            write(6,*) ' current, running, difference: ',vbend,
     &            vtest,vbend-vtest
         endif
         call twh_boxvtorsion(GLB_GET,ibox,vtest)
         if ( (vtg - vtest)**2 .gt. ENERGY_THRESH) then
            lproblem = .true.
            write(6,*) '### problem with Torsion Energy ###',
     &            '  box ',ibox
            write(6,*) ' current, running, difference: ',vtg,vtest
     &           ,vtg-vtest
         endif
         call twh_boxvfield(GLB_GET,ibox,vtest)
         if ( (vext - vtest)**2 .gt. ENERGY_THRESH) then
            lproblem = .true.
            write(6,*) '### problem with External Energy ###',
     &            '  box ',ibox
            write(6,*) ' current, running: ',vext,vtest
         endif
         call twh_boxvelect(GLB_GET,ibox,vtest)
         if ( (velect - vtest)**2 .gt. ENERGY_THRESH) then
            lproblem = .true.
            write(6,*) '### problem with Coulomb Energy ###',
     &            '  box ',ibox
            write(6,*) ' current, running: ',velect,vtest
         endif
         call twh_boxvsolvation(GLB_GET,ibox,vtest)
         if ( (vsolvation - vtest)**2 .gt. ENERGY_THRESH) then
            lproblem = .true.
            write(6,*) '### problem with Solvation Energy ###',
     &            '  box ',ibox
            write(6,*) ' current, running: ',vsolvation,vtest
         endif
         call twh_boxvquantum(GLB_GET,ibox,vtest)
         if ( (vquantum - vtest)**2 .gt. ENERGY_THRESH) then
            lproblem = .true.
            write(6,*) '### problem with Quantum Energy ###',
     &            '  box ',ibox
            write(6,*) ' current, running: ',vquantum,vtest
            write(6,*) 'difference:',vquantum-vtest
         endif
         if ( (dvinter - scp_get_derivative(ibox, LAMBDA_LJ))**2 
     &         .gt. ENERGY_THRESH) then
            lproblem = .true.
            write(6,*) '### problem with dU/d[lambda_lj] ###',
     &            '  box ',ibox
            write(6,*) ' current, running: ',dvinter,
     &            scp_get_derivative(ibox, LAMBDA_LJ)
         endif
         if ( (dvelect - scp_get_derivative(ibox, LAMBDA_C))**2 
     &         .gt. ENERGY_THRESH) then
            lproblem = .true.
            write(6,*) '### problem with dU/d[lambda_c] ###',
     &            '  box ',ibox
            write(6,*) ' current, running: ',dvelect,
     &            scp_get_derivative(ibox, LAMBDA_C)
         endif
         call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
         if (scp_num_foreign_lambda .ne. 0) then
            call scp_get_foreign_energies(ibox, vf)
            lforeign_prob = .false.
            do l=1,scp_num_foreign_lambda
               if ( (vf(l) - vforeign(l))**2 .gt.  ENERGY_THRESH) then
                  lforeign_prob = .true.
                  lproblem = .true.
               endif
            enddo
            if (lforeign_prob) then
               write(6,*) '### problem: mismatch in foreign energies ',
     &            '### box ', ibox
               write(6,*) 'lambda_lj   lambda_c  current     running'
               do l=1,scp_num_foreign_lambda
                  write(6,*) scp_get_foreign_lambda(LAMBDA_LJ, l),
     &               scp_get_foreign_lambda(LAMBDA_C, l), vforeign(l), 
     &               vf(l)
               enddo
            endif
         endif
      enddo

c     --- final warning if there was a problem
      if ( lproblem ) then
         write(6,*) '### WARNING: energies calculated at the end of ',
     &        ' the run do not match running totals,',
     &        ' as described above###'
      endif

      return
      end


      subroutine twh_main_get_IO_redirect_info(lfailure,outputlocation
     &     ,lredirect)
c     ******************************************************************
c     * Obtain directory name for remapping files, as described in     *
c     * mainloop.                                                      *
c     * This subroutine returns with io unit 6 opened, if STDOUT is    *
c     * redirected.                                                    *
c     * Variables:                                                     *
c     *   outputlocation -- used for parallel routines.  The rank of   *
c     *     this node, which becomes the index to the directory listing*
c     *     in towhee_parallel.  As passed to towheemainloop           *
c     *   lredirect -- returned value indicating whether io unit 6 must*
c     *           be closed (if false, output goes to stdout)          *
c     *   lfailure -- set to true if error                             *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "wrapperF.h"
#define FUNCTION_GET_STRING_LENGTH
#define FUNCTION_READ_DIR_STRING
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lredirect,lfailure
      integer outputlocation
c     --- local variables
      character*MAXDIRLENGTH stmp,stdout_fn,io_dir
      logical lreadwrapper
      integer io_dl,iskip

      io_dir = ''
      call twh_lreadwrapper(GLB_GET,lreadwrapper)
      if (outputlocation .eq. 0 .and. .not. lreadwrapper) then
         lredirect = .false.
c        --- set the io_directory
         call twh_io_directory(GLB_SET,io_dir)
         return
      endif

      DEBUG_WRITE('outputlocation: ')
      DEBUG_WRITE(outputlocation)

      if ( outputlocation .ne. 0 ) then
         lredirect = .true.
c            --- first, handle odd case of outputlocation = 0 and 
c            --- lreadwrapper = .true. These two modes are mutually
c            --- exclusive; print an error and die.
         if (lreadwrapper) then
            write(6,*) "outputlocation .ne. 0 and lreadwrapper true."
            write(6,*) "Weird mix of modes is not supported"
            lfailure = .true.
            return
         endif

c        --- scan through the towhee_directories file for the current
c        --- working directory
c        --- NOTE: this is an exception to the normal rules as this
c        --- is the sequence for finding the directory information
         open(21,file='towhee_parallel',form='formatted')
c        --- skip the top two lines (one of these is a comments)
         call twh_skip_lines(lfailure,21,1)
         if ( lfailure ) return
c        --- get the suffix for the standard out filename
         call twh_read_string(lfailure,21,'filename', .false.,stdout_fn)
         if ( lfailure ) return
c        --- get the local working directory
         do iskip = 1,outputlocation-1
            stmp = twh_read_dir_string(lfailure,21,'ignored dir'
     &           , .false.)
            if ( lfailure ) return
         enddo
         io_dir  = twh_read_dir_string(lfailure,21,'working dir'
     &        , .false.)
         if ( lfailure ) return

         rewind(21)
         close(21)
      elseif (lreadwrapper) then
         lredirect = wrap_lredirect_stdout
         io_dl = twh_get_string_length(wrap_output_dir)
         io_dir = wrap_output_dir(1:io_dl)

c     --- construct stdout_fn, where standard out is redirected to.
c     --- for the case of wrappers, we append "towhee.out" to the output
c     --- directory given by wrap_output_dir
         stdout_fn = 'towhee.out'
      endif

c     --- compute the length of the directory characters
      io_dl = twh_get_string_length(io_dir)
c        --- make sure the final character is a '/'
      if ( io_dir(io_dl:io_dl) .ne. '/' ) then
         io_dir(io_dl+1:io_dl+1) = '/'
         io_dl = io_dl + 1
      endif

      if ( io_dl  .gt. MAXDIRLENGTH ) then
         write(6,*) 'TOWHEE: dirlength exceed MAXDIRLENGTH'
         lfailure = .true.
         return
      endif

c     --- set the io_directory
      call twh_io_directory(GLB_SET,io_dir)

      if (lredirect) then
c        --- open standard output.  This is opened with "NEW", to 
c        --- prevent us from overwriting exising output.
         call twh_openfile(lfailure,6,stdout_fn,io_dir,'UNKNOWN'
     &        ,'DEFAULT')
         if ( lfailure ) return
      endif

      return
      end
