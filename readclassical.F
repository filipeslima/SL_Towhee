#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2011 Marcus G. Martin                       *
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
      subroutine twh_readclassical(lfailure,numtypes)
c     ******************************************************************
c     * reads in all input parameters from towhee_input that are       *
c     * related to the classical potential style                       *
c     *                                                                *
c     * split out of readtowhee 11-01-2004 by M.G. Martin              *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_CHECK_LABEL
#define FUNCTION_READ_DIR_STRING
#define FUNCTION_READ_FLOAT
#define FUNCTION_READ_LABELED_FLOAT
#define FUNCTION_READ_LABELED_INTEGER
#define FUNCTION_READ_LABELED_LOGICAL
#include "functions.h"

c     --- variables passed to/from the subroutine
      logical lfailure
      integer numtypes
c     --- local variables
c     --- logical scalars
      logical lreadmixrule,lreadlshift,lreadltailc,lreadrmin
      logical lreadrcut,lreadrcutin,lreadinterpolatestyle,lreadrpd
      logical lreadscaling,scp_scale_coulomb,loutlammps,lshift
      logical ltailc,lewald
c     --- character strings
      character*(FFNAMELEN) cvalue,keynames
      dimension keynames(2)
      character*20 idname,coulombstyle,interpolatestyle
      character*30 cplocal,cmlocal
      character*50 electrostatic_form
      character*(MAXDIRLENGTH) ffname
      integer nbondtypes,nangletypes,ndihedtypes,nimprotypes,numboxes
      integer noftypes,ffnumber
      integer ifile,ibox,naatypes,nbitypes
      integer nhbondtypes,classical_pot_num
      integer icoeff,iadtotal,iadjust,ikey,kmax
c     --- double precision scalars
      double precision dielect,dtemp,ewald_prec,kalp,dvalue
      double precision rcelect,rcut,rcutsq,rminsq,rmin,rcutin

c     --- retrieve constants
      call twh_numboxes(GLB_GET,numboxes)

c     --- At some point, it may be good   
c     --- to strictly separate the reading of user parameters from the
c     --- initialization and calculation of internal parameters.  
                                                                     
c     --- set subroutine idname
      idname = 'readclassical'

c     --- read in number of force fields
      ffnumber = twh_read_labeled_integer(lfailure,4, 'ffnumber'
     &     , .true., idname)
      if ( lfailure ) return
      call twh_ffnumber(GLB_SET,ffnumber)

      if ( ffnumber .lt. 1 ) then
         write(6,*) 
     &        'READCLASSICAL: must specify at least one forcefield file'
         lfailure = .true.
         return
      endif

c     --- read in ff_filename
      call twh_read_label(lfailure,4, 'ff_filename',.true.,idname)
      if ( lfailure ) return
      do ifile = 1,ffnumber
         ffname = twh_read_dir_string(lfailure,4,'ff_filename',.true.)
         if ( lfailure ) return
         call twh_ff_filename(GLB_SET,ifile,ffname)
      enddo

c     --- read potential name
      call twh_read_labeled_string(lfailure,4,'classical_potential'
     &     ,.true. ,idname,cplocal)
      if ( lfailure ) return
      call twh_classical_potential(GLB_SET,cplocal)
c     --- get the classical potential integer
      call twh_cp_number(GLB_GET,classical_pot_num)
      write(6,'(1x,a,1x,i5,1x,a,1x,a)') 'READCLASSICAL: pot_num:'
     &     ,classical_pot_num,'potential name:',cplocal

c     --- determine the read logicals
      call twh_setpotentiallogic(lfailure,lreadmixrule,lreadlshift
     &     ,lreadltailc,lreadrmin,lreadrcut,lreadrcutin
     &     ,lreadinterpolatestyle,lreadrpd,lreadscaling)
      if ( lfailure ) return
      
      if ( lreadmixrule ) then
c        --- read in the mixing rule
         call twh_read_labeled_string(lfailure,4,'classical_mixrule'
     &        , .true. , idname,cmlocal)
         if ( lfailure ) return
         call twh_classical_mixrule(GLB_SET,cmlocal)

         if ( cmlocal .eq. 'LB plus manual' ) then
c           --- get all of the manual scaling information
            iadtotal = twh_read_labeled_integer(lfailure,4
     &           ,'mixrule_adjust_total',.true.,idname)
            if ( lfailure ) return
            call twh_mixrule_adjust_total(GLB_SET,iadtotal)

            do iadjust = 1,iadtotal
               call twh_read_labeled_string(lfailure,4
     &              ,'mixrule_adjust_key',.true.,idname,cvalue)
               if ( lfailure ) return
               call twh_mixrule_adjust_key(GLB_SET,iadjust,cvalue)

               call twh_read_labeled_string_array(lfailure,4
     &              ,'mixrule_adjust_keynames',2,.true.,idname,keynames)
               if ( lfailure ) return
               do ikey = 1,2
                  call twh_mixrule_adjust_keynames(GLB_SET,iadjust,ikey
     &                 ,keynames(ikey))
               enddo

               call twh_read_label(lfailure,4,'mixrule_adjustments'
     &              ,.true. ,idname)
               if ( lfailure ) return
               do icoeff = 1,4
                  dtemp = twh_read_float(lfailure,4
     &                 ,'mixrule_adjustments',.true.)
                  if ( lfailure ) return
                  call twh_mixrule_adjustments(GLB_SET,iadjust,icoeff
     &                 ,dtemp)
               enddo
            enddo
         endif
      endif

c     Read potential scaling parameters for select classical potentials.
      if ( lreadscaling ) then
         call twh_read_scaling(lfailure)
         if ( lfailure ) return
      endif

      if ( lreadinterpolatestyle ) then
c        --- read in the interpolatestyle
         call twh_read_labeled_string(lfailure,4,'interpolatestyle'
     &        , .true.  ,idname,interpolatestyle)
         if ( lfailure ) return
         call twh_interpolatestyle(GLB_SET,interpolatestyle)
c        --- error check
         if ( interpolatestyle .eq. 'cubicspline' ) then
            write(6,*) 'Tables interpolated using cubic splines'
         elseif ( interpolatestyle .eq. 'linear' ) then
            write(6,*) 'Tables interpolated linearly'
         else
            write(6,*) 'unknown interpolatestyle'
            lfailure = .true.
            return
         endif
      endif
      
      if ( lreadlshift ) then
c        --- read in cut-and-shift logical
         lshift = twh_read_labeled_logical(lfailure,4, 'lshift'
     &        , .true., idname)
         if ( lfailure ) return
         call twh_lshift(GLB_SET,lshift)
      endif
      
      if ( lreadltailc ) then
c        --- read in analytical tail correction logical
         ltailc = twh_read_labeled_logical(lfailure,4, 'ltailc'
     &        , .true., idname)
         if ( lfailure ) return
         call twh_ltailc(GLB_SET,ltailc)
      endif
      
      if ( lreadrmin ) then
c        --- read in distance below which we assume hard overlap
         rmin = twh_read_labeled_float(lfailure,4, 'rmin', .true.
     &        , idname)
         if ( lfailure ) return
         call twh_rmin(GLB_SET,rmin)
c        --- error check
         if ( rmin .lt. 0.0d0 ) then
            write(6,*) 'READCLASSICAL: rmin must be non-negative'
            lfailure = .true.
            return
         endif
      endif
c     --- retrieve rmin for use in the rest of this subroutine
      call twh_rmin(GLB_GET,rmin)
      
      if ( lreadrcut ) then
c        --- read in vdw cutoff
         rcut = twh_read_labeled_float(lfailure,4, 'rcut', .true.
     &        , idname)
         if ( lfailure ) return
         call twh_rcut(GLB_SET,rcut)
c        --- error check
         if ( rcut .lt. 0.0d0 ) then
            write(6,*) 'READCLASSICAL: rcut must be non-negative'
            lfailure = .true.
            return
         endif
         if ( rcut .lt. rmin ) then
            write(6,*) 'READCLASSICAL: rcut must be greater than rmin'
            lfailure = .true.
            return
         endif
c        --- set rcutsq
         rcutsq = rcut*rcut
         call twh_rcutsq(GLB_SET,rcutsq)
      endif
      
      if ( lreadrcutin ) then
c        --- read in inner cutoff for configuration-bias move
         rcutin = twh_read_labeled_float(lfailure,4, 'rcutin', .true.
     &        , idname)
         if ( lfailure ) return
         call twh_rcutin(GLB_SET,rcutin)
c        --- error check
         if ( rcutin .lt. 0.0d0 ) then
            write(6,*) 'READCLASSICAL: rcutin must be non-negative'
            lfailure = .true.
            return
         endif
         if ( rmin .gt. rcutin ) then
            write(6,*) 'READCLASSICAL: rcutin must be greater than rmin'
            lfailure = .true.
            return
         endif
      endif
      if ( lreadrpd ) then
c        --- read in the radial_pressure_delta
         dvalue = twh_read_labeled_float(lfailure,4
     &        ,'radial_pressure_delta', .true., idname)
         if ( lfailure ) return
         call twh_radial_pressure_delta(GLB_SET,dvalue)
         if ( dvalue .le. 0.0d0 ) then
            write(6,*) 'radial_pressure_delta must be strictly positive'
            lfailure = .true.
            return
         endif
      endif
      
c     --- set rminsq
      rminsq = rmin*rmin
      call twh_rminsq(GLB_SET,rminsq)
      
      call twh_read_labeled_string(lfailure,4,'electrostatic_form'
     &     ,.true. ,idname,electrostatic_form)
      if ( lfailure ) return
      call twh_electrostatic_form(GLB_SET,electrostatic_form)

      if ( electrostatic_form .eq. 'none' ) then
c        --- no electrostatic interactions
         call twh_lcoulomb(GLB_SET,.false.)
         call twh_lewald(GLB_SET,.false.)
         call twh_ewald_kalp(GLB_SET,0.0d0)
         call twh_ewald_kmax(GLB_SET,1,0)
c        --- reset qqfact to zero
         dvalue = 0.0d0
         call twh_constant_qqfact(GLB_SET,dvalue)
      elseif ( electrostatic_form .eq. 'coulomb' ) then
         call twh_lcoulomb(GLB_SET,.true.)
c        --- read in the coulombstyle
         call twh_read_labeled_string(lfailure,4, 'coulombstyle'
     &        ,.true.,idname,coulombstyle)
         if ( lfailure ) return
         call twh_coulombstyle(GLB_SET,coulombstyle)
         if ( coulombstyle .eq. 'ewald_fixed_kmax' ) then
c           --- coulombic interactions with fixed kmax values
            call twh_lewald(GLB_SET,.true.)
c           --- read in ewald sum kalp and kmax parameters
            kalp = twh_read_labeled_float(lfailure,4, 'kalp', .true.
     &           , idname)
            if ( lfailure ) return
            call twh_ewald_kalp(GLB_SET,kalp)
            if ( kalp .le. 0.0d0 ) then
               write(6,*) 'READCLASSICAL:kalp must be greater than 0.0'
               lfailure = .true.
               return
            endif
            kmax = twh_read_labeled_integer(lfailure,4,'kmax',.true.
     &           ,idname)
            if ( lfailure ) return
            if ( kmax .le. 0 ) then
               write(6,*) 'READCLASSICAL: kmax must be greater than 0'
               lfailure = .true.
               return
            elseif ( kmax .gt. MAXKMAX ) then
               write(6,*) 'READCLASSICAL: kmax exceeds MAXKMAX'
               write(6,*) 'either lower kmax or modify MAXKMAX'
     &              , 'in preproc.h'
               lfailure = .true.
               return
            endif
c           --- set all of the boxes to this initial kmax
            do ibox = 1,numboxes
               call twh_ewald_kmax(GLB_SET,ibox,kmax)
            enddo
         elseif ( coulombstyle .eq. 'ewald_fixed_cutoff' ) then
c           --- coulombic interactions with fixed electrostatic cutoff
            call twh_lewald(GLB_SET,.true.)
c           --- read in the ewald precision
            ewald_prec = twh_read_labeled_float(lfailure,4
     &           , 'ewald_prec', .true., idname)
            if ( lfailure ) return
            call twh_ewald_prec(GLB_SET,ewald_prec)
c           --- safety check
            if ( ewald_prec .le. 0.0d0 ) then
               write(6,*) 'READCLASSICAL: ewald_prec must be positive'
               lfailure = .true.
               return
            endif
c           --- read in the electostatic cutoff
            rcelect = twh_read_labeled_float(lfailure,4,'rcelect'
     &           ,.true., idname)
            if ( lfailure ) return
            if ( rcelect .lt. 0.0d0 ) then
               write(6,*) 'READCLASSICAL: rcelect must be non-negative'
               lfailure = .true.
               return
            endif
c           --- set all of the boxes to this initial rcelect
            call twh_rcelect(GLB_INIT,numboxes,rcelect)
         elseif ( coulombstyle .eq. 'minimum image' ) then
c            --- explicit minimum image coulombic interactions
            call twh_lewald(GLB_SET,.false.)
            kmax = 0
            call twh_ewald_kalp(GLB_SET,0.0d0)
c           --- set all of the boxes to this initial kmax
            do ibox = 1,numboxes
               call twh_ewald_kmax(GLB_SET,ibox,kmax)
            enddo
         else
c           --- unknown coulombstyle
            write(6,*) 'READCLASSICAL: unknown coulombstyle'
            lfailure = .true.
            return
         endif

c        --- read in the dielectric constant
         dielect = twh_read_labeled_float(lfailure,4, 'dielect', 
     &        .true., idname)
         if ( lfailure ) return
c        --- make sure dielectric constant is greater than zero
         if ( dielect .le. 0.0d0 ) then
            write(6,*) 'READCLASSICAL: dielect must be positive '
            lfailure = .true.
            return
         endif
c        --- reset qqfact to account for the dielectric constant
         dvalue = CON_EETOKCAL*CON_KCALTOK/dielect
         call twh_constant_qqfact(GLB_SET,dvalue)
      else
         write(6,*) 'READCLASSICAL: unknown electrostatic_style'
         lfailure = .true.
         return
      endif

c     --- incomplete implementation issue.  How do self and correct
c     --- terms
c     --- scale with cmix_lambda_c?  At this point they don't, and so 
c     --- will yield incorrect corrections.  This should be fixed
c     --- at some point.
      call twh_scp_scale_coulomb(GLB_GET,scp_scale_coulomb)
      call twh_lewald(GLB_GET,lewald)
      if (lewald .and. scp_scale_coulomb) then
         write(6,*) "WARNING: Coulomb scaling is not fully implemented",
     &      "for Ewald."
         write(6,*) "Self and correction energies may be incorrect. "
      endif
      

c ******************************************
      
c     --- read the forcefield from files
      write(6,*) 'Setting up force field parameters from files'
c     --- initialize the parameters
      numtypes = 0
      nbondtypes = 0
      nangletypes = 0
      ndihedtypes = 0
      nimprotypes = 0
      naatypes = 0
      nbitypes = 0
      noftypes = 0
      nhbondtypes = 0
      
c     --- read in the force field parameters from the files
      do ifile = 1,ffnumber
         write(6,'(a25,i5)') ' opening forcefield file:',ifile
        call twh_rwforcefield(lfailure,READ_FLAG,numtypes,nbondtypes
     &        ,nangletypes,ndihedtypes,nimprotypes,naatypes,nbitypes
     &        ,noftypes,nhbondtypes,ifile)
        if ( lfailure ) return
      enddo
      
c     --- set up the nonbonded mixing and convert bonded
      call twh_setmixterms(lfailure,numtypes)
      if ( lfailure ) return

c     --- optional read max_bond_length(default value set in initialize)
      if ( twh_check_label(lfailure,4,'max_bond_length') ) then
         dvalue = twh_read_labeled_float(lfailure,4, 'max_bond_length'
     &        ,.true., idname)
         if ( lfailure ) return
         call twh_max_bond_length(GLB_SET,dvalue)
c        --- error check
         if ( dvalue .lt. 0.0d0 ) then
            write(6,*) 'READCLASSICAL: max_bond_length must be positive'
            write(6,*) dvalue
         endif
      else
c        --- retrieve the default and output it to std_output
         call twh_max_bond_length(GLB_GET,dvalue)
         write(6,'(a24,1x,f5.2)') 'default max_bond_length:',dvalue
      endif
      
c     --- optional read in the number of external fields
      if ( twh_check_label(lfailure,4,'nfield') ) then
         call twh_read_field(lfailure,4)
      endif

      call twh_loutlammps(GLB_GET,loutlammps)
      if ( loutlammps ) then
         call twh_writedynamo(lfailure,numtypes)
         if ( lfailure ) return
      endif
      
      return
      end

      subroutine twh_read_scaling(lfailure)
c     ******************************************************************
c     * Reads potential scaling information, and appropriately sets    *
c     * variables in scaled_potential.F.                               *
c     *                                                                *
c     * Variables:                                                     *
c     *                                                                *
c     * split from twh_readclassical on 04-27-2006 by MAW              *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
c     --- Linear scaling of potentials (other than Scaled Lennard Jones)
c     ---- with lambda is possible, and is mostly implemented internally
c     --- The input format here currently does not support it, but if
c     ---- there is interest, this shouldn't be hard.  
c     --- See vcoulomb and vtwobody for implementation, and 
c     --- scaled_potential general architecture.
c     --- In addition, there is no reason why coulomb and twobody
c     --- scaling need to be coupled; they are logically distinct.
      implicit none
#define FUNCTION_GET_SCALING_STYLE
#define FUNCTION_READ_LABELED_FLOAT
#define FUNCTION_READ_LABELED_LOGICAL
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
      logical scp_eval_derivative
      double precision cmix_lambda_lj,cmix_lambda_c
      double precision cmix_foreign_lambda(MAX_FOREIGN_LAMBDA,2)
      double precision v(MAX_FOREIGN_LAMBDA)
      integer cmix_npair, cmix_num_foreign_lambda, l
      character*10 cmix_pair_list 
      dimension cmix_pair_list(CROSSTYPEMAX,2)
      character*20 idname
      character*30 cvalue

      idname = 'twh_read_scaling'
            
c     --- read in the cmix rescaling style
      call twh_read_labeled_string(lfailure,4,'cmix_rescaling_style'
     &     , .true.,idname,cvalue)
      if ( lfailure ) return
      call scp_set_scaling_style(cvalue)

c     For rescaling_style 'none', do nothing; parameters as set will
c     reproduce 6-12 LJ
      if ( scp_get_scaling_style() .eq. SCP_SC ) then
c     --- soft-core scaling
c           --- cross terms in the style of Beutler et al. 1994
c           --- the following parameters are mandatory
c           ---   cmix_lambda_lj (0..1)
c           ---   cmix_alpha_lj (float) 
c           ---   cmix_lambda_c (0..1)
c           ---   cmix_lprintdudl (boolean)
c           ---   cmix_npair (integer)
c           ---   cmix_pair_list (FF+atom pair, npair entries)

c           --- get cmix_alpha_lj and cmix_alpha_c.  Special function
c           --- called, because these parameters may be passed from
c           --- the C wrapper (towhee.c)
         call twh_get_lambda(lfailure,4, cmix_lambda_lj, cmix_lambda_c, 
     &      .true., idname)
         if ( lfailure ) return
         call scp_set_native_lambda(LAMBDA_LJ, cmix_lambda_lj, 
     &         'cmix_lambda_lj')
         call scp_set_native_lambda(LAMBDA_C, cmix_lambda_c, 
     &         'cmix_lambda_c')

c        read and set alpha_lj in one go.
         call scp_set_native_lambda(ALPHA_LJ, 
     &       twh_read_labeled_float(lfailure,4, 'cmix_alpha_lj'
     &        , .true., idname),'cmix_alpha_lj')
         if ( lfailure ) return

c        Print out dU/d[lambda_lj] and dU/d[lambda_coulomb]?  
         scp_eval_derivative = twh_read_labeled_logical(lfailure,4, 
     &      'cmix_lprintdudl', .true., idname)
         if ( lfailure ) return
         call twh_scp_eval_derivative(GLB_SET,scp_eval_derivative)

         call twh_get_scaled_pairs(lfailure,cmix_npair, cmix_pair_list)
         if ( lfailure ) return
         call scp_set_scaled_atoms(cmix_npair, cmix_pair_list)

         
c       --- get the foreign lambdas.  These may come from towhee_input
c       --- (where they are optional) or from the C wrapper (towhee.c)
         call twh_get_foreign_lambda(lfailure,4
     &        , cmix_num_foreign_lambda,cmix_foreign_lambda, .true.
     &        , idname)
         if ( lfailure ) return
         call scp_set_num_foreign_lambda(cmix_num_foreign_lambda)

         do l=1,cmix_num_foreign_lambda
            v(l) = cmix_foreign_lambda(l,LAMBDA_LJ)
         enddo
         call scp_set_foreign_lambdas(LAMBDA_LJ, v)

         do l=1,cmix_num_foreign_lambda
            v(l) = cmix_foreign_lambda(l,LAMBDA_C)
         enddo
         call scp_set_foreign_lambdas(LAMBDA_C, v)

         call twh_scp_scale_coulomb(GLB_SET,.true.)
         call twh_scp_scale_twobody(GLB_SET,.true.)


c     ---  Weeks-Chandler-Anderson also uses cmix_list
c     ---  to keep track of who is repulsive-only
      elseif ( scp_get_scaling_style() .eq. SCP_WCA) then
         call twh_get_scaled_pairs(lfailure,cmix_npair, cmix_pair_list)
         if ( lfailure ) return
         call scp_set_scaled_atoms(cmix_npair, cmix_pair_list)
         
         call twh_scp_scale_coulomb(GLB_SET,.false.)
         call twh_scp_scale_twobody(GLB_SET,.true.)
         
      else
c        --- shouldn't get here, as scp_set_scaling_style() should catch
c        --- this.
         write(6,*) 'unknown cmix_rescaling_style'
         lfailure = .true.
         return
      endif

      return
      end


      subroutine twh_get_scaled_pairs(lfailure,cmix_npair
     &     ,cmix_pair_list)
c     ******************************************************************
c     * A convenience function to read in cmix_npair and cmix_pair_list*
c     *                                                                *
c     * split from twh_readclassical on 04-27-2006 by MAW              *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer cmix_npair
      character*10 cmix_pair_list 
      dimension cmix_pair_list(CROSSTYPEMAX,2)
c     --- local variables
      character*25 idname
      character*10 satmp(2)
      integer ipair, ientry

      idname = 'twh_get_scaled_pairs'


c           --- get number of single atoms to rescale
c           --- this is stored in cmix_npair as that is the more
c           --- general version of the array
      cmix_npair = twh_read_labeled_integer(lfailure,4, 
     &   'cmix_npair', .true., idname)
      if ( lfailure ) return
      if ( cmix_npair .le. 0 ) then
         write(6,*) 'cmix_npair must be positive'
         lfailure = .true.
         return
      endif
      if ( cmix_npair .gt. CROSSTYPEMAX ) then
         write(6,*) 'cmix_npair must be less than or equal to'
     &        ,' crosstypemax (set in preproc.h) '
     &        ,CROSSTYPEMAX
         lfailure = .true.
         return
      endif

c           --- get the cmix pair list 
      call twh_read_label(lfailure,4,'cmix_pair_list',.false.,idname)
      if ( lfailure ) return
      do ipair = 1,cmix_npair
         call twh_read_string_array(lfailure,4,'cmix_pair_list',2
     &        ,.false.,satmp)
         if ( lfailure ) return
         cmix_pair_list(ipair,1) = satmp(1)
         cmix_pair_list(ipair,2) = satmp(2)
         write(6,'(a5,1x,i5,1x,a15,1x,4(a10,1x))') 
     &        'pair:',ipair,'cmix_pair_list:'
     &        ,(cmix_pair_list(ipair,ientry),ientry=1,2)
      enddo

      return
      end


   
      subroutine twh_get_lambda(lfailure,io_unit,lambda_lj,lambda_c
     &     ,print_value,idname)
c     ******************************************************************
c     * Returns lambda_lj and lambda_c parameters.  These values       *
c     * may come from one of two sources: from the towhee_input file   *
c     * or, if lreadwrapper is set (in a call to towheemainloop) from  *
c     * variables set (by the C wrapper) in the wrapper common block.  *
c     *                                                                *
c     * Variables:                                                     *
c     * io_unit specifies the file number to read from both lambda_lj  *
c     * and lambda_c are set by this subroutine.  Note that these are  *
c     * local vars, *not* the global vars cmix_lambda_                 *
c     * print_value specifies whether the label and value output to    *
c     * screen.                                                        *
c     * idname is the calling subroutine, printed only in case of error*
c     *                                                                *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "wrapperF.h"
#define FUNCTION_CHECK_LABEL
#define FUNCTION_READ_LABELED_FLOAT
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer io_unit
      double precision lambda_lj, lambda_c
      logical lfailure,print_value
      character*(*) idname
c     --- local variables
      logical lreadwrapper
      double precision dptmp

      call twh_lreadwrapper(GLB_GET,lreadwrapper)
      if (lreadwrapper) then
c        --- if lambdas come from wrapper, then these parameters are 
c        --- allowed, but optional, in towhee_input.  See if they exist
c        --- and read them so that they're skipped over.
         if (twh_check_label(lfailure,io_unit,'cmix_lambda_lj')) then
            dptmp = twh_read_labeled_float(lfailure,io_unit
     &           , 'cmix_lambda_lj', .false., idname)
            if ( lfailure ) return
         endif
         if (twh_check_label(lfailure,io_unit,'cmix_lambda_c')) then
            dptmp = twh_read_labeled_float(lfailure,io_unit
     &           , 'cmix_lambda_c', .false., idname)
            if ( lfailure ) return
         endif

c        --- actual values come from wrapper
         lambda_lj = wrap_lambda_lj
         lambda_c = wrap_lambda_c
         if (print_value) then
            write(6,*) 'cmix_lambda_lj (from wrapper): ',lambda_lj
            write(6,*) 'cmix_lambda_c (from wrapper): ',lambda_c
         endif
      else
         lambda_lj = twh_read_labeled_float(lfailure,io_unit
     &        , 'cmix_lambda_lj', print_value, idname)
         if ( lfailure ) return
         lambda_c = twh_read_labeled_float(lfailure,io_unit
     &        , 'cmix_lambda_c', print_value, idname)
         if ( lfailure ) return
      endif
      return
      end

      subroutine twh_get_foreign_lambda(lfailure,io_unit
     &     , num_foreign_lambda,foreign_lambda, print_value, idname)
c     ******************************************************************
c     * Returns the array foreign_lambda, as well as its size,         *
c     * num_foreign_lambda                                             *
c     * these values may come from one of two sources: from the        *
c     * towhee_input file or, if lreadwrapper is set (in a call to     *
c     * towheemainloop) from variables  set (by the C wrapper) in the  *
c     * wrapper common block.                                          *
c     * Note that in towhee_input, the parameter num_foreign_lambda is *
c     * an optional parameter.                                         *
c     * Note also that the parameters foreign_lambda and               *
c     * num_foreign_lambda returned by this subroutine are not the     *
c     * global ones (cmix_*); it is expected that the calling function *
c     * set the global parameters.                                     *
c     *                                                                *
c     * Variables:                                                     *
c     * io_unit specifies the file number to read from                 *
c     * num_foreign_lambda is set by this subroutine                   *
c     * foreign_lambda is set by this subroutine                       *
c     * print_value specifies whether the label and value output to    *
c     * screen.                                                        *
c     * idname is the calling subroutine, printed only in case of error*
c     * lfailure is a flag which is set to true in this subroutine in  *
c     * case of error.                                                 *
c     *                                                                *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "wrapperF.h"
#define FUNCTION_CHECK_LABEL
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"

c     --- variables passed to/from the subroutine
      double precision foreign_lambda(MAX_FOREIGN_LAMBDA,2)
      integer io_unit,num_foreign_lambda
      logical print_value,lfailure
      character*(*) idname
c     --- local variables
      logical pv,lreadwrapper
      integer i
      double precision dvalue
      double precision fatmp
      dimension fatmp(2)

c    --- read the optional argument num_foreign_lambda, and read
c    --- in the lambda pairs if required.  Note that this is done 
c    --- whether or not the final values of foreign_lambda come
c    --- from the wrapper, since we need to advance the file beyond
c    --- this section in either case.
c    --- if we will end up reading from wrapper, however, make sure
c    --- file reading does not output to screen - done with pv

      call twh_lreadwrapper(GLB_GET,lreadwrapper)
      if (lreadwrapper) then
         pv = .false.
      else
         pv = print_value
      endif
      if (twh_check_label(lfailure,io_unit,'num_foreign_lambda')) then
         num_foreign_lambda = twh_read_labeled_integer(lfailure,io_unit
     &        , 'num_foreign_lambda', pv, idname)
         if ( lfailure ) return
c    --- check size of num_foreign_lambda
         if (num_foreign_lambda .gt. MAX_FOREIGN_LAMBDA) then
            write(6,*) 'num_foreign_lambda must be less than',
     &          'or equal to MAX_FOREIGN_LAMBDA (set in preproc.h)'
            lfailure = .true.
            return
         endif

c    --- read in the lambda pairs
         call twh_read_label(lfailure,io_unit,'foreign_lambda',pv
     &        , idname)
         if ( lfailure ) return
         do i=1,num_foreign_lambda
            call twh_read_float_array(lfailure,io_unit
     &           ,'foreign_lambda',2,pv,fatmp)
            if ( lfailure ) return
            foreign_lambda(i,LAMBDA_LJ) = fatmp(LAMBDA_LJ)
            foreign_lambda(i,LAMBDA_C) = fatmp(LAMBDA_C)
         enddo
      endif

      if (lreadwrapper) then
         num_foreign_lambda = wrap_num_foreign_lambda

         if (num_foreign_lambda .gt. MAX_FOREIGN_LAMBDA .or.
     &       num_foreign_lambda .lt. 0) then
c  --- this could happen if values in C wrapper not properly assigned...
            write(6,*) "number foreign lambda, from wrapper, set to ",
     &         num_foreign_lambda
            write(6,*) "must be between 0 and MAX_FOREIGN_LAMBDA"
            lfailure = .true.
            return
         endif
         if (print_value) then
            write(6,*) "Number foreign lambda (from wrapper): ",
     &         num_foreign_lambda
            write(6,*) "Foreign lambda (from wrapper): ",
     &         "(lambda_lj, lambda_c) "
         endif
         do 220,i=1,num_foreign_lambda
            call twh_wrap_foreign_lambda_lj(GLB_GET,i,dvalue)
            foreign_lambda(i,LAMBDA_LJ) = dvalue
            call twh_wrap_foreign_lambda_c(GLB_GET,i,dvalue)
            foreign_lambda(i,LAMBDA_C) = dvalue
            if (print_value) then
               write(6,*) "wrapper: ", foreign_lambda(i,LAMBDA_LJ), 
     &            foreign_lambda(i,LAMBDA_C)
            endif
220      continue
      endif
      return
      end

      subroutine twh_read_field(lfailure,ifile)
c     ******************************************************************
c     * reads all of the information related to the additional fields  *
c     *                                                                *
c     * originally split from readclassical.F 10-26-2008 by M.G. Martin*
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- functions
#define FUNCTION_GET_CLASSICAL_POTENTIAL
#define FUNCTION_GET_ELEMENTNAME
#define FUNCTION_GET_NBNAME
#define FUNCTION_READ_LABELED_FLOAT
#define FUNCTION_READ_LABELED_INTEGER
#define FUNCTION_READ_LABELED_LOGICAL
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ifile
c     --- local variables
c     --- character string scalars
      character*11 fchar
      character*15 idname
      character*20 fieldstring
      character*(FFNAMELEN) atomname
      character*30 cplocal
c     --- logical scalars
      logical ltemp
c     --- integer scalars
      integer nfield,hafmolec,ientry,ifield,itype,ivalue,ntii,type
      integer numboxes
c     --- double precision scalars
      double precision dfield,sigma_sf,epsilon_sf,fthree,ffour
      double precision fnbparamone,fnbparamtwo,fnbparam
      double precision invdist,invdist3,invdist9,ljfsig,ljfeps
      double precision zzz,szone,sztwo,szfour,szten,twopi
c     --- double precision arrays
      double precision fatmp
      dimension fatmp(3)

c     --- retrieve constants
      call twh_numboxes(GLB_GET,numboxes)
      call twh_constant_twopi(GLB_GET,twopi)

c     --- set values
      idname = 'read_field'

      nfield = twh_read_labeled_integer(lfailure,ifile, 'nfield'
     &     , .true., idname)
      if ( lfailure ) return
      call twh_nfield(GLB_SET,nfield)
c     --- make sure nfield is non-negative
      if ( nfield .lt. 0 ) then
         write(6,*) 'READ_FIELD: nfield must be non-negative'
         lfailure = .true.
         return
      endif
      if ( nfield .gt. FLDMAX ) then
         write(6,*) 'READ_FIELD: nfield exceeds FLDMAX'
         write(6,*) 'increase the value of FLDMAX in preproc.h'
         lfailure = .true.
         return
      endif
         
      do ifield = 1,nfield
c        --- read in field type
         call twh_read_labeled_string(lfailure,ifile, 'fieldtype'
     &        ,.false.,idname,fieldstring)
         if ( lfailure ) return
         call twh_fieldtypestring(GLB_SET,ifield,fieldstring)
         write(6,'(a6,1x,i5,1x,a10,1x,a20)') 'Field:',ifield
     &        ,'fieldtype:',fieldstring
         
         if ( fieldstring .eq. 'Hard Wall' ) then
c           --- read hard wall information
c           --- box for hard wall field
            ivalue = twh_read_labeled_integer(lfailure,ifile,'hrdbox'
     &           ,.true.,idname)
            if ( lfailure ) return
            call twh_fieldint(GLB_SET,ifield,1,ivalue)
            if ( ivalue .lt. 1 ) then
               write(6,*) 'hrdbox must be positive'
               lfailure = .true.
               return
            endif
            if ( ivalue .gt. numboxes ) then
               write(6,*) 'hrdbox must be <= numboxes'
               lfailure = .true.
               return
            endif
            
c           --- x,y, or z dimension
            call twh_read_labeled_string(lfailure,ifile, 'hrdxyz',.true.
     &           , idname,fchar)
            if ( lfailure ) return
            call twh_fieldchar(GLB_SET,ifield,1,fchar)
            if ( (fchar .ne. 'x') .and.  (fchar .ne. 'y') .and.
     &           (fchar .ne. 'z') ) then
               write(6,*) 'hrdxyz must be x, y or z'
               lfailure = .true.
               return
            endif

c           --- center of hard wall
            dfield = twh_read_labeled_float(lfailure,ifile, 'hrdcen', 
     &         .true., idname)
            if ( lfailure ) return
            call twh_fielddp(GLB_SET,ifield,1,dfield)
            if ( dfield .lt. 0.0d0 ) then
               write(6,*) 'hrdcen must be non-negative'
               lfailure = .true.
               return
            endif

c           --- radius of hard wall
            dfield = twh_read_labeled_float(lfailure,ifile, 'hrdrad', 
     &         .true., idname)
            if ( lfailure ) return
            call twh_fielddp(GLB_SET,ifield,2,dfield)
            if ( dfield .lt. 0.0d0 ) then
               write(6,*) 'hrdrad must be non-negative'
               lfailure = .true.
               return
            endif

c           --- save the square of the distance
            dfield = dfield**2
            call twh_fielddp(GLB_SET,ifield,3,dfield)

c           --- hard wall repulsion style
            call twh_read_labeled_string(lfailure,ifile
     &           ,'hrd_repulsion_style', .true., idname,fchar)
            if ( lfailure ) return
            call twh_fieldchar(GLB_SET,ifield,2,fchar)
            if ( fchar .eq. 'centers' ) then
               write(6,*) 'Centers of all atoms are excluded from the'
     &              ,' hard wall'
            elseif ( fchar .eq. 'hard radii' ) then
               write(6,*)'Hard Sphere radii excluded from the hard wall'
               cplocal = twh_get_classical_potential()
               if ( cplocal .ne. 'Hard Sphere' .and.
     &              cplocal .ne. 'Repulsive Sphere' .and.
     &              cplocal .ne. 'Lennard-Jones') then
                  write(6,*) 'This hrd_repulsion_style is not valid'
     &                 ,' with this classical_potential'
                  lfailure = .true.
                  return
               endif
            endif

c           --- hard wall energy style
            call twh_read_labeled_string(lfailure,ifile
     &           ,'hrd_energy_type', .true., idname,fchar)
            if ( lfailure ) return
            call twh_fieldchar(GLB_SET,ifield,3,fchar)
            if ( fchar .eq. 'infinite' ) then
               write(6,*) 'Infinitely high energy inside of wall'
            elseif (fchar .eq.'finite' ) then
c              --- get the finite value
               dfield = twh_read_labeled_float(lfailure,ifile, 
     &            'hrd_wall_energy', .true., idname)
               if ( lfailure ) return
               call twh_fielddp(GLB_SET,ifield,ifile,dfield)
            else
               write(6,*) 'READ_FIELD: unknown hrd_energy_type'
               lfailure = .true.
               return
            endif

         elseif ( fieldstring .eq. 'Hard Ball' ) then
c           --- read hard ball information
c           --- box for hard ball
            ivalue = twh_read_labeled_integer(lfailure,ifile,'hrbbox'
     &           ,.true.,idname)
            if ( lfailure ) return
            call twh_fieldint(GLB_SET,ifield,FLD_HRBBOX,ivalue)
            if ( ivalue .lt. 1 ) then
               write(6,*) 'hrbbox must be positive'
               lfailure = .true.
               return
            endif
            if ( ivalue .gt. numboxes ) then
               write(6,*) 'hrbbox must be <= numboxes'
               lfailure = .true.
               return
            endif
c           --- hard ball position
            call twh_read_labeled_float_array(lfailure,ifile,'hrbpos',
     &           3,.true.,idname,fatmp)
            if ( lfailure ) return
            call twh_fielddp(GLB_SET,ifield,FLD_HRBPOSX,fatmp(1))
            call twh_fielddp(GLB_SET,ifield,FLD_HRBPOSY,fatmp(2))
            call twh_fielddp(GLB_SET,ifield,FLD_HRBPOSZ,fatmp(3))
c           --- hard ball radius
            dfield = twh_read_labeled_float(lfailure,ifile,'hrbrad'
     &           ,.true.,idname)
            if ( lfailure ) return
            call twh_fielddp(GLB_SET,ifield,FLD_HRBRAD,dfield)
            if ( dfield .lt. 0.0d0 ) then
               write(6,*) 'hrbrad must be non-negative'
               lfailure = .true.
               return
            endif
c           --- hard ball scale factor
            dfield = twh_read_labeled_float(lfailure,ifile,'hrbfac'
     &           ,.true.,idname)
            if ( lfailure ) return
            call twh_fielddp(GLB_SET,ifield,FLD_HRBFAC,dfield)
            if ( ( dfield .lt. 0.0d0 ) .OR.
     &           ( dfield .gt. 1.0d0 ) ) then
               write(6,*) 'hrbfac must be non-negative, smaller' 
     &              , 'than one'
               lfailure = .true.
               return
            endif
c           --- hard ball energy style
            call twh_read_labeled_string(lfailure,ifile
     &           ,'hrb_energy_type' , .true., idname,fchar)
            if ( lfailure ) return
            call twh_fieldchar(GLB_SET,ifield,1,fchar)
            if ( fchar .eq. 'infinite' ) then
               write(6,*) 'Infinitely high energy inside of wall'
            elseif (fchar .eq.'finite' ) then
c              --- get the finite value
               dfield = twh_read_labeled_float(lfailure,ifile, 
     &              'hrb_ball_energy', .true., idname)
               if ( lfailure ) return
               call twh_fielddp(GLB_SET,ifield,FLD_HRBENERGY,dfield)
            else
               write(6,*) 'READ_FIELD: unknown hrb_energy_type'
               lfailure = .true.
               return
            endif

         elseif ( fieldstring .eq. 'LJ 9-3 Wall' ) then

c           --- box for the lj field
            ivalue = twh_read_labeled_integer(lfailure,ifile,'ljfbox'
     &           ,.true.,idname)
            if ( lfailure ) return
            call twh_fieldint(GLB_SET,ifield,1,ivalue)
            if ( ivalue .gt. numboxes ) then
               write(6,*) 'ljfbox must be <= numboxes'
               lfailure = .true.
               return
            endif
            if ( ivalue .lt. 1 ) then
               write(6,*) 'ljfbox must be positive'
               lfailure = .true.
               return
            endif

c           --- x,y, or z dimension
            call twh_read_labeled_string(lfailure,ifile,'ljfxyz', .true.
     &           , idname,fchar)
            if ( lfailure ) return
            call twh_fieldchar(GLB_SET,ifield,1,fchar)
            if ( (fchar .ne. 'x') .and. (fchar .ne. 'y') .and.
     &           (fchar .ne. 'z') ) then
               write(6,*) 'ljfxyz must be x,y or z'
               lfailure = .true.
               return
            endif

c           --- location of the wall
            dfield = twh_read_labeled_float(lfailure,ifile,'ljfcen'
     &           , .true., idname)
            if ( lfailure ) return
            call twh_fielddp(GLB_SET,ifield,1,dfield)

c           --- direction wall faces
            ivalue = twh_read_labeled_integer(lfailure,ifile,'ljfdir'
     &           ,.true.,idname)
            if ( lfailure ) return
            call twh_fieldint(GLB_SET,ifield,2,ivalue)
            if ( iabs(ivalue) .ne. 1 ) then
               write(6,*) 'ljfdir must be -1 or 1'
               lfailure = .true.
               return
            endif

c           --- cutoff distance for wall interactions
            dfield = twh_read_labeled_float(lfailure,ifile,'ljfcut'
     &           , .true. , idname)
            if ( lfailure ) return
            call twh_fielddp(GLB_SET,ifield,2,dfield)
            if ( dfield .le. 0.0d0 ) then
               write(6,*) 'ljfcut must be positive'
               lfailure = .true.
               return
            endif

c           --- ljfshift: shift to zero at cutoff logical
            ltemp = twh_read_labeled_logical(lfailure,ifile,'ljfshift'
     &           , .true. , idname)
            if ( lfailure ) return
            call twh_fieldlogical(GLB_SET,ifield,ltemp)

c           --- ljfrho
            dfield = twh_read_labeled_float(lfailure,ifile,'ljfrho'
     &           , .true. , idname)
            if ( lfailure ) return
            call twh_fielddp(GLB_SET,ifield,3,dfield)
            if ( dfield .le. 0.0d0 ) then
               write(6,*) 'ljfrho must be positive'
               lfailure = .true.
               return
            endif

c           --- ljfntypes
            ivalue = twh_read_labeled_integer(lfailure,ifile
     &           ,'ljfntypes',.true.,idname)
            if ( lfailure ) return
            call twh_fieldint(GLB_SET,ifield,3,ivalue)
            if ( ivalue .lt. 1 ) then
               write(6,*) 'ljfntypes must be positive'
               lfailure = .true.
               return
            endif
            if ( ivalue .gt. NNTYPE ) then
               write(6,*) 'ljfntypes must be <= NNTYPE'
               lfailure = .true.
               return
            endif

c           --- read each type, resolve into the integer value for this 
c           --- atom name, and store the parameters
            do itype = 1,ivalue
c              --- ljfname
               call twh_read_labeled_string(lfailure,ifile,'ljfname'
     &              , .false. , idname,atomname)
               if ( lfailure ) return
               write(6,'(a5,i4,1x,a8,a10)') 'type:',itype
     &              ,'ljfname:',atomname

c              --- ljfsig
               ljfsig = twh_read_labeled_float(lfailure,ifile,
     &            'ljfsig', .false., idname)
               if ( lfailure ) return
               write(6,'(a5,i4,1x,a7,f10.5)') 'type:',itype
     &              ,'ljfsig:',ljfsig

c              --- ljfeps
               ljfeps = twh_read_labeled_float(lfailure,ifile,
     &            'ljfeps', .false., idname)
               if ( lfailure ) return
               write(6,'(a5,i4,1x,a7,f10.5)') 'type:',itype
     &              ,'ljfeps:',ljfeps

c              --- resolve ljfname into type numbers
               ntii = 0
               type = 0
               do while ( ntii .eq. 0 .and. type .lt. NNTYPE)
                  type = type + 1
                  if ( twh_get_nbname(type) .eq. atomname ) then
                     ntii = type
                  endif
               enddo
               if ( ntii .eq. 0 ) then
                  write(6,*) 'READ_FIELD: match not found'
                  lfailure = .true.
                  return
               endif
c              --- save the ntii value
               call twh_fieldnbint(GLB_SET,ifield,itype,1,ntii)
c              ---- convert ljfeps into ljfatr and ljfrep parts
               call twh_fielddp(GLB_GET,ifield,3,dfield)
               fnbparamone = twopi*ljfeps*dfield*2.0d0*(ljfsig**12)
     &              /(45.0d0)
               call twh_fieldnbparam(GLB_SET,ifield,ntii,1,fnbparamone)
               fnbparamtwo = -twopi*ljfeps*dfield*(ljfsig**6)/3.0d0
               call twh_fieldnbparam(GLB_SET,ifield,ntii,2,fnbparamtwo)
               call twh_fieldlogical(GLB_GET,ifield,ltemp)
               if ( ltemp ) then
c                 --- using a shifted potential, reset ljfadj
c                 --- set coordinates to the lj cutoff
                  call twh_fielddp(GLB_GET,ifield,2,dfield)
                  invdist = 1.0d0/dfield
                  invdist3 = invdist*invdist*invdist
                  invdist9 = invdist3*invdist3*invdist3
                  fnbparam = -(
     &                 + fnbparamone*invdist9 
     &                 + fnbparamtwo*invdist3)
                  call twh_fieldnbparam(GLB_SET,ifield,ntii,3
     &                 ,fnbparam)
               endif
c              --- store ljfsig and ljseps
               call twh_fieldnbparam(GLB_SET,ifield,ntii,4,ljfsig)
               call twh_fieldnbparam(GLB_SET,ifield,ntii,5,ljfeps)
            enddo

         elseif ( fieldstring .eq. 'Hooper Umbrella' ) then
c           --- read in umbbox
            ivalue = twh_read_labeled_integer(lfailure,ifile,'umbbox'
     &           , .true., idname)
            if ( lfailure ) return
            call twh_fieldint(GLB_SET,ifield,1,ivalue)
            if ( ivalue .lt. 1 .or. ivalue .gt. numboxes ) then
               write(6,*) 'umbxyz must range from 1 to numboxes'
               lfailure = .true.
               return
            endif

c           --- read in umbxyz
            call twh_read_labeled_string(lfailure,ifile,'umbxyz',.true.
     &           , idname,fchar)
            if ( lfailure ) return
            call twh_fieldchar(GLB_SET,ifield,1,fchar)
            if ( (fchar .ne. 'x') .and. (fchar .ne. 'y') .and.
     &           (fchar .ne. 'z') ) then
               write(6,*) 'umbxyz must be x, y, or z'
               lfailure = .true.
               return
            endif

c           --- read in umbcenter
            dfield = twh_read_labeled_float(lfailure,ifile,'umbcenter'
     &           ,.true. ,idname)
            if ( lfailure ) return
            call twh_fielddp(GLB_SET,ifield,1,dfield)
            if ( dfield .le. 0.0d0 ) then
               write(6,*) 'umbcenter must be positive'
               lfailure = .true.
               return
            endif

c           --- read in umba
            dfield = twh_read_labeled_float(lfailure,ifile,'umba'
     &           , .true. , idname)
            if ( lfailure ) return
            call twh_fielddp(GLB_SET,ifield,2,dfield)

         elseif ( fieldstring .eq. 'Steele Wall' ) then
c           --- box for the lj field
            ivalue = twh_read_labeled_integer(lfailure,ifile
     &           ,'steele box', .true.,idname)
            if ( lfailure ) return
            call twh_fieldint(GLB_SET,ifield,1,ivalue)
            if ( ivalue .gt. numboxes ) then
               write(6,*) 'steele box must be <= numboxes'
               lfailure = .true.
               return
            endif
            if ( ivalue .lt. 1 ) then
               write(6,*) 'steele box must be positive'
               lfailure = .true.
               return
            endif

c           --- x,y, or z dimension
            call twh_read_labeled_string(lfailure,ifile,'steele xyz'
     &           ,.true. ,idname,fchar)
            if ( lfailure ) return
            call twh_fieldchar(GLB_SET,ifield,1,fchar)
            if ( (fchar .ne. 'x') .and. (fchar .ne. 'y') .and.
     &           (fchar .ne. 'z') ) then
               write(6,*) 'steele xyz must be x,y or z'
               lfailure = .true.
               return
            endif

c           --- location of the wall
            dfield = twh_read_labeled_float(lfailure,ifile
     &           ,'steele surface' ,.true.,idname)
            if ( lfailure ) return
            call twh_fielddp(GLB_SET,ifield,1,dfield)

c           --- direction wall faces
            ivalue = twh_read_labeled_integer(lfailure,ifile
     &           ,'steele dir',.true.,idname)
            if ( lfailure ) return
            call twh_fieldint(GLB_SET,ifield,2,ivalue)
            if ( iabs(ivalue) .ne. 1 ) then
               write(6,*) 'steele dir must be -1 or 1'
               lfailure = .true.
               return
            endif

c           --- cutoff distance for wall interactions
            dfield = twh_read_labeled_float(lfailure,ifile
     &           ,'steele cutoff',.true.  ,idname)
            if ( lfailure ) return
            call twh_fielddp(GLB_SET,ifield,2,dfield)
            if ( dfield .le. 0.0d0 ) then
               write(6,*) 'steele cutoff must be positive'
               lfailure = .true.
               return
            endif

c           --- shift to zero at cutoff logical
            ltemp = twh_read_labeled_logical(lfailure,ifile
     &           ,'steele shift' , .true.,idname)
            if ( lfailure ) return
            call twh_fieldlogical(GLB_SET,ifield,ltemp)

c           --- steele delta
            dfield = twh_read_labeled_float(lfailure,ifile
     &           ,'steele delta', .true. ,idname)
            if ( lfailure ) return
            call twh_fielddp(GLB_SET,ifield,3,dfield)
            if ( dfield .le. 0.0d0 ) then
               write(6,*) 'steele delta must be positive'
               lfailure = .true.
               return
            endif

c           --- steele rho_s
            dfield = twh_read_labeled_float(lfailure,ifile
     &           ,'steele rho_s', .true.  , idname)
            if ( lfailure ) return
            call twh_fielddp(GLB_SET,ifield,4,dfield)
            if ( dfield .le. 0.0d0 ) then
               write(6,*) 'steele rho_s must be positive'
               lfailure = .true.
               return
            endif

c           --- steele ntype
            ivalue = twh_read_labeled_integer(lfailure,ifile
     &           ,'steele ntype', .true. , idname)
            if ( lfailure ) return
            call twh_fieldint(GLB_SET,ifield,3,ivalue)
            if ( ivalue .lt. 1 ) then
               write(6,*) 'steele ntype must be positive'
               lfailure = .true.
               return
            endif
            if ( ivalue .gt. NNTYPE ) then
               write(6,*) 'steele ntype must be <= NNTYPE'
               lfailure = .true.
               return
            endif

c           --- read each type, resolve into the integer value for this 
c           --- atom name, and store the parameters
            do itype = 1,ivalue
c              --- steele name
               call twh_read_labeled_string(lfailure,ifile,
     &            'steele name', .false., idname,atomname)
               if ( lfailure ) return
               write(6,'(a5,i4,1x,a12,a10)') 'type:',itype
     &              ,'steele name:',atomname

c              --- sigma_sf
               sigma_sf = twh_read_labeled_float(lfailure,ifile,
     &            'sigma_sf', .false., idname)
               if ( lfailure ) return
               write(6,'(a5,i4,1x,a9,f10.5)') 'type:',itype
     &              ,'sigma_sf:',sigma_sf

c              --- epsilon_sf
               epsilon_sf = twh_read_labeled_float(lfailure,ifile,
     &            'epsilon_sf', .false., idname)
               if ( lfailure ) return
               write(6,'(a5,i4,1x,a11,f10.5)') 'type:',itype
     &              ,'epsilon_sf:',epsilon_sf

c              --- resolve atomname into type numbers
               ntii = 0
               type = 0
               do while ( ntii .eq. 0 .and. type .lt. NNTYPE)
                  type = type + 1
                  if ( twh_get_nbname(type) .eq. atomname ) then
                     ntii = type
                  endif
               enddo
               if ( ntii .eq. 0 ) then
                  write(6,*) 'READ_FIELD: match not found'
                  lfailure = .true.
                  return
               endif
c              --- save the ntii value
               call twh_fieldnbint(GLB_SET,ifield,itype,1,ntii)
c              ---- convert into epsilon_w
c              --- e_w = 2 Pi e_sf rho_s (sigma_sf)^2 delta
               call twh_fielddp(GLB_GET,ifield,3,fthree)
               call twh_fielddp(GLB_GET,ifield,4,ffour)
               fnbparamone = twopi*epsilon_sf*ffour*(sigma_sf**2)*fthree
               call twh_fieldnbparam(GLB_SET,ifield,ntii,1,fnbparamone)
c              --- store sigma_sf in the nb arrays
               call twh_fieldnbparam(GLB_SET,ifield,ntii,2,sigma_sf)
               fnbparam = sigma_sf**4
               call twh_fieldnbparam(GLB_SET,ifield,ntii,3,fnbparam)
               call twh_fieldlogical(GLB_GET,ifield,ltemp)
               if ( ltemp ) then
c                 --- using a shifted potential, set the adjustment
c                 --- using the steele cutoff
                  call twh_fielddp(GLB_GET,ifield,2,zzz)
                  szone = sigma_sf/zzz
                  sztwo = szone*szone
                  szfour = sztwo*sztwo
                  szten = szfour*szfour*sztwo
                  fnbparam = -(fnbparamone*(
     &                 (2.0d0/5.0d0)*szten - szfour
     &                 - fnbparam/
     &                 ( 3.0d0*fthree*(zzz + 0.61d0*fthree)**3)))
                  call twh_fieldnbparam(GLB_SET,ifield,ntii,4,fnbparam)
               endif
c              --- store epsilon_sf in the nb arrays
               call twh_fieldnbparam(GLB_SET,ifield,ntii,5,epsilon_sf)
            enddo

c        --- cdl: harmonic field
         elseif ( fieldstring .eq. 'Harmonic Attractor' ) then
c           --- read harmonic attractor information
c           --- box for the harmonic attractor
            ivalue = twh_read_labeled_integer(lfailure,ifile,'hafbox'
     &           ,.true.,idname)
            if ( lfailure ) return
            call twh_fieldint(GLB_SET,ifield,1,ivalue)
            if ( ivalue .gt. numboxes) then
               write(6,*) 'hafbox must be <= numboxes'
               lfailure = .true.
               return
            endif
            if ( ivalue .lt. 1) then
               write(6,*) 'hafbox must be positive'
               lfailure = .true.
               return
            endif

c           --- force constant for the harmonic attractor
            dfield = twh_read_labeled_float(lfailure,ifile,'hafk'
     &           , .true. , idname)
            if ( lfailure ) return
            call twh_fielddp(GLB_SET,ifield,1,dfield)

c           --- number of atoms to apply harmonic attractor
            ivalue = twh_read_labeled_integer(lfailure,ifile
     &           ,'hafnentries',.true.,idname)
            if ( lfailure ) return
            call twh_fieldint(GLB_SET,ifield,2,ivalue)
            if ( ivalue .lt. 0) then
               write(6,*) 'hafnentries must be >= 0'
               lfailure = .true.
               return
            endif
            if ( ivalue .gt. (NNTYPE*NTMAX)) then
               write(6,*) 'hafnetries must be <= NNTYPE*NTMAX'
               lfailure = .true.
               return
            endif

c           --- reference position type for atoms
            call twh_read_labeled_string(lfailure,ifile,'hafrefpos'
     &           ,.true. ,idname,fchar)
            if ( lfailure ) return
            call twh_fieldchar(GLB_SET,ifield,1,fchar)
            if ( fchar .eq. 'Global') then

               call twh_read_labeled_float_array(lfailure,ifile
     &              ,'hafglobxyz', 3, .true., idname, fatmp)
               if ( lfailure ) return

               call twh_fielddp(GLB_SET,ifield,2,fatmp(1))
               call twh_fielddp(GLB_SET,ifield,3,fatmp(2))
               call twh_fielddp(GLB_SET,ifield,4,fatmp(3))

               if (fatmp(1) .lt. 0.0d0) then
                  write(6,*) 'hafglobx must be >= 0.0'
                  lfailure = .true.
                  return
               endif
               if (fatmp(2) .lt. 0.0d0) then
                  write(6,*) 'hafgloby must be >= 0.0'
                  lfailure = .true.
                  return
               endif
               if (fatmp(3) .lt. 0.0d0) then
                  write(6,*) 'hafglobz must be >= 0.0'
                  lfailure = .true.
                  return
               endif
            elseif ( fchar .ne. 'Initial') then
               write(6,*) 'READ_FIELD: unknown hafrefpos'
               lfailure = .true.
               return
            endif

c           --- key type for the harmonic attractor
            call twh_read_labeled_string(lfailure,ifile,'hafkey',.true.
     &           ,idname,fchar)
            if ( lfailure ) return
            call twh_fieldchar(GLB_SET,ifield,2,fchar)
            call twh_fieldint(GLB_GET,ifield,2,ivalue)
            do ientry = 1,ivalue
               hafmolec = twh_read_labeled_integer(lfailure,ifile,
     &              'hafmolec', .true., idname)
               if ( lfailure ) return
               if (hafmolec .lt. 0) then
                  write(6,*) 'hafmolec must be >= 0'
                  lfailure = .true.
                  return
               endif
               if (hafmolec .gt. NTMAX) then
                  write(6,*) 'hafmolec must be <= NTMAX'
                  lfailure = .true.
                  return
               endif

               if ( fchar .eq. 'Element') then
c                 --- fix a certain element type in a given molecule
c                 --- read atom name
                  call twh_read_labeled_string(lfailure,ifile,
     &               'hafelement', .false., idname,atomname)
                  if ( lfailure ) return
                  write(6,'(a5,i4,1x,a11,a2)') 'type:',ientry
     &              ,'hafelement:',atomname

c                 --- resolve atomname into type numbers
                  ntii = 0
                  type = 0
                  do while ( ntii .eq. 0 .and. type .lt. NNTYPE)
                     type = type + 1
                     if ( twh_get_elementname(type) 
     &                    .eq. atomname(1:2) ) then
                        ntii = type
                     endif
                  enddo
               elseif ( fchar .eq. 'FFtype') then

c                 --- fix a certain forcefield type in a given molecule
c                 --- read atom name
                  call twh_read_labeled_string(lfailure,ifile,
     &               'hafname', .false., idname,atomname)
                  if ( lfailure ) return
                  write(6,'(a5,i4,1x,a8,a10)') 'type:',ientry
     &              ,'hafname:',atomname

c                 --- resolve atomname into type numbers
                  ntii = 0
                  type = 0
                  do while ( ntii .eq. 0 .and. type .lt. NNTYPE)
                     type = type + 1
                     if ( twh_get_nbname(type) .eq. atomname ) then
                        ntii = type
                     endif
                  enddo
               else
                  write(6,*) 'READ_FIELD: unknown hafkey'
                  lfailure = .true.
                  return
               endif 
               if ( ntii .eq. 0 ) then
                  write(6,*) 'READ_FIELD: match not found'
                  lfailure = .true.
                  return
               endif
c              --- save the ntii values
               call twh_fieldnbint(GLB_SET,ifield,ientry,1,hafmolec)
               call twh_fieldnbint(GLB_SET,ifield,ientry,2,ntii)
            enddo
         else
            write(6,*) 'READ_FIELD: unknown fieldtype'
            lfailure = .true.
            return
         endif
      enddo

      return
      end
