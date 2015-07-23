#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2006-2007 Marcus G. Martin                       *
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

c     ******************************************************************
c     * This file contains functions and variables associated with     *
c     * scaling of the Hamiltonian by a parameter(s) lambda,           *
c     * specifically for the 'Scaled Lennard-Jones' classical potential*
c     *                                                                *
c     * The atoms for which scaling takes place are defined by the     *
c     * cmix_pair_list parameters in towhee_input.                     *
c     *                                                                *
c     * The parameters lambda_lj and alpha_lj scale the lennard-jones  *
c     * terms (see vtwobody.F and wtwobody.F).  The parameter lambda_c *
c     * scales the coulombic terms (see vcoulomb.F and wcoulomb.F).    *
c     *                                                                *
c     * We may evaluate the derivative of the potential energy with    *
c     * respect to lambda_c and lambda_lj.                             *
c     *                                                                *
c     * We may also evaluate the potential energy of the current       *
c     * atomic configuration at a different value of lambda; such      *
c     * a parameter, and corresponding energy, are termed 'foreign',   *
c     * to distinguish between 'native' parameters which are the ones  *
c     * used in move acceptance criteria.                              *
c     *                                                                *
c     * While most of these data structures are specific to the Scaled *
c     * Lennard Jones potential, Weeks-Chandler-Anderson (WCA) uses    *
c     * atom scaling data structures to keep track of repulsive-only   *
c     * atom types                                                     *
c     *                                                                *
c     * To use, first call scp_initialize() during initialization      *
c     * call scp_set_scaled_atoms(...) to establish which atoms will   *
c     *  be scaled.                                                    *
c     * when force field information has been read and after the       *
c     *  above call, call scp_initialize_atom_scaling(...).            *
c     * Thereafter, calling scp_scale_atoms(ntii,ntjj) will tell you   *
c     *  whether interaction between two atom types are scaled.        *
c     *                                                                *
c     *                                                                *
c     * Functions begin with scp_ for SCaled Potential.                *
c     *                                                                *
c     * The structure of the file scaled_potential.F is in the spirit  *
c     * of an 'object' data structure.  There are globally defined     *
c     * variables, locally defined variables for use only by the       *
c     * functions/subroutines in this file, and globally accessible    *
c     * functions/subroutines which operate on these variables (see    *
c     * below).                                                        *
c     ******************************************************************

c     FUNCTIONS AND SUBROUTINES:
c     Initialization:
c      subroutine scp_initialize()
c      subroutine scp_set_scaled_atoms(npair, pair_list)
c      subroutine scp_initialize_atom_scaling(ntypes)
c     The above three must be called in that order before 
c     scp_scale_atoms will work
c
c     Scaling information:
c      function scp_scale_atoms( ntii, ntjj )
c      subroutine scp_get_scaled_atoms(npair, pair_list)
c      subroutine scp_set_scaling_style(style)
c      function scp_get_scaling_style()
c      function scp_get_scaling_style_string()
c     
c     Manipulation of native and foreign parameters 
c      subroutine scp_set_num_foreign_lambda(num_foreign_lambda)
c      subroutine scp_set_foreign_lambdas(param, v)
c      function scp_get_foreign_lambda(param, l)
c      function scp_get_native_lambda(param)
c      subroutine scp_set_native_lambda(param, v, param_name)
c
c     Manipulation of foreign energies
c      subroutine scp_set_foreign_energies(box, v)
c      subroutine scp_increment_foreign_energies(box, v)
c      subroutine scp_get_foreign_energies(box, v)
c
c     Manipulation of derivatives of energy wrt lambdas
c      subroutine scp_set_derivative(box, param, v)
c      subroutine scp_increment_derivative(box, param, v)
c      function scp_get_derivative(box, param)
c
c     Internal utility functions
c      subroutine scp_check_param(param, non_lambda_ok, l, fname)


c  initial implementation of dynamic memory allocation via C calls.
c  The arrays foreign_lambda are stored in
c  memory structures defined in scaled_potential_data.c if the 
c  compiler flag USE_DYNAMIC_MEMORY is defined
c  This is *not* yet implemented


      subroutine scp_initialize()
c     ******************************************************************
c     * initialize various variables internal to the scaled_potential  *
c     * family of functions.                                           *
c     *                                                                *
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 10-13-2008 by M.G. Martin                          *
c     ******************************************************************
      implicit none
c     --- local variables
      integer i

      call twh_scp_scale_coulomb(GLB_SET,.false.)
      call twh_scp_scale_twobody(GLB_SET,.false.)
      call twh_scp_eval_derivative(GLB_SET,.false.)

      call twh_scaling_style(GLB_SET,SCP_NONE)

c     -- defaults for params below reproduce non-scaled potentials
      call twh_native_param(GLB_SET,LAMBDA_LJ,1.0d0)
      call twh_native_param(GLB_SET,LAMBDA_C,1.0d0)
      call twh_native_param(GLB_SET,ALPHA_LJ,0.0d0)

#ifndef USE_DYNAMIC_MEMORY
      call twh_foreign_energy(GLB_INIT,MAXBOX,MAX_FOREIGN_LAMBDA,0.0d0)
      call twh_foreign_lambda(GLB_INIT,MAX_FOREIGN_LAMBDA,LAMBDA_LJ
     &     ,0.0d0)
      call twh_foreign_lambda(GLB_INIT,MAX_FOREIGN_LAMBDA,LAMBDA_C
     &     ,0.0d0)
#endif
      call twh_scp_num_foreign_lambda(GLB_SET,0)

      call twh_scaled_pair_list(GLB_INIT,CROSSTYPEMAX,2,'')
      call twh_scaled_npair(GLB_SET,0)

      
      do i=1,MAXBOX
         call twh_du_dlambda(GLB_SET,i,LAMBDA_LJ,0.0d0)
         call twh_du_dlambda(GLB_SET,i,LAMBDA_C,0.0d0)
      enddo

      call twh_lscaleatom(GLB_INIT,NNTYPE,.false.)

      return
      end


      subroutine scp_set_num_foreign_lambda(num_foreign_lambda)
c     ******************************************************************
c     * Set the number of foreign lambdas.  For the case of non-dynamic*
c     * data structures, this just sets the scp_num_foreign_lambda     *
c     * variable, but for the case of dynamic data, it does memory     *
c     * allocation.  Once called, number of foreign lambdas may not be *
c     * changed.                                                       *
c     * num_foreign_lambdas -- how many foreign lambdas we have        *
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 10-13-2008 by M.G. Martin                          *
c     ******************************************************************
      implicit none
      integer num_foreign_lambda
c     --- local variables
      integer scp_num_foreign_lambda
   
c     Error checks
c     First, make sure called only once
      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
      if (scp_num_foreign_lambda .ne. 0) then
         write(6,*) "scp_set_num_foreign_lambdas may only be called ",
     &         "once"
         stop
      endif

#ifndef USE_DYNAMIC_MEMORY
      if (num_foreign_lambda .gt. MAX_FOREIGN_LAMBDA) then
         write(6,*) 'too many foreign lambdas'
         stop
      endif
#else
c     testing only
c      double precision v

c      i = 123
c      call initializec(i)
c      i = -12
c      call get_number(i)
c      write(6,*) "in fortran, got ", i
c      call initialize_foreign_lambda(3)
c      call get_foreign_lambda(3, v)
c      write(6,*) "in fortran, got ", v
c      stop

#endif

      call twh_scp_num_foreign_lambda(GLB_SET,num_foreign_lambda)

      return
      end


      subroutine scp_set_foreign_lambdas(param, v)
c     ******************************************************************
c     * Set the values of foreign lambdas.                             *
c     * param -- parameter name we are setting: LAMBDA_C or LAMBDA_LJ  *
c     * v -- array of parameter values (scp_num_foreign_lambda of them)*
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 10-13-2008 by M.G. Martin                          *
c     ******************************************************************
      implicit none
      integer param
      double precision v(MAX_FOREIGN_LAMBDA)
c     --- local variables
      integer l,scp_num_foreign_lambda
   
c     --- error checks
      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
      if (scp_num_foreign_lambda .gt. MAX_FOREIGN_LAMBDA) then
         write(6,*) 'too many foreign lambdas'
         stop
      endif
      call scp_check_param(param, .false., -1, 
     &      "scp_set_foreign_lambdas")

      do l=1,scp_num_foreign_lambda
         call twh_foreign_lambda(GLB_SET,l,param,v(l))
      enddo
      return
      end


      double precision function scp_get_foreign_lambda(param, l)
c     ******************************************************************
c     * returns the requested foreign lambda                           *
c     * param -- parameter name we are setting: LAMBDA_C or LAMBDA_LJ  *
c     * l -- the index of the requested parameter                      *
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 12-22-2007 by M.G. Martin                          *
c     ******************************************************************
      implicit none
      integer param,l
      double precision dvalue

      call scp_check_param(param, .false., l, "scp_get_foreign_lambdas")
      call twh_foreign_lambda(GLB_GET,l,param,dvalue)
      scp_get_foreign_lambda = dvalue
      return 
      end


      double precision function scp_get_native_lambda(param)
c     ******************************************************************
c     * Returns the requested lambda/parameter value                   *
c     * even though this function ends with lambda, it can also returns*
c     * the parameter ALPHA_LJ                                         *
c     * param -- LAMBDA_C, LAMBDA_LJ, or ALPHA_LJ                      *
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 12-20-2007 by M.G. Martin                          *
c     ******************************************************************
      implicit none
      integer param
      call scp_check_param(param, .true., -1, "scp_get_native_lambda")
      call twh_native_param(GLB_GET,param,scp_get_native_lambda)
      return
      end


      subroutine scp_set_native_lambda(param, v, param_name)
c     ******************************************************************
c     * set the requested lambda/parameter to given value              *
c     * The values of v must be between 0 and 1                        *
c     * param -- LAMBDA_C, LAMBDA_LJ, or ALPHA_LJ                      *
c     * v -- the value                                                 *
c     * param_name -- text name of parameter, for user output.         *
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 12-20-2007 by M.G. Martin                          *
c     ******************************************************************
      implicit none
      integer param
      double precision v
      character*(*) param_name

      call scp_check_param(param, .true., -1, "scp_set_native_lambda")

      if ( v .lt. 0.0d0 .or. v .gt. 1.0d0) then
         write(6,*) param_name, ' must be between 0 and 1'
         stop
      endif

      call twh_native_param(GLB_SET,param,v)
      return
      end


      subroutine scp_set_foreign_energies(box, v)
c     ******************************************************************
c     * Set the total classical energies of the given box to the given *
c     * values. This subroutine is safe even when there are no foreign *
c     * lambdas.                                                       *
c     * Units for energy are in [K]                                    *
c     * box is the box of these energies                               *
c     * v is a vector containing scp_num_foreign_lambda entries        *
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 10-13-2008 by M.G. Martin                          *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer box
      double precision v(MAX_FOREIGN_LAMBDA)
c     --- local variables
      integer l,scp_num_foreign_lambda

      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
      if (scp_num_foreign_lambda .eq. 0) return
      do l = 1,scp_num_foreign_lambda
         call twh_foreign_energy(GLB_SET,box,l,v(l))
      enddo
      return
      end


      subroutine scp_increment_foreign_energies(box, v)
c     ******************************************************************
c     * Adds the given energies to the total classical energies in the *
c     * given box. This subroutine is safe even when there are no      *
c     * foreign lambdas.                                               *
c     * Units for energy are in [K]                                    *
c     * box is the box of these energies                               *
c     * v is a vector containing scp_num_foreign_lambda entries        *
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 10-13-2008 by M.G. Martin                          *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer box
      double precision v(MAX_FOREIGN_LAMBDA)
c     --- local variables
      integer l,scp_num_foreign_lambda
      
      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
      if (scp_num_foreign_lambda .eq. 0) return
      do l = 1,scp_num_foreign_lambda
         call twh_foreign_energy(GLB_INCR,box,l,v(l))
      enddo
      return
      end


      subroutine scp_get_foreign_energies(box, v)
c     ******************************************************************
c     * Returns the foreign energies for the requested box in the      *
c     * vector v.  Safe to use if no foreign lambdas.                  *
c     * Units for energy are in [K]                                    *
c     * box  is the box of the requested energies                      *
c     * v is a vector, will be filled with scp_num_foreign_lambda value*
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 10-13-2008 by M.G. Martin                          *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer box
      double precision v(MAX_FOREIGN_LAMBDA)
c     --- local variables
      integer l,scp_num_foreign_lambda

      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
      if (scp_num_foreign_lambda .eq. 0) return
      do l = 1,scp_num_foreign_lambda
         call twh_foreign_energy(GLB_GET,box,l,v(l))
      enddo
      return
      end


      subroutine scp_set_derivative(box, param, v)
c     ******************************************************************
c     * Set the value of the derivative of energy with respect to the  *
c     * given parameter (e.g. dU/d[lambda_c]) for the given box        *
c     * box is the box under consideration                             *
c     * Units are in [K]                                               *
c     * param is either LAMBDA_C or LAMBDA_LJ                          *
c     * v is the value of the derivative                               *
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 12-22-2007 by M.G. Martin                          *
c     ******************************************************************
      implicit none
      integer box,param
      double precision v

      call scp_check_param(param, .false., -1, "scp_set_derivative")
      call twh_du_dlambda(GLB_SET,box,param,v)

      return
      end


      subroutine scp_increment_derivative(box, param, v)
c     ******************************************************************
c     * Add v to the the derivative of energy with respect to the      *
c     * given parameter (e.g. dU/d[lambda_c]) for the given box        *
c     * box is the box under consideration                             *
c     * Units are in [K]                                               *
c     * param is either LAMBDA_C or LAMBDA_LJ                          *
c     * v is the to be added                                           *
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 12-22-2007 by M.G. Martin                          *
c     ******************************************************************
      implicit none
      integer box, param
      double precision v

      call scp_check_param(param, .false., -1, 
     &      "scp_increment_derivative")
      call twh_du_dlambda(GLB_INCR,box,param,v)

      return
      end


      double precision function scp_get_derivative(box, param)
c     ******************************************************************
c     * Return the derivative of energy with respect to the given      *
c     * parameter for the given box.                                   *
c     * box is the box under consideration                             *
c     * Units are in [K]                                               *
c     * param is either LAMBDA_C or LAMBDA_LJ                          *
c     * function returns the derivative                                *
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 09-25-2008 by M.G. Martin                          *
c     ******************************************************************
      implicit none
      integer box, param,potentialstyle
c     --- retrieve constants
      call twh_potentialstyle(GLB_GET,potentialstyle)

      if ( potentialstyle .eq. POT_INTERNAL ) then
         call scp_check_param(param, .false., -1, "scp_get_derivative")
         call twh_du_dlambda(GLB_GET,box,param,scp_get_derivative)
      else
c        --- not defined for quantum
         scp_get_derivative = 0.0d0
      endif
      return
      end


      subroutine scp_initialize_atom_scaling(ntypes)
c     ******************************************************************
c     * mark which atoms will be rescaled.  Rescaling takes place      *
c     * (currently) only for the Scaled Lennard-Jones and              *
c     * Weeks-Chandler-Anderson classical potentials.                  *
c     * Atoms are scaled according to their type and forcefield.       *
c     * This function must be called only after forcefield information *
c     * has been read by rwforcefield(), and after calling             *
c     * scp_set_scaled_atoms.                                          *
c     * ntypes is (I think...) the number of atom types we have; set   *
c     *     by rwforcefield                                            *
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 12-22-2007 by M.G. Martin                          *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NBFF
#define FUNCTION_GET_NBNAME
#include "functions.h"

c     arguments passed to subroutine
      integer ntypes
c     local variables
      logical lscale
      integer itype, ipair,nscale
      character*10 scaleone,scaletwo

c     --- Return early if scaled_npair not set; this may happen for
c     --- different potentials, or if scp_set_scaled_atoms has not yet
c     --- been called (the latter case being an error)
      call twh_scaled_npair(GLB_GET,nscale)
      if ( nscale .eq. 0 ) return

      do itype = 1,ntypes
c             --- check to see if we are using special scale factors
c             --- for this atom
         lscale = .false.
         do ipair = 1,nscale
            call twh_scaled_pair_list(GLB_GET,ipair,1,scaleone)
            call twh_scaled_pair_list(GLB_GET,ipair,2,scaletwo)
c           --- check for matches
            if ( scaleone .eq. twh_get_nbff(itype) 
     &           .and.
     &           scaletwo .eq. twh_get_nbname(itype)
     &           ) then
               lscale = .true.
            endif
         enddo
         call twh_lscaleatom(GLB_SET,itype,lscale)
      enddo

      return
      end


      subroutine scp_set_scaled_atoms(npair, pair_list)
c     ******************************************************************
c     * Define the atoms which are to be scaled.                       *
c     * npair (aka cmix_npair) is the number of atom types to be scaled*
c     *     (with Scaled Lennard-Jones potential)                      *
c     * pair_list (aka cmix_pair_list) is as defined in the towhee     *
c     *     manual; an array size (npair,2) with the force field name  *
c     *     and atom type in the first and second fields, respectively.*
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 12-22-2007 by M.G. Martin                          *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer npair,l
      character*10 pair_list
      dimension pair_list(CROSSTYPEMAX,2)
c     --- local variables
      character*10 plist

      if ( npair .gt. CROSSTYPEMAX ) then
c     --- in practice, this ought to be caught before we get here, but 
c     --- be safe anyhow
         write(6,*) 'npair exceeds CROSSTYPEMAX in scp_set_scaled_atoms'
         stop
      endif

      do l=1,npair
         plist = pair_list(l,1)
         call twh_scaled_pair_list(GLB_SET,l,1,plist)
         plist = pair_list(l,2)
         call twh_scaled_pair_list(GLB_SET,l,2,plist)
      enddo
      call twh_scaled_npair(GLB_SET,npair)

      return
      end

      subroutine scp_get_scaled_atoms(npair, pair_list)
c     ******************************************************************
c     * Get the atoms which are to be scaled.  Returned values:        *
c     * npair (aka cmix_npair) is the number of atom types to be scaled*
c     *     (with Scaled Lennard-Jones potential)                      *
c     * pair_list (aka cmix_pair_list) is as defined in the towhee     *
c     *     manual; an array size (npair,2) with the force field name  *
c     *     and atom type in the first and second fields, respectively.*
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 12-22-2007 by M.G. Martin                          *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer npair
      character*10 pair_list
      dimension pair_list(CROSSTYPEMAX,2)
c     --- local variables
      integer l,nscale
      character*10 plist

      call twh_scaled_npair(GLB_GET,nscale)
      do l=1,nscale
         call twh_scaled_pair_list(GLB_GET,l,1,plist)
         pair_list(l,1) = plist
         call twh_scaled_pair_list(GLB_GET,l,2,plist)
         pair_list(l,2) = plist
      enddo
      npair = nscale

      return
      end


      subroutine scp_set_scaling_style(style)
c     ******************************************************************
c     * Set the scaling style.  Currently legal values:                *
c     *   none                                                         *
c     *   soft-core                                                    *
c     *   wca                                                          *
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 12-22-2007 by M.G. Martin                          *
c     ******************************************************************
      implicit none
      character*30 style

c     the scaling_style is an integer taking on the values
c     defined in preproc.h.  The string version is defined in the
c     corresponding _STR macros.
      if (style .eq. SCP_NONE_STR) then
c        'none'
         call twh_scaling_style(GLB_SET,SCP_NONE)
      elseif (style .eq. SCP_SC_STR) then
c        'soft-core'
         call twh_scaling_style(GLB_SET,SCP_SC)
      elseif (style .eq. SCP_WCA_STR) then
c        'WCA'
         call twh_scaling_style(GLB_SET,SCP_WCA)
      else
         write(6,*) 'unknown cmix_rescaling_style: ',style
         stop
      endif

      return
      end


      integer function scp_get_scaling_style()
c     ******************************************************************
c     * returns the integer value of scaling style.  See preproc.h for *
c     * macro definitions.                                             *
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 12-22-2007 by M.G. Martin                          *
c     ******************************************************************
      implicit none
      call twh_scaling_style(GLB_GET,scp_get_scaling_style)
      return
      end


      character*30 function scp_get_scaling_style_string()
c     ******************************************************************
c     * returns the string value of scaling style.  See preproc.h for  *
c     * macro definitions.                                             *
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 12-22-2007 by M.G. Martin                          *
c     ******************************************************************
      implicit none
c     --- local variables
      integer style

      call twh_scaling_style(GLB_GET,style)
      if ( style .eq. SCP_NONE) then
         scp_get_scaling_style_string = SCP_NONE_STR
      elseif ( style .eq. SCP_WCA) then
         scp_get_scaling_style_string = SCP_WCA_STR
      elseif ( style .eq. SCP_SC) then
         scp_get_scaling_style_string = SCP_SC_STR
c     --- weird, something went bad in set_scaling_style().
      else 
         scp_get_scaling_style_string = 'UNKNOWN'
      endif
      return
      end


      logical function scp_scale_atoms( ntii, ntjj )
c     ******************************************************************
c     * returns true if the interactions between atoms ntii and ntjj   *
c     * are to be scaled.                                              *
c     * atom types given by, e.g.                                      *
c     *      ntii = twh_get_ntype(imolty,iunit)                        *
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 12-22-2007 by M.G. Martin                          *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer ntii, ntjj
c     --- local variables
      logical lscaleii,lscalejj

c     --- scaling is done if one, but not both, atoms are scalable
      call twh_lscaleatom(GLB_GET,ntii,lscaleii)
      call twh_lscaleatom(GLB_GET,ntjj,lscalejj)
      scp_scale_atoms = ( lscaleii .neqv. lscalejj ) 
      return
      end


      subroutine scp_check_param(param, non_lambda_ok, l, fname)
c     ******************************************************************
c     * an internal convenience function to make sure that passed      *
c     * params are legal.                                              *
c     * param -- see below for legal values                            *
c     * l -- index of requested parameter; must .le. num_foreign_lambda*
c     *      -1 to skip check                                          *
c     * if non_lambda_ok, [LAMBDA_LJ, LAMBDA_C, ALPHA_LJ] are legal    *
c     * else [LAMBDA_LJ, LAMBDA_C] legal.                              *
c     * prints message and exits on error.                             *
c     * fname is a string specifying calling function/procedure        *
c     *                                                                *
c     * originally written 04-25-2006 by MAW                           *
c     * last update 10-13-2008 by M.G. Martin                          *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the function
      logical non_lambda_ok
      integer param, l
      character*(*) fname
c     --- local variables
      integer scp_num_foreign_lambda

      if (l .ne. -1) then
         call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
         if (l .gt. scp_num_foreign_lambda) then
            write(6,*) "Requested index ",l," greater than number of ",
     &         "foreign lambda (",scp_num_foreign_lambda,")"
            stop
         endif
      endif

      if (param .eq. LAMBDA_C) return
      if (param .eq. LAMBDA_LJ) return
      if (param .eq. ALPHA_LJ .and. non_lambda_ok) return

      write(6,*) 'Illegal parameter passed to ',fname
      stop
      end

