c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
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
c     * last modified 08-16-2011 by M.G. Martin                        *
c     ******************************************************************
c 
c  * definitions for data passed between towhee.c and mainloop.F
c  * This is the Fortran version
c  * See, http://consult.stanford.edu/pub/programming/fortran/fortran-c
c  
c  It is assumed that preproc.h is already included by the time this is
c  called
c
c  The common blocks shared with the C wrapper are /wrap/ and /swrap/
      double precision wrap_lambda_lj,wrap_lambda_c
      integer wrap_nstep,wrap_num_foreign_lambda,wrap_random_seed
      logical wrap_linit, wrap_lredirect_stdout

      common /wrap/ wrap_lambda_lj,wrap_lambda_c
     &,wrap_nstep
     &,wrap_num_foreign_lambda,wrap_random_seed
     &,wrap_linit, wrap_lredirect_stdout

c     --- strings must be in a separate common block
      character*MAXDIRLENGTH wrap_towhee_input_file
      character*MAXDIRLENGTH wrap_output_dir

      common /swrap/ wrap_towhee_input_file,wrap_output_dir

