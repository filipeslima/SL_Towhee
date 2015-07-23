/*     $Id: wrapperC.h,v 1.6 2010/09/04 14:16:50 lperi Exp $ */
/* 
* MCCCS - Towhee: A Monte Carlo molecular simulation program     *
* Copyright (C) 2003-2005 Marcus G. Martin                       *
* see the file license.gpl for the full license information      *
*                                                                *
* This program is free software; you can redistribute it and/or  *
* modify it under the terms of the GNU General Public License    *
* as published by the Free Software Foundation; either version 2 *
* of the License, or (at your option) any later version.         *
*                                                                *
* This program is distributed in the hope that it will be useful,*
* but WITHOUT ANY WARRANTY; without even the implied warranty of *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
* GNU General Public License for more details.                   *
*                                                                *
* You should have received a copy of the GNU General Public      *
* License along with this program; if not, write to the Free     *
* Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,*
* MA  02111-1307, USA.                                           
*/
/* 
   Written 09-21-2005 by MAW
   last modified 07-02-2009 by M.G. Martin
*/
/*
 * definitions for data passed between towhee.c and mainloop.F
 * This is the C version
 * See, http://consult.stanford.edu/pub/programming/fortran/fortran-c
 */

#include "preproc.h"

/* TODO: separate data passed which overrides data in towhee_input (e.g.
 * lambda_lj) from data passed which has no analogy in towhee_input (and whose
 * visibility is not controlled by lreadwrapper in mainloop.F), such as
 * towhee_input_file.
 * last modified 08-16-2011 by M.G. Martin
 */
typedef struct {
    double lambda_lj, lambda_c;
    int num_steps, num_foreign_lambda, random_seed;
    int linit, lredirect_stdout;
} wrapper;

/* strings must be in a separate struct */
typedef struct {
    char towhee_input_file[MAXDIRLENGTH];
    char output_dir[MAXDIRLENGTH];
} str_wrapper;

extern wrapper wrap_;
#ifdef INTEL_VISUAL_FORTRAN
#define swrap_ SWRAP
#endif
extern str_wrapper swrap_;
