#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2014 Marcus G. Martin                       *
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
      subroutine twh_ffdevel(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nimprotypes)
c     ******************************************************************
c     * holding place for some force field parameters that I'm trying  *
c     * out for various projects.  These could change at any time and  *
c     * I don't suggest using them unless you are helping to work on   *
c     * these parameterizations.                                       *
c     *                                                                *
c     * originally written 02-14-2005 by M.G. Martin                   *
c     * last modified 07-15-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nimprotypes
c     --- local variables
      character*(FFNAMELEN) ffname

c     --- set variables for this forcefield
      ffname = 'devel'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'LB or Geometric')

      call twh_ffdevel_nonbond(lfailure,ffname,ntypes)
      if ( lfailure ) return
      call twh_ffdevel_bond(lfailure,ffname,nbondtypes)
      if ( lfailure ) return
      call twh_ffdevel_angle(lfailure,ffname,nangletypes)
      if ( lfailure ) return
      call twh_ffdevel_dihedral(lfailure,ffname,ndihedtypes)
      if ( lfailure ) return
      call twh_ffdevel_improper(lfailure,ffname,nimprotypes)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffdevel_nonbond(lfailure,ffname,ntypes)
c     ******************************************************************
c     * originally written 04-29-2008 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      character*(*) ffname
      integer ntypes
c     --- local variables
      integer ii,type,ielement
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal

      type = 0
c     --- ionic liquids
c     --- 1-Ethyl-3-methylimidazolium
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.80d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.050d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 'sp2'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CPH1) ring carbon bonded to a single N in imidazolium
      nbnlocal ='CPH1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.80d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.050d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 'sp2'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CPH2) ring carbon bonded to two N in imidazolium
      nbnlocal ='CPH2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.85d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.200d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 'sp2'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' N'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NR1) ring nitrogen in imidazolium
      nbnlocal ='NR1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.175d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.055d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CT3I) methyl group on the ethyl branch in imidazolium
      nbnlocal ='CT3I'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.275d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.02d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN7B) side-chain carbon bonded to N in imidazolium
      nbnlocal ='CN7B'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.90d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.046d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' H'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HR1) bonded to CPH2
      nbnlocal ='HR1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.468d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0078d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' H'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HR3) bonded to CPH1
      nbnlocal ='HR3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.32d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.022d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' H'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HA) bonded to CT3I
      nbnlocal ='HA'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.32d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.022d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' H'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HN7) bonded to CN7B
      nbnlocal ='HN7'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- butanesulfonate
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.34d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.024d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' H'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HA3) bonded to CT3
      nbnlocal ='HA3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.040d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.078d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CT3) methyl group carbon
      nbnlocal ='CT3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.010d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.056d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CT2) methylene group carbon not bonded to S
      nbnlocal ='CT2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.34d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.028d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' H'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HA2) bonded to CT2
      nbnlocal ='HA2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.1d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.47d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' S'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (SBS) in sulfonate
      nbnlocal ='SBS'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.7d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.12d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' O'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (OBS) in sulfonate
      nbnlocal ='OBS'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.010d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.056d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CBS) carbon bonded to S in sulfonate
      nbnlocal ='CBS'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.340d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.028d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' H'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HBS) bonded to the carbon that is bonded to S in sulfonate
      nbnlocal ='HBS'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.09900d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.200d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' P'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (P3) P in phosphonium cation
      nbnlocal ='P3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.96431d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (C1P) alkane sidechain carbons bonded to P in phosphonium cation
      nbnlocal ='C1P'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.96431d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (C1PR) alkane ring carbons bonded to P in phosphonium cation
      nbnlocal ='C1PR'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.96431d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (C2R) alkane ring carbons not bonded to P in phosphonium cation
      nbnlocal ='C2R'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.96431d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CT) alkane sidechain carbons not bonded to P in phosphonium
c     cation
      nbnlocal ='CT'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.40308d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.03d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' H'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (H1) H bonded to C bonded to P in phosphonium
c     cation
      nbnlocal ='H1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.40308d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.03d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' H'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HC) H bonded to C not bonded to P in phosphonium
c     cation
      nbnlocal ='HC'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.99d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 'sp2'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CCM) carbon in the center of tricyanomethanide
      nbnlocal ='CCM'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.75d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.20d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 'sp'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN) carbon triple bonded to N in tricyanomethanide
      nbnlocal ='CN'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.85d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.60d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      elelocal = ' N'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_getelement(lfailure,'getmass',masslocal,elelocal
     &     ,ielement )
      if ( lfailure ) return
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NC) N in tricyanomethanide
      nbnlocal ='NC'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffdevel_bond(lfailure,ffname,nbondtypes)
c     ******************************************************************
c     * bond lengths - Harmonic                                        *
c     *                                                                *
c     * originally written 04-29-2008 by M.G. Martin                   *
c     * last modified 07-15-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(*) ffname
      logical lfailure
      integer nbondtypes
c     --- local variables
      integer type,isame
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 

      type = 0
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.361d0
      vclocal(1) = (410.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CPH1 - CPH1
      isame = isame + 1
      vnlocal(isame,1) = 'CPH1'
      vnlocal(isame,2) = 'CPH1'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.337d0
      vclocal(1) = (400.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NR1 - CPH2
      isame = isame + 1
      vnlocal(isame,1) = 'NR1'
      vnlocal(isame,2) = 'CPH2'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.382d0
      vclocal(1) = (400.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NR1 - CPH1
      isame = isame + 1
      vnlocal(isame,1) = 'NR1'
      vnlocal(isame,2) = 'CPH1'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.477d0
      vclocal(1) = (220.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NR1 - CN7B
      isame = isame + 1
      vnlocal(isame,1) = 'NR1'
      vnlocal(isame,2) = 'CN7B'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.077d0
      vclocal(1) = (365.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CPH1 - HR3
      isame = isame + 1
      vnlocal(isame,1) = 'CPH1'
      vnlocal(isame,2) = 'HR3'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.078d0
      vclocal(1) = (340.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CPH2 - HR1
      isame = isame + 1
      vnlocal(isame,1) = 'CPH2'
      vnlocal(isame,2) = 'HR1'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.090d0
      vclocal(1) = (309.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7B - HN7
      isame = isame + 1
      vnlocal(isame,1) = 'CN7B'
      vnlocal(isame,2) = 'HN7'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.524d0
      vclocal(1) = (200.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7B - CT3I
      isame = isame + 1
      vnlocal(isame,1) = 'CN7B'
      vnlocal(isame,2) = 'CT3I'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.092d0
      vclocal(1) = (322.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CT3I - HA
      isame = isame + 1
      vnlocal(isame,1) = 'CT3I'
      vnlocal(isame,2) = 'HA'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.093d0
      vclocal(1) = (346.60d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CBS - HBS
      isame = isame + 1
      vnlocal(isame,1) = 'CBS'
      vnlocal(isame,2) = 'HBS'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.842d0
      vclocal(1) = (178.60d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CBS - SBS
      isame = isame + 1
      vnlocal(isame,1) = 'CBS'
      vnlocal(isame,2) = 'SBS'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.492d0
      vclocal(1) = (528.50d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- SBS - OBS
      isame = isame + 1
      vnlocal(isame,1) = 'SBS'
      vnlocal(isame,2) = 'OBS'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.523d0
      vclocal(1) = (298.50d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CBS - CT2
      isame = isame + 1
      vnlocal(isame,1) = 'CBS'
      vnlocal(isame,2) = 'CT2'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.096d0
      vclocal(1) = (331.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CT2 - HA2
      isame = isame + 1
      vnlocal(isame,1) = 'CT2'
      vnlocal(isame,2) = 'HA2'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.534d0
      vclocal(1) = (287.8d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CT2 - CT2
      isame = isame + 1
      vnlocal(isame,1) = 'CT2'
      vnlocal(isame,2) = 'CT2'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.532d0
      vclocal(1) = (292.5d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CT2 - CT3
      isame = isame + 1
      vnlocal(isame,1) = 'CT2'
      vnlocal(isame,2) = 'CT3'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.095d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CT3 - HA3
      isame = isame + 1
      vnlocal(isame,1) = 'CT3'
      vnlocal(isame,2) = 'HA3'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.823d0
      vclocal(1) = (206.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- P3 - C1P
      isame = isame + 1
      vnlocal(isame,1) = 'P3'
      vnlocal(isame,2) = 'C1P'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.839d0
      vclocal(1) = (206.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- P3 - C1PR
      isame = isame + 1
      vnlocal(isame,1) = 'P3'
      vnlocal(isame,2) = 'C1PR'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.547d0
      vclocal(1) = (281.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- C1PR - C2R
      isame = isame + 1
      vnlocal(isame,1) = 'C1PR'
      vnlocal(isame,2) = 'C2R'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.538d0
      vclocal(1) = (281.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- C2R - C2R
      isame = isame + 1
      vnlocal(isame,1) = 'C2R'
      vnlocal(isame,2) = 'C2R'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.534d0
      vclocal(1) = (281.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- C1P - CT
      isame = isame + 1
      vnlocal(isame,1) = 'C1P'
      vnlocal(isame,2) = 'CT'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.092d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- C1P - H1
      isame = isame + 1
      vnlocal(isame,1) = 'C1P'
      vnlocal(isame,2) = 'H1'
c     --- C1PR - H1
      isame = isame + 1
      vnlocal(isame,1) = 'C1PR'
      vnlocal(isame,2) = 'H1'
c     --- C2R - HC
      isame = isame + 1
      vnlocal(isame,1) = 'C2R'
      vnlocal(isame,2) = 'HC'
c     --- CT - HC
      isame = isame + 1
      vnlocal(isame,1) = 'CT'
      vnlocal(isame,2) = 'HC'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.408d0
      vclocal(1) = (430.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CCM - CN
      isame = isame + 1
      vnlocal(isame,1) = 'CCM'
      vnlocal(isame,2) = 'CN'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.167d0
      vclocal(1) = (1210.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- CN - NC
      isame = isame + 1
      vnlocal(isame,1) = 'CN'
      vnlocal(isame,2) = 'NC'
c     --- set vibnames
      call twh_set_vibnames(type,isame,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffdevel_angle(lfailure,ffname,nangletypes)
c     ******************************************************************
c     * angles - harmonic                                              *
c     *                                                                *
c     * originally written 04-29-2008 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(*) ffname
      logical lfailure
      integer nangletypes
c     --- local variables
      integer type,isame
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)

      type = 0
      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.18d0
      bclocal(1) = (130.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NR1 - CPH1 - CPH1
      isame = isame + 1
      bnlocal(isame,1) = 'NR1'
      bnlocal(isame,2) = 'CPH1'
      bnlocal(isame,3) = 'CPH1'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.28d0
      bclocal(1) = (130.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CPH1 - NR1 - CPH2
      isame = isame + 1
      bnlocal(isame,1) = 'CPH1'
      bnlocal(isame,2) = 'NR1'
      bnlocal(isame,3) = 'CPH2'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.08d0
      bclocal(1) = (130.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NR1 - CPH2 - NR1
      isame = isame + 1
      bnlocal(isame,1) = 'NR1'
      bnlocal(isame,2) = 'CPH2'
      bnlocal(isame,3) = 'NR1'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.88d0
      bclocal(1) = (130.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CPH2 - NR1 - CN7B
      isame = isame + 1
      bnlocal(isame,1) = 'CPH2'
      bnlocal(isame,2) = 'NR1'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.83d0
      bclocal(1) = (130.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CPH1 - NR1 - CN7B
      isame = isame + 1
      bnlocal(isame,1) = 'CPH1'
      bnlocal(isame,2) = 'NR1'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.33d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN7 - CN7B - NR1
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'NR1'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.97d0
      bclocal(1) = (35.5d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN7 - CN7B - HN7
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'HN7'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.46d0
      bclocal(1) = (25.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HR1 - CPH2 - NR1
      isame = isame + 1
      bnlocal(isame,1) = 'HR1'
      bnlocal(isame,2) = 'CPH2'
      bnlocal(isame,3) = 'NR1'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.09d0
      bclocal(1) = (25.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HR3 - CPH1 - NR1
      isame = isame + 1
      bnlocal(isame,1) = 'HR3'
      bnlocal(isame,2) = 'CPH1'
      bnlocal(isame,3) = 'NR1'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 130.73d0
      bclocal(1) = (25.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HR3 - CPH1 - CPH1
      isame = isame + 1
      bnlocal(isame,1) = 'HR3'
      bnlocal(isame,2) = 'CPH1'
      bnlocal(isame,3) = 'CPH1'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.69d0
      bclocal(1) = (33.4d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HN7 - CN7B - CT3I
      isame = isame + 1
      bnlocal(isame,1) = 'HN7'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'CT3I'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.86d0
      bclocal(1) = (33.4d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA - CT3I - CN7B
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT3I'
      bnlocal(isame,3) = 'CN7B'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.03d0
      bclocal(1) = (35.5d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA - CT3I - HA
      isame = isame + 1
      bnlocal(isame,1) = 'HA'
      bnlocal(isame,2) = 'CT3I'
      bnlocal(isame,3) = 'HA'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.44d0
      bclocal(1) = (140.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NR1 - CN7B - CT3I
      isame = isame + 1
      bnlocal(isame,1) = 'NR1'
      bnlocal(isame,2) = 'CN7B'
      bnlocal(isame,3) = 'CT3I'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 105.69d0
      bclocal(1) = (37.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HBS - CBS - SBS
      isame = isame + 1
      bnlocal(isame,1) = 'HBS'
      bnlocal(isame,2) = 'CBS'
      bnlocal(isame,3) = 'SBS'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.55d0
      bclocal(1) = (33.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HBS - CBS - HBS
      isame = isame + 1
      bnlocal(isame,1) = 'HBS'
      bnlocal(isame,2) = 'CBS'
      bnlocal(isame,3) = 'HBS'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.35d0
      bclocal(1) = (74.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OBS - SBS - CBS
      isame = isame + 1
      bnlocal(isame,1) = 'OBS'
      bnlocal(isame,2) = 'SBS'
      bnlocal(isame,3) = 'CBS'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.07d0
      bclocal(1) = (95.6d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OBS - SBS - OBS
      isame = isame + 1
      bnlocal(isame,1) = 'OBS'
      bnlocal(isame,2) = 'SBS'
      bnlocal(isame,3) = 'OBS'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.97d0
      bclocal(1) = (61.6d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- SBS - CBS - CT2
      isame = isame + 1
      bnlocal(isame,1) = 'SBS'
      bnlocal(isame,2) = 'CBS'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.77d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HBS - CBS - CT2
      isame = isame + 1
      bnlocal(isame,1) = 'HBS'
      bnlocal(isame,2) = 'CBS'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.79d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA2 - CT2 - CBS
      isame = isame + 1
      bnlocal(isame,1) = 'HA2'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CBS'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.20d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA2 - CT2 - HA2
      isame = isame + 1
      bnlocal(isame,1) = 'HA2'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'HA2'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.59d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA2 - CT2 - CT2
      isame = isame + 1
      bnlocal(isame,1) = 'HA2'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.23d0
      bclocal(1) = (53.4d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CBS - CT2 - CT2
      isame = isame + 1
      bnlocal(isame,1) = 'CBS'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.33d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA2 - CT2 - CT3
      isame = isame + 1
      bnlocal(isame,1) = 'HA2'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT3'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.26d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA3 - CT3 - CT2
      isame = isame + 1
      bnlocal(isame,1) = 'HA3'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.62d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HA3 - CT3 - HA3
      isame = isame + 1
      bnlocal(isame,1) = 'HA3'
      bnlocal(isame,2) = 'CT3'
      bnlocal(isame,3) = 'HA3'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.14d0
      bclocal(1) = (53.4d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CT3 - CT2 - CT2
      isame = isame + 1
      bnlocal(isame,1) = 'CT3'
      bnlocal(isame,2) = 'CT2'
      bnlocal(isame,3) = 'CT2'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 97.1d0
      bclocal(1) = (58.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- C1PR - P3 - C1PR
      isame = isame + 1
      bnlocal(isame,1) = 'C1PR'
      bnlocal(isame,2) = 'P3'
      bnlocal(isame,3) = 'C1PR'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.6d0
      bclocal(1) = (58.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- C1P - P3 - C1P
      isame = isame + 1
      bnlocal(isame,1) = 'C1P'
      bnlocal(isame,2) = 'P3'
      bnlocal(isame,3) = 'C1P'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.8d0
      bclocal(1) = (58.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- C1PR - P3 - C1P
      isame = isame + 1
      bnlocal(isame,1) = 'C1P'
      bnlocal(isame,2) = 'P3'
      bnlocal(isame,3) = 'C1PR'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.2d0
      bclocal(1) = (55.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- P3 - C1PR - C2R
      isame = isame + 1
      bnlocal(isame,1) = 'P3'
      bnlocal(isame,2) = 'C1PR'
      bnlocal(isame,3) = 'C2R'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.8d0
      bclocal(1) = (42.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- P3 - C1PR - H1
      isame = isame + 1
      bnlocal(isame,1) = 'P3'
      bnlocal(isame,2) = 'C1PR'
      bnlocal(isame,3) = 'H1'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.40d0
      bclocal(1) = (42.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- P3 - C1P - H1
      isame = isame + 1
      bnlocal(isame,1) = 'P3'
      bnlocal(isame,2) = 'C1P'
      bnlocal(isame,3) = 'H1'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.2d0
      bclocal(1) = (55.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- P3 - C1P - CT
      isame = isame + 1
      bnlocal(isame,1) = 'P3'
      bnlocal(isame,2) = 'C1P'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.9d0
      bclocal(1) = (67.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- C1PR - C2R - C2R
      isame = isame + 1
      bnlocal(isame,1) = 'C1PR'
      bnlocal(isame,2) = 'C2R'
      bnlocal(isame,3) = 'C2R'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.6d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- H1 - C1PR - C2R
      isame = isame + 1
      bnlocal(isame,1) = 'H1'
      bnlocal(isame,2) = 'C1PR'
      bnlocal(isame,3) = 'C2R'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.7d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- H1 - C1PR - H1
      isame = isame + 1
      bnlocal(isame,1) = 'H1'
      bnlocal(isame,2) = 'C1PR'
      bnlocal(isame,3) = 'H1'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.5d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- H1 - C1P - H1
      isame = isame + 1
      bnlocal(isame,1) = 'H1'
      bnlocal(isame,2) = 'C1P'
      bnlocal(isame,3) = 'H1'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.1d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HC - C2R - C1PR
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'C2R'
      bnlocal(isame,3) = 'C1PR'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.9d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HC - C2R - C2R
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'C2R'
      bnlocal(isame,3) = 'C2R'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.0d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HC - C2R - HC
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'C2R'
      bnlocal(isame,3) = 'HC'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.0d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- H1 - C1P - CT
      isame = isame + 1
      bnlocal(isame,1) = 'H1'
      bnlocal(isame,2) = 'C1P'
      bnlocal(isame,3) = 'CT'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.1d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HC - CT - C1P
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'C1P'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.8d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- HC - CT - HC
      isame = isame + 1
      bnlocal(isame,1) = 'HC'
      bnlocal(isame,2) = 'CT'
      bnlocal(isame,3) = 'HC'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 179.9999d0
      bclocal(1) = (30.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- NC - CN - CCM
      isame = isame + 1
      bnlocal(isame,1) = 'NC'
      bnlocal(isame,2) = 'CN'
      bnlocal(isame,3) = 'CCM'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
c     --- standard harmonic (no urey-bradly)
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (44.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- CN - CCM - CN
      isame = isame + 1
      bnlocal(isame,1) = 'CN'
      bnlocal(isame,2) = 'CCM'
      bnlocal(isame,3) = 'CN'
c     --- set nbsame
      call twh_set_bendnames(type,isame,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffdevel_dihedral(lfailure,ffname,ndihedtypes)
c     ******************************************************************
c     * torsions - Charmm style                                        *
c     *                                                                *
c     * originally written 04-29-2008 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(*) ffname
      logical lfailure
      integer ndihedtypes
c     --- local variables
      integer type,isame
      double precision onepi
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)

c     --- fetch onepi
      call twh_constant_pi(GLB_GET,onepi)

      type = 0
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (14.0d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- NR1 - CPH1 - CPH1 - NR1
      isame = isame + 1
      tnlocal(isame,1) = 'NR1'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'NR1'
c     --- CPH1 - NR1 - CPH2 - NR1
      isame = isame + 1
      tnlocal(isame,1) = 'CPH1'
      tnlocal(isame,2) = 'NR1'
      tnlocal(isame,3) = 'CPH2'
      tnlocal(isame,4) = 'NR1'
c     --- CPH1 - CPH1 - NR1 - CPH2
      isame = isame + 1
      tnlocal(isame,1) = 'CPH1'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CPH2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (3.0d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HR3 - CPH1 - CPH1 - NR1
      isame = isame + 1
      tnlocal(isame,1) = 'HR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'NR1'
c     --- CPH1 - NR1 - CPH2 - HR1
      isame = isame + 1
      tnlocal(isame,1) = 'CPH1'
      tnlocal(isame,2) = 'NR1'
      tnlocal(isame,3) = 'CPH2'
      tnlocal(isame,4) = 'HR1'
c     --- HR3 - CPH1 - NR1 - CPH2
      isame = isame + 1
      tnlocal(isame,1) = 'HR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CPH2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (2.0d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HR3 - CPH1 - CPH1 - HR3
      isame = isame + 1
      tnlocal(isame,1) = 'HR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'CPH1'
      tnlocal(isame,4) = 'HR3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
c     --- basically a null torsion, plus 1-4 terms
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CN7B - NR1 - CPH2 - NR1
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NR1'
      tnlocal(isame,3) = 'CPH2'
      tnlocal(isame,4) = 'NR1'
c     --- CN7B - NR1 - CPH2 - HR1
      isame = isame + 1
      tnlocal(isame,1) = 'CN7B'
      tnlocal(isame,2) = 'NR1'
      tnlocal(isame,3) = 'CPH2'
      tnlocal(isame,4) = 'HR1'
c     --- CPH1 - CPH1 - NR1 - CN7B
      isame = isame + 1
      tnlocal(isame,1) = 'CPH1'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CN7B'
c     --- HR3 - CPH1 - NR1 - CN7B
      isame = isame + 1
      tnlocal(isame,1) = 'HR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CN7B'
c     --- HR3 - CPH1 - NR1 - CN7B
      isame = isame + 1
      tnlocal(isame,1) = 'HR3'
      tnlocal(isame,2) = 'CPH1'
      tnlocal(isame,3) = 'NR1'
      tnlocal(isame,4) = 'CN7B'
c     --- CPH1 - NR1 - CN7B - HN7
      isame = isame + 1
      tnlocal(isame,1) = 'CPH1'
      tnlocal(isame,2) = 'NR1'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'HN7'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1950d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CPH2 - NR1 - CN7B - HN7
      isame = isame + 1
      tnlocal(isame,1) = 'CPH2'
      tnlocal(isame,2) = 'NR1'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'HN7'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.10d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CPH2 - NR1 - CN7B - CT3I
      isame = isame + 1
      tnlocal(isame,1) = 'CPH2'
      tnlocal(isame,2) = 'NR1'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'CT3I'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.20d0)*(CON_KCALTOK)
      tclocal(2) = 4.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CPH1 - NR1 - CN7B - CT3I
      isame = isame + 1
      tnlocal(isame,1) = 'CPH1'
      tnlocal(isame,2) = 'NR1'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'CT3I'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
c     --- another null with a 3 mult instead of 1.  I do not think
c     --- that could make any difference, but just to be consistent
c     --- with the info here it is
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.00d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HA - CT3I - CN7B - NR1
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CT3I'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'NR1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1950d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HA - CT3I - CN7B - HN7
      isame = isame + 1
      tnlocal(isame,1) = 'HA'
      tnlocal(isame,2) = 'CT3I'
      tnlocal(isame,3) = 'CN7B'
      tnlocal(isame,4) = 'HN7'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1758d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- OBS - SBS - CBS - HBS
      isame = isame + 1
      tnlocal(isame,1) = 'OBS'
      tnlocal(isame,2) = 'SBS'
      tnlocal(isame,3) = 'CBS'
      tnlocal(isame,4) = 'HBS'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1963d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- OBS - SBS - CBS - CT2
      isame = isame + 1
      tnlocal(isame,1) = 'OBS'
      tnlocal(isame,2) = 'SBS'
      tnlocal(isame,3) = 'CBS'
      tnlocal(isame,4) = 'CT2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1540d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HA2 - CT2 - CBS - HBS
      isame = isame + 1
      tnlocal(isame,1) = 'HA2'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CBS'
      tnlocal(isame,4) = 'HBS'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1021d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HA2 - CT2 - CBS - SBS
      isame = isame + 1
      tnlocal(isame,1) = 'HA2'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CBS'
      tnlocal(isame,4) = 'SBS'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1842d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CT2 - CT2 - CBS - HBS
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CBS'
      tnlocal(isame,4) = 'HBS'
c     --- HA3 - CT3 - CT2 - CT2
      isame = isame + 1
      tnlocal(isame,1) = 'HA3'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CT2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = (-1.6288d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(3+1) = (-0.2844d0)*(CON_KCALTOK)
      tclocal(3+2) = 2.0d0
      tclocal(3+3) = (180.0d0)*(onepi/180.0d0)
      tclocal(6+1) = (0.1040d0)*(CON_KCALTOK)
      tclocal(6+2) = 3.0d0
      tclocal(6+3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CT2 - CT2 - CBS - SBS
      isame = isame + 1
      tnlocal(isame,1) = 'CT2'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CBS'
      tnlocal(isame,4) = 'SBS'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1540d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HA2 - CT2 - CT2 - HA2
      isame = isame + 1
      tnlocal(isame,1) = 'HA2'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HA2'
c     --- HA3 - CT3 - CT2 - HA2
      isame = isame + 1
      tnlocal(isame,1) = 'HA3'
      tnlocal(isame,2) = 'CT3'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HA2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1842d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HA2 - CT2 - CT2 - CBS
      isame = isame + 1
      tnlocal(isame,1) = 'HA2'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CBS'
c     --- CT3 - CT2 - CT2 - HA2
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'HA2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,4)
      tclocal(1) = (0.5218d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(3+1) = (-0.0826d0)*(CON_KCALTOK)
      tclocal(3+2) = 2.0d0
      tclocal(3+3) = (180.0d0)*(onepi/180.0d0)
      tclocal(6+1) = (0.2002d0)*(CON_KCALTOK)
      tclocal(6+2) = 3.0d0
      tclocal(6+3) = (0.0d0)*(onepi/180.0d0)
      tclocal(9+1) = (-0.1348d0)*(CON_KCALTOK)
      tclocal(9+2) = 4.0d0
      tclocal(9+3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CT3 - CT2 - CT2 - CBS
      isame = isame + 1
      tnlocal(isame,1) = 'CT3'
      tnlocal(isame,2) = 'CT2'
      tnlocal(isame,3) = 'CT2'
      tnlocal(isame,4) = 'CBS'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1354d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- C2R - C1PR - P3 - C1PR
      isame = isame + 1
      tnlocal(isame,1) = 'C2R'
      tnlocal(isame,2) = 'C1PR'
      tnlocal(isame,3) = 'P3'
      tnlocal(isame,4) = 'C1PR'
c     --- C2R - C1PR - P3 - C1P
      isame = isame + 1
      tnlocal(isame,1) = 'C2R'
      tnlocal(isame,2) = 'C1PR'
      tnlocal(isame,3) = 'P3'
      tnlocal(isame,4) = 'C1P'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1108d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- H1 - C1PR - P3 - C1PR
      isame = isame + 1
      tnlocal(isame,1) = 'H1'
      tnlocal(isame,2) = 'C1PR'
      tnlocal(isame,3) = 'P3'
      tnlocal(isame,4) = 'C1PR'
c     --- H1 - C1PR - P3 - C1P
      isame = isame + 1
      tnlocal(isame,1) = 'H1'
      tnlocal(isame,2) = 'C1PR'
      tnlocal(isame,3) = 'P3'
      tnlocal(isame,4) = 'C1P'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = (-0.3882d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(3+1) = (0.1181d0)*(CON_KCALTOK)
      tclocal(3+2) = 2.0d0
      tclocal(3+3) = (180.0d0)*(onepi/180.0d0)
      tclocal(6+1) = (-0.0854d0)*(CON_KCALTOK)
      tclocal(6+2) = 3.0d0
      tclocal(6+3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- C2R - C2R - C1PR - P3
      isame = isame + 1
      tnlocal(isame,1) = 'C2R'
      tnlocal(isame,2) = 'C2R'
      tnlocal(isame,3) = 'C1PR'
      tnlocal(isame,4) = 'P3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1830d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- C2R - C2R - C1PR - H1
      isame = isame + 1
      tnlocal(isame,1) = 'C2R'
      tnlocal(isame,2) = 'C2R'
      tnlocal(isame,3) = 'C1PR'
      tnlocal(isame,4) = 'H1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0556d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HC - C2R - C1PR - P3
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'C2R'
      tnlocal(isame,3) = 'C1PR'
      tnlocal(isame,4) = 'P3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1590d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HC - C2R - C1PR - H1
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'C2R'
      tnlocal(isame,3) = 'C1PR'
      tnlocal(isame,4) = 'H1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,3)
      tclocal(1) = (0.8700d0)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      tclocal(3+1) = (-0.0785d0)*(CON_KCALTOK)
      tclocal(3+2) = 2.0d0
      tclocal(3+3) = (180.0d0)*(onepi/180.0d0)
      tclocal(6+1) = (0.1395d0)*(CON_KCALTOK)
      tclocal(6+2) = 3.0d0
      tclocal(6+3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- C1PR - C2R - C2R - C1PR
      isame = isame + 1
      tnlocal(isame,1) = 'C1PR'
      tnlocal(isame,2) = 'C2R'
      tnlocal(isame,3) = 'C2R'
      tnlocal(isame,4) = 'C1PR'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1830d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HC - C2R - C2R - C1PR
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'C2R'
      tnlocal(isame,3) = 'C2R'
      tnlocal(isame,4) = 'C1PR'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1590d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HC - C2R - C2R - HC
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'C2R'
      tnlocal(isame,3) = 'C2R'
      tnlocal(isame,4) = 'HC'
c     --- HC - CT - C1P - H1
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'CT'
      tnlocal(isame,3) = 'C1P'
      tnlocal(isame,4) = 'H1'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1354d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CT - C1P - P3 - C1P
      isame = isame + 1
      tnlocal(isame,1) = 'CT'
      tnlocal(isame,2) = 'C1P'
      tnlocal(isame,3) = 'P3'
      tnlocal(isame,4) = 'C1P'
c     --- CT - C1P - P3 - C1PR
      isame = isame + 1
      tnlocal(isame,1) = 'CT'
      tnlocal(isame,2) = 'C1P'
      tnlocal(isame,3) = 'P3'
      tnlocal(isame,4) = 'C1PR'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1108d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- H1 - C1P - P3 - C1P
      isame = isame + 1
      tnlocal(isame,1) = 'H1'
      tnlocal(isame,2) = 'C1P'
      tnlocal(isame,3) = 'P3'
      tnlocal(isame,4) = 'C1P'
c     --- H1 - C1P - P3 - C1PR
      isame = isame + 1
      tnlocal(isame,1) = 'H1'
      tnlocal(isame,2) = 'C1P'
      tnlocal(isame,3) = 'P3'
      tnlocal(isame,4) = 'C1PR'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0556d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- HC - CT - C1P - P3
      isame = isame + 1
      tnlocal(isame,1) = 'HC'
      tnlocal(isame,2) = 'CT'
      tnlocal(isame,3) = 'C1P'
      tnlocal(isame,4) = 'P3'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- NC - CN - CCM - CN
      isame = isame + 1
      tnlocal(isame,1) = 'NC'
      tnlocal(isame,2) = 'CN'
      tnlocal(isame,3) = 'CCM'
      tnlocal(isame,4) = 'CN'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffdevel_improper(lfailure,ffname,nimprotypes)
c     ******************************************************************
c     * improper torsions                                              *
c     *                                                                *
c     * originally written 04-29-2008 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(*) ffname
      logical lfailure
      integer nimprotypes
c     --- local variables
      integer type,isame
      character*(FFNAMELEN) ninlocal 
      dimension ninlocal(MAXSAMEIMP,MAXIMPNAME)
      double precision onepi
      double precision iclocal
      dimension iclocal(MINIMPCOEFF:MAXIMPCOEFF)

c     --- fetch onepi
      call twh_constant_pi(GLB_GET,onepi)

      type = 0
      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (0.50d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- CPH2 - NR1 - NR1 - HR1
      isame = isame + 1
      ninlocal(isame,1) = 'CPH2'
      ninlocal(isame,2) = 'NR1'
      ninlocal(isame,3) = 'NR1'
      ninlocal(isame,4) = 'HR1'
c     --- CPH1 - CPH1 - NR1 - HR3
      isame = isame + 1
      ninlocal(isame,1) = 'CPH1'
      ninlocal(isame,2) = 'CPH1'
      ninlocal(isame,3) = 'NR1'
      ninlocal(isame,4) = 'HR3'
c     --- set nimpsame
      call twh_set_impnames(type,isame,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (0.60d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- NR1 - CPH1 - CPH2 - CN7B
      isame = isame + 1
      ninlocal(isame,1) = 'NR1'
      ninlocal(isame,2) = 'CPH1'
      ninlocal(isame,3) = 'CPH2'
      ninlocal(isame,4) = 'CN7B'
c     --- set nimpsame
      call twh_set_impnames(type,isame,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (26.7d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- CCM - CN - CN - CN
      isame = isame + 1
      ninlocal(isame,1) = 'CCM'
      ninlocal(isame,2) = 'CN'
      ninlocal(isame,3) = 'CN'
      ninlocal(isame,4) = 'CN'
c     --- set nimpsame
      call twh_set_impnames(type,isame,ninlocal)

      nimprotypes = type
      call twh_ffcheck(lfailure,4,type)
      if ( lfailure ) return
      
      return
      end
